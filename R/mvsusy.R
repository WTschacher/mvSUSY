isTRUEorFALSE = function(x) {
  isTRUE(x) || isFALSE(x)
}
omega = function(x, ...) {
  Hact = det(abs(cov(x)))
  Hpot = prod(diag(cov(x)))
  1-(Hact/Hpot)
}
lambda_max = function(x, ...) {
  corx = cor(x)
  corx[is.na(corx)] = 0
  eigenwerte = eigen(corx)
  eigenv = eigenwerte$values
  max(eigenv) / (sum(eigenv)*.01)
}
eigenvalue = function(x) {
  eigenwerte = cor(x)
  eigenwerte[is.na(eigenwerte)] = 0
  eigen(eigenwerte)$values
}
as.mvsusy = function(x) {
  if (inherits(x, "mvsusy"))
    return(x)
  if (!is.list(x))
    stop("only list class objects can be turned into mvsusy objects")
  class(x) = unique(c("mvsusy", class(x)))
  x
}

mvsusy = function(x, segment, Hz,
                  method = c("lambda_max","omega"),
                  max_pseudo = 1000,
                  seed = 1) {
  ## 2.1 Warnings - check parameter settings
  if (!is.data.frame(x))
    stop("'x' must be a data.frame")
  if (is.null(names(x)))
    stop("'x' must have named columns")
  if (!all(vapply(x, is.numeric, FALSE, USE.NAMES=FALSE)))
    stop("'x' must have numeric columns")
  if (!is.numeric(segment) || length(segment)!=1L || is.na(segment))
    stop("'segment' must be scalar non-NA numeric")
  if (!is.numeric(Hz) || length(Hz)!=1L || is.na(Hz))
    stop("'Hz' must be scalar non-NA numeric")
  cols = names(x)
  N = length(cols)
  x = na.omit(x)
  data = x ## return for plot type time series
  nx = ncol(x)
  if (nx < 2L)
    stop("'x' must have at least 2 columns")
  segment = as.integer(segment)
  Hz = as.integer(Hz)
  method = match.arg(method)
  segmentHz = segment*Hz
  nsegment = floor(nrow(x)/segmentHz)
  if (segmentHz > ((nsegment*segmentHz) / ncol(x)))
    stop("Segment size is invalid: maximum segment size is ", (nsegment*segmentHz) / (ncol(x)*Hz))
  #  3.0 Real synchrony
  ## 3.1 Reshape data for segmenting
  if (nrow(x)%%segmentHz) { ## cut
    x = x[seq_len(nsegment*segmentHz),, drop=FALSE]
  }
  x$segment = rep(seq_len(nsegment), each=segmentHz)
  x$row = rep(seq_len(segmentHz), nsegment)
  xx = as.data.frame(dcast(as.data.table(x), row ~ segment, value.var=cols))
  xx$row = NULL
  ## 3.2 Generate list of real segment-combinations
  comb_real = paste0(cols, "_", rep(seq_len(nsegment), each=nx))
  comb_real = split(comb_real, ceiling(seq_along(comb_real)/nx))
  names(comb_real) = paste0('segment_combo_', seq_along(comb_real))
  ## 3.3 Expand data based on list of real segment-combinations
  df_real = lapply(comb_real, function(cols, x) x[cols], x=xx)
  ## 3.4 Define mv-synchrony function (based on method)
  method_fun = if (method=="lambda_max") lambda_max else if (method=="omega") omega else stop("internal error: unsupported method should ba cought by now")
  ## 3.5 Apply synchrony function per real segment-combination
  synchrony_real = sapply(df_real, method_fun)
  ## 3.6 Define eigenvalue function (based on method)
  ### done in global package namespace
  #  4.0 Pseudo synchrony
  ## 4.1 Generate list of pseudo segment-combinations
  err = character()
  tryCatch(comb_pseudo <- permuteSample(
    v = seq_len(nsegment), m = nx, n = max_pseudo, seed = seed,
    FUN = function(i) paste(cols, i, sep="_")
  ), error = function(e) err <<- e$message)
  if (length(err)) {
    if (err=="n exceeds the maximum number of possible results")
      stop("'max_pseudo' argument value exceeds the maximum number of permutations, decrease value of the argument")
    else
      stop("RcppAlgos::permuteSample returned error:\n", err)
  }
  ## 4.2 Expand data based on list of pseudo segment-combinations
  df_pseudo = lapply(comb_pseudo, function(cols, x) x[cols], x=xx)
  names(df_pseudo) = paste0("segment_combo_pseudo_", seq_along(df_pseudo))
  ## 4.3 Apply synchrony function per pseudo segment-combination
  synchrony_pseudo = sapply(df_pseudo, method_fun)
  #  5.0 Summarize synchrony data
  ## 5.1 Reshape synchrony data
  matrix_pseudo = data.frame(variable = "synchrony_pseudo", value=synchrony_pseudo)
  matrix_real = data.frame(variable = "synchrony_real", value=synchrony_real)
  synchrony = rbind(matrix_pseudo, matrix_real)
  rownames(synchrony) = NULL
  synchrony$variable = as.factor(synchrony$variable)
  ### 5.1.1 Parametric and nonparametric tests
  t_tests = t.test(value ~ variable, data = synchrony)
  wilcox_tests = wilcox.test(value ~ variable, data = synchrony)
  ### 5.1.2 Real and pseudo synchrony indices
  real_mean = mean(matrix_real$value)
  real_sd = sd(matrix_real$value)
  pseudo_mean = mean(matrix_pseudo$value)
  pseudo_sd = sd(matrix_pseudo$value)
  ### 5.1.3 Effect size
  ES_synchrony = (real_mean - pseudo_mean) / pseudo_sd
  ## 5.2 Generate eigenvalue data
  if (method=="lambda_max") {
    eigenvalue_real = sapply(df_real, eigenvalue)
    eigenvalue_pseudo = sapply(df_pseudo, eigenvalue)
    eigenvalue_real = data.frame(data="real", segment=rep(colnames(eigenvalue_real), each=nx), value=c(eigenvalue_real))
    eigenvalue_pseudo = data.frame(data="pseudo", segment=rep(colnames(eigenvalue_pseudo), each=nx), value=c(eigenvalue_pseudo))
    ev = rbind(eigenvalue_pseudo, eigenvalue_real)
    value = .N = NULL ## resolve NSE check notes only
    ev = as.data.frame(as.data.table(ev)[order(-value), "segment_num" := seq_len(.N), by=c("data","segment")])
  } else {
    ev = as.data.frame(NULL)
  }
  ans = list(
    method=method, n_col=nx, n_row=nrow(x), seed=seed, n_pseudo=length(comb_pseudo),
    segment_size_s = segment, data_per_segment = segmentHz*nx,
    real_mean=real_mean, real_sd=real_sd, pseudo_mean=pseudo_mean, pseudo_sd=pseudo_sd, ES_synchrony=ES_synchrony, EV=ev,
    t_tests = t_tests, wilcox_tests = wilcox_tests,
    nsegment = nsegment, max_pseudo = max_pseudo,
    synchrony = synchrony, data = data, segmentHz = segmentHz
  )
  as.mvsusy(ans)
}

as.data.frame.mvsusy = function(x, row.names=NULL, optional=FALSE, ...) {
  if (!inherits(x, "mvsusy"))
    stop("'x' must be an object of class 'mvsusy'")
  ans = data.frame(
    x$method, x$n_col, x$n_row, x$seed, x$n_pseudo, x$segment_size_s, x$data_per_segment,
    x$real_mean, x$real_sd, x$pseudo_mean, x$pseudo_sd, x$ES_synchrony,
    x$t_tests$statistic, x$t_tests$p.value,
    x$wilcox_tests$statistic, x$wilcox_tests$p.value
  )
  colnames(ans) = c(
    "method", "ncol", "nrow", "seed", "npseudo", "segment_size_s", "data_per_segment",
    "real_mean", "real_sd", "pseudo_mean", "pseudo_sd", "ES",
    "t_statistic", "p_value", "statistic_nonpar", "p_value_nonpar"
  )
  as.data.frame(ans)
}

print.mvsusy = function(x, ...) {
  if (!inherits(x, "mvsusy"))
    stop("'x' must be an object of class 'mvsusy'")
  df = as.data.frame(x)
  df$real_mean = round(df$real_mean,5)
  df$real_sd = round(df$real_sd,5)
  df$pseudo_mean = round(df$pseudo_mean,5)
  df$pseudo_sd = round(df$pseudo_sd,5)
  df$ES = round(df$ES,5)
  df$t_statistic = round(df$t_statistic,5)
  df$p_value = format(df$p_value, scientific=FALSE)
  df$statistic_nonpar = format(df$statistic_nonpar, scientific=FALSE)
  df$p_value_nonpar = format(df$p_value_nonpar, scientific=FALSE)
  colnames(df) = c(
    "method", "n(col)", "n(row)", "seed", "n(pseudo)", "segment size (s)", "data per segment",
    "mean(synchrony real)", "sd(synchrony real)", "mean(synchrony pseudo)", "sd(synchrony pseudo)", "ES",
    "t statistic", "p-value", "statistic nonpar", "p-value nonpar"
  )
  print(df, ...=...)
  invisible(x)
}

## plot various types of plots using mvsusy object
plot.mvsusy = function(x, type=c("eigenvalue","density","free scale","segment-wise","time series"), ..., plotly) {
  if (!inherits(x, "mvsusy"))
    stop("'x' must be an object of class 'mvsusy'")
  type = match.arg(type)
  if (!missing(plotly)) {
    if (!(identical(plotly, TRUE) || identical(plotly, FALSE)))
      stop("Argument 'plotly' must be TRUE or FALSE")
    if (isTRUE(plotly) && type!="time series")
      stop("Using 'plotly' is only implemented for time series type of mvSUSY plot")
  }
  # check NSE notes
  segment_num = segment = variable = . = mean.var = segment_id= NULL
  if (type=="eigenvalue") {
    if (x$method!="lambda_max")
      stop("plot mvSUSY of type eigenvalue is only for mvSUSY computed using 'lambda_max' method")
    p = ggplot(x$EV, aes(x = as.factor(segment_num), group = segment)) +
      geom_line(aes(y = value, colour = data)) +
      scale_color_manual(values = c('real' = 'chartreuse4', 'pseudo' = 'brown3'))+
      scale_y_continuous(breaks = seq(0, max(x$EV$value), by = 1)) +
      facet_wrap(~data, labeller = as_labeller(c(`pseudo` = paste( "n segments =", x$max_pseudo), `real` = paste( "n segments =", x$nsegment)))) +
      labs(y = "eigenvalue", x = "dimension",
           title = "eigenvalues per segment for real and pseudo data", subtitle = paste( "segmentsize =", x$segment_size_s)) +
      theme_light()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(), text = element_text(vjust = 0, size = 15, family="serif"), strip.text = element_text(size=15))
  } else if (type=="density") {
    p = ggplot(x$synchrony, aes(x = value, colour = variable))+
      geom_density(aes(x = value, fill = variable), alpha = .1, linewidth = 1, show.legend = FALSE)+
      geom_rug(aes(x = value, y = 0))+
      scale_colour_manual('data',
                          label = c('synchrony_pseudo'= "pseudo", 'synchrony_real'="real"),
                          values = c('synchrony_real' = 'chartreuse4', 'synchrony_pseudo' = 'brown3'))+
      scale_fill_manual('data',
                        values = c('synchrony_real' = 'chartreuse4', 'synchrony_pseudo' = 'brown3'))+
      labs(y = "density", x = "synchrony") +
      theme_light()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(vjust = 0, size = 12, family="serif"))+
      labs(title = "Dean")
  } else if (type=="free scale") {
    value = NULL ## fix check NSE notes
    vlines = as.data.frame(as.data.table(x$synchrony)[, .(mean.var=mean(value)), by="variable"])
    rm(value)
    p = ggplot(x$synchrony, aes(x = value, fill = variable))+
      geom_histogram(color="#e9ecef",alpha=0.8, position = 'identity')+ ## bins argument removed vs script as per email
      geom_vline(data = vlines, aes(xintercept = mean.var), linetype = "dashed", size = 0.8, alpha = 0.9)+
      facet_wrap(~variable, scales ="free", labeller = as_labeller(c('synchrony_pseudo'= paste("n segments =", x$max_pseudo), 'synchrony_real'=paste("n segments =", x$nsegment))))+
      theme_light()+
      scale_fill_manual(values = c("brown3","chartreuse4"))+
      theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(text = element_text(vjust = 0, size = 15, family="serif"),  strip.text = element_text(size=15))+
      labs(x="synchrony", title = "histogram of multivariate synchrony", subtitle = paste("segmentsize =", x$segment_size_s, "method = ", x$method, sep = " "))
  } else if (type=="segment-wise") {
    real = x$synchrony[x$synchrony$variable == "synchrony_real", , drop=FALSE]
    real$segment_id = seq_len(nrow(real))
    p = ggplot(real, aes(x=segment_id, y=value))+
      geom_bar(stat = "identity", fill = "chartreuse4", alpha=.8, width=.9)+
      geom_hline(aes(yintercept = x$pseudo_mean,
                     linetype = "mean surrogate synchrony"), color = "brown3", linewidth=1.5, alpha=.5)+
      geom_hline(aes(yintercept = x$pseudo_mean-x$pseudo_sd,
                     linetype = "sd surrogate synchrony"), color = "brown3", linewidth=1, alpha=.5)+
      geom_hline(aes(yintercept = x$pseudo_mean+x$pseudo_sd,
                     linetype = "sd surrogate synchrony"), color = "brown3", linewidth=1, alpha=.5, show.legend=FALSE)+
      labs(title="synchrony per segment", x = "segment", y = "synchrony")+
      theme_minimal()+
      theme(text = element_text(vjust = 0, size = 12, family="serif"),  strip.text = element_text(size=12))+
      theme(legend.position = "bottom", legend.title = element_blank())+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  } else if (type=="time series") {
    data = x$data
    data$counter = seq_len(nrow(data))
    counter = column = NULL ## check NSE notes
    data_ts = as.data.frame(melt(as.data.table(data), id.vars="counter", measure.vars=colnames(x$data), variable.name = "column")[order(counter, column)])
    rm(counter, column)
    sci_palette = pal_npg("nrc", alpha = 0.7)(9)
    sci_palette = colorRampPalette(sci_palette)(ncol(x$data))
    update_geom_defaults("line", list(color = sci_palette))
    p = ggplot(data_ts, aes(x=counter, y=value, group=column))+
      geom_line(linewidth=0.8, alpha=.5)+
      geom_vline(xintercept = seq(0, nrow(data), by = x$segmentHz), alpha=0.25, linetype="longdash")+
      labs(x="time", y="value", colour=NULL)+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(text = element_text(vjust = 0, size = 15, family="serif"))+
      theme(legend.position = "none")
    if (missing(plotly)) plotly = TRUE
    if (plotly && requireNamespace("plotly", quietly=TRUE)) {
      p = plotly::layout(
        plotly::ggplotly(p),
        title=paste("time series with segment size =", x$segment_size_s)
      )
    } else {
      if (plotly) { ## verbose message disabled when plotly=FALSE
        cat("For interactive mvSUSY time series plot install 'plotly' package\n")
      }
      p = p + ggtitle(paste("time series with segment size =", x$segment_size_s))
    }
  }
  p
}
