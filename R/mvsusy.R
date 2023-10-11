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
    t_tests = t_tests, wilcox_tests = wilcox_tests
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
plot.mvsusy = function(x,  ...) {
  if (!inherits(x, "mvsusy"))
    stop("'x' must be an object of class 'mvsusy'")
  invisible()
}
