# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright Wolfgang Tschacher

# Multivariate Surrogate Synchrony (mv-SUSY)
# R-Script by Deborah Meier, Wolfgang Tschacher


# Packages ------------------------------------------------------------------------------------------------------------------

  library(data.table) # Data Transformation - subset, select, and arrange data
  library(dplyr)      # Data Transformation - select, arrange, and mutate data
  library(tidyverse)  # Data Transformation - Pivoting into long and wide data forms

  library(psych)      # Basic descriptive statistics
  library(rstatix)    # Basic statistical tests - t-test, Wilcoxon test  
  library(RcppAlgos)  # Combination and Permutation - Surrogate Analysis

  library(plotly)     # Data Visualization - Interactive Plots
  library(ggplot2)    # Data Visualization - Plots
  library(ggsci)      # Data Visualization - Color palettes inspired by scientific journal

  library(imputeTS)   # Data Visualization and Imputation of missing values - Plots
  library(VIM)        # Data Visualization and Imputation of missing values - Plots

#  1.0 Data import        [import required] ---------------------------------------------------------------------------------------------------------------
### From csv, xlsx, ... [optional]----

 data_table <-read.table(file.choose(),
                        header = TRUE, na.strings=".", 
                        sep = " ", dec = ".")
 data <- data_table
 
## as.numeric(data_table)

 
### Simulation [optional]----

 set.seed(16)
 data_table <- replicate(n = 5, expr = sample(1:10, 500, replace =T)) %>% data.frame()
 data_table$counter <- (1:nrow(data_table))
 colnames(data_table) <- c("t1", "t2", "t3", "t4","t5", "time")
 
 data <- data_table[,1:5]
 
 data_ts1 <- data_table %>% pivot_longer(data_table, cols = 1:5, names_to = "participant" )
 
### Environment [optional]----

 data_table <- data[,1:5]
 
 data_table <- df_Zoom.SPR_end[2:11]
 
 data_table <- motion_12[motion_12$piece == "dean", 1:45] %>% data.frame()
 
 data_table <- motion_12_smooth[motion_12_smooth$piece == "dean",4:48] %>% data.frame()
 
 data_table <- motion_12_noise[motion_12_noise$piece == "dean",1:45] %>% data.frame()
 
 data <- select(motion_12[motion_12$piece == "brahms",1:45], -c(token_zvmuh)) %>% data.frame()
 
 data_table <- df_Kinect_ER10US01R
 
 data_table <- df_Kinect_AL04EL12L
 
 data <- data %>% 
   #filter(piece == "brahms") %>%
   #select(-'mtrkt') %>%
   select(., -c(piece, timestamp)) %>% data.table()
 
 data <- data_10 %>% 
   filter(piece == "beethoven") %>%
   #select(-'mtrkt') %>%
   select(., -c(piece, timestamp)) %>% data.table()
 
 data_table <- data
 
 data <- data[1:100,1:5]
 
 results_global <- results
 
## 1.1 Data overview  [optional] ----
 
 str(data) 
 
 psych::describe(data_table)
 
 Hmisc::describe(data_table)
 
 data_table %>%  VIM::aggr(numbers = TRUE,  axes = TRUE, prop = FALSE)
 
## 1.2 Data processing 

 data <- na_interpolation(data_table[], option = "spline") %>% data.frame()   # use linear, spline or sline
 
 data <- apply(data_table, 2, FUN = function(x){forecast::na.interp(x) }) %>% data.frame()
 
 data <- na.omit(data_table) %>% data.frame()
 
 data <- na.omit(data_table)
 
 data <- as_tibble(data)
 
### 1.2.1 Plot interpolation for a specific column
 
 data_table[,4] %>%
     ggplot_na_imputations(.,na_interpolation(data_table[,4], option = "spline"))+
     theme_minimal()
 

#  1.0 Data import: Datafile must be stored as 'data' in the environment [import required]
 
 data_table <-read.table(file = "/Users/wolfgangtschacher/Documents/Projekte/Konzertprojekt-sneak preview/MEA-Musiker-Auswertungen/SUSY-Berechnungen/MEA-topstage-C3-Beethoven.txt",
                         header = T, na.strings=".", sep = " ", comment.char="", colClasses =  "numeric")
 data <- data_table
 
 #  2.0 Parameter settings [specification required]
 
 name_data      <- "J001_02_RESP_10Hz"  
 seed_number    <-  001            # choose random number: allows to replicate results (surrogates)
 piece          <- ""        # [optional]
 participant    <- ""              # [optional]
 
 fps            <- 10
 segment_size   <- 10
 max_pseudo     <- 1000 
 
 method         <- "lambda_max"       # define method: "lambda_max", "omega", "main_eigenvalue", "range_eigenvalue", "skew_eigenvalue",  "covmax"
 transposition  <- "no"               # yes/ no - only available for 'lambda_max'
 
 delete_memory  <- "yes"              # no - results are attached to previous results/ yes - deletion of previous results
 
## computation
 
 segmentsize    <- segment_size*fps                
 nsegment       <- floor(nrow(data)/segmentsize)
 
 RcppAlgos::permuteCount(nsegment, ncol(data))   # max pseudo segment-combinations
 
## 2.1 Warnings - check parameter settings
 
 if(!method %in% c("lambda_max", "main_eigenvalue", "range_eigenvalue", "skew_eigenvalue", "omega", "covmax") ) 
    warning('method is invalid! choose one of the following computation methods: "lambda_max", "main_eigenvalue", "range_eigenvalue", "skew_eigenvalue", "omega", "covmax" ')
 
 if(segmentsize > (nsegment*segmentsize)/ncol(data))  
    warning(paste('segment size is invalid: Maximum segment size is', ((nsegment)*segmentsize)/ (ncol(data)*fps)))
 
#  3.0 Real synchrony  --------------------------------------------------------------------------------------- 
## 3.1 Reshape data for segmenting
 
 
 if(!is.integer(nrow(data)/segmentsize)) { data_cut <- data[1:(nsegment*segmentsize),]}
 
 data_base <- data_cut
 data_base$segment <- base::rep(paste("segment", 1:c(nsegment), sep = ""), each = segmentsize)
 
 data_base <-  data_base %>%
          group_by(segment) %>%
          mutate(row = row_number()) %>%
          pivot_wider(names_from = segment, values_from = data_base %>% select(-segment) %>% colnames(),
                      names_sep = ".") %>%  
          select(-row) 
 
## 3.2 Generate list of real segment-combinations
 
 comb_real <- paste(colnames(data),sep = ".", "segment",rep(1:(nsegment), each=ncol(data)))
 
 comb_real <- sub(pattern = "segment.", replacement= "segment", comb_real)
 
 comb_real <- split(comb_real, ceiling(seq_along(comb_real)/ncol(data)))
 
 comb_real <- lapply(comb_real, function(i){list(i)})
 
## 3.3 Expand data based on list of real segment-combinations
 
 df_real <- lapply(comb_real[], function(i){ data_base %>% select(c(i[[1]]))})
 
 names(df_real) <- paste0('segment_combo', seq_along(df_real))
 
## 3.4 Define mv-synchrony function (based on method)

   if(method == "lambda_max" & transposition == "no") {
      method_fun <- function(x) {
       eigenwerte <- x %>% cor() %>% replace(is.na(.),0) %>% eigen()
       max(eigenwerte$values)/ (sum(eigenwerte$values)*.01) }
   } else if(method == "lambda_max" & transposition == "yes"){
     method_fun <- function(x) {
       eigenwerte <- x %>% t() %>% cor() %>% replace(is.na(.),0) %>% eigen()
       sum(eigenwerte$values[eigenwerte$values >=1] / (sum(eigenwerte$values)*.01)) }
   } else if(method == "main_eigenvalue"){
      method_fun <- function(x) {
       eigenwerte <- x %>% cor() %>% replace(is.na(.),0) %>% eigen()
       sum(eigenwerte$values[eigenwerte$values >=1] / (sum(eigenwerte$values)*.01)) }
   } else if(method == "range_eigenvalue"){
      method_fun <- function(x) {
       eigenwerte <- x %>% cor() %>% replace(is.na(.),0) %>% eigen()
       max(eigenwerte$values)-min(eigenwerte$values) }
   } else if(method == "skew_eigenvalue"){
      method_fun <- function(x) {
        eigenwerte <- x %>% cor() %>% replace(is.na(.),0) %>% eigen()
        skew(eigenwerte$values) }
   } else if(method == "omega"){
       method_fun <- function(x) {
          Hact <- x %>% cov() %>% abs() %>% det() 
          Hpot <- x %>% cov() %>% diag() %>% prod()
          Omega <- 1-(Hact/Hpot)
          Omega}
   } else if(method == "covmax"){
      method_fun <- function(x) {
        cov_act <- x %>% cov() %>% abs() %>% .[lower.tri(., diag = F)] %>% mean()
        cov_pot <- x %>% cov() %>% diag() %>% mean()
        covmax <- cov_act/cov_pot
        covmax}
    }
 
## 3.5 Apply synchrony function per real segment-combination   

  synchrony_real <- sapply(df_real[], function(i){method_fun(i)})
  
## 3.6 Define eigenvalue function (based on method) 
  
  if(method %in% c("lambda_max", "main_eigenvalue", "range_eigenvalue", "skew_eigenvalue", "covmax")) {
     fun_eigenvalue <- function(x) {
        eigenwerte <- x %>% cor() %>% replace(is.na(.),0) %>% eigen()
        eigenwerte$values }
  }
 
#  4.0 Pseudo synchrony  --------------------------------------------------------------------------------------
## 4.1 Generate list of pseudo segment-combinations

 comb_pseudo <- c(paste("segment", sep ="", rep(1:(nsegment))))
  
 comb_pseudo <- RcppAlgos::permuteSample(v = comb_pseudo, 
                                        m = ncol(data), n = max_pseudo, seed = seed_number,
                                        FUN = function(i){paste(colnames(data),i,sep = ".")})

## 4.2 Expand data based on list of pseudo segment-combinations
 
 df_pseudo <- lapply(comb_pseudo, function(i){ data_base %>% select(all_of(i))})

 names(df_pseudo) <- paste0('segment_combo', seq_along(df_pseudo))
 
## 4.3 Apply synchrony function per pseudo segment-combination

 synchrony_pseudo <- sapply(df_pseudo[], function(i){method_fun(i)})
 
#  5.0 Summarize synchrony data --------------------------------------------------------------------------------------
## 5.1 Reshape synchrony data
 
 matrix_pseudo   <- synchrony_pseudo %>% 
   as.list %>% 
   data.frame() %>%
   pivot_longer(., cols = colnames(.), names_to = "variable", values_to = "value") %>%
   mutate(variable = "synchrony_pseudo")
 
 matrix_real <- synchrony_real %>%
   as.list %>% 
   data.frame() %>%
   pivot_longer(., cols = colnames(.), names_to = "variable", values_to = "value") %>%
   mutate(variable = "synchrony_real")
 
 matrix_synchrony <- rbind(matrix_pseudo, matrix_real)
 
 matrix_synchrony$variable <- matrix_synchrony$variable %>% as.factor()
 
### 5.1.1 Parametric and nonparametric tests
 
 t_tests <- matrix_synchrony %>% t.test(value ~ variable, data = .)
 
 wilcox_tests <- matrix_synchrony %>% wilcox.test(value ~ variable,data= .)
 
### 5.1.2 Real and pseudo synchrony indices
 
 mean_synchrony <- mean(matrix_real$value)
 
 sd_synchrony <- sd(matrix_real$value)
 
 mean_synchrony_pseudo <- mean(matrix_pseudo$value)
 
 sd_synchrony_pseudo <- sd(matrix_pseudo$value)
 
### 5.1.3 Effect size
 
 ES_synchrony <- (mean_synchrony - mean_synchrony_pseudo)/ sd_synchrony_pseudo
 
## 5.2 Generate eigenvalue data
 
 matrix_eigenvalue <- sapply(df_real[], function(i){fun_eigenvalue(i)})
 
 matrix_eigenvalue <- matrix_eigenvalue %>%
    data.frame() %>%
    pivot_longer(., cols = colnames(matrix_eigenvalue), names_to = "segment", values_to = "eigenvalue") %>%
    mutate(data = "real")
 
 matrix_eigenvalue_pseudo <- sapply(df_pseudo[], function(i){fun_eigenvalue(i)})
 
 colnames(matrix_eigenvalue_pseudo) <- paste0('segment_combo_pseudo', seq_along(colnames(matrix_eigenvalue_pseudo)))
 
 matrix_eigenvalue_pseudo <- matrix_eigenvalue_pseudo %>%
    data.frame() %>%
    pivot_longer(., cols = colnames(matrix_eigenvalue_pseudo), names_to = "segment", values_to = "eigenvalue") %>%
    mutate(data = "pseudo")
 
 matrix_eigen <- rbind(matrix_eigenvalue_pseudo, matrix_eigenvalue) %>% 
    group_by(data, segment) %>%
    mutate(segment_num = order(eigenvalue, decreasing = T)) %>%
    ungroup()
 
#  6.0 Output ------------------------------------------------------------------------------------------------------------------------------
## 6.1 Print output
 
 mv_SUSY <- cbind(name_data, method, piece, participant, ncol(data), nrow(data_cut), seed_number, length(comb_pseudo) , segment_size, segmentsize*ncol(data),
                  round(mean_synchrony,5), round(sd_synchrony,5),
                  round(mean_synchrony_pseudo,5), round(sd_synchrony_pseudo,5), round(ES_synchrony,5), 
                  round(t_tests$statistic,5),
                  format(t_tests$p.value,scientific = F), format(wilcox_tests$statistic,scientific = F), format(wilcox_tests$p.value,scientific = F))
 
 colnames(mv_SUSY) <- c("data", "method", "piece", "token", "n(col)", "n(row)", "seed", "n(pseudo)", "segment size (s)", "data per segment",
                        "mean(synchrony real)", "sd(synchrony real)", 
                        "mean(synchrony pseudo)", "sd(synchrony pseudo)", "ES", "t statistic", "p-value", "statistic nonpar", "p-value nonpar") 
 
 if(delete_memory == "yes") { results <- data.frame(matrix(vector(), ncol = 19))}
 if(delete_memory == "yes") { colnames(results) <- colnames(mv_SUSY)}
 
 View(results<- rbind(results[], mv_SUSY))


 
 
# Save output as .txt file [optional]----
 
 write.table(results, 
             "mv-SUSY_sneak-preview_motion_piece.txt", 
             #paste(name_data,"_",fps,"fps_",method,"_",segment_size,"s",".txt",sep = ""),
             append = FALSE, sep = " ", dec = ".",
             row.names = F, col.names = TRUE)
 
 write.table(results[132:134,],                                            # define rows from results to be stored in the txt-file
             "mv-SUSY_sneak-preview_motion_concert.txt", 
             #paste(name_data,"_",fps,"fps_",method,"_",segment_size,"s",".txt",sep = ""),
             append = FALSE, sep = " ", dec = ".",
             row.names = F, col.names = TRUE)
 
#  7.0 Plots              [optional]---- -----------------------------------------------------------------------------------------------------------------------------------
## 7.1 Eigenvalue (screeplot per segment)
 
 p1<- matrix_eigen %>%
    ggplot(., aes(x = as.factor(segment_num), group = segment)) +
    geom_line(aes(y = eigenvalue, colour = data)) +
    scale_color_manual(values = c('real' = 'chartreuse4', 'pseudo' = 'brown3'))+   
    scale_y_continuous(breaks = seq(0, max(matrix_eigen$eigenvalue), by = 1)) +
    #geom_line(aes(y = 1))+
    facet_wrap(~data, labeller = as_labeller(c(`pseudo` = paste( "n segments =", max_pseudo), `real` = paste( "n segments =", nsegment)))) +
    labs(y = "eigenvalue", x = "dimension", 
         title = "eigenvalues per segment for real and pseudo data", subtitle = paste( "segmentsize =", segment_size)) +
    theme_light()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_blank(), text = element_text(vjust = 0, size = 15, family="serif"), strip.text = element_text(size=15))
 p1
 
## 7.2 density histogram
 
 p2 <- matrix_synchrony %>%
    ggplot(., aes(x = value, colour = variable))+
    geom_density(aes(x = value, fill = variable), alpha = .1, size = 1, show.legend = F) +
    geom_rug(aes(x = value, y = 0)) +
    #xlim(5.5,13) +
    scale_colour_manual('data', 
                        label = c('synchrony_pseudo'= "pseudo", 'synchrony_real'="real"),
                        values = c('synchrony_real' = 'chartreuse4', 'synchrony_pseudo' = 'brown3'))+
    scale_fill_manual('data', 
                      values = c('synchrony_real' = 'chartreuse4', 'synchrony_pseudo' = 'brown3'))+
    labs(y = "density", x = "synchrony") +
    theme_light()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), text = element_text(vjust = 0, size = 12, family="serif")) +
    #labs(title = "density histogram of multivariate synchrony", subtitle = paste( "segmentsize =", segment_size,"method = ", method, sep = " ") +
    labs(title = "Dean")      
 p2
 
 sjPlot::save_plot("sneak-preview_11_dean.jpg", width = 13, height = 8)
 
## 7.3 Free scale histogram 
 
 vlines <- matrix_synchrony %>% 
    group_by(variable) %>%
    summarise(mean.var = mean(value), .groups = 'drop')
 
 p3 <- matrix_synchrony %>%
    ggplot(., aes(x = value, fill = variable))+
    geom_histogram(color="#e9ecef",alpha=0.8, position = 'identity', bins = nsegment)+
    #stat_bin(bins = nsegment, geom = "text", aes(label = ..count.., group = variable)) +
    geom_vline(data = vlines, aes(xintercept = mean.var),linetype = "dashed", size = 0.8, alpha = 0.9)+
    facet_wrap(~variable, scales ="free", labeller = as_labeller(c('synchrony_pseudo'= paste("n segments =", max_pseudo), 'synchrony_real'=paste("n segments =", nsegment))))+
    theme_light()+
    scale_fill_manual(values = c("brown3","chartreuse4"))+
    theme(legend.position = "none",panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(text = element_text(vjust = 0, size = 15, family="serif"),  strip.text = element_text(size=15))+
    labs(x="synchrony", title = "histogram of multivariate synchrony", subtitle = paste( "segmentsize =", segment_size,"method = ", method, sep = " ") )
 p3
 
 p3 <- matrix_synchrony %>%
   ggplot(., aes(x = value, fill = variable))+
   geom_histogram(color="#e9ecef",alpha=0.8, position = 'identity', bins = nsegment*.5)+
   #stat_bin(bins = nsegment, geom = "text", aes(label = ..count.., group = variable)) +
   geom_vline(data = vlines, aes(xintercept = mean.var),linetype = "dashed", size = 0.8, alpha = 0.9)+
   facet_wrap(~variable, scales ="free",labeller = as_labeller(c('synchrony_pseudo'= paste("n segments =", max_pseudo), 'synchrony_real'=paste("n segments =", nsegment))))+
   theme_light()+
   scale_fill_manual(values = c("brown3","chartreuse4"))+
   theme(legend.position = "none",panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
   theme(text = element_text(vjust = 0, size = 15, family="serif"),  strip.text = element_text(size=15))+
   labs(x="synchrony", title = "histogram of multivariate synchrony", subtitle = paste( "segmentsize =", segment_size,"method = ", method, sep = " ") )
 p3
 
## 7.4 Segment-wise synchrony
 
 p4 <- matrix_real %>%
    ggplot(., aes(x=rownames(.), y=synchrony_real))+
    geom_bar(stat = "identity", fill = "chartreuse4", alpha=.8, width = .8)+
    coord_flip() +
    scale_x_discrete(limits = as.numeric(rownames(matrix_real))) +
    geom_hline(aes(yintercept = mean_synchrony_pseudo, 
                   linetype = "mean synchrony in pseudo segments"), color = "brown3", size =1.5, alpha=.5)+
    geom_hline(aes(yintercept = mean_synchrony_pseudo-sd_synchrony_pseudo, 
                   linetype = "sd synchrony in pseudo segments"), color = "brown3", size =1, alpha=.5)+
    geom_hline(aes(yintercept = mean_synchrony_pseudo+sd_synchrony_pseudo, 
                   linetype = "sd synchrony in pseudo segments"), color = "brown3", size =1, alpha=.5)+
    scale_x_discrete(labels = c(1:c(nsegment)))+
    #labs(title=paste( "multivariate synchrony per segment"), subtitle = paste( "segmentsize =", segment_size, "method = ", method, sep = " "), x="segment", y =      "measure of synchrony") +
    labs(title=paste( "synchrony per segment"), x="segment", y = "synchrony") + 
    theme(text = element_text(vjust = 0, size = 12, family="serif"),  strip.text = element_text(size=12))+
    theme(legend.position = "bottom")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
 p4
 
 p4 <- matrix_real %>%
   mutate(segment_nr = rownames(.) %>% as.numeric()) %>%
   ggplot(., aes(x=segment_nr, y=value))+
   geom_bar(stat = "identity", fill = "chartreuse4", alpha=.8, width = .9) +
   #coord_flip() +
   #xlim(nrow(matrix_real)+1,0) +
   geom_hline(aes(yintercept = mean_synchrony_pseudo, 
                  linetype = "mean surrogate synchrony"), color = "brown3", size =1.5, alpha=.5)+
   geom_hline(aes(yintercept = mean_synchrony_pseudo-sd_synchrony_pseudo, 
                  linetype = "sd surrogate synchrony"), color = "brown3", size =1, alpha=.5)+
   geom_hline(aes(yintercept = mean_synchrony_pseudo+sd_synchrony_pseudo, 
                  linetype = "sd surrogate synchrony"), color = "brown3", size =1, alpha=.5, show.legend = F)+
   #scale_x_discrete(labels = c(1:c(nsegment)))+
   #labs(title=paste( "multivariate synchrony per segment"), subtitle = paste( "segmentsize =", segment_size, "method = ", method, sep = " "), x="segment", y =      "measure of synchrony") +
   labs(title=paste( "synchrony per segment"), x="segment", y = "synchrony") + 
   theme_minimal()+
   theme(text = element_text(vjust = 0, size = 12, family="serif"),  strip.text = element_text(size=12))+
   theme(legend.position = "bottom", legend.title = element_blank())+
   theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
 p4
 
## 7.5 Time series
 
 data_ts <- data %>%
   mutate(counter = seq(nrow(data)))
   
  data_ts <- data_ts %>%  
   pivot_longer(., cols = data_ts %>% select(-counter) %>% colnames(), names_to = "column", values_to = "value")
 
 sci_palette <- pal_npg("nrc", alpha = 0.7)(9)
 sci_palette <- colorRampPalette(sci_palette)(ncol(data))
 update_geom_defaults("line", list(color = sci_palette))
 
 plot_ts <- data_ts %>% 
    ggplot(., aes(x=counter ,y=value, group=column)) +
    geom_line(size=0.8, alpha = .5)+
    geom_vline(xintercept = seq(0, nrow(data), by = segmentsize), alpha=0.25, linetype = "longdash")+
    labs(x="time", y="value", colour=NULL) +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    theme(text = element_text(vjust = 0, size = 15, family="serif"))+
    theme(legend.position = "none")

 plot_ts <- plot_ts %>%
    ggplotly() %>%
    layout(title=paste("time series with segment size =", segment_size))
 plot_ts
 
