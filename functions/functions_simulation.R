########################
# function simulaition #
########################

add_noise <- function(x, per){
  
  n_days <- dim(x)[2]
  
  nval_sel <- round(dim(x)[1] * per)
  
  out_list <- vector(mode = 'list', length = dim(x)[2])
  
  for(i in 1:n_days){out_list[[i]] <- boxplot(x[, i])$out}
  
  # split values between above and below the mean
  
  av_tr <- colMeans(x, na.rm = TRUE)
  
  out_positions <- vector(mode = 'list', length = n_days)
  
  for(i in 1:n_days){
    
    v_low <- out_list[[i]][out_list[[i]] < av_tr[i]]
    v_high <- out_list[[i]][out_list[[i]] > av_tr[i]]
    v_tot <- c(v_low, v_high)
    
    per_low <- length(v_low)/length(v_tot)
    per_high <- length(v_high)/length(v_tot)
    
    nval_sel_low <- round(nval_sel * per_low) - length(v_low)
    nval_sel_high <- round(nval_sel * per_high) - length(v_high)
    
    # sample values
    
    samp_high <- sample(x = v_high, size = nval_sel_high, replace = TRUE)
    samp_high <- samp_high + abs(rnorm(n = nval_sel_high))
    
    samp_low <- sample(x = v_low, size = nval_sel_low, replace = TRUE)
    samp_low <- samp_low - abs(rnorm(n = nval_sel_low))
    
    extra_out <- c(samp_low, samp_high)
    
    # replace the values
    
    pos_out <- which(x[, i] %in% out_list[[i]])
    samp_pos <- 1:dim(x)[1]
    samp_pos <- samp_pos[-pos_out]
    pos_rep <- sample(x = samp_pos, size = length(extra_out), replace = FALSE)
    
    x[pos_rep, i] <- extra_out
    
    # store outliers positions
    
    out_positions[[i]] <- c(pos_rep, pos_out)
    
    
  }
  
  return(list(pheno_data = x, out_pos = out_positions))
  
}

add_noise2 <- function(x, per, coeff_std = 3){
  
  n_days <- dim(x)[2]
  N <- dim(x)[1]
  
  nval_sel <- round(N * per)
  
  out_list <- vector(mode = 'list', length = dim(x)[2])
  
  for(i in 1:n_days){out_list[[i]] <- boxplot(x[, i], plot = FALSE)$out}
  
  # split values between above and below the mean
  
  av_tr <- colMeans(x, na.rm = TRUE)
  
  out_positions <- vector(mode = 'list', length = n_days)
  
  for(i in 1:n_days){
    
    v_low <- out_list[[i]][out_list[[i]] < av_tr[i]]
    v_high <- out_list[[i]][out_list[[i]] > av_tr[i]]
    v_tot <- c(v_low, v_high)
    
    per_low <- length(v_low)/length(v_tot)
    per_high <- length(v_high)/length(v_tot)
    
    nval_sel_low <- round(nval_sel * per_low) 
    nval_sel_high <- round(nval_sel * per_high)
    
    std_i <- sd(x[, i], na.rm = TRUE)
    
    noise_low <- -abs(rnorm(n = nval_sel_low, mean = 0, sd = coeff_std * std_i))
    noise_high <- abs(rnorm(n = nval_sel_high, mean = 0, sd = coeff_std * std_i))
    
    noise_vect <- c(noise_low, noise_high)
    
    av_pos <- 1:N
    av_pos <- av_pos[!is.na(x[, i])]
    
    pos_samp <- sample(x = av_pos, size = nval_sel_low + nval_sel_high)
    
    n_values <- x[pos_samp, i] + noise_vect
    n_values[n_values <0] <- 0
    
    x[pos_samp, i] <- n_values
    
    # store outliers positions
    
    out_positions[[i]] <- pos_samp
    
  }
  
  return(list(pheno_data = x, out_pos = out_positions))
  
}

add_miss <- function(x, per, out_pos = NULL){
  
  n_days <- dim(x)[2]
  nval_sel <- round(dim(x)[1] * per)
  
  for(i in 1:n_days){
    
    if(is.null(out_pos)){
      
      sel_pos <- sample(x = 1:dim(x)[1], size = nval_sel, replace = FALSE)
      x[sel_pos, i] <- NA
      
    } else {
      
      av_pos <- 1:dim(x)[1]
      av_pos <- av_pos[-out_pos[[i]]]
      sel_pos <- sample(x = av_pos, size = nval_sel, replace = FALSE)
      
      x[sel_pos, i] <- NA
      
    }
    
  }
  
  return(x)
  
}


SimVal_proc <- function(exp_des_data, pheno_data, out_det = TRUE,
                        miss_imp = TRUE, fixed = NULL,
                        random = ~ rep +  rep:block + row_f + col_f, ref_day,
                        single_mixed_model = FALSE, out_p_val = 0.05,
                        G_BLUEs_ref) {
  
  # Check if specified column were included in the experimental design
  # with the right format.
  ################
  
  if(!('genotype' %in% colnames(exp_des_data))){
    
    stop('There is no column labelled genotype in exp_des_data')
  }
  
  if(!('col' %in% colnames(exp_des_data))){
    
    stop('There is no column labelled col in exp_des_data')
  }
  
  if(!is.numeric(exp_des_data$col)){
    
    stop('The col information in exp_des_data must be numeric')
  }
  
  if(!('row' %in% colnames(exp_des_data))){
    
    stop('There is no column labelled row in exp_des_data')
  }
  
  if(!is.numeric(exp_des_data$row)){
    
    stop('The row information in exp_des_data must be numeric')
  }
  
  if(!('col_f' %in% colnames(exp_des_data))){
    
    stop('There is no column labelled col_f in exp_des_data')
  }
  
  if(!is.factor(exp_des_data$col_f)){
    
    stop('The col_f information in exp_des_data must be factor')
  }
  
  if(!('row_f' %in% colnames(exp_des_data))){
    
    stop('There is no column labelled row_f in exp_des_data')
  }
  
  if(!is.factor(exp_des_data$row_f)){
    
    stop('The row_f information in exp_des_data must be factor')
  }
  
  ##############
  
  # transform variable into factor
  
  
  
  geno_id <- unique(exp_des_data$genotype)
  geno_id <- as.character(geno_id)
  n_geno <- length(geno_id)
  
  G_BLUEs_ref <- G_BLUEs_ref[geno_id]
  
  if(!is.null(fixed)){fixed <- as.formula(fixed)}
  
  exp_des_data$genotype <- factor(exp_des_data$genotype)
  
  if(!single_mixed_model){ # regular procedure
    
    # Outliers detection
    if(out_det) {
      
      pheno_data <- outliers_det_boxplot(data = pheno_data, plot = FALSE)
      
    }
    
    # Missing values imputation
    if(miss_imp){
      
      pheno_data <- miss_imp_PMM(data = pheno_data, plot = FALSE)
      
    }
    
    # SpATS model with spatial adjustment
    
    data <- cbind.data.frame(exp_des_data, pheno_data[, ref_day])
    colnames(data)[dim(data)[2]] <- 'pheno'
    
    m <- tryCatch(SpATS(response = "pheno", genotype = "genotype",
                        geno.decomp = NULL, genotype.as.random = FALSE,
                        spatial = ~PSANOVA(col, row, nseg = c(20,20)),
                        fixed = fixed,
                        random = as.formula(random),
                        data = data,
                        control = list(maxit = 50, tolerance = 1e-06, monitoring = 1)),
                  error = function(e) NULL)
    
    if(!is.null(m)){
      
      pred <- predict(m, which = 'genotype')
      BLUE_i <- pred$predicted.values
      names(BLUE_i) <- as.character(pred$genotype)
      
      G_BLUEs_ref_i <- G_BLUEs_ref[as.character(pred$genotype)]
      
      print(plot(G_BLUEs_ref_i, BLUE_i))
      res_i <- rmse(actual = BLUE_i, predicted = G_BLUEs_ref_i)^2
      cor_i <- cor(BLUE_i, G_BLUEs_ref_i, use = 'complete.obs')
      
      return(list(RMSE= res_i, cor = cor_i))
      
    } else {
      
      
      return(list(RMSE= NA, cor = NA))
      
    }
    
  } else { # single step mixed model
    
    data <- cbind.data.frame(exp_des_data, pheno_data[, ref_day])
    colnames(data)[dim(data)[2]] <- 'pheno'
    
    p_val <- 0
    
    while(p_val < out_p_val){
      
      m <- tryCatch(SpATS(response = "pheno", genotype = "genotype",
                          geno.decomp = NULL, genotype.as.random = FALSE,
                          spatial = ~PSANOVA(col, row, nseg = c(20,20)),
                          fixed = fixed,
                          random = as.formula(random),
                          data = data,
                          control = list(maxit = 50, tolerance = 1e-06, monitoring = 1)),
                    error = function(e) NULL)
      
      # test null
      
      if(!is.null(m)){
        
        m_resid <- m$residuals
        
        test <- grubbs.test(m_resid, type = 10, two.sided = TRUE) # Grubb test
        p_val <- test$p.value
        
        if(p_val < out_p_val){ # test if p_val is lower than ... (presence of an outlier)
          
          # remove the outlying value from the data
          
          pos_max_res <- which(abs(m_resid) == max(abs(m_resid), na.rm = TRUE))
          
          data$pheno[pos_max_res] <- NA
          
          
        }
        
      } else{
        
        break()
        
      }
      
    }
    
    
    if(!is.null(m)){
      
      pred <- predict(m, which = 'genotype')
      BLUE_i <- pred$predicted.values
      names(BLUE_i) <- as.character(pred$genotype)
      
      G_BLUEs_ref_i <- G_BLUEs_ref[as.character(pred$genotype)]
      
      print(plot(G_BLUEs_ref_i, BLUE_i))
      res_i <- rmse(actual = BLUE_i, predicted = G_BLUEs_ref_i)^2
      cor_i <- cor(BLUE_i, G_BLUEs_ref_i, use = 'complete.obs')
      
      return(list(RMSE= res_i, cor = cor_i))
      
    } else {
      
      return(list(RMSE= NA, cor = NA))
      
    }
    
  }
  
  
  
}