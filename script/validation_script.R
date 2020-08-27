################################################################################
################################################################################
#
# SpaTemHTP: A Data Analysis Pipeline for Efficient Processing and Utilization
# of Temporal High-Throughput Phenotyping Data
#
# Soumyashree Kar, Vincent Garin, Jana Kholova, Vincent Vadez, Surya S. Durbha,
# Ryokei Tanaka, Hiroyoshi Iwata, Milan O. Urban, J Adinarayana
#
################################################################################
################################################################################

setwd('F:/ICRISAT/Phenotyping/LeasyScan_article/SpaTemHTP_Validation')
# set your working directory here (change later)

# Library

# source('./functions/detect_OL_Boxplot.R') ### package functions
# source('./functions/GapFill_mice.R')
# source('./functions/det_OL_Box_vect.R')
# source('./functions/geno_comp_miss.R')

# download the SpaTemHTP package from the github repository
# https://github.com/ICRISAT-GEMS/SpaTemHTP


# devtools::install_github("ICRISAT-GEMS/SpaTemHTP")

library(SpaTemHTP)


# library(SpaTempHTP) change later when package is ready
library(SpATS)
library(lme4)
library(VIM)
library(mice)
library(plyr)
library(outliers)
library(R.utils)
library(dplyr)
library(drc)


# Cross-validation
##################

# Datasets: 2 populations (chickpea-CP, sorghum-SG), 2 experiments (E1, E2), 2
# traits (leaf area 3D-LA3D, plant height-PH)
data_id <- c("CP_E1_LA3D", "CP_E1_PH", "CP_E2_LA3D", "CP_E2_PH",
            "SG_E1_LA3D", "SG_E1_PH", "SG_E2_LA3D", "SG_E2_PH")

n_days <- c(23, 23, 38, 38, 23, 23, 22, 22)
n_e_d_col <- 5 # Number of experimental design information columns

n_data <- length(data_id)

### Strategy grid

out_detect <- c(FALSE, TRUE)
miss_imp <- c(FALSE, TRUE)
spat_adj <- c(FALSE, TRUE)

strat_tab <- expand.grid(out_detect, miss_imp, spat_adj, stringsAsFactors = FALSE)

colnames(strat_tab) <- c('out_det', 'miss_imp', 'spat_adj')

n_scen <- dim(strat_tab)[1] + 1  ### Check how the scenario 9 is done
n_rep <- 5                      ### Potential change if we select all days and not only 5.
n_fold <- 10

### Space to store the results

pred_res <- vector(mode = 'list', length = n_data) # Precition ability results
he_res <- vector(mode = 'list', length = n_data) # heritability results
names(pred_res) <- names(he_res) <- data_id

for(i in 1:n_data){
  pred_res[[i]] <- matrix(NA, nrow = n_scen, ncol = n_rep*n_fold)
  he_res[[i]] <- matrix(NA, nrow = n_scen, ncol = n_rep*n_fold)
}


### Loop over the different dataset

for(r in 1:n_data){
  
  load(file.path('data', paste0(data_id[r], '.RData')))
  
  # select different date for each CV replication 
  rep_time_id <- sample(1:n_days[r], n_rep)
  # rep_time_id <- 6:n_days[r] ### could be changed
  
  run_ind <- 1 # indicator to fill the table results
  
  for(i in 1:n_rep){  # cross-validation loop
    
    ### Split the observed phenotypic data into n_fold
    
    # select the position of the observed values. We can only predict what is observed
    obs_val_id <- which(!is.na(data[, n_e_d_col + rep_time_id[i]]))
    
    folds <- cut(seq(1, length(obs_val_id)), breaks = n_fold, labels = FALSE)
    folds <- sample(folds)
    
    for(j in 1:n_fold){ # Within a replication iterate over the folds 
      
      data_j <- data # keep data as reference raw data
      
      # split data into training and validation
      vs_id <- obs_val_id[which(folds == j)]
      data_ts <- data_j[-vs_id, ]
      data_vs <- data_j[vs_id, ]
      
      for(z in 1:4){  # iteration over the different strategies (S1-S9)
                      # for each four combinations of out_det and miss_imp
                      # we directly cacluate the spatially and non-spatially
                      # adjusted models.
        
        data_ts_z <- data_ts # Reinitiate the TS and VS dataset at each iteration
        data_vs_z <- data_vs
        
        # Outliers detection
        if(strat_tab[z, ]$out_det) {
          
          data_ts_z <- data.frame(data_ts_z[, 1:5],
                                  outliers_det_boxplot(data = data_ts_z[, 6:dim(data_ts)[2]]))
          
        }
        
        # Missing values imputation
        if(strat_tab[z, ]$miss_imp){
          
          data_ts_z <- data.frame(data_ts_z[, 1:5],
                                  miss_imp_PMM(data = data_ts_z[, 6:dim(data_ts)[2]]))
          
        }
        
        ### Mixed model computation
        
        data_ts_z$col_f <- factor(data_ts_z$col)
        data_ts_z$row_f <- factor(data_ts_z$row)
        data_ts_z$block <- factor(data_ts_z$block)
        data_ts_z$rep <- factor(data_ts_z$rep)
        data_ts_z$genotype <- factor(data_ts_z$genotype)
        
        
        data_vs_z$col_f <- factor(data_vs_z$col)
        data_vs_z$row_f <- factor(data_vs_z$row)
        data_vs_z$block <- factor(data_vs_z$block)
        data_vs_z$rep <- factor(data_vs_z$rep)
        data_vs_z$genotype <- factor(data_vs_z$genotype)
        
        # set the phenotypic value
        
        ph_id <- n_e_d_col + rep_time_id[i]
        colnames(data_ts_z)[ph_id] <- colnames(data_vs_z)[ph_id] <- 'pheno'
        
        
        # check that there are no genotype with complete missing values
        prob_geno <- geno_comp_miss(d = data_ts_z)
        
        if(!is.null(prob_geno)){ 
          
          data_ts_z <- data_ts_z[!(data_ts_z$genotype %in% prob_geno), ]
          
          if(any(data_vs_z$genotype %in% prob_geno)){
            
            MgenoVs_id <- which(data_vs_z$genotype %in% prob_geno)
            data_vs_z <- data_vs_z[-MgenoVs_id, ]
            
          }
          
        }
        
        ### mixed model without spatial adjustment
        
        m_no_Sadj <- tryCatch(lmer(pheno ~ 1 + (1|rep) + (1|rep:block) +
                                     (1|row_f) + (1|col_f) + (1|genotype),
                                   data = data_ts_z), error = function(e) NULL)
        
        if(!is.null(m_no_Sadj)){
          
          pheno_pred <- tryCatch(predict(m_no_Sadj, newdata = data_vs_z,
                                         allow.new.levels = TRUE),
                                 error = function(e) NULL)
          
          var.mat <- as.data.frame(VarCorr(m_no_Sadj))
          Geno_pos <- which(var.mat$grp == 'genotype')
          Res_pos <- which(var.mat$grp == 'Residual')
          he <- var.mat[Geno_pos, 4]/(var.mat[Geno_pos, 4]+var.mat[Res_pos, 4])
          print(paste('h2(LMM) = ', round(he, 2)))
          he_res[data_id[r]][[1]][z, run_ind] <- he
          
          if(!is.null(pheno_pred)){
            
            cor_res <- cor(pheno_pred, data_vs_z$pheno, use = 'complete.obs')
            print(paste('cor(LMM) = ', round(cor_res, 2)))
            pred_res[data_id[r]][[1]][z, run_ind] <- cor_res
            
          }
          
        }
          
          ### mixed model with spatial adjustment (SpATS)
          
          m_Sadj <- tryCatch(SpATS(response = "pheno", genotype = "genotype",
                                   geno.decomp = NULL, genotype.as.random = TRUE,
                                   spatial = ~ PSANOVA(col, row, nseg = c(20,20)),
                                   random = ~ rep +  rep:block + row_f + col_f,
                                   data = data_ts_z,
                                   control = list(maxit = 50, tolerance = 1e-06, monitoring = 1)),
                             error = function(e) NULL)
          
          if(!is.null(m_Sadj)){
            
            pheno_pred <- tryCatch(predict.SpATS(m_Sadj, newdata = data_vs_z),
                                   error = function(e) NULL)
            
            he_res[data_id[r]][[1]][z+4, run_ind] <- he_Sadj <- getHeritability(m_Sadj)
            print(paste('h2(SpATS) = ', round(he_Sadj, 2)))
            
            if(!is.null(pheno_pred)){
              
              cor_res <- cor(pheno_pred$pheno, pheno_pred$predicted.values,
                             use = 'complete.obs')
              print(paste('cor(SpATS) = ', round(cor_res, 2)))
              pred_res[data_id[r]][[1]][z+4, run_ind] <- cor_res
              
            } 
            
          } 
          
          # strategy (S9) single-step mixed model analysis
          
          if (z == 1){ # Start from the 1st strategy with non modified data
            
            print('strategy 9')
            
            p_val <- 0
            count_iter <- 0
            
            while(p_val < 0.05){
              
              # compute the mixed model
              
              m_Sadj <- tryCatch(SpATS(response = "pheno", genotype = "genotype",
                                       geno.decomp = NULL, genotype.as.random = TRUE,
                                       spatial = ~PSANOVA(col, row, nseg = c(20,20)),
                                       random = ~ rep +  rep:block + row_f + col_f,
                                       data = data_ts_z,
                                       control = list(maxit = 50, tolerance = 1e-06,
                                                      monitoring = 1)),
                                 error = function(e) NULL)
              
              # test null
              
              if(!is.null(m_Sadj)){
                
                m_resid <- m_Sadj$residuals
                
                test <- grubbs.test(m_resid, type = 10, two.sided = TRUE) # Grubb test
                p_val <- test$p.value
                
                if(p_val < 0.05){ # test if p_val is lower than 0.05 (presence of an outlier)
                  
                  # remove the outlying value from the data
                  
                  pos_max_res <- which(abs(m_resid) == max(abs(m_resid), na.rm = TRUE))
                  
                  data_ts_z$pheno[pos_max_res] <- NA
                  
                  count_iter <-  count_iter + 1
                  
                  print(paste0('out_iteration_', count_iter))
                  
                }
                
              } else{
                
                break()
                
              }
              
            }
            
            # save the results
            
            if(!is.null(m_Sadj)){
              
              pheno_pred <- tryCatch(predict.SpATS(m_Sadj, newdata = data_vs_z),
                                     error = function(e) NULL)
              
              he_res[data_id[r]][[1]][z+8, run_ind] <- he_Sadj <- getHeritability(m_Sadj)
              print(paste('h2(SpATS) = ', round(he_Sadj, 2)))
              
              if(!is.null(pheno_pred)){
                
                cor_res <- cor(pheno_pred$pheno, pheno_pred$predicted.values,
                               use = 'complete.obs')
                print(paste('cor(SpATS) = ', round(cor_res, 2)))
                pred_res[data_id[r]][[1]][z+8, run_ind] <- cor_res
                
              } 
              
            }
            
          } # end single-step mixed model analysis (S9)
        
      } # End different strategies computation
      
      run_ind <- run_ind + 1
      
    } # End fold loop
    
  } # End CV loop
  
} # End iteration over the different datasets

### results processing

# construction of the average table

# calculation of the t-test

# Final ouput of this section should be the table 3 of the manuscript

###########

# Between experiment comparisons
################################

data_id <- data.frame(E1 = c("CP_E1_LA3D", "CP_E1_PH", "SG_E1_LA3D", "SG_E1_PH"),
                      E2 = c("CP_E2_LA3D", "CP_E2_PH", "SG_E2_LA3D", "SG_E2_PH"))

n_days <- c(23, 23, 38, 38, 23, 23, 22, 22)
n_e_d_col <- 5 # Number of experimental design information columns

n_data <- length(data_id)

out_detect <- c(FALSE, TRUE)
miss_imp <- c(FALSE, TRUE)

strat_tab <- expand.grid(out_detect, miss_imp, stringsAsFactors = FALSE)

colnames(strat_tab) <- c('out_det', 'miss_imp')

# create space to store the results

# Structure of the list:

# a) level 1: four combination of crop and trait (CP_LA3D, CP_PH, SG_LA3D, SG_PH)
# b) level 2: within level 1, two experiments (E1, E2)
# c) level 3: within level 2, four strategies (S1-5)

res <- vector(mode = 'list', length = 4)
names(res) <- c('CP_LA3D', 'CP_PH', 'SG_LA3D', 'SG_PH')

for(i in 1:4){
  
  res[[i]] <- vector(mode = 'list', length = 2)
  names(res[[i]]) <- c('E1', 'E2')
  
  for(j in 1:2){
    
    res[[i]][[j]] <- vector(mode = 'list', length = 5)
    names(res[[i]][[j]]) <- paste0('S', 1:5)
    
  }
  
}


for(i in 1:4){ # Iteration over the different combination of crop and trait
  
  
  for(j in 1:2){ # iteration over the two experiments
    
    
    for(k in 1:5){ # iteration over the different strategies
      
      load(file.path('data', paste0(data_id[i, j], '.RData')))
      
      exp_des_data <- data[, 1:n_e_d_col]
      
      exp_des_data$col_f <- factor(exp_des_data$col)
      exp_des_data$row_f <- factor(exp_des_data$row)
      exp_des_data$rep <- factor(exp_des_data$rep)
      exp_des_data$block <- factor(exp_des_data$block)
      
      if(k != 5){ # Strategies S1-S5
        
        G_BLUEs <- SpaTemHTP_proc(exp_des_data,
                                  pheno_data = data[, (n_e_d_col+1):dim(data)[2]],
                                  out_det = strat_tab$out_det[k],
                                  miss_imp = strat_tab$miss_imp[k], sp_adj = TRUE,
                                  random = ~ rep +  rep:block + row_f + col_f,
                                  plot = TRUE)
        
        res[[i]][[j]][[k]] <- G_BLUEs
        
      } else { # Strategy S5
        
        
        G_BLUES <- SpaTemHTP_proc(exp_des_data,
                                    pheno_data = data[, (n_e_d_col+1):dim(data)[2]],
                                    single_mixed_model = TRUE,
                                    random = ~ rep +  rep:block + row_f + col_f,
                                    plot = TRUE)
        
        res[[i]][[j]][[k]] <- G_BLUEs
        
      }
      
    }
    
  }
  
}

### Save genotype BLUEs TS

res_TS <- res
save(res_TS, file = './results/res_TS.RData')

### results processing

# construction of the average table. You can start from res_TS
# (!For the moment, it does not contain the result of sc 9. I will add
# them later)

# calculation of the t-test

# Final ouput of this section should be the table 3 of the manuscript

#######

# Assessment of genotype growth pattern 
####################################### 

load(file = './results/res_TS.RData')

n_scen <- 5 # change later to 5 when adding S9

### CP experiments

CP_res <- res_TS[1:2]
CP_res <- list(LA3D_E1 = CP_res[[1]][[1]], LA3D_E2 = CP_res[[1]][[2]],
               PH_E1 = CP_res[[2]][[1]], PH_E2 = CP_res[[2]][[2]])

res_tab_lin <- matrix(NA, n_scen, 4)
res_tab_log <- matrix(NA, n_scen, 4)
CP_log_curve <- vector(mode = 'list', length = 4)
names(CP_log_curve) <- c('LA3D_E1', 'LA3D_E2', 'PH_E1', 'PH_E2')
for(i in 1:4){CP_log_curve[[i]] <- vector(mode = 'list', length = n_scen)
names(CP_log_curve[[i]]) <- paste0('S', 1:n_scen)}

for(i in 1:4){
  
  G_TS <- CP_res[[i]]
  
  for(j in 1:n_scen){
    
    trend_j <- TS_trend(data = G_TS[[j]])
    
    res_tab_lin[j, i] <- round(mean(trend_j$lin_R2, na.rm = TRUE), 2)
    res_tab_log[j, i] <- round(mean(trend_j$log_R2, na.rm = TRUE), 2)
    CP_log_curve[[i]][[j]] <- trend_j$log_val
    
  }
  
}

save(CP_log_curve, file = './results/CP_log_curve.RData')

### SG experiments

SG_res <- res_TS[3:4]
SG_res <- list(LA3D_E1 = SG_res[[1]][[1]], LA3D_E2 = SG_res[[1]][[2]],
               PH_E1 = SG_res[[2]][[1]], PH_E2 = SG_res[[2]][[2]])

res_tab_lin <- matrix(NA, n_scen, 4)
res_tab_log <- matrix(NA, n_scen, 4)
SG_log_curve <- vector(mode = 'list', length = 4)
names(SG_log_curve) <- c('LA3D_E1', 'LA3D_E2', 'PH_E1', 'PH_E2')
for(i in 1:4){SG_log_curve[[i]] <- vector(mode = 'list', length = n_scen)
names(SG_log_curve[[i]]) <- paste0('S', 1:n_scen)}

for(i in 1:4){
  
  G_TS <- SG_res[[i]]
  
  for(j in 1:n_scen){
    
    trend_j <- TS_trend(data = G_TS[[j]])
    
    res_tab_lin[j, i] <- round(mean(trend_j$lin_R2, na.rm = TRUE), 2)
    res_tab_log[j, i] <- round(mean(trend_j$log_R2, na.rm = TRUE), 2)
    SG_log_curve[[i]][[j]] <- trend_j$log_val
    
  }
  
}

save(SG_log_curve, file = './results/SG_log_curve.RData')

###########

# Plots of the genotype BLUEs TS (Figure 2)
###########################################

load(file = './results/res_TS.RData')

### CPE2 LA3D

# raw data

load(file.path('data', "CP_E2_LA3D.RData"))

geno_id <- data$genotype
n_days <- dim(data)[2] - 5
n_geno <- 288

g_means <- matrix(NA, n_geno, n_days)

for(i in 1:dim(g_means)[2]){
  
  g_av <- tapply(X = data[, i + 5], INDEX = geno_id,
                 FUN = function(x) mean(x, na.rm = TRUE))
  g_means[, i] <- g_av
  
}

# calculate the biological trend (logistic curve) on the raw data

trend_raw <- TS_trend(g_means)
dt_trend <- colMeans(x = trend_raw$log_val, na.rm = TRUE)
dt_trend <- data.frame(trait = dt_trend, day = 1:n_days, geno = 'G_av')

dt <- data.frame(geno = rep(names(g_av), n_days), day = rep(1:n_days, each = n_geno),
                 trait = c(g_means))

ggplot(data = dt, aes(x = day, y = trait, group = geno)) +
  geom_line(aes(group = geno)) + geom_line(data = dt_trend, color = 'red', size = 1) +
  ggtitle('CPE2-LA3D raw data')


# Strategies (S5-9)

load(file = './results/CP_log_curve.RData')
log_curve <- CP_log_curve$LA3D_E2

data_sc <- res_TS[[1]][[2]]
title_id <- paste('CPE2-LA3D', paste0('S', 5:9))

for(i in 1:5){
  
  data_i <- data_sc[[i]]
  n_days <- dim(data_i)[2]
  
  dt <- data.frame(geno = rep(rownames(data_i), n_days),
                   day = rep(1:n_days, each = dim(data_i)[1]),
                   trait = c(data_i))
  
  dt_trend <- colMeans(x = log_curve[[i]], na.rm = TRUE)
  dt_trend <- data.frame(trait = dt_trend, day = 1:n_days, geno = 'G_av')
  
  plot <- ggplot(data = dt, aes(x = day, y = trait, group = geno)) +
    geom_line(aes(group = geno)) + geom_line(data = dt_trend, color = 'red', size = 1) +
    ggtitle(title_id[i])
  print(plot)
  
}


### SGE2 PH

# raw data

load(file.path('data', "SG_E2_PH.RData"))

geno_id <- data$genotype
n_days <- dim(data)[2] - 5
n_geno <- 384

g_means <- matrix(NA, n_geno, n_days)

for(i in 1:dim(g_means)[2]){
  
  g_av <- tapply(X = data[, i + 5], INDEX = geno_id,
                         FUN = function(x) mean(x, na.rm = TRUE))
  g_means[, i] <- g_av
  
}

# calculate the biological trend (logistic curve) on the raw data

trend_raw <- TS_trend(g_means)
dt_trend <- colMeans(x = trend_raw$log_val, na.rm = TRUE)
dt_trend <- data.frame(trait = dt_trend, day = 1:n_days, geno = 'G_av')

dt <- data.frame(geno = rep(names(g_av), n_days), day = rep(1:n_days, each = n_geno),
                 trait = c(g_means))

ggplot(data = dt, aes(x = day, y = trait, group = geno)) +
  geom_line(aes(group = geno)) + geom_line(data = dt_trend, color = 'red', size = 1) +
  ggtitle('SGE2-PH raw data')

# Strategies (S5-9)

load(file = './results/SG_log_curve.RData')
log_curve <- SG_log_curve$PH_E2

data_sc <- res_TS[[4]][[2]]
title_id <- paste('SG-PH', paste0('S', 5:9))

for(i in 1:5){
  
  data_i <- data_sc[[i]]
  n_days <- dim(data_i)[2]
  
  dt <- data.frame(geno = rep(rownames(data_i), n_days),
                   day = rep(1:n_days, each = dim(data_i)[1]),
                   trait = c(data_i))
  
  dt_trend <- colMeans(x = log_curve[[i]], na.rm = TRUE)
  dt_trend <- data.frame(trait = dt_trend, day = 1:n_days, geno = 'G_av')
  
  plot <- ggplot(data = dt, aes(x = day, y = trait, group = geno)) +
    geom_line(aes(group = geno)) + geom_line(data = dt_trend, color = 'red', size = 1) +
    ggtitle(title_id[i])
  print(plot)
  
}

########


# Clustering
############

### Soumya (Only the results from the manuscript. Wait for choice to keep clustering)

#########

# Determination of optimal time window (OTW)
############################################

### Soumya (Only the manuscript results)

##########

# Gc x TW model
###############

### Soumya (Only the manuscript results)

##########
