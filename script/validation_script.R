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


devtools::install_github("ICRISAT-GEMS/SpaTemHTP")

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
library(kernlab)
library(fpc)
library(tidyverse)
library(cluster)
library(xts)
library(ecp)
library(ggplot2)
library(lubridate)

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
  
  for(i in 1:n_rep){ # cross-validation loop
    
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
          he_res[data_id[i]][[1]][z, run_ind] <- he
          
          if(!is.null(pheno_pred)){
            
            cor_res <- cor(pheno_pred, data_vs_z$pheno, use = 'complete.obs')
            print(paste('cor(LMM) = ', round(cor_res, 2)))
            pred_res[data_id[i]][[1]][z, run_ind] <- cor_res
            
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
            
            he_res[data_id[i]][[1]][z+4, run_ind] <- he_Sadj <- getHeritability(m_Sadj)
            print(paste('h2(SpATS) = ', round(he_Sadj, 2)))
            
            if(!is.null(pheno_pred)){
              
              cor_res <- cor(pheno_pred$pheno, pheno_pred$predicted.values,
                             use = 'complete.obs')
              print(paste('cor(SpATS) = ', round(cor_res, 2)))
              pred_res[data_id[i]][[1]][z+4, run_ind] <- cor_res
              
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
              
              he_res[data_id[i]][[1]][z+8, run_ind] <- he_Sadj <- getHeritability(m_Sadj)
              print(paste('h2(SpATS) = ', round(he_Sadj, 2)))
              
              if(!is.null(pheno_pred)){
                
                cor_res <- cor(pheno_pred$pheno, pheno_pred$predicted.values,
                               use = 'complete.obs')
                print(paste('cor(SpATS) = ', round(cor_res, 2)))
                pred_res[data_id[i]][[1]][z+8, run_ind] <- cor_res
                
              } 
              
            }
            
          } # end single-step mixed model analysis (S9)
          
        }
        
      } # End different strategies computation
      
    } # End fold loop
    
  } # End CV loop
  
} # End iteration over the different datasets



### results processing

n_scen <- 9

res <- cor_res
res <- as.matrix(res)

t.testRES <- as.data.frame(matrix(NA, nrow = 16, ncol = 2))
rownames(t.testRES) <- c("OL_No", "OL_Yes", "OL_Diff", "OL_pVal",
                         "MI_No", "MI_Yes", "MI_Diff", "MI_pVal",
                         "SP_No", "SP_Yes", "SP_Diff", "SP_pVal",
                         "S8", "S9", "Diff", "pVal")
colnames(t.testRES) <- c("CPE2_LA3D_CORR", "CPE2_LA3D_H2") # 2 columns needed to store RMSE and H2 results


### outlier vs no outliers ###

res_yes <- res[scen_tab$out_det %in% 'yes', ]
res_no <- res[scen_tab$out_det %in% 'no', ]

# OL-t.test
no_OL <- colMeans(res_yes)
OL <- colMeans(res_no)
trt_ind <- factor(c(rep('cont', ncol(res)), rep('trt', ncol(res))))
y <- c(no_OL, OL)
aov.m <- anova(lm(y~trt_ind))
aov_t.test <- t.test(x = no_OL, y = OL)
Sum_of_Sqs <- aov.m$`Sum Sq`
t_stat_val <- aov_t.test$statistic
aov_p_val <- aov_t.test$p.value
grp_means <- aov_t.test$estimate
t.testRES[1:4, 1] <- round(c(grp_means[1], grp_means[2], (grp_means[1]-grp_means[2]), aov_p_val), 2)


### miss imp vs no miss imp ###

res_yes <- res[scen_tab$miss_imp %in% 'yes', ]
res_no <- res[scen_tab$miss_imp %in% 'no', ]

# MICE-t.test
no_MI <- colMeans(res_yes)
MI <- colMeans(res_no)
trt_ind <- factor(c(rep('cont', ncol(res)), rep('trt', ncol(res))))
y <- c(no_MI, MI)
aov.m <- anova(lm(y~trt_ind))
aov_t.test <- t.test(x = no_MI, y = MI)
Sum_of_Sqs <- aov.m$`Sum Sq`
t_stat_val <- aov_t.test$statistic
aov_p_val <- aov_t.test$p.value
grp_means <- aov_t.test$estimate
t.testRES[5:8, 1] <- round(c(grp_means[1], grp_means[2], (grp_means[1]-grp_means[2]), aov_p_val), 2)


### spatial corr vs no spatial corr ###

res_yes <- res[scen_tab$spat_adj %in% 'yes', ]
res_no <- res[scen_tab$spat_adj %in% 'no', ]

# SpATS-t.test
no_SP <- colMeans(res_yes)
SP <- colMeans(res_no)
trt_ind <- factor(c(rep('cont', ncol(res)), rep('trt', ncol(res))))
y <- c(no_SP, SP)
aov.m <- anova(lm(y~trt_ind))
aov_t.test <- t.test(x = no_SP, y = SP)
Sum_of_Sqs <- aov.m$`Sum Sq`
t_stat_val <- aov_t.test$statistic
aov_p_val <- aov_t.test$p.value
grp_means <- aov_t.test$estimate
t.testRES[9:12, 1] <- round(c(grp_means[1], grp_means[2], (grp_means[1]-grp_means[2]), aov_p_val), 2)


### SpATS vs. S9 t.test ###
S8 <- res[8, ]
S9 <- res[9, ]
trt_ind <- factor(c(rep('cont', length(na.approx(S8))), rep('trt', length(na.approx(S9)))))
y <- c(na.approx(S8), na.approx(S9))
aov.m <- anova(lm(y~trt_ind))
aov_t.test <- t.test(x = S8, y = S9)
Sum_of_Sqs <- aov.m$`Sum Sq`
t_stat_val <- aov_t.test$statistic
aov_p_val <- aov_t.test$p.value
grp_means <- aov_t.test$estimate
t.testRES[13:16, 1] <- round(c(grp_means[1], grp_means[2], (grp_means[1]-grp_means[2]), aov_p_val), 2)




# heritability
######################
res <- he_res
res <- as.matrix(res)


### outlier vs no outliers ###

res_yes <- res[scen_tab$out_det %in% 'yes', ]
res_no <- res[scen_tab$out_det %in% 'no', ]

# OL-t.test
no_OL <- colMeans(res_yes)
OL <- colMeans(res_no)
trt_ind <- factor(c(rep('cont', ncol(res)), rep('trt', ncol(res))))
y <- c(no_OL, OL)
aov.m <- anova(lm(y~trt_ind))
aov_t.test <- t.test(x = no_OL, y = OL)
Sum_of_Sqs <- aov.m$`Sum Sq`
t_stat_val <- aov_t.test$statistic
aov_p_val <- aov_t.test$p.value
grp_means <- aov_t.test$estimate
t.testRES[1:4, 2] <- round(c(grp_means[1], grp_means[2], (grp_means[1]-grp_means[2]), aov_p_val), 2)


### miss imp vs no miss imp ###

res_yes <- res[scen_tab$miss_imp %in% 'yes', ]
res_no <- res[scen_tab$miss_imp %in% 'no', ]

# MICE-t.test
no_MI <- colMeans(res_yes)
MI <- colMeans(res_no)
trt_ind <- factor(c(rep('cont', ncol(res)), rep('trt', ncol(res))))
y <- c(no_MI, MI)
aov.m <- anova(lm(y~trt_ind))
aov_t.test <- t.test(x = no_MI, y = MI)
Sum_of_Sqs <- aov.m$`Sum Sq`
t_stat_val <- aov_t.test$statistic
aov_p_val <- aov_t.test$p.value
grp_means <- aov_t.test$estimate
t.testRES[5:8, 2] <- round(c(grp_means[1], grp_means[2], (grp_means[1]-grp_means[2]), aov_p_val), 2)


### spatial corr vs no spatial corr ###

res_yes <- res[scen_tab$spat_adj %in% 'yes', ]
res_no <- res[scen_tab$spat_adj %in% 'no', ]

# SpATS-t.test
no_SP <- colMeans(res_yes)
SP <- colMeans(res_no)
trt_ind <- factor(c(rep('cont', ncol(res)), rep('trt', ncol(res))))
y <- c(no_SP, SP)
aov.m <- anova(lm(y~trt_ind))
aov_t.test <- t.test(x = no_SP, y = SP)
Sum_of_Sqs <- aov.m$`Sum Sq`
t_stat_val <- aov_t.test$statistic
aov_p_val <- aov_t.test$p.value
grp_means <- aov_t.test$estimate
t.testRES[9:12, 2] <- round(c(grp_means[1], grp_means[2], (grp_means[1]-grp_means[2]), aov_p_val), 2)


### SpATS vs. S9 t.test ###
S8 <- res[8, ]
S9 <- res[9, ]
trt_ind <- factor(c(rep('cont', length(na.approx(S8))), rep('trt', length(na.approx(S9)))))
y <- c(na.approx(S8), na.approx(S9))
aov.m <- anova(lm(y~trt_ind))
aov_t.test <- t.test(x = S8, y = S9)
Sum_of_Sqs <- aov.m$`Sum Sq`
t_stat_val <- aov_t.test$statistic
aov_p_val <- aov_t.test$p.value
grp_means <- aov_t.test$estimate
t.testRES[13:16, 2] <- round(c(grp_means[1], grp_means[2], (grp_means[1]-grp_means[2]), aov_p_val), 2)

write.csv(t.testRES,  paste0(output.loc, '/CPE2_LA3D_t.testRES.csv'))

# Final ouput of this section is the table 3 of the manuscript


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
# c) level 3: within level 2, four strategies (S1-4)

res <- vector(mode = 'list', length = dim(data_id)[1])
names(res) <- c('CP_LA3D', 'CP_PH', 'SG_LA3D', 'SG_PH')

for(i in 1:4){
  
  res[[i]] <- vector(mode = 'list', length = 2)
  names(res[[i]]) <- c('E1', 'E2')
  
  for(j in 1:2){
    
    res[[i]][[j]] <- vector(mode = 'list', length = 4)
    names(res[[i]][[j]]) <- paste0('S', 1:4)
    
  }
  
}


for(i in 1:4){ # Iteration over the different combination of crop and trait
  
  
  for(j in 1:2){ # iteration over the two experiments
    
    
    for(k in 1:4){ # iteration over the different strategies
      
      load(file.path('data', paste0(data_id[i, j], '.RData')))
      
      exp_des_data <- data[, 1:n_e_d_col]
      
      exp_des_data$col_f <- factor(exp_des_data$col)
      exp_des_data$row_f <- factor(exp_des_data$row)
      exp_des_data$rep <- factor(exp_des_data$rep)
      exp_des_data$block <- factor(exp_des_data$block)
  
      G_BLUEs <- SpaTemHTP_proc(exp_des_data,
                                pheno_data = data[, (n_e_d_col+1):dim(data)[2]],
                                out_det = strat_tab$out_det[k],
                                miss_imp = strat_tab$miss_imp[k], sp_adj = TRUE,
                                random = ~ rep +  rep:block + row_f + col_f,
                                plot = TRUE)
      
      res[[i]][[j]][[k]] <- G_BLUEs
      
          
    }
    
  }
  
}

Sys.time()

### Save genotype BLUEs TS

res_TS <- res
save(res_TS, file = './results/res_TS.RData')


### results processing

n_scen <- 9

betEx_t.testRES <- as.data.frame(matrix(NA, nrow = 16, ncol = 4)
rownames(betEx_t.testRES) <- c("OL_No", "OL_Yes", "OL_Diff", "OL_pVal",
                         "MI_No", "MI_Yes", "MI_Diff", "MI_pVal",
                         "SP_No", "SP_Yes", "SP_Diff", "SP_pVal",
                         "S8", "S9", "Diff", "pVal")
colnames(betEx_t.testRES) <- c('CP_LA3D', 'CP_PH', 'SG_LA3D', 'SG_PH') # 4 columns for each of the dataset

for(i in 1:4)
{

  res <- res_TS[[i]] 
  res <- as.matrix(res)

	### outlier vs no outliers ###

	res_yes <- res[scen_tab$out_det %in% 'yes', ]
	res_no <- res[scen_tab$out_det %in% 'no', ]

	# OL-t.test
	no_OL <- colMeans(res_yes)
	OL <- colMeans(res_no)
	trt_ind <- factor(c(rep('cont', ncol(res)), rep('trt', ncol(res))))
	y <- c(no_OL, OL)
	aov.m <- anova(lm(y~trt_ind))
	aov_t.test <- t.test(x = no_OL, y = OL)
	Sum_of_Sqs <- aov.m$`Sum Sq`
	t_stat_val <- aov_t.test$statistic
	aov_p_val <- aov_t.test$p.value
	grp_means <- aov_t.test$estimate
	betEx_t.testRES[1:4, i] <- round(c(grp_means[1], grp_means[2], (grp_means[1]-grp_means[2]), aov_p_val), 2)


	### miss imp vs no miss imp ###

	res_yes <- res[scen_tab$miss_imp %in% 'yes', ]
	res_no <- res[scen_tab$miss_imp %in% 'no', ]

	# MICE-t.test
	no_MI <- colMeans(res_yes)
	MI <- colMeans(res_no)
	trt_ind <- factor(c(rep('cont', ncol(res)), rep('trt', ncol(res))))
	y <- c(no_MI, MI)
	aov.m <- anova(lm(y~trt_ind))
	aov_t.test <- t.test(x = no_MI, y = MI)
	Sum_of_Sqs <- aov.m$`Sum Sq`
	t_stat_val <- aov_t.test$statistic
	aov_p_val <- aov_t.test$p.value
	grp_means <- aov_t.test$estimate
	betEx_t.testRES[5:8, i] <- round(c(grp_means[1], grp_means[2], (grp_means[1]-grp_means[2]), aov_p_val), 2)


	### spatial corr vs no spatial corr ###

	res_yes <- res[scen_tab$spat_adj %in% 'yes', ]
	res_no <- res[scen_tab$spat_adj %in% 'no', ]

	# SpATS-t.test
	no_SP <- colMeans(res_yes)
	SP <- colMeans(res_no)
	trt_ind <- factor(c(rep('cont', ncol(res)), rep('trt', ncol(res))))
	y <- c(no_SP, SP)
	aov.m <- anova(lm(y~trt_ind))
	aov_t.test <- t.test(x = no_SP, y = SP)
	Sum_of_Sqs <- aov.m$`Sum Sq`
	t_stat_val <- aov_t.test$statistic
	aov_p_val <- aov_t.test$p.value
	grp_means <- aov_t.test$estimate
	betEx_t.testRES[9:12, i] <- round(c(grp_means[1], grp_means[2], (grp_means[1]-grp_means[2]), aov_p_val), 2)


	### SpATS vs. S9 t.test ###
	S8 <- res[8, ]
	S9 <- res[9, ]
	trt_ind <- factor(c(rep('cont', length(na.approx(S8))), rep('trt', length(na.approx(S9)))))
	y <- c(na.approx(S8), na.approx(S9))
	aov.m <- anova(lm(y~trt_ind))
	aov_t.test <- t.test(x = S8, y = S9)
	Sum_of_Sqs <- aov.m$`Sum Sq`
	t_stat_val <- aov_t.test$statistic
	aov_p_val <- aov_t.test$p.value
	grp_means <- aov_t.test$estimate
	betEx_t.testRES[13:16, i] <- round(c(grp_means[1], grp_means[2], (grp_means[1]-grp_means[2]), aov_p_val), 2)

}
save(betEx_t.testRES, file = './results/betEx_t.testRES.RData')

#######

# Plots of the genotype BLUEs TS (Figure 2)
###########################################

# Vincent does that. In any case, there is already a ploting option integrate
# in the package function SpaTemHTP_proc

########

# Assessment of genotype growth pattern (linear trend fit)
#########################################################

### Soumya

########

# Determination of optimal time window (OTW)
############################################

# For each crop-type in an experiment find the OTW, e.g. CP_E2 

cpaRES <- OTW <- list()

for (i in 1:2){

  cpaip <- as.data.frame(res[[i]][[2]][[4]])
  
  deflt <- 3
  
  blueKKmeans<-kkmeans(as.matrix(cpaip), deflt, 
                       kernel="polydot",alg="kkmeans",p=1, na.action=na.omit)
  blue.res<-as.data.frame(centers(blueKKmeans))
  
  colnames(blue.res)<-colnames(cpaip)
  
  clust_diff<-matrix(nrow = 1, ncol = ncol(blue.res))
  
  for(i in 1:ncol(blue.res))
  {
    # clust_diff[1,i]<-(abs(blue.res[1,i]-blue.res[2,i]) + abs(blue.res[1,i]-blue.res[3,i]) + abs(blue.res[3,i]-blue.res[2,i]))
    clust_diff[1,i]<-(sum(dist(blue.res[,i], method = "euclidean")))
  }
  
  colnames(clust_diff) <- colnames(blue.res)
  
  
  ###### Start CPA for OTW identification ######
  
  ip.cpa<-cbind(t(clust_diff), Herit.Res)
  colnames(ip.cpa)[1]<-"BLUE_CD"
  colnames(ip.cpa)[2]<-"h2"
  head(ip.cpa)
  
  ip.cpa.ts <- xts(ip.cpa, order.by=as.Date(rownames(ip.cpa), "%Y-%m-%d"))
    
  ecp.ph<-e.cp3o(Z=ip.cpa.ts, K=4, minsize=3, alpha=1, verbose=FALSE)
  
  E<-ecp.ph$estimates
  
  dates<-rownames(ip.cpa)
  bp.date<-xts(dates[E], order.by = as.Date(dates[E], "%Y-%m-%d"))
  
  cpaRES[[i]] <- list(ip.cpa.ts=ip.cpa.ts, bp.date=bp.date)
  
  range01 <- function(x) {(x-min(x))/(max(x) - min(x))}
  
  # Find OTW
  TWs <- length(bp.date)+1
  TWmetr <- as.data.frame(matrix(NA, nrow = 3, ncol = TWs))
  colnames(TWmetr) <- paste0("TW-", 1:TWs)
  rownames(TWmetr) <- c("medCD", "slpCD", "medH2")
  
  for(i in 1:TWs){
    
    if(i == 1){
      r.ind <- which(dates %in% bp.date[i])
      tmp.tw <- ip.cpa.ts[1:r.ind, ]
      
          } else if (i == TWs){
            r.ind <- which(dates %in% bp.date[i-1])
            tmp.tw <- ip.cpa.ts[r.ind:nrow(ip.cpa.ts), ]
        
             } else {
               r.ind1 <- which(dates %in% bp.date[i-1])
               r.ind2 <- which(dates %in% bp.date[i])
               tmp.tw <- ip.cpa.ts[(r.ind1) : (r.ind2-1), ]
               
               } # end if-else
    
    # get the median of cluster-distance
    TWmetr[1 ,i] <- round(mean(tmp.tw$BLUE_CD), 2)
    # get the slope of cluster-distance
    l.mod <- lm(tmp.tw$BLUE_CD ~ c(1:dim(tmp.tw)[1]))
    l.mod.st <- summary(l.mod)
    TWmetr[2 ,i] <- l.mod.st$coefficients[2, 1]
    # get the median of heritability
    TWmetr[3 ,i] <- round(mean(tmp.tw$h2), 2)
    
  } # end for loop
  
  TWmetr.sc <- as.data.frame(t(apply(TWmetr, 1, range01)))
  OTWid <- which.max(apply(TWmetr.sc[c(1,3), ], 2, sum))
  
  BLUEs <- x$BLUEs
  
  if (OTWid == 1) {
    c2 <- bp.date[(OTWid)]
    col2 <- which(colnames(BLUEs) %in% as.character(c2))
    OTW[[i]] <- BLUEs[ ,1:(col2)]
  
    } else if (OTWid == length(bp.date)) {
      c1 <- bp.date[(OTWid)]
      col1 <- which(colnames(BLUEs) %in% as.character(c1))
      OTW[[i]] <- BLUEs[ ,(col1:ncol(BLUEs))]
    
      } else {
        c1 <- bp.date[(OTWid-1)]
        c2 <- bp.date[(OTWid)]
        
        col1 <- which(colnames(BLUEs) %in% as.character(c1))
        col2 <- which(colnames(BLUEs) %in% as.character(c2))
        
        OTW[[i]] <- BLUEs[ ,col1:(col2-1)]
        
      }

save(OTW, file = './results/OTW.RData')

##########


# Clustering
###############

# Cluster using both the traits, LA3D and PH
clusDF <- rowmeans(OTW[[1]]) + rowmeans(OTW[[2]])
pc <- princomp(clusDF)
plot(pc)

# Scale
data2 <- data.frame(scale(clusDF))
# Verify variance is uniform
plot(sapply(data2, var))

# Proceed with principal components
pc <- princomp(data2)
plot(pc)
plot(pc, type='l')
summary(pc)

cp.res <- pc$scores

# First for principal components
comp <- data.frame(pc$scores[,1:2])
# Plot
plot(comp, pch=16, col=rgb(0,0,0,0.5))
# Apply k-means with k=4
k <- kmeans(comp, 3, nstart=25, iter.max=1000)
palette(alpha(brewer.pal(9,'Set1'), 0.5))
plot(comp, col=k$clust, pch=16)

cp.clusID <- as.data.frame(k$cluster)

save(as.data.frame(cbind(cp.res, cp.clusID, rownames(OTW[[1]]))), file = './results/FinalClusters.RData')


##########



# Gc x TW model
###############
  	
for (i in 1:length(cpaRES)) { # For each CP_E2 trait

allData <- as.data.frame(res[[i]][[2]][[4]])

    bp.date <- cpaRES[[i]]$bp.date
    TWs <- length(bp.date)+1
        
    aovDF <- data.frame(Trait = c(unlist(allData[,1:bp.date[1]])), 
	                Geno = rep(cp.clusID, ncol(allData[,1:bp.date[1]]), 
                        TW = rep(1, nrow(allData)*ncol(allData[,1:bp.date[1]]))
  	  
    for(j in 2:length(bp.date)-1) {
	
	tmp <- data.frame(Trait = c(unlist(allData[ ,bp.date[j]:bp.date[j+1]])), 
	                  Geno = rep(cp.clusID, ncol(allData[ ,bp.date[j]:bp.date[j+1]), 
                          TW = rep(j, nrow(allData)*ncol(allData[ ,bp.date[j]:bp.date[j+1]))

	aovDF <- as.data.frame(rbind(aovDF, tmp)) }

	tmp1 <- data.frame(Trait = c(unlist(allData[ ,bp.date[length(bp.date)]:ncol[allData]])), 
                           Geno = rep(cp.clusID, ncol(allData[ ,bp.date[length(bp.date)]:ncol[allData]], 
                           TW = rep(length(bp.date)+1, nrow(allData)*ncol(allData[ ,bp.date[length(bp.date)]:ncol[allData]))
	aovDF <- as.data.frame(rbind(aovDF, tmp1))

        aov.model<-aov(Trait ~ Geno + TW + Geno*TW, data = aovDF)

summary(aov.model)
}

##########
