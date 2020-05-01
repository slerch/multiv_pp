# code to compute test statistics of DM tests
# results are saved in a specific data frame format to simplify plotting later on

rm(list=ls())

# parameters 
input_rho0 <- c(0.1, 0.25, 0.5, 0.75, 0.9)
input_eps <- c(1)
input_sigma <-  1
input_rho <- c(0.1, 0.25, 0.5, 0.75, 0.9)
input_d <- 5
input_par <- expand.grid(input_rho0, input_eps, input_sigma, input_rho, input_d)
names(input_par) <- c("rho0", "eps", "sigma", "rho", "d")

df_raw <- data.frame(input_par)
df_raw$simID <- 1:nrow(df_raw)

flist <- list.files("/path/to/Rdata-files/")
existing <- as.numeric(sapply(flist, FUN = function(x) as.numeric(strsplit(strsplit(x, "_setting3A_")[[1]][2], ".Rdata"))))

df <- df_raw[which(is.element(df_raw$simID, existing)),] 

load("/path/to/example/Rdata-output-file.Rdata")

input_models <- names(res$es_list)
input_scores <- names(res)

df_use <- as.data.frame(df[1,])
df_use$model <- as.character("a")
df_use$score <- as.character("a")

model_score_grid <- expand.grid(input_models, input_scores)

for(i in 1:nrow(df)){
  df_use[((i-1)*nrow(model_score_grid)+1):(i*nrow(model_score_grid)),] <- df[i,]
  df_use[((i-1)*nrow(model_score_grid)+1):(i*nrow(model_score_grid)),]$model <- as.character(model_score_grid$Var1)
  df_use[((i-1)*nrow(model_score_grid)+1):(i*nrow(model_score_grid)),]$score <- as.character(model_score_grid$Var2)
}

df100 <- data.frame(cbind(zoo::coredata(df_use)[rep(seq(nrow(df_use)),100),]))
df100$value <- NA
head(df100)

library(forecast) # for DM test function

for(filename in flist){
  
  load(paste0("/path/to/Rdata-files/", filename))
  ID <- as.numeric(as.numeric(strsplit(strsplit(filename, "_setting3A_")[[1]][2], ".Rdata")))
  print(which(flist == filename))
  
  for(this_model in input_models){
    for(this_score in input_scores){
      ind <- which(df100$simID == ID & df100$model == this_model & df100$score == this_score)  
      
      # deal with CRPS specifically (use only first dimension, not all 5 recorded ones)
      
      if(this_score == "crps_list"){
        dm_teststat_vec <- rep(NA, 100)
        for(MC_rep in 1:100){
          tmp <- NA
          tryDM <- try(tmp_DM <- dm.test(e1 = res[[which(input_scores == this_score)]][[which(input_models == this_model)]][,,1][MC_rep,],
                                         e2 = res[[which(input_scores == this_score)]][[which(input_models == "ecc.q")]][,,1][MC_rep,],
                                         h = 1, power = 1), silent = TRUE)
          if(class(tryDM) != "try-error"){
            tmp <- tmp_DM$statistic
          } else{
            tmp <- 0
          }
          dm_teststat_vec[MC_rep] <- tmp
        }
      } else{
        dm_teststat_vec <- rep(NA, 100)
        for(MC_rep in 1:100){
          tmp <- NA
          tryDM <- try(tmp_DM <- dm.test(e1 = res[[which(input_scores == this_score)]][[which(input_models == this_model)]][MC_rep,],
                                         e2 = res[[which(input_scores == this_score)]][[which(input_models == "ecc.q")]][MC_rep,],
                                         h = 1, power = 1), silent = TRUE)
          if(class(tryDM) != "try-error"){
            tmp <- tmp_DM$statistic
          } else{
            tmp <- 0
          }
          dm_teststat_vec[MC_rep] <- tmp
        }
      }
      
      
      df100$value[ind] <- dm_teststat_vec
    }
  }
  
}

save(df100, file = "/path/to/save/output/df_mvNormal_StructuralChange.Rdata")

