

# code to compute test statistics of DM tests
# results are saved in a specific data frame format to simplify plotting later on

rm(list=ls())


# parameters to run
# Exemplarily 45 parameter and rho/rho0 combinations, 36 of them presented in the paper

input_par <- as.data.frame(matrix(NA, nrow =45, ncol = 9))
names(input_par) <- c("d", "gev0loc_obs", "gev0shape_obs", "gev0scale_obs",
                      "gev0loc_ens", "gev0shape_ens", "gev0scale_ens", "rho0","rho")


for(i in 1:18){
  input_par[i, which(is.element(names(input_par),c("gev0loc_obs", "gev0shape_obs", "gev0scale_obs")))] <- c(0, -0.1, 1)
}
for(i in 19:36){
  input_par[i, which(is.element(names(input_par),c("gev0loc_obs", "gev0shape_obs", "gev0scale_obs")))] <- c(1, 0.3, 1)
}
for(i in 37:45){
  input_par[i, which(is.element(names(input_par),c("gev0loc_obs", "gev0shape_obs", "gev0scale_obs")))] <- c(0, 0, 1)
}


for(i in 1:9){
  input_par[i, which(is.element(names(input_par),c("gev0loc_ens", "gev0shape_ens", "gev0scale_ens")))] <- c(1, 0, 0.2)
}
for(i in 10:18){
  input_par[i, which(is.element(names(input_par),c("gev0loc_ens", "gev0shape_ens", "gev0scale_ens")))] <- c(0, 0, 2)
}
for(i in 19:27){
  input_par[i, which(is.element(names(input_par),c("gev0loc_ens", "gev0shape_ens", "gev0scale_ens")))] <- c(1, 0, 0.2)
}
for(i in 28:36){
  input_par[i, which(is.element(names(input_par),c("gev0loc_ens", "gev0shape_ens", "gev0scale_ens")))] <- c(0, 0, 2)
}
for(i in 37:45){
  input_par[i, which(is.element(names(input_par),c("gev0loc_ens", "gev0shape_ens", "gev0scale_ens")))] <- c(0, 0, 1)
}


input_par$rho0 <- rep(c(0.25,0.5,0.75),each=3)
input_par$d <- 4
input_par$rho <- rep(c(0.25,0.5,0.75),times=3)



# Use input_par matrix from above now
df_raw <- data.frame(input_par)
df_raw$simID <- 1:nrow(df_raw)

flist <- list.files("/path/to/Rdata-files/")

# Load one of the Rdata files exemplarily to have object res in R workspace
load("/path/to/example/Rdata-file.Rdata")

input_models <- names(res$es_list)
input_scores <- names(res)

df_use <- as.data.frame(df_raw[1,])
df_use$model <- as.character("a")
df_use$score <- as.character("a")

model_score_grid <- expand.grid(input_models, input_scores)

for(i in 1:nrow(df_raw)){
  df_use[((i-1)*nrow(model_score_grid)+1):(i*nrow(model_score_grid)),] <- df_raw[i,]
  df_use[((i-1)*nrow(model_score_grid)+1):(i*nrow(model_score_grid)),]$model <- as.character(model_score_grid$Var1)
  df_use[((i-1)*nrow(model_score_grid)+1):(i*nrow(model_score_grid)),]$score <- as.character(model_score_grid$Var2)
}

df100 <- data.frame(cbind(zoo::coredata(df_use)[rep(seq(nrow(df_use)),100),]))
df100$value <- NA
head(df100)

library(forecast) # for DM test function

for(filename in flist){

  load(paste0("/path/to/Rdata-files/", filename))
  ID <- as.numeric(as.numeric(strsplit(strsplit(filename, "setting_")[[1]][2], ".Rdata")))
  print(which(flist == filename))

  for(this_model in input_models){
    for(this_score in input_scores){
      ind <- which(df100$simID == ID & df100$model == this_model & df100$score == this_score)

      # deal with CRPS specifically (use only first dimension, not all 4 recorded ones)

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


save(df100, file = "/path/to/save/output/df_gev0.Rdata")



