# Note: only the figure for high rho/rho0 and a = 0.7 (last plot below) is part of the main paper, the rest of the figures is presented in the supplemental material.

rm(list=ls())

library(ggplot2)
library(gridExtra)

load("/path/to/save/output/df_mvNormal_StructuralChange_B.Rdata")

df <- df100; rm(df100)
df1 <- subset(df, eps == 1)
input_rho0 <- unique(df1$rho0)
input_rho <- unique(df1$rho)
input_sigma <-  unique(df1$sigma)
input_a <- unique(df1$a)

df1$value <- (-1)*df1$value
df0 <- df1
input_scores <- unique(df1$score)

plot_folder <- "/path/to/save/figures/"

df1_save <- df0 
df1 <- subset(df1_save, model != "ens")
df2 <- subset(df1, model != "emos.q")


##
## low rho0 / rho0s
##

this_rho0 <- 0.25

df2 <- subset(df2, rho0 == this_rho0)
df2 <- subset(df2, rho == 0.2 | rho == 0.25 | rho == 0.3)

mypal <- colorspace::rainbow_hcl(5)
mypal_use <- c("decc.q" = mypal[1],
               "ecc.q" = mypal[2],
               "ecc.s" = mypal[3],
               "gca" = mypal[4],
               "ssh" = mypal[5])

df2$model <- factor(df2$model, levels = c("decc.q", "ecc.q", "ecc.s", "gca", "ssh"))

df2 <- subset(df2, model != "ecc.q")

model_vec <- c("dECC", "ECC-S", "GCA", "SSh")

this_sigma <- input_sigma[2]

for(this_sigma in input_sigma){
  
  this_score <- "es_list"
  
  dfplot <- subset(df2, sigma == this_sigma & score == this_score)
  dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 
  
  p1 <- ggplot(dfplot, aes(model, value, colour = model)) + 
    geom_rect(data = subset(dfplot, equalRhos == 1), color = "black", size = 2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf) 
  p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25) 
  p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
  p1 <- p1 + facet_grid(rows = vars(a), cols = vars(rho), 
                        labeller = label_bquote(rows = a==.(a), 
                                                cols = rho==.(rho)))
  sigval <- round(this_sigma, 2)
  scval <- strsplit(this_score, split = "_")[[1]][1]
  if(scval == "es"){
    p1 <- p1 + ggtitle(bquote(list("Energy Score", sigma == .(sigval), rho[0] ==.(this_rho0))))
  } 
  
  p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
  p1 <- p1 + xlab("Model") + ylab("DM test statistic") 
  p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", label = model_vec) 
  p1 <- p1 + scale_x_discrete(label = model_vec)
  p1save <- p1
  
  ## VS
  
  this_score <- "vs1_list"
  
  dfplot <- subset(df2, sigma == this_sigma & score == this_score)
  dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 
  
  p1 <- ggplot(dfplot, aes(model, value, colour = model)) + 
    geom_rect(data = subset(dfplot, equalRhos == 1), color = "black", size = 2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf) 
  p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25) 
  p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
  p1 <- p1 + facet_grid(rows = vars(a), cols = vars(rho), 
                        labeller = label_bquote(rows = a==.(a), 
                                                cols = rho==.(rho)))
  sigval <- round(this_sigma, 2)
  scval <- strsplit(this_score, split = "_")[[1]][1]
  if(scval == "vs1"){
    p1 <- p1 + ggtitle(bquote(list("Variogram Score", sigma == .(sigval), rho[0] ==.(this_rho0))))
  } 
  
  p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
  p1 <- p1 + xlab("Model") + ylab("DM test statistic") 
  p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", label = model_vec) 
  p1 <- p1 + scale_x_discrete(label = model_vec)
  p2save <- p1
  
  pdf(paste0(plot_folder, "S3B_lowRho0_sigma", round(this_sigma,2), ".pdf"),  width = 7.5, height = 11, pointsize = 12)
  grid.arrange(p1save, p2save, ncol = 1)
  dev.off()
  
}


##
## medium rho0 / rho0s
##

this_rho0 <- 0.5

df2 <- subset(df1, rho0 == this_rho0)
df2 <- subset(df2, rho == 0.4 | rho == 0.45 | rho == 0.5 | rho == 0.55 | rho == 0.6)

mypal <- colorspace::rainbow_hcl(5)
mypal_use <- c("decc.q" = mypal[1],
               "ecc.q" = mypal[2],
               "ecc.s" = mypal[3],
               "gca" = mypal[4],
               "ssh" = mypal[5])

df2$model <- factor(df2$model, levels = c("decc.q", "ecc.q", "ecc.s", "gca", "ssh"))

df2 <- subset(df2, model != "ecc.q")

model_vec <- c("dECC", "ECC-S", "GCA", "SSh")

this_sigma <- input_sigma[2]

for(this_sigma in input_sigma){
  
  this_score <- "es_list"
  
  dfplot <- subset(df2, sigma == this_sigma & score == this_score)
  dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 
  
  p1 <- ggplot(dfplot, aes(model, value, colour = model)) + 
    geom_rect(data = subset(dfplot, equalRhos == 1), color = "black", size = 2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf) 
  p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25) 
  p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
  p1 <- p1 + facet_grid(rows = vars(a), cols = vars(rho), 
                        labeller = label_bquote(rows = a==.(a), 
                                                cols = rho==.(rho)))
  sigval <- round(this_sigma, 2)
  scval <- strsplit(this_score, split = "_")[[1]][1]
  if(scval == "es"){
    p1 <- p1 + ggtitle(bquote(list("Energy Score", sigma == .(sigval), rho[0] ==.(this_rho0))))
  } 
  
  p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
  p1 <- p1 + xlab("Model") + ylab("DM test statistic") 
  p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", label = model_vec) 
  p1 <- p1 + scale_x_discrete(label = model_vec)
  p1save <- p1
  
  ## VS
  
  this_score <- "vs1_list"
  
  dfplot <- subset(df2, sigma == this_sigma & score == this_score)
  dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 
  
  p1 <- ggplot(dfplot, aes(model, value, colour = model)) + 
    geom_rect(data = subset(dfplot, equalRhos == 1), color = "black", size = 2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf) 
  p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25) 
  p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
  p1 <- p1 + facet_grid(rows = vars(a), cols = vars(rho), 
                        labeller = label_bquote(rows = a==.(a), 
                                                cols = rho==.(rho)))
  sigval <- round(this_sigma, 2)
  scval <- strsplit(this_score, split = "_")[[1]][1]
  if(scval == "vs1"){
    p1 <- p1 + ggtitle(bquote(list("Variogram Score", sigma == .(sigval), rho[0] ==.(this_rho0))))
  } 
  
  p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
  p1 <- p1 + xlab("Model") + ylab("DM test statistic") 
  p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", label = model_vec) 
  p1 <- p1 + scale_x_discrete(label = model_vec)
  p2save <- p1
  
  pdf(paste0(plot_folder, "S3B_medRho0_sigma", round(this_sigma,2), ".pdf"),  width = 11, height = 11, pointsize = 12)
  grid.arrange(p1save, p2save, ncol = 1)
  dev.off()
  
}



##
## high rho0 / rho0s
##

this_rho0 <- 0.75

df2 <- subset(df1, rho0 == this_rho0)
df2 <- subset(df2, rho == 0.7 | rho == 0.75 | rho == 0.8)

mypal <- colorspace::rainbow_hcl(5)
mypal_use <- c("decc.q" = mypal[1],
               "ecc.q" = mypal[2],
               "ecc.s" = mypal[3],
               "gca" = mypal[4],
               "ssh" = mypal[5])

df2$model <- factor(df2$model, levels = c("decc.q", "ecc.q", "ecc.s", "gca", "ssh"))

df2 <- subset(df2, model != "ecc.q")

model_vec <- c("dECC", "ECC-S", "GCA", "SSh")

this_sigma <- input_sigma[2]

for(this_sigma in input_sigma){
  
  this_score <- "es_list"
  
  dfplot <- subset(df2, sigma == this_sigma & score == this_score)
  dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 
  
  p1 <- ggplot(dfplot, aes(model, value, colour = model)) + 
    geom_rect(data = subset(dfplot, equalRhos == 1), color = "black", size = 2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf) 
  p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25) 
  p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
  p1 <- p1 + facet_grid(rows = vars(a), cols = vars(rho), 
                        labeller = label_bquote(rows = a==.(a), 
                                                cols = rho==.(rho)))
  sigval <- round(this_sigma, 2)
  scval <- strsplit(this_score, split = "_")[[1]][1]
  if(scval == "es"){
    p1 <- p1 + ggtitle(bquote(list("Energy Score", sigma == .(sigval), rho[0] ==.(this_rho0))))
  } 
  
  p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
  p1 <- p1 + xlab("Model") + ylab("DM test statistic") 
  p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", label = model_vec) 
  p1 <- p1 + scale_x_discrete(label = model_vec)
  p1save <- p1
  
  ## VS
  
  this_score <- "vs1_list"
  
  dfplot <- subset(df2, sigma == this_sigma & score == this_score)
  dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 
  
  p1 <- ggplot(dfplot, aes(model, value, colour = model)) + 
    geom_rect(data = subset(dfplot, equalRhos == 1), color = "black", size = 2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf) 
  p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25) 
  p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
  p1 <- p1 + facet_grid(rows = vars(a), cols = vars(rho), 
                        labeller = label_bquote(rows = a==.(a), 
                                                cols = rho==.(rho)))
  sigval <- round(this_sigma, 2)
  scval <- strsplit(this_score, split = "_")[[1]][1]
  if(scval == "vs1"){
    p1 <- p1 + ggtitle(bquote(list("Variogram Score", sigma == .(sigval), rho[0] ==.(this_rho0))))
  } 
  
  p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
  p1 <- p1 + xlab("Model") + ylab("DM test statistic") 
  p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", label = model_vec) 
  p1 <- p1 + scale_x_discrete(label = model_vec)
  p2save <- p1
  
  pdf(paste0(plot_folder, "S3B_highRho0_sigma", round(this_sigma,2), ".pdf"),  width = 7.5, height = 11, pointsize = 12)
  grid.arrange(p1save, p2save, ncol = 1)
  dev.off()
  
}


##
## high rho/rho0, but only a = 0.7 [shown in the paper, rest in the supplemental material]
##


this_rho0 <- 0.75

df2 <- subset(df1, rho0 == this_rho0)
df2 <- subset(df2, rho == 0.7 | rho == 0.75 | rho == 0.8)
df2 <- subset(df2, a == 0.7)

mypal <- colorspace::rainbow_hcl(5)
mypal_use <- c("decc.q" = mypal[1],
               "ecc.q" = mypal[2],
               "ecc.s" = mypal[3],
               "gca" = mypal[4],
               "ssh" = mypal[5])

df2$model <- factor(df2$model, levels = c("decc.q", "ecc.q", "ecc.s", "gca", "ssh"))

df2 <- subset(df2, model != "ecc.q")

model_vec <- c("dECC", "ECC-S", "GCA", "SSh")

this_sigma <- input_sigma[2]

this_score <- "es_list"

dfplot <- subset(df2, sigma == this_sigma & score == this_score)
dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 

p1 <- ggplot(dfplot, aes(model, value, colour = model)) + 
  geom_rect(data = subset(dfplot, equalRhos == 1), color = "black", size = 2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf) 
p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25) 
p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
p1 <- p1 + facet_grid(rows = vars(a), cols = vars(rho), 
                      labeller = label_bquote(rows = a==.(a), 
                                              cols = rho==.(rho)))
sigval <- round(this_sigma, 2)
scval <- strsplit(this_score, split = "_")[[1]][1]
if(scval == "es"){
  p1 <- p1 + ggtitle(bquote(list("Energy Score", sigma == .(sigval), rho[0] ==.(this_rho0))))
} 

p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
p1 <- p1 + xlab("Model") + ylab("DM test statistic") 
p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", label = model_vec) 
p1 <- p1 + scale_x_discrete(label = model_vec)
p1save <- p1

## VS

this_score <- "vs1_list"

dfplot <- subset(df2, sigma == this_sigma & score == this_score)
dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 

p1 <- ggplot(dfplot, aes(model, value, colour = model)) + 
  geom_rect(data = subset(dfplot, equalRhos == 1), color = "black", size = 2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf) 
p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25) 
p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
p1 <- p1 + facet_grid(rows = vars(a), cols = vars(rho), 
                      labeller = label_bquote(rows = a==.(a), 
                                              cols = rho==.(rho)))
sigval <- round(this_sigma, 2)
scval <- strsplit(this_score, split = "_")[[1]][1]
if(scval == "vs1"){
  p1 <- p1 + ggtitle(bquote(list("Variogram Score", sigma == .(sigval), rho[0] ==.(this_rho0))))
} 

p1 <- p1 + theme_bw() + theme(legend.position = "none")
p1 <- p1 + xlab("Model") + ylab("DM test statistic") 
p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", label = model_vec) 
p1 <- p1 + scale_x_discrete(label = model_vec)
p2save <- p1

pdf(paste0(plot_folder, "S3B_highRho0_a07_sigma", round(this_sigma,2), ".pdf"),  width = 7.5, height = 5.5, pointsize = 12)
grid.arrange(p1save, p2save, ncol = 1)
dev.off()
