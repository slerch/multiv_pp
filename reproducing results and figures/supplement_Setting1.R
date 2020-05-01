rm(list=ls())

library(ggplot2)
library(gridExtra)

load("/path/to/save/output/df_mvNormal.Rdata")

df <- df100; rm(df100)
input_eps <- unique(df$eps)
input_rho0 <- c(0.1, 0.25, 0.5, 0.75, 0.9)
input_sigma <-  unique(df$sigma)
input_rho <- c(0.1, 0.25, 0.5, 0.75, 0.9)
input_scores <- unique(df$score)

df$value <- (-1)*df$value

plot_folder <- "/path/to/save/figures/"

df1 <- subset(df, model != "ens" & model != "emos.q")

mypal <- colorspace::rainbow_hcl(5)
mypal_use <- c("decc.q" = mypal[1],
               "ecc.q" = mypal[2],
               "ecc.s" = mypal[3],
               "gca" = mypal[4],
               "ssh" = mypal[5])

df1$model <- factor(df1$model, levels = c("decc.q", "ecc.q", "ecc.s", "gca", "ssh"))
model_vec <- c("dECC", "ECC-Q", "ECC-S", "GCA", "SSh")

library(gridExtra)

for(this_eps in input_eps){
  print(this_eps)
  for(this_sigma in input_sigma){
    cat("...", this_sigma, "\n")
    
    # top part = ES
    this_score <- "es_list"
    
    dfplot <- subset(df1, eps == this_eps & sigma == this_sigma & score == this_score)
    dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 
    
    p1 <- ggplot(dfplot, aes(model, value, colour = model)) + 
      geom_rect(data = subset(dfplot, equalRhos == 1), color = "black", size = 2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf) 
    p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25) 
    p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
    p1 <- p1 + facet_grid(rows = vars(rho), cols = vars(rho0), 
                          labeller = label_bquote(rows = rho==.(rho), 
                                                  cols = rho[0]==.(rho0)))    
    sigval <- round(this_sigma, 2)
    scval <- strsplit(this_score, split = "_")[[1]][1]
    if(scval == "es"){
      p1 <- p1 + ggtitle(bquote(list("Energy Score", epsilon == .(this_eps), sigma == .(sigval))))
    } 
    if(scval == "vs1"){
      p1 <- p1 + ggtitle(bquote(list("Variogram Score", epsilon == .(this_eps), sigma == .(sigval))))
    }
    p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
    p1 <- p1 + xlab("Model") + ylab("DM test statistic") 
    p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", labels = model_vec) 
    p1 <- p1 + scale_x_discrete(label = model_vec)
    p1save <- p1
    
    # bottom part = VS
    this_score <- "vs1_list"
    
    dfplot <- subset(df1, eps == this_eps & sigma == this_sigma & score == this_score)
    dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 
    
    p1 <- ggplot(dfplot, aes(model, value, colour = model)) + 
      geom_rect(data = subset(dfplot, equalRhos == 1), color = "black", size = 2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf) 
    p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25) 
    p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
    p1 <- p1 + facet_grid(rows = vars(rho), cols = vars(rho0), 
                          labeller = label_bquote(rows = rho==.(rho), 
                                                  cols = rho[0]==.(rho0)))    
    sigval <- round(this_sigma, 2)
    scval <- strsplit(this_score, split = "_")[[1]][1]
    if(scval == "es"){
      p1 <- p1 + ggtitle(bquote(list("Energy Score", epsilon == .(this_eps), sigma == .(sigval))))
    } 
    if(scval == "vs1"){
      p1 <- p1 + ggtitle(bquote(list("Variogram Score", epsilon == .(this_eps), sigma == .(sigval))))
    }
    p1 <- p1 + theme_bw() + theme(legend.position = "none") # no legend here
    p1 <- p1 + xlab("Model") + ylab("DM test statistic") 
    p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", labels = model_vec) 
    p1 <- p1 + scale_x_discrete(label = model_vec)
    p2save <- p1
    
    plotname <- paste0("supplement_S1_eps", this_eps, "_sigma", 100*round(this_sigma,2), ".pdf")
    pdf(paste0(plot_folder, plotname), width = 15, height = 24, pointsize = 11)
    grid.arrange(p1save, p2save, ncol = 1)
    dev.off()
    
  }
}

##
##
##

## During the revision process, additional simulations with varying numbers of members and dimensions were run, Code to reproduce corresponding plots from the supplemental material is provided in the following:

## To obtain the corresponding Rdata files, the simulation code needs to be run with the corresponding choices of m and d

## 
## additional numbers of dimensions, d = 2,3,4,10,20,30,50
##

load("/path/to/Rdata/")

df <- df100; rm(df100)

names(df)
df1 <- subset(df, eps == 1)
input_rho0 <- unique(df1$rho0)
input_sigma <-  unique(df1$sigma)
input_rho <- unique(df1$rho)

df1$value <- (-1)*df1$value

df0 <- df1

input_scores <- unique(df1$score)
plot_folder <- "..."


df1_save <- df0 
df1 <- subset(df1_save, model != "ens")
df1 <- subset(df1, d == this_d)

df2 <- subset(df1, model != "emos.q")

mypal <- colorspace::rainbow_hcl(5)
mypal_use <- c("decc.q" = mypal[1],
               "ecc.q" = mypal[2],
               "ecc.s" = mypal[3],
               "gca" = mypal[4],
               "ssh" = mypal[5])

df2$model <- factor(df2$model, levels = c("decc.q", "ecc.q", "ecc.s", "gca", "ssh"))

df2 <- subset(df2, model != "ecc.q")
df2 <- subset(df2, d == 2 | d == 3 | d == 4 | d == 10 | d == 20 | d == 30 | d == 50)
df2$d <- factor(df2$d)

model_vec <- c("dECC", "ECC-S", "GCA", "SSh")

## ES

this_score <- "es_list"

this_sigma <- input_sigma[1]

dfplot <- subset(df2, sigma == this_sigma & score == this_score)
dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 

p1 <- ggplot(dfplot, aes(model, value, colour = model, fill=model, alpha=d))
p1 <- p1 + geom_rect(data = subset(dfplot, equalRhos == TRUE), color = "black", size=2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf)
p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25)
p1 <- p1 + geom_boxplot(outlier.shape = NA) + guides(color = guide_legend(override.aes = list(fill = NA)), alpha="none") + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
p1 <- p1 + facet_grid(rows = vars(rho), cols = vars(rho0), 
                      labeller = label_bquote(rows = rho==.(rho), 
                                              cols = rho[0]==.(rho0)))
p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
p1 <- p1 + xlab("Model") + ylab("DM test statistic")
p1 <- p1 + scale_x_discrete(label = model_vec)
p1 <- p1 + scale_alpha_discrete("d", range = c(0.1, 0.7))
p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", labels = model_vec)
p1 <- p1 + scale_fill_manual(values = mypal_use, name = "Model", labels = model_vec)
sigval <- round(this_sigma, 2)
scval <- strsplit(this_score, split = "_")[[1]][1]
if(scval == "es"){
  p1 <- p1 + ggtitle(bquote(list("Energy Score", sigma == .(sigval))))
} 
if(scval == "vs1"){
  p1 <- p1 + ggtitle(bquote(list("Variogram Score", sigma == .(sigval))))
}
p1save <- p1

# next sigma 

this_sigma <- input_sigma[3]

dfplot <- subset(df2, sigma == this_sigma & score == this_score)
dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 

p1 <- ggplot(dfplot, aes(model, value, colour = model, fill=model, alpha=d))
p1 <- p1 + geom_rect(data = subset(dfplot, equalRhos == TRUE), color = "black", size=2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf)
p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25)
p1 <- p1 + geom_boxplot(outlier.shape = NA) + guides(color = guide_legend(override.aes = list(fill = NA)), alpha="none") + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
p1 <- p1 + facet_grid(rows = vars(rho), cols = vars(rho0), 
                      labeller = label_bquote(rows = rho==.(rho), 
                                              cols = rho[0]==.(rho0)))
p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
p1 <- p1 + xlab("Model") + ylab("DM test statistic")
p1 <- p1 + scale_x_discrete(label = model_vec)
p1 <- p1 + scale_alpha_discrete("d", range = c(0.1, 0.7))
p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", labels = model_vec)
p1 <- p1 + scale_fill_manual(values = mypal_use, name = "Model", labels = model_vec)
sigval <- round(this_sigma, 2)
scval <- strsplit(this_score, split = "_")[[1]][1]
if(scval == "es"){
  p1 <- p1 + ggtitle(bquote(list("Energy Score", sigma == .(sigval))))
} 
if(scval == "vs1"){
  p1 <- p1 + ggtitle(bquote(list("Variogram Score", sigma == .(sigval))))
}
p2save <- p1

pdf(paste0(plot_folder, "ES_S1rev_dim", ".pdf"), width = 7.5, height = 11, pointsize = 12)
grid.arrange(p1save, p2save, ncol = 1)
dev.off()

## VS

this_score <- "vs1_list"

this_sigma <- input_sigma[1]

dfplot <- subset(df2, sigma == this_sigma & score == this_score)
dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 

p1 <- ggplot(dfplot, aes(model, value, colour = model, fill=model, alpha=d))
p1 <- p1 + geom_rect(data = subset(dfplot, equalRhos == TRUE), color = "black", size=2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf)
p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25)
p1 <- p1 + geom_boxplot(outlier.shape = NA) + guides(color = guide_legend(override.aes = list(fill = NA)), alpha="none") + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
p1 <- p1 + facet_grid(rows = vars(rho), cols = vars(rho0), 
                      labeller = label_bquote(rows = rho==.(rho), 
                                              cols = rho[0]==.(rho0)))
p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
p1 <- p1 + xlab("Model") + ylab("DM test statistic")
p1 <- p1 + scale_x_discrete(label = model_vec)
p1 <- p1 + scale_alpha_discrete("d", range = c(0.1, 0.7))
p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", labels = model_vec)
p1 <- p1 + scale_fill_manual(values = mypal_use, name = "Model", labels = model_vec)
sigval <- round(this_sigma, 2)
scval <- strsplit(this_score, split = "_")[[1]][1]
if(scval == "es"){
  p1 <- p1 + ggtitle(bquote(list("Energy Score", sigma == .(sigval))))
} 
if(scval == "vs1"){
  p1 <- p1 + ggtitle(bquote(list("Variogram Score", sigma == .(sigval))))
}
p1save <- p1

# next sigma 

this_sigma <- input_sigma[3]

dfplot <- subset(df2, sigma == this_sigma & score == this_score)
dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 

p1 <- ggplot(dfplot, aes(model, value, colour = model, fill=model, alpha=d))
p1 <- p1 + geom_rect(data = subset(dfplot, equalRhos == TRUE), color = "black", size=2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf)
p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25)
p1 <- p1 + geom_boxplot(outlier.shape = NA) + guides(color = guide_legend(override.aes = list(fill = NA)), alpha="none") + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
p1 <- p1 + facet_grid(rows = vars(rho), cols = vars(rho0), 
                      labeller = label_bquote(rows = rho==.(rho), 
                                              cols = rho[0]==.(rho0)))
p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
p1 <- p1 + xlab("Model") + ylab("DM test statistic")
p1 <- p1 + scale_x_discrete(label = model_vec)
p1 <- p1 + scale_alpha_discrete("d", range = c(0.1, 0.7))
p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", labels = model_vec)
p1 <- p1 + scale_fill_manual(values = mypal_use, name = "Model", labels = model_vec)
sigval <- round(this_sigma, 2)
scval <- strsplit(this_score, split = "_")[[1]][1]
if(scval == "es"){
  p1 <- p1 + ggtitle(bquote(list("Energy Score", sigma == .(sigval))))
} 
if(scval == "vs1"){
  p1 <- p1 + ggtitle(bquote(list("Variogram Score", sigma == .(sigval))))
}
p2save <- p1

pdf(paste0(plot_folder, "VS_S1rev_dim", ".pdf"), width = 7.5, height = 11, pointsize = 12)
grid.arrange(p1save, p2save, ncol = 1)
dev.off()

## 
## additional numbers of members, m = 5, 10, 20, 35, 50, 100
##

load("/path/to/Rdata/")

df <- df100; rm(df100)

names(df)
df1 <- subset(df, eps == 1)
input_rho0 <- unique(df1$rho0)
input_sigma <-  unique(df1$sigma)
input_rho <- unique(df1$rho)

df1$value <- (-1)*df1$value

df0 <- df1

input_scores <- unique(df1$score)
plot_folder <- "..."

df1_save <- df0 
df1 <- subset(df1_save, model != "ens")
df2 <- subset(df1, model != "emos.q")

mypal <- colorspace::rainbow_hcl(5)
mypal_use <- c("decc.q" = mypal[1],
               "ecc.q" = mypal[2],
               "ecc.s" = mypal[3],
               "gca" = mypal[4],
               "ssh" = mypal[5])

df2$model <- factor(df2$model, levels = c("decc.q", "ecc.q", "ecc.s", "gca", "ssh"))
df2 <- subset(df2, model != "ecc.q")
model_vec <- c("dECC", "ECC-S", "GCA", "SSh")
df2$m <- factor(df2$m)

## ES

this_score <- "es_list"

this_sigma <- input_sigma[1]

dfplot <- subset(df2, sigma == this_sigma & score == this_score)
dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 

p1 <- ggplot(dfplot, aes(model, value, colour = model, fill=model, alpha=m))
p1 <- p1 + geom_rect(data = subset(dfplot, equalRhos == TRUE), color = "black", size=2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf)
p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25)
p1 <- p1 + geom_boxplot(outlier.shape = NA) + guides(color = guide_legend(override.aes = list(fill = NA)), alpha="none") + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
p1 <- p1 + facet_grid(rows = vars(rho), cols = vars(rho0), 
                      labeller = label_bquote(rows = rho==.(rho), 
                                              cols = rho[0]==.(rho0)))
p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
p1 <- p1 + xlab("Model") + ylab("DM test statistic")
p1 <- p1 + scale_x_discrete(label = model_vec)
p1 <- p1 + scale_alpha_discrete("m", range = c(0.1, 0.7))
p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", labels = model_vec)
p1 <- p1 + scale_fill_manual(values = mypal_use, name = "Model", labels = model_vec)
sigval <- round(this_sigma, 2)
scval <- strsplit(this_score, split = "_")[[1]][1]
if(scval == "es"){
  p1 <- p1 + ggtitle(bquote(list("Energy Score", sigma == .(sigval))))
} 
if(scval == "vs1"){
  p1 <- p1 + ggtitle(bquote(list("Variogram Score", sigma == .(sigval))))
}
p1save <- p1

# next sigma 

this_sigma <- input_sigma[4]

dfplot <- subset(df2, sigma == this_sigma & score == this_score)
dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 

p1 <- ggplot(dfplot, aes(model, value, colour = model, fill=model, alpha=m))
p1 <- p1 + geom_rect(data = subset(dfplot, equalRhos == TRUE), color = "black", size=2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf)
p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25)
p1 <- p1 + geom_boxplot(outlier.shape = NA) + guides(color = guide_legend(override.aes = list(fill = NA)), alpha="none") + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
p1 <- p1 + facet_grid(rows = vars(rho), cols = vars(rho0), 
                      labeller = label_bquote(rows = rho==.(rho), 
                                              cols = rho[0]==.(rho0)))
p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
p1 <- p1 + xlab("Model") + ylab("DM test statistic")
p1 <- p1 + scale_x_discrete(label = model_vec)
p1 <- p1 + scale_alpha_discrete("m", range = c(0.1, 0.7))
p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", labels = model_vec)
p1 <- p1 + scale_fill_manual(values = mypal_use, name = "Model", labels = model_vec)
sigval <- round(this_sigma, 2)
scval <- strsplit(this_score, split = "_")[[1]][1]
if(scval == "es"){
  p1 <- p1 + ggtitle(bquote(list("Energy Score", sigma == .(sigval))))
} 
if(scval == "vs1"){
  p1 <- p1 + ggtitle(bquote(list("Variogram Score", sigma == .(sigval))))
}
p2save <- p1

library(gridExtra)

pdf(paste0(plot_folder, "ES_S1rev_members", ".pdf"), width = 7.5, height = 11, pointsize = 12)
grid.arrange(p1save, p2save, ncol = 1)
dev.off()


## VS

this_score <- "vs1_list"

this_sigma <- input_sigma[1]

dfplot <- subset(df2, sigma == this_sigma & score == this_score)
dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 

p1 <- ggplot(dfplot, aes(model, value, colour = model, fill=model, alpha=m))
p1 <- p1 + geom_rect(data = subset(dfplot, equalRhos == TRUE), color = "black", size=2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf)
p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25)
p1 <- p1 + geom_boxplot(outlier.shape = NA) + guides(color = guide_legend(override.aes = list(fill = NA)), alpha="none") + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
p1 <- p1 + facet_grid(rows = vars(rho), cols = vars(rho0), 
                      labeller = label_bquote(rows = rho==.(rho), 
                                              cols = rho[0]==.(rho0)))
p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
p1 <- p1 + xlab("Model") + ylab("DM test statistic")
p1 <- p1 + scale_x_discrete(label = model_vec)
p1 <- p1 + scale_alpha_discrete("m", range = c(0.1, 0.7))
p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", labels = model_vec)
p1 <- p1 + scale_fill_manual(values = mypal_use, name = "Model", labels = model_vec)
sigval <- round(this_sigma, 2)
scval <- strsplit(this_score, split = "_")[[1]][1]
if(scval == "es"){
  p1 <- p1 + ggtitle(bquote(list("Energy Score", sigma == .(sigval))))
} 
if(scval == "vs1"){
  p1 <- p1 + ggtitle(bquote(list("Variogram Score", sigma == .(sigval))))
}
p1save <- p1

# next sigma

this_sigma <- input_sigma[4]

dfplot <- subset(df2, sigma == this_sigma & score == this_score)
dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho) 

p1 <- ggplot(dfplot, aes(model, value, colour = model, fill=model, alpha=m))
p1 <- p1 + geom_rect(data = subset(dfplot, equalRhos == TRUE), color = "black", size=2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf)
p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25)
p1 <- p1 + geom_boxplot(outlier.shape = NA) + guides(color = guide_legend(override.aes = list(fill = NA)), alpha="none") + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
p1 <- p1 + facet_grid(rows = vars(rho), cols = vars(rho0), 
                      labeller = label_bquote(rows = rho==.(rho), 
                                              cols = rho[0]==.(rho0)))
p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
p1 <- p1 + xlab("Model") + ylab("DM test statistic")
p1 <- p1 + scale_x_discrete(label = model_vec)
p1 <- p1 + scale_alpha_discrete("m", range = c(0.1, 0.7))
p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", labels = model_vec)
p1 <- p1 + scale_fill_manual(values = mypal_use, name = "Model", labels = model_vec)
sigval <- round(this_sigma, 2)
scval <- strsplit(this_score, split = "_")[[1]][1]
if(scval == "es"){
  p1 <- p1 + ggtitle(bquote(list("Energy Score", sigma == .(sigval))))
} 
if(scval == "vs1"){
  p1 <- p1 + ggtitle(bquote(list("Variogram Score", sigma == .(sigval))))
}
p2save <- p1

library(gridExtra)

pdf(paste0(plot_folder, "VS1_S1rev_members", ".pdf"), width = 7.5, height = 11, pointsize = 12)
grid.arrange(p1save, p2save, ncol = 1)
dev.off()