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


