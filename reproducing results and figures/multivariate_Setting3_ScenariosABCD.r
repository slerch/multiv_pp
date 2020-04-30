# -------------------------------------------------------------------------- #

# Matrix plots for all rho/rho0 combinations for the Scenarios A, B, C, D
# presented in paper and supplement
# Each figure has 2 subpanels, one for ES and one for VS



library(ggplot2)
library(gridExtra)

mypal <- colorspace::rainbow_hcl(5)

mypal_use <- c("decc.q" = mypal[1],
               "ecc.q" = "blue",
               "ecc.s" = mypal[3],
               "gca" = mypal[4],
               "ssh" = mypal[5],
               "ens" = "black",
               "emos.q" = mypal[2])
               

plot_folder <- "/path/to/save/figures/"

# Function to identify outliers in the data frame model-wise
# for deleting them for plotting purposes
source(file = "/path/to/Rfunctions/identify_outliers.r")

               
load("/path/to/save/output/df_gev0.Rdata")


df <- df100; rm(df100)


# Plots for fixed combination of GEV parameters in obs und ens

names(df)

# Scenarios A, B, C, D obtained by extracting results for the respective parameter combinations of GEV0

# Scenario A  (Supplement)
#df1 <- subset(df, gev0loc_obs == 0 & gev0shape_obs==-0.1 & gev0scale_obs==1 & gev0loc_ens==1 & gev0shape_ens==0 & gev0scale_ens==0.2)
# Scenario B (main paper)
df1 <- subset(df, gev0loc_obs == 0 & gev0shape_obs==-0.1 & gev0scale_obs==1 & gev0loc_ens==0 & gev0shape_ens==0 & gev0scale_ens==2)
# Scenario C (Supplement)
#df1 <- subset(df, gev0loc_obs == 1 & gev0shape_obs==0.3 & gev0scale_obs==1 & gev0loc_ens==0 & gev0shape_ens==0 & gev0scale_ens==2)
# Scenario D (Supplement)
#df1 <- subset(df, gev0loc_obs == 0 & gev0shape_obs==0 & gev0scale_obs==1 & gev0loc_ens==0 & gev0shape_ens==0 & gev0scale_ens==1)

str(df1)
df1_save <- df1

# correct sign of values to have positive values indicate improvements over ECC.Q as in SS plots
df1$value <- (-1)*df1$value

input_scores <- unique(df1$score)


# Remove reference model ECC-Q
df1 <- subset(df1, model != "ecc.q")
df1$model <- factor(df1$model, levels = c("ens", "emos.q", "decc.q", "ecc.s", "gca", "ssh"))
model_vec <- c("ENS", "EMOS-Q", "dECC", "ECC-S", "GCA", "SSh")



    #plotfname <- paste0(plot_folder, "SA_DMtest_ESVS_3x3.pdf")
    plotfname <- paste0(plot_folder, "SB_DMtest_ESVS_3x3.pdf")
    #plotfname <- paste0(plot_folder, "SC_DMtest_ESVS_3x3.pdf")
    #plotfname <- paste0(plot_folder, "SD_DMtest_ESVS_3x3.pdf")


# ES

this_score <- "es_list"
dfplot <- subset(df1, score == this_score)

dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho)


# For the GEV0 results, some extreme outliers can and do appear frequently.
# In order to have a better view of the main bulk of data points and
# not letting the comparison of boxplots suffer from single extreme outliers,
# a removal of outliers outside the standard fences (1.5 times IQR) is
# conducted before plotting
# For this, source the outlier identification function above and
# then apply the following
mod_out <- by(dfplot$value, dfplot$model, hf, coef=1.5)
mod_names <- levels(dfplot$model)
keep_ind <- NULL
for (i in 1:length(mod_names)){
  sel <- mod_names[i]
  add_keep <- which(dfplot$model == sel)[!mod_out[[sel]]]
  keep_ind <- c(keep_ind, add_keep)
}

dfplot <- dfplot[keep_ind,]


    p1 <- ggplot(dfplot, aes(model, value, colour = model))
    p1 <- p1 + geom_rect(data = subset(dfplot, equalRhos == TRUE), color = "black", size=2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf)
    p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25)
    p1 <- p1 + geom_boxplot() + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")

    p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
    p1 <- p1 + xlab("Model") + ylab("DM test statistic")
    p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", label = model_vec)
    p1 <- p1 + scale_x_discrete(label = model_vec)
    p1 <- p1 + facet_grid(rows = vars(rho), cols = vars(rho0),
                       labeller = label_bquote(rows = rho ==.(rho),
                                              cols = rho[0]==.(rho0)))

    # Specify m to be the number of ensemble members you are visualizing
    # The other parameter values of observation and forecast distributions are automatically read from the
    # data frame containing the DM Test results, see to have a data frame with the respective structure
    # or enter the parameter values you plot manually here
    p1 <- p1 + ggtitle(bquote(list("Energy Score", m == 50, mu[0] == .(dfplot$gev0loc_obs), sigma[0] == .(dfplot$gev0scale_obs),
                 xi[0] == .(dfplot$gev0shape_obs), mu == .(dfplot$gev0loc_ens),
                 sigma == .(dfplot$gev0scale_ens), xi == .(dfplot$gev0shape_ens))))


    p1_save <- p1




# VS

this_score <- "vs1_list"
dfplot <- subset(df1, score == this_score)

dfplot$equalRhos <- (dfplot$rho0 == dfplot$rho)


# For the GEV0 results, some extreme outliers can and do appear frequently.
# In order to have a better view of the main bulk of data points and
# not letting the comparison of boxplots suffer from single extreme outliers,
# a removal of outliers outside the standard fences (1.5 times IQR) is
# conducted before plotting
# For this, source the outlier identification function above and
# then apply the following
mod_out <- by(dfplot$value, dfplot$model, hf, coef=1.5)
mod_names <- levels(dfplot$model)
keep_ind <- NULL
for (i in 1:length(mod_names)){
  sel <- mod_names[i]
  add_keep <- which(dfplot$model == sel)[!mod_out[[sel]]]
  keep_ind <- c(keep_ind, add_keep)
}

dfplot <- dfplot[keep_ind,]


    p1 <- ggplot(dfplot, aes(model, value, colour = model))
    p1 <- p1 + geom_rect(data = subset(dfplot, equalRhos == TRUE), color = "black", size=2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf)
    p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25)
    p1 <- p1 + geom_boxplot() + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")

    p1 <- p1 + theme_bw() + theme(legend.position = "none")
    p1 <- p1 + xlab("Model") + ylab("DM test statistic")
    p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", label = model_vec)
    p1 <- p1 + scale_x_discrete(label = model_vec)
    p1 <- p1 + facet_grid(rows = vars(rho), cols = vars(rho0),
                      labeller = label_bquote(rows = rho ==.(rho),
                                              cols = rho[0]==.(rho0)))


    # Specify m to be the number of ensemble members you are visualizing
    # The other parameter values of observation and forecast distributions are automatically read from the
    # data frame containing the DM Test results, see to have a data frame with the respective structure
    # or enter the parameter values you plot manually here
    p1 <- p1 + ggtitle(bquote(list("Variogram Score", m == 50, mu[0] == .(dfplot$gev0loc_obs), sigma[0] == .(dfplot$gev0scale_obs),
              xi[0] == .(dfplot$gev0shape_obs), mu == .(dfplot$gev0loc_ens),
              sigma == .(dfplot$gev0scale_ens), xi == .(dfplot$gev0shape_ens))))



    p2_save <- p1



pdf(plotfname, width = 10, height = 13, pointsize = 11)
                  grid.arrange(p1_save, p2_save, ncol = 1)
dev.off()





