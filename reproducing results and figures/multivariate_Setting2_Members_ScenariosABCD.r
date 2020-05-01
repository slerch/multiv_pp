# -------------------------------------------------------------------------- #

# Matrix plots for all rho/rho0 combinations for the Scenarios A, B, C, D
# presented in paper and supplement
# Within each rho-rho0 panel the number of Members are represented
# in the plot by grouped boxplots for
# each model, with varying shading of the respective model colour
# This plots are separately proudced for ES and VS in order not to have
# too much information in a single plot




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


# In each panel boxplots different numbers of members set side by side, for each of the settings A,B,C,D
# Here example for plotting 5, 20, 50, 100 members

# For the automated plotting, you need to define a large data frame containing results of all number of members you wish to plot
# Load each data frame, define number of members in that data frame as new colum, rbind all of them together,
# define member numbers as factor variable

load("/path/to/save/output/df_gev0_5Members.Rdata")
dfm5 <- df100
dfm5$m <- 5

load("/path/to/save/output/df_gev0_20Members.Rdata")
dfm20 <- df100
dfm20$m <- 20

load("/path/to/save/output/df_gev0_50Members.Rdata")
dfm50 <- df100
dfm50$m <- 50

load("/path/to/save/output/df_gev0_100Members.Rdata")
dfm100 <- df100
dfm100$m <- 100

df <- rbind(dfm5, dfm20, dfm50, dfm100)
df$m <- factor(df$m)


# Plots for fixed combination of GEV parameters in obs und ens

names(df)

# Scenarios A, B, C, D obtained by extracting results for the respective parameter combinations of GEV0

# Scenario A  (Supplement)
#df1 <- subset(df, gev0loc_obs == 0 & gev0shape_obs==-0.1 & gev0scale_obs==1 & gev0loc_ens==1 & gev0shape_ens==0 & gev0scale_ens==0.2)
# Scenario B (Supplement)
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


# As the plot contains a lot of different layers, it may be necessary to save
# as png instead of pdf

    #plotfname <- paste0(plot_folder, "SA_DMtest_ES_Members_3x3.pdf")
    plotfname <- paste0(plot_folder, "SB_DMtest_ES_Members_3x3.pdf")
    #plotfname <- paste0(plot_folder, "SC_DMtest_ES_Members_3x3.pdf")
    #plotfname <- paste0(plot_folder, "SD_DMtest_ES_Members_3x3.pdf")

# The alternative plot names in case you wish to save as png
    #plotfname <- paste0(plot_folder, "SA_DMtest_ES_Members_3x3.png")
    #plotfname <- paste0(plot_folder, "SB_DMtest_ES_Members_3x3.png")
    #plotfname <- paste0(plot_folder, "SC_DMtest_ES_Members_3x3.png")
    #plotfname <- paste0(plot_folder, "SD_DMtest_ES_Members_3x3.png")


# Exemplarily for ES, for plotting VS select the score here

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



    p1 <- ggplot(dfplot, aes(model, value, colour = model, fill=model, alpha=m))
    p1 <- p1 + geom_rect(data = subset(dfplot, equalRhos == TRUE), color = "black", size=2, fill = NA, xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf)
    p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75", alpha=0.25)
    p1 <- p1 + geom_boxplot() + guides(color = guide_legend(override.aes = list(fill = NA)), alpha="none") + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")

    p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
    p1 <- p1 + xlab("Model") + ylab("DM test statistic")
    p1 <- p1 + scale_x_discrete(label = model_vec)
    p1 <- p1 + scale_alpha_discrete("m", range = c(0.1, 0.7))
    p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", labels = model_vec)
    p1 <- p1 + scale_fill_manual(values = mypal_use, name = "Model", labels = model_vec)

    p1 <- p1 + facet_grid(rows = vars(rho), cols = vars(rho0),
                       labeller = label_bquote(rows = rho ==.(rho),
                                              cols = rho[0]==.(rho0)))


    p1 <- p1 + ggtitle(bquote(list("Energy Score", mu[0] == .(dfplot$gev0loc_obs), sigma[0] == .(dfplot$gev0scale_obs),
  xi[0] == .(dfplot$gev0shape_obs), mu == .(dfplot$gev0loc_ens),
  sigma == .(dfplot$gev0scale_ens), xi == .(dfplot$gev0shape_ens))))


    p1_save <- p1


 pdf(plotfname, width = 15, height = 10, pointsize = 11)
 p1_save
 dev.off()
 
 
 # In case you save as png
 #png(plotfname, width = 1000, height = 800, pointsize = 11)
 #p1_save
 #dev.off()
