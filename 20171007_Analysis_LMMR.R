####
# This code is generating a specific spectral disease index for the pathosystem
# Austropuccinia psidii and Backhousia citriodora. It utilizes a 
# spectral dataset containing following columns:
#
# Type= Categorical variables referring to spectral class to be classified.
# Wavelength: Contains class IDs -> Healthy = Australian Botanical Garden Mount 
# Annan, Treated = Fungicide treated plants plantation, Untreated = Untreated 
# plants plantation and 350,351,....,2500: Each column contains values of 
# spectral reflectance at the specified wavelength in [%].
#
# The full manuscript was published in: ...
####


# Loading Packages and Functions------------------------------------------------


install.packages(c("cowplot", "gdata", "glmulti", "hsdar", "plyr",
                   "PresenceAbsence", "prospectr", "rJava", "tidyverse",
                   "VSURF"))

#Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_151') 

library(cowplot)
library(gdata)
library(glmulti)
library(hsdar)
library(plyr)
library(PresenceAbsence)
library(prospectr)
library(rJava)
library(tidyverse)
library(VSURF)
library(reshape2)
library(gridExtra)

source('R/raw2speclib_hsdar.R')
source('R/DropCatVar_July2017.R')
source('R/ModelSelect.R')
source('R/LMMRindex.R')
source('R/LMMRloop.R')
source('R/index2prob_test.R')


# Set Project Structure --------------------------------------------------------

dir.create('data', FALSE, FALSE)
dir.create('R', FALSE, FALSE)
dir.create('doc', FALSE, FALSE)
dir.create('output', FALSE, FALSE)

# Loading and Preparing Data ---------------------------------------------------

ori.data <- read.csv('data/data.wo.out.binned.cut.csv')

levels(ori.data$Type)

data.log <- DropClass(ori.data, ori.data$Type, 'Healthy') # Drop healthy leaves

data.log$Type <-
    as.numeric(data.log$Type == 'Untreated') # 1 = treated; 0 = untreated

data.log[names(data.log)[-1]] <-
    log(data.log[names(data.log)[-1]]) # logs frequencies


# Selecting Subset of Relevant Wavebands ---------------------------------------

    # First step using VSURF

#feature.set <- 
    #VSURF(data.log[,2:2152], data.log[,1], clusterType = "FORK", ntree = 2000, 
          #mtry = 50) Takes a while, therefore saved/loaded as .rds

# saveRDS(feature.set,'output')
feature.set <- readRDS('output/featuresforindex.rds') 

    # Second step using glmulti

model.1 <- GLM.featureselect(data.log, feature.set)

best.bands <-
    row.names(summary(model.1)$coefficients)[c(2, 3, 4, 5)]

df.bands <-
    cbind('Type' = ori.data$Type, ori.data[best.bands])

df.bands <- DropClass(df.bands, df.bands$Type, "Healthy")

index_vals <-
    index4col(df.bands[, 2:5], model.1) # Function that generates LMMR values

result_df <-
    cbind('Type' = df.bands$Type, index_vals) # Attach LMMR values to dataframe

names(result_df)[2] <- c('LMMR') #Rename col name

modelcoeffs <- coef(model.1) # Extracting model coefficients for LMMR design

#saveRDS(modelcoeffs, 'output/coefficients.rds') # Export for future use
#saveRDS(best.bands, 'output/bestbands.rds') # Export for future use


# Attach Indices to Compare to LMMR  -------------------------------------------


tospectra <- read.csv("data/data.wo.out.binned.cut.csv", check.names = FALSE)

tospectra <- DropClass(tospectra, tospectra$Type, "Healthy")

spectra <- raw2speclib(tospectra) # Use hsdar to build spectral library

ID4spectra <- as.numeric(names(tospectra[2:202]))

spectra@ID

# Define spectral vegetation indices to use them in hsdar pkg

#ARI <-  '((R550)^−1) − ((R700)^−1)' #Anthocyanin Reflectance Index
#RVSI <- '(((R712) + (R752))/2) − (R732)' #Red_Edge Vegetation Stress Index
NBNDVI <- '(R850-R680)/(R850+R680)'
#SIPI <- '(R800-R445)/(R800+R680)'
LMMR2 <- '((R545/R555)^(5/3))*(R1485/R2175)'

# ?vegindex OPTIONAL: Call indices known by the hsdar pkg

ind <- c("PRI", "MCARI", NBNDVI, LMMR2) # add more indices HERE (some known by hsdar)


for (i in ind) {
    result_df[[paste(i, 'prob', sep = "_")]] <-
        Index2Prob(result_df, spectra, i)
    
} 

ind.c <- c('PRI', 'MCARI', 'NBNDVI', 'LMMR2') #manipulate if necessary

colnames(result_df)[3:length(result_df)] <- ind.c

write_csv(result_df,'output/result_df.csv')

# Visualize Results ------------------------------------------------------------

plot_list <- list()
indi <- names(result_df[, 2:length(result_df)])

for (i in indi) {
    plot_list[[i]] <-
        ggplot(result_df,
               aes_string(result_df$Type, result_df[, i], colour = result_df$Type)) +
        geom_jitter(height = 0) +
        geom_boxplot(colour = 'black',
                     alpha = 0.5,
                     outlier.alpha = 0) +
        theme(
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 16, face = "bold"),
            legend.position = "none"
        ) +
        xlab(label = paste(i)) +
        ylab(label = 'Index Value')+
        ylim(0, 1)
}



p1 <- plot_grid(
    plot_list[[2]],
    plot_list[[3]],
    plot_list[[4]],
    plot_list[[5]],
    labels = c("a", 'b', 'c', 'd'),
    ncol = 4,
    nrow = 1
)
p1

ggsave(
    "output/LMMR_compare_v2.pdf",
    plot = p1,
    width = 40,
    height = 20,
    units = "cm",
    dpi = 400
)

# Spectra plots

source('R/prepgg_December2017.R')

spectra.gg <- prep.gg(tospectra)

pspec <- ggplot(spectra.gg, aes(Wavelength, Reflectance, colour = Type)) +
  geom_line(size = .5)+
  annotate(
    "rect",
    xmin = 500,
    xmax = 570,
    ymin = -Inf,
    ymax = Inf,
    alpha = .2,
    fill = 'green'
  ) +
  theme_set(theme_bw(base_size = 20))+
  theme(legend.position = c(.90, .88), legend.title = element_blank(), legend.background = element_blank())

p5 <- plot_list[[2]]
p5 <- p5+
  theme(#axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())
p6 <- plot_list[[3]]
p6 <- p6+
  theme(#axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank())
p7 <- plot_list[[4]]
p7 <- p7+
  theme(#axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank())
p8 <- plot_list[[5]]
p8 <- p8+
  theme(#axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank())+
    labs(x = "LMMR")

plot.res <- ggdraw() +
  draw_plot(p5, x = 0, y = .5, width = .25, height = .5) +
  draw_plot(p6, x = .25, y = .5, width = .25, height = .5) +
  draw_plot(p7, x = .5, y = .5, width = .25, height = .5) +
  draw_plot(p8, x = .75, y = .5, width = .25, height = .5) +
  draw_plot(pspec, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 12,
                  x = c(0, 0.25, 0.5, 0.75, 0), y = c(1, 1, 1, 1, 0.5))

ggsave(
  "output/Results.pdf",
  plot = plot.res,
  width = 40,
  height = 20,
  units = "cm",
  dpi = 400
)

# Calculate Accuracy Metrics ---------------------------------------------------


result_df$Type <-
    as.numeric(result_df$Type == 'Untreated') # 1 = treated; 0 = untreated

myDat <- cbind(ID = seq(1, dim(result_df)[1]), result_df)

myCMX <- cmx(myDat) #Select columns here to get confusion matrix

acc <- presence.absence.accuracy(myDat)

acc_fin <- acc[, c(1, 3, 4, 5, 6, 7)]

saveRDS(acc_fin, 'data/Table1.rds')

