####
# This code is generating a specific spectral disease index for the pathosystem
# Austropuccinia psidii and Backhousia citriodora. Therefore, it utilizes a 
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


# install.packages(c("cowplot", "gdata", "glmulti", "hsdar", "plyr", 
#                    "PresenceAbsence", "prospectr", "rJava", "tidyverse", 
#                    "VSURF"))

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

source('functions/raw2speclib_hsdar.R')
source('functions/DropCatVar_July2017.R')
source('functions/ModelSelect.R')
source('functions/LMMRindex.R')
source('functions/LMMRloop.R')
source('functions/index2prob_test.R')


# Set Project Structure --------------------------------------------------------


dir.create('data', FALSE, FALSE)
dir.create('data/raw', FALSE, FALSE)
dir.create('data/temp', FALSE, FALSE)
dir.create('output', FALSE, FALSE)
dir.create('output/tables', FALSE, FALSE)
dir.create('output/figures', FALSE, FALSE)
dir.create('manuscript', FALSE, FALSE)
dir.create('manuscript/templates', FALSE, FALSE)
dir.create('analysis', FALSE, FALSE)
dir.create('functions', FALSE, FALSE)


# Loading and Preparing Data ---------------------------------------------------


ori.data <- read.csv('data/raw/data.wo.out.binned.cut.csv')

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

# saveRDS(feature.set,'data/temp/featuresforindex.rds')
feature.set <- readRDS('data/temp/featuresforindex.rds') 

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

result_df <- rename(result_df, LMMR = index_val)

modelcoeffs <- coef(model.1) # Extracting model coefficients for LMMR design

saveRDS(modelcoeffs, 'data/temp/coefficients.rds') # Export for future use
saveRDS(best.bands, 'data/temp/bestbands.rds') # Export for future use


# Attach Indices to Compare to LMMR  -------------------------------------------


tospectra <-
    read.csv("data/raw/data.wo.out.binned.cut.csv", check.names = FALSE)

tospectra <- DropClass(tospectra, tospectra$Type, "Healthy")

spectra <- raw2speclib(tospectra) # Use hsdar to build spectral library

# Define spectral vegetation indices to use them in hsdar pkg

#ARI <-  '((R550)^−1) − ((R700)^−1)' #Anthocyanin Reflectance Index
#RVSI <- '(((R712) + (R752))/2) − (R732)' #Red_Edge Vegetation Stress Index
NBNDVI <- '(R850-R680)/(R850+R680)'
#SIPI <- '(R800-R445)/(R800+R680)'

# ?vegindex OPTIONAL: Call indices known by the hsdar pkg

ind <- c("PRI", "MCARI", NBNDVI) # add more indices HERE (some known by hsdar)


for (i in ind) {
    result_df[[paste(i, 'prob', sep = "_")]] <-
        Index2Prob(result_df, spectra, i)
    
} # HERE IS A PROBLEM!!!

ind.c <- c('PRI', 'MCARI', 'NBNDVI') #manipulate if necessary

colnames(result_df)[3:length(result_df)] <- ind.c

write_csv(result_df,'output/tables/result_df.csv')

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
        ylab(label = 'Disease Prob.')
}



p1 <- plot_grid(
    plot_list[[1]],
    plot_list[[2]],
    plot_list[[3]],
    plot_list[[4]],
    labels = c("a", 'b', 'c', 'd'),
    ncol = 2,
    nrow = 2
)
p1
ggsave(
    "output/figures/LMMR_compare_v2.pdf",
    plot = p1,
    width = 40,
    height = 20,
    units = "cm",
    dpi = 400
)
####
# 8. Calculate Accuracy Metrics
####

# turning treated and untreated into numbers i.e. 0 and 1

result_df$Type <-
    as.numeric(result_df$Type == 'Untreated') # 1 = treated; 0 = untreated

myDat <- cbind(ID = seq(1, dim(result_df)[1]), result_df)

myCMX <- cmx(myDat)

acc <- presence.absence.accuracy(myDat)

acc_fin <- acc[, c(1, 3, 4, 5, 6, 7)]

saveRDS(acc_fin, 'data/temp/Table1.rds')

