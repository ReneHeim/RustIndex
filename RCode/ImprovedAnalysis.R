####
# 1. Loading packages and functions
####

library(hsdar)
library(gdata)
library(plyr)
library(rJava)
library(glmulti)
library(tidyverse)
library(PresenceAbsence)
library(cowplot)
library(VSURF)

source('RCode/functions/raw2speclib_hsdar.R')
source('RCode/functions/DropCatVar_July2017.R')
source('RCode/functions/ModelSelect.R')
source('RCode/functions/LMMRindex.R')
source('RCode/functions/LMMRloop.R')
source('RCode/functions/index2prob_test.R')

dir.create('data', FALSE, FALSE)
dir.create('data/raw', FALSE, FALSE)
dir.create('data/temp', FALSE, FALSE)
dir.create('output/results', FALSE, FALSE)
dir.create('output/figures', FALSE, FALSE)
dir.create('paper', FALSE, FALSE)
dir.create('paper/templates', FALSE, FALSE)
dir.create('RCode', FALSE, FALSE)
dir.create('RCode/functions', FALSE, FALSE)

####
# 2. Getting raw spectra in and rename factor levels
####

ori.data <- read.csv('data/raw/data.wo.out.binned.cut.csv')

ori.data$Type <-
    revalue(
        ori.data$Type,
        c(
            "Treated" = "Uninfected",
            "Untreated" = "Infected",
            "Healthy" = "Healthy"
        )
    )
levels(ori.data$Type)

feature.set <- readRDS('data/raw/featuresforindex.rds')

####
# 3. Drop healthy leaves; convert treated/untreated to 0/1 coding; logs frequencies
####

data.log <- DropClass(ori.data, ori.data$Type, 'Healthy')
data.log$Type <-
    as.numeric(data.log$Type == 'Infected') # 1 = treated; 0 = untreated
data.log[names(data.log)[-1]] <-
    log(data.log[names(data.log)[-1]]) # logs frequencies from second to last col

#feature.set <- VSURF(data.log[,2:2152], data.log[,1], clusterType = "FORK", ntree = 2000, mtry = 50)

####
# 4. Select best model and select important predictor vars from original data
####

model.1 <- GLM.featureselect(data.log, feature.set)

best.bands <-
    row.names(summary(model.1)$coefficients)[c(2, 3, 4, 5)]
df.bands <-
    cbind('Type' = ori.data$Type, ori.data[best.bands])
df.bands <- DropClass(df.bands, df.bands$Type, "Healthy")

index_vals <-
    index4col(df.bands[, 2:5], model.1) #use function to calculate index values
result_df <-
    cbind('Type' = df.bands$Type, index_vals) #create final dataframe including index values
result_df <- rename(result_df, LMMR = index_val)

modelcoeffs <- coef(model.1)
saveRDS(modelcoeffs, 'coefficients.rds')
saveRDS(best.bands, 'bands.rds')

####
# 5. Add other indices to the result dataframe
####
tospectra <-
    read.csv("data/raw/data.wo.out.binned.cut.csv", check.names = FALSE)
tospectra <- DropClass(tospectra, tospectra$Type, "Healthy")

spectra <- raw2speclib(tospectra)

#rm(tospectra, index_vals)# remove redundant objects

####
# 6. Add indices to result_df and run through logit to get probs
####

#ARI <-  '((R550)^−1) − ((R700)^−1)' #Anthocyanin Reflectance Index
#RVSI <- '(((R712) + (R752))/2) − (R732)' #Red_Edge Vegetation Stress Index
NBNDVI <- '(R850 − R680)/(R850 + R680)'
#SIPI <- '(R800-R445)/(R800+R680)'

ind <- c("PRI", "MCARI", NBNDVI) #add more indices HERE


for (i in ind) {
    result_df[[paste(i, 'prob', sep = "_")]] <-
        Index2Prob(result_df, spectra, i)
    
}

ind.c <- c('PRI', 'MCARI', 'NBNDVI') #manipulate if necessary

colnames(result_df)[3:length(result_df)] <- ind.c

####
# 7. Visualize Results
####

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
    as.numeric(result_df$Type == 'Infected') # 1 = treated; 0 = untreated

myDat <- cbind(ID = seq(1, dim(result_df)[1]), result_df)

myCMX <- cmx(myDat)

acc <- presence.absence.accuracy(myDat)

acc_fin <- acc[, c(1, 3, 4, 5, 6, 7)]

saveRDS(acc_fin, 'Table1.rds')