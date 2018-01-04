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


# 1. Loading Packages and Functions------------------------------------------------


install.packages(c("cowplot", "gdata", "glmulti", "hsdar", "plyr",
                   "PresenceAbsence", "prospectr", "rJava", "tidyverse",
                   "VSURF", "reshape2", "caret"))

#Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_151') #Set path to Java dir for rJava

# Please install correct version of Java (32 or 64 bit) before setting the path to the Java dir.

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
library(caret)


source('R/20171224_FUN_raw2speclibhsdar.R')
source('R/20170601_FUN_DropCatVar.R')
source('R/20171224_FUN_ModelSelect.R')
source('R/2011224_FUN_LMMRindex.R')
source('R/20171224_FUN_LMMRloop.R')
source('R/20171224_FUN_index2prob.R')
source('R/20170601_FUN_exportVSURF.R')
source('R/20171224_FUN_prepggwide2long.R')

#Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_151') #Set path to Java dir for rJava

# 2. Set Project Structure --------------------------------------------------------

dir.create('data', FALSE, FALSE)
dir.create('R', FALSE, FALSE)
dir.create('doc', FALSE, FALSE)
dir.create('output', FALSE, FALSE)


# 3. Loading and Preparing Data ---------------------------------------------------

ori.data <- read.csv('data/data.wo.out.binned.cut.csv') #Get original data

levels(ori.data$Type) #Check levels of categorical response 

data.log <- DropClass(ori.data, ori.data$Type, 'Healthy') #Drop factor level 'Healthy'

data.log$Type <-
    as.numeric(data.log$Type == 'Untreated') #Transform remaining factor levels into 1 and 0 

data.log[names(data.log)[-1]] <-
    log(data.log[names(data.log)[-1]]) #Logs frequencies (R's log() computes natural logarithm)


# 4. Selecting Subset of Relevant Wavebands ---------------------------------------

    # A) VSURF Feature Selection (Runs ~30 hours) ## really? not only 3 see line 84

# feature.set <- list()
# runs <- seq(1,10,1)
# for(i in runs){
#   
#   feature.set[[i]] <-
#          VSURF(data.log[,2:202], data.log[,1], clusterType = "FORK", ntree = 2000,
#                mtry = 50) #Takes ~3h, therefore saved/loaded as .rds
#   
# }
# saveRDS(feature.set, 'output/features.rds')

feature.set <- readRDS('output/features.rds') # B) Load features as A) runs a while

    # C) VSURF outputs column indices (1,2..), therefore turn into waveband names (500, 505...)

runs <- seq(1,10,1) #Create sequence depending on length of feature.set
res <- list() #Create output object

for(i in runs){
  res[[i]] <- export.VSURF(feature.set[[i]]$varselect.pred, data.log[, 2:202])
  }

VSURF.selection <- unlist(res) #Unlist list to create a vector containing all selected bands
VSURF.selection <- sort(unique(VSURF.selection)) #Only select unique bands from vector and sort


# 5. Model Selection ---------------------------------------------------------

model.1 <- GLM.featureselect(data.log, VSURF.selection)

best.bands <-
  row.names(summary(model.1)$coefficients)[c(2, 3, 4, 5)]

df.bands <-
  cbind('Type' = ori.data$Type, ori.data[best.bands]) #Select only best bands cols from df

df.bands <- DropClass(df.bands, df.bands$Type, "Healthy")

index_vals <-
  index4col(df.bands[, 2:5], model.1) # Function that generates LMMR values

result_df <-
  cbind('Type' = df.bands$Type, index_vals) # Attach LMMR values to dataframe

names(result_df)[2] <- c('LMMR') #Rename col name

mod.out <- list()

mod.out[['BestBands']] <- best.bands
mod.out[['Coefficients']] <- coef(model.1)
mod.out[['IndexDF']] <- result_df

saveRDS(mod.out, 'output/modeloutput.rds') # Export for future use

    # A) Do 95% confint for coefficient pairs overlap?

y = coef(model.1)
x = seq_along(y)
ci = confint(model.1)
xlim = range(x) + c(-0.5,0.2)
ylim = range(ci)
ylim = ylim + 0.1*c(-1,+1)*diff(ylim) # extend it a little
ylab = bquote(hat(beta))
xlab = "coefficient"
par(mar=c(4.5,5,1,1), las=1)
plot(y, pch=16, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, xaxt="n", bty="n")
axis(1, at=x, labels=names(y), tick=FALSE)
abline(h=0, lty=3)
arrows(x,ci[,1],x,ci[,2], code=3, angle=90, length=0.05)


# 6. Build spectral library and compute indices ------------------------------

    # A) Build spectral library

tospectra <- read.csv("data/data.wo.out.binned.cut.csv", check.names = FALSE)

tospectra <- DropClass(tospectra, tospectra$Type, "Healthy")

spectra <- raw2speclib(tospectra) # Use hsdar to build spectral library


    # B) Define spectral vegetation indices according to 'hsdar' pkg style


#ARI <-  '((R550)^−1) − ((R700)^−1)' #Anthocyanin Reflectance Index
#RVSI <- '(((R712) + (R752))/2) − (R732)' #Red_Edge Vegetation Stress Index
NBNDVI <- '(R850-R680)/(R850+R680)'
#SIPI <- '(R800-R445)/(R800+R680)'
LMMR2 <- '((R545/R555)^(5/3))*(R1505/R2195)'

ind <- c("PRI", "MCARI", NBNDVI, LMMR2) # Add self-defined or hsdar pkg indices here (?vegindex)

    # C) Build df containing disease probabilities for each vegetation index

for (i in ind) {
    result_df[[paste(i, 'prob', sep = "_")]] <-
        Index2Prob(result_df, spectra, i)
    
} 

ind.c <- c('PRI', 'MCARI', 'NBNDVI', 'LMMRopt') # Choose suitable column names

colnames(result_df)[3:length(result_df)] <- ind.c

write_csv(result_df,'output/vegindex_df.csv')

# 7. Visualize Results ------------------------------------------------------------

plot_list <- list()
indi <- names(result_df[, 2:length(result_df)])

for (i in indi) {
    plot_list[[i]] <-
        ggplot(result_df,
               aes_string(result_df$Type, result_df[, i], colour = result_df$Type)) +
        geom_jitter(aes(shape=Type),height = 0) +
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


# Spectra plots


spectra.gg <- prep.gg(tospectra) # Transforms wide to long for ggplot2 readibility

# Subplot E (Spectra including spectral regions and best bands)

bands4gg <-as.numeric(gsub('X', '', best.bands))

pspec <- ggplot(spectra.gg, aes(Wavelength, Reflectance, colour = Type)) +
  geom_line(aes(linetype=Type), size = 1)+
  geom_point(aes(shape=Type), size = 2)+
  annotate(
    "rect",
    xmin = 500,
    xmax = 570,
    ymin = -Inf,
    ymax = Inf,
    alpha = .2,
    fill = 'green'
  ) +
  annotate(
    "rect",
    xmin = 570,
    xmax = 590,
    ymin = -Inf,
    ymax = Inf,
    alpha = .2,
    fill = 'yellow'
  ) +
  annotate(
    "rect",
    xmin = 590,
    xmax = 610,
    ymin = -Inf,
    ymax = Inf,
    alpha = .2,
    fill = 'orange'
  ) +
  annotate(
    "rect",
    xmin = 610,
    xmax = 700,
    ymin = -Inf,
    ymax = Inf,
    alpha = .2,
    fill = c("red")
  ) +
  annotate(
    "rect",
    xmin = 700,
    xmax = 1300,
    ymin = -Inf,
    ymax = Inf,
    alpha = .2,
    fill = c("lightgrey")
  ) +
  annotate(
    "rect",
    xmin = 1300,
    xmax = 2500,
    ymin = -Inf,
    ymax = Inf,
    alpha = .2,
    fill = 'white'
  ) +
  annotate(
    "text",
    x = 600,
    y = 30,
    label = "VIS",
    fontface = "bold",
    size = 5
  ) +
  annotate(
    "text",
    x = 1000,
    y = 30,
    label = "NIR",
    fontface = "bold",
    size = 5
  ) +
  annotate(
    "text",
    x = 1900,
    y = 30,
    label = "SWIR",
    fontface = "bold",
    size = 5
  ) +
  geom_vline(
    xintercept = bands4gg,
    col = "black",
    linetype = "twodash",
    size = 1,
    alpha = .5
  ) +
  theme_set(theme_bw(base_size = 20))+
  theme(legend.position = c(.90, .88), legend.title = element_blank(), legend.background = element_blank())

pspec

p5 <- plot_list[[2]]
p5 <- p5+
  theme(#axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank())
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
    axis.ticks.x=element_blank())+
    labs(x = "LMMR", y="Disease Prob.")

plot.res <- ggdraw() +
  draw_plot(p8, x = 0, y = .5, width = .25, height = .5) +
  draw_plot(p5, x = .25, y = .5, width = .25, height = .5) +
  draw_plot(p6, x = .5, y = .5, width = .25, height = .5) +
  draw_plot(p7, x = .75, y = .5, width = .25, height = .5) +
  draw_plot(pspec, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 12,
                  x = c(0, 0.25, 0.5, 0.75, 0), y = c(1, 1, 1, 1, 0.5))

ggsave(
  "output/20180103Results.png",
  plot = plot.res,
  width = 40,
  height = 20,
  units = "cm",
  dpi = 400
)

# Calculate Accuracy Metrics ---------------------------------------------------

    # A) Get AIC values from logit models

VILMMRopt <- vegindex(spectra, LMMR2)
VIPRI <- vegindex(spectra, 'PRI')
VIMCARI <- vegindex(spectra, 'MCARI')
VINBNDVI <- vegindex(spectra, NBNDVI)


AICtab <- cbind(result_df, VILMMRopt)
AICtab <- cbind(AICtab, VIPRI)
AICtab <- cbind(AICtab, VIMCARI)
AICtab <- cbind(AICtab, VINBNDVI)

LMMR.logit <- glm(AICtab[,1]~AICtab$VILMMRopt, family=binomial(link=probit))
LMMR.logit$aic

PRI.logit <- glm(AICtab[,1]~AICtab$VIPRI, family=binomial(link=probit))
PRI.logit$aic

MCARI.logit <- glm(AICtab[,1]~AICtab$VIMCARI, family=binomial(link=probit))
MCARI.logit$aic

NBNDVI.logit <- glm(AICtab[,1]~AICtab$VINBNDVI, family=binomial(link=probit))
NBNDVI.logit$aic

aics <- vector()

aics[1] <- LMMR.logit$aic
aics[2] <- PRI.logit$aic
aics[3] <- MCARI.logit$aic
aics[4] <- NBNDVI.logit$aic

DeltaAIC <- round(aics/aics[1], 3)
  
  
    # B) Design accuracy assessment table

result_df$Type <-
    as.numeric(result_df$Type == 'Untreated') # 1 = treated; 0 = untreated

myDat <- cbind(ID = seq(1, dim(result_df)[1]), result_df)

myCMX <- cmx(myDat) #Select columns here to get confusion matrix

acc <- presence.absence.accuracy(myDat)

acc_fin <- acc[, c(1, 3, 4, 5, 6, 7)] 

acc_fin[, c(2:6)] <- round(acc_fin[,c(2:6)], 3)

acc_fin[c(2,3,4,5),] <- acc_fin[c(5,2,3,4),]

acc_fin <- acc_fin[c(2:5),]
acc_fin

Table1 <- cbind(acc_fin,DeltaAIC)

saveRDS(Table1, 'output/Table1.rds')

write.csv(Table1, 'output/Table1.csv')

?cmx
