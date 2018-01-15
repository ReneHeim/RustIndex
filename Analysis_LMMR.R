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


# 1. Install and Load Packages ----------------------------------------------------------------------

# install.packages(c("devtools","cowplot", "gdata", "glmulti", "hsdar", "plyr",
#                    "PresenceAbsence", "prospectr", "rJava", "tidyverse",
#                    "VSURF", "reshape2", "caret"))

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
library(devtools)

devtools::install_github("sachsmc/plotROC")
library(plotROC)


# 2. Loading Functions and set project structure-----------------------------------------------------

source('R/20170601_FUN_DropCatVar.R')
source('R/20170601_FUN_exportVSURF.R')
source('R/20171224_FUN_raw2speclibhsdar.R')
source('R/20171224_FUN_prepggwide2long.R')

-dir.create('data', FALSE, FALSE) # Contains original data only (Do not modify!)
-dir.create('R', FALSE, FALSE) # Contains functions
-dir.create('output', FALSE, FALSE) # Code generated output

# 3. Loading and Preparing Data ---------------------------------------------------------------------

ori.data <- read.csv('data/data.wo.out.binned.cut.csv') #Get original data

levels(ori.data$Type) #Check levels of categorical response 

data <- DropClass(ori.data, ori.data$Type, 'Healthy') #Drop factor level 'Healthy'

Type <- data$Type # Exract for later use

data$Type <-
  as.numeric(data$Type == 'Untreated') #Transform remaining factor levels into 1 and 0 

data[names(data)[-1]] <-
  log(data[names(data)[-1]]) #Logs frequencies (R's log() computes natural logarithm)


# 4. Selecting Subset of Relevant Wavebands ---------------------------------------------------------

set.seed(20180111)

# A) 10x VSURF Feature Selection (Runs ~30 hours)

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

feature.set <- readRDS('output/features.rds') # List conatining 10 feature selection objects

# C) VSURF outputs only important column numbers (1,2..), 
#    therefore turn into waveband names (500, 505...)

runs <- seq(1,10,1) # Create sequence depending on length of feature.set
band.vectors <- list() # Create output object

for(i in runs){
  band.vectors[[i]] <- export.VSURF(feature.set[[i]]$varselect.pred, data[, 2:202])
} # export.VSURF turns list of column numbers into list of according wavebands

VSURF.selection <- unlist(band.vectors) #Unlist list to create a vector containing all selected bands
VSURF.selection <- sort(unique(VSURF.selection)) #Only select unique bands from vector and sort

# 5. Model Selection to find 4 best bands to explain seperation between healthy and treated-----------

multi.model <- glmulti(y=names(data)[1],xr=paste0('X',VSURF.selection),data,maxsize=4,level = 1, family=binomial)

# 5.1 glmulti does not provide model coefficients, therefore add another logistic regression to
#     find coefficients for the spectral index

model.1 <- glm(as.formula(summary(multi.model)$bestmodel),data,family=binomial)

saveRDS(model.1, 'output/LMMRmodel.RDS')
model.1 <- readRDS(file = 'output/LMMRmodel.RDS')

# 5.2 Build LMMR complex equation

coefficients(model.1)
best.bands <- row.names(summary(model.1)$coefficients)[c(2, 3, 4, 5)] # Extract best bands


LMMR.model.eq <- 'log[P/(1 - P)] = 18.387 + 75.382 log[R545] - 78.809 log[R555] + 45.993 log[R1505] - 46.831 log[R2195]'

# LMMR model is the foundation for the simplified LMMR index (see article)

# 5.3 Do absolute 95% confint for coefficient pairs overlap? (Requirement to build ratio index)

confint(model.1)

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


# 6. Compare simplified LMMR and other common indices-----------------------------------------

# A) Build spectral library

tospectra <- read.csv("data/data.wo.out.binned.cut.csv", check.names = FALSE)

tospectra <- DropClass(tospectra, tospectra$Type, "Healthy")

spectra <- raw2speclib(tospectra) # Use hsdar to build spectral library

# B) Define spectral vegetation indices and LMMR 

NBNDVI <- '(R850-R680)/(R850+R680)'
LMMR <- '((R545/R555)^(5/3))*(R1505/R2195)'

index <- c("PRI", "MCARI", NBNDVI, LMMR) # Add self-defined or hsdar pkg indices here (?vegindex)

# C) Compute spectral index values for all indices

index.list <- list()

index.list[['PRI']] <- vegindex(spectra, index[1])
index.list[['MCARI']] <- vegindex(spectra, index[2])
index.list[['NBNDVI']] <- vegindex(spectra, index[3])
index.list[['LMMR']] <- log(vegindex(spectra, index[4])) # Needs to be transformed 
                                                                # to the log scale as it was
                                                                # developed on that scale.
index.df <- do.call(cbind.data.frame, index.list)
index.df$Type <- Type

# D) Logistic regression on each index outcome (index.df) using caret pkg

inTrain <- createDataPartition(y = index.df$Type, p = .75, list = FALSE)
Train.75 <- index.df[inTrain,]
Test.25 <- index.df[-inTrain,]

ctrl <- trainControl(method = "repeatedcv",repeats = 10, classProbs = TRUE,summaryFunction = twoClassSummary)

# PRI
PRI.model <- train(Type ~ PRI, data = Train.75, method = "glm", trControl = ctrl, metric = 'Kappa')
summary(PRI.model)# estimates
PRI.pred <- predict(PRI.model, newdata = Test.25)
summary(PRI.pred)

confmat.PRI <- confusionMatrix(data = PRI.pred, Test.25$Type)

# MCARI
MCARI.model <- train(Type ~ MCARI, data = Train.75, method = "glm", trControl = ctrl, metric = c("Kappa"))
MCARI.model
MCARI.pred <- predict(MCARI.model, newdata = Test.25)
confmat.MCARI <- confusionMatrix(data = MCARI.pred, Test.25$Type)

# NBNDVI
NBNDVI.model <- train(Type ~ NBNDVI, data = Train.75, method = "glm", trControl = ctrl, metric = c("Kappa"))
NBNDVI.model
NBNDVI.pred <- predict(NBNDVI.model, newdata = Test.25)
confmat.NBNDVI <- confusionMatrix(data = NBNDVI.pred, Test.25$Type)

# LMMR
LMMR.model <- train(Type ~ LMMR, data = Train.75, method = "glm", trControl = ctrl, metric = c("Kappa"))
LMMR.model
LMMR.pred <- predict(LMMR.model, newdata = Test.25)
confmat.LMMR <- confusionMatrix(data = LMMR.pred, Test.25$Type)

confmat.LMMR$overall[2]
# E) Create Figure 1 - Training data/ Model evaluation (plotROC pkg)

long.Train75 <- melt_roc(Train.75, "Type", c("PRI", "MCARI", "NBNDVI", "LMMR"))
names(long.Train75) <- c('Type', 'Value', 'Index')

ggplot(long.Train75, aes(d = Type, m = Value, color = Index)) + 
  geom_roc(n.cuts = 50, labels = FALSE, aes(shape=Index)) + 
  style_roc(xlab = "1 - Specificity", ylab="Sensitivity")


# F) Create Table 1 - Test data accuracy metrics

OA <- round(c('PRI'= confmat.PRI$overall[1], 
        'MCARI'= confmat.MCARI$overall[1], 
        'NBNDVI'= confmat.NBNDVI$overall[1],
        'LMMR'= confmat.LMMR$overall[1]),3)*100

Kappa <- round(c('PRI'= confmat.PRI$overall[2], 
        'MCARI'= confmat.MCARI$overall[2], 
        'NBNDVI'= confmat.NBNDVI$overall[2],
        'LMMR'= confmat.LMMR$overall[2]), 3)*100

Sensitivity <- round(c('PRI'= confmat.PRI$byClass[1], 
        'MCARI'= confmat.MCARI$byClass[1], 
        'NBNDVI'= confmat.NBNDVI$byClass[1],
        'LMMR'= confmat.LMMR$byClass[1]), 3)*100

Specificity <- round(c('PRI'= confmat.PRI$byClass[2], 
                 'MCARI'= confmat.MCARI$byClass[2], 
                 'NBNDVI'= confmat.NBNDVI$byClass[2],
                 'LMMR'= confmat.LMMR$byClass[2]), 3)*100

testresult.df <- as.data.frame(cbind(OA, Kappa, Sensitivity, Specificity))

names(testresult.df) <- c('OA[%]', 'Kappa[%]', 'Sensitivity[%]', 'Specificity[%]')
row.names(testresult.df) <- c('PRI', 'MCARI', 'NBNDVI', 'LMMR')
# 7. Visualize Results ------------------------------------------------------------------------------




# B) Plot spectra and show final most important wavebands

spectra.gg <- prep.gg(tospectra) # Transforms wide to long for ggplot2 readibility

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

# C) Build Zoomed version

# D) Build Figure 2 (Figure 1 was designed outside of R using photographs and Inkscape)

plot.res <- ggdraw() +
  draw_plot(p9, x = 0, y = .5, width = .2, height = .5) +
  draw_plot(p8, x = 0.2, y = .5, width = .2, height = .5) +
  draw_plot(p5, x = .4, y = .5, width = .2, height = .5) +
  draw_plot(p6, x = .6, y = .5, width = .2, height = .5) +
  draw_plot(p7, x = .8, y = .5, width = .2, height = .5) +
  draw_plot(pspec, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F"), size = 12,
                  x = c(0, 0.2, 0.4, 0.6,0.8, 0), y = c(1, 1, 1, 1, 1, 0.5))
plot.res

ggsave(
  "output/20180108Results.png",
  plot = plot.res,
  width = 40,
  height = 20,
  units = "cm",
  dpi = 400
)



