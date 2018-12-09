####
# This code is generating the LMMR index (LemonMyrtle-MyrtleRust), a spectral disease index for the 
# pathosystem Austropuccinia psidii and Backhousia citriodora. It utilizes a spectral dataset 
# containing following columns:
#
# Type: Categorical variables referring to spectral classes (Healthy, Treated and Untreated). Data was
# already cleaned and manipulated in Heim et al. (2018a).
#
# X505-X2500: Spectral reflectance [%]at a specific waveband.
#
# Please refer to Heim et al. (2018b) for more information.
#
# Heim, R.H.J., Wright, I.J., Chang, H.-C., Carnegie, A.J., Pegg, G.S., Lancaster, E.K., Falster, D.S.,
# Oldeland, J., 2018a. Detecting myrtle rust ( Austropuccinia psidii ) on lemon myrtle trees using 
# spectral signatures and machine learning. Plant Pathol. https://doi.org/10.1111/ppa.12830
#
# Heim, RHJ; Wright, IJ; Allen, AP; Geedicke, I and Oldeland, J. (2018) 
# Developing a spectral disease index for myrtle rust (*Austropuccinia psidii*), 
# Plant Pathology, XX(xx), pp. XXX. doi: XXX.
#
# The following code is structured according to the Method section (Fig. 2 A and B) in 
# Heim et al. (2018b)
####


# 1. Install and Load Packages ----------------------------------------------------------------------

# install.packages(c("cowplot", "gdata", "glmulti", "hsdar", "plyr",
#                    "PresenceAbsence", "prospectr", "rJava", "tidyverse",
#                    "VSURF", "reshape2", "caret"))

#Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_151') #Set path to Java dir for rJava

# Please install correct version of Java (32 or 64 bit) before setting the path to the Java dir.

library(caret)
library(cowplot)

library(devtools)

library(gdata)
library(ggplot2)
library(glmulti)

library(hsdar)

library(magrittr)

library(plyr)
library(PresenceAbsence)
library(prospectr)

library(rJava)
library(reshape2)

library(VSURF)

# 2. Loading Functions and set project structure-----------------------------------------------------

source('R/FUN_drop_cat_var.R')
source('R/FUN_exportVSURF.R')
source('R/FUN_raw2speclibhsdar.R')
source('R/FUN_prepggwide2long.R')

dir.create('output', FALSE, FALSE) # Code generated output

############################################################################################
# A SECTION - Raw data to linear model
############################################################################################

# A 3. Loading and Preparing Data ---------------------------------------------------------------------

ori.data <- read.csv('data/data.wo.out.binned.cut.csv') #Get original data

levels(ori.data$Type) #Check levels of categorical response 

data <- drop_class(ori.data, ori.data$Type, 'Healthy') #Drop factor level 'Healthy'

Type <- data$Type # Exract for later use

data$Type <-
  as.numeric(data$Type == 'Untreated') #Transform remaining factor levels into 1 and 0 

data[names(data)[-1]] <-
  log(data[names(data)[-1]]) #Logs frequencies (R's log() computes natural logarithm)


# A 4. Feature Selection (202 to 12 bands) -------------------------------------------------

#set.seed(20180111)

# A) Selection was repeated 10x to account for random selection encounters (~ 30h runtime)

# feature.set <- list()
# runs <- seq(1,10,1)
# for(i in runs){
#   
#   feature.set[[i]] <-
#          VSURF(data.log[,2:202], data.log[,1], clusterType = "FORK", ntree = 2000,
#                mtry = 50) #Takes ~3h, therefore saved/loaded as .rds
#   
# }
# saveRDS(feature.set, 'data/features.rds')

feature.set <- readRDS('data/features.rds') # Reload stored feature selection object

# C) VSURF outputs only important column numbers (1,2..), 
#    therefore turn into waveband names (500, 505...)

runs <- seq(1,length(feature.set),1) # Create sequence depending on length of feature.set
band.vectors <- list() # Create output object

for(i in runs){
  band.vectors[[i]] <- export_vsurf(feature.set[[i]]$varselect.pred, data[, 2:202])
} # export.VSURF turns list of column numbers into list of according wavebands

VSURF.selection <- unlist(band.vectors) #Unlist list to create a vector containing all selected bands
VSURF.selection <- sort(unique(VSURF.selection)) #Only select unique bands from vector and sort

# A5. Model Selection to find 4 best bands to explain seperation between healthy and treated-----------
set.seed(20180117)
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

# Based on the LMMR model we designed the LMMR index through mathematical simplification.
# We initially log transformed the spectral reflectance to apply Eq. 2 and 3 
# (Heim et al., 2018b). Additionally, to summarize the model coefficients and yield
# 5/3 as a simplification, the absolute 95% confidence intervals for coefficient pairs 
# should overlap as this indicates products of approximately equal magnitudes.

# 5.3 Plot confidence intervals

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

############################################################################################
# B SECTION - Linear Model to classification report
############################################################################################

# Please refer to the according  article for the required steps yielding the LMMR.

# B 6. Compare simplified LMMR and other common indices-----------------------------------------

# A) Build spectral library

tospectra <- read.csv("data/data.wo.out.binned.cut.csv", check.names = FALSE)

tospectra <- drop_class(tospectra, tospectra$Type, "Healthy")

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
index.df$Type <- Type # index.df is the new df for the classification

# D) Logistic regression (classification) on each index outcome (index.df) using caret pkg

inTrain <- createDataPartition(y = index.df$Type, p = .75, list = FALSE)
Train.75 <- index.df[inTrain,]
Test.25 <- index.df[-inTrain,]

ctrl <- trainControl(method = "boot",
                     number = 1000, 
                     classProbs = FALSE, 
                     savePredictions = TRUE)

# PRI Training
PRI.model <- train(Type ~ PRI, data = Train.75, method = "glm", trControl = ctrl, metric = c('Kappa'))
summary(PRI.model)# estimates
confusionMatrix(PRI.model)
PRIvalues <- PRI.model$finalModel$fitted.values
PRI.model$pred
# PRI Testing
PRI.pred <- predict(PRI.model, newdata = Test.25)

sink("output/PRITest25.txt", append=FALSE, split=FALSE)
confusionMatrix(data = PRI.pred, Test.25$Type)
sink()
confmat.PRI <- confusionMatrix(data = PRI.pred, Test.25$Type)


# MCARI Training
MCARI.model <- train(Type ~ MCARI, data = Train.75, method = "glm", trControl = ctrl, metric = c("Kappa"))
summary(MCARI.model)
confusionMatrix(MCARI.model)
MCARIvalues <- MCARI.model$finalModel$fitted.values
# MCARI Test
MCARI.pred <- predict(MCARI.model, newdata = Test.25)

sink("output/MCARITest25.txt", append=FALSE, split=FALSE)
confusionMatrix(data = MCARI.pred, Test.25$Type)
sink()
confmat.MCARI <- confusionMatrix(data = MCARI.pred, Test.25$Type)

# NBNDVI Training
NBNDVI.model <- train(Type ~ NBNDVI, data = Train.75, method = "glm", trControl = ctrl, metric = c("Kappa"))
summary(NBNDVI.model)
confusionMatrix(NBNDVI.model)
NBNDVIvalues <- NBNDVI.model$finalModel$fitted.values
# NBNDVI Test
NBNDVI.pred <- predict(NBNDVI.model, newdata = Test.25)

sink("output/NBNDVITest25.txt", append=FALSE, split=FALSE)
confusionMatrix(data = NBNDVI.pred, Test.25$Type)
sink()
confmat.NBNDVI <- confusionMatrix(data = NBNDVI.pred, Test.25$Type)

# LMMR Training
LMMR.model <- train(Type ~ LMMR, data = Train.75, method = "glm", trControl = ctrl, metric = c("Kappa"))
summary(LMMR.model)
confusionMatrix(LMMR.model)
LMMRvalues <- LMMR.model$finalModel$fitted.values
# LMMR Test
LMMR.pred <- predict(LMMR.model, newdata = Test.25)

sink("output/LMMRTest25.txt", append=FALSE, split=FALSE)
confusionMatrix(data = LMMR.pred, Test.25$Type)
sink()
confmat.LMMR <- confusionMatrix(data = LMMR.pred, Test.25$Type)
# E) Create df from training fitted values

train.res.df <- as.data.frame(Train.75$Type)
train.res.df$PRIval <- PRIvalues
train.res.df$MCARIval <- MCARIvalues
train.res.df$NBNDVIval <- NBNDVIvalues
train.res.df$LMMRval <- LMMRvalues

names(train.res.df) <- c('Type', 'PRI', 'MCARI', 'NBNDVI', 'LMMR')

# 7. Visualize Results ------------------------------------------------------------------------------
indi <- names(train.res.df[, 2:length(train.res.df)])
plot_list <- list()
for (i in indi) {

plot_list[[i]] <- ggplot(train.res.df, aes_string(train.res.df$Type, train.res.df[, i], colour = train.res.df$Type))+
  geom_jitter(aes(shape=Type),height = 0)+
  geom_boxplot(colour = 'black', alpha = 0.5, outlier.alpha = 0)+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none")+
        xlab(label = paste(i))+
        ylab(label = 'Index Value')+
        ylim(0, 1)
}


p1 <- plot_grid(
  plot_list[[4]],
  plot_list[[2]],
  plot_list[[3]],
  plot_list[[1]],
  labels = c("a", 'b', 'c', 'd'),
  ncol = 4,
  nrow = 1)

p1
  
  
# Spectra plots

# B) Plot spectra and show final most important wavebands
    
spectra.gg <- prep_gg(tospectra) # Transforms wide to long for ggplot2 readibility
  
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
  theme(legend.position = c(.90, .88), 
        legend.title = element_blank(), 
        legend.background = element_blank())+
  labs(x = "Wavelength [nm]", y = "Reflectance [%]")

pspec
  
# C) Get bar/jitter plots from 7A and modify to combine with 7B

p5 <- plot_list[[2]]#MCARI
p5 <- p5+
  theme(#axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y=element_blank())

p6 <- plot_list[[3]]#NBNDVI
p6 <- p6+
    theme(#axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y=element_blank())

p7 <- plot_list[[4]]#LMMR
p7 <- p7+
    theme(#axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y=element_blank())

p9 <- plot_list[[1]]#PRI
p9 <- p9+
    theme(#axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())+
      labs(x = "PRI", y="Disease Prob.")
  
# D) Build Figure 2 (Figure 1 was designed outside of R using photographs and Inkscape)
    
plot.res <- ggdraw() +
    draw_plot(p9, x = 0, y = .5, width = .25, height = .5) +
    draw_plot(p5, x = .25, y = .5, width = .25, height = .5) +
    draw_plot(p6, x = .5, y = .5, width = .25, height = .5) +
    draw_plot(p7, x = .75, y = .5, width = .25, height = .5) +
      draw_plot(pspec, x = 0, y = 0, width = 1, height = 0.5) +
        draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 12,
                        x = c(0, 0.25, 0.5, 0.75, 0), y = c(1, 1, 1, 1, 0.5))
  
plot.res
  
ggsave("output/Figure2.boxspectra.png",
    plot = plot.res,
    width = 40,
    height = 20,
    units = "cm",
    dpi = 400
  )
