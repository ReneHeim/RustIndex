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

# install.packages(c("cowplot", "gdata", "glmulti", "hsdar", "plyr",
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


# 2. Loading Functions and set project structure-----------------------------------------------------

source('R/20170601_FUN_DropCatVar.R')
source('R/20170601_FUN_exportVSURF.R')
source('R/20171224_FUN_raw2speclibhsdar.R')
source('R/20171224_FUN_index2prob.R')
source('R/20171224_FUN_prepggwide2long.R')

-dir.create('data', FALSE, FALSE)
-dir.create('R', FALSE, FALSE)
-dir.create('doc', FALSE, FALSE)
-dir.create('output', FALSE, FALSE)

# 3. Loading and Preparing Data ---------------------------------------------------------------------

ori.data <- read.csv('data/data.wo.out.binned.cut.csv') #Get original data

levels(ori.data$Type) #Check levels of categorical response 

data.log <- DropClass(ori.data, ori.data$Type, 'Healthy') #Drop factor level 'Healthy'

Type <- data.log$Type

data.log$Type <-
  as.numeric(data.log$Type == 'Untreated') #Transform remaining factor levels into 1 and 0 

data.log[names(data.log)[-1]] <-
  log(data.log[names(data.log)[-1]]) #Logs frequencies (R's log() computes natural logarithm)


# 4. Selecting Subset of Relevant Wavebands ---------------------------------------------------------

set.seed(20180108)

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
  band.vectors[[i]] <- export.VSURF(feature.set[[i]]$varselect.pred, data.log[, 2:202])
} # export.VSURF turns list of column numbers into list of according wavebands

VSURF.selection <- unlist(band.vectors) #Unlist list to create a vector containing all selected bands
VSURF.selection <- sort(unique(VSURF.selection)) #Only select unique bands from vector and sort

# 5. Model Selection to find 4 best bands to explain seperation between healthy and treated-----------

multi.model <- glmulti(y=names(data.log)[1],xr=paste0('X',VSURF.selection),data.log,maxsize=4,level = 1, family=binomial)

# 5.1 glmulti does not provide model coefficients, therefore add another logistic regression to
#     find coefficients for the spectral index

model.1 <- glm(as.formula(summary(multi.model)$bestmodel),data.log,family=binomial)

# 5.2 Build LMMR complex equation

coefficients(model.1)
best.bands <- row.names(summary(model.1)$coefficients)[c(2, 3, 4, 5)] # Extract best bands

#P=(exp(coef(model.1)[1]+(coef(model.1)[2]*log(best.bands[1]))+(coef(model.1)[3]*log(best.bands[2]))+(coef(model.1)[4]*log(best.bands[3]))+(coef(model.1)[5]*log(best.bands[4]))))/(1+(exp(coef(model.1)[1]+(coef(model.1)[2]*log(best.bands[1]))+(coef(model.1)[3]*log(best.bands[2]))+(coef(model.1)[4]*log(best.bands[3]))+(coef(model.1)[5]*log(best.bands[4])))))

# 5.3 Do 95% confint for coefficient pairs overlap?

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


# 6. Build spectral library and compute indices ------------------------------------------------------

# A) Build spectral library

tospectra <- read.csv("data/data.wo.out.binned.cut.csv", check.names = FALSE)

tospectra <- DropClass(tospectra, tospectra$Type, "Healthy")

spectra <- raw2speclib(tospectra) # Use hsdar to build spectral library

# B) Define spectral vegetation indices according to 'hsdar' pkg style

LMMR.complex <- '(exp(coef(model.1)[1]+(coef(model.1)[2]*log(R545))+(coef(model.1)[3]*log(R555))+(coef(model.1)[4]*log(R1505))+(coef(model.1)[5]*log(R2195))))/(1+(exp(coef(model.1)[1]+(coef(model.1)[2]*log(R545))+(coef(model.1)[3]*log(R555))+(coef(model.1)[4]*log(R1505))+(coef(model.1)[5]*log(R2195)))))'
LMMR_complex <- vegindex(spectra, LMMR.complex)
# LMMRcomplex values already are probabilities as the equation was built using a logistic regression
# All other indices must be submitted to a logistic regression to convert their values to probabilities

NBNDVI <- '(R850-R680)/(R850+R680)'
LMMR.simple <- '((R545/R555)^(5/3))*(R1505/R2195)'

index <- c("PRI", "MCARI", NBNDVI, LMMR.simple) # Add self-defined or hsdar pkg indices here (?vegindex)

# C) Compute spectral index values for all indices (except LMMR.complex, see 6. B)

index.list <- list()

index.list[['PRI']] <- vegindex(spectra, index[1])
index.list[['MCARI']] <- vegindex(spectra, index[2])
index.list[['NBNDVI']] <- vegindex(spectra, index[3])
index.list[['LMMR.simple']] <- vegindex(spectra, index[4])

index.df <- do.call(cbind.data.frame, index.list)

# D) Each column in 'index.df' must be passed through a logistic regression to yield probability values

index.prob.list <- list() # This list will contain all logistic regression models

index.prob.list[['PRI.model']] <- glm(data.log[,1]~index.df$PRI, family=binomial(link=probit))
index.prob.list[['MCARI.model']] <- glm(data.log[,1]~index.df$MCARI, family=binomial(link=probit))
index.prob.list[['NBNDVI.model']] <- glm(data.log[,1]~index.df$NBNDVI, family=binomial(link=probit))
index.prob.list[['LMMR.simple.model']] <- glm(data.log[,1]~index.df$LMMR.simple, family=binomial(link=probit))

# E) Build a df containing LMMR.complex plus probability values stored under each model in index.prob.list

result.df <- as.data.frame(Type)
result.df <- cbind(result.df, LMMR_complex)
result.df <- cbind(result.df, 'PRI'= index.prob.list$PRI.model$fitted.values)
result.df <- cbind(result.df, 'MCARI'= index.prob.list$MCARI.model$fitted.values)
result.df <- cbind(result.df, 'NBNDVI'= index.prob.list$NBNDVI.model$fitted.values)
result.df <- cbind(result.df, 'LMMR.simple'=index.prob.list$LMMR.simple.model$fitted.values)

# F) Is the simplified LMMR similar to the complex LMMR?

# D) Explore the logistic regression model of the LMMR


coefficients(index.prob.list$LMMR.simple.model) # Find in Equation 4 in article

compare <- (index.prob.list$LMMR.simple.model$aic)-(model.1$aic) # Shows similarity of complex and simplified LMMR

compare


write_csv(result.df,'output/vegindex_df.csv') # This file contains all used indices converted into 
# probability values

# 7. Visualize Results ------------------------------------------------------------------------------

# A) Create box plots/jitter plots to visually compare LMMR.complex, LMMR.simple, PRI, MCARI and NBNDVI

plot_list <- list()
indi <- names(result.df[, 2:length(result.df)])

for (i in indi) {
  plot_list[[i]] <-
    ggplot(result.df,
           aes_string(result.df$Type, result.df[, i], colour = result.df$Type)) +
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
  plot_list[[1]],
  plot_list[[5]],
  plot_list[[2]],
  plot_list[[3]],
  plot_list[[4]],
  labels = c("a", 'b', 'c', 'd', 'e'),
  ncol = 5,
  nrow = 1
)
p1


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

# C) Get bar/jitter plots from 7A and modify to combine with 7B

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
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank())+
  labs(x = "LMMR.s")

p9 <- plot_list[[1]]
p9 <- p9+
  theme(#axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())+
  labs(x = "LMMR.c", y="Disease Prob.")

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

# 8. Calculate Accuracy Metrics -------------------------------------------------------

# A) Get AIC values from logit models

aics <- list()

aics['LMMRcomplex'] <- model.1$aic
aics['LMMRsimple'] <- index.prob.list$LMMR.simple.model$aic
aics['PRI'] <- index.prob.list$PRI.model$aic
aics['MCARI'] <- index.prob.list$MCARI.model$aic
aics['NBNDVI'] <- index.prob.list$NBNDVI.model$aic

aicvec <- unlist(aics)
DeltaAIC <- round(aicvec-aicvec[1], 3)


# B) Extract confusion matrices for Table 1

confmat.LMMRc<- cmx(myDat[,c(1,2,3)]) #Select columns here to get confusion matrix
confmat.LMMRs<- cmx(myDat[,c(1,2,7)])
confmat.PRI<- cmx(myDat[,c(1,2,4)])
confmat.MCARI<- cmx(myDat[,c(1,2,5)])
confmat.NBNDVI<- cmx(myDat[,c(1,2,6)])

write.csv(confmat.LMMRc, 'output/confmatLMMRc.csv', row.names = FALSE)
write.csv(confmat.LMMRc, 'output/confmatLMMRs.csv', row.names = FALSE)
write.csv(confmat.PRI, 'output/confmatPRI.csv', row.names = FALSE)
write.csv(confmat.MCARI, 'output/confmatMCARI.csv', row.names = FALSE)
write.csv(confmat.NBNDVI, 'output/confmatNBNDVI.csv', row.names = FALSE)

# C) Design accuracy assessment for Table 2

result.temp <- result.df

result.temp$Type <- ifelse(result.temp$Type=="Untreated",1,0) # 1 = untreated; 0 = treated

myDat <- cbind(ID = seq(1, dim(result.temp)[1]), result.temp)


acc <- presence.absence.accuracy(myDat)

acc_fin <- acc[, c(1, 3, 4, 5, 6, 7)] 

acc_fin[, c(2:6)] <- round(acc_fin[,c(2:6)], 3)

acc_fin[c(1,2,3,4,5),] <- acc_fin[c(1,5,3,4,2),]

Table2 <- cbind(acc_fin[,2:5],DeltaAIC)

saveRDS(Table2, 'output/Table1.rds')

write.csv(Table2, 'output/Table1.csv')




