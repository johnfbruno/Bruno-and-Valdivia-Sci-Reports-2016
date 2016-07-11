# R Code for SREP29778
# Coral reef degradation is not correlated with local human population density
# Bruno, J. and Valdivia, A. 2016

#Analysis to compare reef state w human population density
#Response variables include coral cover, macroalgae cover, 
#Human population density comes from the SEDAC webpage


#Load World data
  world <- read.csv(file = "World.csv", header = TRUE, stringsAsFactors = TRUE)
  world <- read.csv(file = "World+human.csv", header = TRUE, stringsAsFactors = TRUE)
  
  head(world)
  names(world)
  
#Load global coral cover data 
  global <- read.csv("Global_Cover_1_5_10.csv")
  
#Load extra columns from coral data base
  extra <- read.csv("world_extra.csv")
  names(extra)
  
#Add columns to world data
  world <- cbind(world, extra[,c(2,3,5,6,7)])
  
#Explore world data with histograms and probabilities distributions
  par(mfcol=c(2,2))
  hist(world$MACROALGAE, prob=T, col="grey");lines(density(world$MACROALGAE, adjust=2),col=4, lwd=2)
  hist(world$HARD_CORAL, prob=T, col="grey");lines(density(world$HARD_CORAL, adjust=2),col=4, lwd=2)
  hist(world$PEOPLE,prob=T, col="grey");lines(density(world$PEOPLE, adjust=2),col=4, lwd=2) #this population variable needs to be calculated again
  hist(log(world$PEOPLE+1),prob=T, col="grey");lines(density(log(world$PEOPLE+1), adjust=2),col=4, lwd=2) #this population variable needs to be calculated again
  
#Transform Coral cover with logit anf run hist
  install.packages ('car')
  library(car)
  HARD_CORAL_logit <- logit(world$HARD_CORAL*0.01)
  hist(HARD_CORAL_logit, prob=T, col="grey", breaks=100); lines(density(HARD_CORAL_logit, adjust=2),col=4, lwd=2)
  
#### Get Human population data from SEDAC website ####
# Download data from the SEDAC webpage ans save it in subfolder in wd

# Install potential packages
  install.packages("sp")
    library(sp)
  install.packages("dismo")
    library(dismo)
  install.packages("rgdal")
    library(rgdal)
  install.packages("RArcInfo") # to import ArcGIS data
    library(RArcInfo)
  install.packages("gtools")
    library(gtools)
  install.packages("rworldmap")
      library(rworldmap)
  install.packages('mapdata')
    library(mapdata)
 
#Load raster layer from subfolder
#Load human population data from year 2005
  Pop <- raster(paste0(getwd(), "/Human Population Data/gl_gpwfe_pcount_05_bil_25/glp05ag.bil", 
                    sep = ""))

#Look at the info
  Pop 
  
#Plot data
  plot(Pop)
  
#Reclassify raster and set 0 to NA
  Pop.c <- reclassify(Pop, cbind(0, NA))
  
#Plot again to see results
  plot(Pop.c)

# Define spatial geographical projections,datum WGS84
  crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  
    
# Project for human pop data
  proj4string(Pop) <- crs.geo

# Get World sites with coordinates
  #colum 1 site names, colunm 8 and 9 lat and lon
  w.coord <- world[,c(1,8,9)]
  head(w.coord)
  
# Define coordinates columns
  coordinates(w.coord)<- c("LON","LAT")

# Define spatial projections for sites
  proj4string(w.coord) <- crs.geo
  summary(w.coord)

# Plot Sites on Map (TEST) #####
  
  ## Color code subregions
  

png("./Graphs/Fig 1test.png", 6.5,4,"in",res=600)
    par(oma=c(5,0,0,0))
    map("world", fill=T, col="grey90",  wrap=T, lty=1, lwd=0.4, 
        ylim=c(-90,90), mar=c(0,0,0,0), boundary=T, resolution=0)
    plot(w.coord, pch=16, cex=0.4, col=as.numeric(world$Subregion), add=T)
    legend(-180,-90, pch=16, col=rainbow(21, alpha=0.75), legend=subregions, 
           border="white", title = "SUBREGIONS", bty="n", xpd=NA, ncol=5, cex=0.7)
dev.off()
  
# Load data with coasts and countries
  data("coastsCoarse")
  data("countriesLow")
  data("world2HiresMapEnv")
  
# Now plot sites on map
# Set projections as the rest
  proj4string(countriesLow) <- crs.geo 
  proj4string(coastsCoarse) <- crs.geo
  plot(coastsCoarse, add=T, col="grey50")
  plot(countriesLow, add=T, col="grey95", lwd=0.1)
  
# Extract human population density for the world within 100 km buffer
  world$human100km <- extract(Pop, w.coord, buffer = 100000, fun = sum)

# Extract human population density for the world within 50 km buffer 
  world$human50km <- extract(Pop, w.coord, buffer = 50000, fun = sum)
  
# Extract human population density for the world within 25 km buffer 
  world$human25km <- extract(Pop, w.coord, buffer = 25000, fun = sum)

  write.csv(world, "world+human.csv")
  
# Check for correlation between human pop count within 100 km and within 50 km

  plot(human100km~human50km, world)

#### Check for spatial autocorrelation of raw data with spline correlograms ####

install.packages("ncf")
library(ncf)

#Spatial Autocorrelation for raw data of corals

 CORAL.cor <- spline.correlog(world$LON, world$LAT, world$HARD_CORAL, latlon=TRUE, resamp=1000)
    plot(CORAL.cor)

# Calculate autocorrelation by Ocean basin for corals
    #Atlantic Ocean
    Atl <- subset(world, OCEAN=="Atlantic")
    Coral.Atl.cor <- spline.correlog(Atl$LON, Atl$LAT, Atl$HARD_CORAL, latlon=TRUE, resamp=100)
    plot(Coral.Atl.cor)

    #Indian Ocean    
    Ind <- subset(world, OCEAN=="Indian")
    Coral.Ind.cor <- spline.correlog(Ind$LON, Ind$LAT, Ind$HARD_CORAL, latlon=TRUE, resamp=100)
    plot(Coral.Ind.cor)    

    #Pacific Ocean
    Pac <- subset(world, OCEAN=="Pacific")
    Coral.Pac.cor <- spline.correlog(Pac$LON, Pac$LAT, Pac$HARD_CORAL, latlon=TRUE, resamp=100)
    plot(Coral.Pac.cor)

#Spatial autocorrelation of raw data for macroalgae

  ALGAE.cor <- spline.correlog(world$LON, world$LAT, world$MACROALGAE, latlon=TRUE, resamp=100)
    plot(ALGAE.cor)

# Calculate autocorrelation by Ocean basin for Algae
    
  #Atlantic Ocean
    Atl <- subset(world, OCEAN=="Atlantic")
    Algae.Atl.cor <- spline.correlog(Atl$LON, Atl$LAT, Atl$MACROALGAE, latlon=TRUE, resamp=100)
    plot(Algae.Atl.cor)

   #Indian Ocean    
    Ind <- subset(world, OCEAN=="Indian")
    Algae.Ind.cor <- spline.correlog(Ind$LON, Ind$LAT, Ind$MACROALGAE, latlon=TRUE, resamp=100)
    plot(Algae.Ind.cor)    

  #Pacific Ocean
    Pac <- subset(world, OCEAN=="Pacific")
    Algae.Pac.cor <- spline.correlog(Pac$LON, Pac$LAT, Pac$MACROALGAE, latlon=TRUE, resamp=100)
    plot(Algae.Pac.cor)
    
    
#### Exploratory Models with human population data  ####

# Basic linear model for the world for HARD_CORAL cover and 3 pop densities

  modC1 <- lm(HARD_CORAL ~ human100km, world)
    summary(modC1)
    par(mfcol=c(2,2));plot(modC1)
  modC2<- lm(HARD_CORAL~ human50km, world) 
    summary(modC2);plot(modC2)
  modC3<- lm(HARD_CORAL~ human25km, world)
    summary(modC3); plot(modC1)
    
#the model with human pop within 50km has a slightly lower AIC
   AIC(modC1,modC2, modC3)
  
# Compare models
  anova(modC1,modC2, modC3)
  
##Try log of humans and compare AIC
  modC1log <- lm(HARD_CORAL ~ log(human100km+1), world)
    summary(modC1log);plot(modC1log)
  modC2log<- lm(HARD_CORAL~ log(human50km+1), world) 
    summary(modC2log);plot(modC2log)
  modC3log<- lm(HARD_CORAL~ log(human25km+1), world) 
    summary(modC3log);plot(modC3log)
  AIC(modC1,modC2,modC3,modC1log,modC2log,modC3log)

#Model corals by subregions and compare AIC
  #without log of human pop density
  modC4<- lm(HARD_CORAL~ human50km + SUBREGION, world) #the model with human pop within 50km has a lower AIC
  summary(modC4);
par(mfcol=c(2,2))
  plot(modC4)
  
  #with log of human pop density
  modC5<- lm(HARD_CORAL~ log(human50km+1) + SUBREGION, world) #the model with human pop within 50km has a lower AIC
  summary(modC4);
par(mfcol=c(2,2))
 plot(modC5)
  AIC(modC2,modC4, modC5)
  
# Model Corals by regions as random effects
  install.packages("lme4")
  library(lme4)
  library(car)
  
# Logit transform corals and rerun models 
  #Transform Hard coral cover with logit
  world$HARD_CORAL_logit <- logit(world$HARD_CORAL*0.01, adjust=0.025)
  hist(world$HARD_CORAL_logit)

  world$MACROALGAE_logit <- logit(world$MACROALGAE*0.01)
  hist(world$MACROALGAE_logit)

# Compare histograms for corals and algae
  png("./Logit_Coral_Algae_hist.png", 7,7,"in", res=600)
  par(mfcol=c(3,2))
  hist(world$HARD_CORAL,breaks=100,col="grey", prob=T, main="Hard Corals")
    lines(density(world$HARD_CORAL, adjust=2), col=4)
  hist(log(world$HARD_CORAL+1), breaks=100, col="grey", prob=T, main= "Log of Hard Corals")
    lines(density(log(world$HARD_CORAL+1), adjust=2), col=4)
  hist(world$HARD_CORAL_logit, breaks=100, col="grey", prob=T, main= "Logit of Hard Corals")
    lines(density(world$HARD_CORAL_logit, adjust=2), col=4)
  hist(world$MACROALGAE,breaks=100,col="grey", prob=T, main="Macroalgae")
    lines(density(world$MACROALGAE, adjust=2), col=4)
  hist(log(world$MACROALGAE+1), breaks=100, col="grey", prob=T, main= "Log of Macroalgae")
    lines(density(log(world$MACROALGAE+1), adjust=2), col=4)
  hist(world$MACROALGAE_logit, breaks=100, col="grey", prob=T, main= "Logit of Macroalgae")
    lines(density(world$MACROALGAE_logit, adjust=2), col=4)
  dev.off() 
  
### GLMM Models ####
  
# Run the models with logit transformation
# Cover should be in probabilities so divide by 100
# But first run null model (i.e no predictor) but with random effect
  library(lme4)  
  modC6null <-glmer(HARD_CORAL*0.01~1+(1|SUBREGION), data = world, 
                       family = binomial("logit"), nAGQ=2)
  AIC(modC6null)
  summary(modC6null)
  plot(modC6null)

# Run model with humans
  modC6 <- glmer(HARD_CORAL*0.01~scale(human50km)+(1|SUBREGION), data = world, 
                 family = binomial("logit"), nAGQ=2)
  summary(modC6)
  plot(modC6)
  AIC(modC6null, modC6)

#Check again if using log of human gives us a better fit
  modC7 <- glmer(HARD_CORAL*0.01~scale(log(human50km+1))+(1|SUBREGION), data = world, 
                 family = binomial("logit"), nAGQ=2)
  summary(modC7)
  plot(modC7)
  
  par(mfcol=c(2,2))
  plot(modC6); plot(modC7)
  AIC(modC6null, modC6, modC7)
  
#Compare models with anova 
  anova(modC6null, modC6, modC7)

#Try again the model with human within 100km 
  modC8 <- update(modC6,.~. +scale(human100km)-scale(human50km))
  summary(modC8)
  AIC(modC6null, modC6, modC7, modC8)
  
#Still model modc6 is the best based on AIC  
  
  #Check for spatial autocorrelation of modc6 residuals
  modC6.autocorr <- spline.correlog(world$LON, world$LAT, resid(modC6), latlon=TRUE, resamp=100)
  plot(modC6.autocorr)  

  #Compare Spatial Autocorrelation of the model with the raw data
   par(mfcol=c(2,1), mar=c(2,5,3,2))
   plot(CORAL.cor)
   mtext("raw data", 3, 1)
   plot(modC6.autocorr)
   mtext("model residuals", 3, 1)
   
   #there is still spatial autocorrelation problems
   world$REGION <- world$REGION_PS

   #Try Nested model with all combinations of REGION, SUBREGION, OCEAN and COUNTRY
   modC9 <- glmer(HARD_CORAL*0.01~scale(human50km)+ (1|OCEAN/REGION/SUBREGION/COUNTRY), data = world, 
                  family = binomial("logit"), nAGQ=1)
   summary(modC9)
   plot(modC9)
   
   AIC(modC6null,modC6,modC9)
   anova(modC6null,modC6,modC9)

#Adding region or ocean or subregion as random effect does not improve model
   
#Check for spatial autocorrelation of modc9 residuals
   modC9.autocorr <- spline.correlog(world$LON, world$LAT, resid(modC9), latlon=TRUE, resamp=10)
   
#Compare Spatial Autocorrelation of the model with the raw data
   par(mfcol=c(3,1), mar=c(2,5,1,2))
   plot(CORAL.cor)
   mtext("raw data", 3, 1)
   plot(modC6.autocorr)
   mtext("model residuals", 3, 1)
   plot(modC9.autocorr)

# Run model by Ocean basin

  modCAtl <- glmer(HARD_CORAL*0.01~scale(log(human50km+1))+(1|SUBREGION), data = subset(world, OCEAN=="Atlantic"), 
               family = binomial("logit"), nAGQ=2)
  modCInd <- glmer(HARD_CORAL*0.01~scale(log(human50km+1))+(1|SUBREGION), data = subset(world, OCEAN=="Indian"), 
                 family = binomial("logit"), nAGQ=2)
  modCPac <- glmer(HARD_CORAL*0.01~scale(log(human50km+1))+(1|SUBREGION), data = subset(world, OCEAN=="Pacific"), 
                 family = binomial("logit"), nAGQ=2)

  summary(modCAtl)
  summary(modCInd)
  summary(modCPac)
    AIC(modCAtl)
    AIC(modCInd)
    AIC(modCPac)

## Check model residuals

  plot(modCAtl)
  plot(modCInd)
  plot(modCPac)

# Check for spatial autocorrelation of model residuals
  
  modCAtl.autocorr <- spline.correlog(Atl$LON, Atl$LAT, resid(modCAtl), latlon=TRUE, resamp=100)
  modCInd.autocorr <- spline.correlog(Ind$LON, Ind$LAT, resid(modCInd), latlon=TRUE, resamp=100)
  modCPac.autocorr <- spline.correlog(Pac$LON, Pac$LAT, resid(modCPac), latlon=TRUE, resamp=100)

  plot(modCAtl.autocorr)
  plot(modCInd.autocorr)
  plot(modCPac.autocorr)

png("./spatial_correlation.png", 5,7,"in",res=600)
  par(mfcol=c(3,2), mar=c(2,5,2,2), cex=0.8)
    plot(Coral.Atl.cor, ylim=c(-0.5,0.5))
  mtext("Atlantic Basin Raw Coral data",3,0.5, cex=0.8)
    plot(Coral.Ind.cor, ylim=c(-0.5,0.5))
  mtext("Indian Basin Raw Coral data",3,0.5, cex=0.8)
    plot(Coral.Pac.cor, ylim=c(-0.6,0.6))
  mtext("Pacific Basin Raw Coral data",3,0.5, cex=0.8)
    plot(modCAtl.autocorr, ylim=c(-0.6,0.6))
  mtext("Atlantic model residuals",3,0.5,cex=0.8)
    plot(modCInd.autocorr, ylim=c(-0.5,0.5))
  mtext("Indian model residuals",3,0.5, cex=0.8)
    plot(modCPac.autocorr, ylim=c(-0.5,0.5))
  mtext("Pacific model residuals",3,0.5, cex=0.8)
dev.off()

### GAMM MODELS ####

install.packages('gamm4')
library(gamm4)

install.packages('mgcv')
library(mgcv)

# For Hard Corals at global analysis

  modC10 <- gamm4(HARD_CORAL*0.01~s(log(human50km+1)), random=~(1|SUBREGION), data = world, 
               family = binomial("logit"))

  summary(modC10$mer)
  summary(modC10$gam)
  plot(modC10$mer)
  plot(modC10$gam, residuals=T, cex=2)
  
# Check for spatial autocorrelation 

  modC10.corr <- spline.correlog(world$LON, world$LAT, resid(modC10$mer), latlon=TRUE, resamp=100)
  plot(modC10.corr)
  
# GAMM by Ocean basin

  modC10Atl <- gamm4(HARD_CORAL*0.01~s(log(human50km+1)), random=~(1|SUBREGION), data = subset(world, OCEAN=="Atlantic"), 
                 family = binomial("logit"))
  modC10Ind <- gamm4(HARD_CORAL*0.01~s(log(human50km+1)), random=~(1|SUBREGION), data = subset(world, OCEAN=="Indian"), 
                 family = binomial("logit"))
  modC10Pac <- gamm4(HARD_CORAL*0.01~s(log(human50km+1)), random=~(1|SUBREGION), data = subset(world, OCEAN=="Pacific"), 
                 family = binomial("logit"))
  
  summary(modC10Atl$mer)
  summary(modC10Ind$mer)
  summary(modC10Pac$mer)

  summary(modC10Atl$gam)
  summary(modC10Ind$gam)
  summary(modC10Pac$gam)
  
  plot(modC10Atl$mer)
  plot(modC10Ind$mer)
  plot(modC10Pac$mer)

  plot(modC10Atl$gam)
  plot(modC10Ind$gam)
  plot(modC10Pac$gam)

# Run model with gamm to add correlation structure (gamm4 does not allow correlation structure)

  modC11 <- gamm(HARD_CORAL*0.01~s(log(human50km+1)), random=list(SUBREGION=~1), data = world, 
                family = binomial("logit"))
  
  summary(modC11$lme)
  summary(modC11$gam)
  plot(modC11$lme)
  plot(modC11$gam)

#Check Spatial autocorrelation

  modC11.corr <- spline.correlog(world$LON, world$LAT, resid(modC11$lme), latlon=TRUE, resamp=100)
  plot(modC11.corr)

#Add correlation structure

#Global

  modC11c <- gamm(HARD_CORAL*0.01~s(log(human50km+1), k=100), random=list(OCEAN=~1, SUBREGION=~1), data = world, 
             family = binomial('logit'), correlation=corAR1())

  summary(modC11c$lme)
  summary(modC11c$gam)
  plot(modC11c$lme)
  plot(modC11c$gam, residuals=T)
  
#With Latitude  as variable
  modC11ct <- gamm(HARD_CORAL*0.01~s(log(human50km+1), k=100),  random=list(OCEAN=~1, SUBREGION=~1), 
                   data = subset (world), 
                  family = binomial('logit'), correlation=corAR1(), gamma=1.4)
  
  summary(modC11ct$lme)
  summary(modC11ct$gam)
  plot(modC11ct$lme)
  plot(modC11ct$gam, residuals=T)

#GAM by subregions for Corals 

  modC11cANI <- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Andaman_Nicobar"),
                  family = binomial('logit'), correlation=corAR1())
  summary(modC11cANI)
  plot(modC11cANI)

  modC11cANT <- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Antilles"),
                  family = binomial('logit'), correlation=corAR1())
  summary(modC11cANT)
  plot(modC11cANT)

  modC11cBAH <- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Bahamas"),
                  family = binomial('logit'), correlation=corAR1(), gamma=1.4)
  summary(modC11cBAH)
  plot(modC11cBAH)

  modC11cBRZ <- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Brazil"),
                  family = binomial('logit'), correlation=corAR1())
  summary(modC11cBRZ)
  plot(modC11cBRZ)

  modC11cCAM <- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="C. America"),
                  family = binomial('logit'), correlation=corAR1())
  summary(modC11cCAM)
    plot(modC11cCAM)

  modC11cEIP <- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="E Indonesia_PNG"),
                  family = binomial('logit'), correlation=corAR1())
  summary(modC11cEIP)
  plot(modC11cEIP)

  modC11cECA <- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="E. Caribbean"),
                  family = binomial('logit'), correlation=corAR1())
  summary(modC11cECA)
  plot(modC11cECA)

  modC11cEIO <- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="E. Indian Ocean"),
                  family = binomial('logit'), correlation=corAR1())
  summary(modC11cEIO)
  plot(modC11cEIO)

  modC11cFKY<- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Florida Keys"),
                  family = binomial('logit'), correlation=corAR1())
  summary(modC11cFKY)
  plot(modC11cFKY)

  modC11cGBR <- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="GBR"),
                family = binomial('logit'), correlation=corAR1())
  summary(modC11cGBR)
      plot(modC11cGBR)

  modC11cHAW <- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Hawaiian Islands"),
                  family = binomial('logit'), correlation=corAR1())
  summary(modC11cHAW)
  plot(modC11cHAW)

  modC11cMBR <- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Mesoamerican BR"),
                  family = binomial('logit'), correlation=corAR1())
  summary(modC11cMBR)
  plot(modC11cMBR)

  modC11cPHI <- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Philippines_Spra"),
                  family = binomial('logit'), correlation=corAR1())
  summary(modC11cPHI)
  plot(modC11cPHI)

  modC11cSEP <- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="SE Pacific"),
                  family = binomial('logit'), correlation=corAR1())
  summary(modC11cSEP)
  plot(modC11cSEP)

    modC11cSCS <- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="South China Sea"),
                    family = binomial('logit'), correlation=corAR1())
    summary(modC11cSCS)
    plot(modC11cSCS)

  modC11cSWP<- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="SW Pacific"),
                  family = binomial('logit'), correlation=corAR1())
  summary(modC11cSWP)
  plot(modC11cSWP)

  modC11cTWJ<- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Taiwan_Japan"),
                 family = binomial('logit'), correlation=corAR1())
  summary(modC11cTWJ)
  plot(modC11cTWJ)

  modC11cWIN<- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="W Indonesia"),
                 family = binomial('logit'), correlation=corAR1())
  summary(modC11cWIN)
  plot(modC11cWIN)

  modC11cWEA<- glm(HARD_CORAL*0.01~log(human50km+1), data = subset(world, SUBREGION=="W. Australia"),
                 family = binomial('logit'))
  summary(modC11cWEA)
  plot(modC11cWEA)

  modC11cWEP<- gam(HARD_CORAL*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Western Pacific"),
                   family = binomial('logit'), correlation=corAR1())
  summary(modC11cWEP)
  plot(modC11cWEP)


#GAM by subregions for Algae 

modM11cANI <- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Andaman_Nicobar"),
                  family = binomial('logit'), correlation=corAR1())
summary(modM11cANI)
plot(modM11cANI)

modM11cANT <- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Antilles"),
                  family = binomial('logit'), correlation=corAR1())
summary(modM11cANT)
plot(modM11cANT)

modM11cBAH <- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Bahamas"),
                  family = binomial('logit'), correlation=corAR1())
summary(modM11cBAH)
plot(modM11cBAH)

modM11cBRZ <- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Brazil"),
                  family = binomial('logit'), correlation=corAR1())
summary(modM11cBRZ)
plot(modM11cBRZ)

modM11cCAM <- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="C. America"),
                  family = binomial('logit'), correlation=corAR1())
summary(modM11cCAM)
plot(modM11cCAM)

modM11cEIP <- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="E Indonesia_PNG"),
                  family = binomial('logit'), correlation=corAR1())
summary(modM11cEIP)
plot(modM11cEIP)

modM11cECA <- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="E. Caribbean"),
                  family = binomial('logit'), correlation=corAR1())
summary(modM11cECA)
plot(modM11cECA)

modM11cEIO <- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="E. Indian Ocean"),
                  family = binomial('logit'), correlation=corAR1())
summary(modM11cEIO)
plot(modM11cEIO)
AIC(modM11cEIO)

modM11cFKY<- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Florida Keys"),
                 family = binomial('logit'), correlation=corAR1())
summary(modM11cFKY)
plot(modM11cFKY)

modM11cGBR <- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="GBR"),
                  family = binomial('logit'), correlation=corAR1())
summary(modM11cGBR)
plot(modM11cGBR)

modM11cHAW <- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Hawaiian Islands"),
                  family = binomial('logit'), correlation=corAR1())
summary(modM11cHAW)
plot(modM11cHAW)

modM11cMBR <- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Mesoamerican BR"),
                  family = binomial('logit'), correlation=corAR1())
summary(modM11cMBR)
plot(modM11cMBR)

modM11cPHI <- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Philippines_Spra"),
                  family = binomial('logit'), correlation=corAR1())
summary(modM11cPHI)
plot(modM11cPHI)

modM11cSEP <- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="SE Pacific"),
                  family = binomial('logit'), correlation=corAR1())
summary(modM11cSEP)
plot(modM11cSEP)

modM11cSCS <- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="South China Sea"),
                  family = binomial('logit'), correlation=corAR1())
summary(modM11cSCS)
plot(modM11cSCS)

modM11cSWP<- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="SW Pacific"),
                 family = binomial('logit'), correlation=corAR1())
summary(modM11cSWP)
plot(modM11cSWP)

modM11cTWJ<- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Taiwan_Japan"),
                 family = binomial('logit'), correlation=corAR1())
summary(modM11cTWJ)
plot(modM11cTWJ)

modM11cWIN<- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="W Indonesia"),
                 family = binomial('logit'), correlation=corAR1())
summary(modM11cWIN)
plot(modM11cWIN)

modM11cWEP<- gam(MACROALGAE*0.01~s(log(human50km+1)), data = subset(world, SUBREGION=="Western Pacific"),
                 family = binomial('logit'), correlation=corAR1())
summary(modM11cWEP)
plot(modM11cWEP)


#Check Spatial autocorrelation

  acf(residuals(modC11c$lme))

  modC11c.corr <- spline.correlog(world$LON, world$LAT, resid(modC11c$lme), latlon=TRUE, resamp=100)
  plot(modC11c.corr, ylim=c(-0.5,0.5))
  
# Coral GAMM wtih gamm by Ocean basin

  modC11Atl <- gamm(HARD_CORAL*0.01~s(log(human50km+1)), random=list(SUBREGION=~1), data = subset(world, OCEAN=="Atlantic"), 
                  correlation = corAR1(), family = binomial("logit"))

  modC11Ind <- gamm(HARD_CORAL*0.01~s(log(human50km+1)), random=list(SUBREGION=~1), data = subset(world, OCEAN=="Indian"), 
                  correlation = corAR1(), family = binomial("logit"))

  modC11Pac <- gamm(HARD_CORAL*0.01~s(log(human50km+1)), random=list(SUBREGION=~1), data = subset(world, OCEAN=="Pacific"), 
                  correlation = corAR1(), family = binomial("logit"))

  summary(modC11Atl$lme)
  summary(modC11Ind$lme)
  summary(modC11Pac$lme)

  summary(modC11Atl$gam)
  summary(modC11Ind$gam)
  summary(modC11Pac$gam)

  plot(modC11Atl$lme)
  plot(modC11Ind$lme)
  plot(modC11Pac$lme)

  plot(modC11Atl$gam, residuals=T)
  plot(modC11Ind$gam, residuals=T)
  plot(modC11Pac$gam, residuals=T)

#check for Spatial autocorrelation

  modC11Atl.corr <- spline.correlog(world$LON, world$LAT, resid(modC11Atl$lme), latlon=TRUE, resamp=100)
  plot(modC11Atl.corr)

  modC11Ind.corr <- spline.correlog(world$LON, world$LAT, resid(modC11Ind$lme), latlon=TRUE, resamp=100)
  plot(modC11Ind.corr)

  modC11Pac.corr <- spline.correlog(world$LON, world$LAT, resid(modC11Pac$lme), latlon=TRUE, resamp=100)
  plot(modC11Pac.corr)
  
  
### FIGURE S1. Spatial autocorrelation spline correlograms for Coral and Algae
  
  png("C:/Users/avaldivia/Dropbox/Coral threats-reef isolation paper/Graphs/Fig S1.png", 6,6,"in", res=600)
  
    par(mfcol=c(2,2), mar=c(2,2,3,1), oma=c(3,3,0,0), cex=0.8)
    
      plot(CORAL.cor, ylim=c(-0.5,0.5))
        mtext("Moran's I Correlation",2,3, cex=0.8)
        mtext("Global raw data coral cover",3,1, cex=0.8)
        mtext("A", 3, 1, adj=-0.3, font=2)
      plot(modC11c.corr, ylim=c(-0.5,0.5))
        mtext("Moran's I Correlation",2,3, cex=0.8)
        mtext("Distance (km)",1,3, cex=0.8)
        mtext("Model residuals for coral cover",3,1, cex=0.8)
        mtext("C", 3, 1, adj=-0.3, font=2)
      plot(ALGAE.cor, ylim=c(-0.5,0.5))
        mtext("Global raw data macroalgae cover",3,1, cex=0.8)
        mtext("B", 3, 1, adj=-0.2, font=2)
      plot(modM11c.corr, ylim=c(-0.5,0.5))
        mtext("Model residuals for macroalgae cover",3,1, cex=0.8)
        mtext("Distance (km)",1,3, cex=0.8)
        mtext("D", 3, 1, adj=-0.2, font=2)
  
    dev.off()
  
# Run models for Macroalgae with the same correlation structure
   
#Global
  #with gamm4
  modM10 <- gamm4(MACROALGAE*0.01~s(log(human50km+1)), random=~(1|SUBREGION), data = world, 
                  family = binomial("logit"))
  plot(modM10$gam, residuals=T, cex=2)
  
  #with gamm
  modM11c <- gamm(MACROALGAE*0.01~s(log(human50km+1)), random=list(OCEAN=~1, SUBREGION=~1), data = world, 
                family = binomial("logit"), correlation=corAR1())

  summary(modM11c$lme)
  summary(modM11c$gam)
  plot(modM11c$lme)
  plot(modM11c$gam)

#By Basin
  
  modM11Atl <- gamm(MACROALGAE*0.01~s(log(human50km+1)), random=list(SUBREGION=~1), data = subset(world, OCEAN=="Atlantic"), 
                  correlation = corAR1(), family = binomial("logit"))

  modM11Ind <- gamm(MACROALGAE*0.01~s(log(human50km+1)), random=list(SUBREGION=~1), data = subset(world, OCEAN=="Indian"), 
                  correlation = corAR1(), family = binomial("logit"))

  modM11Pac <- gamm(MACROALGAE*0.01~s(log(human50km+1)), random=list(SUBREGION=~1), data = subset(world, OCEAN=="Pacific"), 
                  correlation = corAR1(), family = binomial("logit"))

  summary(modM11Atl$lme)
  summary(modM11Ind$lme)
  summary(modM11Pac$lme)

  summary(modM11Atl$gam)
  summary(modM11Ind$gam)
  summary(modM11Pac$gam)

  plot(modM11Atl$lme)
  plot(modM11Ind$lme)
  plot(modM11Pac$lme)

  plot(modM11Atl$gam, residuals=T)
  plot(modM11Ind$gam, residuals=T)
  plot(modM11Pac$gam, residuals=T)

#check for Spatial autocorrelation

  modM11c.corr <- spline.correlog(world$LON, world$LAT, resid(modM11c$lme), latlon=TRUE, resamp=100)
  plot(modM11c.corr)

  modM11Atl.corr <- spline.correlog(world$LON, world$LAT, resid(modM11Atl$lme), latlon=TRUE, resamp=100)
  plot(modM11Atl.corr)

  modM11Ind.corr <- spline.correlog(world$LON, world$LAT, resid(modM11Ind$lme), latlon=TRUE, resamp=100)
  plot(modM11Ind.corr)

  modM11Pac.corr <- spline.correlog(world$LON, world$LAT, resid(modM11Pac$lme), latlon=TRUE, resamp=100)
  plot(modM11Pac.corr)
 

# Build figure of coral, algae, vs human pop density in ggplot ####

##Figure 1 ####
  
 ## Color code subregions
  
  subregions=c("Andaman / Nicobar","Antilles","Bahamas", "Brazil", "Central America", "E. Indonesia",
               "E. Caribbean", "E. Indian Ocean", "Florida Keys", "Great Barrier Reef",
               "Hawaiian Islands", "Mesoamerican BR", "N. Central Australia", "Philippines",
               "SE. Pacific", "S. China Sea", "SW. Pacific", "Taiwan / Japan", "W. Indonesia",
               "W. Australia", "W. Pacific")
  
  subregions_colors=rainbow(21)
 
# Recode subregion colors
  
  world$Subregion <- world$SUBREGION
  world$SubregionColor <- as.character(recode(world$SUBREGION, 
         "'Andaman_Nicobar' = '#FF0000BF';
         'Antilles' = '#FF4900BF';
         'Bahamas' = '#FF9200BF';
         'Brazil' = '#FFDB00BF';
         'C. America' = '#DBFF00BF';
         'E Indonesia_PNG' = '#92FF00BF';
         'E. Caribbean' = '#49FF00BF';
         'E. Indian Ocean' = '#00FF00BF';
         'Florida Keys' = '#00FF49BF';
         'GBR' = '#00FF92BF';
         'Hawaiian Islands' = '#00FFDBBF';
         'Mesoamerican BR' = '#00DBFFBF';
         'NCentrral Austra' = '#0092FFBF';
         'Philippines_Spra' = '#0049FFBF';
         'SE Pacific' = '#0000FFBF';
         'South China Sea' = '#4900FFBF';
         'SW Pacific' = '#9200FFBF';
         'Taiwan_Japan' = '#DB00FFBF';
         'W Indonesia' = '#FF00DBBF';
         'W. Australia' = '#FF0092BF';
         'Western Pacific' = '#FF0049BF' "))

  world$SubregionNames <- (recode(world$SUBREGION, 
        "'Andaman_Nicobar' = 'Andaman/Nicobar';
         'E Indonesia_PNG' = 'E. Indonesia';
         'GBR' = 'Great Barrier Reef';
         'Philippines_Spra' = 'Philippines';
         'SE Pacific' = 'SE. Pacific';
         'South China Sea' = 'S. China Sea';
         'SW Pacific' = 'SW. Pacific';
         'Taiwan_Japan' = 'Taiwan/Japan';
         'W Indonesia' = 'W. Indonesia';
         'W. Australia' = 'W. Australia';
         'Western Pacific' = 'W. Pacific'"))
  
  png("./Graphs/Fig 1.png", 6.5,4,"in",res=600)
  par(oma=c(5,0,0,0))
  map("world", fill=T, col="grey90",  wrap=T, lty=1, lwd=0.4, 
      ylim=c(-90,90), mar=c(0,0,0,0), boundary=T, resolution=0)
  plot(w.coord, pch=16, cex=0.4, col=world$SubregionColor, add=T)
  legend(-180,-90, pch=16, col=rainbow(21, alpha=0.75), legend=subregions, 
         border="white", title = "SUBREGIONS", bty="n", xpd=NA, ncol=5, cex=0.7)
  dev.off()  
  

## Figure 2  
   
png("C:/Users/avaldivia/Dropbox/Coral threats-reef isolation paper/Graphs/Fig 2v2.png",8,7,"in",res=600, pointsize=10)
  
    par(mfrow=c(2,2), mar=c(2,2,1,1), oma=c(2,2,1,3))
    
    nf1 <- layout(matrix(c(1:12),2,6,byrow = TRUE), c(4,1.5,0.8,4,1.5,2), c(4,4), TRUE)
    layout.show(nf1)
      
    plot(HARD_CORAL~log(human50km+1), world, pch=16, cex=0.8, col= world$SubregionColor, ylab='')
      mtext("Coral cover (%)", 2, line=2.5)
      mtext("CORAL COVER", 3, line=1)
      mtext("A", 3, line =1, at=-4, font=2)
      #boxplot(humans, horizontal=T, notch=F, log="x")
      #boxplot(world$HARD_CORAL, vertical=T, notch=F)
    boxplot(world$HARD_CORAL)
    plot(1,1, col="white", axes=NULL)
      
    plot(MACROALGAE~log(human50km+1), world, pch=16, cex=0.8, col=world$SubregionColor, ylim=c(0,100), ylab='')
      mtext("Macroalgae cover (%)", 2, line=2.5)
      mtext("MACROALGAE COVER", 3, line=1)
      mtext("B", 3, line =1, at=-4, font=2)
      #boxplot(humans, horizontal=T, notch=F, log="x")
      #boxplot(world$MACROALGAE, vertical=T, notch=F)
    boxplot(world$MACROALGAE)
    plot(1,1, col="white", axes=NULL)
      
    legend("topleft", pch=16, col=rainbow(21, alpha=0.75), legend=subregions, 
           border="white", title = "SUBREGIONS", bty="n", xpd=NA)
      
    plot(modC11c$gam, ylim=c(-2.0,4.5), ylab='', residuals=T, cex=1.5, shade=T, shade.col="grey90")
      mtext("GAMM smooth (coral cover)", 2, line=2.5)
      mtext("C", 3, line =1, at=-4., font=2)
      
    plot(1,1, col=0, axes=NULL)
    plot(1,1, col=0, axes=NULL)
      
    plot(modM11c$gam, ylim=c(-2.0,4.5), ylab='', residuals=T, cex=1.5, shade=T, shade.col="grey90")
      mtext("GAMM smooth (macroalgae cover)", 2, line=2.5)
      mtext("Logarithm of human population density within 50 km", 1, line=2.5, at=-7)
      mtext("D", 3, line =1, at=-4, font=2)
      
    plot(1,1, col=0, axes=NULL)
    plot(1,1, col=0,axes=NULL)
    
  dev.off()        

## FIGURE S2 an S3. Subregional scatter plots ####

library(ggplot2)

#Build layout
  nf <- layout(matrix(c(1,2,3,4),2,2,byrow = TRUE), c(1,1), TRUE)
  nf1 <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(4,1), c(1,4), TRUE)
    layout.show(nf)
    layout.show(nf1)

# Exclude subregions with one location such as N central Australia

dat <- subset(world, SubregionNames != "NCentrral Austra")

  png('C:/Users/Abel Valdivia/Documents/My Dropbox/Coral threats-reef isolation paper/Graphs/Fig S2 Coral.png', 6.5,5.5, "in", res=600, pointsize=8) 

    p1 <- ggplot(dat, aes(log(human50km+1),HARD_CORAL))+ 
      facet_wrap(~SubregionNames)+ 
      geom_smooth(method="gam", fill="grey80")+
      geom_point(pch=16, cex=0.8)+ #scale_x_log10()+
      theme_bw()+theme(legend.title=element_blank())+
      theme(text = element_text(size=10))+
      xlab("Logarithm of human population density within 50 km of reef") + ylab("Coral cover (%)")
  p1

  dev.off()

  png('C:/Users/Abel Valdivia/Documents/My Dropbox/Coral threats-reef isolation paper/Graphs/Fig S3 Algae.png', 6.5,5.5, "in", res=600, pointsize=8)

  p2 <- ggplot(dat, aes(log(human50km+1), MACROALGAE))+ facet_wrap(~SubregionNames)+
    geom_smooth(method="gam",fill="grey80")+
    geom_point(pch=16, cex=0.8)+ ylim(c(-15,100))+
    theme_bw()+ theme(legend.title=element_blank())+
    theme(text = element_text(size=10))+
    xlab("Logarithm of human population density within 50 km of reef") + ylab("Macroalgae cover (%)")
  p2
  
  dev.off()
    
## Figure 3 Histograms of benthic cover with no people ###

  # Extract data with no people
    dat0Human <- subset(world,human50km==0)

  #Make figure

png('C:/Users/Abel Valdivia/Documents/My Dropbox/Coral threats-reef isolation paper/Graphs/Fig 3 Hist.png', 6,3.5, "in", res=600, pointsize=10)

  par(mar=c(2,3,1,1),oma=c(2,1,2,1))
  
  nf <- layout(matrix(c(1,2,3,4),2,2,byrow = TRUE), c(5,5),c(1.5,4), TRUE)
    layout.show(nf)
  
    boxplot(dat0Human$HARD_CORAL, horizontal=T, ylim=c(0,80))
      mtext("HARD CORALS", 3, line=1)
    mtext("A", 2, line=2, adj=1, font=2, las=2)
    boxplot(dat0Human$MACROALGAE, horizontal=T, ylim=c(0,80))
      mtext("MACROALGAE",3, line=1)
    mtext("B", 2, line=2, adj=1, font=2, las=2)
    hist(dat0Human$HARD_CORAL, breaks=25,col="grey", prob=T,  ylim=c(0,0.1), xlim=c(0,80), main=NULL)
      lines(density(dat0Human$HARD_CORAL, adjust=2), col=4)
      mtext("Density", 2, line=3)
    hist(dat0Human$MACROALGAE, breaks=25, col="grey", prob=T,  ylim=c(0,0.1), xlim=c(0,80), main=NULL)
      lines(density(dat0Human$MACROALGAE, adjust=2), col=4)
    mtext("Benthic Cover (%)", 1, line =2.5, adj=-1)

dev.off()


