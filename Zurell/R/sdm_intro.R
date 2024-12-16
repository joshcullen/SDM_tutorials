
### Code from Intro to SDM from Damaris Zurell ###
#https://damariszurell.github.io/SDM-Intro/

library(raster)
library(data.table)
library(randomForest)
library(lattice)
library(PresenceAbsence)
library(RColorBrewer)
library(geodata)
library(mecofun)


### Load data ###

avi_dat <- read.table('Zurell/data/Data_SwissBreedingBirds.csv', header=T, sep=',')

summary(avi_dat)



### Subset data for Ringed ouzel and relevant predictors ###

avi_cols <- c('Turdus_torquatus', 'bio_5', 'bio_2', 'bio_14', 'std', 'rad', 'blockCV_tile')

avi_df <- data.frame(avi_dat)[avi_cols]

summary(avi_df)



### Acquire new bioclimatic data ###

# Please note that you have to set download=T if you haven't downloaded the data before:
bio_curr <- worldclim_country(country = "Switzerland", var='bio', res=0.5, lon=5.5, lat=45.5,
                              path='Zurell/data', download=T)[[c(2,5,14)]]

# Please note that you have to set download=T if you haven't downloaded the data before:
# bio_fut <- getData('CMIP5', var='bio', res=0.5, lon=5.5, lat=45.5, rcp=45, model='NO', year=50, path='data', download=F)[[c(2,5,14)]]
bio_fut <- cmip6_tile(lon = 5, lat = 45, model = "ACCESS-CM2", ssp = 585, time = "2041-2060", var = "bioc",
                      res = 5, path = 'Zurell/data')




### Clip rasters to Switzerland boundaries ###

# A spatial mask of Switzerland in Swiss coordinates
bg <- rast('/vsicurl/https://damariszurell.github.io/SDM-Intro/CH_mask.tif')

# The spatial extent of Switzerland in Lon/Lat coordinates is roughly:
ch_ext <- c(5, 11, 45, 48)

# Crop the climate layers to the extent of Switzerland
bio_curr <- crop(bio_curr, ch_ext)
plot(bio_curr)

# Re-project to Swiss coordinate system and clip to Swiss political bounday
bio_curr <- project(bio_curr, bg)
bio_curr <- resample(bio_curr, bg)
bio_curr <- mask(bio_curr, bg)
names(bio_curr) <- c('bio_2', 'bio_5', 'bio_14')

# For storage reasons the temperature values in worldclim are multiplied by 10. For easier interpretability, we change it back to °C.
# bio_curr[[1]] <- bio_curr[[1]]/10
# bio_curr[[2]] <- bio_curr[[2]]/10

# Repeat above steps for future climate layers
# bio_fut <- crop(bio_fut, ch_ext)
# bio_fut <- projectRaster(bio_fut, bg)
# bio_fut <- resample(bio_fut, bg)
# bio_fut <- mask(bio_fut, bg)
# names(bio_fut) <- c('bio_2', 'bio_5', 'bio_14')
# bio_fut[[1]] <- bio_fut[[1]]/10
# bio_fut[[2]] <- bio_fut[[2]]/10



### Fit SDM via GLM ###

# Fit GLM
m_glm <- glm(Turdus_torquatus ~ bio_2 + I(bio_2^2) + bio_5 + I(bio_5^2) + bio_14 + I(bio_14^2),
              family='binomial', data=avi_df)

summary(m_glm)


## Inspect partial effects

# Names of our variables:
pred <- c('bio_2', 'bio_5', 'bio_14')

# We want three panels next to each other:
par(mfrow=c(1,3)) 

# Plot the partial responses
partial_response(m_glm, predictors = avi_df[,pred])



# We prepare the response surface by making a dummy data set where two predictor variables range from their minimum to maximum value, and the remaining predictor is kept constant at its mean:
xyz <- data.frame(expand.grid(seq(min(avi_df[,pred[1]]),max(avi_df[,pred[1]]),length=50),
                              seq(min(avi_df[,pred[2]]),max(avi_df[,pred[2]]),length=50)),
                  mean(avi_df[,pred[3]]))
names(xyz) <- pred

# Make predictions
xyz$z <- predict(m_glm, xyz, type='response')
summary(xyz)


# Make a colour scale
cls <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(100)

# plot 3D-surface
wireframe(z ~ bio_2 + bio_5, data = xyz, zlab = list("Occurrence prob.", rot=90),
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), zlim = c(0, 1), 
          main='GLM', xlab='bio_2', ylab='bio_5', screen=list(z = 120, x = -70, y = 3))


# Plot inflated response curves:
par(mfrow=c(1,3))
inflated_response(m_glm, predictors = avi_df[,pred], method = "stat6", lwd = 3, main='GLM') 
