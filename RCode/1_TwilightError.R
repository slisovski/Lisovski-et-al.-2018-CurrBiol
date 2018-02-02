#############################################
#### Lisovski et al. 2018 Current Biology ###
#############################################

library(TwGeos)
library(SGAT)
library(geosphere)
library(GeoLight)
library(MASS)

lon.breed <- -84.29
lat.breed <-  36.27


threshold <- 1
offset    <- 16

## Individual ID and colors
names <- c("GCM05", "GCM06", "GCM09", "GCM13", "GCM16") 
cols  <- c("yellow", "orange", "darkgreen", "cornflowerblue", "purple")


## Calibration function
calib <- function(twl_c, lat, zenith.start = 96) {
  z0    <- seq(zenith.start-10, zenith.start+10, by = 0.2)
  crds1 <- lapply(cbind(z0), function(x) thresholdPath(twl_c$Twilight, twl_c$Rise, zenith = x)$x)
  dist1 <- unlist(lapply(crds1, function(x1) median(abs(x1[,2]-lat))))
  z0[which.min(dist1)]
}


## The raw ligth recordings were downloaded from https://conservancy.umn.edu/handle/11299/183086?show=full and safed in folder LightData
## See also:	Kramer, Gunnar, R; Streby, Henry M; Peterson, Sean M; Lehman, Justin A; Buehler, David A; Wood, 
## Petra B; McNeil, Darin J; Larkin, Jeffrey L; Andersen, David E. (2016). Raw Light-Level Geolocator Data from 
## Golden-Winged Warblers Breeding at Three Sites in North America. Retrieved from the Data Repository for the University of Minnesota



###################################
## Calculation of Twiligth Error ##
###################################

tm <- seq(as.POSIXct("2013-05-01", tz = "GMT"), as.POSIXct("2014-05-10", tz = "GMT"), by = "day")
rise <- rep(c(TRUE, FALSE), length(tm))
c.dat <- data.frame(Twilight = twilight(rep(tm, each = 2), lon = lon.breed, lat = lat.breed, 
                                        rise = rise, zenith = 93), Rise = rise)



twilight.dev <- list()
crds.breed   <- list()


### GCM05
GCM05 <- readLig("LightData/GCM052014.lig")
GCM05 <- GCM05[-1,]

    lightImage(GCM05, offset = offset)
    tsimagePoints(c.dat$Twilight, offset = offset, pch = 16, cex = 0.25,
                  col = "orange")

twl <- read.csv("LightData/twl_GCM05.csv")
twl$Twilight <- as.POSIXct(twl$Twilight, tz = "GMT")
twl <- twilightAdjust(twl, 120)

calib.tm  <- c(as.POSIXct("2013-05-10", tz = "GMT"), as.POSIXct("2013-06-10", tz = "GMT"))  
twl_calib <- subset(twl, Twilight>=calib.tm[1] & Twilight<=calib.tm[2])

### get zenith angle and alpha parameters
sun  <- solar(twl_calib[,1])
z    <- refracted(zenith(sun, lon.breed, lat.breed))
# plot(z)

twl_t   <- twilight(twl_calib[,1], lon.breed, lat.breed, rise = twl_calib[,2], zenith = max(z)+0.1)
twl_dev <- ifelse(twl_calib$Rise, as.numeric(difftime(twl_calib[,1], twl_t, units = "mins")),
                  as.numeric(difftime(twl_t, twl_calib[,1], units = "mins")))

hist(twl_dev, freq = F, breaks = 26)
seq <- seq(0, 80, length = 100)
fitml_ng <- fitdistr(twl_dev, "gamma")
lines(seq, dgamma(seq, fitml_ng$estimate[1], fitml_ng$estimate[2]), col = "firebrick", lwd = 3, lty = 2)

twilight.dev[[1]] <- list(zenith = c(median(z), max(z)), twl_dev)

### 
twl_dev_all0 <- twilight(twl[,1], lon.breed, lat.breed, rise = twl[,2], zenith = max(z)+0.1)
twl_dev_all  <- ifelse(twl$Rise, as.numeric(difftime(twl[,1], twl_dev_all0, units = "mins")),
                       as.numeric(difftime(twl_dev_all0, twl[,1], units = "mins")))


zenith <- calib(twl[twl$Twilight<=as.POSIXct("2013-06-26"),], lat.breed, 96)
crds <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith)
plot(crds$x[crds$time<as.POSIXct("2013-06-26", tz = "GMT"),], xlim = c(-95, -70), ylim =c(15, 45))
plot(wrld_simpl, add = T)
points(lon.breed, lat.breed, pch = 22, bg = "white", cex = 2, lwd = 2)

crds.breed[[1]]   <- data.frame(crds$time, crds$x)

out <- data.frame(id1 = "GCM05", id2 = "GTN05", time = crds$time, zenithT = median(z), zenith = max(z), tw_error = twl_dev_all, lon = crds$x[,1], lat =  crds$x[,2])



### GCM06
GCM06 <- readLig("LightData/GCM062014.lig") 
GCM06 <- GCM06[-1,]
  GCM06$Date <- GCM06$Date+(5*60) ## time adjustment

lightImage(GCM06, offset = offset)
tsimagePoints(c.dat$Twilight, offset = offset, pch = 16, cex = 0.25,
              col = "orange")


# twl <- TwGeos::preprocessLight(GCM06[-1,], threshold = threshold, offset = offset, gr.Device = "x11")
# write.csv(twl, "LightData/twl_GCM06.csv", row.names = F)
twl <- read.csv("LightData/twl_GCM06.csv")
  twl$Twilight <- as.POSIXct(twl$Twilight, tz = "GMT")
  twl <- twilightAdjust(twl, 120)

calib.tm  <- c(as.POSIXct("2013-05-09", tz = "GMT"), as.POSIXct("2013-06-09", tz = "GMT"))  
twl_calib <- subset(twl, Twilight>=calib.tm[1] & Twilight<=calib.tm[2])

### get zenith angle and alpha parameters
sun  <- solar(twl_calib[,1])
z    <- refracted(zenith(sun, lon.breed, lat.breed))
plot(z)

twl_t   <- twilight(twl_calib[,1], lon.breed, lat.breed, rise = twl_calib[,2], zenith = max(z)+0.1)
twl_dev <- ifelse(twl_calib$Rise, as.numeric(difftime(twl_calib[,1], twl_t, units = "mins")),
                  as.numeric(difftime(twl_t, twl_calib[,1], units = "mins")))

hist(twl_dev, freq = F, breaks = 26)
seq <- seq(0, 80, length = 100)
fitml_ng <- fitdistr(twl_dev, "gamma")
lines(seq, dgamma(seq, fitml_ng$estimate[1], fitml_ng$estimate[2]), col = "firebrick", lwd = 3, lty = 2)

twilight.dev[[2]] <- list(zenith = c(median(z), max(z)), twl_dev)


twl_dev_all0 <- twilight(twl[,1], lon.breed, lat.breed, rise = twl[,2], zenith = max(z)+0.1)
twl_dev_all  <- ifelse(twl$Rise, as.numeric(difftime(twl[,1], twl_dev_all0, units = "mins")),
                       as.numeric(difftime(twl_dev_all0, twl[,1], units = "mins")))

###
zenith <- calib(twl[twl$Twilight<=as.POSIXct("2013-06-26"),], lat.breed, 96)
crds <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith)
plot(crds$x[crds$time<as.POSIXct("2013-06-26", tz = "GMT"),], xlim = c(-95, -70), ylim =c(15, 45))
plot(wrld_simpl, add = T)
points(lon.breed, lat.breed, pch = 22, bg = "white", cex = 2, lwd = 2)

crds.breed[[2]]   <- data.frame(crds$time, crds$x)

out <- rbind(out,
             data.frame(id1 = "GCM06", id2 = "GTN06", time = crds$time, zenithT = median(z), zenith = max(z), tw_error = twl_dev_all, lon = crds$x[,1], lat =  crds$x[,2]))


### GCM09
GCM09 <- readLig("LightData/GCM092014.lig")
GCM09 <- GCM09[-1,]
  GCM09$Date <- GCM09$Date+(5*60) ## time adjustment

lightImage(GCM09, offset = offset)
tsimagePoints(c.dat$Twilight, offset = offset, pch = 16, cex = 0.25,
              col = "orange")

# twl <- TwGeos::preprocessLight(GCM09, threshold = threshold, offset = offset, gr.Device = "x11")
# write.csv(twl, "LightData/twl_GCM09.csv", row.names = F)
twl <- read.csv("LightData/twl_GCM09.csv")
  twl$Twilight <- as.POSIXct(twl$Twilight, tz = "GMT")
  twl <- twilightAdjust(twl, 120)


calib.tm  <- c(as.POSIXct("2013-05-13", tz = "GMT"), as.POSIXct("2013-06-13", tz = "GMT"))  
twl_calib <- subset(twl, Twilight>=calib.tm[1] & Twilight<=calib.tm[2])

### get zenith angle and alpha parameters
sun  <- solar(twl_calib[,1])
z    <- refracted(zenith(sun, lon.breed, lat.breed))
plot(z)

twl_t   <- twilight(twl_calib[,1], lon.breed, lat.breed, rise = twl_calib[,2], zenith = max(z)+0.1)
twl_dev <- ifelse(twl_calib$Rise, as.numeric(difftime(twl_calib[,1], twl_t, units = "mins")),
                  as.numeric(difftime(twl_t, twl_calib[,1], units = "mins")))

hist(twl_dev, freq = F, breaks = 26)
seq <- seq(0, 80, length = 100)
fitml_ng <- fitdistr(twl_dev, "gamma")
lines(seq, dgamma(seq, fitml_ng$estimate[1], fitml_ng$estimate[2]), col = "firebrick", lwd = 3, lty = 2)

twilight.dev[[3]] <- list(zenith = c(median(z), max(z)), twl_dev)

twl_dev_all0 <- twilight(twl[,1], lon.breed, lat.breed, rise = twl[,2], zenith = max(z)+0.1)
twl_dev_all  <- ifelse(twl$Rise, as.numeric(difftime(twl[,1], twl_dev_all0, units = "mins")),
                       as.numeric(difftime(twl_dev_all0, twl[,1], units = "mins")))


###
zenith <- calib(twl[twl$Twilight<=as.POSIXct("2013-06-26"),], lat.breed, 96)
crds <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith)
plot(crds$x[crds$time<as.POSIXct("2013-06-26", tz = "GMT"),], xlim = c(-95, -70), ylim =c(15, 45))
plot(wrld_simpl, add = T)
points(lon.breed, lat.breed, pch = 22, bg = "white", cex = 2, lwd = 2)

crds.breed[[3]]   <- data.frame(crds$time, crds$x)


out <- rbind(out,
             data.frame(id1 = "GCM09", id2 = "GTN09", time = crds$time, zenithT = median(z), zenith = max(z), tw_error = twl_dev_all, lon = crds$x[,1], lat =  crds$x[,2]))



### GCM13
GCM13 <- readLig("LightData/GCM132014.lig")
GCM13 <- GCM13[-1,]
  GCM13$Date <- GCM13$Date+(5*60) ## time adjustment

lightImage(GCM13, offset = offset)
tsimagePoints(c.dat$Twilight, offset = offset, pch = 16, cex = 0.25,
              col = "orange")

# twl <- TwGeos::preprocessLight(GCM13, threshold = threshold, offset = offset, gr.Device = "x11")
# write.csv(twl, "LightData/twl_GCM13.csv", row.names = F)
twl <- read.csv("LightData/twl_GCM13.csv")
  twl$Twilight <- as.POSIXct(twl$Twilight, tz = "GMT")
  twl <- twilightAdjust(twl, 120)

calib.tm  <- c(as.POSIXct("2013-05-16", tz = "GMT"), as.POSIXct("2013-06-16", tz = "GMT"))  
twl_calib <- subset(twl, Twilight>=calib.tm[1] & Twilight<=calib.tm[2])

### get zenith angle and alpha parameters
sun  <- solar(twl_calib[,1])
z    <- refracted(zenith(sun, lon.breed, lat.breed))
plot(z)

twl_t   <- twilight(twl_calib[,1], lon.breed, lat.breed, rise = twl_calib[,2], zenith = max(z)+0.1)
twl_dev <- ifelse(twl_calib$Rise, as.numeric(difftime(twl_calib[,1], twl_t, units = "mins")),
                  as.numeric(difftime(twl_t, twl_calib[,1], units = "mins")))

hist(twl_dev, freq = F, breaks = 26)
seq <- seq(0, 80, length = 100)
fitml_ng <- fitdistr(twl_dev, "gamma")
lines(seq, dgamma(seq, fitml_ng$estimate[1], fitml_ng$estimate[2]), col = "firebrick", lwd = 3, lty = 2)

twilight.dev[[4]] <- list(zenith = c(median(z), max(z)), twl_dev)


twl_dev_all0 <- twilight(twl[,1], lon.breed, lat.breed, rise = twl[,2], zenith = max(z)+0.1)
twl_dev_all  <- ifelse(twl$Rise, as.numeric(difftime(twl[,1], twl_dev_all0, units = "mins")),
                       as.numeric(difftime(twl_dev_all0, twl[,1], units = "mins")))


###
zenith <- calib(twl[twl$Twilight<=as.POSIXct("2013-06-26"),], lat.breed, 96)
crds <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith)
plot(crds$x[crds$time<as.POSIXct("2013-06-26", tz = "GMT"),], xlim = c(-95, -70), ylim =c(15, 45))
plot(wrld_simpl, add = T)
points(lon.breed, lat.breed, pch = 22, bg = "white", cex = 2, lwd = 2)

crds.breed[[4]]   <- data.frame(crds$time, crds$x)

out <- rbind(out,
             data.frame(id1 = "GCM13", id2 = "GTN13", time = crds$time, zenithT = median(z), zenith = max(z), tw_error = twl_dev_all, lon = crds$x[,1], lat =  crds$x[,2]))



### GCM16
GCM16 <- readLig("LightData/GCM162014.lig")
GCM16 <- GCM16[-1,]
  GCM16$Date <- GCM16$Date+(7.5*60) ## time adjustment

lightImage(GCM16, offset = offset)
tsimagePoints(c.dat$Twilight, offset = offset, pch = 16, cex = 0.25,
              col = "orange")

# twl <- TwGeos::preprocessLight(GCM16, threshold = threshold, offset = offset, gr.Device = "x11")
# write.csv(twl, "LightData/twl_GCM16.csv", row.names = F)
twl <- read.csv("LightData/twl_GCM16.csv")
  twl$Twilight <- as.POSIXct(twl$Twilight, tz = "GMT")
  twl <- twilightAdjust(twl, 120)

calib.tm  <- c(as.POSIXct("2013-05-16", tz = "GMT"), as.POSIXct("2013-06-16", tz = "GMT"))  
twl_calib <- subset(twl, Twilight>=calib.tm[1] & Twilight<=calib.tm[2])

### get zenith angle and alpha parameters
sun  <- solar(twl_calib[,1])
z    <- refracted(zenith(sun, lon.breed, lat.breed))
plot(z)

twl_t   <- twilight(twl_calib[,1], lon.breed, lat.breed, rise = twl_calib[,2], zenith = max(z)+0.1)
twl_dev <- ifelse(twl_calib$Rise, as.numeric(difftime(twl_calib[,1], twl_t, units = "mins")),
                  as.numeric(difftime(twl_t, twl_calib[,1], units = "mins")))

hist(twl_dev, freq = F, breaks = 26)
seq <- seq(0, 80, length = 100)
fitml_ng <- fitdistr(twl_dev, "gamma")
lines(seq, dgamma(seq, fitml_ng$estimate[1], fitml_ng$estimate[2]), col = "firebrick", lwd = 3, lty = 2)

twilight.dev[[5]] <- list(zenith = c(median(z), max(z)), twl_dev)


twl_dev_all0 <- twilight(twl[,1], lon.breed, lat.breed, rise = twl[,2], zenith = max(z)+0.1)
twl_dev_all  <- ifelse(twl$Rise, as.numeric(difftime(twl[,1], twl_dev_all0, units = "mins")),
                       as.numeric(difftime(twl_dev_all0, twl[,1], units = "mins")))

###
zenith <- calib(twl[twl$Twilight<=as.POSIXct("2013-06-26"),], lat.breed, 96)
crds <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith)
plot(crds$x[crds$time<as.POSIXct("2013-06-26", tz = "GMT"),], xlim = c(-95, -70), ylim =c(15, 45))
plot(wrld_simpl, add = T)
points(lon.breed, lat.breed, pch = 22, bg = "white", cex = 2, lwd = 2)

crds.breed[[5]]   <- data.frame(crds$time, crds$x)


out <- rbind(out,
             data.frame(id1 = "GCM16", id2 = "GTN16", time = crds$time, zenithT = median(z), zenith = max(z), tw_error = twl_dev_all, lon = crds$x[,1], lat =  crds$x[,2]))


# save(crds.breed,   file = "ReAnalysis_crds.breed.RData")
# save(twilight.dev, file = "ReAnalysis_twl_dev.RData")
# save(out,   file = "allData.RData")
# load("ReAnalysis_twl_dev.RData")
# load("ReAnalysis_crds.breed.RData")



#############
### Plot ####
#############

dev <- unlist(lapply(twilight.dev, function(x) x[[2]]))
seq <- seq(0, 60, length = 100)

opar <- par(mar = c(6,6,1,1), mgp = c(3.8,0.9,0), cex.lab =1.4)
bp <- hist(dev, xlim = c(0, 60), las = 1, 
           ylab = "Density distribution", xlab = "Twiligth error (min)",
           main = "", col = "grey80",
           breaks = seq(0, 60, by = 1.5), freq = FALSE, ylim = c(0, 0.08))

for(i in 1:length(twilight.dev)) {
  fitml_ng <- fitdistr(twilight.dev[[i]][[2]], "gamma")
  lines(seq, dgamma(seq, fitml_ng$estimate[1], fitml_ng$estimate[2]), 
        col = adjustcolor(cols[[i]], alpha.f = 0.8) , lwd = 2, lty = 1)
}


fitml_ng <- fitdistr(dev, "gamma")
# lines(seq, dgamma(seq, fitml_ng$estimate[1], fitml_ng$estimate[2]), col = "firebrick", lwd = 3, lty = 2)
legend("topright", names, lty = 1, lwd = 2, pch = NA, col = cols, bty = "n", cex = 1.4)
par(opar)