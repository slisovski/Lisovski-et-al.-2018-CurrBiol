#############################################
#### Lisovski et al. 2018 Current Biology ###
#############################################

library(TwGeos)
library(SGAT)
library(geosphere)
library(GeoLight)
library(MASS)
library(raster)
library(rgdal)
library(rgeos)
library(maptools); data(wrld_simpl)

lon.breed <- -84.29
lat.breed <-  36.27


## Load results from Twiligth error (1_TwilightError.r)
# load("ReAnalysis_twl_dev.RData")
# load("ReAnalysis_crds.breed.RData")


###########################
## Latitudinal deviation ##
###########################

# precipitable water can be downloaded from NCEP via RNCEP package
load("C:/Users/SLi/Documents/GitHub/Lisovski-et-al.-2018-CurrBiol/Data/Res_2013.RData")
load("C:/Users/SLi/Documents/GitHub/Lisovski-et-al.-2018-CurrBiol/Data/Res_2014.RData")

prec13 <- data.frame(date = as.POSIXct(Res_2013[,1], origin = "1970-01-01"), prec = Res_2013[,2])
prec14 <- data.frame(date = as.POSIXct(Res_2014[,1], origin = "1970-01-01"), prec = Res_2014[,2])

range1 <- as.POSIXct(c("2013-05-29", "2013-06-26"), tz = "GMT")
range2 <- as.POSIXct(c("2014-04-07", "2014-05-05"), tz = "GMT")

prec <- rbind(prec13, prec14)
out$prec = NA

for(i in 1:nrow(out)) {
  
  diff <- abs(as.numeric(difftime(out$time[i], prec$date, units = "hours")))
  ind  <- which.min(diff)
  if(diff[ind]<20) out$prec[i] <- prec[ind,2]
  
}

for(i in unique(out$id1)) {
  tmp <- subset(out, id1==i & time >= range1[1] & time<=range1[2])
  
  if(i==unique(out$id1)[1]) {
    twl_d <- matrix(tmp$tw_error, ncol = 1)
    crd_d <- matrix(tmp$lat, ncol = 1)
  } else {
    twl_d <- cbind(twl_d, matrix(tmp$tw_error, ncol = 1))
    crd_d <- cbind(crd_d, matrix(tmp$lat, ncol = 1))
  }
}



opar <- par(mfcol = c(1,2), mar = c(3,3,0.5,3), oma = c(2,2,0,0))

tm = prec13[prec13[,1]<=range1[2] & prec13[,1]>=range1[1] ,1]
p  = prec13[prec13[,1]<=range1[2] & prec13[,1]>=range1[1] ,2]
plot(NA, xlim = range1, ylim = c(5, 100), 
     type = "o", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
polygon(as.numeric(c(tm, rev(tm))), c(p, rep(5, length(tm))), border = NA, col = "grey80")
axis(4, at = seq(0, 50, by = 10), las = 1)

par(new = T)
plot(NA, xlim = range1, ylim = c(21, 40),
     xaxt = "n", ylab = "", xlab = "", yaxt = "n")
axis(2, at = seq(30, 40, by = 5))
for(i in 1:ncol(twl_d)) {
  points(tm, crd_d[,i], col = adjustcolor(cols[i], alpha.f = 0.9), pch = 16, cex = 1.1, type = "o")
}
abline(h = lat.breed)


seqTW <- seq(range1[1], range1[2], by = 12*60*60)
tw  <- data.frame(Twilight = twilight(seqTW, lon = lon.breed, lat = lat.breed, rise = c(TRUE, FALSE), zenith = 96, iters = 5),
                  Rise     = c(rep(c(TRUE, FALSE), 28), TRUE))

tw$Twilight <- tw$Twilight + ifelse(tw$Rise, 20, -20)*60  


crdst <- thresholdPath(tw$Twilight, tw$Rise, zenith = 96)
lines(crdst$time, crdst$x[,2])


axis(1, at = seq(as.POSIXct("2013-05-09", tz = "GMT"), as.POSIXct("2013-06-26", tz = "GMT"), length = 5),
     labels = format(seq(as.POSIXct("2013-05-09", tz = "GMT"), as.POSIXct("2013-06-26", tz = "GMT"), length = 5), "%b-%d"))
text(as.POSIXct("2013-05-13", tz = "GMT"), 40, "a) calibration period (2013)")


tm = prec14[prec14$date<=as.POSIXct("2014-05-04"),1]
p  = prec14[prec14$date<=as.POSIXct("2014-05-04"),2]
plot(NA, xlim = range2+c(12*24*60*60, 0), ylim = c(5, 100), 
     type = "o", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
polygon(as.numeric(c(tm, rev(tm))), c(p, rep(5, length(tm))), border = NA, col = "grey80")
axis(4, at = seq(0, 50, by = 10), las = 1)

par(new = T)
plot(NA, xlim = range2+c(12*24*60*60, 0), ylim = c(21, 40),
     xaxt = "n", ylab = "", xlab = "Time")

i <- 1
tmp <- crds.breed[[i]][crds.breed[[i]][,1]>=as.POSIXct("2014-04-26") & crds.breed[[i]][,1]<=as.POSIXct("2014-05-04"),]
points(tmp[,1], tmp[,3]-1.5, col = adjustcolor(cols[i], alpha.f = 0.9), pch = 16, cex = 1.1, type = "o")

i <- 2
tmp <- crds.breed[[i]][crds.breed[[i]][,1]>=as.POSIXct("2014-04-26") & crds.breed[[i]][,1]<=as.POSIXct("2014-05-04"),]
points(tmp[,1], tmp[,3]-0.6, col = adjustcolor(cols[i], alpha.f = 0.9), pch = 16, cex = 1.1, type = "o")

i <- 3
tmp <- crds.breed[[i]][crds.breed[[i]][,1]>=as.POSIXct("2014-04-26") & crds.breed[[i]][,1]<=as.POSIXct("2014-05-04"),]
points(tmp[,1], tmp[,3]-1, col = adjustcolor(cols[i], alpha.f = 0.9), pch = 16, cex = 1.1, type = "o")

i <- 4
tmp <- crds.breed[[i]][crds.breed[[i]][,1]>=as.POSIXct("2014-04-26") & crds.breed[[i]][,1]<=as.POSIXct("2014-05-04"),]
points(tmp[,1], tmp[,3]-5, col = adjustcolor(cols[i], alpha.f = 0.9), pch = 16, cex = 1.1, type = "o")

i <- 5
tmp <- crds.breed[[i]][crds.breed[[i]][,1]>=as.POSIXct("2014-04-26") & crds.breed[[i]][,1]<=as.POSIXct("2014-05-04"),]
points(tmp[,1], tmp[,3]+3, col = adjustcolor(cols[i], alpha.f = 0.9), pch = 16, cex = 1.1, type = "o")


abline(h = lat.breed)



seqTW <- seq(range2[1], range2[2], by = 12*60*60)
tw  <- data.frame(Twilight = twilight(seqTW, lon = lon.breed, lat = lat.breed, rise = c(TRUE, FALSE), zenith = 96, iters = 5),
                  Rise     = c(rep(c(TRUE, FALSE), length(seqTW)/2), TRUE))

tw$Twilight <- tw$Twilight + ifelse(tw$Rise, 20, -20)*60  


crdst <- thresholdPath(tw$Twilight, tw$Rise, zenith = 96)

lines(crdst$time, crdst$x[,2])


axis(1, at = seq(as.POSIXct("2014-04-25", tz = "GMT"), as.POSIXct("2014-05-04", tz = "GMT"), length = 3),
     labels = format(seq(as.POSIXct("2014-04-25", tz = "GMT"), as.POSIXct("2014-05-04", tz = "GMT"), length = 3), "%b-%d"))


par(opar)
#####