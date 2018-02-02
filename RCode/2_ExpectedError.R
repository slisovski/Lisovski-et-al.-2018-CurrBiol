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



##########################
## Expected Error Range ##
##########################

## Map boundaries ###
xlim = c(-95, -70)
ylim = c(15, 45)


## Paths extracted from Streby et al. 2016 ####
dates <- seq(as.POSIXct('2014-04-26 12:00:00', tz = "GMT"), as.POSIXct('2014-05-02 12:00:00', tz = "GMT"), by = "day")

paths <- list(Date = dates,
              GTN05 = cbind(c(-85.29389, -85.29389, -82.48356, -86.83244, -86.83244, -85.37684, -85.29389),
                            c(36.26889, 36.26889, 30.88530, 30.41831, 30.41831, 33.62979, 36.26889)),
              GTN06 = cbind(c(-85.29389, -85.29389, -85.37087, -82.56917, -85.03393, -86.15258, -85.29389),
                            c(36.26889, 36.26889, 29.78579, 27.8102, 22.58259, 34.96984, 36.26889)),
              GTN09 = cbind(c(-93.25813, -85.29389, -90.05831, -86.17697, -86.17697, -85.25954, -85.29389),
                            c(31.50613, 36.26889, 30.90732, 30.36583, 30.36583, 33.00142, 36.26889)),
              GTN13= cbind(c(-85.29389, -82.32519, -83.07325, -86.83768, -86.83768, -83.99534, -85.29389),
                           c(36.26889, 33.48906, 32.22929, 30.46721, 30.46721, 35.78956, 36.26889)),
              GTN16= cbind(c(-85.29389, -85.29389, -81.313, -85.64561, -85.64561, -83.63923, -85.29389),
                           c(36.26889, 36.26889, 34.07701, 30.22698, 30.22698, 36.02401, 36.26889)))


## Simulation of locations based on twilight error (takes time - go for coffee!)
## Could be optimised using parallel computing (e.g. parApply)

loc = 50000 # number of locations
crds.out <- cbind(NA, NA)
for(i in 1:loc) {
  sr <- sunrise(seq(as.POSIXct("2014-04-26", tz = "EST"), as.POSIXct("2014-05-02", tz = "EST"), by = "days"),
                lon = lon.breed, lat = lat.breed, iters = 5)+(rgamma(7, fitml_ng$estimate[1], fitml_ng$estimate[2])*60)
  ss <- sunset(seq(as.POSIXct("2014-04-26", tz = "EST"), as.POSIXct("2014-05-02", tz = "EST"), by = "days"),
               lon = lon.breed, lat = lat.breed, iters = 5)-(rgamma(7, fitml_ng$estimate[1], fitml_ng$estimate[2])*60)
  
  twl <- data.frame(Twilight = as.POSIXct(c(sr, ss), origin = "1970-01-01", tz = "GTM"), 
                    Rise = c(rep(TRUE, length(sr)), rep(FALSE, length(ss))))
  
  crds <- thresholdLocation(twl[order(twl[,1]),][,1], twl[order(twl[,1]),][,2], zenith = 94.81235+runif(1, -1.2, 1.2))
  
  crds.out <- rbind(crds.out, crds$x)
}


r <- raster(xmn = xlim[1], xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], res = c(0.5, 0.5))
rpt <- rasterize(crds.out, r, fun = "sum")


### normalization
plot((rpt[]/cellStats(rpt, "max")) * 100)


scl <- function(x) {
  (x - min(x,na.rm = TRUE)) / diff(range(x, na.rm = TRUE))
}
subx <- setValues(rpt, scl(values(rpt)))


plot(NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "", col = "grey97", border = NA)
plot(subx, add =T, legend = T, col = rev(terrain.colors(200)))
plot(wrld_simpl, col = "transparent", border = "grey20", add = T)
points(lon.breed, lat.breed, pch = 23, bg = "white", lwd = 2, cex = 2)
contour(subx, levels = c(0.05, 10, 40, 60)/100, add = TRUE, drawlabels = FALSE)

for(i in 1:length(paths)) {
  crds <- data.frame(dates, paths[[i+1]])
  lines(crds[2:6,2], crds[2:6,3], type = "p", pch = "x", cex = 1.6,
        col = "grey10")
}

axis(1)
axis(2, las = 1)
mtext("Longitude", 1, line = 3.5, cex = 1.4)
mtext("Latitude", 2, line = 3.5, cex = 1.4)



### Location estimates during the previous breeding season (calibration period - not affected by tornado)

plot(NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "", col = "grey97", border = "grey20")
plot(wrld_simpl, add = T)
for(i in c(1,2,3,4,5)) {
  tmp <- crds.breed[[i]][crds.breed[[i]][,1]<=as.POSIXct("2013-06-30"),-1]
  tmp[,2] <- tmp[,2]
  points(tmp, col = adjustcolor(cols[i], alpha.f = 0.9), pch = 16, cex = 1.1, type = "o")
}
points(lon.breed, lat.breed, pch = 23, bg = "white", lwd = 2, cex = 2)