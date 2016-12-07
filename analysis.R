##########
# This file contains all data process and analysis code for smoothing GPS data
# by Kalman filter.


# import functions
source('~/GitProjects/denoising_GPS_data/kalman_filter_functions.R')


# ------ process data and store data chunks into separated csv files -------

# read data
setwd('~/GitProjects/denoising_GPS_data/')
gps.data <- read.csv('data/campus_gammaradiation.csv', stringsAsFactors=FALSE)
# long keep three columns -- 'timestamp', 'lon', 'lat'
gps.data <- gps.data[ , c('timestamp', 'lon', 'lat')]
# convert `timestamp` from characters to time format
gps.data$timestamp <- as.POSIXct(gps.data$timestamp, format='%Y-%m-%d %H:%M:%S')

#sort data using timestamp
sorted_gpsdata <- gps.data[order(gps.data$timestamp), ]

# check time interval between two sample times.
time.shift <- shift_vec(sorted_gpsdata$timestamp, 1)
time.diff <- as.numeric(difftime(sorted_gpsdata$timestamp, time.shift))
time.diff[is.na(time.diff)] <- 1

# split data based on time interval
# time interval = 30s
index <- c(which (time.diff > 0.5*60), 814458)

length(which((index - shift_vec(index, 1)) >= 100))

# data chunk into csv; only data chunk in which time interval between samples
# is less than 30s, and sample size is larger than 100.
num <- 1
for (i in 1:length(index)) {
    if (i == 1) {
        diff <- index[i] - 1
    } else {
        diff <- index[i] - index[i-1]
    }
    if (diff >= 100) {
        data <- sorted_gpsdata[index[i-1]:(index[i]-1), ]
        # convert longitude and latitude to UTM
        for (j in 1:nrow(data)) {
            UTM.location <- longlatToUTM(data[j, 2], data[j, 3])
            data$UTMlon[j] <- UTM.location[1,1]
            data$UTMlat[j] <- UTM.location[1,2]
        }
        # calculate time difference
        time.shift <- shift_vec(data$timestamp, 1)
        time.diff <- as.numeric(difftime(data$timestamp, time.shift))
        time.diff[is.na(time.diff)] <- 0
        data$timediff <- time.diff
        # calculate speed on longitude and latitude
        # on longitude
        lon.shift <- shift_vec(data$UTMlon, -1)
        lon.diff <- as.numeric(lon.shift - data$UTMlon)
        vLon <- lon.diff[1:(length(lon.diff)-1)] / time.diff[-1]
        data$vLon <- c(vLon, vLon[length(vLon)])
        # on latitude
        lat.shift <- shift_vec(data$UTMlat, -1)
        lat.diff <- as.numeric(lat.shift - data$UTMlat)
        vLat <- lat.diff[1:(length(lat.diff)-1)] / time.diff[-1]
        data$vLat <- c(vLat, vLat[length(vLat)])
        write.csv(data, paste('data/', num, '.csv', sep=''))
        num <- num + 1
    }
}


# ---------------------- Smoothing GPS data -----------------------------

# read data samples
for (i in 379:1150) {
    csv.name <- paste('data/', i, '.csv', sep='')
    data <- read.csv(csv.name)
    pngname <- paste(i, '.png', sep='')
    plotGPS_png(data, long=3, lat=4, pngname)
}


# read data 289
i <- 289
csv.name <- paste('data/', i, '.csv', sep='')
data <- read.csv(csv.name)


lon.shift <- shift_vec(data$lon, -1)
lon.diff <- as.numeric(lon.shift - data$lon)
vLon <- lon.diff[1:(length(lon.diff)-1)] / data$timediff[-1]
data$vLon <- c(0, vLon)
# on latitude
lat.shift <- shift_vec(data$lat, -1)
lat.diff <- as.numeric(lat.shift - data$lat)
vLat <- lat.diff[1:(length(lat.diff)-1)] / data$timediff[-1]
data$vLat <- c(0, vLat)


###################################################
# determine matrix Q and R using empirical method #
###################################################

# transformation matrix
H.matrix <- matrix(0, nrow=2, ncol=4)
H.matrix[1, 1] <- 1
H.matrix[2, 2] <- 1

# determine matrix w and R
w.matrix <- var(data[, c('vLon', 'vLat')]) #* 1000000
R.matrix <- var(data[, c('UTMlon', 'UTMlat')])

# initial Xhat.posterior
Xhat.initial <- matrix(c(data$lon[1], data$lat[1], data$vLon[1], data$vLat[1]), nrow=4)

# covariance matrix of posterior distribution. -- it does not matter what the
# initial is. It will only affect the speed of convergence.
P.initial <- diag(c(var(data$lon), var(data$lat), var(data$vLon), car(data$vLat)), 4)


# matrices storing denoising results
location.results <- matrix(0, nrow=4, ncol=nrow(data))
location.results[, 1] <- Xhat.initial
P.results <- matrix(0, nrow=16, ncol=nrow(data))
P.results[, 1] <- matrix(P.initial, ncol=1)


# run Kalman filter
for (i in 2:nrow(data)) {
    Xhat.initial <- location.results[, i-1]
    P.initial <- matrix(P.results[, i-1], ncol=4)
    locations <- t(data[i, c('timediff', 'lon', 'lat')])
    # run kalman filter
    results <- kalman_update(locations, Xhat.initial, P.initial, H.matrix, w.matrix, R.matrix)
    location.results[, i] <- results$Xhat.posterior
    P.results[, i] <- matrix(results$P.posterior, ncol=1)
}

location.results.t <- t(location.results)

# plot longitude and latitude vs time separately
par(mfrow=c(2,1), mar=c(4,5,2,2))
plot(data[,3], type='l', col='blue', lwd=3, bty='n', ylab='Longitude', xlab='Time')
lines(location.results.t[,1], col='red', lwd=3)
legend('topleft', c('true data', 'denoised data'), lty=c(1,1), lwd=c(3,3), col=c('blue', 'red'))
title(main='Empirical Method')
plot.ts(data[,4], type='l', col='blue', lwd=3, bty='n', ylab='Latitude', xlab='Time')
lines(location.results.t[,2], col='red', lwd=3)

# plot longitude vs latitude
plot(data[,3], data[,4], bty='n', xlab='Longitude', ylab='Latitude')
lines(location.results.t[, 1], location.results.t[, 2], col='red', lwd=3)

# plot longitude and latitude on map
plotGPS_png(location.results.t, long=1, lat=2, '289ss.png', width=1000, height=800)

# plot convergence of variance of both longitude and latitude
par(mfrow=c(2,1), mar=c(4,5,2,2))
plot(P.results[1,], type='l', col='blue', lwd=3, bty='n', ylab='Variance', xlab='Time')
lines(P.results[6, ], col='red', lwd=3)
legend('topright', c('Variance of Longitude', 'Variance of Latitude'), lty=c(1,1), lwd=c(3,3), col=c('blue', 'red'))
title(main='Empirical Method')


#####################################################################
# determine matrix Q and R using Autocovariance Least Square method #
#####################################################################

# state transition matrix F
F.matrix <- diag(1, 4)
F.matrix[1, 3] <- 1
F.matrix[2, 4] <- 1

# positions of delta.time in matrix F
positions <- list(c(1,3), c(2,4))

# transformation matrix
H.matrix <- matrix(0, nrow=2, ncol=4)
H.matrix[1, 1] <- 1
H.matrix[2, 2] <- 1

# initialize matrix K for function QR
K.guess <- matrix(rep(0.5), nrow=4, ncol=2)

# determine matrix Q and R
QR <- ALS(data[, c(7,3,4)], H.matrix, F.matrix, K.guess, N=55, positions=positions)

# process noise covariance matrix
Q.matrix <- QR$Q.matrix #/ 100000000

# measurement noise covariance matrix
R.matrix <- QR$R.matrix

# initial Xhat.posterior
Xhat.initial <- matrix(c(data$lon[1], data$lat[1], data$vLon[1], data$vLat[1]), nrow=4)

# covariance matrix of posterior distribution. -- it does not matter what the
# initial is. It will only affect the speed of convergence.
P.initial <- diag(c(var(data$lon), var(data$lat), var(data$vLon), car(data$vLat)), 4)


# matrices storing denoising results
location.results <- matrix(0, nrow=4, ncol=nrow(data))
location.results[, 1] <- Xhat.initial
P.results <- matrix(0, nrow=16, ncol=nrow(data))
P.results[, 1] <- matrix(P.initial, ncol=1)

# run Kalman filter
for (i in 2:nrow(data)) {
    Xhat.initial <- location.results[, i-1]
    P.initial <- matrix(P.results[, i-1], ncol=4)
    locations <- t(data[i, c('timediff', 'lon', 'lat')])
    # run kalman filter
    results <- kalman_update2(locations, Xhat.initial, P.initial, F.matrix, H.matrix, Q.matrix, R.matrix, positions)
    location.results[, i] <- results$Xhat.posterior
    P.results[, i] <- matrix(results$P.posterior, ncol=1)
}

location.results.t <- t(location.results)

# plot longitude and latitude vs time separately
par(mfrow=c(2,1), mar=c(4,5,2,2))
plot(data[,3], type='l', col='blue', lwd=3, bty='n', ylab='Longitude', xlab='Time')
lines(location.results.t[,1], col='red', lwd=3)
legend('topleft', c('true data', 'denoised data'), lty=c(1,1), lwd=c(3,3), col=c('blue', 'red'))
title(main='Aucovariance Least-Squared Method')
plot.ts(data[,4], type='l', col='blue', lwd=3, bty='n', ylab='Latitude', xlab='Time')
lines(location.results.t[,2], col='red', lwd=3)


# plot longitude vs latitude
plot(data[,3], data[,4], bty='n', xlab='Longitude', ylab='Latitude')
lines(location.results.t[, 1], location.results.t[, 2], col='red', lwd=3)

# plot longitude and latitude on map
plotGPS_png(location.results.t, long=1, lat=2, '289s.png', width=1000, height=800)

# plot convergence of variance of both longitude and latitude
plot(P.results[1,], type='l', col='blue', lwd=3, bty='n', ylab='Variance', xlab='Time', ylim=c(500,4000))
lines(P.results[6, ], col='red', lwd=3)
legend('topright', c('Variance of Longitude', 'Variance of Latitude'), lty=c(1,1), lwd=c(3,3), col=c('blue', 'red'))
title(main='Empirical Method')
