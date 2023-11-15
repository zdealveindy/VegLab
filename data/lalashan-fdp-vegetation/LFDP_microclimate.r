# Function for plotting wind rose ----
plot_windrose <- function(data,
                          spd,
                          dir,
                          spdres = 2,
                          dirres = 30,
                          spdmin = 2,
                          spdmax = 20,
                          spdseq = NULL,
                          palette = "YlGnBu",
                          countmax = NA,
                          debug = 0){
  require(ggplot2)
  require(RColorBrewer)
  
  # Look to see what data was passed in to the function
  if (is.numeric(spd) & is.numeric(dir)){
    # assume that we've been given vectors of the speed and direction vectors
    data <- data.frame(spd = spd,
                       dir = dir)
    spd = "spd"
    dir = "dir"
  } else if (exists("data")){
    # Assume that we've been given a data frame, and the name of the speed 
    # and direction columns. This is the format we want for later use.    
  }  
  
  # Tidy up input data ----
  n.in <- NROW(data)
  dnu <- (is.na(data[[spd]]) | is.na(data[[dir]]))
  data[[spd]][dnu] <- NA
  data[[dir]][dnu] <- NA
  
  # figure out the wind speed bins ----
  if (missing(spdseq)){
    spdseq <- seq(spdmin,spdmax,spdres)
  } else {
    if (debug >0){
      cat("Using custom speed bins \n")
    }
  }
  # get some information about the number of bins, etc.
  n.spd.seq <- length(spdseq)
  n.colors.in.range <- n.spd.seq - 1
  
  # create the color map
  spd.colors <- colorRampPalette(brewer.pal(min(max(3,
                                                    n.colors.in.range),
                                                min(9,
                                                    n.colors.in.range)),                                               
                                            palette))(n.colors.in.range)
  
  if (max(data[[spd]],na.rm = TRUE) > spdmax){    
    spd.breaks <- c(spdseq,
                    max(data[[spd]],na.rm = TRUE))
    spd.labels <- c(paste(c(spdseq[1:n.spd.seq-1]),
                          '-',
                          c(spdseq[2:n.spd.seq])),
                    paste(spdmax,
                          "-",
                          max(data[[spd]],na.rm = TRUE)))
    spd.colors <- c(spd.colors, "grey50")
  } else{
    spd.breaks <- spdseq
    spd.labels <- paste(c(spdseq[1:n.spd.seq-1]),
                        '-',
                        c(spdseq[2:n.spd.seq]))    
  }
  data$spd.binned <- cut(x = data[[spd]],
                         breaks = spd.breaks,
                         labels = spd.labels,
                         ordered_result = TRUE)
  # clean up the data
  data. <- na.omit(data)
  
  # figure out the wind direction bins
  dir.breaks <- c(-dirres/2,
                  seq(dirres/2, 360-dirres/2, by = dirres),
                  360+dirres/2)  
  # dir.labels <- c(paste(360-dirres/2,"-",dirres/2),
  #                 paste(seq(dirres/2, 360-3*dirres/2, by = dirres),
  #                       "-",
  #                       seq(3*dirres/2, 360-dirres/2, by = dirres)),
  #                 paste(360-dirres/2,"-",dirres/2))
  dir.labels <- (dir.breaks[-length(dir.breaks)] + dir.breaks[-1])/2
  dir.labels[length(dir.labels)] <- dir.labels[1]
  # assign each wind direction to a bin
  dir.binned <- cut(data[[dir]],
                    breaks = dir.breaks,
                    ordered_result = TRUE)
  levels(dir.binned) <- dir.labels
  data$dir.binned <- dir.binned
  
  # Run debug if required ----
  if (debug>0){    
    cat(dir.breaks,"\n")
    cat(dir.labels,"\n")
    cat(levels(dir.binned),"\n")       
  }  
  
  # deal with change in ordering introduced somewhere around version 2.2
  if(packageVersion("ggplot2") > "2.2"){    
    cat("Hadley broke my code\n")
    data$spd.binned = with(data, factor(spd.binned, levels = rev(levels(spd.binned))))
    spd.colors = rev(spd.colors)
  }
  
  # create the plot ----
  p.windrose <- ggplot(data = data,
                       aes(x = dir.binned,
                           fill = spd.binned)) +
    geom_bar() + 
    scale_x_discrete(drop = FALSE,
                     labels = waiver()) +
    coord_polar(start = -((dirres/2)/360) * 2*pi) +
    scale_fill_manual(name = "Wind speed (m/s)", 
                      values = spd.colors,
                      drop = FALSE) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          #panel.border = element_rect(colour = "blank"),
          panel.grid.major = element_line(colour="grey65"))
  
  # adjust axes if required
  if (!is.na(countmax)){
    p.windrose <- p.windrose +
      ylim(c(0,countmax))
  }
  
  # print the plot
  print(p.windrose)  
  
  # return the handle to the wind rose
  return(p.windrose)
}

# Temperature and precipitation from weather station ----
# setwd("p:\\VegLab\\Weather_Station-HOBO-TMS\\Weather station data\\R\\")
library(lubridate)
library(ggplot2)

## temperature -----
TEM <- read.csv("p:\\VegLab\\Weather_Station-HOBO-TMS\\Weather station data\\R\\data&code_for paper\\temperature.csv")
TEM$time <- mdy_hm(TEM$time)

# Mean maximum daily temperature 2
# Mean minimum daily temperature 3
# Absolute monthly minimum temperature 4
# Mean daily temperature 5
TEM$floor_day <- floor_date(TEM$time, "day")

TEM_cal <- aggregate(temp_roll ~ floor_day, data = TEM, function(x) c(min(x), max(x), mean(x)))
t <- cbind(TEM_cal$temp_roll[,1], TEM_cal$temp_roll[,2], TEM_cal$temp_roll[,3])
TEM_cal[,2:4] <- t
colnames(TEM_cal) <- c("floor_day", "Min_t", "Max_t", "Mean_t")

TEM_cal$floor_month <- floor_date(TEM_cal$floor_day, "month")
t <- cbind(aggregate(cbind(Min_t, Max_t, Mean_t) ~ floor_month, data = TEM_cal, FUN = mean), aggregate(Min_t~ floor_month, data = TEM_cal, FUN = min)[,2])
TEM_cal <- t[, c(1,2,3,5,4)]
colnames(TEM_cal)[4] <- "Ab_min_t"


## precipitation -----
PRE <- read.csv("p:\\VegLab\\Weather_Station-HOBO-TMS\\Weather station data\\R\\data&code_for paper\\weather_station.csv")
PRE$time <- mdy_hm(PRE$time)
PRE$floor_month <- floor_date(PRE$time, "month")
PRE_cal <- aggregate(rain ~ floor_month, data = PRE, FUN = sum)

# plot
TEM_PRE <- cbind(TEM_cal[, c(1,5)], PRE_cal$rain)
colnames(TEM_PRE) <- c("month", "tem", "pre")
TEM_PRE$month <- 1:12

ggplot(TEM_PRE) + geom_line(aes(x = month, y = tem), size = 1.2, colour = "#ff6666") + geom_bar(aes(x = month, y = pre/60), width = 0.8, stat = "identity", fill = "#005ac8", colour = "#005ac8") + scale_y_continuous(name = "Temperature (\u00B0C)", sec.axis = sec_axis(~ .*60, name = "Precipitation (mm)")) + scale_x_continuous(name = "Month", labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), breaks = 1:12) + theme_light() + theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14))

# ggsave("p:\\VegLab\\Weather_Station-HOBO-TMS\\Weather station data\\R\\figure\\temperature_precipitation_WS.tiff", width = 150, height = 150, units = 'mm', dpi = 200)

# Wind rose diagram ----
source("p:\\VegLab\\Weather_Station-HOBO-TMS\\Weather station data\\R\\plot_windrose.r")
library(lubridate)

WS <- read.csv("p:\\VegLab\\Weather_Station-HOBO-TMS\\Weather station data\\R\\data&code_for paper\\weather_station.csv")
WS$time <- mdy_hm(WS$time)
WS <- WS[!is.na(WS$wind_speed) & WS$wind_speed != 0 & !is.na(WS$wind_direction),]
windrose <- plot_windrose(spd = WS$wind_speed, dir = WS$wind_direction, spdres = 0.5, dirres = 10, spdmin = 0, spdmax = 6)

windrose
# ggsave("figure\\annual_windrose_2021.tiff", width = 100, height = 80, units = 'mm', dpi = 600, compression = 'lzw', scale = 1.5)

# Foggy time ----
library(lubridate)

WS <- read.csv("p:\\VegLab\\Weather_Station-HOBO-TMS\\Weather station data\\R\\data&code_for paper\\weather_station.csv")
WS$time <- mdy_hm(WS$time)

# visibility: visibility, m
# Light: 500-1000 m, medium: 100-500 m, heavy: < 100 m
FOG <- WS[WS$visibility <= 1000, c("time", "visibility")]
FOG$level <- cut(FOG$visibility, c(0, 100, 500, 1000), labels = c("H", "M", "L"))

month_min <- vector(mode = "numeric", length = 12)
month_min[c(1, 3, 5, 7, 8, 10, 12)] <- 31*24*60
month_min[c(4, 6, 9, 11)] <- 30*24*60
month_min[2] <- 28*24*60

FOG$month <- month(FOG$time)
FOG_count <- aggregate(level ~ month, data = FOG, table)
FOG_count <- cbind(FOG_count$month, FOG_count$level[,1], FOG_count$level[,2], FOG_count$level[,3])
colnames(FOG_count) <- c("month", "heavy", "medium", "light")

# percentage
FOG_count[,-1] <- (FOG_count[,-1])*10/month_min*100
FOG_count <- t(FOG_count[,-1])
colnames(FOG_count) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# length(unique(floor_date(FOG$time, "day"))) # 329

# tiff("figure\\foggy_time_2021.tiff", width = 150, height = 100, units = 'mm', res = 600, compression = 'lzw', pointsize = 10)
barplot(FOG_count[3:1,], xlab = list ("Month", cex = 1.2), ylab = list ("Foggy time (%)", cex = 1.2), las = 1, col = gray.colors(3)[3:1])
legend("top", title = 'Fog intensity', legend = c("heavy", "medium", "light"), col = gray.colors(3), pch = 15, box.lty = 0, cex = 1.2)
# dev.off()


# Temperature and relative humidity in LFDP ----
library(lubridate)
library(tidyverse)
library(tidyquant)

LFDP_1 <- read.csv("p:\\VegLab\\Weather_Station-HOBO-TMS\\Weather station data\\R\\data&code_for paper\\LFPD_1.csv")
LFDP_2 <- read.csv("p:\\VegLab\\Weather_Station-HOBO-TMS\\Weather station data\\R\\data&code_for paper\\LFPD_2.csv")
LFDP_3 <- read.csv("p:\\VegLab\\Weather_Station-HOBO-TMS\\Weather station data\\R\\data&code_for paper\\LFPD_3.csv")
LFDP_1$time <- mdy_hm(LFDP_1$time)
LFDP_2$time <- mdy_hm(LFDP_2$time)
LFDP_3$time <- mdy_hm(LFDP_3$time)

## windows average ----
width_of_window = 5
LFDP_1 <- tq_mutate(LFDP_1, select = c(temp, RH), mutate_fun = rollapply, width = width_of_window, FUN = function(y) mean(y, na.rm = T), align = 'center', col_rename = c('temp_roll', 'RH_roll'))[c(-1,-2,-35038,-35039),]

LFDP_2 <- tq_mutate(LFDP_2, select = c(temp, RH), mutate_fun = rollapply, width = width_of_window, FUN = function(y) mean(y, na.rm = T), align = 'center', col_rename = c('temp_roll', 'RH_roll'))[c(-1,-2,-35036,-35037),]

LFDP_3 <- tq_mutate(LFDP_3, select = c(temp, RH), mutate_fun = rollapply, width = width_of_window, FUN = function(y) mean(y, na.rm = T), align = 'center', col_rename = c('temp_roll', 'RH_roll'))[c(-1,-2,-35037,-35038),]

times <- intersect(LFDP_1$time, LFDP_2$time)
times <- intersect(times, LFDP_3$time)

LFDP_1 <- LFDP_1[LFDP_1$time %in% times,]
LFDP_2 <- LFDP_2[LFDP_2$time %in% times,]
LFDP_3 <- LFDP_3[LFDP_3$time %in% times,]

## Temperature mean ----
temp <- cbind.data.frame(LFDP_1$time, LFDP_1$temp_roll, LFDP_2$temp_roll, LFDP_3$temp_roll)
colnames(temp) <- c("time", "LFDP_1", "LFDP_2", "LFDP_3")
temp$Day <- floor_date(temp$time, unit = "days")

temp <- aggregate(cbind(LFDP_1, LFDP_2, LFDP_3) ~ Day, data = temp, FUN = mean)

temp$Month <- floor_date(temp$Day, unit = "months")
BKUP <- temp_devi3 <- temp
t <- aggregate(cbind(LFDP_1, LFDP_2, LFDP_3) ~ Month, data = temp, FUN = function(x) c(mean(x), sd(x)))
t <- cbind(as.character(t$Month), t$LFDP_1[,1], t$LFDP_2[,1], t$LFDP_3[,1], t$LFDP_1[,2], t$LFDP_2[,2], t$LFDP_3[,2])
colnames(t) <- c("Month", "LFDP_1M", "LFDP_2M", "LFDP_3M", "LFDP_1S", "LFDP_2S", "LFDP_3S")
temp <- as.data.frame(t)
temp$Month <- ymd(temp$Month)
temp[, 2:7] <- sapply(temp[, 2:7], FUN = as.numeric)

## Temperature deviation ----
# monthly mean daily differences (mean daily Temperature differences)
temp_devi3$LFDP_12 <- temp_devi3$LFDP_1 - temp_devi3$LFDP_2
temp_devi3$LFDP_32 <- temp_devi3$LFDP_3 - temp_devi3$LFDP_2
t <- aggregate(cbind(LFDP_12, LFDP_32) ~ Month, data = temp_devi3, FUN = function(x) c(mean(x), sd(x)))
t <- cbind(as.character(t$Month), t$LFDP_12[,1], t$LFDP_32[,1], t$LFDP_12[,2], t$LFDP_32[,2])
colnames(t) <- c("Month", "LFDP_12M", "LFDP_32M", "LFDP_12S", "LFDP_32S")
temp_devi3 <- as.data.frame(t)
temp_devi3$Month <- ymd(temp_devi3$Month)
temp_devi3[, 2:5] <- sapply(temp_devi3[, 2:5], FUN = as.numeric)


# Temperature, mean, deviation
# tiff("figure\\Temperature_mean_deviation_LFDP.tiff", width = 240, height = 240, units = 'mm', res = 600, compression = 'lzw', pointsize = 11)

par(mfrow = c(2,1))
par(mar = c(2, 4.1, 2, 2))

# mean
plot(LFDP_1M ~ rownames(temp), data = temp, ylim = c(1.5, 21), xlab = "Month", ylab = "Temperature (\u00B0C)", pch = 16, col = "#1b9e77CC", cex = 1.35, las = 1, xaxt = "n")
axis(side = 1, at = 1:24, labels = paste(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"), rep (c("21", "22", "23"), each = 12), sep = '\n')[8:31], padj = 0.3)
lines(LFDP_1M ~ rownames(temp), data = temp, col = "#1b9e77CC", lwd = 2)
for (i in 1:24){
  lines(x = c(i-0.05, i-0.05), y = c(temp$LFDP_1M[i] + temp$LFDP_1S[i], temp$LFDP_1M[i] - temp$LFDP_1S[i]), col = "#1b9e77CC", lwd = 1)
  lines(x = c(i-0.05-0.1, i-0.05+0.1), y = c(temp$LFDP_1M[i] + temp$LFDP_1S[i], temp$LFDP_1M[i] + temp$LFDP_1S[i]), col = "#1b9e77CC", lwd = 1)
  lines(x = c(i-0.05-0.1, i-0.05+0.1), y = c(temp$LFDP_1M[i] - temp$LFDP_1S[i], temp$LFDP_1M[i] - temp$LFDP_1S[i]), col = "#1b9e77CC", lwd = 1)
}

points(LFDP_2M ~ rownames(temp), data = temp, pch = 16, col = "#d95f02CC", cex = 1.35)
lines(LFDP_2M ~ rownames(temp), data = temp, col = "#d95f02CC", lwd = 2)
for (i in 1:24){
  lines(x = c(i, i), y = c(temp$LFDP_2M[i] + temp$LFDP_2S[i], temp$LFDP_2M[i] - temp$LFDP_2S[i]), col = "#d95f02CC", lwd = 1)
  lines(x = c(i-0.1, i+0.1), y = c(temp$LFDP_2M[i] + temp$LFDP_2S[i], temp$LFDP_2M[i] + temp$LFDP_2S[i]), col = "#d95f02CC", lwd = 1)
  lines(x = c(i-0.1, i+0.1), y = c(temp$LFDP_2M[i] - temp$LFDP_2S[i], temp$LFDP_2M[i] - temp$LFDP_2S[i]), col = "#d95f02CC", lwd = 1)
}

points(LFDP_3M ~ rownames(temp), data = temp, pch = 16, col = "#7570b3CC", cex = 1.35)
lines(LFDP_3M ~ rownames(temp), data = temp, col = "#7570b3CC", lwd = 2)
for (i in 1:24){
  lines(x = c(i+0.05, i+0.05), y = c(temp$LFDP_3M[i] + temp$LFDP_3S[i], temp$LFDP_3M[i] - temp$LFDP_3S[i]), col = "#7570b3CC", lwd = 1)
  lines(x = c(i+0.05-0.1, i+0.05+0.1), y = c(temp$LFDP_3M[i] + temp$LFDP_3S[i], temp$LFDP_3M[i] + temp$LFDP_3S[i]), col = "#7570b3CC", lwd = 1)
  lines(x = c(i+0.05-0.1, i+0.05+0.1), y = c(temp$LFDP_3M[i] - temp$LFDP_3S[i], temp$LFDP_3M[i] - temp$LFDP_3S[i]), col = "#7570b3CC", lwd = 1)
}

legend(x = 16.5, y = 21.5, legend = c("valley", "ridge", "windward"), col = c("#1b9e77CC", "#d95f02CC", "#7570b3CC"), pch = 16, lwd = 1.35, bty = "n", cex = 1.2)
text(x = 0.4, y = 21, "A", xpd = NA, cex = 1.5)

# deviation
par(mar = c(3, 4.1, 2, 2))
plot(rep(0, 24) ~ rownames(temp_devi3), ylim = c(-0.65, 0.4), xlab = "Month", ylab = "Temperature deviation", pch = 16, col = "#d95f02CC", cex = 1.35, las = 1, xaxt = "n")
axis(side = 1, at = 1:24, labels = paste(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"), rep (c("21", "22", "23"), each = 12), sep = '\n')[8:31], padj = 0.3)
lines(x = c(1,24), y = c(0,0), col = "#d95f02CC", lwd = 2)

points(LFDP_12M ~ rownames(temp_devi3), data = temp_devi3, pch = 16, col = "#1b9e77CC", cex = 1.35)
lines(LFDP_12M ~ rownames(temp_devi3), data = temp_devi3, col = "#1b9e77CC", lwd = 2)
for (i in 1:24){
  lines(x = c(i-0.05, i-0.05), y = c(temp_devi3$LFDP_12M[i] + temp_devi3$LFDP_12S[i], temp_devi3$LFDP_12M[i] - temp_devi3$LFDP_12S[i]), col = "#1b9e77CC", lwd = 1)
  lines(x = c(i-0.05-0.1, i-0.05+0.1), y = c(temp_devi3$LFDP_12M[i] + temp_devi3$LFDP_12S[i], temp_devi3$LFDP_12M[i] + temp_devi3$LFDP_12S[i]), col = "#1b9e77CC", lwd = 1)
  lines(x = c(i-0.05-0.1, i-0.05+0.1), y = c(temp_devi3$LFDP_12M[i] - temp_devi3$LFDP_12S[i], temp_devi3$LFDP_12M[i] - temp_devi3$LFDP_12S[i]), col = "#1b9e77CC", lwd = 1)
}

points(LFDP_32M ~ rownames(temp_devi3), data = temp_devi3, pch = 16, col = "#7570b3CC", cex = 1.35)
lines(LFDP_32M ~ rownames(temp_devi3), data = temp_devi3, col = "#7570b3CC", lwd = 2)
for (i in 1:24){
  lines(x = c(i+0.05, i+0.05), y = c(temp_devi3$LFDP_32M[i] + temp_devi3$LFDP_32S[i], temp_devi3$LFDP_32M[i] - temp_devi3$LFDP_32S[i]), col = "#7570b3CC", lwd = 1)
  lines(x = c(i+0.05-0.1, i+0.05+0.1), y = c(temp_devi3$LFDP_32M[i] + temp_devi3$LFDP_32S[i], temp_devi3$LFDP_32M[i] + temp_devi3$LFDP_32S[i]), col = "#7570b3CC", lwd = 1)
  lines(x = c(i+0.05-0.1, i+0.05+0.1), y = c(temp_devi3$LFDP_32M[i] - temp_devi3$LFDP_32S[i], temp_devi3$LFDP_32M[i] - temp_devi3$LFDP_32S[i]), col = "#7570b3CC", lwd = 1)
}

legend("topright", legend = c("reference", "valley-ridge", "windward-ridge"), col = c("#d95f02CC", "#1b9e77CC", "#7570b3CC"), pch = 16, lwd = 1.35, bty = "n", cex = 1.2)
text(x = 0.4, y = 0.4, "B", xpd = NA, cex = 1.5)

# dev.off()



## Relative humidity ----
RH <- cbind.data.frame(LFDP_1$time, LFDP_1$RH_roll, LFDP_2$RH_roll, LFDP_3$RH_roll)
colnames(RH) <- c("time", "LFDP_1", "LFDP_2", "LFDP_3")
RH$Day <- floor_date(RH$time, unit = "days")

RH <- aggregate(cbind(LFDP_1, LFDP_2, LFDP_3) ~ Day, data = RH, FUN = mean)
RH$Month <- floor_date(RH$Day, unit = "months")
t <- aggregate(cbind(LFDP_1, LFDP_2, LFDP_3) ~ Month, data = RH, FUN = function(x) c(sum(x >= 99), sum(x <= 75)))
t <- cbind.data.frame(t$Month, t$LFDP_1[,1], t$LFDP_2[,1], t$LFDP_3[,1], t$LFDP_1[,2], t$LFDP_2[,2], t$LFDP_3[,2])
colnames(t) <- c("Month", "LFDP_1W", "LFDP_2W", "LFDP_3W", "LFDP_1D", "LFDP_2D", "LFDP_3D")
RH <- t

mdays <- rep(c(31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30, 31), 2) # from Aug. to Jul.

# wet, dry days
# tiff("figure\\wet_dry_day_LFDP.tiff", width = 240, height = 240, units = 'mm', res = 600, compression = 'lzw', pointsize = 11)

par(mfrow = c(2,1))

# W
par(mar = c(2, 4.1, 2, 2))
plot(LFDP_1W/mdays ~ rownames(RH), data = RH, ylim = c(0.15,1.1), xlab = "Month", ylab = "Proportion of wet days", pch = 16, col = "#1b9e77CC", cex = 1.35, las = 1, xaxt = "n")
axis(side = 1, at = 1:24, labels = paste(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"), rep (c("21", "22", "23"), each = 12), sep = '\n')[8:31], padj = 0.3)
lines(LFDP_1W/mdays ~ rownames(RH), data = RH, col = "#1b9e77CC", lwd = 2)

points(LFDP_2W/mdays ~ rownames(RH), data = RH, pch = 16, col = "#d95f02CC", cex = 1.35)
lines(LFDP_2W/mdays ~ rownames(RH), data = RH, col = "#d95f02CC", lwd = 2)

points(LFDP_3W/mdays ~ rownames(RH), data = RH, pch = 16, col = "#7570b3CC", cex = 1.35)
lines(LFDP_3W/mdays ~ rownames(RH), data = RH, col = "#7570b3CC", lwd = 2)

legend("topright", legend = c("valley", "ridge", "windward"), col = c("#1b9e77CC", "#d95f02CC", "#7570b3CC"), pch = 16, lwd = 1.35, bty = "n", cex = 1.2)
text(x = 0.4, y = 1.1, "A", xpd = NA, cex = 1.5)

# D
par(mar = c(3, 4.1, 2, 2))
plot(LFDP_1D/mdays ~ rownames(RH), data = RH, ylim = c(0,0.23), xlab = "Month", ylab = "Proportion of dry days", pch = 16, col = "#1b9e77CC", cex = 1.35, las = 1, xaxt = "n")
axis(side = 1, at = 1:24, labels = paste(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"), rep (c("21", "22", "23"), each = 12), sep = '\n')[8:31], padj = 0.3)
lines(LFDP_1D/mdays ~ rownames(RH), data = RH, col = "#1b9e77CC", lwd = 2)

points(LFDP_2D/mdays ~ rownames(RH), data = RH, pch = 16, col = "#d95f02CC", cex = 1.35)
lines(LFDP_2D/mdays ~ rownames(RH), data = RH, col = "#d95f02CC", lwd = 2)

points(LFDP_3D/mdays ~ rownames(RH), data = RH, pch = 16, col = "#7570b3CC", cex = 1.35)
lines(LFDP_3D/mdays ~ rownames(RH), data = RH, col = "#7570b3CC", lwd = 2)

legend("top", legend = c("valley", "ridge", "windward"), col = c("#1b9e77CC", "#d95f02CC", "#7570b3CC"), pch = 16, lwd = 1.35, bty = "n", cex = 1.2)
text(x = 0.4, y = 0.23, "B", xpd = NA, cex = 1.5)

# dev.off()

# Soil moisture in LFDP ----
library(lubridate)
library(tidyquant)

TMS353 <- read.csv("p:\\VegLab\\Weather_Station-HOBO-TMS\\Weather station data\\R\\data&code_for paper\\TMS353.csv")
TMS355 <- read.csv("p:\\VegLab\\Weather_Station-HOBO-TMS\\Weather station data\\R\\data&code_for paper\\TMS355.csv")
TMS357 <- read.csv("p:\\VegLab\\Weather_Station-HOBO-TMS\\Weather station data\\R\\data&code_for paper\\TMS357.csv")
TMS359 <- read.csv("p:\\VegLab\\Weather_Station-HOBO-TMS\\Weather station data\\R\\data&code_for paper\\TMS359.csv")
TMS353$time <- mdy_hm(TMS353$time)
TMS355$time <- mdy_hm(TMS355$time)
TMS357$time <- mdy_hm(TMS357$time)
TMS359$time <- mdy_hm(TMS359$time)

## window average ----
width_of_window = 9
TMS353 <- tq_mutate(TMS353, select = soil_moisture, mutate_fun = rollapply, width = width_of_window, FUN = function(y) mean(y, na.rm = T), align = 'center', col_rename = 'soil_moisture_roll')
TMS355 <- tq_mutate(TMS355, select = soil_moisture, mutate_fun = rollapply, width = width_of_window, FUN = function(y) mean(y, na.rm = T), align = 'center', col_rename = 'soil_moisture_roll')
TMS357 <- tq_mutate(TMS357, select = soil_moisture, mutate_fun = rollapply, width = width_of_window, FUN = function(y) mean(y, na.rm = T), align = 'center', col_rename = 'soil_moisture_roll')
TMS359 <- tq_mutate(TMS359, select = soil_moisture, mutate_fun = rollapply, width = width_of_window, FUN = function(y) mean(y, na.rm = T), align = 'center', col_rename = 'soil_moisture_roll')

SM <- cbind.data.frame(TMS353$time, TMS353$soil_moisture_roll, TMS355$soil_moisture_roll, TMS357$soil_moisture_roll, TMS359$soil_moisture_roll)
SM <- SM[c(-1:-4, -2973:-2976),]
colnames(SM) <- c("Time", "TMS353", "TMS355", "TMS357", "TMS359")
rownames(SM) <- 1:nrow(SM)

SM$Day <- floor_date(SM$Time, unit = "days")
SM <- aggregate(cbind(TMS353, TMS355, TMS357, TMS359) ~ Day, data = SM, FUN = mean)

# tiff("p:\\VegLab\\Weather_Station-HOBO-TMS\\Weather station data\\R\\figure\\soil_moisture_LFDP.tiff", width = 160, height = 120, units = 'mm', res = 600, compression = 'lzw')
boxplot(SM[,2:5], ylab = "Soil moisture", las = 1)
# dev.off()

#
VWC <- SM
VWC[,2:5] <- (-0.0000000134)*(VWC[,2:5])^2 + (0.000249622)*(VWC[,2:5]) - 0.157889

# tiff("p:\\VegLab\\Weather_Station-HOBO-TMS\\Weather station data\\R\\figure\\vwc_LFDP.tiff", width = 160, height = 120, units = 'mm', res = 600, compression = 'lzw')
boxplot(VWC[,2:5], ylab = "Volumetric water content", las = 1)
# dev.off()

# tiff("p:\\VegLab\\Weather_Station-HOBO-TMS\\Weather station data\\R\\figure\\soil_moisture_vwc_LFDP.tiff", width = 150, height = 75, units = 'mm', res = 600, compression = 'lzw', pointsize = 7)

par(mfrow = c(1,2))

boxplot(SM[,2:5], ylab = "Soil moisture", las = 1)
text(x = 0.45, y = 3580, "A", xpd = NA, cex = 1.1)

boxplot(VWC[,2:5], ylab = "Volumetric water content", las = 1)
text(x = 0.45, y = 0.563, "B", xpd = NA, cex = 1.1)

# dev.off()




