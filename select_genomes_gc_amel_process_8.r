setwd('~/deming_lab/cold_random_HGT/select_genomes')

## moving average function
ma <- function(x, n){
  filter(x, rep(1 / n, n), sides = 2)}

#gc1 <- read.table('master_all_out_clean_gc1.txt.gz', sep = '\t', header = F, row.names = NULL)
#gc2 <- read.table('master_all_out_clean_gc2.txt.gz', sep = '\t', header = F, row.names = NULL)

gc1 <- read.table('master_all_out_clean_2_gc1.txt.gz', sep = '\t', header = F, row.names = NULL)
gc2 <- read.table('master_all_out_clean_2_gc2.txt.gz', sep = '\t', header = F, row.names = NULL)

## for testing only
g <- levels(gc1[,1])[2]
l = 1000
##

l = length(levels(gc1[,1]))
good_means <- matrix(ncol = 8, nrow = l)
n <- 0

for(g in levels(gc1[,1])){

  n <- n + 1
  gc1_temp <- gc1[which(gc1[,1] == g),]
  gc2_temp <- gc2[which(gc2[,1] == g),]
  
  pfam <- as.character(gc1_temp[1,4])
  group <- as.character(gc1_temp[1,3])
  strain <- as.character(gc1_temp[1,2])
  
  ## reinitialize variables so that it will break if something doesn't calculate
  
  gc1_mean <- NULL
  gc2_mean <- NULL  
  new_mean <- NA  
  compare <- NULL
  confidence <- NA
  gc1_smooth <- NA
  gc2_smooth <- NA
  
  ## calculate new values
  
  gc1_smooth <- ma(gc1_temp[,5], 10)
  gc2_smooth <- ma(gc2_temp[,5], 10)
  
  gc1_smooth <- na.omit(gc1_smooth)
  gc2_smooth <- na.omit(gc2_smooth)
    
  gc1_mean <- mean(gc1_temp[,5])
  gc2_mean <- mean(gc2_temp[,5])    
  compare <- ks.test(gc1_smooth, gc2_smooth)
  
  ## welcome to the conditional labryinth!
  
  if(gc1_mean > 0 && gc2_mean > 0){
    if(compare$p.value > 0.05){
      confidence <- 'high' ## populations are similar
      new_mean <- mean(c(gc1_mean, gc2_mean))
    }else{
      if(abs(gc1_mean - gc2_mean) < 100){
        confidence <- 'med' ## populations are different but means are within 100 My
        new_mean <- mean(c(gc1_mean, gc2_mean))
      }
    }
  }else{
    if(gc1_mean > 0){
      confidence <- 'low' ## gc1 mean is only positive
      new_mean <- gc1_mean
    }else{
      if(gc2_mean > 0){
        confidence <- 'low' ## gc2 mean is only positive
        new_mean <- gc2_mean
      }
    }
  }
  
  print(paste(n, l, gc1_mean, gc2_mean, new_mean, confidence))
  out <- c(g, strain, pfam, group, new_mean, confidence, gc1_mean, gc2_mean)
  good_means[n,] <- out
  
}

write.table(good_means, 'aquisition_means_2_v8.txt', sep = '\t', quote = F, row.names = F, col.names = F)

#### below here is work with means ####

good_means <- read.table('aquisition_means_2_v8.txt',
                         sep = '\t',
                         colClasses = c('factor', 'factor', 'factor', 'factor', 'numeric', 'factor', 'numeric', 'numeric'))

good_means <- read.table('aquisition_means_2_v8.txt',
                         sep = '\t')
good_means <- good_means[which(good_means[,6] %in% c('high', 'med', 'low')),]
good_means[,6] <- as.numeric(good_means[,6])
cold_good_means <- good_means[which(good_means[,4] == 'cold'),]
control_good_means <- good_means[which(good_means[,4] == 'control'),]

table(good_means[,6])

m <- 650 # extent of glacial record, alternatively x intercept of polynomial fit

#### plot prob range control ####

control_range_mids <- matrix(ncol = 2, nrow = length(control_good_means[,1]))

control_mean_range <- vector(length = 10 ** 6)
n <- 0

for(i in 1:length(control_good_means[,7])){
  if(control_good_means[i,7] > 0 && control_good_means[i,8] > 0){
    for(d in c(control_good_means[i,7]:control_good_means[i,8])){
      n <- n + 1
      control_mean_range[n] <- d
      print(paste(d,control_good_means[i,7], control_good_means[i,8]))
    }
    out <- c(mean(c(control_good_means[i,7], control_good_means[i,8])), length(c(control_good_means[i,7]:control_good_means[i,8])))
    control_range_mids[i,] <- out
  }else{
    if(control_good_means[i,7] > 0){
      for(d in c((control_good_means[i,7] - 5):(control_good_means[i,7] + 5))){
        n <- n + 1
        control_mean_range[n] <- d
        print(paste(d, control_good_means[i,7]))
      }
      out <- c(control_good_means[i,7], 10)
      control_range_mids[i,] <- out
  }else{
      if(control_good_means[i,8] > 0){
        for(d in c((control_good_means[i,8] - 5):(control_good_means[i,8] + 5))){
          n <- n + 1
          control_mean_range[n] <- d
          print(paste(d, control_good_means[i,8]))
        }
        out <- c(control_good_means[i,8], 10)
        control_range_mids[i,] <- out
      }
    }
  }
}

## make sure that there is no correlation between time and range, or analysis is invalid!

plot(control_range_mids[,1],
     control_range_mids[,2],
     xlab = 'midpoint',
     ylab = 'range')

abline(lm(control_range_mids[,2] ~ control_range_mids[,1]),
       col = 'red')

summary(lm(control_range_mids[,2] ~ control_range_mids[,1]))

#### plot cold mean range ####

cold_range_mids <- matrix(ncol = 2, nrow = length(cold_good_means[,1]))

cold_mean_range <- vector(length = 10 ** 6)
n <- 0

for(i in 1:length(cold_good_means[,7])){
  if(cold_good_means[i,7] > 0 && cold_good_means[i,8] > 0){
    for(d in c(cold_good_means[i,7]:cold_good_means[i,8])){
      n <- n + 1
      cold_mean_range[n] <- d
      print(paste(d,cold_good_means[i,7], cold_good_means[i,8]))
    }
    out <- c(mean(c(cold_good_means[i,7], cold_good_means[i,8])), length(c(cold_good_means[i,7]:cold_good_means[i,8])))
    cold_range_mids[i,] <- out
  }else{
    if(cold_good_means[i,7] > 0){
      for(d in c((cold_good_means[i,7] - 5):(cold_good_means[i,7] + 5))){
        n <- n + 1
        cold_mean_range[n] <- d
        print(paste(d, cold_good_means[i,7]))
      }
      out <- c(cold_good_means[i,7], 10)
      cold_range_mids[i,] <- out
    }else{
      if(cold_good_means[i,8] > 0){
        for(d in c((cold_good_means[i,8] - 5):(cold_good_means[i,8] + 5))){
          n <- n + 1
          cold_mean_range[n] <- d
          print(paste(d, cold_good_means[i,8]))
        }
        out <- c(cold_good_means[i,8], 10)
        cold_range_mids[i,] <- out
      }
    }
  }
}

## make sure that there is no correlation between time and range, or analysis is invalid!

plot(cold_range_mids[,1],
     cold_range_mids[,2],
     xlab = 'midpoint',
     ylab = 'range')

abline(lm(cold_range_mids[,2] ~ cold_range_mids[,1]),
       col = 'red')

summary(lm(cold_range_mids[,2] ~ cold_range_mids[,1]))

hist(cold_range_mids[,2], breaks = 100)

#### plot both raw ####

library(dplR)
control_mean_range_hist <- hist(control_mean_range[intersect(which(control_mean_range > 10), which(control_mean_range < m))], breaks = m)
control_mean_range_spline <- ffcsaps(control_mean_range_hist$counts, control_mean_range_hist$mids, nyrs = 0.67 * m)

cold_mean_range_hist <- hist(cold_mean_range[intersect(which(cold_mean_range > 10), which(cold_mean_range < m))], breaks = m)
cold_mean_range_spline <- ffcsaps(cold_mean_range_hist$counts, cold_mean_range_hist$mids, nyrs = 0.67 * m)

pdf('raw_date_spline.pdf', width = 8, height = 6)

plot(control_mean_range_hist$counts ~ control_mean_range_hist$mids,
     type = 'n',
     ylab = 'Relative HGT event probability',
     xlab = 'Ma',
     ylim = c(0,max(control_mean_range_hist$counts) + 10)
)

points(control_mean_range_hist$counts ~ control_mean_range_hist$mids,
      col = 'orange',
      pch = 19,
      cex = 0.6)

points(cold_mean_range_hist$counts ~ cold_mean_range_hist$mids,
      col = 'blue',
      pch = 19,
      cex = 0.6)

points(control_mean_range_spline ~ control_mean_range_hist$mids,
       type = 'l',
       col = 'orange',
       lwd = 2)

points(cold_mean_range_spline ~ cold_mean_range_hist$mids,
       type = 'l',
       col = 'blue',
       lwd = 2)

legend('topright',
       legend = c('Cold', 'Control'),
       pch = c(20,20),
       col = c('blue', 'orange'),
       bg = 'white')

dev.off()

#### plot control norm ####

control_mean_range_norm <- control_mean_range_hist$counts - control_mean_range_spline
control_mean_range_norm_smooth <- ma(control_mean_range_norm, 50)

pdf('normalized_HGT_occurrence_control.pdf', width = 8, height = 6)

plot(control_mean_range_norm ~ control_mean_range_hist$mids,
     type = 'h',
     col = 'orange',
     ylab = 'HGT event anomaly',
     xlab = 'Ma',
     ylim = c(min(control_mean_range_norm) - 10,max(control_mean_range_norm) + 10),
     xlim = c(10, max(control_mean_range_hist$mids))
)

points(control_mean_range_norm_smooth ~ control_mean_range_hist$mids,
       type = 'l',
       lwd = 2,
       col = 'black')

rect(c(1,105,270,390,635),
     c(min(control_mean_range_norm) - 10),
     c(50,220,370,475,650),
     c(min(control_mean_range_norm) - 5),
     border = F,
     col = 'grey')


lines(c(65,65),
      c(min(control_mean_range_norm) - 10,max(control_mean_range_norm) + 10),
      lty = 2)

lines(c(200,200),
      c(min(control_mean_range_norm) - 10,max(control_mean_range_norm) + 10),
      lty = 2)

lines(c(250,250),
      c(min(control_mean_range_norm) - 10,max(control_mean_range_norm) + 10),
      lty = 2)

lines(c(375,375),
      c(min(control_mean_range_norm) - 10,max(control_mean_range_norm) + 10),
      lty = 2)

lines(c(450,450),
      c(min(control_mean_range_norm) - 10,max(control_mean_range_norm) + 10),
      lty = 2)

dev.off()

#### plot cold norm ####

cold_mean_range_norm <- cold_mean_range_hist$counts - cold_mean_range_spline
cold_mean_range_norm_smooth <- ma(cold_mean_range_norm, 50)

pdf('normalized_HGT_occurrence_cold.pdf', width = 8, height = 6)

plot(cold_mean_range_norm ~ cold_mean_range_hist$mids,
     type = 'h',
     col = 'blue',
     ylab = 'HGT event anomaly',
     xlab = 'Ma',
     ylim = c(min(cold_mean_range_norm) - 10,max(cold_mean_range_norm) + 10),
     xlim = c(10, max(cold_mean_range_hist$mids))
)

points(cold_mean_range_norm_smooth ~ cold_mean_range_hist$mids,
       type = 'l',
       lwd = 2,
       col = 'black')

rect(c(1,105,270,390,635),
     c(min(cold_mean_range_norm) - 10),
     c(50,220,370,475,650),
     c(min(cold_mean_range_norm) - 5),
     border = F,
     col = 'grey')


lines(c(65,65),
      c(min(cold_mean_range_norm) - 10,max(cold_mean_range_norm) + 10),
      lty = 2)

lines(c(200,200),
      c(min(cold_mean_range_norm) - 10,max(cold_mean_range_norm) + 10),
      lty = 2)

lines(c(250,250),
      c(min(cold_mean_range_norm) - 10,max(cold_mean_range_norm) + 10),
      lty = 2)

lines(c(375,375),
      c(min(cold_mean_range_norm) - 10,max(cold_mean_range_norm) + 10),
      lty = 2)

lines(c(450,450),
      c(min(cold_mean_range_norm) - 10,max(cold_mean_range_norm) + 10),
      lty = 2)

dev.off()

#### analyze timescale ####

glacial <- c(1:50,105:220,270:370,390:475, 635:650) # Veizer et al. 2000
f_glacial <- length(glacial) / m # fraction of time that is glacial

cold_glacial_mids <- which(ceiling(cold_mean_range_hist$mids) %in% glacial)
glacial_hgt_norm_cold <- sum(cold_mean_range_norm[cold_glacial_mids])

control_glacial_mids <- which(ceiling(control_mean_range_hist$mids) %in% glacial)
glacial_hgt_norm_control <- sum(control_mean_range_norm[control_glacial_mids])

#### MC simulation ####

## cold

s <- 1000000 ## number of simulations to run
g <- 0
mcs_cold <- vector("numeric",length = s)

for(i in seq(1,s)){
  rn <- sample(1:length(cold_mean_range_norm), length(cold_glacial_mids))
  mc_glacial_hgt_norm_cold <- sum(cold_mean_range_norm[rn])
  if(mc_glacial_hgt_norm_cold > glacial_hgt_norm_cold){ g <- g + 1 }
  mcs_cold[i] <- mc_glacial_hgt_norm_cold
  print(paste(i, mc_glacial_hgt_norm_cold))
}

print(paste(c('cold', g, s, g / s)))

pdf('cold_monte_carlo.pdf')

hist_mcs_cold <- hist(mcs_cold, breaks = 100,
                         xlab = 'HGT anomaly during cold period',
                         main = 'cold',
                         col = 'black',
                         cex.axis = 1,
                         cex.lab = 1)

lines(c(glacial_hgt_norm_cold, glacial_hgt_norm_cold), c(0, 100000), col = 'orange', lwd = 2)

box()

dev.off()

## control

s <- 1000000 ## number of simulations to run
g <- 0
mcs_control <- vector("numeric",length = s)

for(i in seq(1,s)){
  rn <- sample(1:length(control_mean_range_norm), length(control_glacial_mids))
  mc_glacial_hgt_norm_control <- sum(control_mean_range_norm[rn])
  if(mc_glacial_hgt_norm_control > glacial_hgt_norm_control){ g <- g + 1 }
  mcs_control[i] <- mc_glacial_hgt_norm_control
  print(paste(i, mc_glacial_hgt_norm_control))
}

print(paste(c('control', g, s, g / s)))

pdf('control_monte_carlo.pdf')

hist_mcs_control <- hist(mcs_control, breaks = 100,
                      xlab = 'HGT anomaly during cold period',
                      main = 'control',
                      col = 'black',
                      cex.axis = 1,
                      cex.lab = 1)

lines(c(glacial_hgt_norm_control, glacial_hgt_norm_control), c(0, 100000), col = 'orange', lwd = 2)

box()

dev.off()

#### explore specific region ####

temp_period <- control_good_means[intersect(which(control_good_means[,5] < 150), which(control_good_means[,5] > 100)),]
