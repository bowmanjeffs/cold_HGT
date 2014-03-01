setwd('~/deming_lab/cold_random_HGT/select_genomes')

gc1 <- read.table('master_all_out_clean_gc1.txt.gz', sep = '\t', header = F, row.names = NULL)
gc2 <- read.table('master_all_out_clean_gc2.txt.gz', sep = '\t', header = F, row.names = NULL)

## for testing only
g <- levels(gc1[,1])[2]
##

l = length(levels(gc1[,1]))

good_means <- matrix(ncol = 14, nrow = l)
library(nortest)
n <- 0

for(g in levels(gc1[,1])){
  n <- n + 1
  print(paste(n, 'out of', l))
  
  ##gc1
  
  r <- which(gc1[,1] == g)
  gc1_mins <- gc1[which(gc1[r,5] == 1),]
  gc1_maxs <- gc1[which(gc1[r,5] == 1000),]
  gc1_temp <- gc1[r,]
  pfam <- gc1_temp[1,4]
  group <- gc1_temp[1,3]
  strain <- gc1_temp[1,2]
  #gc1_temp <- gc1_temp[which(gc1_temp[,5] != 1000),]
  gc1_temp <- gc1_temp[which(gc1_temp[,5] != 1),]
  
  gc1_mean <- mean(gc1_temp[,5])
  gc1_sd <- sd(gc1_temp[,5])
  gc1_n <- length(gc1_temp[,5])
  error <- qnorm(0.975) * gc1_sd/sqrt(gc1_n)
  gc1_lower <- gc1_mean - error
  gc1_upper <- gc1_mean + error
  
  ##gc2
  
  r <- which(gc2[,1] == g)
  gc2_mins <- gc2[which(gc2[r,5] == 1),]
  gc2_maxs <- gc2[which(gc2[r,5] == 1000),]
  gc2_temp <- gc2[r,]
  #gc2_temp <- gc2_temp[which(gc2_temp[,2] != 1000),]
  gc2_temp <- gc2_temp[which(gc2_temp[,5] != 1),]
  
  gc2_mean <- mean(gc2_temp[,5])
  gc2_sd <- sd(gc2_temp[,5])
  gc2_n <- length(gc2_temp[,5])
  error <- qnorm(0.975) * gc2_sd/sqrt(gc2_n)
  gc2_lower <- gc2_mean - error
  gc2_upper <- gc2_mean + error
  
  ##both
  
  both_raw_mean <- mean(c(gc1[r,5], gc2[r,5]))
  both_raw_sd <- sd(c(gc1[r,5], gc2[r,5]))
  both_raw_n <- length(r)
  both_raw_error <- qnorm(0.975) * both_raw_sd/sqrt(both_raw_n)
  both_norm <- NULL

  try(both_norm <- ad.test(c(gc1_temp[,5], gc2_temp[,5])))
  
  try(if(both_norm$p.value < 0.05){
    both_mean <- mean(c(gc1_temp[,5], gc2_temp[,5]))
    both_sd <- sd(c(gc1_temp[,5], gc2_temp[,5]))
    both_n <- gc1_n + gc2_n
    both_error <- qnorm(0.975) * both_sd/sqrt(both_n)
    both_mins <- length(gc1_mins[,1]) + length(gc2_mins[,1])
    both_maxs <- length(gc1_maxs[,1]) + length(gc2_maxs[,1])
    both_norm_mean <- both_mean * (1 - (length(gc1_mins[,1]) + length(gc2_mins[,1]))/2000)
    out <- c(g, paste(strain), paste(pfam), paste(group), both_norm_mean, both_mean, both_raw_mean, both_error, both_raw_error, both_sd, both_n, both_mins, both_maxs, both_norm$p.value)
    #good_means <- rbind(good_means, out)
    good_means[n,] <- out
    print(paste(strain, pfam, both_norm_mean, both_mean, both_raw_mean, both_error, both_raw_error, both_mins), sep = ' ')
  }, silent = F
  )
   
}

plot(good_means[,5] ~ good_means[,7],
     xlab = 'raw mean',
     ylab = 'normalized mean'
     )

plot(good_means[,5] ~ good_means[,6],
     xlab = 'mean no ones',
     ylab = 'normalized mean'
)

write.table(good_means, 'aquisition_means.txt', sep = '\t', quote = F, row.names = F, col.names = F)

##### Analyze timescale #####

sepkoski <- read.table('sepkoski_curve.txt', header = T, fill = T, sep = '\t')

aquisition_dates <- read.table('aquisition_means.txt', sep = '\t')
aquisition_dates_cold <- aquisition_dates[which(aquisition_dates[,4] == 'cold'),]
aquisition_dates_control <- aquisition_dates[which(aquisition_dates[,4] == 'control'),]

## need to change these two lines to reflect the population you want to analyze
numeric_means <- as.numeric(aquisition_dates_control[,5])
numeric_errors <- as.numeric(aquisition_dates_control[,8])
##

glacial <- c(1:50,105:220,270:370,390:475, 635:650) # Veizer et al. 2000
glacial_hgt <- numeric_means[which(ceiling(numeric_means) %in% glacial)]

m <- 650 # extent of glacial record, alternatively x intercept of polynomial fit
f_glacial <- length(glacial) / m # fraction of time that is glacial
f_glacial_hgt <- length(glacial_hgt) / length(numeric_means[which(numeric_means < m)])

## bin data and export for model fit in pydavis

master_hist <- hist(numeric_means[which(numeric_means < m)], breaks = seq(0,m + 5,10))
write.table(data.frame(master_hist$counts, master_hist$mids), 'master_dates_histogram.txt', quote = F, sep = '\t', col.names = F, row.names = F)

##### correcting for more glaciation in recent history - when HGT is more detectable #####

## internal detrend in place of external model fit using pydavis ##

library('pracma')

master_hist_norm <- pracma::detrend(master_hist$counts[which(master_hist$mids <= m)], tt = 'linear')
master_hist_lm <- lm(master_hist$counts[which(master_hist$mids <= m)] ~ master_hist$mids[which(master_hist$mids <= m)])
master_hist_anom <- master_hist_norm[which(master_hist_norm > 0)]
master_hist_anom_time <- master_hist$mids[which(master_hist_norm > 0)]

pdf('not_normalized_HGT_occurrence.pdf', width = 8, height = 6)

plot(master_hist$counts ~ master_hist$mids[which(master_hist$mids <= m)],
     type = 'n',
     ylab = 'HGT events',
     xlab = 'Ma',
     ylim = c(0,max(master_hist$counts) + 10)
)

rect(c(1,105,270,390,635),
     c(0),
     c(50,220,370,475,650),
     c(max(master_hist$counts) + 10),
     border = F,
     col = 'grey')

points(master_hist$counts ~ master_hist$mids[which(master_hist$mids <= m)], pch = 20)
abline(master_hist_lm, col = 'red')

dev.off()

## plot normalized HGT occurrence

pdf('normalized_HGT_occurrence.pdf', width = 8, height = 6)

plot(master_hist_norm ~ master_hist$mids[which(master_hist$mids <= m)],
     type = 'n',
     ylim = c(min(master_hist_norm) - 1,max(master_hist_norm) + 1),
     ylab = 'Normalized HGT events',
     xlab = 'Ma',
     xlim = c(0,m)
     )

rect(c(1,105,270,390,635),
     c(min(master_hist_norm) - 1),
     c(50,220,370,475,650),
     c(max(master_hist_norm) + 1),
     border = F,
     col = 'grey')

barplot(master_hist_norm[,1],
        width = 10,
        space = 0,
        add = T,
        col = 'black'
        )

lines(sepkoski$p.Lmy * 100 ~ sepkoski$Date,
     type = 'l',
      col = 'red')

dev.off()

## plotting complete

f_glacial <- length(glacial) / m 
glacial_hgt_norm <- sum(master_hist_anom[(ceiling(master_hist_anom_time) %in% glacial)])
f_glacial_hgt_norm <- glacial_hgt_norm / sum(master_hist_anom, na.rm = T)

sepkoski_spline <- spline(sepkoski$Date, sepkoski$q.Lmy, xout = master_hist$mids[which(master_hist$mids <= max(sepkoski$Date))])
## species origination = p.Lmy according to Foote, 2000

plot(master_hist_norm[1:length(sepkoski_spline$y)] ~ sepkoski_spline$y)

sepkoski_lm <- lm(master_hist_norm[1:length(sepkoski_spline$y)] ~ sepkoski_spline$y)

abline(sepkoski_lm)

summary(sepkoski_lm)

#### MC simulation ####

n <- 1000000 ## number of simulations to run
l <- 0
g <- 0
mcs <- vector("numeric",length = n)

master_hist_mids <- master_hist_anom_time
glacial_mids <- which(master_hist_anom_time %in% glacial)
master_hist_anom_sum <- sum(master_hist_anom)

for(i in seq(1,n)){
  print(i)
  rn <- sample(master_hist_anom, length(master_hist_anom))
  mc_glacial_hgt_norm <- sum(rn[glacial_mids])
  f_mc_glacial_hgt_norm <- mc_glacial_hgt_norm / master_hist_anom_sum
  if(f_mc_glacial_hgt_norm < f_glacial_hgt_norm){ l <- l + 1 } else { g <- g + 1 }
  mcs[i] <- f_mc_glacial_hgt_norm
}

hist_mcs <- hist(mcs, breaks = 100,
                 xlab = 'Fraction of events occurring during cold period',
                 main = NULL,
                 col = 'black',
                 cex.axis = 1,
                 cex.lab = 1)

lines(c(f_glacial_hgt_norm, f_glacial_hgt_norm), c(0, 100000), col = 'orange', lwd = 2)

box()

print(g / n)


#### which HGT pfams occur more during cold period than during warm periods? ####

aquisition_dates_cold <- aquisition_dates[which(aquisition_dates[,4] == 'cold'),]

pfam_glacial <- NULL
pfam_warm <- NULL

for(p in unique(as.character(aquisition_dates_cold$V3))){
  temp <- aquisition_dates_cold[which(aquisition_dates_cold$V3 == p),]
  temp_glacial <- length(which(round(temp$V5) %in% glacial))
  temp_total <- length(temp$V3)
  temp_warm <- temp_total - temp_glacial
  pfam_glacial <- append(pfam_glacial, temp_glacial)
  pfam_warm <- append(pfam_warm, temp_warm)
  print(c(p, temp_glacial, temp_warm))
}

cold_pfam_glacial_warm <- data.frame(unique(as.character(aquisition_dates_cold$V3)), pfam_glacial, pfam_warm)

cold_pfam_glacial_warm <- cold_pfam_glacial_warm[order(cold_pfam_glacial_warm[,2] - cold_pfam_glacial_warm[,3], decreasing = T),]

write.table(cold_pfam_glacial_warm[1:50,], 'cold_pfam_glacial_warm.txt', row.names = F, col.names = F, quote = F, sep = '\t')

## go through same exercise for control group, at least for those in top 50 of cold

aquisition_dates_control <- aquisition_dates[which(aquisition_dates[,4] == 'control'),]

pfam_glacial <- NULL
pfam_warm <- NULL

for(p in cold_pfam_glacial_warm[1:50,1]){
  temp <- aquisition_dates_control[which(aquisition_dates_control$V3 == p),]
  temp_glacial <- length(which(round(temp$V5) %in% glacial))
  temp_total <- length(temp$V3)
  temp_warm <- temp_total - temp_glacial
  pfam_glacial <- append(pfam_glacial, temp_glacial)
  pfam_warm <- append(pfam_warm, temp_warm)
  print(c(p, temp_glacial, temp_warm))
}

control_pfam_glacial_warm <- data.frame(cold_pfam_glacial_warm[1:50,1], pfam_glacial, pfam_warm)

write.table(control_pfam_glacial_warm, 'control_pfam_glacial_warm.txt', row.names = F, col.names = F, quote = F, sep = '\t')
