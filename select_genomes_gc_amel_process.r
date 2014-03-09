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

numeric_means_control <- as.numeric(aquisition_dates_control[,5])
numeric_errors_control <- as.numeric(aquisition_dates_control[,8])

numeric_means_cold <- as.numeric(aquisition_dates_cold[,5])
numeric_errors_cold <- as.numeric(aquisition_dates_cold[,8])

glacial <- c(1:50,105:220,270:370,390:475, 635:650) # Veizer et al. 2000
#glacial <- c(1:55, 105:170, 250:330, 418:460, 635:650)
glacial_hgt_cold <- numeric_means_cold[which(ceiling(numeric_means_cold) %in% glacial)]
glacial_hgt_control <- numeric_means_control[which(ceiling(numeric_means_control) %in% glacial)]

m <- 650 # extent of glacial record, alternatively x intercept of polynomial fit
f_glacial <- length(glacial) / m # fraction of time that is glacial
#f_glacial_hgt_cold <- length(glacial_hgt) / length(numeric_means[which(numeric_means_cold < m)])

## bin data and export for model fit in pydavis

master_hist_cold <- hist(numeric_means_cold[which(numeric_means_cold < m)], breaks = seq(0,m + 5,10))
master_hist_control <- hist(numeric_means_control[which(numeric_means_control < m)], breaks = seq(0,m + 5,10))

#write.table(data.frame(master_hist$counts, master_hist$mids), 'master_dates_histogram.txt', quote = F, sep = '\t', col.names = F, row.names = F)

##### correcting for more glaciation in recent history - when HGT is more detectable #####

## internal detrend in place of external model fit using pydavis ##

library('pracma')

master_hist_norm_cold <- pracma::detrend(master_hist_cold$counts[which(master_hist_cold$mids <= m)], tt = 'linear')
master_hist_lm_cold <- lm(master_hist_cold$counts[which(master_hist_cold$mids <= m)] ~ master_hist_cold$mids[which(master_hist_cold$mids <= m)])
master_hist_anom_cold <- master_hist_norm_cold[which(master_hist_norm_cold > 0)]
master_hist_anom_time_cold <- master_hist_cold$mids[which(master_hist_norm_cold > 0)]

master_hist_norm_control <- pracma::detrend(master_hist_control$counts[which(master_hist_control$mids <= m)], tt = 'linear')
master_hist_lm_control <- lm(master_hist_control$counts[which(master_hist_control$mids <= m)] ~ master_hist_control$mids[which(master_hist_control$mids <= m)])
master_hist_anom_control <- master_hist_norm_control[which(master_hist_norm_control > 0)]
master_hist_anom_time_control <- master_hist_control$mids[which(master_hist_norm_control > 0)]

bin_comp_col <- colorRampPalette(c('black', 'violet', 'blue', 'green', 'yellow', 'orange', 'red'))(length(master_hist_norm_cold))

plot(master_hist_norm_cold ~ master_hist_norm_control,
     pch = 19,
     col = bin_comp_col)

summary(lm(master_hist_norm_cold ~ master_hist_norm_control))

pdf('not_normalized_HGT_occurrence_cold.pdf', width = 8, height = 6)

plot(master_hist_cold$counts ~ master_hist_cold$mids[which(master_hist_cold$mids <= m)],
     type = 'n',
     ylab = 'HGT events',
     xlab = 'Ma',
     ylim = c(0,max(master_hist_cold$counts) + 10)
)

rect(c(1,105,270,390,635),
     c(0),
     c(50,220,370,475,650),
     c(max(master_hist_cold$counts) + 10),
     border = F,
     col = 'grey')

points(master_hist_cold$counts ~ master_hist_cold$mids[which(master_hist_cold$mids <= m)], pch = 20)
abline(master_hist_lm_cold, col = 'red')

dev.off()

## plot cold normalized HGT occurrence

pdf('normalized_HGT_occurrence_cold.pdf', width = 8, height = 6)

plot(master_hist_norm_cold ~ master_hist_cold$mids[which(master_hist_cold$mids <= m)],
     type = 'n',
     ylim = c(min(master_hist_norm_cold) - 1,max(master_hist_norm_cold) + 1),
     ylab = 'Normalized HGT events',
     xlab = 'Ma',
     xlim = c(0,m)
     )

rect(c(1,105,270,390,635),
     c(min(master_hist_norm_cold) - 1),
     c(50,220,370,475,650),
     c(max(master_hist_norm_cold) + 1),
     border = F,
     col = 'grey')

barplot(master_hist_norm_cold[,1],
        width = 10,
        space = 0,
        add = T,
        col = 'black'
        )

lines(sepkoski$p.Lmy * 100 ~ sepkoski$Date,
     type = 'l',
      col = 'red')

dev.off()

## plot control normalized HGT occurrence

pdf('normalized_HGT_occurrence_control.pdf', width = 8, height = 6)

plot(master_hist_norm_control ~ master_hist_control$mids[which(master_hist_control$mids <= m)],
     type = 'n',
     ylim = c(min(master_hist_norm_control) - 1,max(master_hist_norm_control) + 1),
     ylab = 'Normalized HGT events',
     xlab = 'Ma',
     xlim = c(0,m)
)

rect(c(1,105,270,390,635),
     c(min(master_hist_norm_control) - 1),
     c(50,220,370,475,650),
     c(max(master_hist_norm_control) + 1),
     border = F,
     col = 'grey')

barplot(master_hist_norm_control[,1],
        width = 10,
        space = 0,
        add = T,
        col = 'black'
)

lines(sepkoski$p.Lmy * 100 ~ sepkoski$Date,
      type = 'l',
      col = 'red')

dev.off()

## difference between these two

pdf('diff_in_normalized_HGT_occurrence.pdf', width = 8, height = 6)

plot(master_hist_norm_control ~ master_hist_control$mids[which(master_hist_control$mids <= m)],
     type = 'n',
     ylim = c(min(master_hist_norm_control) - 1,max(master_hist_norm_control) + 1),
     ylab = 'Cold - control normalized HGT events',
     xlab = 'Ma',
     xlim = c(0,m)
)

rect(c(1,105,270,390,635),
     c(min(master_hist_norm_control) - 1),
     c(50,220,370,475,650),
     c(max(master_hist_norm_control) + 1),
     border = F,
     col = 'grey')

barplot(master_hist_norm_cold[,1] - master_hist_norm_control[,1],
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
glacial_hgt_norm_cold <- sum(master_hist_anom_cold[(ceiling(master_hist_anom_time_cold) %in% glacial)])
f_glacial_hgt_norm_cold <- glacial_hgt_norm_cold / sum(master_hist_anom_cold, na.rm = T)

f_glacial <- length(glacial) / m 
glacial_hgt_norm_control <- sum(master_hist_anom_control[(ceiling(master_hist_anom_time_control) %in% glacial)])
f_glacial_hgt_norm_control <- glacial_hgt_norm_control / sum(master_hist_anom_control, na.rm = T)

#### look for correlations with sepkoski curve ####

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
mcs_cold <- vector("numeric",length = n)

master_hist_mids_cold <- master_hist_anom_time_cold # time periods with positive anomalies
glacial_mids_cold <- which(master_hist_anom_time_cold %in% glacial)
master_hist_anom_sum_cold <- sum(master_hist_anom_cold)

for(i in seq(1,n)){
  print(i)
  rn <- sample(master_hist_anom_cold, length(master_hist_anom_cold))
  mc_glacial_hgt_norm_cold <- sum(rn[glacial_mids_cold])
  f_mc_glacial_hgt_norm_cold <- mc_glacial_hgt_norm_cold / master_hist_anom_sum_cold
  if(f_mc_glacial_hgt_norm_cold < f_glacial_hgt_norm_cold){ l <- l + 1 } else { g <- g + 1 }
  mcs_cold[i] <- f_mc_glacial_hgt_norm_cold
}

hist_mcs_cold <- hist(mcs_cold, breaks = 100,
                 xlab = 'Fraction of events occurring during cold period',
                 main = 'cold',
                 col = 'black',
                 cex.axis = 1,
                 cex.lab = 1)

lines(c(f_glacial_hgt_norm_cold, f_glacial_hgt_norm_cold), c(0, 100000), col = 'orange', lwd = 2)

box()

print(paste(c('cold', g, n, g / n)))


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
