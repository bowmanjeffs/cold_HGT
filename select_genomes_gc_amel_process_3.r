setwd('~/deming_lab/cold_random_HGT/select_genomes')

## moving average function
ma <- function(x, n){
  filter(x, rep(1 / n, n), sides = 2)}

## establish a temp dataset
get_temp <- function(level, table){
  temp <- table[which(table[,1] == level),]
  mins <- length(which(temp[,5] == 1))
  temp_reduced <- temp[which(temp[,5] != 1),]
  if(length(temp_reduced[,5]) > 500){
    temp_reduced <- temp_reduced[1:(length(temp_reduced[,5]) - mins),]
  }else temp_reduced <- temp_reduced[1,]
  #return(temp_reduced)
  #print(length(temp_reduced[,5]))
}

## if distribution is approximately normal determine mean
get_mean <- function(table){
  gamma <- fitdist(table[,5], 'gamma', method = 'mle')
  norm <- fitdist(table[,5], 'norm', method = 'mle')
  if(gamma$loglik < norm$loglik){
    table_mean <- mean(table[,5])
  }else table_mean <- NULL
  return(table_mean)
}

gc1 <- read.table('master_all_out_clean_gc1.txt.gz', sep = '\t', header = F, row.names = NULL)
gc2 <- read.table('master_all_out_clean_gc2.txt.gz', sep = '\t', header = F, row.names = NULL)

library(nortest)

## for testing only
g <- levels(gc1[,1])[2]
l = 1000
##

l = length(levels(gc1[,1]))
good_means <- matrix(ncol = 8, nrow = l)
n <- 0
k <- 0

for(g in levels(gc1[,1])){
  n <- n + 1
  gc1_temp <- gc1[which(gc1[,1] == g),]
  gc2_temp <- gc2[which(gc2[,1] == g),]
  
  pfam <- as.character(gc1_temp[1,4])
  group <- as.character(gc1_temp[1,3])
  strain <- as.character(gc1_temp[1,2])
  
  ## reinitialize variables so that it will break if something doesn't calculate
  
  gc1_norm <- list(p.value = 0)
  gc2_norm <- list(p.value = 0)
  
  gc1_smooth <- NULL
  gc2_smooth <- NULL
  
  gc1_mean <- NULL
  gc2_mean <- NULL  
  new_mean <- NULL
  
  compare <- NULL
  
  gc1_smooth <- ma(gc1_temp[,5], 10)
  gc2_smooth <- ma(gc2_temp[,5], 10)
  
  gc1_smooth <- na.omit(gc1_smooth)
  gc2_smooth <- na.omit(gc2_smooth)
  
  try(gc1_norm <- ad.test(gc1_smooth))
  try(gc2_norm <- ad.test(gc2_smooth))
  
  if(is.na(gc1_norm$p.value) == T){
    gc1_norm <- list(p.value = 0)
  }
  
  if(is.na(gc2_norm$p.value) == T){
    gc2_norm <- list(p.value = 0)
  }
  
  compare <- ks.test(gc1_smooth, gc2_smooth)
  
  ## welcome to the conditional labryinth!
  
  if(gc1_norm$p.value > 0.05 && gc2_norm$p.value > 0.05){
    if(compare$p.value > 0.05){
      confidence <- 'high' # distributions are normal and populations are similar
      gc1_mean <- mean(gc1_smooth)
      gc2_mean <- mean(gc2_smooth)
      new_mean <- mean(c(gc1_mean, gc2_mean))
    }else{
      confidence <- 'low' # distributions are normal but populations are dissimilar
      gc1_mean <- mean(gc1_smooth)
      gc2_mean <- mean(gc2_smooth)
      new_mean <- mean(c(gc1_mean, gc2_mean))
    }
  }
  
  if(gc1_norm$p.value <= 0.05 && gc2_norm$p.value > 0.05){
    if(compare$p.value > 0.05){
      confidence <- 'med-high' # distributions are similar, only one is normal
      gc1_mean <- mean(gc1_smooth)
      gc2_mean <- mean(gc2_smooth)
      new_mean <- mean(c(gc1_mean, gc2_mean))
    }else{
      confidence <- 'med' # distributions aren't similar and only one is normal
      gc1_mean <- NA
      gc2_mean <- mean(gc2_smooth)
      new_mean <- gc2_mean
    }  
  }
  
  if(gc1_norm$p.value > 0.05 && gc2_norm$p.value <= 0.05){
    if(compare$p.value > 0.05){
      confidence <- 'med-high' # distributions are similar, only one is normal
      gc1_mean <- mean(gc1_smooth)
      gc2_mean <- mean(gc2_smooth)
      new_mean <- mean(c(gc1_mean, gc2_mean))
    }else{
      confidence <- 'med' # distributions aren't similar and only one is normal
      gc2_mean <- NA
      gc1_mean <- mean(gc1_smooth)
      new_mean <- gc1_mean
    }
  }
  
  if(gc1_norm$p.value <= 0.05 && gc2_norm$p.value <= 0.05){
    if(compare$p.value > 0.05){
      confidence <- 'med' # neither distribution is normal but populations are similar
      gc1_mean <- mean(gc1_smooth)
      gc2_mean <- mean(gc2_smooth)
      new_mean <- mean(c(gc1_mean, gc2_mean))
    }else{
      confidence <- 'very-low' # neither distribution is normal and populations aren't similar
      gc2_mean <- mean(gc2_smooth)
      gc1_mean <- mean(gc1_smooth)
      pdf(paste0('gcamel_hist\\', n, '_dist.pdf'))
      hist(gc1_smooth, breaks = 1000, main = paste('gc1', gc1_norm$p.value), sub = g)
      box()
      hist(gc2_smooth, breaks = 1000, main = paste('gc2', gc2_norm$p.value), sub = g)
      box()
      dev.off()
      try(new_mean <- mean(cbind(gc1_smooth, gc2_smooth)[,which.max(c(gc1_norm$p.value, gc2_norm$p.value))]))
    }
  }
  
  print(paste(n, l, gc1_mean, gc2_mean, new_mean, confidence))
  out <- c(g, strain, pfam, group, new_mean, confidence, gc1_mean, gc2_mean)
  good_means[n,] <- out
  
}

write.table(good_means, 'aquisition_means.txt', sep = '\t', quote = F, row.names = F, col.names = F)

table(good_means[,6])
hist(as.numeric(good_means[,5]), breaks = 100)

#### analyze timescale ####

good_means <- read.table('aquisition_means.txt', sep = '\t')

positive_good_means <- good_means[which(as.numeric(good_means[,5]) > 0),]
cold_good_means <- positive_good_means[which(positive_good_means[,4] == 'cold'),]
control_good_means <- positive_good_means[which(positive_good_means[,4] == 'control'),]
cold_good_dates <- as.numeric(cold_good_means[,5])
control_good_dates <- as.numeric(control_good_means[,5])

glacial <- c(1:50,105:220,270:370,390:475, 635:650) # Veizer et al. 2000
m <- 650 # extent of glacial record, alternatively x intercept of polynomial fit
f_glacial <- length(glacial) / m # fraction of time that is glacial

master_hist_cold <- hist(cold_good_dates[which(cold_good_dates < m)], breaks = seq(0,m + 5,10))
master_hist_control <- hist(control_good_dates[which(control_good_dates < m)], breaks = seq(0,m + 5,10))

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
abline(lm(master_hist_cold$counts ~ master_hist_cold$mids[which(master_hist_cold$mids <= m)]), col = 'red')

dev.off()

##### correcting for more glaciation in recent history - when HGT is more detectable #####

library('pracma')

master_hist_norm_cold <- pracma::detrend(master_hist_cold$counts[which(master_hist_cold$mids <= m)], tt = 'linear', bp = which(master_hist_cold$mids %in% c(105, 150, 215, m)))
master_hist_lm_cold <- lm(master_hist_cold$counts[which(master_hist_cold$mids <= m)] ~ master_hist_cold$mids[which(master_hist_cold$mids <= m)])
master_hist_anom_cold <- master_hist_norm_cold[which(master_hist_norm_cold > 0)]
master_hist_anom_time_cold <- master_hist_cold$mids[which(master_hist_norm_cold > 0)]

master_hist_norm_control <- pracma::detrend(master_hist_control$counts[which(master_hist_control$mids <= m)], tt = 'linear')
master_hist_lm_control <- lm(master_hist_control$counts[which(master_hist_control$mids <= m)] ~ master_hist_control$mids[which(master_hist_control$mids <= m)])
master_hist_anom_control <- master_hist_norm_control[which(master_hist_norm_control > 0)]
master_hist_anom_time_control <- master_hist_control$mids[which(master_hist_norm_control > 0)]

## plot normalized cold HGT events

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

dev.off()

## write cold normalized HGT occurrence

write.table(cbind(master_hist_cold$mids, master_hist_norm_cold[,1]), 'cold_normalized.txt', quote = F, row.names = F, sep = '\t')

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

dev.off()

## plotting complete

f_glacial <- length(glacial) / m 
glacial_hgt_norm_cold <- sum(master_hist_anom_cold[(ceiling(master_hist_anom_time_cold) %in% glacial)])
f_glacial_hgt_norm_cold <- glacial_hgt_norm_cold / sum(master_hist_anom_cold, na.rm = T)

f_glacial <- length(glacial) / m 
glacial_hgt_norm_control <- sum(master_hist_anom_control[(ceiling(master_hist_anom_time_control) %in% glacial)])
f_glacial_hgt_norm_control <- glacial_hgt_norm_control / sum(master_hist_anom_control, na.rm = T)

#### MC simulation ####

## cold

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

## control

n <- 1000000 ## number of simulations to run
l <- 0
g <- 0
mcs_control <- vector("numeric",length = n)

master_hist_mids_control <- master_hist_anom_time_control # time periods with positive anomalies
glacial_mids_control <- which(master_hist_anom_time_control %in% glacial)
master_hist_anom_sum_control <- sum(master_hist_anom_control)

for(i in seq(1,n)){
  print(i)
  rn <- sample(master_hist_anom_control, length(master_hist_anom_control))
  mc_glacial_hgt_norm_control <- sum(rn[glacial_mids_control])
  f_mc_glacial_hgt_norm_control <- mc_glacial_hgt_norm_control / master_hist_anom_sum_control
  if(f_mc_glacial_hgt_norm_control < f_glacial_hgt_norm_control){ l <- l + 1 } else { g <- g + 1 }
  mcs_control[i] <- f_mc_glacial_hgt_norm_control
}

hist_mcs_control <- hist(mcs_control, breaks = 100,
                         xlab = 'Fraction of events occurring during control period',
                         main = 'control',
                         col = 'black',
                         cex.axis = 1,
                         cex.lab = 1)

lines(c(f_glacial_hgt_norm_control, f_glacial_hgt_norm_control), c(0, 100000), col = 'orange', lwd = 2)

box()

print(paste(c('control', g, n, g / n)))

