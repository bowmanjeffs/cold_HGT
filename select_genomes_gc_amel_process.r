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

aquisition_dates <- read.table('aquisition_means.txt', sep = '\t')
aquisition_dates_cold <- aquisition_dates[which(aquisition_dates[,4] == 'cold'),]
aquisition_dates_control <- aquisition_dates[which(aquisition_dates[,4] == 'control'),]

numeric_means <- as.numeric(aquisition_dates[,5])
numeric_errors <- as.numeric(aquisition_dates[,8])

#glacial <- c(1:55,105:183,253:333,421:458) # "climate modes of the phanerozoic"
#glacial <- c(1:34, 65:72, 89:98, 112:125, 144:150, 136:141, 159:168, 183:184, 267:326, 349:361, 443:446) # Royer, 2005
#glacial <- c(112:125, 144:150, 136:141, 159:168, 183:184, 267:326, 349:361, 443:446, 542, 580, 635:650, 660:720, 730:775) # Royer, 2005 and Koll 2011 excepting < 100 Ma
#glacial_start <- c(112, 144, 136, 159, 183, 267, 349, 443, 542, 580, 635, 660, 730)
#glacial_end <- c(125, 150, 141, 168, 184, 326, 361, 446, 542, 580, 650, 720, 775)
glacial <- c(1:50,105:220,270:370,390:475, 635:650, 660:720, 730:775) # Veizer et al. 2000 + snowballs
#glacial <- c(105:220,270:370,390:475,625:775) # Veizer et al. 2000, excluding most recent 100 yrs
glacial_hgt <- numeric_means[which(ceiling(numeric_means) %in% glacial)]

m <- 775 # extent of glacial record, alternatively x intercept of polynomial fit
f_glacial <- length(glacial) / m # fraction of time that is glacial
f_glacial_hgt <- length(glacial_hgt) / length(numeric_means[which(numeric_means < m)])

## bin data and export for model fit in pydavis
master_hist <- hist(numeric_means[which(numeric_means < m)], breaks = seq(0,m + 5,10))

write.table(data.frame(master_hist$counts, master_hist$mids), 'master_dates_histogram.txt', quote = F, sep = '\t', col.names = F, row.names = F)

##### correcting for more glaciation in recent history - when HGT is more detectable #####

## want deviation from a constant rate of HGT, as identified by best fit line pydavis
## representing as fraction of anticipated HGT events - 1

master_hist_norm <- c()

a0 <- 207.220087607263
a1 <- -0.55722808836733
a2 <-  0.000374254659397898
#a3 <- -5.85310237922088e-07
#a4 <- 2.02779268632868e-10

for(i in seq(1:max(which(master_hist$mids < m)))){
  x <- master_hist$mids[i]
  temp <- a0+a1*x+a2*x^2 ## curve fit with pydavis
  if(temp < 1){temp <- 1}
  if(temp <= master_hist$count[i]){
    y_norm <- (master_hist$count[i] / temp) - 1
  }else{y_norm <- 0}
    print(c(y_norm, temp, master_hist$count[i]))
    master_hist_norm <- append(master_hist_norm, y_norm)
}

## plot HGT occurence

plot(master_hist$counts ~ master_hist$mids,
     type = 'n',
     ylab = 'HGT events',
     xlab = 'Ma',
     ylim = c(0,300)
)

rect(c(1,105,270,390,635,660),
     c(0),
     c(50,220,370,475,650,720),
     c(300),
     border = F,
     col = 'grey')

points(master_hist$counts ~ master_hist$mids, pch = 20)

curve(a0+a1*x+a2*x^2, add = T, col = 'red', lwd = 2)


## plot normalized HGT occurrence

plot(master_hist_norm ~ master_hist$mids[which(master_hist$mids < m)],
     type = 'n',
     ylim = c(0,1),
     ylab = 'Normalized HGT events',
     xlab = 'Ma',
     xlim = c(0,m)
     )

rect(c(1,105,270,390,635,660),
     c(0),
     c(50,220,370,475,650,720),
     c(1),
     border = F,
     col = 'grey')

barplot(master_hist_norm,
        width = 10,
        space = 0,
        add = T,
        col = 'black'
        )

## plotting complete

f_glacial <- length(glacial) / m 
glacial_hgt_norm <- sum(master_hist_norm[(ceiling(master_hist$mids) %in% glacial)])
f_glacial_hgt_norm <- glacial_hgt_norm / sum(master_hist_norm, na.rm = T)

## MC simulation
n <- 1000000 ## number of simulations to run
l <- 0
g <- 0
mcs <- vector("numeric",length = n)

master_hist_mids <- master_hist$mids[which(master_hist$mids < m)]
glacial_mids <- which(master_hist_mids %in% glacial)

for(i in seq(1,n)){
  print(i)
  rn <- sample(master_hist_norm, length(master_hist_norm))
  mc_glacial_hgt_norm <- sum(rn[glacial_mids])
  f_mc_glacial_hgt_norm <- mc_glacial_hgt_norm / sum(master_hist_norm)
  if(f_mc_glacial_hgt_norm < f_glacial_hgt_norm){ l <- l + 1 } else { g <- g + 1 }
  mcs[i] <- f_mc_glacial_hgt_norm
}

hist_mcs <- hist(mcs, breaks = 100,
                 xlab = 'Fraction of events occurring during cold period',
                 main = NULL,
                 col = 'black',
                 cex.axis = 1.4,
                 cex.lab = 1.4)

lines(c(f_glacial_hgt_norm, f_glacial_hgt_norm), c(0, 100000), col = 'orange', lwd = 2)

box()

mc_error <- qnorm(0.9875) * sd(mcs)
mc_left <- mean(mcs) - mc_error
mc_right <- mean(mcs) + mc_error