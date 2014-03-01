setwd("~/deming_lab/cold_random_HGT/select_genomes")

#oldpar <- par(no.readonly = TRUE)

##### all genomes #####

gc <- read.table('select_genome_gc.txt.gz',
                         header = F)

anomaly <- abs(gc$V6 - gc$V5) - gc$V7
#values above 0 reflect an anomaly larger than the standard deviation

anomaly_2 <- abs(gc$V6 - gc$V5) - 2*gc$V7

anomaly_3 <- abs(gc$V6 - gc$V5) - 3*gc$V7

gc_anomaly <- cbind(gc,
                    anomaly,
                    anomaly_2,
                    anomaly_3)

cold_gc_anomaly <- gc_anomaly[which(gc_anomaly[,2] == 'cold'),]
cold_gc_anomaly <- rbind(cold_gc_anomaly, gc_anomaly[which(gc_anomaly[,2] == 'special'),])

control_gc_anomaly <- gc_anomaly[which(gc_anomaly[,2] == 'control'),]

row_num <- seq(1, length(cold_gc_anomaly$anomaly))
cold_gc_anomaly <- cbind(cold_gc_anomaly, row_num)

row_num <- seq(1, length(control_gc_anomaly$anomaly))
control_gc_anomaly <- cbind(control_gc_anomaly, row_num)

cold_gc_anomaly_large <- cold_gc_anomaly[which(cold_gc_anomaly[,9] > 0),]
control_gc_anomaly_large <- control_gc_anomaly[which(control_gc_anomaly[,9] > 0),]

length(control_gc_anomaly_large[,1]) / length(control_gc_anomaly[,1])
length(cold_gc_anomaly_large[,1]) / length(cold_gc_anomaly[,1])

## read groups file

gc_groups <- read.table('select_genomes.final.groups')

gc_groups[,1] <- sub('.combined.fna', '', gc_groups[,1])

#gc_groups[,1] <- sub('-', '_', gc_groups[,1], perl = TRUE)

#gc_groups[,1] <- sub('\.', '_', gc_groups[,1], perl = TRUE)

gc_groups <- gc_groups[order(gc_groups$V2),]

gc_cold_group <- gc_groups[which(gc_groups$V2 == 'cold'),]

gc_control_group <- gc_groups[which(gc_groups$V2 == 'control'),]

gc_special_group <- gc_groups[which(gc_groups$V2 == 'special'),]

## plot cold

pdf('gc_anomaly_all_orf_cold.pdf', width = 8, height = 6)

par(mar=c(14,5,4,2))

plot(cold_gc_anomaly$anomaly ~ cold_gc_anomaly$row_num,
     type='n',
     ylab = 'abs(strain mean - gene) - strain sd',
     xaxt = 'n',
     xlab = "",
     main = 'cold genomes gc anomalies',
     cex.main = 0.8,
     cex.lab = 0.8)

cold_strain_levels <- levels(factor(cold_gc_anomaly$V1))

row_strain_means <- tapply(cold_gc_anomaly$row_num,
                           factor(cold_gc_anomaly$V1),
                            mean
                            )

axis(1,
     at = row_strain_means,
     labels = cold_strain_levels,
     las = 3,
     cex.axis = 0.6
     )

points(cold_gc_anomaly$row_num[which(cold_gc_anomaly$anomaly < 0)],
       cold_gc_anomaly$anomaly[which(cold_gc_anomaly$anomaly < 0)],
       col='blue',
       pch = 20,
       cex = 0.5
       )

points(cold_gc_anomaly$row_num[intersect(which(cold_gc_anomaly$anomaly_2 < 0), which(cold_gc_anomaly$anomaly > 0))],
       cold_gc_anomaly$anomaly[intersect(which(cold_gc_anomaly$anomaly_2 < 0), which(cold_gc_anomaly$anomaly > 0))],
       col='green',
       pch = 20,
       cex = 0.5
       )

points(cold_gc_anomaly$row_num[intersect(which(cold_gc_anomaly$anomaly_3 < 0), which(cold_gc_anomaly$anomaly_2 > 0))],
       cold_gc_anomaly$anomaly[intersect(which(cold_gc_anomaly$anomaly_3 < 0), which(cold_gc_anomaly$anomaly_2 > 0))],
       col='orange',
       pch = 20,
       cex = 0.5
       )

points(cold_gc_anomaly$row_num[which(cold_gc_anomaly$anomaly_3 > 0)],
       cold_gc_anomaly$anomaly[which(cold_gc_anomaly$anomaly_3 > 0)],
       col='red',
       pch = 20,
       cex = 0.5
       )

#y is the same for all strains
y = c(min(cold_gc_anomaly$anomaly) - 3, max(cold_gc_anomaly$anomaly) + 3)

x <- 0
for(i in cold_strain_levels){
  x <- x + length(which(gc$V1 == i))
  lines(x = c(x,x),
        y,
        col = 'gray')
} 

legend('topright',
       legend = c('0 - 1 sd','1 - 2 sd','2 - 3 sd','> 3 sd'),
       pch = 20,
       col = c('blue','green','orange','red'),
       bg = 'white',
       cex = 0.8,
       y.intersp = 1.0
       )

dev.off()

write.table(cold_gc_anomaly[which(cold_gc_anomaly$anomaly_2 > 0),1],
            'large_gc_anomaly_cold.txt',
            quote = F,
            row.names = F,
            col.names = F
            )

## plot control

pdf('gc_anomaly_all_orf_control.pdf', width = 8, height = 6)

par(mar=c(14,5,4,2))

plot(control_gc_anomaly$anomaly ~ control_gc_anomaly$row_num,
     type='n',
     ylab = 'abs(strain mean - gene) - strain sd',
     xaxt = 'n',
     xlab = "",
     main = 'control genomes gc anomalies',
     cex.main = 0.8,
     cex.lab = 0.8)

control_strain_levels <- levels(factor(control_gc_anomaly$V1))

row_strain_means <- tapply(control_gc_anomaly$row_num,
                           factor(control_gc_anomaly$V1),
                           mean
)

axis(1,
     at = row_strain_means,
     labels = control_strain_levels,
     las = 3,
     cex.axis = 0.6
)

points(control_gc_anomaly$row_num[which(control_gc_anomaly$anomaly < 0)],
       control_gc_anomaly$anomaly[which(control_gc_anomaly$anomaly < 0)],
       col='blue',
       pch = 20,
       cex = 0.5
)

points(control_gc_anomaly$row_num[intersect(which(control_gc_anomaly$anomaly_2 < 0), which(control_gc_anomaly$anomaly > 0))],
       control_gc_anomaly$anomaly[intersect(which(control_gc_anomaly$anomaly_2 < 0), which(control_gc_anomaly$anomaly > 0))],
       col='green',
       pch = 20,
       cex = 0.5
)

points(control_gc_anomaly$row_num[intersect(which(control_gc_anomaly$anomaly_3 < 0), which(control_gc_anomaly$anomaly_2 > 0))],
       control_gc_anomaly$anomaly[intersect(which(control_gc_anomaly$anomaly_3 < 0), which(control_gc_anomaly$anomaly_2 > 0))],
       col='orange',
       pch = 20,
       cex = 0.5
)

points(control_gc_anomaly$row_num[which(control_gc_anomaly$anomaly_3 > 0)],
       control_gc_anomaly$anomaly[which(control_gc_anomaly$anomaly_3 > 0)],
       col='red',
       pch = 20,
       cex = 0.5
)

#y is the same for all strains
y = c(min(control_gc_anomaly$anomaly) - 3, max(control_gc_anomaly$anomaly) + 3)

x <- 0
for(i in control_strain_levels){
  x <- x + length(which(gc$V1 == i))
  lines(x = c(x,x),
        y,
        col = 'gray')
} 

legend('topright',
       legend = c('0 - 1 sd','1 - 2 sd','2 - 3 sd','> 3 sd'),
       pch = 20,
       col = c('blue','green','orange','red'),
       bg = 'white',
       cex = 0.8,
       y.intersp = 1.0
)

dev.off()

write.table(control_gc_anomaly[which(control_gc_anomaly$anomaly_2 > 0),1],
            'large_gc_anomaly_control.txt',
            quote = F,
            row.names = F,
            col.names = F
)


##### comparison #####

cold_anomaly_strain <- matrix(nrow = length(cold_strain_levels), ncol = 2)
n <- 0
cold_length <- NULL

for (i in cold_strain_levels){
  n = n + 1
  temp_anomaly <- cold_gc_anomaly_large[which(cold_gc_anomaly_large$V1 == i),1]
  temp_all <- cold_gc_anomaly[which(cold_gc_anomaly$V1 == i),1]
  #temp_row <- c(i, length(temp_anomaly)/length(temp_all))
  temp_row <- c(i, length(temp_anomaly))
  cold_anomaly_strain[n,] <- temp_row
  print(c(i,length(temp_anomaly)/length(temp_all)), quote=F)
  cold_length <- append(cold_length, length(temp_all))
}

control_anomaly_strain <- matrix(nrow = length(control_strain_levels), ncol = 2)
n <- 0
control_length <- NULL

for (i in control_strain_levels){
  n = n + 1
  temp_anomaly <- control_gc_anomaly_large[which(control_gc_anomaly_large$V1 == i),1]
  temp_all <- control_gc_anomaly[which(control_gc_anomaly$V1 == i),1]
  #temp_row <- c(i, length(temp_anomaly)/length(temp_all))
  temp_row <- c(i, length(temp_anomaly))
  control_anomaly_strain[n,] <- temp_row
  print(c(i,length(temp_anomaly)/length(temp_all)), quote=F)
  control_length <- append(control_length, length(temp_all))
}

boxplot(as.numeric(cold_anomaly_strain[,2]), as.numeric(control_anomaly_strain[,2]),
        notch = T)

t.test(as.numeric(cold_anomaly_strain[,2]), as.numeric(control_anomaly_strain[,2]))
hist(as.numeric(control_anomaly_strain[,2]), breaks = 50)
hist(as.numeric(cold_anomaly_strain[,2]), breaks = 50)

## distributions not normal

library(nortest)

ad.test(as.numeric(cold_anomaly_strain[,2]))
ad.test(as.numeric(control_anomaly_strain[,2]))

qqnorm(as.numeric(control_anomaly_strain[,2]))

wilcox.test(as.numeric(cold_anomaly_strain[,2]), as.numeric(control_anomaly_strain[,2]))

## one to one comparison

hgt_compare <- NULL

for(x in gc_cold_group[,3]){
  control_hgt_temp <- gc_control_group[which(gc_control_group[,3] == x),1]
  cold_hgt_temp <- gc_cold_group[which(gc_cold_group[,3] == x),1]
  hgt_cold <- as.numeric(cold_anomaly_strain[which(cold_anomaly_strain[,1] == cold_hgt_temp),2])
  hgt_control <- as.numeric(control_anomaly_strain[which(control_anomaly_strain[,1] == control_hgt_temp),2])
  delta_hgt_temp <- hgt_cold - hgt_control
  print(c(cold_hgt_temp, control_hgt_temp, delta_hgt_temp))
  hgt_compare <- rbind(hgt_compare, c(cold_hgt_temp, hgt_cold, control_hgt_temp, hgt_control, delta_hgt_temp))
}

write.table(hgt_compare, 'hgt_compare.txt', sep = '\t', quote = F)

## compare pfams between populations

cold_pfams <- data.frame(table(cold_gc_anomaly_large$V4))
control_pfams <- data.frame(table(control_gc_anomaly_large$V4))

diff_pfams <- data.frame(as.character(cold_pfams[,1]), (cold_pfams[,2] - control_pfams[,2]))
diff_pfams <- diff_pfams[-which(diff_pfams[,2] == 0),]
diff_pfams <- diff_pfams[order(abs(diff_pfams[,2]), decreasing = T),]
diff_pfams_cold_elevated <- diff_pfams[-which(diff_pfams[,2] < 0),]

write.table(diff_pfams, 'diff_pfams.txt', sep = '\t', quote = F)

cold_elevated_pfams <- cold_gc_anomaly_large[which(as.character(cold_gc_anomaly_large$V4) %in% as.character(diff_pfams_cold_elevated[1:50,1])),]

## which HGT pfam is most broadly distributed across cold genomes?

cold_elevated_pfams_distrib <- NULL

for(p in unique(as.character(cold_elevated_pfams$V4))){
  temp <- c(p,length(unique(cold_elevated_pfams[which(cold_elevated_pfams$V4 == p),1])))
  cold_elevated_pfams_distrib <- rbind(cold_elevated_pfams_distrib, temp)
}

cold_elevated_pfams_distrib <- cold_elevated_pfams_distrib[order(cold_elevated_pfams_distrib[,1]),]



