setwd('~/deming_lab/cold_random_HGT/select_genomes')

library(vegan)

##### combined #####

## read groups file

groups <- read.table('select_genomes.final.groups')

groups[,1] <- sub('.combined.fna', '', groups[,1])

groups[,1] <- sub('-', '_', groups[,1], perl = TRUE)

groups[,1] <- sub('\.', '_', groups[,1], perl = TRUE)

groups <- groups[order(groups$V2),]

cold_group <- groups[which(groups$V2 == 'cold'),]

control_group <- groups[which(groups$V2 == 'control'),]

special_group <- groups[which(groups$V2 == 'special'),]

## calculate compositional vector distance - can skip if you don't want to recalculate

name <- '5mer_normalized_phylogeny_vector_output.txt.gz'

d <- as.matrix(read.table(name,
                          header = T))

d_trans <- d[which(rowSums(d) != 0),]

#d_trans <- d_trans + abs(min(d_trans)) + 1

d_dist_raw <- as.matrix(vegdist(t(d_trans), "euclidean"))

d_dist <- d_dist_raw

#d_dist <- d_dist / max(d_dist, na.rm = T) # no longer necessary given downstream transformation and scaling

row.names(d_dist) <- sub('-', '_', row.names(d_dist), perl = TRUE)

row.names(d_dist) <- sub('.combined', '', row.names(d_dist), perl = TRUE)

row.names(d_dist) <- sub('\\.', '_', row.names(d_dist), perl = TRUE)

colnames(d_dist) <- row.names(d_dist)

row.names(d_dist) <- colnames(d_dist)

## set order

row_order <- c()

for(n in groups$V1){
  print(n)
  print(n %in% row.names(d_dist))
  row_order <- append(row_order, which(row.names(d_dist) == n))
}

d_dist <- d_dist[row_order, row_order]


write.table(data.frame(d_dist), file = 'euclidean_dist.txt', quote = F)

## if you skipped the recalculation load earlier calc

d_dist <- as.matrix(read.table('euclidean_dist.txt'))

##### 16S #####

## read 16S distance matrix and get in correct order

ss <- as.matrix(read.table('combined_16S.good.filter.square.dist', row.names = 1, skip = 1))

row.names(ss) <- sub('-', '_', row.names(ss), perl = TRUE)

row.names(ss) <- sub('.combined', '', row.names(ss), perl = TRUE)

row.names(ss) <- sub('\\.', '_', row.names(ss), perl = TRUE)

colnames(ss) <- row.names(ss)

row_order <- c()

for(n in groups$V1){
  row_order <- append(row_order, which(row.names(ss) == n))
}

ss <- ss[row_order, row_order]

##### normalize mean to 0 and var to 1 for both matrices #####

d_dist[d_dist == 0] <- NA
ss[ss == 0] <- NA

ss_vector <- NULL

for(r in seq(1, length(ss[,1]))){
  for(c in seq(1, length(ss[1,]))){
    ss_vector <- append(ss_vector, ss[c,r])
  }
}

ss_mean <- mean(ss_vector, na.rm = T)
ss_sd <- sd(ss_vector, na.rm = T)
ss <- (ss - ss_mean) / ss_sd
ss <- ss + abs(min(ss, na.rm = T))
ss <- ss / max(ss, na.rm = T)

## d_dist

d_dist_vector <- NULL

for(r in seq(1, length(d_dist[,1]))){
  for(c in seq(1, length(d_dist[1,]))){
    d_dist_vector <- append(d_dist_vector, d_dist[c,r])
  }
}

d_dist_mean <- mean(d_dist_vector, na.rm = T)
d_dist_sd <- sd(d_dist_vector, na.rm = T)
d_dist <- (d_dist - d_dist_mean) / d_dist_sd
d_dist <- d_dist + abs(min(d_dist, na.rm = T))
d_dist <- d_dist / max(d_dist, na.rm = T)

plot(ss ~ d_dist)

##### generate plots #####

color <- colorRampPalette(c('white', 'blue', 'green', 'yellow', 'orange', 'red'))(100)

## divg

divg <- as.matrix(abs(d_dist - ss))

pdf('divg_heatmap.pdf', width = 10, height = 10)

ratio_heat <- heatmap(divg,
                      margins = c(25,25),
                      revC = T,
                      Rowv = NA,
                      Colv = NA,
                      col = color,
                      symm = T,
                      cexRow = 0.9,
                      cexCol = 0.9,
                      main = 'normalized compositional vector distance',
                      ColSideColors = c(rep('blue', 20), rep('black', 20), 'red'),
                      RowSideColors = c(rep('blue', 20), rep('black', 20), 'red')
                      )

dev.off()

## 16S

pdf('16S_heatmap.pdf', width = 10, height = 10)

ratio_heat <- heatmap(ss,
                      margins = c(25,25),
                      revC = T,
                      Rowv = NA,
                      Colv = NA,
                      col = color,
                      symm = T,
                      cexRow = 0.9,
                      cexCol = 0.9,
                      main = 'normalized 16S rRNA distance',
                      ColSideColors = c(rep('blue', 20), rep('black', 20), 'red'),
                      RowSideColors = c(rep('blue', 20), rep('black', 20), 'red')
)

dev.off()

## cv

pdf('cv_heatmap.pdf', width = 10, height = 10)

ratio_heat <- heatmap(d_dist,
                      margins = c(25,25),
                      revC = T,
                      Rowv = NA,
                      Colv = NA,
                      col = color,
                      symm = T,
                      cexRow = 0.9,
                      cexCol = 0.9,
                      main = 'abs(genome divergence - 16S gene divergence)',
                      ColSideColors = c(rep('blue', 20), rep('black', 20), 'red'),
                      RowSideColors = c(rep('blue', 20), rep('black', 20), 'red')
)

dev.off()

#### population comparison ####

ss_cold <- ss[which(row.names(ss) %in% cold_group[,1]), which(row.names(ss) %in% cold_group[,1])]

ss_control <- ss[which(row.names(ss) %in% control_group[,1]), which(row.names(ss) %in% control_group[,1])]

d_dist_cold <- d_dist[which(row.names(d_dist) %in% cold_group[,1]), which(row.names(d_dist) %in% cold_group[,1])]

d_dist_control <- d_dist[which(row.names(d_dist) %in% control_group[,1]), which(row.names(d_dist) %in% control_group[,1])]

divg_cold <- divg[which(row.names(divg) %in% cold_group[,1]), which(row.names(divg) %in% cold_group[,1])]

divg_control <- divg[which(row.names(divg) %in% control_group[,1]), which(row.names(divg) %in% control_group[,1])]

## mean, sd

cold_row_mean <- apply(divg_cold, 1, mean, na.rm = TRUE)
cold_row_sd <- apply(divg_cold, 1, sd, na.rm = TRUE)
cold_row_sum <- apply(divg_cold, 1, sum, na.rm = TRUE)

control_row_mean <- apply(divg_control, 1, mean, na.rm = TRUE)
control_row_sd <- apply(divg_control, 1, sd, na.rm = TRUE)
control_row_sum <- apply(divg_control, 1, sum, na.rm = TRUE)

library(nortest)

ad.test(c(cold_row_mean, control_row_mean)) # looks normal, I think - use t.test

t.test(cold_row_mean, control_row_mean) # not signif
wilcox.test(cold_row_mean, control_row_mean) # signif

mean(control_row_mean)
mean(cold_row_mean)

boxplot(cold_row_mean, control_row_mean, notch = T)

## box plots of ss and d_dist and divG

divg_cold_vector <- NULL
divg_control_vector <- NULL

for(r in seq(1, length(divg_cold[,1]))){
  for(c in seq(1, length(divg_cold[1,]))){
    divg_cold_vector <- append(divg_cold_vector, divg_cold[c,r])
  }
}

for(r in seq(1, length(divg_control[,1]))){
  for(c in seq(1, length(divg_control[1,]))){
    divg_control_vector <- append(divg_control_vector, divg_control[c,r])
  }
}

ss_cold_norm_vector <- NULL
ss_control_norm_vector <- NULL

for(r in seq(1, length(ss_cold[,1]))){
  for(c in seq(1, length(ss_cold[1,]))){
    ss_cold_norm_vector <- append(ss_cold_norm_vector, ss_cold[c,r])
  }
}

for(r in seq(1, length(ss_control[,1]))){
  for(c in seq(1, length(ss_control[1,]))){
    ss_control_norm_vector <- append(ss_control_norm_vector, ss_control[c,r])
  }
}

d_dist_cold_norm_vector <- NULL
d_dist_control_norm_vector <- NULL

for(r in seq(1, length(d_dist_cold[,1]))){
  for(c in seq(1, length(d_dist_cold[1,]))){
    d_dist_cold_norm_vector <- append(d_dist_cold_norm_vector, d_dist_cold[c,r])
  }
}

for(r in seq(1, length(d_dist_control[,1]))){
  for(c in seq(1, length(d_dist_control[1,]))){
    d_dist_control_norm_vector <- append(d_dist_control_norm_vector, d_dist_control[c,r])
  }
}

boxplot(cold_row_mean,
        control_row_mean,
        divg_cold_vector,
        divg_control_vector,
        d_dist_cold_norm_vector,
        d_dist_control_norm_vector,
        notch = T,
        ss_cold_norm_vector,
        ss_control_norm_vector,
        col = c('blue', 'red', 'blue', 'red', 'blue', 'red', 'blue', 'red'),
        names = c('mdivG', 'mdivG', 'divG', 'divG', 'CV', 'CV', '16S', '16S'),
        ylab = 'Normalized distance')

t.test(d_dist_cold_norm_vector, d_dist_control_norm_vector) # very signif
t.test(ss_cold_norm_vector, ss_control_norm_vector) # not signif
t.test(divg_cold_vector, divg_control_vector) # signif

## 1 to 1 comparisons ##

cold_temp <- NULL
control_temp <- NULL
divg_temp <- NULL
cold_divg <- NULL
control_divg <- NULL

for(x in cold_group[,3]){
  control_temp <- append(control_temp, control_group[which(control_group[,3] == x),1])
  cold_temp <- append(cold_temp, cold_group[which(cold_group[,3] == x),1])
  divg_temp <- append(divg_temp, sum(divg_cold[which(row.names(divg_cold) == tail(cold_temp, n = 1)),], na.rm = T) - sum(divg_control[which(row.names(divg_control) == tail(control_temp, n = 1)),], na.rm = T))
  cold_divg <- append(cold_divg, sum(divg_cold[which(row.names(divg_cold) == tail(cold_temp, n = 1)),], na.rm = T))
  control_divg <- append(control_divg, sum(divg_control[which(row.names(divg_control) == tail(control_temp, n = 1)),], na.rm = T))
}

col_sum_compare <- data.frame(cold_temp, cold_divg, control_temp, control_divg, divg_temp)

write.table(col_sum_compare, 'divg_compare.txt', quote = F, sep = '\t')

##### shared genera  - below here not in current manuscript #####

select_cold_groups <- c('Aeromonas_salmonicida_A449_uid16723',
                 'Flavobacterium_psychrophilum_JIP02_86_uid19979',
                'Glaciecola_psychrophila_170_uid174842',
                 'Pseudoalteromonas_haloplanktis_TAC125_uid15713',
                 'Shewanella_halifaxensis_HAW_EB4_uid20241',
                 'Shewanella_sediminis_HAW-EB3_uid18789',
                 'Shewanella_violacea_DSS12_uid34739',
                 'Terriglobus_saanensis_SP1PR4_uid48971')

select_control_groups <- c('Aeromonas_veronii_B565_uid63671',
                           'Flavobacterium_branchiophilum_FL_15_uid67123',
                           'Glaciecola_agarilytica_4H_3_7_YE_5_uid62887',
                           'Pseudoalteromonas_atlantica_T6c_uid13454',
                           'Shewanella_MR-7_uid13903',
                           'Shewanella_denitrificans_OS217_uid13390',
                           'Shewanella_oneidensis_uid335',
                           'Shewanella_putrefaciens_200_uid13392',
                           'Terriglobus_saanensis_SP1PR4_uid48971')

select_groups <- c(select_cold_groups, select_control_groups)

select_divg <- divg[which(row.names(divg) %in% select_groups),which(colnames(divg) %in% select_groups)]

## individual ##

select_divg_cold <- select_divg[which(row.names(select_divg) %in% cold_group[,1]), which(row.names(select_divg) %in% cold_group[,1])]

select_divg_control <- select_divg[which(row.names(select_divg) %in% control_group[,1]), which(row.names(select_divg) %in% control_group[,1])]

## mean, sd

select_cold_row_mean <- apply(select_divg_cold, 1, mean, na.rm = TRUE)
select_cold_row_sd <- apply(select_divg_cold, 1, sd, na.rm = TRUE)
select_cold_row_sum <- apply(select_divg_cold, 1, sum, na.rm = TRUE)

select_control_row_mean <- apply(select_divg_control, 1, mean, na.rm = TRUE)
select_control_row_sd <- apply(select_divg_control, 1, sd, na.rm = TRUE)
select_control_row_sum <- apply(select_divg_control, 1, sum, na.rm = TRUE)

t.test(select_cold_row_mean, select_control_row_mean)

mean(select_control_row_mean)
mean(select_cold_row_mean)

boxplot(select_cold_row_mean, select_control_row_mean)

##### sea ice #####

seaice_groups_cold <- c('Glaciecola_psychrophila_170_uid174842',
                        'Octadecabacter_arcticus_238_uid19331',
                        'Psychroflexus_torquis_ATCC_700755_uid13542',
                        'Psychromonas_ingrahamii_37_uid16187')

seaice_groups_control <- c('Glaciecola_agarilytica_4H_3_7_YE_5_uid62887',
                           'Ketogulonicigenium_vulgare_Y25_uid51787',
                           'Flavobacteriales_bacterium_HTCC2170_uid13595',
                           'Marinobacter_aquaeolei_VT8_uid13239')

seaice_groups <- c(seaice_groups_cold, seaice_groups_control)

seaice_divg <- divg[which(row.names(divg) %in% seaice_groups),which(colnames(divg) %in% seaice_groups)]

seaice_divg_cold <- seaice_divg[which(row.names(seaice_divg) %in% seaice_groups_cold), which(row.names(seaice_divg) %in% seaice_groups_cold)]

seaice_divg_control <- seaice_divg[which(row.names(seaice_divg) %in% seaice_groups_control), which(row.names(seaice_divg) %in% seaice_groups_control)]

seaice_cold_row_mean <- apply(seaice_divg_cold, 1, mean, na.rm = TRUE)

seaice_control_row_mean <- apply(seaice_divg_control, 1, mean, na.rm = TRUE)

t.test(seaice_cold_row_mean, seaice_control_row_mean)



