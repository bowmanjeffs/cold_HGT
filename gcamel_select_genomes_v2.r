#setwd('~/deming_lab/cold_random_HGT/select_genomes/gc_amel_develop')

args <- commandArgs() # working directory
wd <- args[3]

setwd(wd)

rm(list=ls())
# parameters needed:
#  1) GC content of genome at each position (.ffn file)
#  2) GC content of xene at each position (.fasta file)
#  3) transition/transversion ratio of genome (assume 2:1)
#  4) substitution rate of genome (assume 0.0245% and 0.455% per Myr for nonsynonymous and synonymous)
#  5) linear (?) relationships between {GC1,GC2} and {GC3} from Muta & Osawa
# 
# wanted:
#  1) GC content of xene at introduction
#  2) time since introduction of xene

# method to get 1) and 2)
#  a) calculate delGC for 1Myr time interval
#     delGC = S * (tstv+0.5)/(tstv+1) * (gc_genome - gc_xene)
#  b) add delGC to gc_xene
#  c) calculate correlation between GC1,GC2~GC3 and Muta&Osawa
#  d) repeat a-c for 1000Myr
#  e) find minimum correlation for all c) -- this point is initial GC and time since transfer

#### define functions ####

#zunzun.com
#finds gc3 given gc1
fit.sigmoid.gc1 <- function(x_in) {

	a <- 1.4180794193642083E+00
	b <- 1.2403719287318339E+02
	c <- 7.5494898937434289E+00

	y = a / (1.0 + b*exp(-1.0 * c * x_in))
}

#zunzun.com
#finds gc3 given gc2
fit.sigmoid.gc2 <- function(x_in) {

	a = 1.0734349192542731E+00;
	b = 2.7879639168098345E+02;
	c = 1.3906208449024092E+01;

	y = a / (1.0 + b*exp(-1.0 * c * x_in));
}

#### run parameters ####

reps <- 1000 # bootstrap reps
mya_max <- 1000 # thousands of years to ameliorate
mya_min <- 1
tstv <- 2 # ratio of transition to transversion
Srate <- c(0.123,0.045,0.668)/100/2 # weighted substitution rates, see lawrence and ochman, 1997, table 1 - 3rd position mutates fastest because mutations are synonymous

file.genome <- 'genome.ffn'

gcs.host <- as.numeric(unlist(strsplit(system("python gc_by_pos_3.py 0 1 genome.ffn",intern=TRUE),"\t")))

list.xenes <- list.files(pattern='.fasta')
tminmax <- matrix(ncol = 4, nrow = reps * length(list.xenes))
oldgc <- matrix(ncol = 4, nrow = reps * length(list.xenes))
aquisition.dates <- matrix(ncol = 4, nrow = reps * length(list.xenes))

xs <- 1:100/100 # gc1 for fit.sigmoid.gc1, gc2 for fit.sigmoid.gc2
gc1.fit <- fit.sigmoid.gc1(xs) # gc3 estimate based on gc1
gc2.fit <- fit.sigmoid.gc2(xs) # gc3 estimate based on gc2

#plot(gc1.fit ~ xs)
#points(gc2.fit ~ xs, pch = 19)

#### execute main loop ####

myfile <- list.xenes[1] # for testing only

n <- 0
f <- 0

for (myfile in list.xenes) {
  f <- f + 1
  file.xene <- myfile
  print(paste(file.xene, f, 'out of', length(list.xenes)))

  ## jacknife codons in xene [# reps] times
  gcs.xene.get <- as.numeric(unlist(strsplit(system(paste("python gc_by_pos_3.py 0 ",reps,paste('"',file.xene,'"',sep=''),sep=" "),intern=TRUE),"\t")))
  gcs.xene <- matrix(gcs.xene.get,ncol=3,byrow=TRUE)
  
#   plot(c(seq(-mya_max, mya_max, length = length(c(-5:5)))),
#        c(-5:5),
#        type = 'n',
#        xlab = 'mya',
#        ylab = 'GC3 amel ~ GC1,GC2 amel per million years')
#   
#   legend('bottomleft',
#          legend = c('GC3 amelioration prediction', 'GC3 amelioration predicted from GC1', 'GC3 amelioration predicted from GC2'),
#          col = c('black', 'blue', 'red'),
#          lty = 1)

 for (rep in 1:reps) {
  n = n + 1
  rev.xene.save <- matrix(ncol = 3, nrow = length(mya_min:mya_max)) ## nrows = number of time steps
  fwd.xene.save <- matrix(ncol = 3, nrow = length(mya_min:mya_max)) ## nrows = number of time steps

  delGC.rev <- c(0,0,0)
  delGC.fwd <-c(0,0,0)
  
  fwd.xene <- gcs.xene[rep,]
  rev.xene <- gcs.xene[rep,]

  dist.gc1 <- rep.int(0,mya_max + abs(mya_min))
  dist.gc2 <- rep.int(0,mya_max + abs(mya_min))

  ## predict gc content of each codon position at each time point (actual reverse amerlioration step)
  for (i in mya_min:mya_max) {
    delGC.rev <- Srate * (tstv+0.5)/(tstv+1) * (gcs.host - rev.xene) # going further back, difference between host and gene decreases
    rev.xene <- rev.xene - delGC.rev
    rev.xene.save[i,] <- rev.xene
    
    delGC.fwd <- Srate * (tstv+0.5)/(tstv+1) * (gcs.host - fwd.xene)
    fwd.xene <- fwd.xene + delGC.fwd;
    fwd.xene.save[i,] <- fwd.xene
    
  }

#   points(rev.xene.save[,3] ~ c(mya_min:mya_max),
#          type = 'l',
#          col = 'black')
#   
#   points(fwd.xene.save[,3] ~ c(-mya_min:-mya_max),
#          type = 'l',
#          col = 'black')
#   
#   points(fit.sigmoid.gc1(rev.xene.save[,1]) ~ c(mya_min:mya_max),
#          type = 'l',
#          col = 'red')
#   
#   points(fit.sigmoid.gc1(fwd.xene.save[,1]) ~ c(-mya_min:-mya_max),
#          type = 'l',
#          col = 'red')
#          
#   points(fit.sigmoid.gc2(rev.xene.save[,2]) ~ c(mya_min:mya_max),
#          type = 'l',
#          col = 'blue')
#   
#   points(fit.sigmoid.gc2(fwd.xene.save[,2]) ~ c(-mya_min:-mya_max),
#          type = 'l',
#          col = 'blue')
                
  ## calculate distance between predicted gc3 (based on input of gc1, gc2) and actual gc3
  dist.gc1.rev <- abs(fit.sigmoid.gc1(rev.xene.save[,1]) - rev.xene.save[,3])
  min.dist.gc1.rev <- dist.gc1.rev[which.min(dist.gc1.rev)]
  
  dist.gc1.fwd <- abs(fit.sigmoid.gc1(fwd.xene.save[,1]) - fwd.xene.save[,3])
  min.dist.gc1.fwd <- dist.gc1.fwd[which.min(dist.gc1.fwd)] 
  
  if(min.dist.gc1.fwd < min.dist.gc1.rev){
    min.dist.gc1 <- -1 * which.min(dist.gc1.fwd)
  }else{
    min.dist.gc1 <- which.min(dist.gc1.rev)
  }

  dist.gc2.rev <- abs(fit.sigmoid.gc2(rev.xene.save[,2]) - rev.xene.save[,3])
  min.dist.gc2.rev <- dist.gc2.rev[which.min(dist.gc2.rev)]
  
  dist.gc2.fwd <- abs(fit.sigmoid.gc2(fwd.xene.save[,2]) - fwd.xene.save[,3])
  min.dist.gc2.fwd <- dist.gc2.fwd[which.min(dist.gc2.fwd)] 
  
  if(min.dist.gc2.fwd < min.dist.gc2.rev){
    min.dist.gc2 <- -1 * which.min(dist.gc2.fwd)
  }else{
    min.dist.gc2 <- which.min(dist.gc2.rev)
  }
  
  ## timepoint with minimum distance between calculated and actual gc3 is taken as intersect
  tminmax[n,] <- c(file.xene,rep,min.dist.gc2,min.dist.gc1)
  oldgc[n,] <- c(file.xene,rep,rev.xene.save[which.min(dist.gc2)],rev.xene.save[which.min(dist.gc1)])

 }
}

tempfiles <- tempfile("gcamel-",tmpdir=".")
write.csv(tminmax,file=paste(tempfiles[1],"-tminmax.csv",sep=""))
write.csv(oldgc,file=paste(tempfiles[1],"-oldgc.csv",sep=""))
