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


fit.sigmoid.gc1 <- function(x_in) {
#zunzun.com
#finds gc3 given gc1

	a <- 1.4180794193642083E+00
	b <- 1.2403719287318339E+02
	c <- 7.5494898937434289E+00

	y = a / (1.0 + b*exp(-1.0 * c * x_in))
}

fit.sigmoid.gc2 <- function(x_in) {
#zunzun.com
#finds gc3 given gc2

	a = 1.0734349192542731E+00;
	b = 2.7879639168098345E+02;
	c = 1.3906208449024092E+01;

	y = a / (1.0 + b*exp(-1.0 * c * x_in));
}

reps <- 1000 #bootstrap reps
mya <- 1000 #millions of years to ameliorate
tstv <- 2
Srate <- c(0.123,0.045,0.668)/100/2 #weighted substitution rates

file.genome <- 'genome.ffn'

gcs.host <- as.numeric(unlist(strsplit(system("python gc_by_pos_3.py 0 1 genome.ffn",intern=TRUE),"\t")))

list.xenes <- list.files(pattern='.fasta')
tminmax <- matrix(ncol = 4, nrow = reps * length(list.xenes))
oldgc <- matrix(ncol = 4, nrow = reps * length(list.xenes))

xs <- 1:100/100
gc1.fit <- fit.sigmoid.gc1(xs)
gc2.fit <- fit.sigmoid.gc2(xs)

n = 0

for (myfile in list.xenes) {
  file.xene <- myfile
  print(file.xene)

  gcs.xene.get <- as.numeric(unlist(strsplit(system(paste("python gc_by_pos_3.py 0 ",reps,paste('"',file.xene,'"',sep=''),sep=" "),intern=TRUE),"\t")))
  gcs.xene <- matrix(gcs.xene.get,ncol=3,byrow=TRUE)

 for (rep in 1:reps) {
  n = n + 1
  print(rep);
  rev.xene.save <- matrix(ncol = 3, nrow = reps)

  delGC.rev <-c(0,0,0)
  rev.xene <- gcs.xene[rep,]

  dist.gc1 <- rep.int(0,mya)
  dist.gc2 <- rep.int(0,mya)
  dist.both <- rep.int(0,mya)

  for (i in 1:mya) {
    delGC.rev <- Srate * (tstv+0.5)/(tstv+1) * (gcs.host - rev.xene)
    rev.xene <- rev.xene - delGC.rev;
    rev.xene.save[rep,] <- rev.xene
  }

  dist.gc1 <- abs(fit.sigmoid.gc1(rev.xene.save[,1]) - rev.xene.save[,3])
  dist.gc2 <- abs(fit.sigmoid.gc2(rev.xene.save[,2]) - rev.xene.save[,3])

  tminmax[n,] <- c(file.xene,rep,which.min(dist.gc2),which.min(dist.gc1))
  oldgc[n,] <- c(file.xene,rep,rev.xene.save[which.min(dist.gc2)],rev.xene.save[which.min(dist.gc1)])

 }
}

tempfiles <- tempfile("gcamel-",tmpdir=".")
write.csv(tminmax,file=paste(tempfiles[1],"-tminmax.csv",sep=""))
write.csv(oldgc,file=paste(tempfiles[1],"-oldgc.csv",sep=""))




