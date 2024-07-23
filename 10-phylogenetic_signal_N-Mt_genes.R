library(ape)
library(data.table)

prop.weight <- function(raw){raw/apply(raw, 1, sum)}

setwd("Mito-nuclear_genes/")

chr.map <- read.table("chromosomes_map.tsv", sep="\t", header=T)
chr.map <- chr.map[1:46,c(4,7,10)]
chr.map$Chromosome.name <- paste("SUPER_", chr.map$Chromosome.name, sep="")

genes.gff <- read.gff("mito_nuclear_genes_2.gff")
genes.gff <- genes.gff[genes.gff$type == "gene", ]
genes.gff <- genes.gff[grep("NC_", genes.gff$seqid), ]

chr <- c()
for(i in 1:nrow(genes.gff)){
  chr <- c(chr, chr.map$Chromosome.name[which(chr.map$RefSeq.seq.accession == genes.gff[i,1])])
}

genes.gff <- cbind(genes.gff, chr)
genes.gff.aut <- genes.gff[genes.gff$chr != "SUPER_Z", ]

results <- matrix(ncol=27, nrow = 0)

for(i in unique(genes.gff.aut$chr)){
  curr.gff <- as.data.table(genes.gff.aut[genes.gff.aut$chr == i, c(4,5,10), ])
  curr.w <- read.table(paste("../TWISST_weights/",i,"_all_W", sep=""), header=T)
  curr.w <- prop.weight(curr.w)
  curr.stats <- read.table(paste("../Neighbour_Joining_trees/",i,"_all_500S_maxmiss60_windows_stats.tsv", sep=""), header=T)
  curr.stats <- curr.stats[curr.stats$TREE == "YES", ]
  curr.stats <- curr.stats[curr.stats$NTIPS == 45, ]
  curr.stats <- curr.stats[which(curr.stats$PROP.MISS < 0.4 & curr.stats$CHR.SIZE < 50000), ]
  curr.stats.w <- as.data.table(cbind(curr.stats, curr.w))
  setkey(curr.gff, start, end)
  setkey(curr.stats.w, CHR.START, CHR.END)
  results <- rbind(results, foverlaps(curr.gff, curr.stats.w), use.names=F)
  print(i)
}

colnames(results) <- colnames(foverlaps(curr.gff, curr.stats.w))

mean.weights <- apply(na.omit(results[,10:24]), 2, mean)
names(mean.weights) <- gsub("topo", "T", names(mean.weights))

color.topo <- c("gray", "gray", "gray", "gray", "gray", "white", "red", "red", "black", "gray",
                "white", "yellow", "yellow", "white", "gray")

barplot(mean.weights[order(mean.weights, decreasing = T)], ylab="Mean weight", xlab="Topology", col=color.topo[order(mean.weights, decreasing=T)])
legend("topright", legend=c("geography", "color (red)", "color (yellow)", "color (both)", "others"), 
       fill=c("gray", "red", "yellow", "black","white"))

rft.mt <- mean(na.omit(apply(results[,16:18], 1, sum)))

w.files.all <- list.files("../TWISST_weights/", pattern = "SUPER_[0-9]+_all_W")
w.all <- list()
w.all$w <- lapply(paste("../TWISST_weights/", w.files.all, sep=""), read.table, header=TRUE)
weights.combined.all <- rbindlist(w.all$w, use.names = T)
weights.combined.all <- na.omit(weights.combined.all)
weights.combined.all <- prop.weight(weights.combined.all)

av.w.r <- c()
seeds <- c()
for(i in 1:5000){
  seed <- runif(1,0,6000)
  seeds <- c(seeds, seed)
  set.seed(seed)
  curr.w <- weights.combined.all[sample(nrow(weights.combined.all),283), ]
  av.w.r <- c(av.w.r, mean(apply(curr.w[,7:9], 1, sum)))
  print(i)
}

plot(density(av.w.r), lwd=2, col="red", xlab="Average weight RFT topologies", main="", xlim=c(0,0.2))
abline(v=rft.mt, lwd=2)

mean(av.w.r)
quantile(av.w.r, 0.95)

# cline width SA

setwd("../HZAR_single_SNPs_South_Africa/")
Results.files <- list.files(".", pattern="Estimates")
lbl.results <- do.call(rbind, lapply(Results.files, read.table))
lbl.results <- na.omit(lbl.results)

fixed.loci.names <- colnames(lbl.results)

loci.coord <- read.table("SNPs_South_Africa.plink.map")
fixed.loci.coord <- loci.coord[loci.coord$V2 %in% lbl.results$locus, ]
fixed.loci.coord <- fixed.loci.coord[,c(1,4,2)]
colnames(fixed.loci.coord) <- c("CHR", "pos", "locus")
lbl.loci.table <- merge(lbl.results, fixed.loci.coord)
lbl.loci.table <- lbl.loci.table[grep("SUPER", lbl.loci.table$CHR), ]
lbl.loci.table <- lbl.loci.table[lbl.loci.table$CHR != "SUPER_Z", ]
lbl.loci.table <- lbl.loci.table[lbl.loci.table$CHR != "SUPER_W", ]

mt.wi <- c()
for(i in unique(lbl.loci.table$CHR)){
  curr.gff <- genes.gff.aut[genes.gff.aut$chr == i, c(4,5,10), ]
  curr.gff$start <- curr.gff$start - 50000
  curr.gff$end <- curr.gff$end + 50000
  curr.lbl <- lbl.loci.table[lbl.loci.table$CHR == i, ]
  for(j in 1:nrow(curr.lbl)){
    if(length(which(between(curr.lbl$pos[j], curr.gff$start, curr.gff$end))) > 0){
      mt.wi <- c(mt.wi, curr.lbl$width[j])
    }
  }
}

av.wi.r <- c()
seeds <- c()
for(i in 1:5000){
  seed <- runif(1,0,6000)
  seeds <- c(seeds, seed)
  set.seed(seed)
  av.wi.r <- c(av.wi.r, mean(lbl.loci.table$width[sample(nrow(lbl.loci.table),38)]))
  print(i)
}

plot(density(av.wi.r), lwd=2, col="red", xlab="Mean cline width (km)", main="")
abline(v=mean(mt.wi), lwd=2)

# cline width U

setwd("../HZAR_single_SNPs_Uganda-Kenya/")
Results.files <- list.files(".", pattern="Estimates")
lbl.results <- do.call(rbind, lapply(Results.files, read.table))
lbl.results <- na.omit(lbl.results)

fixed.loci.names <- colnames(lbl.results)

loci.coord <- read.table("SNPs_Uganda-Kenya.plink.map")
fixed.loci.coord <- loci.coord[loci.coord$V2 %in% lbl.results$locus, ]
fixed.loci.coord <- fixed.loci.coord[,c(1,4,2)]
colnames(fixed.loci.coord) <- c("CHR", "pos", "locus")
lbl.loci.table <- merge(lbl.results, fixed.loci.coord)
lbl.loci.table <- lbl.loci.table[grep("SUPER", lbl.loci.table$CHR), ]
lbl.loci.table <- lbl.loci.table[lbl.loci.table$CHR != "SUPER_Z", ]
lbl.loci.table <- lbl.loci.table[lbl.loci.table$CHR != "SUPER_W", ]

mt.wi <- c()
for(i in unique(lbl.loci.table$CHR)){
  curr.gff <- genes.gff.aut[genes.gff.aut$chr == i, c(4,5,10), ]
  curr.gff$start <- curr.gff$start - 50000
  curr.gff$end <- curr.gff$end + 50000
  curr.lbl <- lbl.loci.table[lbl.loci.table$CHR == i, ]
  for(j in 1:nrow(curr.lbl)){
    if(length(which(between(curr.lbl$pos[j], curr.gff$start, curr.gff$end))) > 0){
      mt.wi <- c(mt.wi, curr.lbl$width[j])
    }
  }
}

av.wi.r <- c()
seeds <- c()
for(i in 1:5000){
  seed <- runif(1,0,6000)
  seeds <- c(seeds, seed)
  set.seed(seed)
  av.wi.r <- c(av.wi.r, mean(lbl.loci.table$width[sample(nrow(lbl.loci.table),80)]))
  print(i)
}

hist(av.wi.r)
plot(density(av.wi.r), lwd=2, col="red", xlab="Mean cline width (km)", main="")
abline(v=mean(mt.wi), lwd=2)

