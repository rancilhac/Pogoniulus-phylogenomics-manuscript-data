library(stringr)
library(ggplot2)

setwd("HZAR_single_SNPs_Uganda-Kenya/")

Results.files <- list.files(".", pattern="Estimates")
lbl.results <- do.call(rbind, lapply(Results.files, read.table))
lbl.results <- na.omit(lbl.results)

fixed.loci.names <- colnames(lbl.results)

loci.coord <- read.table("./SNPs_Uganda-Kenya.plink.map")
fixed.loci.coord <- loci.coord[loci.coord$V2 %in% lbl.results$locus, ]
fixed.loci.coord <- fixed.loci.coord[,c(1,4,2)]
colnames(fixed.loci.coord) <- c("CHR", "pos", "locus")
lbl.loci.table <- merge(lbl.results, fixed.loci.coord)
lbl.loci.table <- lbl.loci.table[grep("SUPER", lbl.loci.table$CHR), ]
lbl.loci.table <- lbl.loci.table[lbl.loci.table$CHR != "SUPER_Z", ]
lbl.loci.table <- lbl.loci.table[lbl.loci.table$CHR != "SUPER_W", ]
lbl.loci.table <- lbl.loci.table[lbl.loci.table$CHR != "SUPER_W_unloc_1", ]

## Execute the next two lines to remove SNPs with outlier clines
lbl.loci.table <- lbl.loci.table[abs(lbl.loci.table$center) < 100, ]
lbl.loci.table <- lbl.loci.table[lbl.loci.table$width < 300, ]

ggplot(lbl.loci.table, aes(x=width)) + geom_density(fill="red", colour="red") + theme_classic()
ggplot(lbl.loci.table, aes(x=width)) + geom_density(fill="red", colour="red") + xlim(0,100)+ theme_classic()

prop.weight <- function(raw){raw/apply(raw, 1, sum)}

table.rec <- read.table("../TWISST_weights/recombination_500S_windows.txt", header=T)

results <- matrix(nrow = 0, ncol=25)

for(i in unique(lbl.loci.table$CHR)){
  curr.loci <- lbl.loci.table[lbl.loci.table$CHR == i, ]
  curr.rec <- table.rec[table.rec$CHR == i, ]
  curr.TW <- paste("~/Work/Tinkerbirds_Cyprus/Phylogeny_chrysopus/Trees_sliding_windows/windows_500S/new_trees/",i, "_all_W", sep="")
  TW <- read.table(curr.TW, header=T)
  TW <- prop.weight(TW)
  curr.stats <- paste("~/Work/Tinkerbirds_Cyprus/Phylogeny_chrysopus/Trees_sliding_windows/windows_500S/new_trees/",i, "_all_500S_maxmiss60_windows_stats.tsv", sep="")
  stats <- read.table(curr.stats, header=T)
  stats <- stats[stats$PROP.MISS < 0.4 & stats$CHR.SIZE < 50000 & stats$NTIPS == 45, ]
  stats <- na.omit(stats)
  pos <- c()
  rec <- c()
  w <- c()
  for(l in 1:length(curr.loci$pos)){
    w <- c(w, which(stats$CHR.START <= curr.loci$pos[l] & stats$CHR.END >= curr.loci$pos[l]))
    if(length(which(stats$CHR.START <= curr.loci$pos[l] & stats$CHR.END >= curr.loci$pos[l])) > 0){
      rec <- c(rec, curr.rec[which(curr.rec$END.R > curr.loci$pos[l]), ][1,4])
      pos <- c(pos, l)}
  }
  
  tab <- cbind(curr.loci[pos, ], stats[w,1:3], TW[w, ], rec)
  results <- rbind(results, tab)
  
}

plot(results$rec, results$width, pch=19, xlab="Rho", ylab="Cline width", main="chrysoconus/affinis")
reg <- lm(results$width~results$rec)
abline(reg, lwd=2)
cor.test(results$rec, results$width, method="spearman")


boxplot(results$width[results$topo1 > 0.8], results$width[results$topo2 > 0.8], results$width[results$topo3 > 0.8], 
        results$width[results$topo7 > 0.8], results$width[results$topo8 > 0.8], results$width[results$topo9 > 0.8],
        results$width[results$topo4 > 0.8],results$width[results$topo5 > 0.8],
        results$width[results$topo6 > 0.8], results$width[results$topo10 > 0.8], 
        results$width[results$topo11 > 0.8], results$width[results$topo12 > 0.8], results$width[results$topo13 > 0.8],
        results$width[results$topo14 > 0.8], results$width[results$topo15 > 0.8],
        col=c("blue", "blue", "blue", "red", "red", "red", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey"),
        xax="", xlab="Topology", ylab="Cline width", pch=19, cex=0.8, ylim=c(0,100))
axis(1, at=c(1:15), labels = c("T1", "T2", "T3", "T7", "T8", "T9", "T4", "T5", "T6", "T10", "T11", "T12", "T13", "T14", "T15"))
legend("topright", c("chrysoconus/affinis", "affinis/pusillus", "others"), fill=c("blue", "red", "grey"))


rest <- c(results$width[results$topo7 > 0.8], results$width[results$topo8 > 0.8], results$width[results$topo9 > 0.8],
          results$width[results$topo4 > 0.8],results$width[results$topo5 > 0.8],
          results$width[results$topo6 > 0.8], results$width[results$topo10 > 0.8], 
          results$width[results$topo11 > 0.8], results$width[results$topo12 > 0.8], results$width[results$topo13 > 0.8],
          results$width[results$topo14 > 0.8], results$width[results$topo15 > 0.8])
chry.aff <- c(results$width[results$topo1 > 0.8], results$width[results$topo2 > 0.8], results$width[results$topo3 > 0.8])
pus.aff <- c(results$width[results$topo7 > 0.8], results$width[results$topo8 > 0.8], results$width[results$topo9 > 0.8])

wil.r <- wilcox.test(chry.aff, rest, alternative="great")
wil.re <- wilcox.test(chry.aff, pus.aff, alternative="great")
boxplot(chry.aff, rest, pus.aff, xax="", ylab="Cline width", pch=19, ylim=c(0,100))
axis(1, at=c(1:3), labels = c("extoni/pusillus", "all others", "pusillus/affinis"))
