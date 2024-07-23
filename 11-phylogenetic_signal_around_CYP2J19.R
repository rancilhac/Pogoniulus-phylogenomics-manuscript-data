library(ggplot2)
library(ape)
library(data.table)

prop.weight <- function(raw){raw/apply(raw, 1, sum)}

chr8.w <- read.table("TWISST_weights/SUPER_8_all_W", header=T)
chr8.w <- prop.weight(chr8.w)

stats.chr8 <- read.table("Neighbour_Joining_trees/SUPER_8_all_500S_maxmiss60_windows_stats.tsv", header=T)
stats.chr8 <- stats.chr8[stats.chr8$NTIPS == 45, ]
stats.chr8 <- stats.chr8[stats.chr8$PROP.MISS < 0.4 & stats.chr8$CHR.SIZE < 50000, ]

trees.chr8 <- read.tree("Neighbour_Joining_trees/SUPER_8_all_500S_maxmiss60_NJ_TWISST.trees")

#GWAS <- read.csv("East and SouthernAfrica_cordblasts_12Feb24Bonferroni.csv", header=T, sep=",")
#GWAS <- GWAS[GWAS$chr == "SUPER_8", ]

CYP.w <- chr8.w[stats.chr8$CHR.START > 2.2e+6 & stats.chr8$CHR.END < 2.6e+6, ]
stats.CYP <- stats.chr8[stats.chr8$CHR.START > 2.2e+6 & stats.chr8$CHR.END < 2.6e+6, ]
CYP.w.red <- apply(CYP.w[,7:9], 1, sum)
stats.CYP$W.r <- CYP.w.red
stats.CYP$CHR.MID <- stats.CYP$CHR.START + (stats.CYP$CHR.END - stats.CYP$CHR.START)/2

par(mfrow=c(1,1))
plot(stats.CYP$CHR.MID, CYP.w.red, pch=19, cex=1, xlab="Position (chr8)", ylab="Cumulative weight T7,8,9", col="white")
rect(xleft=2377713, ybottom=-0.01, xright=2392689, ytop=1.01, col="gray",border=F)
segments(x0 = stats.CYP$CHR.START, y0 = stats.CYP$W.r, x1 = stats.CYP$CHR.END, y1 = stats.CYP$W.r, col = "black", lwd=2)



