library(ape)

## define accessory fuctions

prop.weight <- function(raw){raw/apply(raw, 1, sum)}

weight.threshold <- function(weights){
  length(which(weights == 1))
}

## read Topology weighting data

setwd("Neighbour_Joining_trees/")

stats.files <- list.files(".", pattern = "all_.*.tsv")
stats <- list()
stats$stats <- lapply(stats.files, read.table, header=TRUE)
stats.combined <- rbindlist(stats$stats, use.names = T)
stats.combined <- stats.combined[!stats.combined$CHR == "SUPER_Z", ]

stats.combined <- stats.combined[stats.combined$NTIPS == 45, ]
stats.combined <- stats.combined[which(stats.combined$PROP.MISS < 0.8 & stats.combined$CHR.SIZE < 50000), ]
stats.combined <- na.omit(stats.combined)

setwd("../TWISST_weights/")

w.files.all <- list.files(".", pattern = "SUPER_[0-9]+_all_W")
w.all <- list()
w.all$w <- lapply(w.files.all, read.table, header=TRUE)
weights.combined.all <- rbindlist(w.all$w, use.names = T)
weights.combined.all <- na.omit(weights.combined.all)
weights.combined.all <- prop.weight(weights.combined.all)

## read linkage map

linkage.map <- read.table("../Linkage_map_YFT_LDhat_100kb.txt")
colnames(linkage.map) <- c("chr", "start", "median", "end", "n_variants", "rho")
linkage.map$end <- linkage.map$end*1000
linkage.map$start <- linkage.map$start*1000
linkage.map$median <- linkage.map$median*1000

chromosome <- end.wind <- med.pos <- rho.wind <- c()

for(c in unique(stats.combined$CHR)){
  curr.stats <- stats.combined[stats.combined$CHR == c, ]
  curr.linkmap <- linkage.map[linkage.map$chr == c, ]
  for(i in 1:nrow(curr.stats)){
    MED <- (curr.stats$CHR.END[i]+curr.stats$CHR.START[i])/2
    rec <- which(curr.linkmap$end > MED)[1]
    chromosome <- c(chromosome, c)
    med.pos <- c(med.pos, MED)
    rho.wind <- c(rho.wind, curr.linkmap$rho[rec])
    end.wind <- c(end.wind, curr.linkmap$end[rec])
  }
  print(c)
}

rec.table <- as.data.frame(cbind(chromosome, med.pos, end.wind, rho.wind))
colnames(rec.table) <- c("CHR", "WIN.TW", "END.R", "rho")
rec.table$rho <- as.numeric(rec.table$rho)
write.table(rec.table,"recombination_500S_windows.txt", quote = F, row.names = F)


q005 <- quantile(rec.table$rho, 0.05)
q02 <- quantile(rec.table$rho, 0.2)
q08 <- quantile(rec.table$rho, 0.8)
q095 <- quantile(rec.table$rho, 0.95)

## extract the proportion of windows with weight = 1.0 in different categories of recombination rate
weights.q005 <- weights.combined.all[rec.table$rho < q005, ]
w1.q005 <- apply(weights.q005, 2, weight.threshold)/sum(apply(weights.q005, 2, weight.threshold))
weights.q02 <- weights.combined.all[rec.table$rho > q005 & rec.table$rho < q02, ]
w1.q02 <- apply(weights.q02, 2, weight.threshold)/sum(apply(weights.q02, 2, weight.threshold))
weights.q08 <- weights.combined.all[rec.table$rho > q02 & rec.table$rho < q08, ]
w1.q08 <- apply(weights.q08, 2, weight.threshold)/sum(apply(weights.q08, 2, weight.threshold))
weights.q095 <- weights.combined.all[rec.table$rho > q08 & rec.table$rho < q095, ]
w1.q095 <- apply(weights.q095, 2, weight.threshold)/sum(apply(weights.q095, 2, weight.threshold))
weights.q1 <- weights.combined.all[rec.table$rho > q095, ]
w1.q1 <- apply(weights.q1, 2, weight.threshold)/sum(apply(weights.q1, 2, weight.threshold))

w1.rec <- rbind(w1.q005, w1.q02, w1.q08, w1.q095)
colnames(w1.rec) <- 1:15

## plot weights

par(mfrow=c(1,1))
barplot(height=w1.rec, beside = T, xlab="Topology", ylab="Average weight")



## Extract trees from non-adjacent windows with low recombination rate (< 5% quantile of rho distribution)

nadj.windows <- stats.combined[seq(1, nrow(stats.combined), by=2), ]
nadj.rec <- rec.table[seq(1, nrow(stats.combined), by=2), ]
nadj.LR <- nadj.rec[nadj.rec$rho < q005, ]

trees <- read.tree("all_autosomes_NJ_500S.trees")

trees.LR <- trees[as.numeric(row.names(nadj.LR))]
write.tree(trees.LR, "all_autosomes_NJ_500S_low_recombination.trees")
