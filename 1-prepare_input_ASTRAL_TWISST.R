library(ape)

setwd("Neighbour_Joining_trees")

thin.pos <- function(pos, start, interval){
  thinned.pos <- c(start)
  p <- pos[start]
  i <- 2
  while(i <= length(pos)){
    if(pos[i]-p >= interval){ thinned.pos <- c(thinned.pos, i)
    p <- pos[i]
    i <- i+1 }
    else if(pos[i]-p < interval){ i <- i+1 }
  }
  return(thinned.pos)
}

thin.ranges <- function(ranges, increment){
  thinned.pos <- c(1)
  i <- 2
  end.pos <- ranges[1,2]
  while(i <= nrow(ranges)){
    if(ranges[i,1] >= end.pos+increment){ thinned.pos <- c(thinned.pos, i)
    i <- i+1
    end.pos <- ranges[i-1,2]}
    else {i <- i+1}
  }
  return(thinned.pos)
}

## Autosomes

for(i in c(1:44)){
  
  trees.name.all <- paste("SUPER_",i,"_all_500S_maxmiss60_NJ_trees.trees", sep="")
  trees.name.allo <- paste("SUPER_",i,"_allopatric_500S_maxmiss60_NJ_trees.trees", sep="")
  trees.name.sym <- paste("SUPER_",i,"_sympatric_500S_maxmiss60_NJ_trees.trees", sep="")
  trees.name.ran <- paste("SUPER_",i,"_control2_500S_maxmiss60_NJ_trees.trees", sep="")
  stats.name.all <- paste("SUPER_",i,"_all_500S_maxmiss60_windows_stats.tsv", sep="")
  stats.name.allo <- paste("SUPER_", i, "_allopatric_500S_maxmiss60_windows_stats.tsv", sep="")
  stats.name.sym <- paste("SUPER_", i, "_sympatric_500S_maxmiss60_windows_stats.tsv", sep="")
  stats.name.ran <- paste("SUPER_", i, "_control2_500S_maxmiss60_windows_stats.tsv", sep="")
  
  
  trees.all <- read.tree(trees.name.all)
  trees.allo <- read.tree(trees.name.allo)
  trees.sym <- read.tree(trees.name.sym)
  trees.ran <- read.tree(trees.name.ran)
  stats.all <- read.table(stats.name.all, header=T)
  stats.allo <- read.table(stats.name.allo, header=T)
  stats.sym <- read.table(stats.name.sym, header=T)
  stats.ran <- read.table(stats.name.ran, header=T)
  
  stats.all <- stats.all[!is.na(stats.all$TREE), ]
  row.names(stats.all) <- c(1:nrow(stats.all))
  stats.miss <- stats.all[which(stats.all$PROP.MISS < 0.4 & stats.all$CHR.SIZE < 50000), ]
  stats.astral <- stats.miss[thin.ranges(stats.miss[,2:3], 10000), ]
  stats.twisst <- stats.miss[stats.miss$NTIPS == 45, ]
  trees.all.astral <- trees.all[as.numeric(row.names(stats.astral))]
  trees.all.twisst <- trees.all[as.numeric(row.names(stats.twisst))]
  write.tree(trees.all.astral, paste("SUPER_",i,"_all_500S_maxmiss60_NJ_ASTRAL.trees", sep=""))
  write.tree(trees.all.twisst, paste("SUPER_",i,"_all_500S_maxmiss60_NJ_TWISST.trees", sep=""))
  
  stats.allo <- stats.allo[!is.na(stats.allo$TREE), ]
  row.names(stats.allo) <- c(1:nrow(stats.allo))
  stats.miss <- stats.allo[which(stats.allo$PROP.MISS < 0.4 & stats.allo$CHR.SIZE < 50000), ]
  stats.astral <- stats.miss[thin.ranges(stats.miss[,2:3], 10000), ]
  stats.twisst <- stats.miss[stats.miss$NTIPS == 13, ]
  trees.allo.astral <- trees.allo[as.numeric(row.names(stats.astral))]
  trees.allo.twisst <- trees.allo[as.numeric(row.names(stats.twisst))]
  write.tree(trees.allo.astral, paste("SUPER_",i,"_allopatric_500S_maxmiss60_NJ_ASTRAL.trees", sep=""))
  write.tree(trees.allo.twisst, paste("SUPER_",i,"_allopatric_500S_maxmiss60_NJ_TWISST.trees", sep=""))
  
  stats.sym <- stats.sym[!is.na(stats.sym$TREE), ]
  row.names(stats.sym) <- c(1:nrow(stats.sym))
  stats.miss <- stats.sym[which(stats.sym$PROP.MISS < 0.4 & stats.sym$CHR.SIZE < 50000), ]
  stats.astral <- stats.miss[thin.ranges(stats.miss[,2:3], 10000), ]
  stats.twisst <- stats.miss[stats.miss$NTIPS == 34, ]
  trees.sym.astral <- trees.sym[as.numeric(row.names(stats.astral))]
  trees.sym.twisst <- trees.sym[as.numeric(row.names(stats.twisst))]
  write.tree(trees.sym.astral, paste("SUPER_",i,"_sympatric_500S_maxmiss60_NJ_ASTRAL.trees", sep=""))
  write.tree(trees.sym.twisst, paste("SUPER_",i,"_sympatric_500S_maxmiss60_NJ_TWISST.trees", sep=""))
  
  stats.ran <- stats.ran[!is.na(stats.ran$TREE), ]
  row.names(stats.ran) <- c(1:nrow(stats.ran))
  stats.miss <- stats.ran[which(stats.ran$PROP.MISS < 0.4 & stats.ran$CHR.SIZE < 50000), ]
  stats.astral <- stats.miss[thin.ranges(stats.miss[,2:3], 10000), ]
  stats.twisst <- stats.miss[stats.miss$NTIPS == 13, ]
  trees.ran.astral <- trees.ran[as.numeric(row.names(stats.astral))]
  trees.ran.twisst <- trees.ran[as.numeric(row.names(stats.twisst))]
  write.tree(trees.ran.astral, paste("SUPER_",i,"_control2_500S_maxmiss60_NJ_ASTRAL.trees", sep=""))
  write.tree(trees.ran.twisst, paste("SUPER_",i,"_control2_500S_maxmiss60_NJ_TWISST.trees", sep=""))
  
  
  print(i)
}

## Z chromosome

trees.name.all <- "SUPER_Z_mf_all_500S_maxmiss60_NJ_trees.trees"
trees.name.allo <- "SUPER_Z_mf_allopatric_500S_maxmiss60_NJ_trees.trees"
trees.name.sym <- "SUPER_Z_mf_sympatric_500S_maxmiss60_NJ_trees.trees"
trees.name.ran <- "SUPER_Z_mf_control2_500S_maxmiss60_NJ_trees.trees"
stats.name.all <- "SUPER_Z_mf_all_500S_maxmiss60_windows_stats.tsv"
stats.name.allo <- "SUPER_Z_mf_allopatric_500S_maxmiss60_windows_stats.tsv"
stats.name.sym <- "SUPER_Z_mf_sympatric_500S_maxmiss60_windows_stats.tsv"
stats.name.ran <- "SUPER_Z_mf_control2_500S_maxmiss60_windows_stats.tsv"


trees.all <- read.tree(trees.name.all)
trees.allo <- read.tree(trees.name.allo)
trees.sym <- read.tree(trees.name.sym)
trees.ran <- read.tree(trees.name.ran)
stats.all <- read.table(stats.name.all, header=T)
stats.allo <- read.table(stats.name.allo, header=T)
stats.sym <- read.table(stats.name.sym, header=T)
stats.ran <- read.table(stats.name.ran, header=T)


stats.all <- stats.all[!is.na(stats.all$TREE), ]
row.names(stats.all) <- c(1:nrow(stats.all))
stats.miss <- stats.all[which(stats.all$PROP.MISS < 0.4 & stats.all$CHR.SIZE < 50000), ]
stats.astral <- stats.miss[thin.ranges(stats.miss[,2:3], 10000), ]
stats.twisst <- stats.miss[stats.miss$NTIPS == 45, ]
trees.all.astral <- trees.all[as.numeric(row.names(stats.astral))]
trees.all.twisst <- trees.all[as.numeric(row.names(stats.twisst))]
write.tree(trees.all.astral, paste("SUPER_",i,"_mf_all_500S_maxmiss60_NJ_ASTRAL.trees", sep=""))
write.tree(trees.all.twisst, paste("SUPER_",i,"_mf_all_500S_maxmiss60_NJ_TWISST.trees", sep=""))

stats.allo <- stats.allo[!is.na(stats.allo$TREE), ]
row.names(stats.allo) <- c(1:nrow(stats.allo))
stats.miss <- stats.allo[which(stats.allo$PROP.MISS < 0.4 & stats.allo$CHR.SIZE < 50000), ]
stats.astral <- stats.miss[thin.ranges(stats.miss[,2:3], 10000), ]
stats.twisst <- stats.miss[stats.miss$NTIPS == 13, ]
trees.allo.astral <- trees.allo[as.numeric(row.names(stats.astral))]
trees.allo.twisst <- trees.allo[as.numeric(row.names(stats.twisst))]
write.tree(trees.allo.astral, paste("SUPER_",i,"_mf_allopatric_500S_maxmiss60_NJ_ASTRAL.trees", sep=""))
write.tree(trees.allo.twisst, paste("SUPER_",i,"_mf_allopatric_500S_maxmiss60_NJ_TWISST.trees", sep=""))

stats.sym <- stats.sym[!is.na(stats.sym$TREE), ]
row.names(stats.sym) <- c(1:nrow(stats.sym))
stats.miss <- stats.sym[which(stats.sym$PROP.MISS < 0.4 & stats.sym$CHR.SIZE < 50000), ]
stats.astral <- stats.miss[thin.ranges(stats.miss[,2:3], 10000), ]
stats.twisst <- stats.miss[stats.miss$NTIPS == 34, ]
trees.sym.astral <- trees.sym[as.numeric(row.names(stats.astral))]
trees.sym.twisst <- trees.sym[as.numeric(row.names(stats.twisst))]
write.tree(trees.sym.astral, paste("SUPER_",i,"_mf_sympatric_500S_maxmiss60_NJ_ASTRAL.trees", sep=""))
write.tree(trees.sym.twisst, paste("SUPER_",i,"_mf_sympatric_500S_maxmiss60_NJ_TWISST.trees", sep=""))

stats.ran <- stats.ran[!is.na(stats.ran$TREE), ]
row.names(stats.ran) <- c(1:nrow(stats.ran))
stats.miss <- stats.ran[which(stats.ran$PROP.MISS < 0.4 & stats.ran$CHR.SIZE < 50000), ]
stats.astral <- stats.miss[thin.ranges(stats.miss[,2:3], 10000), ]
stats.twisst <- stats.miss[stats.miss$NTIPS == 13, ]
trees.ran.astral <- trees.ran[as.numeric(row.names(stats.astral))]
trees.ran.twisst <- trees.ran[as.numeric(row.names(stats.twisst))]
write.tree(trees.ran.astral, paste("SUPER_",i,"_mf_control2_500S_maxmiss60_NJ_ASTRAL.trees", sep=""))
write.tree(trees.ran.twisst, paste("SUPER_",i,"_mf_control2_500S_maxmiss60_NJ_TWISST.trees", sep=""))