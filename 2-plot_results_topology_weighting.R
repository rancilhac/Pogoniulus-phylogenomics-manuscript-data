library(ggplot2)
library(ape)
library(data.table)

setwd("TWISST_weights")

prop.weight <- function(raw){raw/apply(raw, 1, sum)}

## read autosomal data
w.files.all <- list.files(".", pattern = "SUPER_[0-9]+_all_W")
w.all <- list()
w.all$w <- lapply(w.files.all, read.table, header=TRUE)
weights.combined.all <- rbindlist(w.all$w, use.names = T)
weights.combined.all <- na.omit(weights.combined.all)
weights.combined.all <- prop.weight(weights.combined.all)
colnames(weights.combined.all) <- gsub("topo", "T", colnames(weights.combined.all))

## read Z chromosome data
weights.Z <- read.table("SUPER_Z_mf_all_W", header=T)
weights.Z <- na.omit(weights.Z)
weights.Z <- prop.weight(weights.Z)
colnames(weights.Z) <- gsub("topo", "T", colnames(weights.Z))

## Plot the proportion of windows with w=1 on autosomes and Z chromosome (Figure 3a) 
weight.threshold <- function(weights){
  length(which(weights == 1))
}

prop.1.aut <- apply(weights.combined.all, 2, weight.threshold)/sum(apply(weights.combined.all, 2, weight.threshold))
prop.1.Z <- apply(weights.Z, 2, weight.threshold)/sum(apply(weights.Z, 2, weight.threshold))

aut.tab <- data.frame(group=rep("autosomes", each=15)[order(prop.1.aut, decreasing=T)],
                      topo=c(1:15)[order(prop.1.aut, decreasing=T)],
                      weights=prop.1.aut[order(prop.1.aut, decreasing=T)],
                      prop=prop.1.aut[order(prop.1.aut, decreasing=T)])

Z.tab <- data.frame(group=rep("chr Z", each=15)[order(prop.1.aut, decreasing=T)],
                      topo=c(1:15)[order(prop.1.aut, decreasing=T)],
                      weights=prop.1.Z[order(prop.1.aut, decreasing=T)],
                      prop=prop.1.Z[order(prop.1.aut, decreasing=T)])

Z.tab$prop <- -Z.tab$prop

plot.tab <- rbind(aut.tab, Z.tab)

ggplot(plot.tab, aes(x=factor(topo, level=order(prop.1.aut, decreasing=T)), y=prop, fill=group)) +
  geom_bar(stat="identity", position="identity") +
  xlab("topology") +
  ylab("proportion of windows with w=1") +
  geom_text(aes(label=round(weights, digits=3)), vjust=0) + scale_fill_manual(values=c("brown1", "cadetblue"))

## Plot the proportion of windows with w=1 on autosomes for four sample subsets (Figure 3b and S2)

## Data "core range" samples
w.files.allo <- list.files(".", pattern = "SUPER_[0-9]+_allopatric_W")
w.allo <- list()
w.allo$w <- lapply(w.files.allo, read.table, header=TRUE)
weights.combined.allo <- rbindlist(w.allo$w, use.names = T)
weights.combined.allo <- na.omit(weights.combined.allo)
weights.combined.allo <- prop.weight(weights.combined.allo)

prop.1.allo <- apply(weights.combined.allo, 2, weight.threshold)/sum(apply(weights.combined.allo, 2, weight.threshold))
names(prop.1.allo) <- c(1:15)

## Data "sympatric" samples
w.files.sym <- list.files(".", pattern = "SUPER_[0-9]+_sympatric_W")
w.sym <- list()
w.sym$w <- lapply(w.files.sym, read.table, header=TRUE)
weights.combined.sym <- rbindlist(w.sym$w, use.names = T)
weights.combined.sym <- na.omit(weights.combined.sym)
weights.combined.sym <- prop.weight(weights.combined.sym)

prop.1.sym <- apply(weights.combined.sym, 2, weight.threshold)/sum(apply(weights.combined.sym, 2, weight.threshold))
names(prop.1.sym) <- c(1:15)

## Data "control" samples
w.files.control <- list.files(".", pattern = "SUPER_[0-9]+_control_W")
w.control <- list()
w.control$w <- lapply(w.files.control, read.table, header=TRUE)
weights.combined.control <- rbindlist(w.control$w, use.names = T)
weights.combined.control <- na.omit(weights.combined.control)
weights.combined.control <- prop.weight(weights.combined.control)

prop.1.control <- apply(weights.combined.control, 2, weight.threshold)/sum(apply(weights.combined.control, 2, weight.threshold))
names(prop.1.control) <- c(1:15)

## Plot

color.topo <- c("gray", "gray", "gray", "gray", "gray", "gray", "red", "red", "black", "gray",
                "white", "yellow", "yellow", "white", "gray")

par(mfrow=c(2,2))

names(prop.1.aut) <- c(1:15)
barplot(prop.1.aut[order(prop.1.aut, decreasing=T)], xlab="topology", ylab="proportion of windows with w=1", 
        col=color.topo[order(prop.1.aut, decreasing=T)], main="all samples")

#legend("topright", legend=c("geography", "color (red)", "color (yellow)", "color (both)", "others"), 
       #fill=c("gray", "red", "yellow", "black","white"))

barplot(prop.1.sym[order(prop.1.sym, decreasing=T)], xlab="topology", ylab="proportion of windows with w=1", 
        col=color.topo[order(prop.1.sym, decreasing=T)], main="'sympatric' samples")

barplot(prop.1.allo[order(prop.1.allo, decreasing=T)], xlab="topology", ylab="proportion of windows with w=1", 
        col=color.topo[order(prop.1.allo, decreasing=T)], main="'core range' samples")

barplot(prop.1.control[order(prop.1.control, decreasing=T)], xlab="topology", ylab="proportion of windows with w=1", 
        col=color.topo[order(prop.1.control, decreasing=T)], main="'control' samples")

## Plot the proportion of windows with w=1 on Z chromosome for four sample subsets (Figure S3)

## Data 'sympatric' samples

weights.Z.sym <- read.table("SUPER_Z_mf_sympatric_W", header=T)
weights.Z.sym <- na.omit(weights.Z.sym)
weights.Z.sym <- prop.weight(weights.Z.sym)

prop.1.Z.sym <- apply(weights.Z.sym, 2, weight.threshold)/sum(apply(weights.Z.sym, 2, weight.threshold))
names(prop.1.Z.sym) <- c(1:15)

## Data 'core range' samples
weights.Z.allo <- read.table("SUPER_Z_mf_allopatric_W", header=T)
weights.Z.allo <- na.omit(weights.Z.allo)
weights.Z.allo <- prop.weight(weights.Z.allo)

prop.1.Z.allo <- apply(weights.Z.allo, 2, weight.threshold)/sum(apply(weights.Z.allo, 2, weight.threshold))
names(prop.1.Z.allo) <- c(1:15)

## Data 'control' samples
weights.Z.control <- read.table("SUPER_Z_mf_control_W", header=T)
weights.Z.control <- na.omit(weights.Z.control)
weights.Z.control <- prop.weight(weights.Z.control)

prop.1.Z.control <- apply(weights.Z.control, 2, weight.threshold)/sum(apply(weights.Z.control, 2, weight.threshold))
names(prop.1.Z.control) <- c(1:15)

par(mfrow=c(2,2))
names(prop.1.Z) <- c(1:15)
barplot(prop.1.Z[order(prop.1.Z, decreasing=T)], col=color.topo[order(prop.1.Z, decreasing=T)],
        ylab="proportion of windows with weight=1", xlab="topology", main="all samples")

#legend("topright", legend=c("geography", "color (red)", "color (yellow)", "color (both)", "others"), 
       #fill=c("gray", "red", "yellow", "black","white"))

barplot(prop.1.Z.sym[order(prop.1.Z.sym, decreasing=T)], col=color.topo[order(prop.1.Z.sym, decreasing=T)],
        ylab="proportion of windows with weight=1", xlab="topology", main="'sympatric' samples")

barplot(prop.1.Z.allo[order(prop.1.Z.allo, decreasing=T)], col=color.topo[order(prop.1.Z.allo, decreasing=T)],
        ylab="proportion of windows with weight=1", xlab="topology", main="'core range' samples")

barplot(prop.1.Z.control[order(prop.1.Z.control, decreasing=T)], col=color.topo[order(prop.1.Z.control, decreasing=T)],
        ylab="proportion of windows with weight=1", xlab="topology", main="'control' samples")
