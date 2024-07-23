library(hzar)
library("data.table")

## The script iterates on loci. l.min is the number of the first locus to be processed in the input file and l.max the last one
l.min <- 1
l.max <- 1852

## The next lines can be used to split the SNPs and run different batches in parallel. Here I use batches of 50 SNPs. The number of the start locus (l.min) must be provided as a command line argument
#args = commandArgs(trailingOnly=T)

#l.min <- as.numeric(args[1])
#if(l.min+49 < 1438){ l.max <- l.min+49
#} else { l.max <- 1438 }

setwd("HZAR_single_SNPs_South_Africa")

diff.freq <- function(freqs, p1, p2){
  diff <- freqs[p1]-freqs[p2]
  return(diff)
}

# read the master table (used to extract distance information)

data.SA <- read.csv("../ddRADseq_samples_metadata_South_Africa.csv", sep=";")
data.SA <- na.omit(data.SA)
colnames(data.SA)[1:2] <- c("Sample", "Subspecies")


data.SA$dist_SA[data.SA$Subspecies == "pusillus"] <- -data.SA$dist_SA[data.SA$Subspecies == "pusillus"]


# read and prepare the data
loci <- fread("Allele_frequencies_South_Africa_HZAR.csv", header=T)
indices <- grep(".A", names(loci))
loci.a <- loci[,..indices]
indices.N <- grep(".N", names(loci))
ind <- loci[, ..indices.N]

loci.a <- as.data.frame(loci.a)
ind <- as.data.frame(ind)

# keep loci only for pops with >1 samples

samp.pop <- function(ind){
  if(min(ind) >= 2){ return(TRUE) }
  else{ return(FALSE) }
}

index <- which(apply(ind, 2, samp.pop))
ind <- ind[,index]
loci.a <- loci.a[,index]

#Extract diagnostic SNPs
pops <- loci$Population

transect <- c()
for(i in pops){
  transect <- c(transect, mean(data.SA[data.SA$locality == i,7]))
}

allo.e <- which(pops == "Buffelsdrift")
allo.p <- which(pops == "Lake_Eland")

diff.freqs <- apply(loci.a, 2, diff.freq, p1=allo.e, p2=allo.p)
fixed <- which(diff.freqs == 1)

fixed.loci <- loci.a[,fixed]
fixed.loci.ind <- ind[,fixed]



## Vizualize the data locus by locus
for(i in 1:ncol(fixed.loci)){
  locus <- as.data.frame(cbind(fixed.loci[,i], pops, transect))
  colnames(locus) <- c("freq", "pop", "transect")
  locus$freq <- as.numeric(locus$freq)
  plot(locus$transect, locus$freq, pch=19, main=paste("locus", i))
  readline(prompt="Press [enter] to continue")
}

infer.cline.single.locus <- function(dist, freq, ind){
  analysis <- list()
  analysis$clines <- list()
  analysis$clines$data <- list()
  analysis$clines$models <- list()
  analysis$clines$fit <- list()
  analysis$clines$mcmc <- list()
  analysis$clines$results <- list()
  
  analysis$clines$data <- hzar.doMolecularData1DPops(dist, freq, ind)
  #hzar.plot.obsData(analysis$clines$data)
  
  # model 1: observed pMin & pMax, no tail
  analysis$clines$models[["model1"]] <- hzar.makeCline1DFreq(analysis$clines$data, scaling="fixed",tails="none")
  # model 2: observed pMin & pMax, two tails
  analysis$clines$models[["model2"]] <- hzar.makeCline1DFreq(analysis$clines$data, scaling="fixed",tails="both")
  # model 3: observed pMin & pMax, left tail
  analysis$clines$models[["model3"]] <- hzar.makeCline1DFreq(analysis$clines$data, scaling="fixed",tails="left")
  # model 4: observed pMin & pMax, right tail
  analysis$clines$models[["model4"]] <- hzar.makeCline1DFreq(analysis$clines$data, scaling="fixed",tails="right")
  
  analysis$clines$models <- sapply(analysis$clines$models, hzar.model.addBoxReq,-500, 500,
                                   simplify = F)
  
  analysis$clines$fit$init <- sapply(analysis$clines$models, hzar.first.fitRequest.old.ML,
                                     obsData=analysis$clines$data,verbose=F, simplify=F)
  
  
  seeds <- list(s1=c(54,34,65,38,45,89),
                s2=c(26,23,27,17,54,78),
                s3=c(395,1254,6985,7895,3652,7459),
                s4=c(98,36,59,36,789,12))
  analysis$clines$fit$init$model1$mcmcParam$chainLength <- 5e4
  analysis$clines$fit$init$model1$mcmcParam$burnin <- 1e4
  analysis$clines$fit$init$model1$mcmcParam$seed[[1]] <- seeds$s1
  analysis$clines$fit$init$model2$mcmcParam$chainLength <- 5e4
  analysis$clines$fit$init$model2$mcmcParam$burnin <- 1e4
  analysis$clines$fit$init$model2$mcmcParam$seed[[1]] <- seeds$s2
  analysis$clines$fit$init$model3$mcmcParam$chainLength <- 5e4
  analysis$clines$fit$init$model3$mcmcParam$burnin <- 1e4
  analysis$clines$fit$init$model3$mcmcParam$seed[[1]] <- seeds$s3
  analysis$clines$fit$init$model4$mcmcParam$chainLength <- 5e4
  analysis$clines$fit$init$model4$mcmcParam$burnin <- 1e4
  analysis$clines$fit$init$model4$mcmcParam$seed[[1]] <- seeds$s4
  
  analysis$clines$mcmc$init <- list()
  analysis$clines$mcmc$init$model1 <- hzar.doFit(analysis$clines$fit$init$model1)
  #plot(hzar.mcmc.bindLL(analysis$clines$mcmc$init$model1))
  analysis$clines$mcmc$init$model2 <- hzar.doFit(analysis$clines$fit$init$model2)
  #plot(hzar.mcmc.bindLL(analysis$clines$mcmc$init$model2))
  analysis$clines$mcmc$init$model3 <- hzar.doFit(analysis$clines$fit$init$model3)
  #plot(hzar.mcmc.bindLL(analysis$clines$mcmc$init$model3))
  analysis$clines$mcmc$init$model4 <- hzar.doFit(analysis$clines$fit$init$model4)
  #plot(hzar.mcmc.bindLL(analysis$clines$mcmc$init$model4))
  
  analysis$clines$fit$chains <- lapply(analysis$clines$mcmc$init, hzar.next.fitRequest)
  analysis$clines$fit$chains <- hzar.multiFitRequest(analysis$clines$fit$chains,
                                                     each=3,
                                                     baseSeed=NULL)
  
  center.init <- runif(12,-30,600)
  width.init <- runif(12,0,630)
  deltaL.init <- runif(12,0,630)
  tauL.init <- runif(12,0,1)
  deltaR.init <- runif(12,0,630)
  tauR.init <- runif(12,0,1)
  
  #init model 1 (no tails)
  for(i in 1:3){
    analysis$clines$fit$chains[[i]]$modelParam$init["center"]= center.init[i]
    analysis$clines$fit$chains[[i]]$modelParam$init["width"]= width.init[i]
  }
  # init model 2 (both tails)
  
  for(i in 4:6){
    analysis$clines$fit$chains[[i]]$modelParam$init["center"]= center.init[i]
    analysis$clines$fit$chains[[i]]$modelParam$init["width"]= width.init[i]
    analysis$clines$fit$chains[[i]]$modelParam$init["deltaL"]= deltaL.init[i]
    analysis$clines$fit$chains[[i]]$modelParam$init["tauL"]= tauL.init[i]
    analysis$clines$fit$chains[[i]]$modelParam$init["deltaR"]= deltaR.init[i]
    analysis$clines$fit$chains[[i]]$modelParam$init["tauR"]= tauR.init[i]
  }
  
  # init model 3 (left tail)
  for(i in 7:9){
    analysis$clines$fit$chains[[i]]$modelParam$init["center"]= center.init[i]
    analysis$clines$fit$chains[[i]]$modelParam$init["width"]= width.init[i]
    analysis$clines$fit$chains[[i]]$modelParam$init["deltaL"]= deltaL.init[i]
    analysis$clines$fit$chains[[i]]$modelParam$init["tauL"]= tauL.init[i]
  }
  
  # init model 4 (right tail)
  for(i in 10:12){
    analysis$clines$fit$chains[[i]]$modelParam$init["center"]= center.init[i]
    analysis$clines$fit$chains[[i]]$modelParam$init["width"]= width.init[i]
    analysis$clines$fit$chains[[i]]$modelParam$init["deltaR"]= deltaR.init[i]
    analysis$clines$fit$chains[[i]]$modelParam$init["tauR"]= tauR.init[i]
  }
  
  analysis$clines$mcmc$chains <- hzar.doChain.multi(analysis$clines$fit$chains,
                                                    doPar=F,
                                                    inOrder=FALSE,
                                                    count=3)
  
  analysis$clines$results$initDGs <- list(nullModel = hzar.dataGroup.null(analysis$clines$data))
  analysis$clines$results$initDGs$model1 <- hzar.dataGroup.add(analysis$clines$mcmc$init$model1)
  analysis$clines$results$initDGs$model2 <- hzar.dataGroup.add(analysis$clines$mcmc$init$model2)
  analysis$clines$results$initDGs$model3 <- hzar.dataGroup.add(analysis$clines$mcmc$init$model3)
  analysis$clines$results$initDGs$model4 <- hzar.dataGroup.add(analysis$clines$mcmc$init$model4)
  
  analysis$clines$results$oDG <- hzar.make.obsDataGroup(analysis$clines$results$initDGs)
  analysis$clines$results$oDG <- hzar.copyModelLabels(analysis$clines$results$initDGs,
                                                      analysis$clines$results$oDG)
  analysis$clines$results$oDG <- hzar.make.obsDataGroup(lapply(analysis$clines$mcmc$chains,
                                                               hzar.dataGroup.add),
                                                        analysis$clines$results$oDG)
  
  analysis$clines$results$AICcTable <- hzar.AICc.hzar.obsDataGroup(analysis$clines$results$oDG)
  
  analysis$clines$results$model.name <- rownames(analysis$clines$results$AICcTable)[[ which.min(analysis$clines$results$AICcTable$AICc)]]
  
  analysis$clines$results$model.selected <- analysis$clines$results$oDG$data.groups[[analysis$clines$results$model.name]]
  
  results <- list()
  results[["AIC"]] <- hzar.AICc.hzar.obsDataGroup(analysis$clines$results$oDG)
  results[["Best.Cline"]] <- analysis$clines$results$model.selected
  
  return(results)
}

AIC <- matrix(ncol=6, nrow=50)
best.models <- list()
results <- matrix(ncol=7, nrow=50)

inc <- 1

for(l in l.min:l.max){
  
  curr.name <- unlist(strsplit(colnames(fixed.loci)[l], "\\."))[1]
  AIC[inc,1] <- curr.name
  results[inc,1] <- curr.name
  
  curr.cline <- try(infer.cline.single.locus(transect, fixed.loci[,l], fixed.loci.ind[,l]), silent = T)
  
  if(class(curr.cline) == "try-error"){
    AIC[inc, 2:6] <- rep("NA", 5)
    results[inc, 2:7] <- rep("NA", 6)
    best.models[inc] <- "NA"
  } else{
    AIC[inc, 2:6] <- curr.cline$AIC$AICc
    width <- round(hzar.get.ML.cline(curr.cline$Best.Cline)$param.all$width, 1)
    center <- round(hzar.get.ML.cline(curr.cline$Best.Cline)$param.all$center, 1)
    widthLL <- round(hzar.getLLCutParam(curr.cline$Best.Cline, names(curr.cline$Best.Cline$data.param))[3:4], 1)
    centerLL <- round(hzar.getLLCutParam(curr.cline$Best.Cline, names(curr.cline$Best.Cline$data.param))[1:2], 1)
    results[inc, 2:7] <- unlist(c(width, widthLL, center, centerLL))
    best.models[inc] <- curr.cline$Best.Cline
  }
  
  cat(paste("done ", inc," loci out of ", ncol(fixed.loci), "\n", sep="") ,file=paste("monitor-progress_", l.min,"-", l.max,".txt", sep=""),append=TRUE)
  inc <- inc+1
  
  
}

AIC <- as.data.frame(AIC)
colnames(AIC) <- c("locus", "Null_model", "Model1_no_tails", "Model2_2tails", "Model3_Ltail", "Model4_Rtail")
write.table(AIC, paste("AIC_hzar_SA_lbl_", l.min,"-", l.max,".txt", sep=""))
results <- as.data.frame(results)
colnames(results) <- c("locus","width", "width_LLlow", "width_LLhigh", "center", "center_LLlow", "center_LLhigh")
write.table(results, paste("Estimates_hzar_SA_lbl_", l.min, "-", l.max,".txt", sep=""))
save.image(file=paste("hzar_lbl_SA_", l.min, "-", l.max, "env.RData", sep=""))
