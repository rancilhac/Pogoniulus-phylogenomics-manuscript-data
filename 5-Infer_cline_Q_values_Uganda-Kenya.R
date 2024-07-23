library(hzar)
library("data.table")
library(gridExtra)

# Read in the data (might need to change path to input table)
data.U <- read.csv("ddRADseq_samples_metadata_Uganda-Kenya.csv", sep=",", header=T)
data.U <- na.omit(data.U)
colnames(data.U)[1:2] <- c("Sample", "Subspecies")

# Transform RFT distances to negative values
data.U$dist_UgandaKenya[data.U$Subspecies == "affinis"] <- -data.U$dist_UgandaKenya[data.U$Subspecies == "affinis"]

# Plot the data
plot(data.U$dist_UgandaKenya, data.U$Q_RFT, pch="+", xlab="distance", ylab="affinis ancestry", main="Uganda contact zone (individuals)", xlim=c(-700, 600))

# Set up structure to store Hzar runs and results
ind.analysis <- list()
ind.analysis$clines <- list()
ind.analysis$clines$data <- list()
ind.analysis$clines$models <- list()
ind.analysis$clines$fit <- list()
ind.analysis$clines$mcmc <- list()
ind.analysis$clines$results <- list()

ind.analysis$clines$data <- hzar.doMolecularData1DPops(data.U$dist_UgandaKenya, data.U$Q_RFT, rep(2, nrow(data.U)))
hzar.plot.obsData(ind.analysis$clines$data)

#Set up the four models to be tested
# model 1: observed pMin & pMax, no tail
ind.analysis$clines$models[["model1"]] <- hzar.makeCline1DFreq(ind.analysis$clines$data, scaling="fixed",tails="none")
# model 2: observed pMin & pMax, two tails
ind.analysis$clines$models[["model2"]] <- hzar.makeCline1DFreq(ind.analysis$clines$data, scaling="fixed",tails="both")
# model 3: observed pMin & pMax, left tail
ind.analysis$clines$models[["model3"]] <- hzar.makeCline1DFreq(ind.analysis$clines$data, scaling="fixed",tails="left")
# model 4: observed pMin & pMax, right tail
ind.analysis$clines$models[["model4"]] <- hzar.makeCline1DFreq(ind.analysis$clines$data, scaling="fixed",tails="right")

ind.analysis$clines$models <- sapply(ind.analysis$clines$models, hzar.model.addBoxReq,-620, 510,
                                     simplify = F)

ind.analysis$clines$fit$init <- sapply(ind.analysis$clines$models, hzar.first.fitRequest.old.ML,
                                       obsData=ind.analysis$clines$data,verbose=F, simplify=F)

# Set up MCMC and seeds (import for reproducibility). MCMC and burn-in lengths were chosen following visual inspection of chain convergence
seeds <- list(s1=c(54,34,65,38,45,89),
              s2=c(26,23,27,17,54,78),
              s3=c(395,1254,6985,7895,3652,7459),
              s4=c(98,36,59,36,789,12))
ind.analysis$clines$fit$init$model1$mcmcParam$chainLength <- 5e5
ind.analysis$clines$fit$init$model1$mcmcParam$burnin <- 5e4
ind.analysis$clines$fit$init$model1$mcmcParam$seed[[1]] <- seeds$s1
ind.analysis$clines$fit$init$model2$mcmcParam$chainLength <- 5e5
ind.analysis$clines$fit$init$model2$mcmcParam$burnin <- 5e4
ind.analysis$clines$fit$init$model2$mcmcParam$seed[[1]] <- seeds$s2
ind.analysis$clines$fit$init$model3$mcmcParam$chainLength <- 5e5
ind.analysis$clines$fit$init$model3$mcmcParam$burnin <- 5e4
ind.analysis$clines$fit$init$model3$mcmcParam$seed[[1]] <- seeds$s3
ind.analysis$clines$fit$init$model4$mcmcParam$chainLength <- 5e5
ind.analysis$clines$fit$init$model4$mcmcParam$burnin <- 5e4
ind.analysis$clines$fit$init$model4$mcmcParam$seed[[1]] <- seeds$s4

# Run optimizer for each model. Remove '#' in front of plot commands to vizualize results
ind.analysis$clines$mcmc$init <- list()
ind.analysis$clines$mcmc$init$model1 <- hzar.doFit(ind.analysis$clines$fit$init$model1)
#plot(hzar.mcmc.bindLL(ind.analysis$clines$mcmc$init$model1))
ind.analysis$clines$mcmc$init$model2 <- hzar.doFit(ind.analysis$clines$fit$init$model2)
#plot(hzar.mcmc.bindLL(ind.analysis$clines$mcmc$init$model2))
ind.analysis$clines$mcmc$init$model3 <- hzar.doFit(ind.analysis$clines$fit$init$model3)
#plot(hzar.mcmc.bindLL(ind.analysis$clines$mcmc$init$model3))
ind.analysis$clines$mcmc$init$model4 <- hzar.doFit(ind.analysis$clines$fit$init$model4)
#plot(hzar.mcmc.bindLL(ind.analysis$clines$mcmc$init$model4))

ind.analysis$clines$fit$chains <- lapply(ind.analysis$clines$mcmc$init, hzar.next.fitRequest)
ind.analysis$clines$fit$chains <- hzar.multiFitRequest(ind.analysis$clines$fit$chains, 
                                                       each=3, 
                                                       baseSeed=NULL)

# Set up random initial parameter values for all models
center.init <- runif(12,-30,600)
width.init <- runif(12,0,630)
deltaL.init <- runif(12,0,630)
tauL.init <- runif(12,0,1)
deltaR.init <- runif(12,0,630)
tauR.init <- runif(12,0,1)

#init model 1 (no tails)
for(i in 1:3){
  ind.analysis$clines$fit$chains[[i]]$modelParam$init["center"]= center.init[i]
  ind.analysis$clines$fit$chains[[i]]$modelParam$init["width"]= width.init[i]
}
  
# init model 2 (both tails)
for(i in 4:6){
  ind.analysis$clines$fit$chains[[i]]$modelParam$init["center"]= center.init[i]
  ind.analysis$clines$fit$chains[[i]]$modelParam$init["width"]= width.init[i]
  ind.analysis$clines$fit$chains[[i]]$modelParam$init["deltaL"]= deltaL.init[i]
  ind.analysis$clines$fit$chains[[i]]$modelParam$init["tauL"]= tauL.init[i]
  ind.analysis$clines$fit$chains[[i]]$modelParam$init["deltaR"]= deltaR.init[i]
  ind.analysis$clines$fit$chains[[i]]$modelParam$init["tauR"]= tauR.init[i]
}

# init model 3 (left tail)
for(i in 7:9){
  ind.analysis$clines$fit$chains[[i]]$modelParam$init["center"]= center.init[i]
  ind.analysis$clines$fit$chains[[i]]$modelParam$init["width"]= width.init[i]
  ind.analysis$clines$fit$chains[[i]]$modelParam$init["deltaL"]= deltaL.init[i]
  ind.analysis$clines$fit$chains[[i]]$modelParam$init["tauL"]= tauL.init[i]
}

# init model 4 (right tail)
for(i in 10:12){
  ind.analysis$clines$fit$chains[[i]]$modelParam$init["center"]= center.init[i]
  ind.analysis$clines$fit$chains[[i]]$modelParam$init["width"]= width.init[i]
  ind.analysis$clines$fit$chains[[i]]$modelParam$init["deltaR"]= deltaR.init[i]
  ind.analysis$clines$fit$chains[[i]]$modelParam$init["tauR"]= tauR.init[i]
}

# Write the initial parameter values to a files (might be useful to reproduce results)
init.matrix <- rbind(center=center.init, width=width.init, deltaL=c(rep("NA", 3), deltaL.init[4:9], rep("NA", 3)), tauL=c(rep("NA", 3), tauL.init[4:9], rep("NA", 3)),
                     deltaR=c(rep("NA", 3), deltaR.init[4:6], rep("NA", 3), deltaR.init[10:12]), tauR=c(rep("NA", 3), tauR.init[4:6], rep("NA", 3), tauR.init[10:12]))
write.table(init.matrix, "parameters_initialization_hzar.txt")

# Run the actual analyses (3 chain for each model)
ind.analysis$clines$mcmc$chains <- hzar.doChain.multi(ind.analysis$clines$fit$chains,
                                                      doPar=F,
                                                      inOrder=FALSE,
                                                      count=3)

# Summarize results
summary(do.call(mcmc.list,
                lapply(ind.analysis$clines$mcmc$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list,
                lapply(ind.analysis$clines$mcmc$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list,
                lapply(ind.analysis$clines$mcmc$chains[7:9],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list,
                lapply(ind.analysis$clines$mcmc$chains[10:12],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )


ind.analysis$clines$results$initDGs <- list(nullModel = hzar.dataGroup.null(ind.analysis$clines$data))
ind.analysis$clines$results$initDGs$model1 <- hzar.dataGroup.add(ind.analysis$clines$mcmc$init$model1)
ind.analysis$clines$results$initDGs$model2 <- hzar.dataGroup.add(ind.analysis$clines$mcmc$init$model2)
ind.analysis$clines$results$initDGs$model3 <- hzar.dataGroup.add(ind.analysis$clines$mcmc$init$model3)
ind.analysis$clines$results$initDGs$model4 <- hzar.dataGroup.add(ind.analysis$clines$mcmc$init$model4)

ind.analysis$clines$results$oDG <- hzar.make.obsDataGroup(ind.analysis$clines$results$initDGs)
ind.analysis$clines$results$oDG <- hzar.copyModelLabels(ind.analysis$clines$results$initDGs,
                                                        ind.analysis$clines$results$oDG)
ind.analysis$clines$results$oDG <- hzar.make.obsDataGroup(lapply(ind.analysis$clines$mcmc$chains,
                                                                 hzar.dataGroup.add),
                                                          ind.analysis$clines$results$oDG)

# Print results summary
print(summary(ind.analysis$clines$results$oDG$data.groups))

# Plot results
hzar.plot.cline(ind.analysis$clines$results$oDG)

# Write AICc for the four models + null model into a file
AICcTable <- hzar.AICc.hzar.obsDataGroup(ind.analysis$clines$results$oDG)
rownames(AICcTable) <- c("Null model", "Model1-no.tails", "Model2-both.tails", "Model3-left.tail", "Model4-right.tail")
write.table(AICcTable, "Hzar_AIC.txt", sep=" ", col.names = T, row.names = T, quote = F)

print(ind.analysis$clines$results$AICcTable <- hzar.AICc.hzar.obsDataGroup(ind.analysis$clines$results$oDG))

print(ind.analysis$clines$results$model.name <- rownames(ind.analysis$clines$results$AICcTable)[[ which.min(AICcTable$AICc)]])

ind.analysis$clines$results$model.selected <- ind.analysis$clines$results$oDG$data.groups[[ind.analysis$clines$results$model.name]]

print(hzar.getLLCutParam(ind.analysis$clines$results$model.selected ,
                         names(ind.analysis$clines$results$model.selected $data.param)))
print(hzar.get.ML.cline(ind.analysis$clines$results$model.selected))

# Get parameter estimates and 95% CI
width <- round(hzar.get.ML.cline(ind.analysis$clines$results$model.selected)$param.all$width, 1)
center <- round(hzar.get.ML.cline(ind.analysis$clines$results$model.selected)$param.all$center, 1)

fz.cline <- hzar.getCredParamRed(ind.analysis$clines$results$model.selected)
width.CI <- center.CI <- c()
for(i in 1:length(fz.cline$clines)){
  width.CI <- c(width.CI, (fz.cline$clines[[i]]$param.free$width))
  center.CI <- c(center.CI, (fz.cline$clines[[i]]$param.free$center))}
cat("center = ", center, " [ ", round(range(center.CI)[1], 1), ";", round(range(center.CI)[2], 1), " ]", sep="")
cat("width = ", width, " [ ", round(range(width.CI)[1],1), ";", round(range(width.CI)[2], 1), " ]", sep="")

# Plot ML cline for best model with 95% CI and print to a file
jpeg("HZAR_cline_best_model.jpg", width=1000, height=1000)
hzar.plot.fzCline(ind.analysis$clines$results$model.selected, pch = "+", main="Cline chrysoconus/affinis",
                  ylab="affinis ancestry", xlab="Distance to contact zone", xlim=c(-620,520), cex=2)

legend("topright", legend=c(rownames(AICcTable)[which.min(AICcTable$AICc)], 
                            paste("width = ", width, " [ ", round(range(width.CI)[1],1), ";", round(range(width.CI)[2], 1), " ]", sep=""), 
                            paste("center = ", center, " [ ", round(range(center.CI)[1], 1), ";", round(range(center.CI)[2], 1), " ]", sep="")), bty="n")
dev.off()

#Save all the variables to an R environment. Can be loaded with load("HZAR_U_env.RData")
save.image(file="HZAR_U_env.RData")


