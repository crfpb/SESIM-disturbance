##############################################
#              output data frames            #
#    SESIM model - The American Naturalist   #
##############################################

# experimental design
experiment <- data.frame(treatment=rep(1:(length(rate)*length(eP)*length(eH)*length(disturb))),
                         dispersal=NA, environment=NA, heterog=NA, disturbance=NA)
# sampled abundances
X_save <- array(data=NA, dim=c(numCom, nSp, length(sampleV)))
# diversity results
results <- data.frame(treatment=rep(NA, dim(experiment)[1]*numCom*reps),
                      dispersal=NA, environment=NA, heterog=NA, disturbance=NA,
                      replicate=NA, patch=rep(1:numCom),
                      alpha=NA, gamma=NA, beta=NA,
                      comrich=NA, metarich=NA, 
                      comeve=NA, metaeve=NA, 
                      comprod=NA, metaprod=NA)
# synchrony & variability results
sync_cv_results <- data.frame(dispersal=as.numeric(), environment=as.numeric(), 
                              heterog=as.numeric(), disturbance=as.numeric(),
                              replicates=as.numeric(), patch=as.numeric(), species=as.character(),
                              scale=as.character(), synchrony=as.numeric(), variability=as.numeric())
