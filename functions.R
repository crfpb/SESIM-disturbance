##############################################
#                 functions                  #
#    SESIM model - The American Naturalist   #
##############################################

# loop id function
ijloop <- function(idloop=loop_id, Comnum=numCom){
  ((idloop*Comnum)-(Comnum-1)):(idloop*Comnum)
}

# beta diversity function
# reference: Ricotta 2017
betaeve <- function(dataset){
  library(vegan)
  n_col <- ncol(dataset)
  n_row <- nrow(dataset)
  total <- apply(dataset, 1, sum)
  species_weights <- total/sum(dataset)
  rel_abu <- sweep(dataset, 1, total, "/")
  rel_abu[is.na(rel_abu)] <- 0
  h_shannon <- vegan::diversity(rel_abu, index = "shannon", MARGIN = 1)
  pielou <- h_shannon/log(n_col)
  pielou_com <- 1-sum(pielou/n_row)
  weighted_pielou <- 1-sum(pielou*species_weights)
  output <- list("Average Beta"=pielou_com,"Species weights"=species_weights, "Weighted Beta" = weighted_pielou)
  return(output)
}

# synchrony function
# reference: function "synch_onerep" modified from codyn package
synchrony <- function(df, time.var, species.var, abundance.var, metric = "Gross"){
  metric = match.arg(metric, choices = c("Gross")) # removed "Loreau" from choices
  # remove any species that were never present
  df <- subset(df, abundance.var > 0)
  # fill in 0
  spplist <- unique(df[species.var])
  # added dimension of species lists #
  sppsize <- dim(spplist)[1]
  # change ends here #
  yearlist <- unique(df[time.var])
  fulllist <- expand.grid(species.var = spplist[[species.var]],
                          time.var = yearlist[[time.var]])
  # recapture original names
  names(fulllist) = c(species.var, time.var)
  df2 <- merge(df[c(species.var, time.var, abundance.var)], fulllist, all.y = T)
  df2[is.na(df2)] <- 0
  # removed if statement for metric="Loreau" #
  corout <- as.data.frame(cbind(species.var = as.character(),
                                "sppcor" =  as.numeric()))
  # check to see if there are species which do not vary within a subplot
  nonvary <- apply(transpose_community(df, time.var, species.var, abundance.var), 2, sd)
  # added NA to 0 #
  #nonvary[is.na(nonvary)] <- 0
  # change ends here #
  if(any(nonvary == 0 | is.na(nonvary))) {
    warning("One or more species has non-varying abundance within a subplot and has been omitted")
    # remove non-varying species from "spplist" for this subplot
    spplist <- data.frame(species =
                            spplist[[1]][is.na(match(as.character(unlist(spplist)), names(nonvary)[which(nonvary == 0)]))]
    )
  }
  # added #
  # return NA
  if(dim(spplist)[1]==0){
    corout <- as.data.frame(cbind(1:sppsize, NA))
  } else {
    # change ends here #  
    for (i in 1:nrow(spplist)){
      myspp <- as.character(spplist[[1]][i])
      focalspp <- df2[which(df2[species.var] == myspp),]
      com.focalspp <- transpose_community(focalspp, time.var, species.var, abundance.var)
      otherspp <- df2[which(df2[species.var] != myspp),]
      com.otherspp <- transpose_community(otherspp, time.var, species.var, abundance.var)
      agg.otherspp <- rowSums(com.otherspp)
      sppcor <- stats::cor(agg.otherspp, com.focalspp)
      subout <- as.data.frame(cbind(myspp, sppcor))
      names(subout) = c(species.var, "sppcor")
      subout$sppcor <- as.numeric(as.character(subout$sppcor))
      corout <- rbind(corout, subout)
    }
  }
  # average correlation for the community
  synchrony <- mean(corout$sppcor)
  return(synchrony)
}
# transpose community
transpose_community <- function(df, time.var, species.var, abundance.var){
  df <- as.data.frame(df)
  # remove unused levels if species is a factor
  df[species.var] <- if(is.factor(df[[species.var]]) == TRUE){factor(df[[species.var]])} else {df[species.var]}
  # sort by time and species
  df <- df[order(df[[time.var]], df[[species.var]]),]
  # cast as a species x time data frame; insert NA for 0
  comdat <- tapply(df[[abundance.var]], list(df[[time.var]], as.vector(df[[species.var]])), sum)
  comdat[is.na(comdat)] <- 0
  comdat <- as.data.frame(comdat)
  # results
  return(comdat)
}

# variability function
# reference: Wang et al. 2019
var.partition <- function(metacomm_tsdata){
  ts_metacom <- apply(metacomm_tsdata,2,sum)
  ts_patch <-  apply(metacomm_tsdata,c(2,3),sum)
  ts_species <- apply(metacomm_tsdata,c(1,2),sum)
  sd_metacom <- sd(ts_metacom)
  sd_patch_k <- apply(ts_patch,2,sd)
  sd_species_i <- apply(ts_species,1,sd)
  sd_species_patch_ik <- apply(metacomm_tsdata,c(1,3),sd)
  mean_metacom <- mean(ts_metacom)
  CV_S_L <- sum(sd_species_patch_ik)/mean_metacom
  CV_C_L <- sum(sd_patch_k)/mean_metacom
  CV_S_R <- sum(sd_species_i)/mean_metacom
  CV_C_R <- sd_metacom/mean_metacom
  phi_S_L2R <- CV_S_R/CV_S_L
  phi_C_L2R <- CV_C_R/CV_C_L
  phi_S2C_L <- CV_C_L/CV_S_L
  phi_S2C_R <- CV_C_R/CV_S_R
  partition_3level <- c(CV_S_L=CV_S_L, CV_C_L=CV_C_L, CV_S_R=CV_S_R, CV_C_R=CV_C_R,
                        phi_S_L2R=phi_S_L2R, phi_C_L2R=phi_C_L2R, phi_S2C_L=phi_S2C_L,
                        phi_S2C_R=phi_S2C_R)
  return(partition_3level)
}

# weighted means
weighted_mean <- function(x, w, ..., na.rm = FALSE){
  if(na.rm){
    df_omit <- na.omit(data.frame(x, w))
    return(weighted.mean(df_omit$x, df_omit$w, ...))
  } 
  weighted.mean(x, w, ...)
}

# Gross' synchrony & variability partitioning function 
syn_CV_partitioning <- function(X_saved){
  # data frame for all results
  synchrony_data <- data.frame(dispersal=as.numeric(), environment=as.numeric(), 
                               heterog=as.numeric(), disturbance=as.numeric(),
                               replicates=as.numeric(), patch=as.numeric(), species=as.character(),
                               scale=as.character(), synchrony=as.numeric(), variability=as.numeric())
  ### local species synchrony ###
  # note:
  # local species needs to be weighted
  # parameters
  patches <- numCom
  # data
  # all species
  all_data <- X_saved
  # data frame
  localspp_out <- data.frame(dispersal=as.numeric(), environment=as.numeric(), 
                             heterog=as.numeric(), disturbance=as.numeric(),
                             replicates=as.numeric(), patch=as.numeric(), species=as.character(),
                             scale=as.character(), synchrony=as.numeric(), variability=as.numeric())
  # loop
  for(p in 1:patches){
    # define patch
    pp <- as.numeric(p)
    # subset data
    d1 <- subset(all_data, patch==pp) 
    data <- subset(d1, select=c(time, replicate, (9:16)))
    # define community
    community_all <- data.frame(time=rep(NA,(dim(data)[1]*nSp)), species=NA, abundances=NA, replicates=NA)
    community_all$time <- data$time
    community_all$replicates <- data$replicate                   
    community_all$species <- rep(c("V1","V2","V3","V4","V5","V6","V7","V8"), each=dim(data)[1])
    community_all$abundances <- c(data$V1,data$V2,data$V3,data$V4,data$V5,data$V6,data$V7,data$V8)
    # calculate synchrony
    Synchrony_all <- synchrony(community_all, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
    # output
    localspp <- data.frame(dispersal=a, environment=ePeriod,
                           heterog=eHeterog, disturbance=pdistr,
                           replicates=r, patch=pp, species=NA,
                           scale="localspp", synchrony=Synchrony_all, variability=NA)
    localspp_out <- rbind(localspp_out, localspp)
  }
  # calculate weights for synchrony
  localspp_all_weights_S <- as.vector(unlist(tapply(localspp_out[,9], localspp_out[,6], function(x) (abs(x)/sum(abs(x))))))
  # weight data for synchrony
  localspp_all_weighted_S <- weighted_mean(localspp_out[,9], localspp_all_weights_S, na.rm=T)
  # final output for local species results
  localspp_data <- data.frame(dispersal=a, environment=ePeriod, 
                              heterog=eHeterog, disturbance=pdistr,
                              replicates=r, patch=NA, species=NA,
                              scale="localspp", synchrony=localspp_all_weighted_S, variability=NA)
  # combine to return object
  synchrony_data <- rbind(synchrony_data, localspp_data)
  
  
  ### regional species synchrony ###
  # create metacommunity abundance data
  reg_list <- split(X_saved[,9:16], X_saved$time)
  reg_sum <- as.data.frame(do.call(rbind,(lapply(names(reg_list), function(x) apply(reg_list[[x]], 2, sum)))))
  treats <- data.frame(str_split_fixed(names(reg_list), "[.]", 1))
  reg_abund <- cbind(treats, reg_sum)
  colnames(reg_abund)[1] <- "time"
  # data
  data_reg <- reg_abund
  # "loop"
  # define community
  community_all <- data.frame(time=rep(NA,dim(data_reg)[1]*nSp), species=NA, abundances=NA, replicates=NA)
  community_all$time <- as.numeric(data_reg$time)
  community_all$replicates <- unique(X_saved$replicate)               
  community_all$species    <- rep(c("V1","V2","V3","V4","V5","V6","V7","V8"), each=dim(data_reg)[1])
  community_all$abundances <- c(data_reg$V1,data_reg$V2,data_reg$V3,data_reg$V4,data_reg$V5,data_reg$V6,data_reg$V7,data_reg$V8)
  # calculate synchrony
  Synchrony_all <- synchrony(community_all, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  # outputs
  out_all <- data.frame(dispersal=a, environment=ePeriod, 
                        heterog=eHeterog, disturbance=pdistr,
                        replicates=r, patch=NA, species=NA,
                        scale="regionalspp", synchrony=Synchrony_all, variability=NA)
  # combine to return object
  synchrony_data <- rbind(synchrony_data, out_all)
  
  
  ### species spatial synchrony ###
  # note:
  # species spatial needs to be weighted
  # setting parameters
  species <- c("V1","V2","V3","V4","V5","V6","V7","V8")
  # dataframe
  sppspatial_out <- data.frame(dispersal=as.numeric(), environment=as.numeric(), 
                               heterog=as.numeric(), disturbance=as.numeric(),
                               replicates=as.numeric(), patch=as.numeric(), species=as.character(),
                               scale=as.character(), synchrony=as.numeric(), variability=as.numeric())
  # loop
  for(spp in 9:16){
    # community dataframe
    community <- data.frame(time=rep(NA,(dim(data)[1]*5)), species=NA, abundances=NA, replicates=NA)
    # subset data
    data <- subset(all_data, select=c(time, patch, replicate, spp))
    # define community
    community$time       <- data$time
    community$species    <- data$patch
    community$abundances <- data[,4]
    community$replicates <- data$replicate                   
    # calculate synchrony
    Synchrony_spp <- synchrony(community, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
    # output
    out_spp <- data.frame(dispersal=a, environment=ePeriod, 
                          heterog=eHeterog, disturbance=pdistr,
                          replicates=r, patch=NA, species=spp, 
                          scale="sppspatial", synchrony=Synchrony_spp, variability=NA)
    # combine outputs
    sppspatial_out <- rbind(sppspatial_out, out_spp)
  }
  # calculate weights for synchrony 4=species
  sppspatial_all_weights_S <- as.vector(unlist(tapply(sppspatial_out[,9], sppspatial_out[,7], function(x) (abs(x)/sum(abs(x))))))
  # weight data for synchrony
  sppspatial_all_weighted_S <- weighted_mean(sppspatial_out[,9], sppspatial_all_weights_S, na.rm=T)
  # final output for local species results
  sppspatial_data <- data.frame(dispersal=a, environment=ePeriod, 
                                heterog=eHeterog, disturbance=pdistr,
                                replicates=r, patch=NA, species=NA,
                                scale="sppspatial", synchrony=sppspatial_all_weighted_S, variability=NA)
  # combine to return object
  synchrony_data <- rbind(synchrony_data, sppspatial_data)
  
  
  ### community spatial synchrony ###
  # data
  all_data$comm <- apply(X_saved[,9:16], 1, sum)
  # community data frame
  community_all <- data.frame(time=rep(NA,(dim(all_data)[1])), species=NA, abundances=NA, replicates=NA)
  # subset data
  data_all <- subset(all_data, select=c(time, patch, replicate, comm))
  # define community
  community_all$time <- data_all$time
  community_all$species <- data_all$patch
  community_all$replicates <- data_all$replicate
  community_all$abundances <- data_all[,4]
  # calculate synchrony
  Synchrony_all <- synchrony(community_all, time.var="time", species.var="species", abundance.var="abundances", metric="Gross")
  # outputs
  out_all <- data.frame(dispersal=a, environment=ePeriod, 
                        heterog=eHeterog, disturbance=pdistr,
                        replicates=r, patch=NA, species=NA,
                        scale="commspatial", synchrony=Synchrony_all, variability=NA)
  # combine to return object
  synchrony_data <- rbind(synchrony_data, out_all)
  
  ###
  
  # variability
  
  # transform abundance into biomass
  X_biom <- X_saved
  
  # create empty array of dimensions
  # row=spp, column=time(10), array=patch(5)
  empty_array_all <- array(NA, dim=c(nSp, length(unique(X_biom$time)), 5))
  
  # populate array with community data
  for(pp in 1:5){
    comm_all <- t(subset(X_biom, X_biom$patch==pp)[,9:16])
    empty_array_all[,,pp] <- comm_all
  }
  
  # calculate synchrony and CV partitioning (Wang et al. 2019)
  part_calc_all <- as.data.frame(var.partition(empty_array_all))
  
  # organize output
  out_CV_all <- data.frame(dispersal=a, environment=ePeriod, 
                           heterog=eHeterog, disturbance=pdistr,
                           replicates=r, patch=NA, species=NA,
                           scale=c("Population","Community","Metapopulation","Metacommunity"), 
                           synchrony=NA, 
                           variability=c(part_calc_all[1,1], part_calc_all[2,1],
                                         part_calc_all[3,1], part_calc_all[4,1]))
  
  # combine all results 
  synchrony_data <- rbind(synchrony_data, out_CV_all)
  
  # return
  return(synchrony_data)
}




