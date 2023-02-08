library(ade4)
library(FactoMineR)
library(dendextend)
library(stats)
library(coRanking)
library(usethis)
library(StatMatch)
library(usedist)

#' Compute elbow
elbow <- function(data, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL) {
  
  if (is.null(xmin) && !is.null(ymin)) {
    stop("`xmin` cannot be NULL if `ymin` is not null.")
  }
  
  if (!is.null(xmin) && is.null(ymin)) {
    stop("`ymin` cannot be NULL if `xmin` is not null.")
  }
  
  if (is.null(xmax) && !is.null(ymax)) {
    stop("`xmax` cannot be NULL if `ymax` is not null.")
  }
  
  if (!is.null(xmax) && is.null(ymin)) {
    stop("`ymax` cannot be NULL if `xmax` is not null.")
  }
  
  data <- data[ , 1:2]
  
  if (!is.null(xmin)) {
    
    xdat <- data.frame(xmin, ymin)
    colnames(xdat) <- colnames(data)
    
    data <- rbind(xdat, data)
  }
  
  if (!is.null(xmax)) {
    
    xdat <- data.frame(xmax, ymax)
    colnames(xdat) <- colnames(data)
    
    data <- rbind(xdat, data)
  }
  
  data <- data[order(data[ , 1]), ]
  
  constant <- data[c(1, nrow(data)), ]
  colnames(constant) <- c("x", "y")
  
  mod <- stats::lm(y ~ x, data = constant)
  
  for (i in 1:nrow(data)) {
    data[i, "constant"] <- round(mod$coef[[1]] + mod$coef[[2]] * data[i, 1], 3)
  }
  
  data[ , "benefits"] <- round(data[ , 2] - data[ , "constant"], 3)
  
  
  data[ , "SelectedorNot"] <- NA
  
  for (i in 1:nrow(data)) {
    if (data$Axe[i] <= data[which.max(data[ , "benefits"]), ]$"Axe") {
      data[i, "SelectedorNot"] <- "Selected"
    } else {
      data[i, "SelectedorNot"] <- "Not Selected"
    }
  }
  
  xxx <- data[-1, ]
  rownames(xxx) <- xxx[ ,1]
  
  return(xxx)
}





#' Compute the Number of Dimensions function
#@@@@trait_df = base de donn?es de traits centr?s-r?duits, uniquement ceux dans l'espace ph?notypique
#@@@@dim_pcoa = nombre d'axes dans l'ACP
#@@@@rep = nombre de r?plicats
#@@@@cores = parall?lisation sous linux, nombre de coeur mobilis? dans l'ordinateur, on ne peut en mettre qu'un sous windows
#@@@@Nb_genotype = permet de r?cup?rer le nombre de g?notype qui sera ins?rer en tant que nom de colonne et de ligne dans la matrice de distance sym?trique

##J'ai laiss? en #ce qu'il avait fait avant dans leur script initial


dimension_funct_taina <- function(trait_df, dim_pcoa = 10, rep = 99, cores = 1,metric_scaled = TRUE) { # , nb_genotype = 1:300
  # function modified by Taïna so that:
  # - PCA instead of PCOA
  # - 1 core instead of 3 by default
  # Previous script is commented when replaced
  # CAPITAL letters to indicate changes made
  
  # MODIFIED TO HAVE GENOTYPES NAMES
  
  ##Computing distance matrix with euclidian distance
  #sp_trdist <- calc_dist(trait_df, trait_category_df, colnames(trait_df), classical_gower)
  #sp_trdist <- ade4::quasieuclid(sp_trdist)

  # MATRICE DE DISTANCE QU'ON POURRAIT EXTRAIRE
  sp_trdist <- dist(trait_df, method = "euclidean", upper=T, diag=T)
  # sp_trdist <- dist_setNames(sp_trdist, nb_genotype) # get column and row names
  
  
  # Computing PCoA-based functional spaces ----
  
  #pcoa_trdist <- ape::pcoa(sp_trdist)
  pcoa_trdist <- PCA(X=trait_df, ncp=dim_pcoa, graph=F) # PCA instead of PCoA
  
  
  # Number of dimensions to keep given the input from user and number of PC
  # with positive eigenvalues
  
  #nbdim <- min(dim_pcoa, ncol(pcoa_trdist$"vectors"))
  nbdim <- min(dim_pcoa, ncol(pcoa_trdist$ind$coord))
  
  # Keeping species coordinates on the 'nbdim' axes and renaming axes
  
  #sp_coord <- pcoa_trdist$"vectors"[ , 1:nbdim]
  sp_coord <- pcoa_trdist$ind$coord[ , 1:nbdim]
  colnames(sp_coord) <- paste0("PC", 1:nbdim)
  
  # Computing quality of multidimensional spaces: storing trait-based distance
  # (=input) in a 3-variables dataframe with 1 row for each pair of species
  # (names stored in the 2 first columns)
  
  distsp_df <- dendextend::dist_long(sp_trdist)
  names(distsp_df) <- c("sp.x", "sp.y", "distsp_df")
  
  fspaces_nm <- vector()
  
  # Increase Nnmber of PCoA dimensions
  
  for (k in 1:nbdim) {
    
    # Computing Euclidean distance between species
    dist_sp_k <- dist(sp_coord[ , 1:k])
    
    # Storing these distances as additional column to previous data frame
    value <- dendextend::dist_long(dist_sp_k)$"distance"
    distsp_df[ , paste0("pcoa_", k, "dim")] <- value
    
    # Adding name of funct space
    fspaces_nm <- c(fspaces_nm, paste0("pcoa_", k, "dim"))
  }
  
  
  # Compute coranking index ----
  
  computeAUCandNullmod <- lapply(1:nbdim, function(i) {
    
    #pcoa_axes <- pcoa_trdist$"vectors"[ , 1:i]
    pcoa_axes <- pcoa_trdist$ind$coord[ , 1:i] # PCA NAMES
    
    D_dimen   <- stats::dist(pcoa_axes)
    
    Co_rank   <- coRanking::coranking(sp_trdist, D_dimen, input_Xi = "dist")
    NX        <- coRanking::R_NX(Co_rank)
    AUC       <- coRanking::AUC_ln_K(NX)
    
    dat <- data.frame(Axe = i, AUC = round(AUC, 3))
    
    if (i == 1) rand_table <- NA
    
    if (i != 1) {
      
      usethis::ui_done(
        paste("Compute NULL model for AUC, nbdim =", i)
      )
      
      if (i == nbdim) usethis::ui_done(paste0("Last one!"))
      
      rand_table <- do.call(rbind, pbmcapply::pbmclapply(1:rep, function (j) {
        
        real_axes <- pcoa_axes[ , -i]
        rand_axis <- sample(pcoa_axes[ , i])
        new_axes  <- data.frame(real_axes, rand_axis)
        
        D_dimen_rand <- dist(new_axes)
        
        Co_rank_rand <- coRanking::coranking(sp_trdist, D_dimen_rand,
                                             input_Xi = "dist")
        NX_rand      <- coRanking::R_NX(Co_rank_rand)
        AUCrand      <- coRanking::AUC_ln_K(NX_rand)
        
        return(AUCrand)
        
      }, mc.cores = cores))
    }
    
    res <- list(dat, rand_table)
  })
  
  
  AUC_SES_pval <- do.call(cbind, pbmcapply::pbmclapply(1:nbdim, function(k) {
    
    AUC_obs <- computeAUCandNullmod[[k]][[1]][ , 2]
    
    if (k == 1) {
      
      SES  <- NA
      pval <- NA
      
    } else {
      
      SES <- (computeAUCandNullmod[[k]][[1]][ , 2] -
                mean(computeAUCandNullmod[[k]][[2]])) /
        sd(computeAUCandNullmod[[k]][[2]])
      
      # p-value = proportion of null value inferior to obs Beta (+1)
      pval <- length(
        which(computeAUCandNullmod[[k]][[2]] >
                computeAUCandNullmod[[k]][[1]]$"AUC")) /
        (length(computeAUCandNullmod[[k]][[2]]) + 1)
    }
    
    res_final <- rbind(AUC_obs, SES, pval)
    return(res_final)
    
  }, mc.cores = cores))
  


  # Detect inflexion point ----
  
  axes_table <- data.frame(Axe = 1:nbdim, AUC = AUC_SES_pval[1, ])
  ebow_meth  <- elbow(axes_table, xmin = 0, ymin = 0)
  
  
  # Compute the choosen quality metric for each functional space ----
  
  if (!metric_scaled) {
    
    # Compute deviance between distance in each functional space and
    # trait-based distance
    
    dev_distsp <- data.frame(
      distsp_df[ , c("sp.x", "sp.y")],
      distsp_df[ , fspaces_nm] - distsp_df[ , "distsp_df"]
    )
    
    
    # Compute absolute deviance ----
    
    abs_dev_distsp <- data.frame(
      dev_distsp[ , c("sp.x", "sp.y")],
      abs(dev_distsp[, fspaces_nm])
    )
    
    # mean absolute deviation:
    
    qual_metric_mad <- apply(abs_dev_distsp[ , fspaces_nm], 2, mean)
    
  } else {
    
    # Compute deviance between distance in each functional space and
    # trait-based distance:
    
    scdistsp <- apply(distsp_df[ , fspaces_nm], 2, function(x) {
      x / max(x) * max(distsp_df[ , "distsp_df"])
    })
    
    
    # Compute deviance ----
    
    dev_scdistsp <- data.frame(
      distsp_df[ , c("sp.x", "sp.y")],
      scdistsp[, fspaces_nm] - distsp_df[, "distsp_df"]
    )
    
    # compute squared deviance:
    sqr_dev_scdistsp <- data.frame (dev_scdistsp[, c("sp.x", "sp.y")],
    (dev_scdistsp[, fspaces_nm])^2)
    
    #compute mean squared deviation:
    qual_metric_msd <- apply(sqr_dev_scdistsp[, fspaces_nm], 2, mean)
    
    
    # Compute absolute deviance ----
    
    abs_dev_scdistsp <- data.frame(
      dev_scdistsp[ , c("sp.x", "sp.y")], abs(dev_scdistsp[ , fspaces_nm])
    )
    
    
    # Compute mean absolute deviation ----
    qual_metric_mad <- apply(abs_dev_scdistsp[ , fspaces_nm], 2, mean)
  }
  
  
  # Return results ----
  
  res <- data.frame(
    dim                  = 1:nbdim,
    quality_fspaces_msd  = round(qual_metric_msd,3),
    MAD                  = round(qual_metric_mad,3),
    AUC                  = ebow_meth$"AUC",
    benefit_AUCebow_meth = ebow_meth$"benefits",
    SelectedbyAUCelbow   = ebow_meth$"SelectedorNot",
    SES                  = AUC_SES_pval[2, ],
    pvalue               = AUC_SES_pval[3, ]
  )
  
  return(res)
}

