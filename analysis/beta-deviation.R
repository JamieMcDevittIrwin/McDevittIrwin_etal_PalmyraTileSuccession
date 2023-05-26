# Dependencies for beta deviation code


# code from Vannette, R. L., and T. Fukami. 2017. Dispersal enhances beta diversity in nectar microbes. Ecology letters 20:901–910.


#### vanette & fukami.r #####
########################################
# code for the function beta_ses #######
########################################

beta.ses <-function(compo, null.matrices)
{
  ## Opens required packages
  require(vegetarian)
  require(vegan)
  
  ## Makes sure 'compo' is in the right format
  dim.compo <- length(dim(compo))
  if(dim.compo==0) 
    stop("'compo' does not have more than 1 dimension")
  if(dim.compo>2) 
    stop("This function does not know what to do when 'compo' has more than 2 dimensions")
  if(dim.compo==2) 
  {
    if(is.null(rownames(compo))) 
      rownames(compo) <- paste("site_", 1:nrow(compo), sep="")
    if(is.null(colnames(compo))) 
      colnames(compo) <- paste("sp_", 1:nrow(compo), sep="")
  } 
  
  ## Checks that values in the cells are all positive integers
  if(sum((compo - round(compo, 0)) != 0) > 0) 
    stop("This function cannot use species abundances that are not integers")
  if(sum(compo<0) > 0)    
    stop("This function cannot use negative species abundances")
  if(sum(is.na(compo)) > 0)    
    stop("This function cannot use NAs")
  
  
  fit <- vegdist(compo, method="bray")
  
  emp.bdisp.dist <- betadisper(fit, group=rep(1, times=nrow(compo)),type="centroid")
  emp.dist <- emp.bdisp.dist$distances
  
  ## Calculates density of individuals per site
  site.densities <- rowSums(compo)
  if(min(site.densities)==0) 
    warning("Some empirical sites seem to be empty of individuals")
  
  ## Calculates regional abundances per species - the regional SAD
  spp.abund <- colSums(compo)
  if(min(spp.abund)<=0) 
    warning("Some species have abundances of less than or equal to zero")
  
  ## Calculates species richness across all sites
  regional.richness <- length(which(spp.abund>0))
  
  ## Calculates species richness at each site
  site.richness <- rowSums(compo>0)
  mean.site.richness <- mean(site.richness, na.rm=TRUE)
  
  
  ## Calculates the number of null matrices to use
  null.N <- length(null.matrices)
  
  
  ## Creates empty objects to hold results from the null matrices
  # need a matrix length (empirical) x width(# of null communities)
  rand.distances <- as.data.frame(matrix(NA, nrow=nrow(compo), ncol=null.N))
  rownames(rand.distances) <- rownames(compo)
  colnames(rand.distances) <- paste("null_", 1:null.N, sep="")
  
  
  for (i in 1:null.N)
  {
    
    spp.abund <- colSums(compo)
    
    ## Defines the focal null composition matrix
    null.compo <- as.matrix(null.matrices[[i]])
    
    ## Calculates null regional abundances per species - the regional SAD
    null.spp.abund <- colSums(null.compo)
    if(identical(null.spp.abund, spp.abund)==FALSE) 
      warning("Simulated and observed species total abundances are not identical")
    
    ## Calculates species richness across all sites
    null.regional.richness <- length(which(null.spp.abund>0))
    
    ## Calculates species richness at each site
    null.site.richness <- rowSums(null.compo>0)
    null.mean.site.richness <- mean(null.site.richness)
    
    ## Calculates null composition distances among all possible pairs of sites
    
    null.dist.details <- vegdist(null.compo, method="bray")   
    rand.bdisper <- betadisper(null.dist.details, group=rep(1, times=nrow(null.compo)),type="centroid")
    rand.distances[,i] <- rand.bdisper$distances
    
  }
  
  
  # now use this to calculate the ses for each value 
  
  ES <- as.data.frame(matrix(NA, nrow=nrow(compo)))
  SES <- as.data.frame(matrix(NA, nrow=nrow(compo)), ncol=1)
  
  ES <- (emp.dist-rowMeans(rand.distances))
  
  for(i in 1:length(ES)){
    SES[i,1]<- ES[i]/sd(rand.distances[i,])
  }
  
  SES
}






##### Tello etal S2 code 2.r ######
### FUNCTION CODE ##################################################################################
assemblages.from.pool.randA <- function(compo, rand.N=999, fix.local.abund=TRUE, fix.rSAD=TRUE, 
                                        save.output=FALSE, save.format="matrices", path.to.save, show.progress=FALSE) 
{
  ## Makes sure there is a path to save files if needed
  if(save.output==TRUE & missing(path.to.save)) path.to.save <- getwd()
  
  ## Makes sure 'compo' is in the right format
  compo <- as.matrix(compo, ncol=ncol(compo))
  if(is.null(rownames(compo))) rownames(compo) <- paste("site_", 1:nrow(compoT), sep="")
  if(is.null(colnames(compo))) colnames(compo) <- paste("sp_", 1:nrow(compoT), sep="")
  
  ## Calculates density of individuals per site
  site.densities <- rowSums(compo)
  if(min(site.densities)<=0) warning("Some sites have no species")
  
  ## Calculates total regional abundance
  regional.abundance <- sum(site.densities)
  
  ## Calculates regional abundances per species - the regional SAD
  spp.abund <- colSums(compo)
  
  ## Checks for a couple of potential problems
  if(length(which((compo - round(compo, 0)) != 0)) > 0) 
    stop("This function can not randomize a matrix with species abundances that are not integers")
  if(min(spp.abund)<=0) 
    warning("Some species have abundances of less than 1")
  
  ## Finds the species with at least one individual
  spp.more.than.0.indices <- which(spp.abund>0)
  spp.more.than.0.names <- colnames(compo)[spp.more.than.0.indices]
  
  ## Produces a vector of species identities for each individual
  individual.id <- as.vector(rep(colnames(compo), spp.abund))
  
  ## Produces an empty list to save null composition matrices
  rand.datasets <- sapply(rep(NA, rand.N), list)
  names(rand.datasets) <- paste("RandDataset_", 1:rand.N, sep="")
  
  for(i in 1:rand.N)
  {
    if(show.progress==TRUE) 
      print(i)
    
    ## When the regional SAD is NOT fixed, produces a null SAD to be used in analyses by 
    ## randomly assigning individuals to species 
    if(fix.rSAD==FALSE)
    {
      # Assigns individuals to species at random
      null.sp.to.ind.assignation.1 <- spp.more.than.0.names # IMPORTANT - This ensures that all 
      # species in the region receive at least 
      # one individual in the null SAD
      null.sp.to.ind.assignation.2 <- character()
      
      if((regional.abundance-length(spp.more.than.0.names))>0) 
        null.sp.to.ind.assignation.2 <- sample(spp.more.than.0.names, 
                                               regional.abundance-length(spp.more.than.0.names), replace=TRUE)
      
      null.sp.to.ind.assignation <- c(null.sp.to.ind.assignation.1, null.sp.to.ind.assignation.2)
      
      # Calculates null abundances for species present with at least one individual
      pre.null.spp.abund <- as.numeric(table(null.sp.to.ind.assignation))
      pre.null.spp.abund <- pre.null.spp.abund[sample(1:length(pre.null.spp.abund))]
      
      # Inserts the null abundances into a vector that might contain some zeroes (i.e., if 
      # there are species present in the empirical composition table that do not have any 
      # individuals)
      null.spp.abund <- spp.abund
      null.spp.abund[spp.more.than.0.indices] <- pre.null.spp.abund
      
      # Produces a vector of null species identities for each individual
      null.individual.id <- as.vector(rep(colnames(compo), null.spp.abund))
    }
    
    ## When the regional SAD is fixed, makes the null SAD identical to the empirical SAD  
    if(fix.rSAD==TRUE)
    {
      null.spp.abund <- spp.abund
      null.individual.id <- individual.id
    }
    
    ## Assigns individuals to local sites at random
    if(fix.local.abund==TRUE) 
      null.site.assignation <- sample(rep(rownames(compo), times=site.densities), 
                                      regional.abundance, replace=FALSE)
    
    if(fix.local.abund==FALSE) 
      null.site.assignation <- sample(rownames(compo), regional.abundance, replace=TRUE)
    
    ## Creates a null species composition matrix using the null assignation of individuals 
    ## to sites
    null.compo <- tapply(rep(1, regional.abundance), list(null.site.assignation, 
                                                          null.individual.id), sum)
    null.compo[is.na(null.compo)] <- 0
    
    ## Reshapes the null composition matrix in case there is only 1 species
    if(ncol(compo)==1) 
    {
      null.compo <- matrix(null.compo, ncol=1)
      colnames(null.compo) <- colnames(compo)
      rownames(null.compo) <- unique(null.site.assignation)
    }
    
    ## Inserts empty plots into the null composition matrix in case there are any
    if(nrow(null.compo)<nrow(compo))
    {
      missing.plots <- setdiff(rownames(compo), rownames(null.compo))
      
      empty.plots <- matrix(0, ncol=ncol(null.compo), nrow=length(missing.plots))
      colnames(empty.plots) <- colnames(null.compo)
      rownames(empty.plots) <- missing.plots
      
      null.compo <- rbind(null.compo, empty.plots)
      
      if(ncol(compo)==1) 
      {
        null.compo <- matrix(null.compo, ncol=1)
        colnames(null.compo) <- colnames(compo)
        rownames(null.compo) <- c(unique(null.site.assignation), missing.plots)
      }
    }
    
    ## Matches the column and row names in the null and empirical composition matrices  
    if(ncol(compo)>1) null.compo <- null.compo[,match(colnames(compo), colnames(null.compo))]
    null.compo <- null.compo[match(rownames(compo), rownames(null.compo)),]
    
    ## Reshapes the null composition matrix in case there is only 1 species
    if(ncol(compo)==1)
    {
      null.compo <- matrix(null.compo, ncol=1)
      colnames(null.compo) <- colnames(compo)
      rownames(null.compo) <- rownames(compo)
    }
    
    ## Adds names to the species with zero abundances (i.e., empty columns)
    if(min(null.spp.abund)<=0)
    {
      colnames(null.compo)[null.spp.abund<=0] <- colnames(compo)[null.spp.abund<=0]
      null.compo[is.na(null.compo)] <- 0
    }
    
    ## Adds the null composition matrix to the list that compiles the results
    if(save.output==FALSE | save.format=="list") 
      rand.datasets[[i]] <- null.compo
    
    ## If requested, the null composition matrix is saved as a file
    if(save.output==TRUE & save.format=="matrices") 
      write.table(null.compo, file=paste(path.to.save, "\\", "RandDataset_", i, ".txt", sep=""), 
                  quote=FALSE, sep="\t", na="NA", dec=".", row.names=TRUE, col.names=TRUE)
  }
  
  ## If saving a list is requested, the full list of null composition matrices is saved as a file 
  if(save.format=="list")
    save(rand.datasets, file=paste(path.to.save, "\\", "RandDatasets", sep="")) 
  
  ## Makes a list of the parameters used in the randomization
  rand.parameters <- c(fix.local.abund, fix.rSAD, rand.N)
  names(rand.parameters) <- c("fix.local.abund", "fix.rSAD", "rand.N")
  
  ## Creates the output to return, depending on whether the main output was saved or not
  if(save.output==TRUE) 
    output <- list(rand.parameters, path.to.save)
  else 
    output <- list(rand.parameters, rand.datasets)
  
  names(output) <- c("rand.parameters", "rand.datasets")
  
  output  
}

