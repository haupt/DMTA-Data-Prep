# Bios 301
# Fall, 2011
# Ryan Haupt

# This code was written as a final project for the course Bios301 at Vanderbilt University.

# Concept: Create a series of functions that can proccess raw DMTA data
# For original excel procedure, refer to Dr. Ungar's DMTA manual, Chapter 4

# DISCLAIMER: This code is not intended to be the most efficient system for organizing DMTA data,
# rather it is the method that most closely follows the steps outlined in the DTMA manual
# future versions of this software would do well to implement more succinct and efficent methods

# Clear the workspace
rm(list=ls())

# Set working directory to wherever results files are
setwd("~/DMTA-Data-Prep")

# Import results spreadsheets
# NOTE: SSFA software returns .csv files not properly formatted for the built-in read.csv function
# So double check after importing to ensure proper formatting before proceeding

# The below code implements using a sample data set provided in the repository and explained in the vignette

batch_results <- read.csv("../DMTA-Data-Prep/Aimaginotherium_batch1.csv")
volume_fill_results <- read.csv("../DMTA-Data-Prep/Aimaginotherium_batch1VolumeFillResults.csv")
Hasfc_results <- read.csv("../DMTA-Data-Prep/Aimaginotherium_batch1_AS_auto-split-stats.csv")

# Correct the units column
batch_results[[9]] <- "mcg"

Asfc_epLsar_Smc_prep <- function(batch_results){
  # Given a .csv file of results, this function returns the median values for Asfc, epLsar, & Smc_prep as a table

  n <- nrow(batch_results)
  n_results <- n/4
  ABCD <- c("A","B","C","D")
  place <- rep(ABCD, n/4)
  
  # Check for appropriate number of rows
  if (n%%4 != 0){
    cat("WARNING: Invalid number of scans may have been entered, function may not run to completion.")
  }
  
  #Create a vector of the filenames
  filename <- batch_results[1]
  
  #next we need to break the filename into it's parts
  filename <- unlist(filename)
  filename <- as.vector(filename)
  filename <- strsplit(filename, "_")
    
  species <- rep(0,n)
  i <- 1
  while (i <= n){        
    species[i] <- filename[[i]][1]
    i <- i+1
  }
  
  museumid <- rep(0,n)
  i <- 1
  while (i <= n){
    museumid[i] <- filename[[i]][2]
    i <- i+1
  }
  
  facet <- rep(0,n)
  i <- 1
  while (i <= n){
    facet[i] <- filename[[i]][3]
    i <- i+1
  }
  
  scan <- rep(0,n)
  i <- 1
  while(i <= n){
    scan[i] <-filename[[i]][4]
    # remove .sur extension
    # In the original manual these letters are left as lower-case, but R is case sensitive
    # And these values need to be compared to the Place column
    # Thus they are converted to upper case
    if (scan[i] == "a.sur"){
      scan[i] <- "A"
    } else if (scan[i] == "b.sur"){
      scan[i] <- "B"
    } else if (scan[i] == "c.sur"){
      scan[i] <- "C"
    } else {
      scan[i] <- "D"
    }
        
    i <- i+1
  }
  
  #recombine split vectors into a matrix
  splitfile <- cbind(species, museumid, facet, scan, batch_results[1], place, batch_results[2:12])
  
  # Need to make absolutely sure that the Scan column is the same as the Place column
  # Going through one by one in order to pinpoint location of errors
  i <- 1
  while (i <= n){
    if(splitfile[[4]][i] != splitfile[[6]][i]){
      return(cat("Error found in sample: ", splitfile[[2]][i], ". Please recheck original .csv file and run function again."))
    } 
    # The else statement needs to be vectorized to output properly
    # But it isn't necessary for the function to work
    #else {
      #cat("For sample", splitfile[[2]][i], ": Facet ", splitfile[[4]][i], " = Place ", splitfile[[6]][i], ".")
   # }
    i <- i+1
  }
  
  # Sort by place
  splitfile <- splitfile[order(place) , ]
  
  
  # Determine range of each scan
  a1 <- 1
  a2 <- n/4
  b1 <- a2 + 1
  b2 <- n/2
  c1 <- b2 + 1
  c2 <- n/4 * 3
  d1 <- c2 + 1
  d2 <- as.numeric(nrow(splitfile))
  
  # Create new data frames for each scan
  a_scans <- splitfile[a1:a2 ,]
  b_scans <- splitfile[b1:b2 ,]
  c_scans <- splitfile[c1:c2 ,]
  d_scans <- splitfile[d1:d2 ,]
  
  # Create empty vectrs to build dataframe of median results
  # (based on "_template for median calculations.xls" in DMTA Manual)
  taxon <- as.vector(a_scans[,1])
  specimen <- as.vector(a_scans[,2])
  facet_median <- as.vector(a_scans[,3])
  check_data <- rep(0,n_results) # Will aggregate the checks below and be displayed in the final data frame
  asfc_median <- rep(0, n_results)
  eplsar_median <- rep(0, n_results)
  smc_median <- rep(0, n_results)
  
  # need to check that samples in a_scans match b, c, & d scans respectively
  check1 <- rep(0, n_results)
  check2 <- rep(0, n_results)
  check3 <- rep(0, n_results)
  
  
  i <- 1
  while (i <= n_results){
    if(a_scans[i,2] != b_scans[i,2]){
      check1[i] <- 1
    }
    if(a_scans[i,2] != c_scans[i,2]){
      check2[i] <- 1
    }
    if(c_scans[i,2] != d_scans[i,2]){
      check3[i] <- 1
    }
    check_data[i] <- check1[i] + check2[i] + check3[i]
    i <- i+1
  }
  
  
  # Need to preserve original filenames to double check against filenames in next function
  filename_original_a <- as.vector(a_scans[,5])
  filename_original_b <- as.vector(b_scans[,5])
  filename_original_c <- as.vector(c_scans[,5])
  filename_original_d <- as.vector(d_scans[,5])
  
  # Create data frame to add result to
  medians <- cbind(taxon, specimen, facet_median, check_data, asfc_median, eplsar_median, smc_median, check1, check2, check3, filename_original_a, filename_original_b, filename_original_c, filename_original_d)
  
  # Calculate medians for complexity
  i <- 1
  while (i <= n_results){
    asfc_temp <- c(a_scans[i,7], b_scans[i,7], c_scans[i, 7], d_scans[i, 7])
    medians[i,5] <- median(asfc_temp)
    i <- i+1
  }
  
  # Calculate medians for anisotropy
  i <- 1
  while (i <= n_results){
    eplsar_temp <- c(a_scans[i,15], b_scans[i,15], c_scans[i, 15], d_scans[i, 15])
    medians[i,6] <- median(eplsar_temp)
    i <- i+1
  }
  
  # Calculate medians for scale of maximal complexity
  i <- 1
  while (i <= n_results){
    smc_temp <- c(a_scans[i,9], b_scans[i,9], c_scans[i, 9], d_scans[i, 9])
    medians[i,7] <- median(smc_temp)
    i <- i+1
  }
  
  return(medians)
}

medians <- Asfc_epLsar_Smc_prep(batch_results)

# Before running Tfv_prep, make sure the check_data column of medians does not contain any '1's as this would indicate mismatched data

Tfv_prep <- function(volume_fill_results, medians){
  # This function receives the results of the 'Asfc_epLsar_Smc_prep' function and the 'volume_fill_results'
  # And returns a single data frame containing compiled results
  
  #Need to remove original filenames from medians to use for comparisons later
  filename_original_a <- as.vector(medians[,11])
  filename_original_b <- as.vector(medians[,12])
  filename_original_c <- as.vector(medians[,13])
  filename_original_d <- as.vector(medians[,14])
  
  # Remake medians without original filenames
  medians <- medians[, 1:10]
  
  
  # Initialize n for fixing header row
  n <- nrow(volume_fill_results)
  
  # first need to recreate table with correct headers
  headers <- c("Species", "Filename", "Sfv(mcg^3)", "Ctfv(mcg^3)", "Ftfv(mcg^3)")
  data <- volume_fill_results[2:n, ]
  volume_fill_results <- data
  names(volume_fill_results) <- headers
  
  # resize n for new table
  n <- nrow(volume_fill_results)
  n_results <- n/4
  ABCD <- c("A","B","C","D")
  place <- rep(ABCD, n/4)
  
  # Check for appropriate number of rows
  if (n%%4 != 0){
    cat("WARNING: Invalid number of scans may have been entered, function may not run to completion.")
  }
  
  volume_sort <- cbind(volume_fill_results[1:2], place, volume_fill_results[3:5])
  
  # Need to diverge from the manual to create a "scans" column for comparison to the "place" column
  #Create a vector of the filenames
  filename <- volume_sort[2]
  filename <- unlist(filename)
  filename <- as.vector(filename)
  filename <- strsplit(filename, "_")
  
  scan <- rep(0,n)
  i <- 1
  while(i <= n){
    scan[i] <-filename[[i]][4]
    # remove .sur extension
    # In the original manual these letters are left as lower-case, but R is case sensitive
    # And these values need to be compared to the Place column
    if (scan[i] == "a.sur"){
      scan[i] <- "A"
    } else if (scan[i] == "b.sur"){
      scan[i] <- "B"
    } else if (scan[i] == "c.sur"){
      scan[i] <- "C"
    } else {
      scan[i] <- "D"
    }
        
    i <- i+1
  }
  
  volume_sort <- cbind(volume_sort, scan)
  
  # Need to make absolutely sure that the Scan column is the same as the Place column
  # Going through one by one in order to pinpoint location of errors
  i <- 1
  while (i <= n){
    if(volume_sort[[3]][i] != volume_sort[[7]][i]){
      return(cat("Errror found in sample: ", volume_sort[[2]][i], ". Please recheck original .csv file and run function again."))
    } 
    # Need to vectorize for proper output but the display is not critical for proper functioning
    #else {
    #  cat("For sample", volume_sort[[2]][i], ": Place ", volume_sort[[3]][i], " = Scan ", volume_sort[[7]][i], ".")
    #}
    i <- i+1
  }
  
  # sort data fram by scan (sorting by place should produce the same results)
  volume_sort <- volume_sort[order(scan) , ]
  
  # Determine range of each scan
  a1 <- 1
  a2 <- n/4
  b1 <- a2 + 1
  b2 <- n/2
  c1 <- b2 + 1
  c2 <- n/4 * 3
  d1 <- c2 + 1
  d2 <- as.numeric(nrow(volume_sort))
  
  # Create new data frames for each scan
  a_scans <- volume_sort[a1:a2 ,]
  b_scans <- volume_sort[b1:b2 ,]
  c_scans <- volume_sort[c1:c2 ,]
  d_scans <- volume_sort[d1:d2 ,]
  
  
  # Now to incorporate volume_sort into medians
  tfv <- rep(0, n_results)
  ftfv <- rep(0, n_results)
  
  #Add these columns to medians
  medians2 <- cbind(medians, tfv, ftfv)
  
  
  # Check volume filenames against originals
  i <- 1
  while (i <= n_results){
    if (a_scans[i,2] != filename_original_a[i]){
      return(cat("ERROR: Mismathced filenames, check original files then run function again."))
    }
    i <- i + 1
  }
  
  #Cacluate median for Tfv
  i <- 1
  temp_a <- as.vector(a_scans[,5])
  temp_b <- as.vector(b_scans[,5])
  temp_c <- as.vector(c_scans[,5])
  temp_d <- as.vector(d_scans[,5])
  while (i <= n_results){
    tfv_temp <- c(as.numeric(temp_a[i]), as.numeric(temp_b[i]), as.numeric(temp_c[i]), as.numeric(temp_d[i]))
    medians2[i,11] <- median(tfv_temp)
    i <- i+1
  }
  
  # Calculate median for FTfv
  i <- 1
  temp_a <- as.vector(a_scans[,6])
  temp_b <- as.vector(b_scans[,6])
  temp_c <- as.vector(c_scans[,6])
  temp_d <- as.vector(d_scans[,6])
  while (i <= n_results){
    ftfv_temp <- c(as.numeric(temp_a[i]), as.numeric(temp_b[i]), as.numeric(temp_c[i]), as.numeric(temp_d[i]))
    medians2[i,12] <- median(ftfv_temp)
    i <- i+1
  }
  
  # Output finalized data frame
  return(medians2)
}

medians <- Tfv_prep(volume_fill_results, medians)

HAsfc_prep <- function(Hasfc_results){
  # Initialize size variables
  n <- nrow(Hasfc_results)
  n_results <- n/4
  ABCD <- c("A","B","C","D")
  place <- rep(ABCD, n_results)
  
  # Check for appropriate number of rows
  if (n%%4 != 0){
    cat("WARNING: Invalid number of scans may have been entered, function may not run to completion.")
  }
  
  #Create a vector of the filenames
  filename <- Hasfc_results[1]
  
  #next we need to break the filename into it's parts
  filename <- unlist(filename)
  filename <- as.vector(filename)
  filename <- strsplit(filename, "_")
  
  species <- rep(0,n)
  i <- 1
  while (i <= n){        
    species[i] <- filename[[i]][1]
    i <- i+1
  }
  
  museumid <- rep(0,n)
  i <- 1
  while (i <= n){
    museumid[i] <- filename[[i]][2]
    i <- i+1
  }
  
  facet <- rep(0,n)
  i <- 1
  while (i <= n){
    facet[i] <- filename[[i]][3]
    i <- i+1
  }
  
  scan <- rep(0,n)
  i <- 1
  while(i <= n){
    scan[i] <-filename[[i]][4]
    # remove .sur extension
    # In the original manual these letters are left as lower-case, but R is case sensitive
    # And these values need to be compared to the Place column
    if (scan[i] == "a.sur"){
      scan[i] <- "A"
    } else if (scan[i] == "b.sur"){
      scan[i] <- "B"
    } else if (scan[i] == "c.sur"){
      scan[i] <- "C"
    } else {
      scan[i] <- "D"
    }
        
    i <- i+1
  }
  
  #recombine split vectors into a matrix
  splitfile <- cbind(species, museumid, facet, scan, Hasfc_results[1], place, Hasfc_results[2:89])
  
  # Need to make absolutely sure that the Scan column is the same as the Place column
  # Going through one by one in order to pinpoint location of errors
  i <- 1
  while (i <= n){
    if(splitfile[[4]][i] != splitfile[[6]][i]){
      return(cat("Errror found in sample: ", splitfile[[2]][i], ". Please recheck original .csv file and run function again."))
    } else {
      cat("For sample", splitfile[[2]][i], ": Facet ", splitfile[[4]][i], " = Place ", splitfile[[6]][i], ".")
    }
    i <- i+1
  }
  
  # Sort by place
  splitfile <- splitfile[order(place) , ]
  
  
  # Determine range of each scan
  a1 <- 1
  a2 <- n/4
  b1 <- a2 + 1
  b2 <- n/2
  c1 <- b2 + 1
  c2 <- n/4 * 3
  d1 <- c2 + 1
  d2 <- as.numeric(nrow(splitfile))
 
  # Separate splitfile into scans
  a_scans <- splitfile[a1:a2 ,]
  b_scans <- splitfile[b1:b2 ,]
  c_scans <- splitfile[c1:c2 ,]
  d_scans <- splitfile[d1:d2 ,]
  
  
  Hasfc <- function(scans){
    # Returns heterogeneity values based on formatted scan dataframes for a particular facet
    n <- as.numeric(nrow(scans)) # Reinitialize necessary sizing variable
   
    # Create blank matrix for results
    output <- matrix(0, nrow = n, ncol = 10)    
    
    # Counting variables and column identifies
    i <- 1
    j <- 19
    k <- 8
    
    # Calculate medians for 'scans' and add results to 'output'
    while(i <= 150){
      raw_dev <- scans[,j]
      raw_asfc <- scans[,k]
      temp <- raw_dev/raw_asfc      
      
      output[i:(i+14)] <- temp
      
      j <- j+1
      k <- k+1
      i <- i+15
      
    }
    # Write table of outputs
    names <- c("2x2", "3x3", "4x4", "5x5", "6x6","7x7","8x8","9x9","10x10","11x11")
    hasfc_facet <- as.data.frame(output)
    names(hasfc_facet) <- names
    
    hasfc_facet <- cbind(scans[1:6], hasfc_facet)
    
    return(hasfc_facet)
  }
  
  Hasfc_a <- Hasfc(a_scans)
  Hasfc_b <- Hasfc(b_scans)
  Hasfc_c <- Hasfc(c_scans)
  Hasfc_d <- Hasfc(d_scans)
  
  
  #need to check that samples in a_scans match b, c, & d scans respectively
  check1 <- rep(0, n_results)
  check2 <- rep(0, n_results)
  check3 <- rep(0, n_results)
  check_data <- rep(0, n_results)
  
  i <- 1
  while (i <= n_results){
    if(Hasfc_a[i,2] != Hasfc_b[i,2]){
      check1[i] <- 1
    }
    if(Hasfc_a[i,2] != Hasfc_c[i,2]){
      check2[i] <- 1
    }
    if(Hasfc_c[i,2] != Hasfc_d[i,2]){
      check3[i] <- 1
    }
    check_data[i] <- check1[i] + check2[i] + check3[i]
    i <- i+1
  }
  
  #Cacluate median for HAsfc
  
  Hasfc_median <- function(y, n_results, A, B, C, D){
    # 'y' is the desired size of the Hasfc square (yxy)
    # n is the number of samples
    # A, B, C, D are four data.frames returned by the Hasfc function
    # Function returns a vector of median values at whatever yxy was specified
    i <- 1
    y <- y + 5
    
    vecA <- A[,y]
    vecB <- B[,y]
    vecC <- C[,y]
    vecD <- D[,y]
    
    median_vector <- rep(0, n_results)
    
    while(i <= n_results)
    {
      temp <- c(vecA[i], vecB[i], vecC[i], vecD[i])
      median_vector[i] <- median(temp)
      i <- i+1
    }
    
    return(median_vector)
  }
  
  #initialize the output dataframe
  hasfc_medians <- Hasfc_median(2, n_results,Hasfc_a, Hasfc_b, Hasfc_c, Hasfc_d)
  
  # Run the function enough times to get completed data.frame
  i <- 3
  while (i <= 11){
    temp_median <- Hasfc_median(i, n_results, Hasfc_a, Hasfc_b, Hasfc_c, Hasfc_d)
    hasfc_medians <- cbind(hasfc_medians, temp_median)
    
    i <- i+1
  }
  
  # Convert from matrix to data frame
  Hasfc_median <- as.data.frame(hasfc_medians)
  
  # Aggregate sample ids, check_data, and results
  Hasfc_median <- cbind(a_scans[1:3], check_data, Hasfc_median, check1, check2, check3 )
  
  # Add column labels
  names <- c("Taxon", "Spec", "Facet", "Check Data", "2x2", "3x3", "4x4", "5x5", "6x6","7x7","8x8","9x9","10x10","11x11")
  names(Hasfc_median) <- names
  
  return(Hasfc_median) 
}

Hasfc_medians <- HAsfc_prep(Hasfc_results)

# Set date for when you originally scanned the teeth
date <- "May-2011"

Master_prep <- function(medians, Hasfc_medians, date){
  # This function prepare a finalized spreadsheet based on the the outputs of the previous functions
  # and date is the date the samples were run
  # This function is NOT based on Ungar's method, but should be more genearlly applicable
  
  # First need to check that the samples from each input match
  n_results1 <- nrow(medians)
  n_results2 <- nrow(Hasfc_medians)
  
  if (n_results1 != n_results2){
    return(cat("ERROR: Sample size mismatch, function cannot continue."))
    }
  
  i <- 1
  while (i <= n_results1){
    if(medians[i,2] != Hasfc_medians[i,2]){
      return(cat("ERROR: Specimen mismatch, double check inputs and run function again."))
    }    
    i <- i+1
  }
  
  # Combine all available data into single data frame  
  final <- cbind(medians[,1:7], medians[,11:12], Hasfc_medians[,5:14], Hasfc_medians[,1:4], date)
  
  return(final)
}

DMTA_Finalized <- Master_prep(medians, Hasfc_medians, date)
