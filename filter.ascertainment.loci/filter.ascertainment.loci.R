################################################################################
##     Function to filter ascertainment-bias loci from a genlight object      ##
##                                                                            ##
##  Author: Diana A Robledo-Ruiz. Research Fellow, Monash University          ##
##  Date: 2022-08-18                                                          ##
##                                                                            ##
##  This function requires:                                                   ##
##   - Input: a genlight object with min 2 pop in gl@other$ind.metrics$pop    ##
##   - User specified parameters:                                             ##
##       - seed = interger, 1 by default                                      ##
##       - n = maximum pop size, size from smallest pop by default            ##
##                                                                            ##
##  Output:                                                                   ##
##    $filtered.gl   - Genlight without ascertainment bias (loci removed)     ##
##    $asc.inds      - Vector with names of ascertainment individuals         ##
##    $removed.loci  - Vector with names of the removed loci                  ##
##    $results.table - Per pop sample size, no. polymorphic loci before       ##
##                     filtering, and no. polymorphic loci after filtering    ##
##                                                                            ##
##  Index:                                                                    ##
##    Line 26:  Function filter.ascertainment.bias                            ##
##    Line 202: Example of use for filter.highly.het                          ##
################################################################################


############################## Defining function ###############################
filter.ascertainment.bias <- function(gl, seed = 1, n = NULL){
  
  ############################ 1. Check pops sizes
  
  # Make sure there are no remnant levels in pop
  gl@other$ind.metrics$pop <- droplevels(gl@other$ind.metrics$pop)
  
  # Rename individual metrics dataframe
  df <- gl@other$ind.metrics
  
  # Check how many pops (at least 2)
  if(nlevels(df$pop) < 2){
    stop("There should be at least 2 populations in gl@other$ind.metrics.")
  }
  
  n.pop <- vector()
  for (i in 1:nlevels(df$pop)) {
    n.pop[i] <- nrow(df[df$pop == levels(df$pop)[i], ])
  }
  
  # Check if pop have different sizes
  if(min(n.pop) == max(n.pop)){
    stop("Population sizes are already equal, there is no ascertainment bias.")
  }
  
  # Make plot of population sizes
  barplot(n.pop,
          names.arg = levels(df$pop),
          las = 2,  # vertical names 
          xlab = "Populations",
          ylab = "Sample size",
          col = rainbow(nlevels(df$pop)))
  
  ########################### 2. Id ind to subsample
  
  # If n was not specified, make it min size
  if(is.null(n)){
    n <- min(n.pop)
  }
  
  message("Populations were subsampled to be max n = ", n)
  
  df$keep <- "drop"
  
  for (i in 1:nlevels(df$pop)) {
    
    # Subset by pop
    tmp <- df[df$pop == levels(df$pop)[i], ]
    
    if (nrow(tmp) <= n) {
      # If it's the required size or smaller pop, add "keep"
      df[df$pop == levels(df$pop)[i], "keep"] <- "keep"
      
    } else {
      # Otherwise randomly subsample rows by rownames
      set.seed(seed)
      x <- sample(row.names(tmp), n)
      df[x, "keep"] <- "keep"
    }
  }
  
  ##################### 3. Function to id monomorphic loci
  
  id.mono <- function(gen.df){    
    
    # Identify monomorphic loci number to be dropped
    mono.index <- vector()  # list of indices of loci to drop
    mono.names <- vector()  # list of names of loci to drop
    
    for (i in 1:ncol(gen.df)) {  # Number of column is the index
      locus <- gen.df[, i] # Vector for a locus
      if (all(locus == 0, na.rm = TRUE) |
          all(locus == 2, na.rm = TRUE) | all(is.na(locus))){
        mono.index <- c(mono.index, i)
        mono.names <- c(mono.names, colnames(gen.df)[i])
      }
    }
    return(list(index = mono.index, 
                names = mono.names))
  }
  
  asc.inds <- df[df$keep == "keep", "id"]
  
  ############### 4. Remove ind, id monomorphic loci, remove loci
  
  # Ind in rows, loci in columns
  gen <- as.data.frame(as.matrix(gl))
  
  # Keep in gen only inds (rows) to keep
  gen <- gen[rownames(df[df$keep == "keep", ]), ]
  
  # Identify monomorphic loci number to be dropped
  loci2drop <- id.mono(gen)
  
  message(length(loci2drop$index), 
          " loci out of ",
          ncol(gen), 
          " (",
          round((length(loci2drop$index)*100)/ncol(gen), 2),
          "%) will be removed from gl")
  
  # Create new genlight from which to drop loci
  gl2 <- gl[ , -loci2drop$index]
  
  # Update the loci metadata to drop loci there too
  gl2@other$loc.metrics <- gl@other$loc.metrics[-loci2drop$index, ]
  
  ############################ 5. Create BEFORE plot
  
  # Ind in rows, loci in columns
  gen <- as.data.frame(as.matrix(gl))
  
  before.poly.pop <- vector()
  
  # Calculate polymorphic sites per pop
  for (i in 1:nlevels(df$pop)){
    # Subset per pop 
    subgen <- gen[rownames(df[df$pop == levels(df$pop)[i], ]), ]
    # Substract monomorphic loci to all loci
    before.poly.pop <- c(before.poly.pop, ncol(gen)-length(id.mono(subgen)$index))
  }
  
  # Make BEFORE plot
  barplot(before.poly.pop,
          names.arg = levels(df$pop),
          las = 2,  # vertical names 
          xlab = "Populations",
          ylab = "# polymorphic loci",
          main = "BEFORE filtering",
          col = rainbow(nlevels(df$pop)))
  
  ############################ 6. Create AFTER plot
  
  # Ind in rows, loci in columns
  gen <- as.data.frame(as.matrix(gl2))
  
  after.poly.pop <- vector()
  
  # Calculate polymorphic sites per pop
  for (i in 1:nlevels(df$pop)){
    # Subset per pop 
    subgen <- gen[rownames(df[df$pop == levels(df$pop)[i], ]), ]
    # Substract monomorphic loci to all loci
    after.poly.pop <- c(after.poly.pop, ncol(gen)-length(id.mono(subgen)$index))
  }
  
  # Make AFTER plot
  barplot(after.poly.pop,
          names.arg = levels(df$pop),
          las = 2,  # vertical names 
          xlab = "Populations",
          ylab = "# polymorphic loci",
          main = "AFTER filtering",
          col = rainbow(nlevels(df$pop)))
  
  
  ############################## 7. Create table 
  results.table <- as.data.frame(matrix(c(n.pop, before.poly.pop, after.poly.pop),
                                        nrow = 3,
                                        byrow = TRUE))
  
  colnames(results.table) <- levels(df$pop)
  rownames(results.table) <- c("Sample size", 
                               "# poly before", 
                               "# poly after")
  
  ############################## 8. Return output 
  return(list("filtered.gl"   = gl2,
              "asc.inds"      = asc.inds,
              "removed.loci"  = loci2drop$names,
              "results.table" = results.table))
}
################################################################################


################################ Example of use ################################
##   new.data <- filter.ascertainment.bias(gl = my.gl,                        ##
##                                       seed = 100,                          ##
##                                       n = 50)                              ##
##   new.data$filtered.gl                                                     ##
##   new.data$asc.inds                                                        ##
##   new.data$removed.loci                                                    ##
##   new.data$results.table                                                   ##
################################################################################
