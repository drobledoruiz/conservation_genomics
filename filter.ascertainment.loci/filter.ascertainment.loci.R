################################################################################
##     Function to filter ascertainment-bias loci from a genlight object      ##
##                                                                            ##
##  Authors: Jesús Castrejón-Figueroa. R developer, Monash University         ##
##           Diana A Robledo-Ruiz. Research Fellow, Monash University         ##
##  Date: 2022-08-18                                                          ##
##                                                                            ##
##  This function requires:                                                   ##
##   - Input: a genlight object with min 2 pop in gl@other$ind.metrics$pop    ##
##   - User specified parameters:                                             ##
##       - seed = interger, 1 by default                                      ##
##       - n = maximum pop size, size from smallest pop by default            ##                                      ##
##                                                                            ##
##  Output:                                                                   ##
##    $filtered.gl   - Genlight object without ascertainment loci (removed)   ##
##    $removed.loci  - Vector with names of the removed (ascertainment) loci  ##
##                                                                            ##
##  Index:                                                                    ##
##    Line 24:  Function filter.ascertainment.loci                            ##
##    Line 183: Example of use for filter.highly.het                          ##
################################################################################


############################## Defining function ###############################
filter.ascertainment.loci <- function(gl, seed = 1, n = NULL){
  
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
          xlab = "Populations",
          ylab = "Population size",
          col = rainbow(nlevels(df$pop)))
  
  ########################### 2. Id ind to subsample
  
  # If n was not specified, make it min size
  if(is.null(n)){
    n <- min(n.pop)
  }
  
  message("Populations were subsampled to be max n = ", n)
  
  df$keep <- "drop"
  
  for (i in 1:length(levels(df$pop))) {
    
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
  
  # Make BEFORE plot
  barplot(after.poly.pop,
          names.arg = levels(df$pop),
          xlab = "Populations",
          ylab = "# polymorphic loci",
          main = "AFTER filtering",
          col = rainbow(nlevels(df$pop)))
  
  
  ############################## 7. Return output 
  return(list("filtered.gl" = gl2,
              "removed.loci" = loci2drop$names))
}
################################################################################


################################ Example of use ################################
## new.data <- filter.ascertainment.loci(gl = my.gl,                          ##
##                                       seed = 100,                          ##
##                                       n = 50)                              ##
## new.data$filtered.gl                                                       ##
## new.data$removed.loci                                                      ##
################################################################################