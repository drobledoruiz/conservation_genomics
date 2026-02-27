################################################################################
##         Function to filter sex-linked SNPs from a genlight object          ##
##                                                                            ##
##  Author: Diana A. Robledo-Ruiz. Research Fellow, Monash University         ##
##  Date: 2023-03-08                                                          ##
##                                                                            ##
##  This function requires:                                                   ##
##   - Input: a genlight object with column 'sex' ('F' or 'M') in ind.metrics ##
##   - User specified parameters:                                             ##
##       - system = sex determination system ('zw' or 'xy').                  ##
##       - plots = produce output plots (default is 'TRUE').                  ##
##       - parallel = run in parallel (default is 'FALSE').                   ##
##                                                                            ##
##  Output:                                                                   ##
##   A list with 6 elements:                                                  ##
##    - $results.table -> Table with statistics (columns) for each loci (rows)##
##    - $w.linked      -> Genlight object with w-linked/y-linked loci         ##
##    - $sex.biased    -> Genlight object with sex-biased scoring rate loci   ##
##    - $z.linked      -> Genlight object with z-linked/x-linked loci         ##
##    - $gametolog     -> Genlight object with zw-gametolog/xy-gametolog loci ##
##    - $autosomal     -> Genlight object with autosomal loci                 ##
##                                                                            ##
## Index:                                                                     ##
##   Line 29: Function                                                        ##
##   Line 586: Example of use for filter.sex.linked                           ##
################################################################################


############################## Defining function ###############################
filter.sex.linked <- function(gl, system = NULL, plots = TRUE, parallel = FALSE) {

  start_time0 <- Sys.time()
  if(parallel){
    library(doParallel)
    ncores <- detectCores()-1
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
  }

  if(is.null(system)){
    stop("You must specify the sex-determination system with the parameter 'system' ('zw' or 'xy').")
  } else {
    if(!(system == 'zw' | system == 'xy')){
      stop("Parameter 'system' must be 'zw' or 'xy'.")
    }
  }
  
  # Transform genotypes to matrix and transpose
  gen <- as.data.frame(t(as.matrix(gl)))

  # Extract IDs per sex (FUNCTION IGNORES ALL UNSEXED INDS!)
  if(!("F" %in% gl@other$ind.metrics$sex | "M" %in% gl@other$ind.metrics$sex)){
    stop("Females and males in @other$ind.metrics$sex must be 'F' or 'M', respectively.")
  }
  
  ids.F <- gl@other$ind.metrics[gl@other$ind.metrics$sex == 'F' & !is.na(gl@other$ind.metrics$sex), 'id']
  ids.M <- gl@other$ind.metrics[gl@other$ind.metrics$sex == 'M' & !is.na(gl@other$ind.metrics$sex), 'id']

  message(paste("Detected ", length(ids.F), " females and ", length(ids.M), " males.", sep = ""))

  # Subset genotypes by sex
  gen.F <- gen[ , (colnames(gen) %in% ids.F)]
  gen.M <- gen[ , (colnames(gen) %in% ids.M)]



  ##################### 1. Sex-linked loci by scoring rate

  # Create a results table with loci as row names and index number per locus
  table <- data.frame(index = c(1:nrow(gen)),
                      row.names = row.names(gen))

  # Count missing (NA) and add as column to table
  table$count.F.miss <- rowSums(is.na(gen.F))
  table$count.M.miss <- rowSums(is.na(gen.M))

  # Count scored ("0", "1" or "2") and add as column to table
  table$count.F.scored <- rowSums(!is.na(gen.F))
  table$count.M.scored <- rowSums(!is.na(gen.M))

  if(parallel){
    message("Starting phase 1. Working in parallel...")
  } else {
    message("Starting phase 1. May take a while...")
  }
  
  # Apply Fisher's exact test (because there are observations with less than 5)
  if(parallel){
    xfisher <- foreach(i = 1:nrow(table), .combine = rbind)%dopar%{
      # Make vector of observed values
      obs <- matrix(c(table[i, "count.F.miss"],
                      table[i, "count.M.miss"],
                      table[i, "count.F.scored"],
                      table[i, "count.M.scored"]),
                    nrow = 2,
                    ncol = 2,
                    dimnames = list(c("F", "M"),
                                    c("miss", "scored")))

      # See if it's possible to use chisq test
      if (sum(obs) >= 1000) {

        # Convert zeros to 1 to not obtain an error with chisq-test
        obs[obs == 0] <- 1

        # Chisq-test
        chisq.res <- chisq.test(obs, correct = FALSE)

        # Add to results table
        return(data.frame(ratio = chisq.res$statistic, p.value = chisq.res$p.value))

      } else {
        # Convert zeros to 1 to not obtain an error with Fisher's test
        obs[obs == 0] <- 1

        # Run Fisher's exact test
        F.test <- fisher.test(obs)

        # Add to results table
        return(data.frame(ratio = F.test$estimate, p.value = F.test$p.value))
      }
    }
    table <- cbind(table, xfisher)
  } else {
    # Add empty columns for chisq-statistic and corresponding p-value
    table$ratio   <- NA
    table$p.value <- NA
    
    # Test for independece of sex and missingness
    for (i in 1:nrow(table)) {

      # Make matrix of observed values
      obs <- matrix(c(table[i, "count.F.miss"],
                      table[i, "count.M.miss"],
                      table[i, "count.F.scored"],
                      table[i, "count.M.scored"]),
                    nrow = 2,
                    ncol = 2,
                    dimnames = list(c("F", "M"),
                                    c("miss", "scored")))


      # See if it's possible to use chisq test
      if (sum(obs) >= 1000) {

        # Convert zeros to 1 to not obtain an error with chisq-test
        obs[obs == 0] <- 1

        # Chisq-test
        chisq.res <- chisq.test(obs, correct = FALSE)

        # Add to results table
        table[i, "ratio"]   <- chisq.res$statistic
        table[i, "p.value"] <- chisq.res$p.value

      } else {

        # Convert zeros to 1 to not obtain an error with Fisher's test
        obs[obs == 0] <- 1

        # Run Fisher's exact test
        F.test <- fisher.test(obs)

        # Add to results table
        table[i, "ratio"]   <- F.test$estimate
        table[i, "p.value"] <- F.test$p.value
      }
    }
  }
  

  # Adjust p-values for multiple comparisons (False discovery rate)
  table$p.adjusted <- p.adjust(table$p.value, method = "fdr")

  # Calculate scoring rate for females and males and add to results table
  table$scoringRate.F <- table$count.F.scored/(table$count.F.scored+table$count.F.miss)

  table$scoringRate.M <- table$count.M.scored/(table$count.M.scored+table$count.M.miss)

  ##### 1.1 W-linked or Y-linked loci
  # For zw sex-determination system
  if(system == "zw") {
    table$w.linked <- NA

    for (i in 1:nrow(table)) {
      if (table[i, "scoringRate.M"] <= 0.1 && table[i, "p.adjusted"] <= 0.01) {
        table[i, "w.linked"] <- TRUE
      } else {
        table[i, "w.linked"] <- FALSE
      }
    }
  }
  
  # For xy sex-determination system
  if(system == "xy") {
    table$y.linked <- NA

    for (i in 1:nrow(table)) {
      if (table[i, "scoringRate.F"] <= 0.1 && table[i, "p.adjusted"] <= 0.01) {
        table[i, "y.linked"] <- TRUE
      } else {
        table[i, "y.linked"] <- FALSE
      }
    }
  }

  
  ##### 1.2 Loci with sex-biased scoring rate
  table$sex.biased <- NA

  for (i in 1:nrow(table)) {
    if (table[i, "p.adjusted"] <= 0.01 && table[i, 11] == FALSE) {
      table[i, "sex.biased"] <- TRUE
    } else {
      table[i, "sex.biased"] <- FALSE
    }
  }
  
    ##### 1.3 Plot BEFORE vs AFTER
  if(plots) {
    message("Building call rate plots.")
    
    # For zw sex-determination system
    if(system == "zw") {
      BEF.mis <- plot(x = table$scoringRate.F,
                      y = table$scoringRate.M,
                      main = "BEFORE",
                      xlab = "Call rate Females",
                      ylab = "Call rate Males",
                      xlim = c(0, 1),
                      ylim = c(0, 1))
  
      AFT.mis <- plot(x = table[table$w.linked == FALSE & table$sex.biased == FALSE,
                                "scoringRate.F"],
                      y = table[table$w.linked == FALSE & table$sex.biased == FALSE,
                                "scoringRate.M"],
                      main = "AFTER",
                      xlab = "Call rate Females",
                      ylab = "Call rate Males",
                      xlim = c(0, 1),
                      ylim = c(0, 1))
    }
  
    # For xy sex-determination system
    if(system == "xy") {
      BEF.mis <- plot(x = table$scoringRate.F,
                      y = table$scoringRate.M,
                      main = "BEFORE",
                      xlab = "Call rate Females",
                      ylab = "Call rate Males",
                      xlim = c(0, 1),
                      ylim = c(0, 1))
  
      AFT.mis <- plot(x = table[table$y.linked == FALSE & table$sex.biased == FALSE,
                                "scoringRate.F"],
                      y = table[table$y.linked == FALSE & table$sex.biased == FALSE,
                                "scoringRate.M"],
                      main = "AFTER",
                      xlab = "Call rate Females",
                      ylab = "Call rate Males",
                      xlim = c(0, 1),
                      ylim = c(0, 1))
    }
  }


  #################### 2. Sex-linked loci by heterozygosity
  # Count heterozygotes ("1") and add as column to results table
  table$count.F.het <- rowSums(gen.F == 1, na.rm = TRUE)
  table$count.M.het <- rowSums(gen.M == 1, na.rm = TRUE)

  # Count homozygotes ("0" or "2") and add as column to results table
  table$count.F.hom <- rowSums(gen.F != 1,  # Ignores NAs
                               na.rm = TRUE)
  table$count.M.hom <- rowSums(gen.M != 1,  # Ignores NAs
                               na.rm = TRUE)

  # Add empty columns for statistic and corresponding p-value
  #  table$stat         <- NA
  #  table$stat.p.value <- NA


  message("Done. Starting phase 2.")
  
  if(parallel){
    xstat <- foreach(i = 1:nrow(table), .combine = rbind)%dopar%{

      # Exclude w.y-linked loci and loci with sex-biased score
      if (table[i, 11] == TRUE | table[i, "sex.biased"] == TRUE) {
        stat.value <- NA
        stat.p.value <- NA

      } else {

        # Make contingency table
        contingency <- matrix(c(table[i, "count.F.het"],
                                table[i, "count.M.het"],
                                table[i, "count.F.hom"],
                                table[i, "count.M.hom"]),
                              nrow = 2,
                              ncol = 2,
                              dimnames = list(c("F", "M"), c("het", "hom")))

        # Check if Yate's correction is necessary (n =< 20, Sokhal & Rohlf 1995)
        if (sum(contingency) >= 1000) {

          # Convert zeros to 1 to not obtain an error with chisq-test
          contingency[contingency == 0] <- 1

          # Chisq-test
          chisq.res <- chisq.test(contingency, correct = FALSE)

          # Add to results table
          stat.value   <- chisq.res$statistic
          stat.p.value <- chisq.res$p.value

        } else {

          # Convert zeros to 1 to not obtain an error with chisq-test
          contingency[contingency == 0] <- 1

          # Run Fisher's exact test
          F.test <- fisher.test(contingency)

          # Add to results table
          stat.value   <- F.test$estimate
          stat.p.value <- F.test$p.value
        }
      }
      return(data.frame(stat = stat.value, stat.p.value = stat.p.value))
    }
    table <- cbind(table,xstat)
  } else {
    
    # Add empty columns for statistic and corresponding p-value
    table$stat         <- NA
    table$stat.p.value <- NA
    
    # Apply test for independence of sex and heterozygosity
    for (i in 1:nrow(table)) {

      # Exclude w.y-linked loci and loci with sex-biased score
      if (table[i, 11] == TRUE | table[i, "sex.biased"] == TRUE) {
        table[i, "stat"]         <- NA
        table[i, "stat.p.value"] <- NA

      } else {

        # Make contingency table
        contingency <- matrix(c(table[i, "count.F.het"],
                                table[i, "count.M.het"],
                                table[i, "count.F.hom"],
                                table[i, "count.M.hom"]),
                              nrow = 2,
                              ncol = 2,
                              dimnames = list(c("F", "M"), c("het", "hom")))

        # See if it's possible to use chisq test
        if (sum(contingency) >= 1000) {

          # Convert zeros to 1 to not obtain an error with chisq-test
          contingency[contingency == 0] <- 1

          # Chisq-test
          chisq.res <- chisq.test(contingency, correct = FALSE)

          # Add to results table
          table[i, "stat"]         <- chisq.res$statistic
          table[i, "stat.p.value"] <- chisq.res$p.value

        } else {

          # Convert zeros to 1 to not obtain an error with Fisher's test
          contingency[contingency == 0] <- 1

          # Run Fisher's exact test
          F.test <- fisher.test(contingency)

          # Add to results table
          table[i, "stat"]         <- F.test$estimate
          table[i, "stat.p.value"] <- F.test$p.value
        }
      }
    }
  }

  # Adjust p-values for multiple comparisons (False discovery rate, least conservative)
  table$stat.p.adjusted <- p.adjust(table$stat.p.value, method = "fdr")

  # Calculate for heterozygosity per sex and add to results table
  table$heterozygosity.F <- table$count.F.het/(table$count.F.het+table$count.F.hom)

  table$heterozygosity.M <- table$count.M.het/(table$count.M.het+table$count.M.hom)

  
  ##### 2.1 Z-linked or X-linked loci AND gametologs
  # For zw sex-determination system
  if(system == "zw") {
    table$z.linked     <- FALSE
    table$gametolog <- FALSE

    for (i in 1:nrow(table)) {
      # Exclude w-linked loci and loci with sex-biased score
      if (!is.na(table[i, "stat.p.adjusted"])) {
        # Exclude autosomal
        if (table[i, "stat.p.adjusted"] <= 0.01) {
          # Identify if heterozygosity is larger in males
          if (table[i, "heterozygosity.M"] > table[i, "heterozygosity.F"]) {
            table[i, "z.linked"] <- TRUE
          } else {
            table[i, "gametolog"] <- TRUE
          }
        }
      }
    }
  }
  
  # For xy sex-determination system
  if(system == "xy") {
    table$x.linked     <- FALSE
    table$gametolog <- FALSE

    for (i in 1:nrow(table)) {
      # Exclude y-linked loci and loci with sex-biased score
      if (!is.na(table[i, "stat.p.adjusted"])) {
        # Exclude autosomal
        if (table[i, "stat.p.adjusted"] <= 0.01) {
          # Identify if heterozygosity is larger in females
          if (table[i, "heterozygosity.F"] > table[i, "heterozygosity.M"]) {
            table[i, "x.linked"] <- TRUE
          } else {
            table[i, "gametolog"] <- TRUE
          }
        }
      }
    }
  }


  ##### 2.2 Plot BEFORE vs AFTER
  if(plots) {
    message("Building heterozygosity plots.")
    
    # For zw sex-determination system
    if(system == "zw") {
      BEF.het <- plot(x = table$heterozygosity.F,
                      y = table$heterozygosity.M,
                      xlab = "% Heterozygous Females",
                      ylab = "% Heterozygous Males",
                      main = "BEFORE",
                      xlim = c(0, 1),
                      ylim = c(0, 1))
  
      AFT.het <- plot(x = table[table$w.linked       == FALSE &
                                  table$sex.biased   == FALSE &
                                  table$z.linked     == FALSE &
                                  table$gametolog == FALSE, "heterozygosity.F"],
                      y = table[table$w.linked       == FALSE &
                                  table$sex.biased   == FALSE &
                                  table$z.linked     == FALSE &
                                  table$gametolog == FALSE, "heterozygosity.M"],
                      main = "AFTER",
                      xlab = "% Heterozygous Females",
                      ylab = "% Heterozygous Males",
                      xlim = c(0, 1),
                      ylim = c(0, 1))
    }
  
    # For xy sex-determination system
    if(system == "xy") {
      BEF.het <- plot(x = table$heterozygosity.F,
                      y = table$heterozygosity.M,
                      xlab = "% Heterozygous Females",
                      ylab = "% Heterozygous Males",
                      main = "BEFORE",
                      xlim = c(0, 1),
                      ylim = c(0, 1))
  
      AFT.het <- plot(x = table[table$y.linked     == FALSE &
                                  table$sex.biased   == FALSE &
                                  table$x.linked     == FALSE &
                                  table$gametolog == FALSE, "heterozygosity.F"],
                      y = table[table$y.linked     == FALSE &
                                  table$sex.biased   == FALSE &
                                  table$x.linked     == FALSE &
                                  table$gametolog == FALSE, "heterozygosity.M"],
                      main = "AFTER",
                      xlab = "% Heterozygous Females",
                      ylab = "% Heterozygous Males",
                      xlim = c(0, 1),
                      ylim = c(0, 1))
    }
    message("Done building heterozygosity plots.")
  }

  #################### 3. Create output of function
  ##### 3.1 Save the indices of each category of loci to later subset gl
  # For zw sex-determination system
  if(system == "zw") {
    a <- table[table$w.linked     == TRUE, "index"]
    b <- table[table$sex.biased   == TRUE, "index"]
    c <- table[table$z.linked     == TRUE, "index"]
    d <- table[table$gametolog == TRUE, "index"]

    autosomal <- table[table$w.linked     == FALSE &
                       table$sex.biased   == FALSE &
                       table$z.linked     == FALSE &
                       table$gametolog == FALSE, "index"]

    message("**FINISHED** Total of analyzed loci: ", nrow(table), ".\n",
            "Found ", length(a)+length(b)+length(c)+length(d), " sex-linked loci:\n",
            "   ",    length(a), " W-linked loci\n",
            "   ",    length(b), " sex-biased loci\n",
            "   ",    length(c), " Z-linked loci\n",
            "   ",    length(d), " ZW gametologs.\n",
            "And ",   length(autosomal), " autosomal loci.")
  }

  if(system == "xy") {
    a <- table[table$y.linked     == TRUE, "index"]
    b <- table[table$sex.biased   == TRUE, "index"]
    c <- table[table$x.linked     == TRUE, "index"]
    d <- table[table$gametolog == TRUE, "index"]

    autosomal <- table[table$y.linked     == FALSE &
                       table$sex.biased   == FALSE &
                       table$x.linked     == FALSE &
                       table$gametolog == FALSE, "index"]

    message("**FINISHED** Total of analyzed loci: ", nrow(table), ".\n",
            "Found ", length(a)+length(b)+length(c)+length(d), " sex-linked loci:\n",
            "   ",    length(a), " Y-linked loci\n",
            "   ",    length(b), " sex-biased loci\n",
            "   ",    length(c), " X-linked loci\n",
            "   ",    length(d), " XY gametologs.\n",
            "And ",   length(autosomal), " autosomal loci.")
  }


  ##### 3.2 Subset gl object (according to IntroTutorial_dartR.pdf, page 32)
  if (length(a) > 0){
    A <- x[, a]  # Loci are columns
  } else {
    A <- NULL
  }
  
  if (length(b) > 0){
    B <- x[, b]
  } else {
    B <- NULL
  }
  
  if (length(c) > 0){
    C <- x[, c]
  } else {
    C <- NULL
  }
  
  if (length(d) > 0){
    D <- x[, d]
  } else {
    D <- NULL
  }
  
  gl.autosomal <- gl[ , autosomal]



  #################### 4. Output
  if(system == "xy"){
    rlist <- list(  "results.table" = table,
                    "y.linked"      = A,
                    "sex.biased"    = B,
                    "x.linked"      = C,
                    "gametolog"     = D,
                    "autosomal"     = gl.autosomal)}

  if(system == "zw"){
    rlist <- list(  "results.table" = table,
                    "w.linked"      = A,
                    "sex.biased"    = B,
                    "z.linked"      = C,
                    "gametolog"     = D,
                    "autosomal"     = gl.autosomal)}
  end_time0 <- Sys.time()
  print(sprintf("Run time: %s", end_time0-start_time0))
  
  if(parallel){
    stopCluster(cl)
  }
  
  return(rlist)
}
################################################################################


################################ Example of use ################################
## filtered.data <- filter.sex.linked(gl = my.genlight,                       ##
##                                    system = "xy",                          ##
##                                    plots = FALSE,                          ##
##                                    parallel = TRUE)                        ##
################################################################################


