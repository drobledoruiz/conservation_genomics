################################################################################
##     Function to filter highly-heterozygous loci from a genlight object     ##
##                                                                            ##
##  Authors: Jesús Castrejón-Figueroa. R developer, Monash University         ##
##           Diana A Robledo-Ruiz. Research Fellow, Monash University         ##
##  Date: 2022-07-21                                                          ##
##                                                                            ##
##  This function requires:                                                   ##
##   - Input: a genlight object                                               ##
##   - User specified parameter: Boolean for Yates's continuity correction.   ##
##                                                                            ##
##  Output:                                                                   ##
##    $results.table - Dataframe with information on highly-heterozygous loci ##
##    $filtered.gl   - Genlight object without highly-heterozygous loci       ##
##                                                                            ##
##  Index:                                                                    ##
##    Line 21: Function filter.highly.het                                     ##
##    Line 126: Example of use for filter.highly.het                          ##
################################################################################

#################### Define function filter.highly.het ####################
filter.highly.het <- function(gl, Yates = FALSE){

  if(!Yates) {
    cc = 0
  } else {
    cc = 0.5
  }

  ################# 1. Start BEFORE plot and results table with observed data
  # Plot
  message("Building BEFORE-filtering plot")

  gen <- as.data.frame(t(as.matrix(gl)))
  n0 <- rowSums(gen == 0, na.rm = TRUE)
  n1 <- rowSums(gen == 1, na.rm = TRUE)
  n2 <- rowSums(gen == 2, na.rm = TRUE)

  fhe <- n1/(n0 + n1 + n2)

  plt.BEF <- plot(fhe,
                  main = "BEFORE",
                  ylab = "Locus heterozygosity",
                  ylim = c(0, 1),
                  xlim = c(0, gl@n.loc))

  # Results per population
  populations <- as.vector(unique(gl@other$ind.metrics$pop))
  n0   <- vector()
  n1   <- vector()
  n2   <- vector()
  pops <- vector()
  loci <- vector()

  for(pop in populations) {
    pop.gl <- gl[gl@other$ind.metrics$pop == pop, ]
    gen    <- as.data.frame(t(as.matrix(pop.gl)))
    n0     <- c(n0, rowSums(gen == 0, na.rm = TRUE))
    n1     <- c(n1, rowSums(gen == 1, na.rm = TRUE))
    n2     <- c(n2, rowSums(gen == 2, na.rm = TRUE))
    loci   <- c(loci, rownames(gen))
    pops   <- c(pops, rep(pop, dim(gen)[[1]]))
  }

  table <- data.frame(loci = loci,
                      pop  = pops,
                      n0   = n0,
                      n1   = n1,
                      n2   = n2,
                      Hobs = n1/(n0 + n1 + n2))

  ################# 2. Function to test HW with chisq
  HWchsq <- function(n0, n1, n2, cc) {
    n    <- n0 + n1 + n2
    p    <- (2*n0 + n1)/(2*n)
    q    <- 1 - p
    en0  <- n*p*p
    en1  <- 2*n*p*q
    en2  <- n*q*q
    Hexp <- en1/(en0 + en1 + en2)
    chsq <- (abs(n0 - en0) - cc)**2/en0 + (abs(n1 - en1) - cc)**2/en1 + (abs(n2 - en2) - cc)**2/en2
    return(c(en0, en1, en2, Hexp, chsq))
  }

  ################# 3. Function to test HW for each locus
  HW.chsqtest <- function(table, cc = cc){
    table$En0     <- NA
    table$En1     <- NA
    table$En2     <- NA
    table$Hexp    <- NA
    table$chsq    <- NA
    table$p.value <- NA

    for(i in 1:nrow(table)){
      n0 <- table$n0[i]
      n1 <- table$n1[i]
      n2 <- table$n2[i]

      chsq <- HWchsq(n0, n1, n2, cc)

      table[i,'En0']     <- chsq[1]
      table[i,'En1']     <- chsq[2]
      table[i,'En2']     <- chsq[3]
      table[i,'Hexp']    <- chsq[4]
      table[i,'chsq']    <- chsq[5]
      table[i,'p.value'] <- pchisq(chsq[5], 1, lower.tail = FALSE)
    }
    return(table)
  }

  ################## 4. Apply chsq test only to loci with Het > 50%
  message("Testing loci for high heterozygosity...")

  table.filter <- na.omit(table[table$Hobs >= 0.5, ] )
  table.filter <- HW.chsqtest(table.filter, cc = cc)

  # Adjust p-values
  table.filter$p.adjusted <- p.adjust(table.filter$p.value, method = "fdr")

  # Keep in table only loci p.adj <= 0.05 and ObsHet >= ExpHet
  table.filter <- table.filter[table.filter$p.adjusted<= 0.05 & table.filter$n1>=table.filter$En1,]
  rownames(table.filter) <- 1:nrow(table.filter)

  # Remove highly-het loci from new filtered gl
  gl.filter <- gl[, !(gl$loc.names %in% table.filter$loci )]

  # Remove highly-het loci from new gl loci metadata
  gl.filter@other$loc.metrics <- gl@other$loc.metrics[!(gl$loc.names %in% table.filter$loci ), ]

  ################## 5. AFTER plot with filtered gl
  message("Building AFTER-filtering plot")

  gen <- as.data.frame(t(as.matrix(gl.filter)))
  n0 <- rowSums(gen == 0, na.rm = TRUE)
  n1 <- rowSums(gen == 1, na.rm = TRUE)
  n2 <- rowSums(gen == 2, na.rm = TRUE)

  fhe <- n1/(n0 + n1 + n2)

  plt.AFT <- plot(fhe,
                  main = "AFTER",
                  ylab = "Locus heterozygosity",
                  ylim = c(0, 1),
                  xlim = c(0, gl@n.loc))

  ################## 6. Output
  message("**FINISHED**. Removed ", gl@n.loc - gl.filter@n.loc, " highly-heterozygote loci.")

  return(list('filtered.gl'   = gl.filter,
              'results.table' = table.filter,
              'removed.loci'  = unique(table.filter$loci)))
}
################################################################################

################################ Example of use ################################
##  filtered.data <- filter.highly.het(gl = my.genlight,                      ##
##                                     Yates = TRUE)                          ##
##  filtered.data$results.table                                               ##
##  filtered.data$filtered.gl                                                 ##
################################################################################
