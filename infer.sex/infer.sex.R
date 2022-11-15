################################################################################
##                 Function to infer sex from sex-linked loci                 ##
##                                                                            ##
##  Authors: Jesús Castrejón-Figueroa. R developer, Monash University         ##
##           Diana A Robledo-Ruiz. Research Fellow, Monash University         ##
##  Date: 2022-07-01                                                          ##
##                                                                            ##
##  This function requires:                                                   ##
##   - Input: the output of function filter.sex.linked (complete list with 6  ##
##            elements).                                                      ##
##   - User-specified parameters:                                             ##
##       - system = sex determination system ('zw' or 'xy').                  ##
##       - seed = interger, chosen randomly by default.                       ##
##                                                                            ##
##  Output:                                                                   ##
##   A dataframe with 11 columns, the last one 'agreed.sex' being the genetic ##
##     sex assigned by at least two types of sex-linked loci.                 ##
##                                                                            ##
##  Index:                                                                    ##
##    Line 25: Function infer.sex                                             ##
##    Line 177: Example of use for filter.highly.het                          ##
################################################################################


########################## Define function infer.sex ###########################
library(plyr)

infer.sex <- function(gl_sex_filtered, system = 'zw', seed = NULL) {

  # Random seed if not specified by user
  if(is.null(seed)) {
    seed <- sample.int(65535, 1)
  }

  if(system == "xy") {
    gl1 <- gl_sex_filtered$y.linked
    gl2 <- gl_sex_filtered$x.linked
  }

  if(system == "zw") {
    gl1 <- gl_sex_filtered$w.linked
    gl2 <- gl_sex_filtered$z.linked
  }

  gl3 <- gl_sex_filtered$gametolog

  # Functions declared below
  w   <-  W.sex(  gl1, system = system)
  z   <-  Z.sex(  gl2, system = system, seed = seed)
  zwg <-  ZWg.sex(gl3, system = system, seed = seed)
  A   <-  data.frame(w, z, zwg)

  Fun <- function(x, y, z){
    d  <- c(x, y, z)
    yy <- count(d)
    yy <- yy[order(-yy$freq),]
    value <- if(yy[1, 2] == 3) yy[1, 1] else sprintf('*%s', yy[1, 1])
    return(value)
  }

  A$agreed.sex <- mapply(Fun, A$W.sex, A$Z.sex, A$ZWg.sex)

  if(system == 'xy'){
    names <- c('y.linked.sex',  '#missing', '#called',
               'x.linked.sex',  '#Het.x',   '#Hom.x',
               'gametolog.sex', '#Het.g',   '#Hom.g', 'agreed.sex')
  } else {
    names <- c('w.linked.sex',  '#called', '#missing',
               'z.linked.sex',  '#Hom.z',  '#Het.z',
               'gametolog.sex', '#Hom.g',  '#Het.g', 'agreed.sex')
  }

  colnames(A) <- names

  A <- cbind(row.names(A), A)

  colnames(A)[1] <- "id"

  return(A)
}

############################### 1. W.sex function
### Map NAs (missing) and scored (called) to 1Dim in [-1,1], IF x<0, F else M

W.sex <- function(gl, system = 'zw'){
  w <- as.matrix(gl)
  w[is.na(w)] <- 3

  n0.w <- rowSums(w == 0 | w == 2 | w == 1, na.rm = TRUE)
  n1.w <- rowSums(w == 3, na.rm = TRUE)

  # Calculate proportion
  sex.score <- function(f, m){
    return( (-f+m)/(f+m) )
  }

  c2 <- sex.score(n0.w, n1.w)

  if(system == 'xy'){
    lab0 <- 'M'
    lab1 <-'F'
  } else {
    lab0 <- 'F'
    lab1 <- 'M'
  }

  W.sex <- ifelse(c2 < 0, lab0, lab1)

  Y <- data.frame(W.sex, n0.w, n1.w)
  return(Y)
}

############################### 2. Z.sex function
### Map Hom and Het to 2Dim and apply kmeans. Choose the label from maximum Hom

Z.sex <- function(gl, system = 'zw', seed = 42){
  z_unclean <- as.matrix(gl)   
  zna <- z_unclean[rowSums(is.na(z_unclean)) > 0,]
  z <- na.omit(z_unclean)

  n0.z = rowSums(z == 0 | z == 2, na.rm = TRUE)
  n1.z = rowSums(z == 1, na.rm = TRUE)

  Z <- t(apply(data.frame(n0.z, n1.z), 1, function(x) x / sum(x) ))

  # Apply k-means
  set.seed(seed)
  km <- kmeans(Z, 2)

  i <- which.max(n1.z) # Largest proportion of '1'
  label <- km$cluster[i]

  if(system == 'xy'){
    lab0 <- 'M'
    lab1 <- 'F'
  } else {
    lab0 <- 'F'
    lab1 <- 'M'
  }

  Z.sex <- ifelse( km$cluster  == label, lab1, lab0)
  Y <- data.frame(Z.sex, n1.z, n0.z)

  Nna <- dim(zna)[1]
  if(Nna == 0){
    return(Y)
  } else{
    rnames <- row.names(zna)
    cnames <- names(zna)
    Yna <- data.frame(rep(NA,nna), rep(0,nna), rep(0,nna))
    row.names(Yna) <- rnames
    names(Yna) <- cnames
    return( rbind( Y , Yna) )
  }
}

############################### 3. ZWg.sex function
### Map Hom and Het to 2Dim and apply kmeans. Choose the label from maximum Het

ZWg.sex <-  function(gl, system = 'zw', seed = 42) {
  zwg_unclean = as.matrix(gl)


  zna <- zwg_unclean[rowSums(is.na(zwg_unclean)) > 0,]
  zwg <- na.omit(zwg_unclean)

  n0.zw = rowSums(zwg == 0 | zwg == 2, na.rm = TRUE)
  n1.zw = rowSums(zwg == 1, na.rm = TRUE)

  Z <- t(apply(data.frame(n0.zw, n1.zw), 1, function(x) x / sum(x) ))

  # Apply k-means
  set.seed(seed)
  km <- kmeans(Z, 2)

  i <- which.max(n1.zw) # Largest proportion of '1'
  label <- km$cluster[i]

  if(system == 'xy') {
    lab0 = 'M'
    lab1 = 'F'
  } else {
    lab0 = 'F'
    lab1 = 'M'
  }

  ZWg.sex <- ifelse( km$cluster  == label, lab0, lab1)
  Y <- data.frame(ZWg.sex, n0.zw, n1.zw)
  
  Nna <- dim(zna)[1]
  if(Nna == 0){
    return(Y)
  } else{
    rnames <- row.names(zna)
    cnames <- names(zna)
    Yna <- data.frame(rep(NA,nna), rep(0,nna), rep(0,nna))
    row.names(Yna) <- rnames
    names(Yna) <- cnames
    return( rbind( Y , Yna) )
  }
}
################################################################################


################################ Example of use ################################
##  inferred.sex <- infer.sex(filtered_gl,                                    ##
##                            system = "xy")                                  ##
################################################################################
