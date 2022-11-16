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
##   A dataframe with 11 columns, the last one 'agreed.sex' being the final   ##
##     genetic sex assigned using sex-linked loci information.                ##
##                                                                            ##
##  Index:                                                                    ##
##    Line 25: Function infer.sex                                             ##
##    Line 253: Example of use for infer.sex                                  ##
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

  # Make sex assignment per type of sex-linked loci (Functions declared below)
  # W/Y-linked
  if (gl1@n.loc >= 3){
    w   <-  W.sex(  gl1, system = system)
  } else {
    message("Not enough W-linked/Y-linked loci. Assigning NA...")
    w <- data.frame(W.sex = rep(NA, length(gl1@ind.names)),
                    n0.w  = rep(NA, length(gl1@ind.names)),
                    n1.w  = rep(NA, length(gl1@ind.names)))
  }
  # Z/X-linked
  if (gl2@n.loc >= 2){
    z   <-  Z.sex(  gl2, system = system, seed = seed)
  } else {
    message("Not enough Z-linked/X-linked loci. Assigning NA...")
    z <- data.frame(Z.sex = rep(NA, length(gl2@ind.names)),
                    n1.z  = rep(NA, length(gl2@ind.names)),
                    n0.z  = rep(NA, length(gl2@ind.names)))
  }
  # Gametologues
  if (gl3@n.loc >= 2){
    g   <-  g.sex(  gl3, system = system, seed = seed)
  } else {
    message("Not enough gametologues. Assigning NA...")
    g <- data.frame(g.sex = rep(NA, length(gl3@ind.names)),
                    n1.g  = rep(NA, length(gl3@ind.names)),
                    n0.g  = rep(NA, length(gl3@ind.names)))
  }

  # Put them all together
  A   <-  data.frame(w, z, g)

  # Function to conciliate assignments
  Fun <- function(x, y, z){
    d  <- c(x, y, z)
    yy <- data.frame(x = c("F", "M"),
                     freq = c(sum(d == "F", na.rm = TRUE),
                              sum(d == "M", na.rm = TRUE)))
    yy <- yy[order(-yy$freq),]
    value <- if(length(unique(na.omit(d))) == 1) unique(na.omit(d)) else sprintf('*%s', yy[1, 1])
    return(value)
  }

  A$agreed.sex <- mapply(Fun, A$W.sex, A$Z.sex, A$g.sex)

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

  message("***FINISHED***")

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
  z <- as.matrix(gl)

  n0.z = rowSums(z == 0 | z == 2, na.rm = TRUE)
  n1.z = rowSums(z == 1, na.rm = TRUE)

  Z_unclean <- t(apply(data.frame(n0.z, n1.z), 1, function(x) x / sum(x) ))
  Z <- na.omit(Z_unclean)
  Zna <- Z_unclean[rowSums(is.na(Z_unclean)) > 0,]

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

  # Assign sex
  Z.sex <- ifelse( km$cluster  == label, lab1, lab0)

  # Assign NA to inds with NA
  if (nrow(Zna) > 0){
    for (i in 1:nrow(Zna)) {
      Z.sex[length(Z.sex)+1] <- NA
      names(Z.sex)[length(names(Z.sex))] <- rownames(Zna)[i]
    }
  }

  # Fuse them
  Y <- data.frame(n1.z, n0.z)

  # Add NAs in appropriate row
  Y$Z.sex <- "STOP"
  for (i in 1:nrow(Y)){
    Y[i, "Z.sex"] <- Z.sex[rownames(Y)[i]]
  }

  Y <- Y[, c("Z.sex", "n1.z", "n0.z")]

  return(Y)
}

############################### 3. ZWg.sex function
### Map Hom and Het to 2Dim and apply kmeans. Choose the label from maximum Het

g.sex <-  function(gl, system = 'zw', seed = 42) {
  z <- as.matrix(gl)

  n0.g = rowSums(z == 0 | z == 2, na.rm = TRUE)
  n1.g = rowSums(z == 1, na.rm = TRUE)

  Z_unclean <- t(apply(data.frame(n0.g, n1.g), 1, function(x) x / sum(x) ))
  Z <- na.omit(Z_unclean)
  Zna <- Z_unclean[rowSums(is.na(Z_unclean)) > 0,]

  # Apply k-means
  set.seed(seed)
  km <- kmeans(Z, 2)

  i <- which.max(n1.g) # Largest proportion of '1'
  label <- km$cluster[i]

  if(system == 'xy'){
    lab0 <- 'M'
    lab1 <- 'F'
  } else {
    lab0 <- 'F'
    lab1 <- 'M'
  }

  # Assign sex (HERE IS THE OPPOSITE)
  Z.sex <- ifelse( km$cluster  == label, lab0, lab1)

  # Assign NA to inds with NA
  if (nrow(Zna) > 0){
    for (i in 1:nrow(Zna)) {
      Z.sex[length(Z.sex)+1] <- NA
      names(Z.sex)[length(names(Z.sex))] <- rownames(Zna)[i]
    }
  }

  # Fuse them
  Y <- data.frame(n1.g, n0.g)

  # Add NAs in appropriate row
  Y$g.sex <- "STOP"
  for (i in 1:nrow(Y)){
    Y[i, "g.sex"] <- Z.sex[rownames(Y)[i]]
  }

  Y <- Y[, c("g.sex", "n1.g", "n0.g")]

  return(Y)

}
################################################################################


################################ Example of use ################################
##  inferred.sex <- infer.sex(filtered_gl,                                    ##
##                            system = "xy")                                  ##
################################################################################
