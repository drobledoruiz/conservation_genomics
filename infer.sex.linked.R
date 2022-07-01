###########################################################################
# Infer sex from the three tests (W,Z,ZW) or (X,Y,XY). 
#    Recieves as argument:
#       * gl_sex_filtered <-  Sex-filtered genlight dataframe.
#       * link            <-  Linkage 'zw' for (W,Z,ZW) or 'xy' for (X,Y,XY).
#       * seed            <-  Random seed for K-Means clustering.
#
#    Output:
#       * A dataframe with 10 coumns. Includes the result for the sex on
#         the W,Z, ZW.glog test, together with the Number of homozygous (#Hom)
#         and heterozygous (#Het) on the Z and ZW.glog-tests. For the Z-test
#         includes the Number of score (#Scored) and missing (#NA) loci.
#              
########################################################################
library(plyr)

infer.sex <- function(gl_sex_filtered,link='zw',seed=42){

  w   <-  W.sex(gl_sex_filtered[[2]],   link=link)
  z   <-  Z.sex(gl_sex_filtered[[4]],   link=link, seed=seed)
  zwg <-  ZWg.sex(gl_sex_filtered[[5]], link=link, seed=seed)  
  A <- data.frame(w,z,zwg)

Fun<- function(x,y,z){
  d <- c(x,y,z) 
  yy<-count(d)
  yy<-yy[order(-yy$freq),]  
  value <- if(yy[1,2]==3) yy[1,1] else sprintf('*%s',yy[1,1])
  return(value)
}

  A$SEX <- mapply(Fun, A$W.sex, A$Z.sex, A$ZWg.sex)
return(A)
  if(link=='xy'){
    names <- c( 'Y.link', '#NA', '#Scored', 
                'X.link', '#Het',   '#Hom', 
                'XY.glog','#Het',   '#Hom', 'SEX')
  }else{
    names <- c( 'W.link', '#Scored','#NA',
                'Z.link', '#Hom','#Het', 
                'ZW.glog','#Hom','#Het', 'SEX')                 
  }

  colnames(A) <- names
  return(A)
}


#-------------------------------------------------------------------
# Map NAs and Scores to 1Dim in [-1,1], IF x<0, F else M
#-------------------------------------------------------------------

W.sex <- function(gl,link='zw'){

  w <- as.matrix(gl)
  w[is.na(w)] <- 1

  n0.w <- rowSums(w == 0 | w == 2, na.rm=TRUE)
  n1.w <- rowSums(w == 1, na.rm=TRUE)

  sex.score <- function(f,m){
  return( (-f+m)/(f+m) )}

  c2 <- sex.score(n0.w,n1.w)
  if(link=='xy'){
    lab0 <- 'M'
    lab1 <-'F'
  }else{
    lab0 <- 'F'
    lab1 <- 'M'
  }
  W.sex <- ifelse( c2 < 0, lab0, lab1)

  Y <- data.frame(W.sex,n0.w,n1.w)
  return(Y)
}

#---------------------------------------------------------------------------
# Map Hom and Het to 2Dim and apply kmeans. Choose the label from maximum Hom
#---------------------------------------------------------------------------
Z.sex <- function(gl,link='zw',seed=42){
  z <- as.matrix(gl)

  n0.z = rowSums(z == 0 | z == 2, na.rm=TRUE)
  n1.z = rowSums(z == 1, na.rm=TRUE)

  Z <- t(apply(data.frame(n0.z,n1.z), 1, function(x) x / sum(x) ))

  set.seed(seed)
  km <- kmeans(Z,2)
  
  i <- which.max(n1.z) # Largest proportion of '1'
  label<- km$cluster[i]

  if(link=='xy'){
    lab0 <- 'M'
    lab1 <- 'F'
  }else{
    lab0 <- 'F'
    lab1 <- 'M'
  }

  Z.sex <-  ifelse( km$cluster  == label, lab1, lab0)
  Y <- data.frame(Z.sex,n1.z,n0.z)
  return(Y)
}

#-----------------------------------------------------------------------------
# Map Hom and Het to 2Dim and apply kmeans. Choose the label from maximum Het
#-----------------------------------------------------------------------------
ZWg.sex <-  function(gl,link='zw',seed=42){
  zwg = as.matrix(gl)

  n0.zw = rowSums(zwg == 0 | zwg == 2, na.rm=TRUE)
  n1.zw = rowSums(zwg== 1, na.rm=TRUE)

  Z <- t(apply(data.frame(n0.zw,n1.zw), 1, function(x) x / sum(x) ))

  set.seed(seed)
  km <- kmeans(Z,2)
  
  i <- which.max(n1.zw) # Largest proportion of '1'
  label<- km$cluster[i]

  if(link=='xy'){
    lab0 = 'M'
    lab1 = 'F'
  }else{
    lab0 = 'F'
    lab1 = 'M'
  }

  ZWg.sex <-  ifelse( km$cluster  == label, lab0, lab1)
  Y <- data.frame(ZWg.sex,n0.zw,n1.zw)
  return(Y)
}
############################## END  ##################################