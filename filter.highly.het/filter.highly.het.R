########################################################################################
##  Function to filter highly heterozygous loci from an autosomal genlight object     ##
##  based on the Chi-Square tests for Hardy-Weinberg equilibrium.                     ##
##                                                                                    ##
##     Recieves as argument:                                                          ##
##       * gl_autosom_filtered <-   Sex-filtered automal genlight object.             ##
##       * Yates (optional)    <-   Boolean for Yates's continuity correction.        ##
##                                                                                    ##
##     Output:                                                                        ##
##       $filtered.gl     :  Genlight object without highly heterozygous loci         ##
##       $highly.het.loci :  DataFrame with information about the removed loci        ##
##                                                                                    ##
########################################################################################

###################### Define MAIN function filter.highly.het ##########################

filter.highly.het <- function(gl,Yates=FALSE){
#---------------------------------------------------
  if(!Yates){
    cc = 0
  }else{
    cc = 0.5
  }
#----------------------------------------------------
gen <- as.data.frame(t(as.matrix(gl)))
n0 <- rowSums(gen == 0, na.rm = TRUE)
n1 <- rowSums(gen == 1, na.rm = TRUE)
n2 <- rowSums(gen == 2, na.rm = TRUE)
table <- data.frame(n0=n0,n1=n1,n2=n2,Hobs=n1/(n0+n1+n2))
#-----------------------------------------------------
HWchsq <- function(n0,n1,n2,cc){
  n <- n0+n1+n2
  p <- (2*n0 + n1)/(2*n)
  q <- 1-p
  en0 <- n*p*p
  en1 <- 2*n*p*q
  en2 <- n*q*q
  Hexp <- en1/(en0+en1+en2)
  chsq <- (abs(n0-en0) - cc)**2/en0 + (abs(n1-en1) - cc)**2/en1 + (abs(n2-en2) - cc)**2/en2
  return(c(en0,en1,en2,Hexp,chsq))
}
#-----------------------------------------------------
HW.chsqtest <- function(table,cc=cc){
  table$En0 <- NA
  table$En1 <- NA
  table$En2 <- NA
  table$Hexp <- NA
  table$chsq <- NA
  table$p.value <- NA

  for(i in 1:nrow(table)){
      n0 <- table$n0[i]
      n1 <- table$n1[i]
      n2 <- table$n2[i]
  
      chsq <- HWchsq(n0,n1,n2,cc)
      
      table[i,'En0'] <- chsq[1]
      table[i,'En1'] <- chsq[2]
      table[i,'En2'] <- chsq[3]
      table[i,'Hexp'] <- chsq[4]
      table[i,'chsq']<- chsq[5]
      table[i,'p.value'] <- pchisq(chsq[5], 1, lower.tail=FALSE)
  }
  return(table)
}
#----------------------------------------------------
# Apply chsq test only to the loci with Het > 50%
table.filter <- table[table$Hobs >= 0.5,] 
table.filter <- HW.chsqtest(table.filter,cc=cc)
table.filter$p.adjusted <- p.adjust(table.filter$p.value, method = "fdr")

# Keep in HetTable only values with p.adjusted <= 0.05 and ObsHet >= ExpHet
table.filter <- table.filter[table.filter$p.adjusted<= 0.05 & table.filter$n1>=table.filter$En1,]
#-----------------------------------------------------
return(list('filtered.gl'=gl[,!(gl$loc.names %in% rownames(table.filter))],
            'highly.het.loci'=table.filter))
}
######################################################################


############################ Example of use ##########################
##                                                                  ##
##  het.filtered.gl  <- filter.highly.het(  gl_autosom,             ##
##                                          Yates = TRUE)           ##
##                                                                  ##
######################################################################