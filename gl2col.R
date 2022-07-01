###########################################################################
# Exports a Genlight file to Colony Format.
#    Recieves as argument:
#       * genlight_file     <-  Genlight dataframe.
#       * filename_out      <-  Path of the output Colony file.
#       * Options           <-  Options for Colony. Set to a default value
#                               (check Colony User-Manual for further information).
#
#    Output:
#       * A Colony-fomated text file.
#
#   Comments:
#       The information of offspring and parental information must be contained 
#       in genlight_file@other$ind.metrics within the columns 'mother', 'father',
#       and 'offspring' taking values 'yes'/'no' (indistinctly of upper/lowercases).
#
########################################################################

library(gdata)

gl2col <- function( genlight_file = NULL,           
                    filename_out = NULL,
                    project_name = 'my_project',
                    output_name = 'my_project',
                    probability_male = 0.5,
                    probability_female = 0.5,
                    seed = 1234,
                    update_allele_freq = 0,
                    di_mono_ecious = 2,
                    inbreed = 1,
                    diploid = 1,
                    polygamy_male = 0,
                    polygamy_female = 0,
                    clone_inference = 1,
                    scale_shibship = 1,
                    sibship_prior = 0,
                    known_allele_freq = 0,
                    num_runs = 1,
                    length_run = 2,
                    monitor_method = 0,
                    monitor_interval = 10000,
                    windows_gui = 0,
                    likelihood = 0,
                    precision_fl = 2,
                    marker_id = 'mk@',
                    marker_type = '0@',
                    allelic_dropout = '0.000@',
                    other_typ_err = '0.05@',
                    paternity_exclusion_threshold = '0 0',
                    maternity_exclusion_threshold = '0 0',
                    paternal_sibship = 0,
                    maternal_sibship = 0,
                    excluded_paternity = 0,
                    excluded_maternity = 0,
                    excluded_paternal_sibships = 0,
                    excluded_maternity_sibships = 0) {


if(is.null(genlight_file)){
    stop('Missing GenLight file.')
}

if(is.null(filename_out)){
    stop('Missing output filename.')
}

loci = dim(genlight_file)[2]

if(is.null(loci)){
    stop('Number of loci must be non-zero.')
} 

x <- sex.parental.ids(genlight_file)

offspring_ids <- x$offs
dad_ids <- x$dad
mum_ids <- x$mum

n_offspring <- length(offspring_ids)
n_males <- length(dad_ids)
n_females <- length(mum_ids)

message(sprintf('%d Mothers detected. \n%d Fathers detected. \n%d Offsprings detected.\n', n_females, n_males, n_offspring))

if( !length(dad_ids) ){
    stop('Missing parenthal IDs.')
}

if( !length(mum_ids) ){
    stop('Missing maternal IDs.')
}

if( !length(offspring_ids) ){
    stop('Missing offspring IDs.')
}

# -----------------------------------------------
message("Exporting Genlight object to Colony format...")
str_1row_with0s <- gl2struc(genlight_file)

# -----------------------------------------------
# Subset dataset with offspring
offspring_geno_to_keep <- match(offspring_ids, rownames(str_1row_with0s))  
offspring_gen <- str_1row_with0s[offspring_geno_to_keep,]  # keep only offspring genotypes

# Subset dataset with only females
mum_geno_to_keep <- match(mum_ids, rownames(str_1row_with0s))  
mum_gen <- str_1row_with0s[mum_geno_to_keep,]  # keep only female genotypes

# Subset dataset with only males
dad_geno_to_keep <- match(dad_ids, rownames(str_1row_with0s))  
dad_gen <- str_1row_with0s[dad_geno_to_keep,]  # keep only male genotypes
# -----------------------------------------------

head_comments <- list('! No. offspring',
'! No. of loci',
'! Seed for random number generator',
'! 0/1=Not updating/updating allele frequency',
'! 2/1=Dioecious/Monoecious species',
'! 0/1=no inbreeding/inbreeding',
'! 0/1=Diploid species/HaploDiploid species',
'! 0/1=Polygamy/Monogamy for males & females',
'! 0/1=Clone/duplicates inference =No/Yes',
'! 0/1=Scale full sibship=No/Yes',
'! 0/1/2/3=No/Weak/Medium/Strong sibship prior',
'! 0/1=Unknown/Known population allele frequency',
'! Number of runs',
'! 1/2/3/4=short/medium/long/very long run',
'! 0/1=Monitor method by Iterate#/Time in second',
'! Monitor interval in Iterate# / in seconds',
'! 0/1=No/Yes for run with Windows GUI',
'! 0/1/2=PairLikelihood score/Fulllikelihood/FPLS',
'! 0/1/2/3=Low/Medium/High/Very high precision FL',
'',
'! Marker Ids (consecutive for all)',
'! Marker types, 0/1=Codominant/Dominant',
'! Allelic dropout rate for all loci',
'! Other typing error rate for all loci')

polygamy <- paste(polygamy_male,polygamy_female,sep=' ')

head.list <- list(n_offspring,loci,seed,update_allele_freq,di_mono_ecious,inbreed,diploid,
                    polygamy,clone_inference,scale_shibship,sibship_prior,known_allele_freq,
                    num_runs,length_run,monitor_method,monitor_interval,windows_gui,likelihood,precision_fl,'',marker_id,
                    marker_type,allelic_dropout,other_typ_err
                    )
# -----------------------------------------------
sink(filename_out)
cat(project_name,'\n')
cat(output_name,'\n')
for(i in 1:length(head.list)){
cat(head.list[[i]],'\t\t',head_comments[[i]],'\n')
}
sink()
# -----------------------------------------------
write.table(offspring_gen,
            file =filename_out,
            append = TRUE,
            quote = FALSE, 
            col.names = FALSE)

sink(filename_out,append =TRUE)
#-------------------------------------------
probabilities <- paste(probability_male,probability_female,sep=' ')
n_indv <- paste(n_males,n_females,sep=' ')
cat('\n')
cat(probabilities,'\t\t','! Probabilities that the father and mother of an offspring are included in candidates','\n')
cat(n_indv,'\t\t','! Numbers of candidate males and females','\n')
cat('\n')
sink()
# -----------------------------------------------
write.table(dad_gen,
            file =filename_out,
            append = TRUE,
            quote = FALSE, 
            col.names = FALSE)
# -----------------------------------------------
sink(filename_out,append =TRUE)
cat('\n')
sink()
# -----------------------------------------------
write.table(mum_gen,
            file =filename_out,
            append = TRUE,
            quote = FALSE, 
            col.names = FALSE)
# -----------------------------------------------

last_comments <- list('! Number of offspring with known paternity, exclusion threshold',
                      '! Number of offspring with known maternity, exclusion threshold',
                      '',
                      '! Number of known paternal sibship',
                      '! Number of known maternal sibship',
                      '',
                      '! Number of offspring with known excluded paternity',
                      '! Number of offspring with known excluded maternity',
                      '',
                      '! Number of offspring with known excluded paternal sibships',
                      '! Number of offspring with known excluded maternal sibships')

last_values <- list(paternity_exclusion_threshold,
                    maternity_exclusion_threshold,
                    '',
                    paternal_sibship,
                    maternal_sibship,
                    '',
                    excluded_paternity,
                    excluded_maternity,
                    '',
                    excluded_paternal_sibships,
                    excluded_maternity_sibships)

sink(filename_out,append =TRUE)  
cat('\n')
for(i in 1:length(last_values)){
cat(last_values[[i]],'\t\t',last_comments[[i]],'\n')
}
sink()
# -----------------------------------------------
cat(crayon::green$bold('GenLight file succesfully exported!'))
}


#------------------------------------------------------------------------------------
# Define function to extract sex parental information 
sex.parental.ids <- function(gen_data){

# read metadata and convert to lowercase
indv.metadata <- gen_data@other$ind.metrics
names(indv.metadata) <- tolower(names(indv.metadata))

# remove leading/trailing white spaces
indv.metadata$mother <- tolower(trim(indv.metadata$mother))
indv.metadata$father <- tolower(trim(indv.metadata$father))
indv.metadata$offspring <- tolower(trim(indv.metadata$offspring))

mum_ids <- indv.metadata[indv.metadata$mother == 'yes', 'id']
dad_ids <- indv.metadata[indv.metadata$father == 'yes', 'id']
offs_ids <- indv.metadata[indv.metadata$offspring == 'yes', 'id']

mum_ids <- as.vector(na.omit(mum_ids))
dad_ids <- as.vector(na.omit(dad_ids))
offs_ids <- as.vector(na.omit(offs_ids))

x = list(offs = offs_ids, dad = dad_ids, mum = mum_ids)
return(x)
}

#------------------------------------------------------------------------------------
# Define function to convert gl matrix to Struct format and convert 2-row-per-ind structure format to 1-row-per-ind
gl2struc <- function(x,
                    addtlColumns = NULL, 
                    ploidy = 2,
                    exportMarkerNames = FALSE){    

    genmat <- as.matrix(x)
    indNames <- dimnames(genmat)[[1]]
    nInd <- dim(genmat)[1] # number of individuals
    # make sets of possible genotypes
    G <- list()
    for(i in 0:ploidy){
        G[[i + 1]] <- c(rep(1, ploidy - i), rep(2, i))
    }
    #G[[ploidy + 2]] <- rep(-9, ploidy) # for missing data
    G[[ploidy + 2]] <- rep(0, ploidy) # for missing data

    # set up data frame for Structure
    StructTab <- data.frame(ind = rep(indNames, each = ploidy))
    # add any additional columns
    if(!is.null(addtlColumns)){
        for(i in 1:dim(addtlColumns)[2]){
            StructTab <- data.frame(StructTab, rep(addtlColumns[,i], each = ploidy))
            if(!is.null(dimnames(addtlColumns)[[2]])){
                names(StructTab)[i + 1] <- dimnames(addtlColumns)[[2]][i]
            } else {
                names(StructTab)[i + 1] <- paste("X", i, sep = "")
            }
        }
    }

    # add genetic data
    for(i in 1:dim(genmat)[2]){
        thesegen <- genmat[,i] + 1
        thesegen[is.na(thesegen)] <- ploidy + 2
        StructTab[[dimnames(genmat)[[2]][i]]] <- unlist(G[thesegen])
    }
    
    # return(StructTab)  # for returning the value of gl2struct dartR function 
    
    data <- StructTab
    # Define dimensions of the matrix (only genotypes, not Ids)
    out <- matrix(NA, nrow = (nrow(data) / 2),        # no. of rows divided by 2
                    ncol = (2 * (ncol(data) - 1)))  # no. of columns minus Ids column times 2
  
    # Select first row per ind, leaving behind first column (Ids), then assign as first column per ind
    out[, seq(1, ncol(out), by = 2)] <- as.matrix(data[seq(1, nrow(data), by = 2), -1])
    # Select second row per ind, leaving behind first column (Ids), then assign as second column per ind
    out[, seq(2, ncol(out), by = 2)] <- as.matrix(data[seq(2, nrow(data), by = 2), -1])
  
    # Select Id column (only first row per ind) and make it rownames for matrix
    rownames(out) <- data[seq(1, nrow(data), by = 2), 1]  
    return(out) 
}
################################## END ######################################