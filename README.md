# Easy-to-use R functions for conservation genomics

A toolset written in *R* for quick and easy handling of genomic data to do common analyses for the genetic management of endangered populations.

**If you use our functions, PLEASE CITE THIS ARTICLE:**. 
Robledo-Ruiz et al. (2023) Easy-to-use R functions to separate reduced-representation genomic datasets into sex-linked and autosomal loci, and conduct sex-assignment. _Molecular Ecology Resources_. (https://doi.org/10.1111/1755-0998.13844).

You can watch a short video abstract of the paper here: https://vimeo.com/840300860 

Do check each function directory for further details on usage.

## Functions

1. **filter.sex.linked:** filters sex-linked SNPs from a genlight object. 
2. **infer.sex:** infers the sex of individuals using sex-linked SNPs.
3. **filter.excess.het:** filters excessively-heterozygous loci.
4. **gl2colony:** exports a COLONY2 file from a genlight object.

## Usage

The functions require as input a genlight object. Genlight objects were originally created by the *adegenet* package, and have been adopted by *dartR* package. A genlight object can be considered a matrix of individuals (rows) and SNPs (columns) in which the genotypes are coded as '0' for homozygous reference, '1' for heterozygous, and '2' for homozygous alternate. Genlight objects can contain metadata for individuals in the form of a dataframe that is attached to the genlight object (slot '@other$ind.metrics'). In this dataframe, individuals are in rows and their attributes in columns. Some of our functions require 'ind.metrics' to have a column named 'sex' or a column 'pop'. 

For a great explanation on how to obtain a genlight object **from DArTseq data** check: http://georges.biomatix.org/storage/app/media/uploaded-files/tutorial3adartrdatastructuresandinput22-dec-21-2.pdf. 

If you have DArTseq data, the easiest way to obtain a genlight object is:
```
# Install dartR
install.packages("dartR")
library(dartR)

# Import data as genlight object
gl <- gl.read.dart(filename = "sample_data_2Row.csv", ind.metafile = "individual_metadata.csv")

# Call genlight object
gl

# Check individual's metafile
View(gl@other$ind.metrics)
```

For an explanation on how to obtain genlight objects from **non-DArTseq data (e.g., RADseq or WGRe-seq)** check page 13 of: http://georges.biomatix.org/storage/app/media/uploaded-files/tutorial3bdartrdatastructuresandinputfromsourcesotherthandartlmagv2-2.pdf. 

Alternatively, if you have a vcf file, the easiest way to obtain a genlight object is:
```
# Install vcfR
install.packages("vcfR")
library(vcfR)

# Import vcf
vcf <- read.vcfR("/path_to/file.vcf")

# Transform vcf to genlight object
gl <- vcfR2genlight(vcf)

# Call genlight object
gl

# Import individual's metafile
metafile <- read.csv(file = "/path_to/individual_metafile.csv")

# Add individual's metafile to genlight object
gl@other$ind.metrics <- metafile
```



You can test our functions on the small datasets that we include in 'data' directory. Download the data files and functions files to your computer and load them to *R*:
```
# Load testing dataset
load(file = "/path_to_data/gl_EYR.R")

# Call it
gl_EYR

# Load function
source("/path_to_function/filter.sex.linked.R")

# Apply function to gl
filtered_data <- filter.sex.linked(gl_EYR, system = "zw")
```

---------------------------------------------------------------------------
## Contact
Do not hesitate to write to us asking any questions. We are happy to help!
- Diana Robledo-Ruiz, diana.robledoruiz1@monash.edu
- Jesus Castrejon-Figueroa, j.castrejon@unsw.edu.au
