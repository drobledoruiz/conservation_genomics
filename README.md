# Easy-to-use R functions for conservation genomics

A toolset written in *R* for quick and easy handling of genomic data to do common analyses for the genetic management of endangered populations. These functions are introduced in the preprint: Robledo-Ruiz et al. **Easy-to-use R functions to separate reduced-representation genomic datasets into sex-linked and autosomal loci, and conduct sex-assignment**. Authorea. December 02, 2022 (https://www.authorea.com/doi/full/10.22541/au.166998603.34117027/v1). **Do check each function directory for further details on usage.**

## Functions

1. **filter.sex.linked:** filters sex-linked SNPs from a genlight object. 
2. **infer.sex:** infers the sex of individuals using sex-linked SNPs.
3. **filter.excess.het:** filters excessively-heterozygous loci.
4. **gl2colony:** exports a COLONY2 file from a genlight object.

## Usage

The functions require as input a genlight object. Genlight objects were originally created by the *adegenet* package, and have been adopted by *dartR* package. A genlight object can be considered a matrix of individuals (rows) and SNPs (columns) in which the genotypes are coded as '0' for homozygous reference, '1' for heterozygous, and '2' for homozygous alternate. Genlight objects can contain metadata for individuals in the form of a dataframe that is attached to the genlight object (slot '@other$ind.metrics'). In this dataframe, individuals are in rows and their attributes in columns. Some of our functions require 'ind.metrics' to have a column named 'sex' or a column 'pop'. 

For a great explanation on how to obtain a genlight object **from DArTseq data** check: http://georges.biomatix.org/storage/app/media/uploaded-files/tutorial3adartrdatastructuresandinput22-dec-21-2.pdf. 

For an explanation on how to obtain genlight objects from **non-DArTseq data (e.g., RADseq)** check page 13 of: http://georges.biomatix.org/storage/app/media/uploaded-files/tutorial3bdartrdatastructuresandinputfromsourcesotherthandartlmagv2-2.pdf. 

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
Contact:
- Diana Robledo-Ruiz, diana.robledoruiz1@monash.edu
- Jesus Castrejon-Figueroa, jcastrejon@ciencias.unam.mx
