# Bioinformatics tools for conservation genomics

A tool set written in R for quick and easy handling of genomic data to do common analyses for the genetic management of endangered populations. These functions are introduced in the preprint Robledo-Ruiz, et al. **Easy-to-use R functions to separate reduced-representation genomic datasets into sex-linked and autosomal loci, and conduct sex-assignment**. Authorea. December 02, 2022 (DOI: 10.22541/au.166998603.34117027/v1). **Do check each function directory for further details on usage.**

## Functions

1. **filter.sex.linked:** filters sex-linked SNPs from a genlight object. 
2. **infer.sex:** infers the sex of individuals using sex-linked SNPs.
3. **filter.excess.het:** filters excessively-heterozygous loci.
4. **gl2colony:** exports a COLONY2 file from a genlight object.

## Usage

The functions require as input a genlight object. A genlight object can be considered a matrix of individuals (rows) and SNPs (columns) in which the data is coded as '0' for homozygous reference, '1' for heterozygous, and '2' for homozygous alternate. Genlight objects were originally created by the adegenet R package, and have been adopted by dartR package. Genlight objects can contain individuals metadata as an attachment in slot @other$ind.metrics. Individuals metadata come in a dataframe with individuals in rows and their attributes in columns. Some of our functions require ind.metrics to have a column named 'sex' or a column 'pop'. For a great explanation on genlight objects and how dartR handles them check http://georges.biomatix.org/storage/app/media/uploaded-files/tutorial3adartrdatastructuresandinput22-dec-21-2.pdf. If you got DArTseq data, the easiest way to obtain a genlight object is:

```
# Install dartR
install.packages("dartR")
library(dartR)

# Import data as genlight object
gl <- gl.read.dart(filename = "sample_data_2Row.csv", ind.metafile = "individual_metadata.csv")
```

---------------------------------------------------------------------------
Contact:
- Diana Robledo-Ruiz, diana.robledoruiz1@monash.edu
- Jesus Castrejon-Figueroa, jcastrejon@ciencias.unam.mx
