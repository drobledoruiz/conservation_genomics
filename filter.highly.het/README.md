# filter.highly.het

A function to filter highly heterozygous loci based on Chi-Square tests for Hardy-Weinberg equilibrium. 

This function requires as input:

- **gl_autosom_filtered** --   A sex-filtered automal genlight object. This is the sixth element of the list created by the function filter.sex.linked.

- **Yates** -- Boolean (optional), to use Yates's continuity correction. *Set to FALSE by default.*


## Intructions

To use the function *infer.sex* it is necesary to load the file *ifilter.highly.het.R*:

```
source('/path_to_function/filter.highly.het.R')
```

then the function can be called. Here we show an example of its use: 

```
het.filtered.gl  <- filter.highly.het(  gl_autosom,           
                                        Yates = TRUE) 
```

The object *het.filtered.gl* contains two elements:

- **het.filtered.gl$filtered.gl**  --  The original (gl_autosom) genlight object without the highly heterozygous loci.         
- **het.filtered.gl$highly.het.loci** --  DataFrame with information about the removed loci.

Tha DataFrame *het.filtered.gl$highly.het.loci* has the folowing structure:

```
           n0   n1  n2   Hobs   en0     en1     en2    Hexp  chsq   p.value   p.adjusted
locus     284  524  77   0.59  336.85  418.29  129.85  0.47  56.52  5.56-14   1.15e-13
```

wiht the columns meaning:

- **n0** -- Number of **AA** genotypes.
- **n1** -- Number of **Aa** genotypes.
- **n2** -- Number of **aa** genotypes.
- **Hobs** -- Observed proportion of heterozygous.
- **en0** -- Expected value for the number of **AA** genotypes.
- **en1** -- Expected value for the number of **Aa** genotypes.
- **en2** -- Expected value for the number of **aa** genotypes.
- **Hexp** -- Expected proportion of heterozygous.
- **chsq** -- Value of chi-square for these genotype proportions. 
- **p.value** -- P-value corresponding to the chi-square value.
- **p.adjusted** -- Adjusted p-value.

In this case the adjusted p-value is significantly smaller than 0.05. Hence this locus would be considered highly heterozygous and would be removed from the genlight object.


---------------------------------------------------------------------------
Contact:
- Diana Robledo-Ruiz, diana.robledoruiz1@monash.edu
- Jesus Castrejon-Figueroa, jcastrejon@ciencias.unam.mx
