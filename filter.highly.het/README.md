# filter.highly.het

A function to filter out highly-heterozygous loci. This function considers a locus as highly-heterozygous if it (i) presents an heterozygosity > 0.5, AND (ii) presents a significant excess of observed heterozygous individuals from those expected according to Hardy-Weinberg equilibrium. 

This function requires as input:

- **gl** - A genlight object.

- **Yates** - Boolean (optional), to use Yates's continuity correction. Recommended for sample sizes < 20. *Set to FALSE by default.*

This function produces as output:

- A dataframe with information about highly-heterozygous loci (filtered out loci).

- A genlight object without highly-heterozygous loci.


## Usage

To use the function it is necesary to load the file *filter.highly.het.R*:

```
source('/path_to_function/filter.highly.het.R')
```

Then the function can be called:

```
filtered.data <- filter.highly.het(gl = my.genlight,           
                                   Yates = TRUE) 
```

The output *filtered.data* contains two elements that can be called as:

- **filtered.data$results.table** - Dataframe with information about highly-heterozygous loci (filtered out loci).
- **filtered.data$filtered.gl** - Genlight object without highly-heterozygous loci.

The results table has the folowing structure:

```
           n0   n1  n2   Hobs   en0     en1     en2    Hexp  chsq   p.value   p.adjusted
locus1    284  524  77   0.59  336.85  418.29  129.85  0.47  56.52  5.56-14   1.15e-13
```

With columns meaning:

- **pop** -- Name of the population.
- **n0** -- Observed number of homozygous reference individuals.
- **n1** -- Observed number of heterozygous individuals.
- **n2** -- Observed number of homozygous alternate individuals.
- **Hobs** -- Observed heterozygosity (proportion of heterozygous individuals).
- **en0** -- Expected number of homozygous reference individuals according to HW equilibrium.
- **en1** -- Expected number of heterozygous individuals according to HW equilibrium.
- **en2** -- Expected number of homozygous alternate individuals according to HW equilibrium.
- **Hexp** -- Expected heterozygosity (proportion of heterozygous individuals).
- **chsq** -- Chi-square statistic for the deviation of observed genotypes from the expected genotypes. 
- **p.value** -- P-value for the chi-square statistic.
- **p.adjusted** -- P-value adjusted for false discovery rate.

In the example above, locus1 exhibits (i) an heterozygosity > 0.5, and (ii) an adjusted p-value <= 0.05, which means that the excess of observed heterozygous individuals is significant (n1 > en1). Hence, this locus is considered highly-heterozygous and is removed from the filtered genlight object.


---------------------------------------------------------------------------
Contact:
- Diana Robledo-Ruiz, diana.robledoruiz1@monash.edu
- Jesus Castrejon-Figueroa, jcastrejon@ciencias.unam.mx
