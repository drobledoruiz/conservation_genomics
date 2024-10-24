# filter.excess.het

A function to filter out excessively-heterozygous SNPs that are likely to be spurious multilocus SNPs (i.e., products of erroneously fusing multiple physically separate loci that are very similar). This function considers a locus as excessively-heterozygous if it (i) presents an heterozygosity > 0.5, *and* (ii) presents a significant excess of observed heterozygous individuals from those expected according to Hardy-Weinberg equilibrium. These criteria are evaluated per population. This is, if a locus is found to fullfill both (i) and (ii) in *any* population, it is considered as excessively-heterozygous and is removed.

This function requires as input:

- **gl** - A genlight object in which ‘ind.metrics’ contains a column named ‘pop’ and each individual is assigned to one population.
- **Yates** - Boolean (optional), to use Yates's continuity correction. Recommended for sample sizes < 20. *Set to FALSE by default.*

This function produces as output:

- A dataframe with information about excessively-heterozygous loci (removed loci).
- A genlight object without excessively-heterozygous loci.
- A vector with the names of the removed loci (i.e., excessively-heterozygous loci).
- Two plots: one BEFORE plot with the heterozygosity of the loci present in the *input* genlight, and one AFTER plot with the heterozygosity of the loci present in the *output* genlight (i.e., without excessively-heterozygous loci).


## Usage

To use the function it is necesary to save the file *filter.excess.het.R*, and load it to *R*:

```
source('/path_to_function/filter.excess.het.R')
```

Then the function can be used:

```
filtered.data <- filter.excess.het(gl = my.genlight,           
                                   Yates = TRUE) 
```

The output *filtered.data* contains three elements that can be called as:

- **filtered.data$results.table** - Dataframe with information about excessively-heterozygous loci (removed loci).
- **filtered.data$filtered.gl** - Genlight object without excessively-heterozygous loci.
- **filtered.data$removed.loci** - Vector with the names of the excessively-heterozygous loci (i.e. removed loci).

The results table has the folowing structure:

```
loci      pop   n0   n1  n2   Hobs   en0     en1     en2    Hexp  chsq   p.value   p.adjusted
locus1    pop1  284  524  77   0.59  336.85  418.29  129.85  0.47  56.52  5.56-14   1.15e-13
```

With columns meaning:
- **loci** -- Name of the locus.
- **pop** -- Name of the population in which the locus was found to be excessively-heterozygous.
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

In the example above, locus1 exhibits (i) an heterozygosity > 0.5, and (ii) an adjusted p-value <= 0.05 in population "pop1", which means that the excess of observed heterozygous individuals is significant (n1 > en1). Hence, this locus is considered excessively-heterozygous and is removed from the filtered genlight object.


---------------------------------------------------------------------------
Contact:
- Diana Robledo-Ruiz, diana.robledoruiz1@monash.edu
- Jesus Castrejon-Figueroa, jcastrejon@ciencias.unam.mx
