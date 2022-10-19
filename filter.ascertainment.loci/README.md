# filter.ascertainment.bias

This function reduces the ascertainment bias product of unequal sampling size for different populations. By definition, reduced-representation sequencing (RRS) 
identifies variable loci present in a sample of individuals. If the sample contains individuals from more than one population, and the sample size for one 
population is larger (e.g. n1 = 20, n2 = 40), there will be loci that are polymorphic only in the larger population simply because rare alleles are more
likely to be found the more individuals are included in a sample (i.e. ascertainment bias). This function reduces this form of ascertainment bias by 
removing loci that become monomorphic after equalizing sample sizes across all populations (e.g. n1 = 20, n2 = 20). The resultant dataset contains
all individuals (e.g. n1 = 20, n2 = 40) genotyped for only the retained loci.

This function requires as input:
  - A genlight object in which genlight@others$ind.metrics has a column named 'pop' and individuals assigned to at least 2 populations.
  - A seed number. Interger used to randomly subsample larger populations to equalize sample sizes. The use of a seed makes results repeatable. *Set to 1 by default*.
  - Maximum sample size (n). Interger used as the maximum sample size for the subsampling of larger populations. *Set to the smallest population size by default*. This can be set to be larger than the smallest sample size if willing to tolerate some ascertainment bias (e.g. instead of equalizing n1 = 20, n2 = 20; set n = 30 to obtain n1 = 20, n2 = 30).
  
This function produces as output:
  - **$filtered.gl** - A genlight object without ascertainment bias (loci removed).
  - **$asc.inds** - A vector with the names of the ascertainment individuals after equalization.
  - **$removed.loci** - A vector with the names of the removed loci.
  - **$results.table** - A table with per population (i) sample size, (ii) number of polymorphic loci *before* filtering out ascertainment bias, and (iii) number of polymorphic loci *after* filtering out ascertainment bias. 

The function also produces 3 plots (based on the information in $results.table):
  - Barplot of sample size per population.
  - Barplot of the number of polymorphic loci present per population *before* filtering out ascertainment bias.
  - Barplot of the number of polymorphic loci present per population *after* filtering out ascertainment bias.


## Usage
```
filtered.data <- filter.ascertainment.bias(gl = my.genlight,
                                           seed = 100,
                                           n = 30)
```

---------------------------------------------------------------------------
Contact:
- Diana Robledo-Ruiz, diana.robledoruiz1@monash.edu
- Jesus Castrejon-Figueroa, jcastrejon@ciencias.unam.mx
