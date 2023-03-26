# filter.sex.linked

This function identifies four types of sex-linked loci using individuals with known sex. 

This function requires as input:
  - A genlight object in which genlight@others$ind.metrics has a column named 'sex' and individuals are assigned 'F' or 'M'. The function ignores individuals that are assigned anything else or nothing at all.
  - The sex determination system of the species: 'zw' or 'xy'.

This function has two parameters that can be modified by the user:
  - 'plots', specifies whether the function should produce the four output plots ('TRUE' or 'FALSE'). *Set to 'TRUE' by default.*
  - 'parallel', specifies whether the function should run in parallel, which saves time ('TRUE' or 'FALSE'). *Set to 'FALSE' by default.*

This function produces as output:
  - A list with 6 elements:
      1. A Results Table (see description below).
      2. Genlight object with **w-linked/y-linked** loci.
      3. Genlight object with **sex-biased** call rate loci.
      4. Genlight object with **z-linked/x-linked** loci.
      5. Genlight object with **gametologous** loci.
      6. Genlight object with **autosomal** loci.
  - Four plots:
      1. A plot BEFORE filtering sex-linked loci by call rate.
      2. A plot AFTER filtering sex-linked loci by call rate.
      3. A plot BEFORE filtering sex-linked loci by heterozygosity.
      4. A plot AFTER filtering sex-linked loci by heterozygosity.

## Dependencies
- [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)



## Usage

To use the function it is necesary to save the file *filter.sex.linked.R*, and load it to *R*:

```
source('/path_to_function/filter.sex.linked.R')
```

Then the function can be called:
```
filtered.data <- filter.sex.linked(gl = my.genlight,
                                   system = "xy",
                                   plots = TRUE,
                                   parallel = TRUE)
```

Setting parameter 'plots = FALSE' saves a little bit of running time for very large datasets (> 50,000 SNPs), especially if not run in parallel. However, we **strongly** encourage to always inspect the output plots at least once to make sure everything is working properly.

The output *filtered.data* contains 6 elements that can be called as:

   - **filtered.data$results.table** - Dataframe with information about highly-heterozygous loci (filtered out loci).
   - **filtered.data$w.linked** or **filtered.data$y.linked** - Genlight object with w-linked/y-linked loci.
   - **filtered.data$sex.biased**    - Genlight object with sex-biased call rate loci.
   - **filtered.data$z.linked** or **filtered.data$x.linked**     - Genlight object with z-linked/x-linked loci.
   - **filtered.data$gametolog**     - Genlight object with gametologous loci.
   - **filtered.data$autosomal**     - Genlight object with autosomal loci.

The results table has the folowing structure:
```
                index  count.F.miss  count.M.miss  count.F.scored  count.M.scored  ratio    p.value  p.adjusted   scoringRate.F  scoringRate.M  w.linked  sex.biased  count.F.het  count.M.het  count.F.hom  count.M.hom       stat   stat.p.value   stat.p.adjusted  heterozygosity.F  heterozygosity.M  z.linked  zw.gametolog
28689726-20-T/A     1             0             0             173             224      0  1.0000000           1       1.0000000      1.0000000     FALSE       FALSE          119          157           54           67  0.9405887     0.82619338         1.0000000         0.6878613         0.7008929     FALSE         FALSE
```

With loci being in rows and columns meaning:

- **index** - Index number to identify loci.
- **count.F.miss** - Count of females that have this locus as missing data (NA).
- **count.M.miss** - Count of males that have this locus as missing data (NA).
- **count.F.scored** - Count of females that have this locus scored (0, 1 or 2; i.e. non-missing).
- **count.M.scored** - Count of males that have this locus scored (0, 1 or 2; i.e. non-missing).
- **ratio** - Fisher's exact test estimate testing for the independence of call rate and sex for this locus.
- **p.value** - P-value for the Fisher's exact test estimate.
- **p.adjusted** - P-value adjusted for false discovery rate.
- **scoringRate.F** - Female call rate (proportion of females that were scored for this locus). This is the x-axis in the 1st and 2nd plot.
- **scoringRate.M** - Male call rate (proportion of males that were scored for this locus). This is the y-axis in the 1st and 2nd plot.
- **w.linked** - Boolean for this locus being w-linked.
- **sex.biased** - Boolean for this locus having sex-biased call rate.

- **count.F.het** - Count of females that are heterozygous for this locus.
- **count.M.het** - Count of males that are heterozygous for this locus.
- **count.F.hom** - Count of females that are homozygous for this locus.
- **count.M.hom** - Count of males that are homozygous for this locus.
- **stat** - Fisher's exact test estimate testing for the independence of heterozygosity and sex for this locus.
- **stat.p.value** - P-value for the Fisher's exact test estimate.
- **stat.p.adjusted** - P-value adjusted for false discovery rate.
- **heterozygosity.F** - Female heterozygosity (proportion of females that were heterozygotes for this locus). This is the x-axis in the 3rd and 4th plot.
- **heterozygosity.M** - Male heterozygosity (proportion of males that were heterozygotes for this locus). This is the y-axis in the 3rd and 4th plot.
- **z.linked** - Boolean for this locus being z-linked.
- **zw.gametolog** - Boolean for this locus being a zw-gametolog.


## How the function works
The function works in 2 phases:
1. Use loci call rate per sex to identify w-linked/y-linked loci and loci with sex-biased call rate.
2. Use the proportion of heterozygous males and females per loci to identify z-linked/x-linked loci and gametologs.

**Phase 1:**
  1. It creates a Results Table with loci in rows.
  2. It adds columns 'count.F.miss' and 'count.M.miss' with counts of the number of females and males with NA (missing data) for each locus, respectively.
  3. It adds columns 'count.F.scored' and 'count.M.scored' with counts of the number of females and males scored (0, 1 or 2) for each locus, respectively.
  4. It builds a contingency table and performs a Fisher's exact test to test for the independence of call rate and sex per locus. It then adds to the Results Table a column with the Fisher's exact test estimate (column 'ratio') and its respective p-value (column 'p.value'). The rationale is that autosomal loci should present no difference in call rate between the sexes, and therefore, a locus in which call rate is biased by sex (i.e. individuals of one sex have significantly more or fewer missing data than expected compared to the other sex) is likely to be sex-linked.
  5. It adjusts p-values to control for the false discovery rate and adds column 'p.adjusted' to the Results Table.
  6. It adds columns 'scoringRate.F' and 'scoringRate.M' with the proportion of females and males scored (0, 1 or 2) for each locus, respectively. It outputs a plot with 'Female call rate' in the x axis and 'Male call rate' in the y axis in which each point is a locus. This is the BEFORE filtering plot. Autosomal loci should have roughly the same call rate for males and females, forming a cloud of points in a diagonal line.
  7. It adds column 'w.linked' in which a locus is signalled as w-linked (TRUE) if its call rate for males is smaller or equal to 0.1 (because males have no W chromosome) and its adjusted p-value is smaller or equal to 0.01. Or it adds column 'y.linked' in which a locus is signalled as y-linked (TRUE) if its call rate for females is smaller or equal to 0.1 (because females have no Y chromosome) and its adjusted p-value is smaller or equal to 0.01. All other loci take value FALSE.
  8. It adds column 'sex.biased' in which a locus is signalled as having sex-biased call rate (TRUE) if its adjusted p-value is smaller or equal to 0.01 and it has not been signalled as w-linked/y-linked.
  9. It outputs the same plot as step 6 but removing w-linked/y-linked and sex-biased call rate loci. This is the AFTER filtering plot and should have only loci (points) that are roughly in the diagonal line.
  
**Phase 2:**
  1. It adds columns 'count.F.het' and 'count.M.het' with counts of the number of females and males scored as heterozygotes (scored '1') for each locus, respectively.
  2. It adds columns 'count.F.hom' and 'count.M.hom' with counts of the number of females and males scored as homozygotes ('0' or '2') for each locus, respectively.
  3. It builds a contingency table and performs a Fisher's exact test to test for the independence of heterozygosity and sex per locus. It then adds to the Results Table a column with the Fisher's exact test estimate (column 'stat') and its respective p-value (column 'stat.p.value'). The rationale is that autosomal loci should present no difference in heterozygosity rate between the sexes, and therefore, a locus in which heterozygosity is biased by sex (i.e. there are significantly more or fewer heterozygote individuals from one sex than the other sex) is likely to be sex-linked.
  4. It adjusts p-values to control for the false discovery rate and adds column 'stat.p.adjusted' to the Results Table.
  5. It adds columns 'heterozygosity.F' and 'heterozygosity.M' with the proportion of females and males that are heterozygous ('1') for each locus, respectively. It outputs a plot with 'Proportion of heterozygous females' in the x axis and 'Proportion of heterozygous males' in the y axis in which each point is a locus. This is the BEFORE filtering plot. Autosomal loci should have roughly the same proportion of heterozygous males and females, forming a cloud of points in a diagonal line.
  6. It adds column 'z.linked' in which a locus is signalled as z-linked (TRUE) if its adjusted p-value is < 0.01 and the proportion of heterozygous males is greater than the proportion of heterozygous females (because females only have one Z chromosome). Or it adds column 'x.linked' in which a locus is signalled as x-linked (TRUE) if its adjusted p-value is < 0.01 and the proportion of heterozygous females is greater than the proportion of heterozygous males (because males only have one X chromosome).
  7. It adds column 'zw.gametolog' in which a locus is signalled as zw-gametolog (TRUE) if its adjusted p-value is < 0.01 and the proportion of heterozygous males is smaller than the proportion of heterozygous females. Or it adds column 'xy.gametolog' in which a locus is signalled as xy-gametolog (TRUE) if its adjusted p-value is < 0.01 and the proportion of heterozygous females is smaller than the proportion of heterozygous males.
  8. It outputs the same plot as step 5 but removing z-linked/x-linked and gametolog loci. This is the AFTER filtering plot and should have only loci (points) that are roughly in the diagonal line.
  
  
---------------------------------------------------------------------------
Contact:
- Diana Robledo-Ruiz, diana.robledoruiz1@monash.edu
- Jesus Castrejon-Figueroa, jcastrejon@ciencias.unam.mx
