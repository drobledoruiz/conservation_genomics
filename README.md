# swiss_army_knife for genetic management

A tool set for quick and easy handling of genetic data to do common analyses for the genetic management of endangered populations.

Contact: Diana Robledo-Ruiz, diana.robledoruiz1@monash.edu

-------------------------------------------------------------------------------------
**Function _filter.sex.linked_**

It identifies four types of sex-linked loci using individuals with assigned sex. 

This function requires as input:
  - A genlight object in which genlight@others$ind.metrics has a column named 'sex' and individuals are assigned 'F' or 'M'. The function ignores individuals that are assigned anything else or nothing at all.
  - The sex determination system of the species. Default is 'zw' but can be set to 'xy'.
  
This function produces as output:
  - A list with 6 elements:
      1. A Results Table (see description below).
      2. Genlight object with w-linked/y-linked loci.
      3. Genlight object with sex-biased scoring rate loci.
      4. Genlight object with z-linked/x-linked loci.
      5. Genlight object with zw-gametolog/xy-gametolog loci.
      6. Genlight object with **autosomal loci**.
  - Four plots:
      1. A plot BEFORE filtering sex-linked loci by scoring rate.
      2. A plot AFTER filtering sex-linked loci by scoring rate.
      3. A plot BEFORE filtering sex-linked loci by heterozygosity.
      4. A plot AFTER filtering sex-linked loci by heterozygosity.

Usage:
```
filtered.data <- filter.sex.linked(gl = my.genlight,
                                   system = "xy")
```

The function works in 2 phases:
1. Use loci scoring rate per sex to identify w-linked/y-linked loci and loci with sex-biased scoring rate.
2. Use the proportion of heterozygous males and females per loci to identify z-linked/x-linked loci and zw-gametologs.

Phase 1:
  1. It creates a Results Table with loci in rows.
  2. It adds columns 'count.F.miss' and 'count.M.miss' with counts of the number of females and males with NA (missing data) for each locus, respectively.
  3. It adds columns 'count.F.scored' and 'count.M.scored' with counts of the number of females and males scored (0, 1 or 2) for each locus, respectively.
  4. It builds a contingency table and performs a Fisher's exact test to test for the independence of scoring rate and sex per locus. It then adds to the Results Table a column with the Fisher's exact test estimate (column 'ratio') and its respective p-value (column 'p.value'). The rationale is that autosomal loci should present no difference in scoring rate between the sexes, and therefore, a locus in which scoring rate is biased by sex (i.e. individuals of one sex have significantly more or fewer missing data than expected compared to the other sex) is likely to be sex-linked.
  5. It adjusts p-values to control for the false discovery rate and adds column 'p.adjusted' to the Results Table.
  6. It adds columns 'scoringRate.F' and 'scoringRate.M' with the proportion of females and males scored (0, 1 or 2) for each locus, respectively. It outputs a plot with 'Female scoring rate' in the x axis and 'Male scoring rate' in the y axis in which each point is a locus. This is the BEFORE filtering plot. Autosomal loci should have roughly the same scoring rate for males and females, forming a cloud of points in a diagonal line.
  7. It adds column 'w.linked' in which a locus is signalled as w-linked (TRUE) if its scoring rate for males is smaller or equal to 0.1 (because males have no W chromosome) and its adjusted p-value is smaller or equal to 0.01. Or it adds column 'y.linked' in which a locus is signalled as y-linked (TRUE) if its scoring rate for females is smaller or equal to 0.1 (because females have no Y chromosome) and its adjusted p-value is smaller or equal to 0.01. All other loci take value FALSE.
  8. It adds column 'sex.biased' in which a locus is signalled as having sex-biased scoring rate (TRUE) if its adjusted p-value is smaller or equal to 0.01 and it has not been signalled as w-linked/y-linked.
  9. It outputs the same plot as step 6 but removing w-linked/y-linked and sex-biased scoring rate loci. This is the AFTER filtering plot and should have only loci (points) that are roughly in the diagonal line.
  
  Phase 2:
  1. It adds columns 'count.F.het' and 'count.M.het' with counts of the number of females and males scored as heterozygotes (scored '1') for each locus, respectively.
  2. It adds columns 'count.F.hom' and 'count.M.hom' with counts of the number of females and males scored as homozygotes ('0' or '2') for each locus, respectively.
  3. It builds a contingency table and performs a Fisher's exact test to test for the independence of heterozygosity and sex per locus. It then adds to the Results Table a column with the Fisher's exact test estimate (column 'stat') and its respective p-value (column 'stat.p.value'). The rationale is that autosomal loci should present no difference in heterozygosity rate between the sexes, and therefore, a locus in which heterozygosity is biased by sex (i.e. there are significantly more or fewer heterozygote individuals from one sex than the other sex) is likely to be sex-linked.
  4. It adjusts p-values to control for the false discovery rate and adds column 'stat.p.adjusted' to the Results Table.
  5. It adds columns 'heterozygosity.F' and 'heterozygosity.M' with the proportion of females and males that are heterozygous ('1') for each locus, respectively. It outputs a plot with 'Proportion of heterozygous females' in the x axis and 'Proportion of heterozygous males' in the y axis in which each point is a locus. This is the BEFORE filtering plot. Autosomal loci should have roughly the same proportion of heterozygous males and females, forming a cloud of points in a diagonal line.
  6. It adds column 'z.linked' in which a locus is signalled as z-linked (TRUE) if its adjusted p-value is < 0.05 and the proportion of heterozygous males is greater than the proportion of heterozygous females (because females only have one Z chromosome). Or it adds column 'x.linked' in which a locus is signalled as x-linked (TRUE) if its adjusted p-value is < 0.05 and the proportion of heterozygous females is greater than the proportion of heterozygous males (because males only have one X chromosome).
  7. It adds column 'zw.gametolog' in which a locus is signalled as zw-gametolog (TRUE) if its adjusted p-value is < 0.05 and the proportion of heterozygous males is smaller than the proportion of heterozygous females. Or it adds column 'xy.gametolog' in which a locus is signalled as xy-gametolog (TRUE) if its adjusted p-value is < 0.05 and the proportion of heterozygous females is smaller than the proportion of heterozygous males.
  8. It outputs the same plot as step 5 but removing z-linked/x-linked and gametolog loci. This is the AFTER filtering plot and should have only loci (points) that are roughly in the diagonal line.
