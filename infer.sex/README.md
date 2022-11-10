# infer.sex

This function uses the output of function *filter.sex.linked* (list of 6 genlight objects) to infer the sex of all individuals. It uses 3 types of sex-linked loci (W-linked/Y-linked, Z-linked/X-linked, and gametologs), assigns a genetic sex for each type, and outputs an agreed sex assigned by at least 2 types of sex-linked loci. **We created this function with the explicit intent that a human checks the evidence for the sex assignments that do NOT agree for all 3 types of sex-linked loci** (denoted as '*M' or '*F'). This human can then use their criterion to validate the assignment.

The function requires as input:
- **gl_sex_filtered** - List of 6 genlight objects created by the function *filter.sex.linked*.
- **system** - The sex determination system of the species. Default is 'zw' but can be set to 'xy'.
- **seed** - Integer, seed for random number generator (< 2³¹) used on the KMeans algorithm. *Set randomly by default.*

It produces as output a dataframe with 11 columns:
- **id** - Individuals ID.
- **w.linked.sex** or **y.linked.sex** - Sex inferred using W-linked or Y-linked loci.
- **#called** - Number of W-linked or Y-linked loci for which the individual had a called genotype (cf. missing genotype).
- **#missing** - Number of W-linked or Y-linked loci for which the individual had a missing genotype (cf. called genotype).
- **z.linked.sex** or **x.linked.sex** - Sex inferred using Z-linked or X-linked loci.
- **#Hom.z** or **#Hom.x** - Number of Z-linked or X-linked loci for which the individual is homozygous.
- **#Het.z** or **#Het.x** - Number of Z-linked or X-linked loci for which the individual is heterozygous.
- **gametolog.sex** - Sex inferred using ZW-gametologs or XY-gametologs.
- **#Hom.g** - Number of gametolog loci for which the individual is homozygous.
- **#Het.g** - Number of gametolog loci for which the individual is heterozygous.
- **agreed.sex** - Genetic sex assigned by at least two types of sex-linked loci.

## Dependencies

- [plyr](https://cran.r-project.org/web/packages/plyr/index.html)

## Usage

To use the function *infer.sex* it is necesary to load the file *infer.sex.R*:

```
source('/path_to_function/infer.sex.R')
```

Then the function can be called: 

```
inferred_sex <- infer.sex(gl_sex_filtered,
		      	  system='xy')
```

The object 'inferred_sex' is a dataframe whose first 4 rows read as: 

```
id       y.linked.sex  #called  #missing  x.linked.sex  #Het.x  #Hom.x  gametolog.sex  #Het.g  #Hom.g  agreed.sex
ind1     F      	1      	179       F     	208   	583      F     		410    	9     	F
ind2     M     		158     2         M      	3    	750      M     		184   	231    	M
ind3     F      	0     	180       F     	203   	579      M     		180   	239   	*F
ind4     F     		12     	114       M      	7    	776      M     		156   	248   	*M
```
The 3 types of sex-linked loci agree that ind1 is a female, and ind2 is a male. However, for ind3 Y- and X-linked loci suggest a female but XY-gametologs suggest a male. Similarly for ind4 X-linked loci and XY-gametologs suggest a male but Y-linked loci suggest a female. In these cases the sex is reported as probable female ('*F') and probable male ('*M'), respectively. **A human must check the evidence for these probable sex assignments and use their criterion to validate them or not.**

## How the function works
It assigns a genetic sex for 3 types of sex-linked loci by:
- **W-linked**: assigns 'M' to an individual if its genotypes for these loci are mostly missing (cf. called genotype), 'F' otherwise. Or **Y-linked**: assigns 'F' to an individual if its genotypes for these loci are mostly missing (cf. called genotype), 'M' otherwise.
- **Z-linked/X-linked**: it uses these loci to perform a k-means clustering of individuals in two clusters (F and M).
- **Gametologs**: it uses these loci to perform a k-means clustering of individuals in two clusters (F and M).
	
The function then outputs an agreed sex:
- 'F' or 'M' if the 3 types of sex-linked loci assigned the same sex.
- '*F' or '*M' if only 2 types of sex-linked loci assigned the same sex.


---------------------------------------------------------------------------
Contact:
- Diana Robledo-Ruiz, diana.robledoruiz1@monash.edu
- Jesus Castrejon-Figueroa, jcastrejon@ciencias.unam.mx
