# infer.sex

Given a list of 6 genlight sex-filtered objects created as the output of the function filter.sex.linked, this function creates a dataframe with the infered sexes of the individuals.

This function requires as input:

- **gl_sex_filtered** --  List of 6 genlight objects created from the function filter.sex.linked.
- **link** --  String, sex determination system "xy" or "zw". *Set to "zw"by default.*
- **seed** --  Integer, seed for random number generator (< 2³¹) used on the KMeans algorithm. *Set a random seed by default.*

## Dependencies

- [plyr](https://cran.r-project.org/web/packages/plyr/index.html)

## Intructions

To use the function *infer.sex* it is necesary to load the file *infer.sex.R*:

```
source('/path_to_function/infer.sex.R')
```

then the function can be called. Here we show an example of its use: 

```
sex_fish <- infer.sex(	gl_sex_filtered,
			system='xy'      )
```

creating the dataframe "sex_fish" whose first 4 columns read as: 

```
         Y.link  #NA  #Scored  X.link  #Het  #Hom  XY.glog  #Het  #Hom  SEX
ind1       F      1     179      F     208   583      F     410    9     F
ind2       M     158     2       M      3    750      M     184   231    M
ind3       F      0     180      F     203   579      M     180   239   *F
ind4       F     12     114      M      7    776      M     156   248   *M
```

The proportion of #NA and #scores in the Y-test suggests ind1 is a female. Because the X-test and XY-test also suggest ind1 is a female the sex of the individual is reported as female (F) on column "SEX". Analogously for ind2 the three test suggest is a male. However for ind3 Y-test and X-test suggest a female but XY-test suggest a male. Similarly for ind4 X-test and XY-test suggest a male but Y-test suggest a female. In these cases the sex is reported as probable female (*F) for ind3 and probable male (*M) for ind4.

---------------------------------------------------------------------------
Contact:
- Diana Robledo-Ruiz, diana.robledoruiz1@monash.edu
- Jesus Castrejon-Figueroa, jcastrejon@ciencias.unam.mx
