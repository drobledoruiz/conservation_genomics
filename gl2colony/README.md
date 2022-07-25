# gl2colony

Export a [COLONY2](https://www.zsl.org/science/software/colony) file from a genlight object. 

This function requires as input:

 - **gl** -- A genlight object in which genlight@others$ind.metrics has three columns named ["father","mother","offspring"] with boolean values ("yes" or "no") indicating if the individual is considered to be "father","mother", and/or "offspring". The function ignores uppercases in column names and values.
 
- **filename_out** -- String, name of the output colony2-formated file.

- Colony2 project-specific parameters: 
	- **project_name** -- String, name of the project, should be a string containing less than 40 letters and numbers. *Set to 'my_project' by default.* 
	
	- **output_name** -- String, output file name which should be a string containing less than 40 letters and numbers. All output files of Colony will use the same output file name but different extension names. *Set to 'my_project' by default.*

	- **probability_father** -- Real, estimated probability that a father of an offspring is included in the CMS samples. *Set to 0.5 by default.*
	
	-  **probability_mother** -- Real, estimated probability that a mother of an offspring is included in the CFS samples. *Set to 0.5 by default.*
	
	- **seed** -- Integer, seed for random number generator (< 2³¹). *Set a random seed by default.*

	- **update_allele_freq** -- Boolean (1 or 0), instructs
Colony to update allele frequencies during the simulated annealing process in searching for the ML configuration. *Set to 0 by default.*

	- **di_mono_ecious** -- Boolean (2 or 1), indicate dioecious
species or monoecious species. *Set to 2 by default.*

	- **inbreed** -- Boolean (0 or 1), indicate the absence or
presence of inbreeding. *Set to 0 by default.*

	- **haplodiploid** -- Boolean (0 or 1), indicate diploid species or haplodiploid species. For monoecious, the indicator value should always be 0 to indicate diploid. *Set to 0 by default.*

	- **polygamy_male** -- Boolean (0 or 1), specify whether males are  polygamous or monogamous. *Set to 0 by default.*
	
	- **polygamy_female** -- Boolean (0 or 1), specify whether females are  polygamous or monogamous. Set to 0 by default.*

	- **clone_inference** -- Boolean (1 or 0), specify whether clones
(or duplicated individuals) are to be inferred (=1) or not (=0). *Set to 1 by default.*
	- **scale_shibship** -- Boolean (1 or 0), specify whether full
sibship size is to be scaled or not. *Set to 1 by default.*

	- **sibship_prior** -- Integer, indicates average paternal sibship size; no sibship prior (0), weak (1), medium (2) or strong (3) sibship prior. Otional: (real) average maternal sibship size. *Set to 0 by default.*
	
	- **known_allele_freq** -- Boolean (1 or 0), indicates whether population allele frequencies for each locus are known and are to be provided or not.

	- **num_runs** -- Integer, number of replicate runs for the dataset. *Set to 1 by default.*

	- **length_run** -- Integer (1, 2, 3, 4), indicates short, medium, long, very long run. *Set to 2 by default.*

	- **monitor_method** -- - Boolean (1 or 0), indicates monitoring the intermediate results by iterate number or running time. Always choose value 0 for run without Windows GUI. *Set to 0 by default.*

- Other colony2 parameters and their default values:
	- **monitor_interval** -- *Set to 10000 by default.*
	- **windows_gui** -- *Set to 0 by default.*
	- **likelihood** -- *Set to 0 by default.*
	- **precision_fl** -- *Set to 2 by default.*
	- **marker_id** -- *Set to 'mk@' by default.*
	- **marker_type** -- *Set to '0@' by default.*
	- **allelic_dropout** -- *Set to '0.000@' by default.*
	- **other_typ_err** -- *Set to '0.05@' by default.*
	- **paternity_exclusion_threshold** -- *Set to '0 0' by default.*
	- **maternity_exclusion_threshold** -- *Set to '0 0' by default.*
	- **paternal_sibship** -- *Set to 0 by default.*
	- **maternal_sibship** -- *Set to 0 by default.*
	- **excluded_paternity** -- *Set to 0 by default.*
	- **excluded_maternity** -- *Set to 0 by default.*
	- **excluded_paternal_sibships** -- *Set to 0 by default.*
	- **excluded_maternity_sibships** -- *Set to 0 by default.*

We refer the user to the [Colony2 user manual](https://usermanual.wiki/Document/ColonyUserGuide.68067402) for more details on these parameters.
 

## Dependencies

- [crayon](https://cran.r-project.org/web/packages/crayon/index.html)


## Intructions

To use the function *gl2colony* it is necesary to load the file *gl2colony.R*:

```
source('/path_to_function/gl2colony.R')
```

then the function can be called. Here we show an example of its use with some non-default parameters: 

```
gl2colony(  gl = my.genlight,                                               
            filename_out = "colony2.dat",                                     
            project_name = "parentage_fish_2022",                             
            output_name =  "parentage_fish_jul_2022",                       
            seed = 1234,                                                    
            probability_father = 0.6,                                       
            probability_mother = 0.4,                                       
            update_allele_freq = 1,                                         
            allelic_dropout = '0.01@',                                      
            other_typ_err = '0.001@'  )
```

creating the file "colony2.dat" whose first 26 columns read as:

```
parentage_fish_2022 
parentage_fish_jul_2022 
711 		 ! No. offspring 
10403 		 ! No. of loci 
1234 		 ! Seed for random number generator 
1 		 ! 0/1=Not updating/updating allele frequency 
2 		 ! 2/1=Dioecious/Monoecious species 
0 		 ! 0/1=no inbreeding/inbreeding 
0 		 ! 0/1=Diploid species/HaploDiploid species 
0 0 		 ! 0/1=Polygamy/Monogamy for males & females 
1 		 ! 0/1=Clone/duplicates inference =No/Yes 
1 		 ! 0/1=Scale full sibship=No/Yes 
0 		 ! 0/1/2/3=No/Weak/Medium/Strong sibship prior 
0 		 ! 0/1=Unknown/Known population allele frequency 
1 		 ! Number of runs 
2 		 ! 1/2/3/4=short/medium/long/very long run 
0 		 ! 0/1=Monitor method by Iterate#/Time in second 
10000 		 ! Monitor interval in Iterate# / in seconds 
0 		 ! 0/1=No/Yes for run with Windows GUI 
0 		 ! 0/1/2=PairLikelihood score/Fulllikelihood/FPLS 
2 		 ! 0/1/2/3=Low/Medium/High/Very high precision FL 
 		  
mk@ 		 ! Marker Ids (consecutive for all) 
0@ 		 ! Marker types, 0/1=Codominant/Dominant 
0.01@ 		 ! Allelic dropout rate for all loci 
0.001@ 		 ! Other typing error rate for all loci 
```


---------------------------------------------------------------------------
Contact:
- Diana Robledo-Ruiz, diana.robledoruiz1@monash.edu
- Jesus Castrejon-Figueroa, jcastrejon@ciencias.unam.mx
