# PLACO
Created by: Yazdan Asgari<br>
Creation date: 24 Dec 2021<br>
Update: 24 Dec 2021<br>
https://cesp.inserm.fr/en/equipe/exposome-and-heredity
<br>
<br>
**Ref:** Ray et al., "A powerful method for pleiotropic analysis under composite null hypothesis identifies novel shared loci between Type 2 Diabetes and Prostate Cancer",
PLoS Genet. 2020 Dec 8;16(12):e1009218, doi: 10.1371/journal.pgen.1009218 [Paper_link](https://pubmed.ncbi.nlm.nih.gov/33290408/)
<br>
<br>

Here is an example of running PLACO for two traits
<br>
**NOTE 1:** For running PLACO, first it is needed to extract common SNPs between two traits.
<br>
**NOTE 2:** PLACO needs **Z** column for running. So, if that column is not available in the GWAS data, a user should calculate it.
<br>
So, there are two parts: 1) finding common SNPs and calculation of Z column, 2) running PLACO function
## First Step
DEFINITION SECTION which should be changed by a user
```r
# directory in which input data exist
path_input_data_trait1 = "~/PLACO/"
path_input_data_trait2 = "~/PLACO/"

# the file names of input data
traits_1 <- "TRAIT_1.txt"
traits_2 <- "TRAIT_2.txt"

# directory in which an output data would be written
path_output_data = "~/PLACO/"

# the file names of the output data 
output_names <- c(  "TRAIT_1_Common_inc_Z",
                    "TRAIT_2_Common_inc_Z")

# used for filtering GWAS data              
info_threshold <- 0.8
MAF_threshold <- 0.01
```
libraries used in the code
```r
library(vroom)
library(dplyr)
library(data.table)
```
RUNNING SECTIONS
<br>
Summary: Reading GWAS data, filtering some SNPs, calculation of Z based on "beta" and "se", extraction of common SNPs between pair of GWAS data
```r
writeLines("\n\n")
print(traits_1)
writeLines("\n\n")
print(traits_2)
writeLines("\n\n")

# reading treat_1
gwas_traits_1 <- vroom(file=paste0(path_input_data_trait1, traits_1))
gwas_traits_1 <- as.data.frame(gwas_traits_1)

# reading treat_2
gwas_traits_2 <- vroom(file=paste0(path_input_data_trait2, traits_2))
gwas_traits_2 <- as.data.frame(gwas_traits_2)
```
**IMPORTANT NOTE**: The following columns MUST be available in the GWAS data. If not, a user could rename the columns in the GWAS data or makes changes in the script.  
- snp (RS ID)
- chr (chromosome)
- bp_hg19 (base pair position)
- Effect_A (Effect Allele)
- nonEffect_A (non-Effect Allele)
- beta (beta value)
- se (standard error)
- pval (P-value)
- info (r2)
- EAF (Effect Allele Frequency)
- MAF (Minor Allele Frequency)

Now a user could continue to run the following codes:
```r
# filtering step
gwas_traits_1 <- select(filter(gwas_traits_1, gwas_traits_1$info > info_threshold & gwas_traits_1$MAF > MAF_threshold),
                        c("snp","chr","bp_hg19","Effect_A","nonEffect_A","beta","se","pval","info","EAF","MAF"))

gwas_traits_2 <- select(filter(gwas_traits_2, gwas_traits_2$info > info_threshold & gwas_traits_2$MAF > MAF_threshold),
                        c("snp","chr","bp_hg19","Effect_A","nonEffect_A","beta","se","pval","info","EAF","MAF"))
                        
# calculation of Z
# if se==0, put a VERY LARGE VALUE for Z (here 1.0e+08)
gwas_traits_1 %>% mutate(Z = case_when(gwas_traits_1$se == 0   ~ 1.0e+08 , gwas_traits_1$se != 0 ~ (gwas_traits_1$beta/gwas_traits_1$se) )) -> gwas_traits_1

gwas_traits_2 %>% mutate(Z = case_when(gwas_traits_2$se == 0   ~ 1.0e+08 , gwas_traits_2$se != 0 ~ (gwas_traits_2$beta/gwas_traits_2$se) )) -> gwas_traits_2

# merging two GWAS data based on "snp", "chr", and "bp_hg19" columns
gwas_merge <- inner_join(gwas_traits_1,gwas_traits_2, by = c("snp", "chr", "bp_hg19"))

# checking number of total, unique, and duplicated rsids
length(gwas_merge$snp)
length(unique(gwas_merge$snp))
dim(gwas_merge[duplicated(gwas_merge$snp), ])[1]

# removing potential duplicated entries
gwas_merge <- gwas_merge[!duplicated(gwas_merge$snp), ]

# checking data after removing the duplicated entries
length(gwas_merge$snp)
length(unique(gwas_merge$snp))
dim(gwas_merge[duplicated(gwas_merge$snp), ])[1]

# extracting common SNPs for the first trait
gwas_traits_1 <- gwas_merge[, c("snp", "chr", "bp_hg19", 
                                "Effect_A.x", "nonEffect_A.x", "beta.x", "se.x", 
                                "pval.x", "info.x", "EAF.x", "MAF.x", "Z.x") ]
# rename the column names
colnames(gwas_traits_1) <- c("snp", "chr", "bp_hg19", 
                             "Effect_A", "nonEffect_A", "beta", "se", 
                             "pval", "info", "EAF", "MAF", "Z") 

# sorting based on the ID column
gwas_traits_1 <- as.data.table(gwas_traits_1)
gwas_traits_1 <- setorder(gwas_traits_1, snp)

# extracting common SNPs for the second  trait
gwas_traits_2 <- gwas_merge[, c("snp", "chr", "bp_hg19", 
                                "Effect_A.y", "nonEffect_A.y", "beta.y", "se.y", 
                                "pval.y", "info.y", "EAF.y", "MAF.y", "Z.y") ]
# rename the column names
colnames(gwas_traits_2) <- c("snp", "chr", "bp_hg19", 
                             "Effect_A", "nonEffect_A", "beta", "se", 
                             "pval", "info", "EAF", "MAF", "Z") 

# sorting based on the ID column
gwas_traits_2 <- as.data.table(gwas_traits_2)
gwas_traits_2 <- setorder(gwas_traits_2, snp)

# writing the output files
vroom_write(gwas_traits_1, file = paste0(path_output_data, output_names[1], ".txt"))
vroom_write(gwas_traits_2, file = paste0(path_output_data, output_names[2], ".txt"))

writeLines("\n\n")
print(paste0('End of Running'))
```



