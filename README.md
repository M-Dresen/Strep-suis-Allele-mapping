Read.me R script for Streptococcus suis allele mapping 

Credits:
R: Core Team R. R: A language and environment for statistical computing [Software]. Vienna, Austria: R Foundation for Statistical Computing. Retrieved from https://www.R-project.org/; 2018.
Ape: Paradis E, Schliep K. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics. 2019;35(3):526-8.

This R script extracts SNPs from the following resistance genes from core genome reference alignment to P1/7: pbp2B, dhfr, dhfr promoter, pbp2x, mraY, parC, and gyrA
The script was written in R v4.2.3

Required packages: library(ape) 
Install package with install.packages("ape") and activate with library(ape) 

Set working directory to directory containing P1/7 reference mapped alignment and R script
-	Either rename alignment to “Alignment.fna” or change the alignment name in the R script

Commands in the manuscript (example dhfr – I102L)
-	dhfr<-as.matrix(trans(x[,813422:813934]))
	#Extracts and translates defined region of the gene
-	(trans(complement(x[…]))
	##Extracts and translates defined region of gene on complementary string)
-	mutations_dhfr<-(102)
	#Stores position (102) of SNP in the variabel mutations_dhfr
-	dhfrm<-dhfr[,mutations_dhfr]
	#Extracts amino acid at this position and stores it in a matrix called dhfrm
-	M_dhfr<- as.character(dhfrm[,1])=='I' 
	#Screens position one of matrix dhfrm for the amino acid indicated (I) and reports either true or false
-	colnames(M_dhfr) <-c("dhfr_I102L_Susceptible")
	Adds column name to M_dhfr
-	M_dhfr_2 <- as.character(dhfrm[,1]) 
	#Checks which allele is at the indicated position in the matrix 
-	colnames(M_dhfr_2) <-c("dhfr_I102L_Allele")
 	#Adds column name to matrix M_dhfr_2
-	M_dhfr_All <- cbind(M_dhfr[,1], M_dhfr_2[,1])
	#Combine results of M_dhfr and M_dhfr_2 into one matrix (True/False + present allele)
-	colnames(M_dhfr_All) <-c("dhfr_I102L_Susceptible", "dhfr_I102L_Allele")
	#Adds column name to matrix M_dhfr_All
-	M_dhfr_3<-as.character(dhfrm)=='L' 
	#Reports whether the mutation relevant for the resistance is present
-	colnames(M_dhfr_3) <-c("dhfr_I102L_Resistant") 
	#Adds column name to  matrix M_dhfr_3
-	write.dna(dhfr, "dhfr.fas", 'fasta') 
	#Save protein alignment of dhfr

Summary and export of results
-	Complete<- cbind (M_Pbp2b_All, M_Pbp2b_3, M_dhfr_All, M_dhfr_3, M_dhfr_prom_All, M_dhfr_prom_3, M_Pbp2X_All, M_Pbp2X_3, M_mraY_All, M_mraY_3, M_parC_All, M_parC_3, M_gyrA_All, M_gyrA_3)
	#Combines all results into one dataframe, includes resistance phenotype and columns indicating the present alleles
-	write.csv (Complete,file="Complete_Allele_mapping.csv")
	#Saves the combined dataframe
-	Complete_1<- cbind (M_Pbp2b_3, M_dhfr_3, M_dhfr_prom_3, M_Pbp2X_3, M_mraY_3, M_parC_3, M_gyrA_3)
	#Only combines results for the resistance testing, reports true when needed resistance allele is present and false when allele is absent
-	write.csv (Complete_1,file="Complete_Allele_mapping_Resistance.csv")
	Saves the file to a csv file

Additional note
-	dhfr promoter: base exchange not allele exchange, no translation of gene region necessary, bases are indicated with small letters
