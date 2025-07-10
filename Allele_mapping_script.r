#R script for Streptococcus suis allele mapping 
#This R script extracts SNPs from the following resistance genes from a core genome reference alignment to S. suis P1/7 (AM946016.1): pbp2B, dhfr, dhfr promoter, pbp2x, mraY, parC, and gyrA
#The script was written in for R v4.2.3 
#Install package ape if necessary and activate library ape
#install.packages("ape")
library(ape)

#Read in alignment mapped against P1/7, exchange file name if input alignment file is named differently
x <- read.dna ("Alignment.fna", 'fasta')

#Detailed description of pbp2b extraction
#pbp2B contains allele exchanges at four possible positions (K479T/A, D512E/Q/K/A, K513E,/D, T515S)
#Other genes follow the same procedure
#pbp2B/penA - K479T/A, D512E/Q/K/A, K513E,/D, T515S
#Extract pbp2b sequence and translate into protein sequence
pbp2b<-as.matrix(trans(complement(x[,1207241:1209316])))
#Store SNP postions in variabel called mutations_pbp2b
mutations_pbp2b<-c(479, 512, 513, 515)
#Store alleles at these SNP positions (479, 512, 513, 515) in a matrix called pbp2bm
pbp2bm<-pbp2b[,mutations_pbp2b]
#Check whether alleles at these positions are a K/D/K/T (normal phenotype) and add column names
M_pbp2b<-cbind(as.character(pbp2bm[,1])=='K', as.character(pbp2bm[,2])=='D', as.character(pbp2bm[,3])=='K', as.character(pbp2bm[,4])=='T' )
colnames(M_pbp2b) <-c("pbp2b_K479T/A_Susceptible", "pbp2b_D512E/Q/K/A_Susceptible", "pbp2b_K513E_Susceptible", "pbp2b_T515S_Susceptible")
#Check which allele is at the position of the possible SNP and add column names
M_pbp2b_2<-cbind(as.character(pbp2bm[,1]), as.character(pbp2bm[,2]), as.character(pbp2bm[,3]), as.character(pbp2bm[,4]))
colnames(M_pbp2b_2) <-c("pbp2b_K479T/A_Allele", "pbp2b_D512E/Q/K/A_Allele", "pbp2b_K513E_Allele", "pbp2b_T515S_Allele")
#Combine all results into matrix called M_pbp2b_All and add column names
M_pbp2b_All <- cbind(M_pbp2b[,1], M_pbp2b_2[,1], M_pbp2b[,2], M_pbp2b_2[,2], M_pbp2b[,3], M_pbp2b_2[,3], M_pbp2b[,4], M_pbp2b_2[,4])
colnames(M_pbp2b_All) <-c("pbp2b_K479T/A_Susceptible", "pbp2b_K479T/A_Allele", "pbp2b_D512E/Q/K/A_Susceptible", "pbp2b_D512E/Q/K/A_Allele", "pbp2b_K513E", "pbp2b_K513E_Allele", "pbp2b_T515S_Susceptible", "pbp2b_T515S_Allele")
#Report whether the alleles relevant for the resistance are present and add column names
M_pbp2b_3<-cbind(grepl("T|A",(pbp2bm[,1])), grepl("E|Q|K|A",(pbp2bm[,2])), grepl("E|D",(pbp2bm[,3])), grepl("S",(pbp2bm[,4])))
colnames(M_pbp2b_3) <-c("pbp2b_K479T/A_Resistant", "pbp2b_D512E/Q/K/A_Resistant", "pbp2b_K513E_Resistant", "pbp2b_T515S_Resistant")
#Save protein alignment of pbp2b/penA gene
write.dna(pbp2b, "pbp2b.fas", 'fasta')

#dhfr - I102L
dhfr<-as.matrix(trans(x[,813422:813934]))
mutations_dhfr<-(102)
dhfrm<-dhfr[,mutations_dhfr]
#Check whether allele at this position is an I 
M_dhfr<-as.character(dhfrm[,1])=='I'
colnames(M_dhfr) <-c("dhfr_I102L_Susceptible")
#Check which allele is at the position of possible SNP
M_dhfr_2<-as.character(dhfrm[,1])
colnames(M_dhfr_2) <-c("dhfr_I102L_Allele")
M_dhfr_All <- cbind(M_dhfr[,1], M_dhfr_2[,1])
colnames(M_dhfr_All) <-c("dhfr_I102L_Susceptible", "dhfr_I102L_Allele")
#Report whether the allele relevant for the resistance is present 
M_dhfr_3<-as.character(dhfrm)=='L'
colnames(M_dhfr_3) <-c("dhfr_I102L_Resistant")
#Save protein alignment of dhfr gene
write.dna(dhfr, "dhfr.fas", 'fasta')

#dhfr promoter - A5G
dhfr_promoter<-as.matrix(x)[,813392:813421]
mutations_dhfr_promoter<-(5)
dhfr_promoterm<-dhfr_promoter[,mutations_dhfr_promoter]
#Check whether allele at this position is an a
M_dhfr_prom<-as.character(dhfr_promoterm[,1])=='a'
colnames(M_dhfr_prom) <-c("dhfr_prom_A5G_Susceptible")
#Check which allele is at the position of the possible SNP
M_dhfr_prom_2<-as.character(dhfr_promoterm[,1])
colnames(M_dhfr_prom_2) <-c("dhfr_prom_A5G_Allele")
M_dhfr_prom_All <- cbind(M_dhfr_prom[,1], M_dhfr_prom_2[,1])
colnames(M_dhfr_prom_All) <-c("dhfr_prom_A5G_Susceptible", "dhfr_prom_A5G_Allele")
#Report whether the alleles relevant for the resistance is present 
M_dhfr_prom_3<-as.character(dhfr_promoterm)=='g'
colnames(M_dhfr_prom_3) <-c("dhfr_prom_A5G_Resistant")
#Save nucleotide alignment of dhfr promoter
write.dna(dhfr_promoter, "dhfr_promoter.fas", 'fasta')

#pbp2X (pbpX) - M437L_C, S445T_C, T467S_C, Y525F_C, T551S_P1, L594Y/F_P2, V596G_P2, N569Q_P3
pbp2X<-as.matrix(trans(complement(x[,1556162:1558411])))
mutations_pbp2X<-c(437, 445, 467, 525, 551, 594, 596, 569)
pbp2Xm<-pbp2X[,mutations_pbp2X]
#Check whether alleles at these positions are an M/S/T/Y/T/L/V/N
M_pbp2X<-cbind(as.character(pbp2Xm[,1])=='M', as.character(pbp2Xm[,2])=='S', as.character(pbp2Xm[,3])=='T', as.character(pbp2Xm[,4])=='Y', as.character(pbp2Xm[,5])=='T', as.character(pbp2Xm[,6])=='L', as.character(pbp2Xm[,7])=='V', as.character(pbp2Xm[,8])=='N')
colnames(M_pbp2X) <-c("pbp2X_M437L_C_Susceptible", "pbp2X_S445T_C_Susceptible", "pbp2X_T467S_C_Susceptible", "pbp2X_Y525F_C_Susceptible", "pbp2X_T551S_P1_Susceptible", "pbp2X_L594Y/F_P2_Susceptible", "pbp2X_V596G_P2_Susceptible", "pbp2X_N569Q_P3_Susceptible")
#Check which allele is at the position of possible SNP
M_pbp2X_2<-cbind(as.character(pbp2Xm[,1]), as.character(pbp2Xm[,2]), as.character(pbp2Xm[,3]), as.character(pbp2Xm[,4]), as.character(pbp2Xm[,5]), as.character(pbp2Xm[,6]), as.character(pbp2Xm[,7]), as.character(pbp2Xm[,8]))
colnames(M_pbp2X_2) <-c("pbp2X_M437L_C_Allele", "pbp2X_S445T_C_Allele", "pbp2X_T467S_C_Allele", "pbp2X_Y525F_C_Allele", "pbp2X_T551S_P1_Allele", "pbp2X_L594Y/F_P2_Allele", "pbp2X_V596G_P2_Allele", "pbp2X_N569Q_P3_Allele")
M_pbp2X_All <- cbind(M_pbp2X[,1], M_pbp2X_2[,1],M_pbp2X[,2], M_pbp2X_2[,2],M_pbp2X[,3], M_pbp2X_2[,3], M_pbp2X[,4], M_pbp2X_2[,4],M_pbp2X[,5], M_pbp2X_2[,5],M_pbp2X[,6], M_pbp2X_2[,6], M_pbp2X[,7], M_pbp2X_2[,7],M_pbp2X[,8], M_pbp2X_2[,8])
colnames(M_pbp2X_All) <-c("pbp2X_M437L_C_Susceptible", "pbp2X_M437L_C_Allele", "pbp2X_S445T_C_Susceptible", "pbp2X_S445T_C_Allele", "pbp2X_T467S_C_Susceptible", "pbp2X_T467S_C_Allele", "pbp2X_Y525F_C_Susceptible", "pbp2X_Y525F_C_Allele", "pbp2X_T551S_P1_Susceptible", "pbp2X_T551_S_P1_Allele", "pbp2X_L594Y/F_P2_Susceptible", "pbp2X_L594Y/F_P2_Allele", "pbp2X_V596G_P2_Susceptible", "pbp2X_V596G_P2_Allele", "pbp2X_N569Q_P3_Susceptible", "pbp2X_N569Q_P3_Allele")
#Report whether the alleles relevant for the resistance are present 
M_pbp2X_3<-cbind(grepl("L|C",(pbp2Xm[,1])), grepl("T|C",(pbp2Xm[,2])), grepl("S|C",(pbp2Xm[,3])), grepl("F|C",(pbp2Xm[,4])), grepl("S",(pbp2Xm[,5])), grepl("Y|F",(pbp2Xm[,6])), grepl("G",(pbp2Xm[,7])), grepl("Q",(pbp2Xm[,8])))
colnames(M_pbp2X_3) <-c("pbp2X_M437L_C_Resistant", "pbp2X_S445T_C_Resistant", "pbp2X_T467S_C_Resistant", "pbp2X_Y525F_C_Resistant", "pbp2X_T551S_P1_Resistant", "pbp2X_L594Y/F_P2_Resistant", "pbp2X_V596G_P2_Resistant", "pbp2X_N569Q_P3_Resistant")
#Save protein alignment of pbp2x gene
write.dna(pbp2X, "pbp2X.fas", 'fasta')

#mraY - M6I/L and/or A4S/T or G8S
mraY<-as.matrix(trans(complement(x[,1555162:1556160])))
mutations_mraY<-c(6, 4, 8)
mraYm<-mraY[,mutations_mraY]
#Check whether alleles at these positions are an M/A/G
M_mraY<-cbind(as.character(mraYm[,1])=='M', as.character(mraYm[,2])=='A', as.character(mraYm[,3])=='G')
colnames(M_mraY) <-c("mraY_M6I/L_Susceptible", "mraY_A4S/T_Susceptible", "mraY_G8S_Susceptible")
#Check which allele is at the position of possible SNP
M_mraY_2<-cbind(as.character(mraYm[,1]), as.character(mraYm[,2]), as.character(mraYm[,3]))
colnames(M_mraY_2) <-c("mraY_M6I/L_Allele", "mraY_A4S/T_Allele", "mraY_G8S_Allele")
M_mraY_All <- cbind(M_mraY[,1], M_mraY_2[,1],M_mraY[,2], M_mraY_2[,2],M_mraY[,3], M_mraY_2[,3])
colnames(M_mraY_All) <-c("mraY_M6I/L_Susceptible", "mraY_M6I/L_Allele", "mraY_A4S/T_Susceptible", "mraY_A4S/T_Allele", "mraY_G8S_Susceptible", "mraY_G8S_Allele")
#Report whether the alleles relevant for the resistance are present 
M_mraY_3<-cbind(grepl("I|L",(mraYm[,1])), grepl("S|T",(mraYm[,2])), grepl("SL",(mraYm[,3])))
colnames(M_mraY_3) <-c("mraY_M6I/L_Resistant", "mraY_A4S/T_Resistant", "mraY_G8S_Resistant")
#Save protein alignment of mraY gene
write.dna(mraY, "mraY.fas", 'fasta')

#parC - S79
parC<-as.matrix(trans(x[,749400:751817]))
mutations_parC<-(79)
parCm<-parC[,mutations_parC]
#Check whether allele at this position is an S
M_parC<-as.character(parCm[,1])=='S'
colnames(M_parC) <-c("parC_S79_Susceptible")
#Check which allele is at the position of possible SNP
M_parC_2<-as.character(parCm[,1])
colnames(M_parC_2) <-c("parC_S79_Allele")
M_parC_All <- cbind(M_parC[,1], M_parC_2[,1])
colnames(M_parC_All) <-c("parC_S79_Susceptible", "parC_S79_Allele")
#Report whether the alleles relevant for the resistance is present 
M_parC_3<-as.character(parCm)!='S'
colnames(M_parC_3) <-c("parC_S79_Resistant")
#Save protein alignment of parC gene
write.dna(parC, "parC.fas", 'fasta')

#gyrA - S81, E85
gyrA<-as.matrix(trans(complement(x[,953534:955978])))
mutations_gyrA<-c(81, 85)
gyrAm<-gyrA[,mutations_gyrA]
#Check whether alleles at these positions are an S/E
M_gyrA<-cbind(as.character(gyrAm[,1])=='S', as.character(gyrAm[,2])=='E')
colnames(M_gyrA) <-c("gyrA_S81_Susceptible", "gyrA_E85_Susceptible")
#Check which allele is at the position of possible SNP
M_gyrA_2<-cbind(as.character(gyrAm[,1]), as.character(gyrAm[,2]))
colnames(M_gyrA_2) <-c("gyrA_S81_Allele", "gyrA_E85_Allele")
M_gyrA_All <- cbind(M_gyrA[,1], M_gyrA_2[,1], M_gyrA[,2], M_gyrA_2[,2])
colnames(M_gyrA_All) <-c("gyrA_S81_Susceptible", "gyrA_S81_Allele", "gyrA_E85_Susceptible", "gyrA_E85_Allele")
#Report whether the alleles relevant for the resistance are present 
M_gyrA_3<-cbind(as.character(gyrAm[,1])!='S', as.character(gyrAm[,2])!='E')
colnames(M_gyrA_3) <-c("gyrA_S81_Resistant", "gyrA_E85_Resistant")
#Save protein alignment of gyrA gene
write.dna(gyrA, "gyrA.fas", 'fasta')

#Combine all the results in one table, includes resistance phenotype and columns indicating the present alleles
Complete<- cbind (M_pbp2b_All, M_pbp2b_3, M_dhfr_All, M_dhfr_3, M_dhfr_prom_All, M_dhfr_prom_3, M_pbp2X_All, M_pbp2X_3, M_mraY_All, M_mraY_3, M_parC_All, M_parC_3, M_gyrA_All, M_gyrA_3)
write.csv (Complete,file="Complete_Allele_mapping.csv")

#Combine only resistance results in one table, True indicates presence of the allele, False absence
Complete_1<- cbind (M_pbp2b_3, M_dhfr_3, M_dhfr_prom_3, M_pbp2X_3, M_mraY_3, M_parC_3, M_gyrA_3)
write.csv (Complete_1,file="Complete_Allele_mapping_Resistance.csv")