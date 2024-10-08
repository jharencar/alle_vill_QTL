### Generating lep-map-3 input files
### K. Uckele
### March 12, 2024

## clear workspace
rm(list = ls())

#library(readr)
library(dplyr)

setwd("/Users/juliaharencar/Documents/Github/alle_vill_QTL/aHMM/full_final_posteriors/concatenated_files_thinned_50000")

files <- list.files(".", pattern = "*posterior")
#files <- c(files[grep("HYBR", files)], "19_117_LASI.posterior", "19_211_BRAC.posterior")

df <- c()

for (i in files) {
  id <- strsplit(i, ".", fixed = TRUE)[[1]][1]
  sample_id <- substr(id, 1, nchar(id) - 3)
  #pop <- substr(id, nchar(id) - 3, nchar(id))
  post <- read.table(i, sep = "\t", header = FALSE)
  # rename cols 
  names(post) <- c("chrom", "position", "ref_homo", "het", "alt_homo")
  # assign genetic codes
  post <- post %>% 
    mutate(gencode = case_when(
      ref_homo > 0.9  ~ "A" ,
      het > 0.9  ~ "H" ,
      alt_homo > 0.9  ~ "B",
      ref_homo <= 0.9 & het <= 0.9 & alt_homo < 0.9 ~ "-"),
      dist_cM = position / 1e6
    )
  
  if(exists("head1") == "FALSE") {
    head1 <- c("ID", paste0(post$chrom, "_", post$position))
    df <- rbind(df, head1)
  }
  #if(exists("head2") == "FALSE") {
  #  head2 <- c("", rep(1, nrow(post)))
  #  df <- rbind(df, head2)
  #}
  if(exists("head2") == "FALSE") {
    head2 <- c("", substr(post$chrom, nchar(post$chrom), nchar(post$chrom)))
    df <- rbind(df, head2)
  }
  if(exists("head3") == "FALSE") {
    head3 <- c("", post$dist_cM)
    df <- rbind(df, head3)
  }
  df <- rbind(df, c(sample_id, post$gencode))
}

rownames(df) <- NULL

dim(df)
#KAU   502 16082
#JGH   

colnames(df) <- df[1,]
df <- df[-1,]

View(df)


write.csv(df, "/Users/juliaharencar/Documents/Github/alle_vill_QTL/f3gen.csv", row.names = FALSE, quote = FALSE)

#write.csv(df[,1:100], "~/Documents/R:QTL/F2gen_subset.csv", row.names = FALSE, quote = FALSE)

## Test to see if there are any errors when loading into R/qtl
library(qtl)
cross <- read.cross("csv", dir="/Users/juliaharencar/Documents/Github/alle_vill_QTL", genfile = "preClean_F3gen.csv", phefile = "F2phe.csv", estimate.map=FALSE)



################################################################################
## create input file for lep map 3
################################################################################

df <- df[-1,]
df.t <- t(df)

#View(df.t)

rownames(df.t)[1] <- "SNP"

write.table(df.t, "/Users/juliaharencar/Documents/Github/alle_vill_QTL/lepmap3/affx_genotypes_thinned_50k_240316.txt", quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")
df.t <- read.table("/Users/juliaharencar/Documents/Github/alle_vill_QTL/lepmap3/affx_genotypes_thinned_50k.txt")
#################################################################################
## create pedigree file for lep map 3
#################################################################################
## dimensions should be: 6 rows, 499+2 columns

#CHR POS F   F   F   F   F   F
#CHR POS female  male    progeny_1   progeny_2   progeny_3   progeny_4
#CHR POS 0   0   male    male    male    male
#CHR POS 0   0   female  female  female  female
#CHR POS 2   1   0   0   0   0
#CHR POS 0   0   0   0   0   0

#First line is the family name
# second line is individual name
# third and fourth are the father and mother
# Line 5 contains the sex of each individual (1 male, 2 female, 0 unknown)
# sixth line is the phenotype (can be 0 for all individuals, this is not currently used). 

firstline <- c("CHR", "POS", rep("F", 420))
secondline <- c("CHR", "POS", df.t[1,-1])
thirdline <- c("CHR", "POS","0","0", "19_467","19_467", rep("19_615M", 26), "21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_238M","21_238M","21_232M","21_238M","21_238M","21_238M","21_232M","21_238M","21_238M","21_238M","21_238M","21_232M","21_238M","21_238M","21_232M","21_238M","21_232M","21_232M","21_238M","21_232M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_232M","21_238M","21_238M","21_259M","21_259M","21_284M","21_284M","21_259M","21_284M","21_259M","21_259M","21_259M","21_259M","21_259M","21_284M","21_284M","21_284M","21_259M","21_259M","21_284M","21_259M","21_259M","21_259M","21_259M","21_284M","21_259M","21_284M","21_284M","21_259M","21_284M","21_284M","21_284M","21_284M","21_259M","21_259M","21_284M","21_284M","21_259M","21_284M","21_259M","21_259M","21_259M","21_259M","21_284M","21_259M","21_259M","21_259M","21_284M","21_284M","21_259M","21_284M","21_284M","21_284M","21_284M","21_284M","21_284M","21_232M","21_232M","21_232M","21_238M","21_232M","21_238M","21_238M","21_238M","21_238M","21_232M","21_238M","21_238M","21_232M","21_238M","21_232M","21_238M","21_238M","21_238M","21_259M","21_238M","21_259M","21_238M","21_238M","21_238M","21_238M","21_238M","21_284M","21_232M","21_232M","21_232M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_232M","21_284M","21_284M","21_259M","21_259M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_284M","21_284M","21_232M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_232M","21_232M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_284M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_232M","21_232M","21_238M","21_238M","21_238M","21_238M","21_232M","21_241M","21_241M","21_284M","21_241M","21_232M","21_241M","21_241M","21_241M","21_238M","21_238M","21_241M","21_238M","21_232M","21_238M","21_238M","21_241M","21_241M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_232M","21_232M","21_232M","21_232M","21_238M","21_232M","21_238M","21_238M","21_238M","21_232M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_232M","21_238M","21_241M","21_259M","21_238M","21_232M","21_232M","21_238M","21_238M","21_259M","21_232M","21_232M","21_232M","21_232M","21_238M","21_232M","21_284M") # not actually sure if this is what I should do here
fourthline <- c("CHR", "POS", "0", "0", "19_500", "19_500", rep("19_615F", 26),"21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_238F","21_238F","21_232F","21_238F","21_238F","21_238F","21_232F","21_238F","21_238F","21_238F","21_238F","21_232F","21_238F","21_238F","21_232F","21_238F","21_232F","21_232F","21_238F","21_232F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_232F","21_238F","21_238F","21_259F","21_259F","21_284F","21_284F","21_259F","21_284F","21_259F","21_259F","21_259F","21_259F","21_259F","21_284F","21_284F","21_284F","21_259F","21_259F","21_284F","21_259F","21_259F","21_259F","21_259F","21_284F","21_259F","21_284F","21_284F","21_259F","21_284F","21_284F","21_284F","21_284F","21_259F","21_259F","21_284F","21_284F","21_259F","21_284F","21_259F","21_259F","21_259F","21_259F","21_284F","21_259F","21_259F","21_259F","21_284F","21_284F","21_259F","21_284F","21_284F","21_284F","21_284F","21_284F","21_284F","21_232F","21_232F","21_232F","21_238F","21_232F","21_238F","21_238F","21_238F","21_238F","21_232F","21_238F","21_238F","21_232F","21_238F","21_232F","21_238F","21_238F","21_238F","21_259F","21_238F","21_259F","21_238F","21_238F","21_238F","21_238F","21_238F","21_284F","21_232F","21_232F","21_232F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_232F","21_284F","21_284F","21_259F","21_259F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_284F","21_284F","21_232F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_232F","21_232F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_284F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_232F","21_232F","21_238F","21_238F","21_238F","21_238F","21_232F","21_241F","21_241F","21_284F","21_241F","21_232F","21_241F","21_241F","21_241F","21_238F","21_238F","21_241F","21_238F","21_232F","21_238F","21_238F","21_241F","21_241F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_232F","21_232F","21_232F","21_232F","21_238F","21_232F","21_238F","21_238F","21_238F","21_232F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_232F","21_238F","21_241F","21_259F","21_238F","21_232F","21_232F","21_238F","21_238F","21_259F","21_232F","21_232F","21_232F","21_232F","21_238F","21_232F","21_284F") 
fifthline <- c("CHR", "POS", "1", "2", "1", "2", rep("0", 3), "1", "2", rep("0", 3), "1", "2", rep("0", 2), "1", "2", "0", "1","2", rep("0", 7), "1", "2", rep("0", 390)) # Is this correct? 1 = vill; 2 = alle (best guess is that alle was mom to wild, selfed F1)
sixthline <- c("CHR", "POS", rep("0", 420))

ped <- rbind(firstline, secondline, thirdline, fourthline, fifthline, sixthline)
View(ped)

#write.table(ped, "/Users/juliaharencar/Documents/Github/alle_vill_QTL/lepmap3/pedigree.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

genotypes <- read.table("/Users/juliaharencar/Documents/Github/alle_vill_QTL/lepmap3/affx_genotypes_thinned_50k.txt", header = TRUE)
View(genotypes)
nrow(genotypes)
ch1 <- sum(substr(genotypes$SNP, 6, 6) == 1)
#[1] 2195 JGH 2167
ch2 <- sum(substr(genotypes$SNP, 6, 6) == 2)
#[1] 2381 JGH 2318
ch3 <- sum(substr(genotypes$SNP, 6, 6) == 3)
#[1] 1766 JGH 1736
ch4 <- sum(substr(genotypes$SNP, 6, 6) == 4)
#[1] 1588 JGH 1509
ch5 <- sum(substr(genotypes$SNP, 6, 6) == 5)
# 2080 JGH 2113
ch6 <- sum(substr(genotypes$SNP, 6, 6) == 6)
# 1702 JGH 1699
ch7 <- sum(substr(genotypes$SNP, 6, 6) == 7)
# 1419 JGH 424
ch8 <- sum(substr(genotypes$SNP, 6, 6) == 8)
# 1564 JGH 1567
ch9 <- sum(substr(genotypes$SNP, 6, 6) == 9)
# 1586 JGH 1548

linkage_groups <- substr(genotypes$SNP, 6, 6)

write.table(linkage_groups, "/Users/juliaharencar/Documents/Github/alle_vill_QTL/lepmap3/manual_map1.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

map1 <- read.table("/Users/juliaharencar/Documents/Github/alle_vill_QTL/lepmap3/map1.txt")
head(map1)
unique(map1$V1)
