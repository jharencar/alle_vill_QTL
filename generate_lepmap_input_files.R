### Generating lep-map-3 input files
### K. Uckele
### March 12, 2024

## clear workspace
rm(list = ls())

#library(readr)
library(dplyr)

setwd("/Users/juliaharencar/Documents/Github/alle_vill_QTL/aHMM/full_final_posteriors/concatenated_files_unthinned")

files <- list.files(".", pattern = "*posterior")
#files <- files[grep("HYBR", files)]

df <- c()

for (i in files) {
  id <- strsplit(i, ".", fixed = TRUE)[[1]][1]
  sample_id <- substr(id, 1, nchar(id) - 3)
  #pop <- substr(id, nchar(id) - 3, nchar(id))
  post <- read.table(i, sep = "\t", header = FALSE)
  # rename cols 
  names(post) <- c("chrom", "position", "ref_homo", "het", "alt_homo")
  # assign genetic codes
  post <- post %>% mutate(gencode = case_when(
    ref_homo > 0.9  ~ "AA" ,
    het > 0.9  ~ "AB" ,
    alt_homo > 0.9  ~ "BB",
    ref_homo <= 0.9 & het <= 0.9 & alt_homo < 0.9 ~ "NoCall"
  ))
  
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
  df <- rbind(df, c(sample_id, post$gencode))
}

rownames(df) <- NULL

dim(df)
#[1]   414 94695 - MINE
#[1]   502 16082 - KAU

colnames(df) <- df[1,]
df <- df[-1,]

View(df)


################################################################################
## create input file for lep map 3
################################################################################

df <- df[-1,]
df.t <- t(df)


rownames(df.t)[1] <- "SNP"

View(df.t)

dim(df.t)

df.t.2 <- cbind(df.t, df.t[,ncol(df.t)])

View(df.t.2)

#df.t.2[1,500] <- "F1_SIM"

write.table(df.t.2, "/Users/kathrynuckele/Documents/lepmap3/affx_genotypes.txt", quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")

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

firstline <- c("CHR", "POS", rep("F", 500))
secondline <- c("CHR", "POS", df.t.2[1,])
thirdline <- c("CHR", "POS", rep("22_194", 499), "0")
fourthline <- c("CHR", "POS", rep("F1_SIM", 499), "0")
fifthline <- c("CHR", "POS", rep("0", 498), "1", "2")
sixthline <- c("CHR", "POS", rep("0", 500))

ped <- rbind(firstline, secondline, thirdline, fourthline, fifthline, sixthline)
View(ped)

write.table(ped, "/Users/kathrynuckele/Documents/lepmap3/pedigree.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")





