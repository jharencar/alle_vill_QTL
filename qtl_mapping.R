## We used the scanone function of R/qtl to perform single QTL model standard 
## interval mapping using the EM algorithm.33 Recombination fraction was estimated 
## using the est.rf() function and markers missing genotype data were excluded 
## using the drop.nullmarker() function

setwd("/Users/juliaharencar/Documents/Github/alle_vill_QTL/")

## clear workspace
rm(list = ls())

library(qtl)
#library(snow)
#library(qtl2)

########### F2s ###########
################################################################################
## Load the data
## Note the use of estimate.map=FALSE. Markers will be assigned dummy locations, 
## with no attempt to estimate inter-marker distances.

CostusF2 <- read.cross("csvs", ".", genfile = "F2gen.csv", phefile = "F3phe.csv", estimate.map=FALSE)

## Omit individuals and markers with lots of missing data
summary(CostusF2)

## look at the pattern of missing data (can take a few minutes to load)
#plotMissing(Costus)

## Plot of number of genotyped markers for each individual (left panel) and 
## number of genotyped individuals for each marker (right panel).
# par(mfrow=c(1,2), las=1)
# plot(ntyped(Costus), ylab="No. typed markers", main="No. genotypes by individual") 
# plot(ntyped(Costus, "mar"), ylab="No. typed individuals", main="No. genotypes by marker")

## Omit the individuals with lots of missing genotype data
CostusF2 <- subset(CostusF2, ind=(ntyped(CostusF2)>15000))
# brings total individuals 391 --> 387

## identify names of markers to drop, based on genotypes-per-marker plot
# drop all markers with less than ~80% genotyped individuals (21*0.8 = 17)
nt.bymar <- ntyped(CostusF2, "mar")
todrop <- names(nt.bymar[nt.bymar < 17]) 
## drop those markers
CostusF2 <- drop.markers(CostusF2, todrop)

# summary(Costus) # brings total markers 16281 --> 16201, % genotyped from 98.3 --> 98.8%

###### REMOVE MARKERS WITH SIGNAL OF DISTORTED SEGREGATION ######

## Look for markers with distorted segregation patterns
## We expect the genotypes to appear with the frequencies 1:2:1. 
## Moderate departures from these frequencies are not unusual and may indicate 
## the presence of partially lethal alleles. Gross departures from these 
## frequencies often indicate problematic markers that should be omitted, 
## at least initially

gt <- geno.table(CostusF2)
gt_sig <- gt[gt$P.value < 0.05/totmar(CostusF2),]

## histogram of p-values
hist(gt_sig$P.value, breaks = seq(0,3.0e-6, 0.00000001))
gt_sig_low_sub <- gt_sig %>% 
  filter(P.value < 2.0e-10)
hist(gt_sig_low_sub$P.value, breaks = seq(0,2.0e-11, 1.0e-20))

## omit the worst of these markers
todrop <- rownames(gt[gt$P.value < 1e-10,])
Costus <- drop.markers(Costus, todrop)

###### PLOT GENOTYPE FREQUENCIES ######

## pull out average genotype frequency across markers by individual
g <- pull.geno(Costus)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))

## plot avg individual genotype frequency for each genotype
par(mfrow=c(1,3), las=1)
for(i in 1:3)
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i], ylim=c(0,1))


########### F3s ###########
################################################################################
## Load the data
## Note the use of estimate.map=FALSE. Markers will be assigned dummy locations, 
## with no attempt to estimate inter-marker distances.

Costus <- read.cross("csvs", ".", genfile = "F3gen.csv", phefile = "F3phe_onehotcovar.csv", estimate.map=FALSE)

#Warning message:
#  In summary.cross(cross) :
#  Some chromosomes > 1000 cM in length; there may be a problem with the genetic map.
#(Perhaps it is in basepairs?)

## Omit individuals and markers with lots of missing data
# summary(Costus)

# F2 intercross
# 
# No. individuals:    391 
# 
# No. phenotypes:     8 
# Percent phenotyped: 100 98.2 96.7 87.2 88.5 88.5 88.7 97.2 
# 
# No. chromosomes:    9 
# Autosomes:      1 2 3 4 5 6 7 8 9 
# 
# Total markers:      16281 
# No. markers:        2195 2381 1766 1588 2080 1702 1419 1564 1586 
# Percent genotyped:  98.3 
# Genotypes (%):      AA:25.4  AB:26.7  BB:47.9  not BB:0.0  not AA:0.0 

## look at the pattern of missing data (can take a few minutes to load)
#plotMissing(Costus)

## Plot of number of genotyped markers for each individual (left panel) and 
## number of genotyped individuals for each marker (right panel).
# par(mfrow=c(1,2), las=1)
# plot(ntyped(Costus), ylab="No. typed markers", main="No. genotypes by individual") 
# plot(ntyped(Costus, "mar"), ylab="No. typed individuals", main="No. genotypes by marker")

## Omit the individuals with lots of missing genotype data
Costus <- subset(Costus, ind=(ntyped(Costus)>15000))
# brings total individuals 391 --> 387

## identify names of markers to drop, based on genotypes-per-marker plot
# drop all markers with less than ~80% genotyped individuals (391*0.8 = 313)
nt.bymar <- ntyped(Costus, "mar")
todrop <- names(nt.bymar[nt.bymar < 313]) 
## drop those markers
Costus <- drop.markers(Costus, todrop)

# summary(Costus) # brings total markers 16281 --> 16201, % genotyped from 98.3 --> 98.8%

## Identify duplicate individuals
# cg <- comparegeno(Costus)
# par(mfrow=c(1,1), las=1)
# hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
# rug(cg[lower.tri(cg)])

## identify pairs with over 90% matching genotypes
# wh <- which(cg > 0.9, arr=TRUE)
# wh <- wh[wh[,1] < wh[,2],]
# wh
# 
# ## inspect the genotype matches for these pairs
# g1 <- pull.geno(Costus, chr = 1)
# g2 <- pull.geno(Costus, chr = 2)
# g3 <- pull.geno(Costus, chr = 3)
# g4 <- pull.geno(Costus, chr = 4)
# g5 <- pull.geno(Costus, chr = 5)
# g6 <- pull.geno(Costus, chr = 6)
# g7 <- pull.geno(Costus, chr = 7)
# g8 <- pull.geno(Costus, chr = 8)
# g9 <- pull.geno(Costus, chr = 9)
# 
# 
# ## Omit the genotypes that mismatch, as these are indicated to be errors in one 
# ## or the other individual (or both)
# for(i in 1:nrow(wh)) {
#   tozero <- !is.na(g1[wh[i,1],]) & !is.na(g1[wh[i,2],]) & g1[wh[i,1],] != g1[wh[i,2],]
#   Costus$geno[[1]]$data[wh[i,1],tozero] <- NA 
# }
# for(i in 1:nrow(wh)) {
#   tozero <- !is.na(g2[wh[i,1],]) & !is.na(g2[wh[i,2],]) & g2[wh[i,1],] != g2[wh[i,2],]
#   Costus$geno[[2]]$data[wh[i,1],tozero] <- NA 
# }
# for(i in 1:nrow(wh)) {
#   tozero <- !is.na(g3[wh[i,1],]) & !is.na(g3[wh[i,2],]) & g3[wh[i,1],] != g3[wh[i,2],]
#   Costus$geno[[3]]$data[wh[i,1],tozero] <- NA 
# }
# for(i in 1:nrow(wh)) {
#   tozero <- !is.na(g4[wh[i,1],]) & !is.na(g4[wh[i,2],]) & g4[wh[i,1],] != g4[wh[i,2],]
#   Costus$geno[[4]]$data[wh[i,1],tozero] <- NA 
# }
# for(i in 1:nrow(wh)) {
#   tozero <- !is.na(g5[wh[i,1],]) & !is.na(g5[wh[i,2],]) & g5[wh[i,1],] != g5[wh[i,2],]
#   Costus$geno[[5]]$data[wh[i,1],tozero] <- NA 
# }
# for(i in 1:nrow(wh)) {
#   tozero <- !is.na(g6[wh[i,1],]) & !is.na(g6[wh[i,2],]) & g6[wh[i,1],] != g6[wh[i,2],]
#   Costus$geno[[6]]$data[wh[i,1],tozero] <- NA 
# }
# for(i in 1:nrow(wh)) {
#   tozero <- !is.na(g7[wh[i,1],]) & !is.na(g7[wh[i,2],]) & g7[wh[i,1],] != g7[wh[i,2],]
#   Costus$geno[[7]]$data[wh[i,1],tozero] <- NA 
# }
# for(i in 1:nrow(wh)) {
#   tozero <- !is.na(g8[wh[i,1],]) & !is.na(g8[wh[i,2],]) & g8[wh[i,1],] != g8[wh[i,2],]
#   Costus$geno[[8]]$data[wh[i,1],tozero] <- NA 
# }
# for(i in 1:nrow(wh)) {
#   tozero <- !is.na(g9[wh[i,1],]) & !is.na(g9[wh[i,2],]) & g9[wh[i,1],] != g9[wh[i,2],]
#   Costus$geno[[9]]$data[wh[i,1],tozero] <- NA 
# }

# ## Omit one individual from each pair
# Costus <- subset(Costus, ind=-wh[,2])
# 
# ## look for duplicate markers (i.e., markers with identical genotypes)
# dup <- findDupMarkers(Costus, exact.only=FALSE, adjacent.only=TRUE)
# Costus <- drop.markers(Costus, unlist(dup))

###### REMOVE MARKERS WITH SIGNAL OF DISTORTED SEGREGATION ######

## Look for markers with distorted segregation patterns
## We expect the genotypes to appear with the frequencies 1:2:1. 
## Moderate departures from these frequencies are not unusual and may indicate 
## the presence of partially lethal alleles. Gross departures from these 
## frequencies often indicate problematic markers that should be omitted, 
## at least initially

# gt <- geno.table(Costus)
# gt_sig <- gt[gt$P.value < 0.05/totmar(Costus),]
# 
# ## histogram of p-values
# hist(gt_sig$P.value, breaks = seq(0,3.0e-6, 0.00000001))
# gt_sig_low_sub <- gt_sig %>% 
#   filter(P.value < 2.0e-10)
# hist(gt_sig_low_sub$P.value, breaks = seq(0,2.0e-11, 1.0e-20))
# 
# ## omit the worst of these markers
# todrop <- rownames(gt[gt$P.value < 1e-10,])
# Costus <- drop.markers(Costus, todrop)

###### PLOT GENOTYPE FREQUENCIES ######

## pull out average genotype frequency across markers by individual
g <- pull.geno(Costus)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))

## plot avg individual genotype frequency for each genotype
par(mfrow=c(1,3), las=1)
for(i in 1:3)
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i], ylim=c(0,1))

###### ESTIMATE RECOMBINATION ######

## Study pairwise marker linkages; look for switched alleles
## The function est.rf() is used to estimate the recombination fraction between 
## each pair and to calculate a LOD score for a test of rec. freq. = 1/2.
## markerlrt() behaves just like est.rf(), but uses a general likelihood ratio 
## test in place of the usual test of pairwise linkage.
# Note: this operation requires high memory and time allocation, recommend running on a cluster/server,
# after figuring out the above parameters

Costus <- est.rf(Costus)

## BELOW HERE IS ALL CHEYENNE'S CODE: 
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 


###### CALCULATE GENOTYPE PROBABILITIES ######

## load est.rf results after running on cluster/server
# 1st scan
#data_sub.ds <- readRDS(file = "./ltreb-qtl-data_sub.ds-est.rf.rds")

## calculate conditional genotype probabilities given multipoint marker data
Costus_prob <- calc.genoprob(Costus)

# 1st scan
write.cross(Costus_prob,file="first_hi.cohort.fam.hotones_covars_AVmap-prob.csv")

saveRDS(Costus_prob,file="first_hi.cohort.fam.hotones_covars_AVmap-prob-prob.rds")

###### SELECT A MODEL ######

# ## determine the appropriate model to use/covariates to include
# days_to_radicle <- as.numeric(pull.pheno(Costus_prob, "days_to_radicle"))
GR_cm.dy <- as.numeric(pull.pheno(Costus_prob, "GR_cm.dy"))
# ldmc <- as.numeric(pull.pheno(Costus_prob, "ldmc"))
# ave_thick <- as.numeric(pull.pheno(Costus_prob, "ave_thick"))
# ave_tough <- as.numeric(pull.pheno(Costus_prob, "ave_tough"))
# N_reabsorp <- as.numeric(pull.pheno(Costus_prob, "N_reabsorp"))
# chloro_reabsorp <- as.numeric(pull.pheno(Costus_prob, "chloro_reabsorp"))
# 
# # potential covariates, per sample; JULIA NEEDS TO calculate and include % heterozygous genotypes
hi <- as.numeric(pull.pheno(Costus_prob, "hybrid_index"))    # C.allenii ancestry proportion
# #het <- as.numeric(pull.pheno(Costus_prob, "heterozygosity")) # % heterozygous genotypes across markers
cohort <- as.factor(pull.pheno(Costus_prob, "cohort"))           # growth cohort (RERUN after converting date to cohort number and including in phenotype data)
fam <- as.factor(pull.pheno(Costus_prob, "family"))           #  F2 parent (RERUN after adding to phenos) 

## select a model by calculating AIC stepwise
#   site.tank is the only recovered covariate
model_all<-lm(GR_cm.dy~hi+cohort+fam)
selectedMod <- step(model_all)

# pull covariates 
#covariates <- pull.pheno(Costus_prob, c("hybrid_index","F21_238","F21_259","F21_241","F21_284","F21_232","coh27","coh24","coh19","coh25","coh26","coh1","coh28","coh31","coh21","coh29","coh30","coh34","coh7","coh8","coh22","coh32","coh33","coh5","coh9","coh6","coh10","coh12","coh15","coh14","coh11","coh13","coh16","coh17","coh20","coh18","coh23")) # ADD fam
#covariates <- pull.pheno(Costus_prob, c("hybrid_index","F21_238","F21_259","F21_241","F21_284","F21_232")) 
# not including all cohorts; for now just including those identified as really problematic for growth rate

## now step through all 31 one-hot encoded cohorts, and families, to remove unneccessary cohorts
coh1 <- as.factor(pull.pheno(Costus_prob, "coh1"))
coh5 <- as.factor(pull.pheno(Costus_prob, "coh5"))
coh6 <- as.factor(pull.pheno(Costus_prob, "coh6"))
coh7 <- as.factor(pull.pheno(Costus_prob, "coh7"))
coh8 <- as.factor(pull.pheno(Costus_prob, "coh8"))
coh9 <- as.factor(pull.pheno(Costus_prob, "coh9"))
coh10 <- as.factor(pull.pheno(Costus_prob, "coh10"))
coh11 <- as.factor(pull.pheno(Costus_prob, "coh11"))
coh12 <- as.factor(pull.pheno(Costus_prob, "coh12"))
coh13 <- as.factor(pull.pheno(Costus_prob, "coh13"))
coh14 <- as.factor(pull.pheno(Costus_prob, "coh14"))
coh15 <- as.factor(pull.pheno(Costus_prob, "coh15"))
coh16 <- as.factor(pull.pheno(Costus_prob, "coh16"))
coh17 <- as.factor(pull.pheno(Costus_prob, "coh17"))
coh18 <- as.factor(pull.pheno(Costus_prob, "coh18"))
coh19 <- as.factor(pull.pheno(Costus_prob, "coh19"))
coh20 <- as.factor(pull.pheno(Costus_prob, "coh20"))
coh21 <- as.factor(pull.pheno(Costus_prob, "coh21"))
coh22 <- as.factor(pull.pheno(Costus_prob, "coh22"))
coh23 <- as.factor(pull.pheno(Costus_prob, "coh23"))
coh24 <- as.factor(pull.pheno(Costus_prob, "coh24"))
coh25 <- as.factor(pull.pheno(Costus_prob, "coh25"))
coh26 <- as.factor(pull.pheno(Costus_prob, "coh26"))
coh27 <- as.factor(pull.pheno(Costus_prob, "coh27"))
coh28 <- as.factor(pull.pheno(Costus_prob, "coh28"))
coh29 <- as.factor(pull.pheno(Costus_prob, "coh29"))
coh30 <- as.factor(pull.pheno(Costus_prob, "coh30"))
coh31 <- as.factor(pull.pheno(Costus_prob, "coh31"))
coh32 <- as.factor(pull.pheno(Costus_prob, "coh32"))
coh33 <- as.factor(pull.pheno(Costus_prob, "coh33"))
coh34 <- as.factor(pull.pheno(Costus_prob, "coh34"))

F21_238 <- as.factor(pull.pheno(Costus_prob, "F21_238"))
F21_259 <- as.factor(pull.pheno(Costus_prob, "F21_259"))
F21_241 <- as.factor(pull.pheno(Costus_prob, "F21_241"))
F21_284 <- as.factor(pull.pheno(Costus_prob, "F21_284"))
F21_232 <- as.factor(pull.pheno(Costus_prob, "F21_232"))

## select a model with cohorts by calculating AIC stepwise
# 1st scan

model_all<-lm(GR_cm.dy~hi+coh5+coh6+coh7+coh8+coh9+coh10+coh11+coh12+coh13+coh14+coh15+coh16+coh17+coh18+coh19+coh20+coh21+coh22+coh23+coh24+coh25+coh26+coh27+coh28+coh29+coh30+coh31+coh32+coh33+coh34+F21_238+F21_259+F21_241+F21_284+F21_232)
selectedMod <- step(model_all)
summary(model_all)

# pull selected covariates (plus hybrid_index)
#  formula = GR_cm.dy ~ coh5 + coh6 + coh7 + coh8 + coh9 + coh10 + 
#coh11 + coh12 + coh13 + coh14 + coh15 + coh16 + coh18 + coh19 + 
#  coh20 + coh21 + coh22 + coh23 + coh26 + coh27 + coh29 + coh30 + 
#  coh31 + F21_259
covariates_GRselect <- pull.pheno(Costus_prob, c("hybrid_index", "coh5","coh6","coh7","coh8","coh9","coh10","coh11","coh12","coh13","coh14","coh15","coh16","coh18","coh19","coh20","coh21","coh22","coh23","coh26","coh27","coh29","coh30","coh31","F21_259"))

###### KK code - 2-dimensional scan #####

## CANT run scantwo locally, could try on humm but for now just do scanone - Error: vector memory exhausted (limit reached?)
#scantwo: Two-dimensional genome scan with a two-QTL model
#scantwo to do permutations for likelihood penalties
scn2_1kperm.hk <- scantwo(Costus_prob, method = "hk", pheno.col = 41:47, addcovar = covariates_GRselect,  n.perm=1000, verbose = TRUE)
plot(scn2_1kperm.hk)

###### RUN SINGLE-QTL SCAN ######
#remove cohort ocovar to look at dormancy bc assoc with hybrid index which is already included. - cohorts really fuck things up...
#PICKUP!!! Eventually - do model selection like Cheyenne for each pheno to see if same cohorts pop up/to include only those that come up across phenos (all from each pheno)
#covariates <- pull.pheno(Costus_prob, c("hybrid_index","F21_238","F21_259","F21_241","F21_284","F21_232")) # ADD fam

## perform a single-qtl scan 
scanone.hk <- scanone(Costus_prob, method="hk", pheno.col = 41:47, addcovar=covariates_GRselect)
plot(scanone.hk,lodcolumn=1, ylab="LOD score - seed dormancy")
plot(scanone.hk,lodcolumn=2, ylab="LOD score - growth rate")
plot(scanone.hk,lodcolumn=3, ylab="LOD score - ldmc")
plot(scanone.hk,lodcolumn=4, ylab="LOD score - thickness")
plot(scanone.hk,lodcolumn=5, ylab="LOD score - toughness")
plot(scanone.hk,lodcolumn=6, ylab="LOD score - %N reabsorbed in drought")
plot(scanone.hk,lodcolumn=7, ylab="LOD score - %chlorophyll reabsorbed in drought")
summary(scanone.hk)
write.table(scanone.hk,file="scanone-hk_hotone_covar_based.on.GR.tsv",sep="\t",quote=FALSE)
#write.table(scanone.hk,file="ltreb-only-CTmax_22-1hot-site-tanks.scanone-hk_second-scan-w-qtl-covariate.tsv",sep="\t",quote=FALSE)
saveRDS(scanone.hk, file="scanone-hk_hotone_covar_based.on.GR.rds")

## run permutations to define the LOD threshold of significance
perm.hk <- scanone(Costus_prob, method="hk", pheno.col = 41:47, n.perm=1000) # ADD BACK  addcovar=covariates,
saveRDS(perm.hk, "scanone-hk_hotone_covar_based.on.GR_1kperm.rds")
#saveRDS(perm.hk, "LTREB-only-ctmax_perm_site-tank.hk.no-xtreme_1hot-site-tank_perm.hk_second-scan-w-qtl-covariate.rds")

## get genome-wide LOD thresholds usi genome-scan-adjusted p-values
# seed dormancy
summary(perm.hk,lodcolumn=1,alpha=c(0.05, 0.1, 0.2))
## 
#lod
# 5%  4.25
# 10% 3.57
# 20% 3.12
cutoff_lod = 3.57 # based on 10%

## look at highest LOD peak and the corresponding p-value on each chromosome
summary(scanone.hk, lodcolumn=1, perms=perm.hk, pvalues=TRUE)

## look at peaks that surpass the 10% LOD threshold
summary(scanone.hk, lodcolumn=1, perms=perm.hk, alpha=0.1, pvalues=TRUE)
summary(scanone.hk[scanone.hk$lod > cutoff_lod,], lodcolumn=1,)

# growth rate
summary(perm.hk,lodcolumn=2,alpha=c(0.05, 0.1, 0.2))
## 
#lod
# 5%  4.25
# 10% 3.57
# 20% 3.12
cutoff_lod = 3.57 # based on 10%

## look at highest LOD peak and the corresponding p-value on each chromosome
summary(scanone.hk, perms=perm.hk, pvalues=TRUE)

## look at peaks that surpass the 10% LOD threshold
summary(scanone.hk, perms=perm.hk, alpha=0.1, pvalues=TRUE)
summary(scanone.hk[scanone.hk$lod > cutoff_lod,])

# ldmc
summary(perm.hk,alpha=c(0.05, 0.1, 0.2))
## 
#lod
# 5%  4.25
# 10% 3.57
# 20% 3.12
cutoff_lod = 3.57 # based on 10%

## look at highest LOD peak and the corresponding p-value on each chromosome
summary(scanone.hk, perms=perm.hk, pvalues=TRUE)

## look at peaks that surpass the 10% LOD threshold
summary(scanone.hk, perms=perm.hk, alpha=0.1, pvalues=TRUE)
summary(scanone.hk[scanone.hk$lod > cutoff_lod,])

# thickness
summary(perm.hk,alpha=c(0.05, 0.1, 0.2))
## 
#lod
# 5%  4.25
# 10% 3.57
# 20% 3.12
cutoff_lod = 3.57 # based on 10%

## look at highest LOD peak and the corresponding p-value on each chromosome
summary(scanone.hk, perms=perm.hk, pvalues=TRUE)

## look at peaks that surpass the 10% LOD threshold
summary(scanone.hk, perms=perm.hk, alpha=0.1, pvalues=TRUE)
summary(scanone.hk[scanone.hk$lod > cutoff_lod,])

# toughness
summary(perm.hk,alpha=c(0.05, 0.1, 0.2))
## 
#lod
# 5%  4.25
# 10% 3.57
# 20% 3.12
cutoff_lod = 3.57 # based on 10%

## look at highest LOD peak and the corresponding p-value on each chromosome
summary(scanone.hk, perms=perm.hk, pvalues=TRUE)

## look at peaks that surpass the 10% LOD threshold
summary(scanone.hk, perms=perm.hk, alpha=0.1, pvalues=TRUE)
summary(scanone.hk[scanone.hk$lod > cutoff_lod,])

# %N reabsorbed in drought
summary(perm.hk,alpha=c(0.05, 0.1, 0.2))
## 
#lod
# 5%  4.25
# 10% 3.57
# 20% 3.12
cutoff_lod = 3.57 # based on 10%

## look at highest LOD peak and the corresponding p-value on each chromosome
summary(scanone.hk, perms=perm.hk, pvalues=TRUE)

## look at peaks that surpass the 10% LOD threshold
summary(scanone.hk, perms=perm.hk, alpha=0.1, pvalues=TRUE)
summary(scanone.hk[scanone.hk$lod > cutoff_lod,])

# %chlorophyll reabsorbed in drought
summary(perm.hk,alpha=c(0.05, 0.1, 0.2))
## 
#lod
# 5%  4.25
# 10% 3.57
# 20% 3.12
cutoff_lod = 3.57 # based on 10%

## look at highest LOD peak and the corresponding p-value on each chromosome
summary(scanone.hk, perms=perm.hk, pvalues=TRUE)

## look at peaks that surpass the 10% LOD threshold
summary(scanone.hk, perms=perm.hk, alpha=0.1, pvalues=TRUE)
summary(scanone.hk[scanone.hk$lod > cutoff_lod,])

# 1st scan: chr22 ScyDAA6-2113-HRSCAF-2539:8811359 ScyDAA6-2113-HRSCAF-2539 8.81 4.37 0.097
# output all markers with LOD score above 10% LOD threshold
write.table(scanone.hk[scanone.hk$lod > cutoff_lod,],file="GR_cm.dy_scanone-hk_hotone_covar-0.1LODcutoff.tsv",sep="\t",quote=FALSE)

## perform a single-qtl scan 
scanone.hk <- scanone(Costus_prob, method="hk", pheno.col = "GR_cm.dy", addcovar=covariates_GRselect)
plot(scanone.hk,col="black")
summary(scanone.hk)
write.table(scanone.hk,file="GR_cm.dy_scanone-hk_hotone_covar.tsv",sep="\t",quote=FALSE)
#write.table(scanone.hk,file="ltreb-only-CTmax_22-1hot-site-tanks.scanone-hk_second-scan-w-qtl-covariate.tsv",sep="\t",quote=FALSE)
saveRDS(scanone.hk, file="GR_cm.dy_scanone-hk_hotone_covar.rds")

## run permutations to define the LOD threshold of significance
perm.hk <- scanone(Costus_prob, method="hk", pheno.col = "days_to_radicle", n.perm=1000) # ADD BACK  addcovar=covariates,
saveRDS(perm.hk, "GR_cm.dy_scanone-hk_hotone_covar.rds")
#saveRDS(perm.hk, "LTREB-only-ctmax_perm_site-tank.hk.no-xtreme_1hot-site-tank_perm.hk_second-scan-w-qtl-covariate.rds")

## get genome-wide LOD thresholds using genome-scan-adjusted p-values
summary(perm.hk,alpha=c(0.05, 0.1, 0.2))
## ctmax ~ hi+17sig.site.tanks <- figure out meaning.. make my equivalent
#lod
# 5%  4.25
# 10% 3.57
# 20% 3.12
cutoff_lod = 3.57 # based on 10%

## look at highest LOD peak and the corresponding p-value on each chromosome
summary(scanone.hk, perms=perm.hk, pvalues=TRUE)

## look at peaks that surpass the 10% LOD threshold
summary(scanone.hk, perms=perm.hk, alpha=0.1, pvalues=TRUE)
summary(scanone.hk[scanone.hk$lod > cutoff_lod,])

# 1st scan: chr22 ScyDAA6-2113-HRSCAF-2539:8811359 ScyDAA6-2113-HRSCAF-2539 8.81 4.37 0.097
# output all markers with LOD score above 10% LOD threshold
write.table(scanone.hk[scanone.hk$lod > cutoff_lod,],file="GR_cm.dy_scanone-hk_hotone_covar-0.1LODcutoff.tsv",sep="\t",quote=FALSE)

## plot chr22 peak
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl.pdf",8.4,8.9)
plot(scanone.hk,chr="ScyDAA6-2113-HRSCAF-2539") # chr22
abline(h=cutoff_lod, col="red", lwd=3)
title("22@8.81")
dev.off()


###### LOD peak support ######

## get the 1.5-LOD support interval for significant peaks
lint = cutoff_lod - 1.5   # 2.83
lodint(scanone.hk,"ScyDAA6-2113-HRSCAF-2539",drop=1.5,expandtomarkers = TRUE)
#ScyDAA6-2113-HRSCAF-2539:7006520  ScyDAA6-2113-HRSCAF-2539  7.006520 2.713022
#ScyDAA6-2113-HRSCAF-2539:8811359  ScyDAA6-2113-HRSCAF-2539  8.811359 4.366771
#ScyDAA6-2113-HRSCAF-2539:10301012 ScyDAA6-2113-HRSCAF-2539 10.301012 2.856064

## get the 2.0-LOD support interval for significant peaks
lint2 = cutoff_lod - 2   # 2.33
lodint(scanone.hk,"ScyDAA6-2113-HRSCAF-2539",drop=2,expandtomarkers = TRUE)
#ScyDAA6-2113-HRSCAF-2539:6963205  ScyDAA6-2113-HRSCAF-2539  6.963205 2.257862
#ScyDAA6-2113-HRSCAF-2539:8811359  ScyDAA6-2113-HRSCAF-2539  8.811359 4.366771
#ScyDAA6-2113-HRSCAF-2539:12213851 ScyDAA6-2113-HRSCAF-2539 12.213851 2.163837

## get 95% Bayes credible interval
bayesint(scanone.hk,"ScyDAA6-2113-HRSCAF-2539",0.95,expandtomarkers = TRUE)
#ScyDAA6-2113-HRSCAF-2539:6984383  ScyDAA6-2113-HRSCAF-2539  6.984383 2.412447
#ScyDAA6-2113-HRSCAF-2539:8811359  ScyDAA6-2113-HRSCAF-2539  8.811359 4.366771
#ScyDAA6-2113-HRSCAF-2539:12192947 ScyDAA6-2113-HRSCAF-2539 12.192947 2.479510

#out.boot <- scanoneboot(data_prob,chr="ScyDAA6-2113-HRSCAF-2539",n.boot=1000,prob=0.95)
#summary(out.boot)
#plot(out.boot)


###### MAKE EFFECT PLOTS ######

### Payne et al 2021 effect plots

##chr22 qtl effect plot
setwd("~/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/LTREB_qtl/CTmax-QTL_first-scan_data")
data_prob.df <- read.csv("ltreb-ctmax-qtl-filtered_22-1hot-site-tank_data-prob.csv")
qtl_genos <- data_prob.df$ScyDAA6.2113.HRSCAF.2539.8811359
qtl_genos[qtl_genos == "-"] <- NA
ctmax <- data_prob.df$ctmax
qtl22_peak_data <- data.frame(ctmax=ctmax,qtl22_geno=qtl_genos)
qtl22_peak_data <- na.omit(qtl22_peak_data)

BB<-qtl22_peak_data[qtl22_peak_data$qtl22_geno=="BB",]$ctmax
MB<-qtl22_peak_data[qtl22_peak_data$qtl22_geno=="AB",]$ctmax
MM<-qtl22_peak_data[qtl22_peak_data$qtl22_geno=="AA",]$ctmax

malcol_22=rgb(0/255,0/255,175/255)
hetcol_22=rgb(100/255,0/255,175/255)
bircol_22=rgb(150/255,0/255,0/255)

## error bar function (plot 1 SD)
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

se <- function(x) 2*(sd(x)/sqrt(length(x)))

pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl-effect-plot-w-points_2se.pdf",5.5,7)

plot(1:3,c(mean(MM),mean(MB),mean(BB)),col=c(malcol_22,hetcol_22,bircol_22),ylim=c(32,38),xlab="Genotype",ylab="Normalized counts",xaxt="n",pch="-",cex=3,xlim=c(0.5,3.5))
error.bar(1:3,c(mean(MM),mean(MB),mean(BB)),c(se(MM),se(MB),se(BB)),col=c(malcol_22,hetcol_22,bircol_22),lwd=2)

malcol=rgb(0/255,0/255,175/255,alpha=0.3)
hetcol=rgb(100/255,0/255,175/255,alpha=0.6)
bircol=rgb(150/255,0/255,0/255,alpha=0.6)

noise<-runif(length(BB),0.2,0.35)
points(rep(3,length(BB))+noise,BB,pch=20,cex=1.8,col=bircol)

noise<-runif(length(MB),0.2,0.35)
points(rep(2,length(MB))+noise,MB,pch=20,cex=1.8,col=hetcol)

noise<-runif(length(MM),0.2,0.35)
points(rep(1,length(MM))+noise,MM,pch=20,cex=1.8,col=malcol)

mtext(c("MM","MB","BB"),at=1:3,side=1)

dev.off()


### other effect plots
## Load data_prob object
data_prob <- readRDS("ltreb-ctmax-qtl-filtered_22-1hot-site-tank_data-prob.rds")

## impute missing genotypes
data_sim <- sim.geno(data_prob,step=1,n.draws=16)

## plot estimated effect of QTL: phenotype vs marker genotype
malcol=rgb(0/255,0/255,175/255)
hetcol=rgb(100/255,0/255,175/255)
bircol=rgb(150/255,0/255,0/255)
effect_colors = c(malcol,hetcol,bircol)

## plot estimated effect of QTL: phenotype vs marker genotype
# 1st scan
mar_chr22 <- find.marker(data_prob, chr='ScyDAA6-2113-HRSCAF-2539', pos=8.811359)
# simple effect plot
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl-effect-plot.pdf",8.4,8.9)
effectplot(data_sim,mname1=mar_chr22,main='22@8.81')
dev.off()
# effect plots with individuals as points
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl-effect-plot-w-points.pdf",5.5,7)
plotPXG(data_prob,marker=mar_chr22,main='22@8.81',infer=FALSE, ylab = expression('CT'['max']*' ('*~degree*C*')'), col = effect_colors)
dev.off()

## get QTL effect on CTmax per genotype at peak
eff <- effectplot(data_sim,mname1=mar_chr22)
eff
# ScyDAA6-2113-HRSCAF-2539:8811359
#      mean      SE
#  AA: 35.98525, 0.1777719
#  AB: 35.65896, 0.1180658
#  BB: 36.00340, 0.1713042

## plot estimated QTL effects along chr22 (shows additive and dominance effects)
# effect summary for 22@8.81:
# a           d             se.a        se.d
# 0.008494465	-0.337265149	0.124726314	0.172637572
effect_table<-effectscan(data_sim, pheno.col="ctmax", c("ScyDAA6-2113-HRSCAF-2539"), draw=F, get.se=T)
write.table(effect_table,file="ltreb-only-CTmax_17-1hot-site-tanks.chr22-qtl_effect-table.tsv",sep="\t",quote=FALSE)
# plot effect across chr22
pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-additive-dominance-effects-plot.pdf",9.2,6.8)
effectscan(data_sim, pheno.col="ctmax", c("ScyDAA6-2113-HRSCAF-2539"), draw=T, get.se=T)
dev.off()


###### MULTI-QTL ANALYSIS ######

# create QTL object with both loci
# what="prob" will out the QTL genotype probabilities for use in HK regression
qtl <- makeqtl(data_prob, chr=c("ScyDAA6-2113-HRSCAF-2539","ScyDAA6-5984-HRSCAF-6694"), pos=c(8.81, 4.41), what="prob")

# set covariates
covariates <- pull.pheno(data_prob, c("hybrid_index","STH.1","STH.2","STH.3","STH.4","STH.7","STH.8","STL.1","STL.2","STL.4","STL.5","STL.6","STL.7","STM.2",
                                      "STM.4","STM.5","STM.6","STM.7"))

# fit the two locus additive model
# “drop one term at a time” table compares the fit of the two-QTL model to the reduced model where one QTL is omitted.
out.fq <- fitqtl(data_prob, pheno.col=1, qtl=qtl, method="hk", get.ests=TRUE, covar=covariates, formula=y~Q1+Q2 + hybrid_index + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + STM.6 + STM.7)
summary(out.fq)

# determine whether there is an interaction between the two QTL by fitting the model with the interaction
out.fqi <- fitqtl(data_prob, pheno.col=1, qtl=qtl, method="hk", get.ests=TRUE, covar=covariates, formula=y~Q1+Q2+Q1:Q2 + hybrid_index + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + STM.6 + STM.7)
summary(out.fqi)

# also assess interaction with addint, which adds one interaction at a time
addint(data_prob, qtl=qtl, method="hk")

# refine the qtl positions
# each QTL is moved to the position giving the highest likelihood,
# and the entire process is repeated until no further improvement in likelihood can be obtained
rqtl <- refineqtl(data_prob, qtl=qtl, method="hk")
rqtl
#                             name                      chr     pos n.gen
# Q1 ScyDAA6-2113-HRSCAF-2539@10.0 ScyDAA6-2113-HRSCAF-2539 10.0490     3
# Q2  ScyDAA6-5984-HRSCAF-6694@5.7 ScyDAA6-5984-HRSCAF-6694  5.7399     3

# look for additional qtl, using the refined qtl positions
out.aq_1qtl <- addqtl(data_prob, qtl=rqtl, method="hk", covar=covariates, formula=y~Q1 + hybrid_index + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + STM.6 + STM.7)
out.aq_2qtl <- addqtl(data_prob, qtl=rqtl, method="hk", covar=covariates, formula=y~Q1+Q2+Q1:Q2 + hybrid_index + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + STM.6 + STM.7)

plot(out.aq_1qtl)
summary(out.aq_1qtl)
out.aq_1qtl[out.aq_1qtl$lod > 3,]
# as expected, recovers chr15 qtl
# ScyDAA6-5984-HRSCAF-6694:5847498  ScyDAA6-5984-HRSCAF-6694 5.85e+00 3.503
plot(out.aq_2qtl)
summary(out.aq_2qtl)
out.aq_2qtl[out.aq_2qtl$lod > 3,]
#                                                        chr      pos   lod
# ScyDAA6-1854-HRSCAF-2213:18020066 ScyDAA6-1854-HRSCAF-2213 1.80e+01 3.512
# ScyDAA6-2469-HRSCAF-2980:1032301  ScyDAA6-2469-HRSCAF-2980 1.03e+00 3.360

# grab the 1.5-LOD for both of these
lod_cutoff = 3.512-1.5
write.csv(out.aq_2qtl[out.aq_2qtl$chr == "ScyDAA6-1854-HRSCAF-2213" & out.aq_2qtl$lod > lod_cutoff,], "CTmax-chr22-chr15-add-qtl_ScyDAA6-1854-HRSCAF-2213_1.5LOD.csv")
lod_cutoff = 3.36-1.5
write.csv(out.aq_2qtl[out.aq_2qtl$chr == "ScyDAA6-2469-HRSCAF-2980" & out.aq_2qtl$lod > lod_cutoff,], "CTmax-chr22-chr15-add-qtl_ScyDAA6-1854-HRSCAF-2213_1.5LOD.csv")


###### ROUGHLY ESTIMATE EFFECT OF QTL WITH AN LM ######

### get a very rough estimate of the QTL effect size with R^2
## load infile as data.frame and do correlations/lm with that
data_prob.df <- read.csv("ltreb-ctmax-qtl-filtered_22-1hot-site-tank_data-prob.csv")

## grab all relevant variables
ctmax <- data_prob.df$ctmax
hi <- data_prob.df$hybrid_index
STH.1<-data_prob.df$STH.1
STH.2<-data_prob.df$STH.2
STH.3<-data_prob.df$STH.3
STH.4<-data_prob.df$STH.4
STH.7<-data_prob.df$STH.7
STH.8<-data_prob.df$STH.8
STL.1<-data_prob.df$STL.1
STL.2<-data_prob.df$STL.2
STL.4<-data_prob.df$STL.4
STL.5<-data_prob.df$STL.5
STL.6<-data_prob.df$STL.6
STL.7<-data_prob.df$STL.7
STM.2<-data_prob.df$STM.2
STM.4<-data_prob.df$STM.4
STM.5<-data_prob.df$STM.5
STM.6<-data_prob.df$STM.6
STM.7<-data_prob.df$STM.7

# grab the genotypes of the locus at the QTL peak
chr22_geno <- data_prob.df$ScyDAA6.2113.HRSCAF.2539.8811359
chr15_geno <- data_prob.df$ScyDAA6.5984.HRSCAF.6694.4408039

## optional: output genotypes for markers that fall under the QTL / within the 2-LOD interval: 6.963205 - 12.213851
firstcol = which(colnames(data_prob.df)=="ScyDAA6.2113.HRSCAF.2539.6963205")
lastcol = which(colnames(data_prob.df)=="ScyDAA6.2113.HRSCAF.2539.12213851")
qtl_genos <- data_prob.df[c(firstcol:lastcol)]
qtl_geno_outdf <- data.frame(ctmax,hi,STH.1,STH.2,STH.3,STH.4,STH.7,STH.8,STL.1,STL.2,STL.4,STL.5,STL.6,STL.7,STM.2,STM.4,STM.5,STM.6,STM.7,qtl_genos)
write.csv(qtl_geno_outdf,'./ltreb-ctmax-qtl_genos-all-markers-under-2LOD-chr22-qtl.csv',row.names=F,quote=F)

### chr22 QTL LM effect size estimate and p-values
## organize variables of interest into new data frame, remove missing data
dp_noNA<-data.frame(ctmax,hi,chr22_geno,STH.1,STH.2,STH.3,STH.4,STH.7,STH.8,STL.1,STL.2,STL.4,STL.5,STL.6,STL.7,STM.2,STM.4,STM.5,STM.6,STM.7)
dp_noNA[dp_noNA == '-']<-NA
dp_noNA<-na.omit(dp_noNA)
dim(dp_noNA)

## perform a linear regression with the full model and get the R^2 value
# full model adjusted R^2: 0.6731 (pval < 2e-16)
lm.res<-lm(ctmax ~ chr22_geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 +
             STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 +
             STM.6 + STM.7,data=dp_noNA)
summary(lm.res)

## get the adjusted R^2 from the null model (i.e. without genotype as a coefficient)
# null adj R^2: 0.6289 (pval < 2e-16)
null.lm.res<-lm(ctmax ~ hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 +
                  STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 +
                  STM.6 + STM.7,data=dp_noNA)
summary(null.lm.res)

## subtract the adjusted R^2 of the null model from the full model to get effect estimate
effect = 0.6731 - 0.6289 # 0.0442 (4.42%)


## get log likelihood difference between focal and null models
summary(lm.res)
model2<-glm(as.numeric(ctmax)~as.numeric(hi)+
              as.factor(STH.1) + as.factor(STH.2) + as.factor(STH.3) +
              as.factor(STH.4) + as.factor(STH.7) + as.factor(STH.8) +
              as.factor(STL.1) + as.factor(STL.2) + as.factor(STL.4) +
              as.factor(STL.5) + as.factor(STL.6) + as.factor(STL.7) +
              as.factor(STM.2) + as.factor(STM.4) + as.factor(STM.5) +
              as.factor(STM.6) + as.factor(STM.7),data=dp_noNA,family="gaussian")
null<-logLik(model2)[1]

model1<-glm(as.numeric(ctmax)~as.numeric(hi)+
              as.factor(STH.1) + as.factor(STH.2) + as.factor(STH.3) +
              as.factor(STH.4) + as.factor(STH.7) + as.factor(STH.8) +
              as.factor(STL.1) + as.factor(STL.2) + as.factor(STL.4) +
              as.factor(STL.5) + as.factor(STL.6) + as.factor(STL.7) +
              as.factor(STM.2) + as.factor(STM.4) + as.factor(STM.5) +
              as.factor(STM.6) + as.factor(STM.7) + as.factor(geno),data=dp_noNA,family="gaussian")
focal<-logLik(model1)[1]

# calculate log likelihood difference
like_diff<-focal-null
like_diff # 10.03532


### chr15 QTL LM effect size estimate
## organize variables of interest into new data frame, remove missing data
dp_noNA<-data.frame(ctmax,hi,chr15_geno,STH.1,STH.2,STH.3,STH.4,STH.7,STH.8,STL.1,STL.2,STL.4,STL.5,STL.6,STL.7,STM.2,STM.4,STM.5,STM.6,STM.7)
dp_noNA[dp_noNA == '-']<-NA
dp_noNA<-na.omit(dp_noNA)
dim(dp_noNA)
# full model adjusted R^2: 0.6417 (pval < 2e-16)
lm.res<-lm(ctmax ~ chr15_geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 +
             STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 +
             STM.6 + STM.7,data=dp_noNA)
summary(lm.res)

## get the adjusted R^2 from the null model (i.e. without genotype as a coefficient)
# null adj R^2: 0.6203 (pval < 2e-16)
null.lm.res<-lm(ctmax ~ hi +STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 +
                  STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 +
                  STM.6 + STM.7,data=dp_noNA)
summary(null.lm.res)

## subtract the adjusted R^2 of the null model from the full model to get effect estimate
effect = 0.6417 - 0.6203 # 0.214 (2.14%)

### both QTL LM effect size estimate
## organize variables of interest into new data frame, remove missing data
dp_noNA<-data.frame(ctmax,hi,chr22_geno,chr15_geno,STH.1,STH.2,STH.3,STH.4,STH.7,STH.8,STL.1,STL.2,STL.4,STL.5,STL.6,STL.7,STM.2,STM.4,STM.5,STM.6,STM.7)
dp_noNA[dp_noNA == '-']<-NA
dp_noNA<-na.omit(dp_noNA)
dim(dp_noNA)
# full model adjusted R^2: 0.7126 (pval < 2e-16)
lm.res<-lm(ctmax ~ chr22_geno + chr15_geno + chr22_geno:chr15_geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 +
             STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 +
             STM.6 + STM.7,data=dp_noNA)
summary(lm.res)

## get the adjusted R^2 from the null model (i.e. without genotype as a coefficient)
# null adj R^2: 0.6306 (pval < 2e-16)
null.lm.res<-lm(ctmax ~ hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 +
                  STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 +
                  STM.6 + STM.7,data=dp_noNA)
summary(null.lm.res)

## subtract the adjusted R^2 of the null model from the full model to get effect estimate
effect = 0.7126 - 0.6306 # 0.082 (8.2%)


### chr22:15 QTL interaction LM effect size estimate
## organize variables of interest into new data frame, remove missing data
library(mltools)
library(data.table)
library(rstatix)
dp_noNA<-data.frame(ctmax,hi,chr22_geno,chr15_geno,STH.1,STH.2,STH.3,STH.4,STH.7,STH.8,STL.1,STL.2,STL.4,STL.5,STL.6,STL.7,STM.2,STM.4,STM.5,STM.6,STM.7)
dp_noNA[dp_noNA == '-']<-NA
dp_noNA<-na.omit(dp_noNA)
dim(dp_noNA)

# full model adjusted R^2: 0.7126 (pval < 2e-16)
lm.res<-lm(ctmax ~ chr22_geno + chr15_geno + chr22_geno:chr15_geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 +
             STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 +
             STM.6 + STM.7,data=dp_noNA)
summary(lm.res)

## get the adjusted R^2 from the null model (i.e. without genotype as a coefficient)
# null adj R^2: 0.687 (pval < 2e-16)
null.lm.res<-lm(ctmax ~ chr22_geno + chr15_geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 +
                  STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 +
                  STM.6 + STM.7,data=dp_noNA)
summary(null.lm.res)

## subtract the adjusted R^2 of the null model from the full model to get effect estimate
effect = 0.7126 - 0.687 # 0.0256 (2.6%)

## eta squared to get effect size per variable
eta_squared(aov(ctmax ~ chr22_geno + chr15_geno + chr22_geno:chr15_geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 + STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 + STM.6 + STM.7, dp_noNA))


### chr22:chr15 QTL comparison of means w anova/tukey
lm.res<-lm(ctmax ~ chr22_geno + chr15_geno + chr22_geno:chr15_geno + hi + STH.1 + STH.2 + STH.3 + STH.4 + STH.7 + STH.8 + STL.1 +
             STL.2 + STL.4 + STL.5 + STL.6 + STL.7 + STM.2 + STM.4 + STM.5 +
             STM.6 + STM.7,data=dp_noNA)
anov <- aov(lm.res)
TukeyHSD(x=anov, 'chr22_geno', conf.level=0.95)
TukeyHSD(x=anov, 'chr15_geno', conf.level=0.95)
TukeyHSD(x=anov, 'chr22_geno:chr15_geno', conf.level=0.95)


###### CREATE A NICE MANHATTAN PLOT OF GENOME-WIDE LODs ######

### To generate an organized genome-wide manhattan plot, convert chromosome names to numbers using
#     ./replace-chrom-name.py ./match_xmac_chroms_xbir-10x_chroms.txt <infile>
#   in ./Scripts/
scanone.hk.df_rename <- read.csv("ltreb-only-CTmax_17-1hot-site-tanks.scanone-hk_clean.tsv_chr-renamed.csv")

## load the qqman manhattan plot function
source("Scripts/input_files/adapted_qqman.R")

## Function: prepares data frame with chr,pos,lod columns to create manhattan plot
#  make sure 'cutoff_lod' is set to the LOD threshold
mh_plot<-function(chr_pos_lod){
  cpl.df <- as.data.frame(chr_pos_lod)
  data_trim<-{}
  # if using non-integer chrom names (i.e. not group numbers),
  # use this instead of following line
  # data_trim$CHR<-chr_to_group(as.character.factor(cpl.df$chr))
  data_trim$CHR<-cpl.df$chr
  data_trim$BP<-(cpl.df$pos)*1e6
  data_trim$P<-cpl.df$lod
  data_trim<-as.data.frame(data_trim)
  data_trim<-na.omit(data_trim)
  # data_trim$CHR <- as.character.factor(data_trim$CHR)
  data_trim$CHR<-data_trim$CHR
  manhattan(data_trim,genomewideline=cutoff_lod,suggestiveline=-1, cex.axis=1.5, cex.lab=4)
  
  return(data_trim)
}

## look at all of the peaks that surpass the genome-wide LOD threshold
scanone.hk.df_rename[scanone.hk.df_rename$lod>cutoff_lod,]
# ScyDAA6-2113-HRSCAF-2539:8811359  22 8.811359 4.366771
# ScyDAA6-2113-HRSCAF-2539:8831590  22 8.831590 4.361586
# ScyDAA6-2113-HRSCAF-2539:8851788  22 8.851788 4.356347

#  and all those under the 2-LOD interval
scanone.hk.df_rename[scanone.hk.df_rename$lod>lint2 & scanone.hk.df_rename$chr == "22" ,]
write.table(scanone.hk.df_rename[scanone.hk.df_rename$lod>lint2 & scanone.hk.df_rename$chr == "22" ,],'ltreb-ctmax_grp22-all-markers-under-2LOD_chr-pos-lod.tsv',sep="\t",row.names=F,quote=F)

## plot
pdf("chr22-ctmax-qtl_first-scan_genome-wide-manhattan.pdf",10.5,5.5)
mh_plot(scanone.hk.df_rename)                                         # genome-wide
title(main="LTREB CTmax QTL: QTL map ctmax~geno+hi+17sig.site.tanks")
dev.off()
pdf("chr22-ctmax-qtl_first-scan_chr22-manhattan.pdf",5.3,5.3)
mh_plot(scanone.hk.df_rename[scanone.hk.df_rename$chr == "22" ,])     # particular chromosome
dev.off()


###### CALCULATE CORRELATION BETWEEN CTMAX AND GENOME-WIDE ANCESTRY ######

## load phenotype data (individuals with excess hybrid_index and heterozygosity have already been removed)
pheno_data <- read.table('LTREB-only-CTmax-phenos-w-hi-het_1hot-site-tank_clean.tsv',header=T)
pheno_data <- pheno_data[pheno_data$hybrid_index>0.15 & pheno_data$hybrid_index<0.85,] # only keep individuals with hybrid index between 0.15 and 0.85
hi<-pheno_data$hybrid_index     # genome-wide X.malinche ancestry proportion
het<-pheno_data$heterozygosity  # genome-wide heterozygosity proportion
ctmax<-pheno_data$ctmax         # CTmax

## look at correlation between CTmax and genome-wide X. malinche ancestry (hybrid index)
model <- lm(ctmax~hi)
summary(model)
cor.test(ctmax,hi,method="pearson") # ctmax,hi: r=+0.1185137, pval=0.1472
# plot correlation
pdf("CTmax-vs-genome-wide-ancestry_corr-plot.pdf",6,6)
plot(hi,ctmax,pch = 16, cex = 1.3, col = "light grey", main = "CTmax and genome-wide X. malinche ancestry",xlab="hybrid index",ylab="CTmax (Celsius)")
abline(34.7295,2.1980, lwd = 3, col="dark blue")
dev.off()

## look at correlation between CTmax and genome-wide heterozygosity)
model <- lm(ctmax~het)
summary(model)
cor.test(ctmax,het,method="pearson") # ctmax,het: r=-0.02285036, pval=0.7806
# plot correlation
pdf("CTmax-vs-genome-wide-heterozygosity_corr-plot.pdf",6,6)
plot(het,ctmax,pch = 16, cex = 1.3, col = "light salmon", main = "CTmax and genome-wide heterozygosity",xlab="heterozygosity",ylab="CTmax (Celsius)")
abline(35.9760,-0.2317, lwd = 3, col="dark blue")
dev.off()









































pdf("~/Documents/R:QTL/plotRF_preformLinkageGroups.pdf", height = 70)
par(mfrow=c(9,1))
plotRF(Costus, chr = 1, alternate.chrid=TRUE)
plotRF(Costus, chr = 2, alternate.chrid=TRUE)
plotRF(Costus, chr = 3, alternate.chrid=TRUE)
plotRF(Costus, chr = 4, alternate.chrid=TRUE)
plotRF(Costus, chr = 5, alternate.chrid=TRUE)
plotRF(Costus, chr = 6, alternate.chrid=TRUE)
plotRF(Costus, chr = 7, alternate.chrid=TRUE)
plotRF(Costus, chr = 8, alternate.chrid=TRUE)
plotRF(Costus, chr = 9, alternate.chrid=TRUE)
dev.off()

## infer linkage groups
## Two markers will be placed in the same linkage groups if they have estimated 
## recombination fraction ≤ max.rf and LOD score ≥ min.lod

#lg <- formLinkageGroups(Costus, max.rf=0.35, min.lod=75) 
lg <- formLinkageGroups(Costus, max.rf=0.35, min.lod=37) 
table(lg[,2])
chromosomes <- substr(rownames(lg), 6, 6)

lg$LG[chromosomes == "1"] ## chromosome 1 is LG 8 and 10
lg$LG[chromosomes == "2"] ## chromosome 2 is LG 1 and 4
lg$LG[chromosomes == "3"] ## chromosome 3 is LG 3
lg$LG[chromosomes == "4"] ## chromosome 4 is LG 4
lg$LG[chromosomes == "5"] ## chromosome 5 is LG 2 and 7 
lg$LG[chromosomes == "6"] ## chromosome 6 is LG 5
lg$LG[chromosomes == "7"] ## chromosome 7 is LG 2
lg$LG[chromosomes == "8"] ## chromosome 8 is LG 6
lg$LG[chromosomes == "9"] ## chromosome 9 is LG 9

## reorganize the markers into these inferred linkage groups
Costus <- formLinkageGroups(Costus, max.rf=0.35, min.lod=40, reorgMarkers=TRUE)

## A plot of the pairwise recombination fractions and LOD scores may indicate 
## how well this worked
pdf("~/Documents/R:QTL/plotRF_preorderMarkers.pdf", height = 80)
par(mfrow=c(9,1))
plotRF(Costus, chr = 1, alternate.chrid=TRUE)
plotRF(Costus, chr = 2, alternate.chrid=TRUE)
plotRF(Costus, chr = 3, alternate.chrid=TRUE)
plotRF(Costus, chr = 4, alternate.chrid=TRUE)
plotRF(Costus, chr = 5, alternate.chrid=TRUE)
plotRF(Costus, chr = 6, alternate.chrid=TRUE)
plotRF(Costus, chr = 7, alternate.chrid=TRUE)
plotRF(Costus, chr = 8, alternate.chrid=TRUE)
plotRF(Costus, chr = 9, alternate.chrid=TRUE)
dev.off()

## plot of LOD scores versus recombination fractions for all pairs (can't do this if using markerlrt)
rf <- pull.rf(Costus)
lod <- pull.rf(Costus, what="lod")
sample <- sample(1:length(rf), 1000)
plot(as.numeric(rf)[sample], as.numeric(lod)[sample], xlab="Recombination fraction", ylab="LOD score")


## Order markers on chromosome 1
Costus <- orderMarkers(Costus, chr=1, window = 2) 
## Order markers on chromosome 2
#Costus <- orderMarkers(Costus, chr=2, window = 2) 
## Order markers on chromosome 8
Costus <- orderMarkers(Costus, chr=8, window = 2) 



pdf("~/Documents/R:QTL/plotRF_postorderMarkers.pdf", height = 80)
par(mfrow=c(9,1))
plotRF(Costus, chr = 1, alternate.chrid=TRUE)
plotRF(Costus, chr = 2, alternate.chrid=TRUE)
plotRF(Costus, chr = 3, alternate.chrid=TRUE)
plotRF(Costus, chr = 4, alternate.chrid=TRUE)
plotRF(Costus, chr = 5, alternate.chrid=TRUE)
plotRF(Costus, chr = 6, alternate.chrid=TRUE)
plotRF(Costus, chr = 7, alternate.chrid=TRUE)
plotRF(Costus, chr = 8, alternate.chrid=TRUE)
plotRF(Costus, chr = 9, alternate.chrid=TRUE)
dev.off()

pdf("~/Documents/R:QTL/plotRF_postorderMarkers_all.pdf")
par(mfrow=c(1,1))
plotRF(Costus, alternate.chrid=TRUE)
dev.off()

summaryMap(Costus)

pdf("~/Documents/R:QTL/plotMap.pdf", width = 12)
plotMap(Costus, show.marker.names=FALSE)
dev.off()

newmap <- est.map(Costus)
logliks <- sapply(newmap, attr, "loglik")
plotMap(Costus, newmap)

##### KATHLEEN'S CODE BELOW
#Calculate conditional genotype probabilities
Costus <- calc.genoprob(Costus, step=2)

#scantwo: Two-dimensional genome scan with a two-QTL model
#scantwo to do permutations for likelihood penalties
operm2.hk <- scantwo(Costus, method = "hk", pheno.col = 2:24, n.cluster= 3, n.perm=1000, verbose = TRUE)
summary(operm2.hk)
summary(out2.hk, perms=operm2.hk, pvalues=TRUE,
        alphas=c(0.05, 0.05, 0.05, 0.05, 0.05), pheno.col = 2)


operm2.hk <- scantwo(Costus, chr=1, method = "hk", pheno.col = 2, n.cluster= 3, n.perm=1000, verbose = TRUE)
summary(operm2.hk)

#calculate penalties to be used in stepwiseqtl based on scantwo permutations
Costus_penalties <- calc.penalties(operm2.hk, alpha = 0.05)

#this step took forever, so writing the results to csv files
write.csv(Costus_penalties, file = "Costus_penalties.csv")

Costus_penalties <- read.csv(file = "Costus_penalties.csv", header = TRUE, row.names = 1)

#stepwiseqtl: Stepwise selection for multiple QTL
#INFA
INFAqtl <- stepwiseqtl(Costus, pheno.col = "INFA", method = "hk", model = "normal", 
                      penalties = as.vector(Costus_penalties[1:3]), max.qtl = 8, keeplodprofile = TRUE,
                      verbose = TRUE, refine.locations = TRUE)
summary(INFAqtl)
plot(INFAqtl)
#fitqtl: Fit a multiple QTL model
INFA.fit <- fitqtl(Costus, pheno.col = "INFA", method = "hk", model = "normal", get.ests = TRUE,
                  qtl = INFAqtl, formula = attr(INFAqtl,"formula"), dropone = TRUE)
summary(INFA.fit)

#Calculate an approximate Bayes credible interval for location
bayesint(INFAqtl, qtl.index = 1)
bayesint(INFAqtl, qtl.index = 2)
bayesint(INFAqtl, qtl.index = 3)




pdf(file="INFA_LOD.pdf", width = 7.5, height = 2.3, pointsize = 8) 
# Plot 1-dimensional LOD profiles for a multiple QTL model
plotLodProfile(INFAqtl, showallchr = TRUE, qtl.labels = FALSE, lwd=1, xlab = "", ylab = "", 
               main ="INFA", frame.plot=FALSE)
dev.off()


