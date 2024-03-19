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
# No. phenotypes:     47 
# Percent phenotyped: 97.9 
# 
# No. chromosomes:    9 
# Autosomes:      1 2 3 4 5 6 7 8 9 
# 
# Total markers:      14509 
# No. markers:        1974 2124 1564 1393 1879 1492 1267 1412 1404 
# Percent genotyped:  98.4 
# Genotypes (%):      AA:25.8  AB:27.0  BB:47.2  not BB:0.0  not AA:0.0 

## look at the pattern of missing data (can take a few minutes to load)
#plotMissing(Costus)

# Plot of number of genotyped markers for each individual (left panel) and
# number of genotyped individuals for each marker (right panel).
par(mfrow=c(1,2), las=1)
plot(ntyped(Costus), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(Costus, "mar"), ylab="No. typed individuals", main="No. genotypes by marker")

## Omit individuals with less than 80% of markers (14509*.80 = 11607)
Costus <- subset(Costus, ind=(ntyped(Costus)>11607))
# none removed

## identify names of markers to drop, based on genotypes-per-marker plot
# drop all markers with less than ~80% genotyped individuals (391*0.8 = 313)
nt.bymar <- ntyped(Costus, "mar")
todrop <- names(nt.bymar[nt.bymar < 313]) 
## drop those markers
Costus <- drop.markers(Costus, todrop)

# summary(Costus) # brings total markers 14509 --> 14416 , % genotyped from 98.4 --> 98.6%

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
days_to_radicle <- as.numeric(pull.pheno(Costus_prob, "days_to_radicle"))
GR_cm.dy <- as.numeric(pull.pheno(Costus_prob, "GR_cm.dy"))
# ldmc <- as.numeric(pull.pheno(Costus_prob, "ldmc"))
ave_thick <- as.numeric(pull.pheno(Costus_prob, "ave_thick"))
ave_tough <- as.numeric(pull.pheno(Costus_prob, "ave_tough"))
N_reabsorp <- as.numeric(pull.pheno(Costus_prob, "N_reabsorp"))
chloro_reabsorp <- as.numeric(pull.pheno(Costus_prob, "chloro_reabsorp"))
# 
# # potential covariates, per sample; JULIA NEEDS TO calculate and include % heterozygous genotypes
hi <- as.numeric(pull.pheno(Costus_prob, "hybrid_index"))    # C.allenii ancestry proportion
# #het <- as.numeric(pull.pheno(Costus_prob, "heterozygosity")) # % heterozygous genotypes across markers
cohort <- as.factor(pull.pheno(Costus_prob, "cohort"))           # growth cohort (RERUN after converting date to cohort number and including in phenotype data)
fam <- as.factor(pull.pheno(Costus_prob, "family"))           #  F2 parent (RERUN after adding to phenos) 

## select a model by calculating AIC stepwise
model_all<-lm(days_to_radicle~hi+cohort+fam)
selectedMod <- step(model_all)
summary(selectedMod)
# keep all (didn't include hi but I will anyway)

model_all<-lm(GR_cm.dy~hi+cohort+fam)
selectedMod <- step(model_all)
summary(selectedMod)
# keep all (didn't include hi but I will anyway)

model_all<-lm(ave_thick~hi+cohort+fam)
selectedMod <- step(model_all)
summary(selectedMod)
# keep all (didn't include hi but I will anyway)

model_all<-lm(ave_tough~hi+cohort+fam)
selectedMod <- step(model_all)
summary(selectedMod)
# keep all (didn't include hi but I will anyway)

model_all<-lm(N_reabsorp~hi+cohort+fam)
selectedMod <- step(model_all)
summary(selectedMod)
# keep all (didn't include hi but I will anyway)

model_all<-lm(chloro_reabsorp~hi+cohort+fam)
selectedMod <- step(model_all)
summary(selectedMod)
# keep all (didn't include hi but I will anyway)


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

## select models for each phenotype by calculating AIC stepwise
### dormancy - not including cohort bc tightly correlated with dormancy (later cohorts = more dormancy)
model_all<-lm(days_to_radicle~hi+F21_238+F21_259+F21_241+F21_284+F21_232)
selectedMod <- step(model_all)
summary(selectedMod)
# pull selected covariates (plus hybrid_index)
#  formula = days_to_radicle ~ F21_259 + F21_241
# makes sense all families are retained here! Maternal effects huge in seeds
covariates_SD_select <- pull.pheno(Costus_prob, c("hybrid_index","F21_259","F21_241"))
                                                 
### growth rate
model_all<-lm(GR_cm.dy~hi+coh5+coh6+coh7+coh8+coh9+coh10+coh11+coh12+coh13+coh14+coh15+coh16+coh17+coh18+coh19+coh20+coh21+coh22+coh23+coh24+coh25+coh26+coh27+coh28+coh29+coh30+coh31+coh32+coh33+coh34+F21_238+F21_259+F21_241+F21_284+F21_232)
selectedMod <- step(model_all)
summary(selectedMod)
# pull selected covariates (plus hybrid_index)
# formula = GR_cm.dy ~ coh5 + coh6 + coh7 + coh8 + coh9 + coh10 + 
#   coh11 + coh12 + coh13 + coh14 + coh15 + coh16 + coh18 + coh19 + 
#   coh20 + coh21 + coh22 + coh23 + coh26 + coh27 + coh29 + coh30 + 
#   coh31 + F21_259
covariates_GR_select <- pull.pheno(Costus_prob, c("hybrid_index", "coh5","coh6","coh7","coh8","coh9","coh10","coh11","coh12","coh13","coh14","coh15","coh16","coh18","coh19","coh20","coh21","coh22","coh23","coh26","coh27","coh29","coh30","coh31","F21_259"))

### thickness (missing cohorts: coh7+coh34+coh1+coh8+coh33)
model_all<-lm(ave_thick~hi+coh1+coh5+coh6+coh9+coh10+coh11+coh12+coh13+coh14+coh15+coh16+coh17+coh18+coh19+coh20+coh21+coh22+coh23+coh24+coh25+coh26+coh27+coh28+coh29+coh30+coh31+coh32+F21_238+F21_259+F21_241+F21_284+F21_232)
selectedMod <- step(model_all)
summary(selectedMod)
# pull selected covariates (plus hybrid_index)
# formula = ave_thick ~ coh1 + coh5 + coh9 + coh10 + coh11 + 
#   coh12 + coh13 + coh14 + coh15 + coh16 + coh17 + coh18 + coh19 + 
#   coh20 + coh21 + coh22 + coh23 + coh24 + coh25 + coh26 + coh30 + 
#   F21_238 + F21_284
covariates_thick_select <- pull.pheno(Costus_prob, c("hybrid_index","coh1","coh5","coh9","coh10","coh11","coh12","coh13","coh14","coh15","coh16","coh17","coh18","coh19","coh20","coh21","coh22","coh23","coh24","coh25","coh26","coh30","F21_238","F21_284"))

# toughness (missing cohorts: coh7+coh34+coh1+coh8+coh33)
model_all<-lm(ave_tough~hi+coh1+coh5+coh6+coh9+coh10+coh11+coh12+coh13+coh14+coh15+coh16+coh17+coh18+coh19+coh20+coh21+coh22+coh23+coh24+coh25+coh26+coh27+coh28+coh29+coh30+coh31+coh32+F21_238+F21_259+F21_241+F21_284+F21_232)
selectedMod <- step(model_all)
summary(selectedMod)
# pull selected covariates (plus hybrid_index)
# formula = ave_tough ~ coh5 + coh6 + coh9 + coh10 + coh11 + 
#   coh13 + coh16 + coh18 + coh21 + coh23 + coh26 + coh31 + F21_238 + 
#   F21_259 + F21_284
covariates_tough_select <- pull.pheno(Costus_prob, c("hybrid_index","coh1","coh5","coh6","coh9","coh10","coh11","coh12","coh13","coh14","coh15","coh16","coh17","coh18","coh19","coh20","coh21","coh22","coh23","coh24","coh25","coh26","coh27","coh28","coh29","coh30","coh31","coh32","F21_238","F21_259","F21_241","F21_284","F21_232"))

# Percent N reabsorbed during drought
model_all<-lm(N_reabsorp~hi+coh1+coh5+coh6+coh7+coh8+coh9+coh10+coh11+coh12+coh13+coh14+coh15+coh16+coh17+coh18+coh19+coh20+coh21+coh22+coh23+coh24+coh25+coh26+coh27+coh28+coh29+coh30+coh31+coh32+coh33+coh34+F21_238+F21_259+F21_241+F21_284+F21_232)
selectedMod <- step(model_all)
summary(selectedMod)
# pull selected covariates (plus hybrid_index)
# formula = N_reabsorp ~ coh1 + coh5 + coh7 + coh8 + coh9 + 
#     coh10 + coh11 + coh12 + coh13 + coh14 + coh15 + coh16 + coh17 + 
#     coh18 + coh19 + coh20 + coh21 + coh22 + coh23 + coh24 + coh25 +  
#     coh26 + coh27 + coh28 + coh29 + coh30 + coh31 + coh32 + coh33 +  
#     F21_238 + F21_259 + F21_241 + F21_284
covariates_Nreabs_select <- pull.pheno(Costus_prob, c("hybrid_index", "coh1","coh5","coh7","coh8","coh9","coh10","coh11","coh12","coh13","coh14","coh15","coh16","coh17","coh18","coh19","coh20","coh21","coh22","coh23","coh24","coh25","coh26","coh27","coh28","coh29","coh30","coh31","coh32","coh33","F21_238","F21_259","F21_241","F21_284"))

# Percent chlorophyll reabsorbed during drought
model_all<-lm(chloro_reabsorp~hi+coh1+coh5+coh6+coh7+coh8+coh9+coh10+coh11+coh12+coh13+coh14+coh15+coh16+coh17+coh18+coh19+coh20+coh21+coh22+coh23+coh24+coh25+coh26+coh27+coh28+coh29+coh30+coh31+coh32+coh33+coh34+F21_238+F21_259+F21_241+F21_284+F21_232)
selectedMod <- step(model_all)
summary(selectedMod)
# pull selected covariates (plus hybrid_index)
# formula = chloro_reabsorp ~ coh5 + coh12 + coh13 + coh14 + 
#   coh15 + coh16 + coh19 + coh20 + coh21 + coh22 + coh24 + coh25 + 
#   coh27 + coh28 + coh30 + coh31 + coh32 + F21_259 + F21_241
covariates_Chlreabs_select <- pull.pheno(Costus_prob, c("hybrid_index","coh5","coh12","coh13","coh14","coh15","coh16","coh19","coh20","coh21","coh22","coh24","coh25","coh27","coh28","coh30","coh31","coh32","F21_259","F21_241"))

###### KK code - 2-dimensional scan #####
## CANT run scantwo locally, could try on humm but for now just do scanone - Error: vector memory exhausted (limit reached?)
#scantwo: Two-dimensional genome scan with a two-QTL model
#scantwo to do permutations for likelihood penalties
scn2_1kperm.hk <- scantwo(Costus_prob, method = "hk", pheno.col = "GR_cm.dy", addcovar = covariates_GRselect,  n.perm=1000, verbose = TRUE)
plot(scn2_1kperm.hk)

###### RUN SINGLE-QTL SCANS ######
#remove cohort covar to look at dormancy bc likely assoc with hybrid index (more vill acestry, slower germ, different cohort) which is already included. 
#- cohorts really fuck things up for germination, I think this is becasue it is overfitting the model since cohort is explaining the same variability as the germination value itself. 
# PICKUP - what I think I should do is give seed germination its own whole set of cohorts based on the batch of seeds (the b,c, d, e)
#covariates <- pull.pheno(Costus_prob, c("hybrid_index","F21_238","F21_259","F21_241","F21_284","F21_232"))

### Seed dormancy
seed.dormancy.scanone.hk <- scanone(Costus_prob, method="hk", pheno.col = "days_to_radicle", addcovar=covariates_SD_select)
plot(seed.dormancy.scanone.hk,lodcolumn=1, ylab="LOD score - seed dormancy")
summary(seed.dormancy.scanone.hk)
write.table(seed.dormancy.scanone.hk,file="seed.dormancy.scanone-hk_hotone_covar.tsv",sep="\t",quote=FALSE)
saveRDS(seed.dormancy.scanone.hk, file="seed.dormancy.scanone-hk_hotone_covar.rds")
seed.dormancy.perm.hk <- scanone(Costus_prob, method="hk", pheno.col ="days_to_radicle", n.perm=1000, addcovar=covariates_SD_select)
saveRDS(seed.dormancy.perm.hk, "seed.dormancy.scanone-hk_hotone_covar_1kperm.rds")
plot(seed.dormancy.perm.hk,lodcolumn=1)
summary(seed.dormancy.scanone.hk, perms=seed.dormancy.perm.hk, alpha=0.1, pvalues=TRUE)
# lod cutoffs
summary(seed.dormancy.perm.hk,alpha=c(0.05, 0.1, 0.2))
# LOD thresholds (1000 permutations)
# lod
# 5%  4.45
# 10% 3.82
# 20% 3.15
SD_10cutoff_lod = 3.74
summary(seed.dormancy.scanone.hk[seed.dormancy.scanone.hk$lod > SD_10cutoff_lod,])
write.table(summary(seed.dormancy.scanone.hk[seed.dormancy.scanone.hk$lod > SD_10cutoff_lod,]), "seed.dormancy_sig.peaks_scanone-hk_hotone_covar_1kperm.txt", quote = F)

#plot
pdf("./QTL_plots/seed.dormancy-qtl.pdf",24,7)
par(mar = c(5, 5, 4, 2))
plot(seed.dormancy.scanone.hk,lodcolumn=1, 
     ylab="LOD score",
     main = "seed dormancy",
     cex.lab = 2,
     cex.main = 2.5)
abline(h=SD_10cutoff_lod, col="red", lwd=3)
dev.off()

### Growth rate
growth.rate.scanone.hk <- scanone(Costus_prob, method="hk", pheno.col = "GR_cm.dy", addcovar=covariates_GR_select)
plot(growth.rate.scanone.hk,lodcolumn=1, ylab="LOD score - growth rate")
summary(growth.rate.scanone.hk)
write.table(growth.rate.scanone.hk,file="growth.rate.scanone-hk_hotone_covar.tsv",sep="\t",quote=FALSE)
saveRDS(growth.rate.scanone.hk, file="growth.rate.scanone-hk_hotone_covar.rds")
growth.rate.perm.hk <- scanone(Costus_prob, method="hk", pheno.col ="GR_cm.dy", n.perm=1000, addcovar=covariates_GR_select)
saveRDS(growth.rate.perm.hk, "growth.rate.scanone-hk_hotone_covar_1kperm.rds") 
plot(growth.rate.perm.hk,lodcolumn=1, ylab="LOD score - growth rate")
summary(growth.rate.scanone.hk, perms=growth.rate.perm.hk, alpha=0.1, pvalues=TRUE)
# lod cutoffs
summary(growth.rate.perm.hk,alpha=c(0.05, 0.1, 0.2))
# LOD thresholds (1000 permutations)
#      lod
# 5%  3.91
# 10% 3.50
# 20% 3.12
GR_10cutoff_lod = 3.50
summary(growth.rate.scanone.hk[growth.rate.scanone.hk$lod > GR_10cutoff_lod,])
write.table(summary(growth.rate.scanone.hk[growth.rate.scanone.hk$lod > GR_10cutoff_lod,]), "growth.rate_sig.peaks_scanone-hk_hotone_covar_1kperm.txt", quote = F)

#plot
pdf("./QTL_plots/growth.rate-qtl.pdf",24,7)
par(mar = c(5, 5, 4, 2))
plot(growth.rate.scanone.hk,lodcolumn=1,
     ylab="LOD score",
     main = "growth rate",
     cex.lab = 2,
     cex.main = 2.5)
abline(h=GR_10cutoff_lod, col="red", lwd=3)
dev.off()


### Leaf thickness
thickness.scanone.hk <- scanone(Costus_prob, method="hk", pheno.col = "ave_thick", addcovar=covariates_thick_select)
plot(thickness.scanone.hk,lodcolumn=1, ylab="LOD score - leaf thickness")
summary(thickness.scanone.hk)
write.table(thickness.scanone.hk,file="thickness.scanone-hk_hotone_covar.tsv",sep="\t",quote=FALSE)
saveRDS(thickness.scanone.hk, file="thickness.scanone-hk_hotone_covar.rds")
thickness.perm.hk <- scanone(Costus_prob, method="hk", pheno.col ="ave_thick", n.perm=1000, addcovar=covariates_thick_select)
saveRDS(thickness.perm.hk, "thickness.scanone-hk_hotone_covar_1kperm.rds")
plot(thickness.perm.hk,lodcolumn=1, ylab="LOD score - leaf thickness")
summary(thickness.scanone.hk, perms=thickness.perm.hk, alpha=0.1, pvalues=TRUE)
# lod cutoffs
summary(thickness.perm.hk,alpha=c(0.05, 0.1, 0.2))
# LOD thresholds (1000 permutations)
# lod
# 5%  3.87
# 10% 3.57
# 20% 3.16
thick_10cutoff_lod = 3.57
summary(thickness.scanone.hk[thickness.scanone.hk$lod > thick_10cutoff_lod,])
write.table(summary(thickness.scanone.hk[thickness.scanone.hk$lod > thick_10cutoff_lod,]), "leaf.thickness_sig.peaks_scanone-hk_hotone_covar_1kperm.txt", quote = F)

#plot
pdf("lf.thickness-qtl.pdf",24,7)
par(mar = c(5, 5, 4, 2))
plot(thickness.scanone.hk,lodcolumn=1, 
     ylab="LOD score",
     main = "leaf thickness",
     cex.lab = 2,
     cex.main = 2.5)
abline(h=thick_10cutoff_lod, col="red", lwd=3)
dev.off()

### Leaf toughness
toughness.scanone.hk <- scanone(Costus_prob, method="hk", pheno.col = "ave_tough", addcovar=covariates_tough_select)
plot(toughness.scanone.hk,lodcolumn=1, ylab="LOD score - leaf toughness")
summary(toughness.scanone.hk)
write.table(toughness.scanone.hk,file="toughness.scanone-hk_hotone_covar.tsv",sep="\t",quote=FALSE)
saveRDS(toughness.scanone.hk, file="toughness.scanone-hk_hotone_covar.rds")
toughness.perm.hk <- scanone(Costus_prob, method="hk", pheno.col ="ave_tough", n.perm=1000, addcovar=covariates_tough_select)
saveRDS(toughness.perm.hk, "toughness.scanone-hk_hotone_covar_1kperm.rds")
plot(toughness.perm.hk,lodcolumn=1, ylab="LOD score - leaf toughness")
summary(toughness.scanone.hk, perms=toughness.perm.hk, alpha=0.1, pvalues=TRUE)
# lod cutoffs
summary(toughness.perm.hk,alpha=c(0.05, 0.1, 0.2))
# LOD thresholds (1000 permutations)
# lod
# 5%  4.11
# 10% 3.79
# 20% 3.28
tough_10cutoff_lod = 3.79
summary(toughness.scanone.hk[toughness.scanone.hk$lod > tough_10cutoff_lod,])
write.table(summary(toughness.scanone.hk[toughness.scanone.hk$lod > tough_10cutoff_lod,]), "leaf.toughness_sig.peaks_scanone-hk_hotone_covar_1kperm.txt", quote = F)

#plot
pdf("lf.toughness-qtl.pdf",24,7)
par(mar = c(5, 5, 4, 2))
plot(toughness.scanone.hk,lodcolumn=1, 
     ylab="LOD score",
     main = "leaf toughness",
     cex.lab = 2,
     cex.main = 2.5)
abline(h=tough_10cutoff_lod, col="red", lwd=3)
dev.off()


### Percent N reabsorbed during drought
Nreabs.scanone.hk <- scanone(Costus_prob, method="hk", pheno.col = "N_reabsorp", addcovar=covariates_Nreabs_select)
plot(Nreabs.scanone.hk,lodcolumn=1, ylab="LOD score - %N reabsorbed in drought")
summary(Nreabs.scanone.hk)
write.table(Nreabs.scanone.hk,file="Nreabs.scanone-hk_hotone_covar.tsv",sep="\t",quote=FALSE)
saveRDS(Nreabs.scanone.hk, file="Nreabs.scanone-hk_hotone_covar.rds")
Nreabs.perm.hk <- scanone(Costus_prob, method="hk", pheno.col ="N_reabsorp", n.perm=1000, addcovar=covariates_Nreabs_select) 
saveRDS(Nreabs.perm.hk, "Nreabs.scanone-hk_hotone_covar_1kperm.rds")
plot(Nreabs.perm.hk,lodcolumn=1, ylab="LOD score - leaf N reabsorption")
summary(Nreabs.scanone.hk, perms=Nreabs.perm.hk, alpha=0.1, pvalues=TRUE)
# lod cutoffs
summary(Nreabs.perm.hk,alpha=c(0.05, 0.1, 0.2))
# LOD thresholds (1000 permutations)
# lod
# 5%  3.85
# 10% 3.49
# 20% 3.12
Nreabs_10cutoff_lod =3.49
summary(Nreabs.scanone.hk[Nreabs.scanone.hk$lod > Nreabs_10cutoff_lod,])
write.table(summary(Nreabs.scanone.hk[Nreabs.scanone.hk$lod > Nreabs_10cutoff_lod,]), "N.reabsorption_sig.peaks_scanone-hk_hotone_covar_1kperm.txt", quote = F)

#plot
pdf("Nreabs-qtl.pdf",24,7)
par(mar = c(5, 5, 4, 2))
plot(Nreabs.scanone.hk,lodcolumn=1, 
     ylab="LOD score",
     main = "% leaf nitrogen reabsorption",
     cex.lab = 2,
     cex.main = 2.5)
abline(h=Nreabs_10cutoff_lod, col="red", lwd=3)
dev.off()


### Percent chlorophyll reabsorbed during drought
Chlreabs.scanone.hk <- scanone(Costus_prob, method="hk", pheno.col = "chloro_reabsorp", addcovar=covariates_Chlreabs_select)
plot(Chlreabs.scanone.hk,lodcolumn=1, ylab="LOD score - Chlorophyll reabsorbed in drought")
summary(Chlreabs.scanone.hk)
write.table(Chlreabs.scanone.hk,file="Chlreabs.scanone-hk_hotone_covar.tsv",sep="\t",quote=FALSE)
saveRDS(Chlreabs.scanone.hk, file="Chlreabs.scanone-hk_hotone_covar.rds")
Chlreabs.perm.hk <- scanone(Costus_prob, method="hk", pheno.col ="chloro_reabsorp", n.perm=1000, addcovar=covariates_Chlreabs_select)
saveRDS(Chlreabs.perm.hk, "Chlreabs.scanone-hk_hotone_covar_1kperm.rds")
plot(Chlreabs.perm.hk,lodcolumn=1, ylab="LOD score - leaf Chlorophyll reabsorption")
summary(Chlreabs.scanone.hk, perms=Chlreabs.perm.hk, alpha=0.1, pvalues=TRUE)
# lod cutoffs
summary(Chlreabs.perm.hk,alpha=c(0.05, 0.1, 0.2))
# LOD thresholds (1000 permutations)
# lod
# 5%  3.99
# 10% 3.49
# 20% 3.08
Chlreabs_10cutoff_lod = 3.49
summary(Chlreabs.scanone.hk[Chlreabs.scanone.hk$lod > Chlreabs_10cutoff_lod,])
write.table(summary(Chlreabs.scanone.hk[Chlreabs.scanone.hk$lod > Chlreabs_10cutoff_lod,]), "chloro.reabsorption_sig.peaks_scanone-hk_hotone_covar_1kperm.txt", quote = F)

#plot
pdf("Chlreabs-qtl.pdf",24,7)
par(mar = c(5, 5, 4, 2))
plot(Chlreabs.scanone.hk,lodcolumn=1, 
     ylab="LOD score",
     main = "% chlorophyll reabsorption",
     cex.lab = 2,
     cex.main = 2.5)
abline(h=Chlreabs_10cutoff_lod, col="red", lwd=3)
dev.off()


###### LOD peak support ######
## get 95% Bayes credible interval

# seed dormancy
bayesint(seed.dormancy.scanone.hk, 0.95,expandtomarkers = TRUE)

# growth rate
bayesint(growth.rate.scanone.hk, 0.95,expandtomarkers = TRUE)

# leaf thickness
bayesint(thickness.scanone.hk, 0.95,expandtomarkers = TRUE)

# leaf toughness
bayesint(toughness.scanone.hk, 0.95,expandtomarkers = TRUE)

# N reabsorption
bayesint(Nreabs.scanone.hk, 0.95,expandtomarkers = TRUE)

# chlorophyll reabsorption
bayesint(Chlreabs.scanone.hk, 0.95,expandtomarkers = TRUE)

#out.boot <- scanoneboot(data_prob,chr="ScyDAA6-2113-HRSCAF-2539",n.boot=1000,prob=0.95)
#summary(out.boot)
#plot(out.boot)


# ###### MAKE EFFECT PLOTS ######
# 
# ##chr22 qtl effect plot
# setwd("~/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/LTREB_qtl/CTmax-QTL_first-scan_data")
# data_prob.df <- read.csv("ltreb-ctmax-qtl-filtered_22-1hot-site-tank_data-prob.csv")
# qtl_genos <- data_prob.df$ScyDAA6.2113.HRSCAF.2539.8811359
# qtl_genos[qtl_genos == "-"] <- NA
# ctmax <- data_prob.df$ctmax
# qtl22_peak_data <- data.frame(ctmax=ctmax,qtl22_geno=qtl_genos)
# qtl22_peak_data <- na.omit(qtl22_peak_data)
# 
# BB<-qtl22_peak_data[qtl22_peak_data$qtl22_geno=="BB",]$ctmax
# MB<-qtl22_peak_data[qtl22_peak_data$qtl22_geno=="AB",]$ctmax
# MM<-qtl22_peak_data[qtl22_peak_data$qtl22_geno=="AA",]$ctmax
# 
# malcol_22=rgb(0/255,0/255,175/255)
# hetcol_22=rgb(100/255,0/255,175/255)
# bircol_22=rgb(150/255,0/255,0/255)
# 
# ## error bar function (plot 1 SD)
# error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
#   if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
#     stop("vectors must be same length")
#   arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
# }
# 
# se <- function(x) 2*(sd(x)/sqrt(length(x)))
# 
# pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl-effect-plot-w-points_2se.pdf",5.5,7)
# 
# plot(1:3,c(mean(MM),mean(MB),mean(BB)),col=c(malcol_22,hetcol_22,bircol_22),ylim=c(32,38),xlab="Genotype",ylab="Normalized counts",xaxt="n",pch="-",cex=3,xlim=c(0.5,3.5))
# error.bar(1:3,c(mean(MM),mean(MB),mean(BB)),c(se(MM),se(MB),se(BB)),col=c(malcol_22,hetcol_22,bircol_22),lwd=2)
# 
# malcol=rgb(0/255,0/255,175/255,alpha=0.3)
# hetcol=rgb(100/255,0/255,175/255,alpha=0.6)
# bircol=rgb(150/255,0/255,0/255,alpha=0.6)
# 
# noise<-runif(length(BB),0.2,0.35)
# points(rep(3,length(BB))+noise,BB,pch=20,cex=1.8,col=bircol)
# 
# noise<-runif(length(MB),0.2,0.35)
# points(rep(2,length(MB))+noise,MB,pch=20,cex=1.8,col=hetcol)
# 
# noise<-runif(length(MM),0.2,0.35)
# points(rep(1,length(MM))+noise,MM,pch=20,cex=1.8,col=malcol)
# 
# mtext(c("MM","MB","BB"),at=1:3,side=1)
# 
# dev.off()
# 
# 
# ### other effect plots
# ## Load data_prob object
# data_prob <- readRDS("ltreb-ctmax-qtl-filtered_22-1hot-site-tank_data-prob.rds")
# 
# ## impute missing genotypes
# data_sim <- sim.geno(data_prob,step=1,n.draws=16)
# 
# ## plot estimated effect of QTL: phenotype vs marker genotype
# malcol=rgb(0/255,0/255,175/255)
# hetcol=rgb(100/255,0/255,175/255)
# bircol=rgb(150/255,0/255,0/255)
# effect_colors = c(malcol,hetcol,bircol)
# 
# ## plot estimated effect of QTL: phenotype vs marker genotype
# # 1st scan
# mar_chr22 <- find.marker(data_prob, chr='ScyDAA6-2113-HRSCAF-2539', pos=8.811359)
# # simple effect plot
# pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl-effect-plot.pdf",8.4,8.9)
# effectplot(data_sim,mname1=mar_chr22,main='22@8.81')
# dev.off()
# # effect plots with individuals as points
# pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-qtl-effect-plot-w-points.pdf",5.5,7)
# plotPXG(data_prob,marker=mar_chr22,main='22@8.81',infer=FALSE, ylab = expression('CT'['max']*' ('*~degree*C*')'), col = effect_colors)
# dev.off()
# 
# ## get QTL effect on CTmax per genotype at peak
# eff <- effectplot(data_sim,mname1=mar_chr22)
# eff
# # ScyDAA6-2113-HRSCAF-2539:8811359
# #      mean      SE
# #  AA: 35.98525, 0.1777719
# #  AB: 35.65896, 0.1180658
# #  BB: 36.00340, 0.1713042
# 
# ## plot estimated QTL effects along chr22 (shows additive and dominance effects)
# # effect summary for 22@8.81:
# # a           d             se.a        se.d
# # 0.008494465	-0.337265149	0.124726314	0.172637572
# effect_table<-effectscan(data_sim, pheno.col="ctmax", c("ScyDAA6-2113-HRSCAF-2539"), draw=F, get.se=T)
# write.table(effect_table,file="ltreb-only-CTmax_17-1hot-site-tanks.chr22-qtl_effect-table.tsv",sep="\t",quote=FALSE)
# # plot effect across chr22
# pdf("ltreb-only-CTmax_17-1hot-site-tanks_chr22-additive-dominance-effects-plot.pdf",9.2,6.8)
# effectscan(data_sim, pheno.col="ctmax", c("ScyDAA6-2113-HRSCAF-2539"), draw=T, get.se=T)
# dev.off()




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


