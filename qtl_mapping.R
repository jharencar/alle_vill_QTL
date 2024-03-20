## We used the scanone function of R/qtl to perform single QTL model standard 
## interval mapping using the EM algorithm.33 Recombination fraction was estimated 
## using the est.rf() function and markers missing genotype data were excluded 
## using the drop.nullmarker() function

setwd("/Users/juliaharencar/Documents/Github/alle_vill_QTL/")

## clear workspace
rm(list = ls())
library(snow)
library(qtl)

#library(snow)
#library(qtl2)

########### F2s - for geno probs ###########
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

## look for duplicate markers (i.e., markers with identical genotypes)
dup <- findDupMarkers(Costus, exact.only=FALSE, adjacent.only=TRUE)
Costus <- drop.markers(Costus, unlist(dup))
# brings total markers 16046 --> 3805, % genotyped from 99.3% --> 97.2%

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


###### CALCULATE GENOTYPE PROBABILITIES ######

## load est.rf results after running on cluster/server
#data_sub.ds <- readRDS(file = "./ltreb-qtl-data_sub.ds-est.rf.rds")

## calculate conditional genotype probabilities given multipoint marker data
Costus_prob <- calc.genoprob(Costus)

# 1st scan
write.cross(Costus_prob,file="first_hi.cohort.fam.hotones_covars_AVmap-prob_drop.dup.markers.csv")

saveRDS(Costus_prob,file="first_hi.cohort.fam.hotones_covars_AVmap-prob-drop.dup.markers.rds")
#Costus_prob <- readRDS(file="first_hi.cohort.fam.hotones_covars_AVmap-prob-prob.rds")
#Costus_prob <- read.cross(file="first_hi.cohort.fam.hotones_covars_AVmap-prob.csv")
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
### dormancy ####
#- not including cohort bc tightly correlated with dormancy (later cohorts = more dormancy)
model_all<-lm(days_to_radicle~hi+F21_238+F21_259+F21_241+F21_284+F21_232)
selectedMod <- step(model_all)
summary(selectedMod)
# pull selected covariates (plus hybrid_index)
#  formula = days_to_radicle ~ F21_259 + F21_241
# makes sense all families are retained here! Maternal effects huge in seeds
covariates_SD_select <- pull.pheno(Costus_prob, c("hybrid_index","F21_259","F21_241"))

### growth rate ####
model_all<-lm(GR_cm.dy~hi+coh5+coh6+coh7+coh8+coh9+coh10+coh11+coh12+coh13+coh14+coh15+coh16+coh17+coh18+coh19+coh20+coh21+coh22+coh23+coh24+coh25+coh26+coh27+coh28+coh29+coh30+coh31+coh32+coh33+coh34+F21_238+F21_259+F21_241+F21_284+F21_232)
selectedMod <- step(model_all)
summary(selectedMod)
# pull selected covariates (plus hybrid_index)
# formula = GR_cm.dy ~ coh5 + coh6 + coh7 + coh8 + coh9 + coh10 + 
#   coh11 + coh12 + coh13 + coh14 + coh15 + coh16 + coh18 + coh19 + 
#   coh20 + coh21 + coh22 + coh23 + coh26 + coh27 + coh29 + coh30 + 
#   coh31 + F21_259
covariates_GR_select <- pull.pheno(Costus_prob, c("hybrid_index", "coh5","coh6","coh7","coh8","coh9","coh10","coh11","coh12","coh13","coh14","coh15","coh16","coh18","coh19","coh20","coh21","coh22","coh23","coh26","coh27","coh29","coh30","coh31","F21_259"))

### thickness ####
#(missing cohorts: coh7+coh34+coh1+coh8+coh33)
model_all<-lm(ave_thick~hi+coh1+coh5+coh6+coh9+coh10+coh11+coh12+coh13+coh14+coh15+coh16+coh17+coh18+coh19+coh20+coh21+coh22+coh23+coh24+coh25+coh26+coh27+coh28+coh29+coh30+coh31+coh32+F21_238+F21_259+F21_241+F21_284+F21_232)
selectedMod <- step(model_all)
summary(selectedMod)
# pull selected covariates (plus hybrid_index)
# formula = ave_thick ~ coh1 + coh5 + coh9 + coh10 + coh11 + 
#   coh12 + coh13 + coh14 + coh15 + coh16 + coh17 + coh18 + coh19 + 
#   coh20 + coh21 + coh22 + coh23 + coh24 + coh25 + coh26 + coh30 + 
#   F21_238 + F21_284
covariates_thick_select <- pull.pheno(Costus_prob, c("hybrid_index","coh1","coh5","coh9","coh10","coh11","coh12","coh13","coh14","coh15","coh16","coh17","coh18","coh19","coh20","coh21","coh22","coh23","coh24","coh25","coh26","coh30","F21_238","F21_284"))

# toughness ####
#(missing cohorts: coh7+coh34+coh1+coh8+coh33)
model_all<-lm(ave_tough~hi+coh1+coh5+coh6+coh9+coh10+coh11+coh12+coh13+coh14+coh15+coh16+coh17+coh18+coh19+coh20+coh21+coh22+coh23+coh24+coh25+coh26+coh27+coh28+coh29+coh30+coh31+coh32+F21_238+F21_259+F21_241+F21_284+F21_232)
selectedMod <- step(model_all)
summary(selectedMod)
# pull selected covariates (plus hybrid_index)
# formula = ave_tough ~ coh5 + coh6 + coh9 + coh10 + coh11 + 
#     coh13 + coh16 + coh18 + coh21 + coh23 + coh26 + coh31 + F21_238 + 
#     F21_259 + F21_284
covariates_tough_select <- pull.pheno(Costus_prob, c("hybrid_index","coh5","coh6","coh9","coh10","coh11","coh13","coh16","coh18","coh21","coh23","coh26","coh31","F21_238","F21_259","F21_284"))

# Percent N reabsorbed during drought ####
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

# Percent chlorophyll reabsorbed during drought ####
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
# scn2_1kperm.hk <- scantwo(Costus_prob, method = "hk", pheno.col = "GR_cm.dy", addcovar = covariates_GRselect,  n.perm=1000, verbose = TRUE)
# plot(scn2_1kperm.hk)

###### RUN SINGLE-QTL SCANS ######
#remove cohort covar to look at dormancy bc likely assoc with hybrid index (more vill acestry, slower germ, different cohort) which is already included. 
#- cohorts really fuck things up for germination, I think this is becasue it is overfitting the model since cohort is explaining the same variability as the germination value itself. 
# PICKUP - what I think I should do is give seed germination its own whole set of cohorts based on the batch of seeds (the b,c, d, e)
#covariates <- pull.pheno(Costus_prob, c("hybrid_index","F21_238","F21_259","F21_241","F21_284","F21_232"))

#### Seed dormancy ####
seed.dormancy.scanone.hk <- scanone(Costus_prob, method="hk", pheno.col = "days_to_radicle", addcovar=covariates_SD_select)
#  Dropping 8 individuals with missing phenotypes.
plot(seed.dormancy.scanone.hk,lodcolumn=1, ylab="LOD score - seed dormancy")
summary(seed.dormancy.scanone.hk)
write.table(seed.dormancy.scanone.hk,file="seed.dormancy.scanone-hk_hotone_covar.tsv",sep="\t",quote=FALSE)
saveRDS(seed.dormancy.scanone.hk, file="seed.dormancy.scanone-hk_hotone_covar.rds")
# permutations
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
#pdf("./QTL_plots/seed.dormancy-qtl.pdf",24,7)
#par(mar = c(5, 5, 4, 2))
#plot all
pdf("./QTL_plots/all_qtl.pdf",40,18)
par(mar = c(5, 5, 4, 2), mfrow = c(3,2))
plot(seed.dormancy.scanone.hk,lodcolumn=1, 
     ylab="LOD score",
     main = "seed dormancy",
     cex.lab = 2,
     cex.main = 2.5)
abline(h=SD_10cutoff_lod, col="red", lwd=3)
dev.off()

#try calculating effect sizes with scanone results:
seed.dormancy.perm.hk <- readRDS("seed.dormancy.scanone-hk_hotone_covar_1kperm.rds")


# # #try scantwo again 
# set.seed(3)
# require(snow)
# SD.scantwo.hk.set3 <- scantwo(Costus_prob, method="hk", pheno.col = "days_to_radicle", n.perm=10, n.cluster = 3, addcovar=covariates_SD_select)
# summary(SD.scantwo.hk)
# #calculate penalties to be used in stepwiseqtl based on scantwo permutations
# seed.dormancy_penalties <- calc.penalties(SD.scantwo.hk.set3, alpha = 0.05)
# 
# #this step took forever, so writing the results to csv files
# write.csv(seed.dormancy_penalties, file = "seed.dormancy_penalties.csv")
# 
# seed.dormancy_penalties <- read.csv(file = "seed.dormancy_penalties.csv", header = TRUE, row.names = 1)
# 
# #stepwiseqtl: Stepwise selection for multiple QTL
# #seed.dormancy
# seed.dormancy.qtl <- stepwiseqtl(Costus_prob, pheno.col = "days_to_radicle", method = "hk", model = "normal",
#                                  penalties = as.vector(Costus_penalties["days_to_radicle",1:3]), max.qtl = 8, keeplodprofile = TRUE,
#                                  verbose = TRUE, refine.locations = TRUE)
# summary(seed.dormancy.qtl)
# plot(seed.dormancy.qtl)
# #fitqtl: Fit a multiple QTL model
# seed.dormancy.fit <- fitqtl(Costus_prob, pheno.col = "days_to_radicle", method = "hk", model = "normal", get.ests = TRUE,
#                             qtl = seed.dormancy.qtl, formula = attr(seed.dormancy.qtl,"formula"), dropone = TRUE)
# summary(seed.dormancy.fit)


## figure out which of two versions bayes
# Calculate an approximate Bayes credible interval for location
bayesint(seed.dormancy.qtl, qtl.index = 1)
bayesint(APHqtl, qtl.index = 2)
bayesint(APHqtl, qtl.index = 3)
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



pdf(file="APH_LOD.pdf", width = 7.5, height = 2.3, pointsize = 8) 
# Plot 1-dimensional LOD profiles for a multiple QTL model
plotLodProfile(APHqtl, showallchr = TRUE, qtl.labels = FALSE, lwd=1, xlab = "", ylab = "", 
               main ="APH", frame.plot=FALSE)
dev.off()


#### Growth rate ####
growth.rate.scanone.hk <- scanone(Costus_prob, method="hk", pheno.col = "GR_cm.dy", addcovar=covariates_GR_select)
# Dropping 14 individuals with missing phenotypes.
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


#### Leaf thickness ####
thickness.scanone.hk <- scanone(Costus_prob, method="hk", pheno.col = "ave_thick", addcovar=covariates_thick_select)
# Dropping 46 individuals with missing phenotypes.
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
pdf("./QTL_plots/lf.thickness-qtl.pdf",24,7)
par(mar = c(5, 5, 4, 2))
plot(thickness.scanone.hk,lodcolumn=1, 
     ylab="LOD score",
     main = "leaf thickness",
     cex.lab = 2,
     cex.main = 2.5)
abline(h=thick_10cutoff_lod, col="red", lwd=3)
dev.off()

#### Leaf toughness ####
toughness.scanone.hk <- scanone(Costus_prob, method="hk", pheno.col = "ave_tough", addcovar=covariates_tough_select)
#  Dropping 46 individuals with missing phenotypes.
# addcovar appears to be over-specified; consider dropping columns
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
#      lod
# 5%  3.77
# 10% 3.44
# 20% 3.07
tough_10cutoff_lod = 3.44
summary(toughness.scanone.hk[toughness.scanone.hk$lod > tough_10cutoff_lod,])
write.table(summary(toughness.scanone.hk[toughness.scanone.hk$lod > tough_10cutoff_lod,]), "leaf.toughness_sig.peaks_scanone-hk_hotone_covar_1kperm.txt", quote = F)

#plot
pdf("./QTL_plots/lf.toughness-qtl.pdf",24,7)
par(mar = c(5, 5, 4, 2))
plot(toughness.scanone.hk,lodcolumn=1, 
     ylab="LOD score",
     main = "leaf toughness",
     cex.lab = 2,
     cex.main = 2.5)
abline(h=tough_10cutoff_lod, col="red", lwd=3)
dev.off()


#### Percent N reabsorbed during drought ####
Nreabs.scanone.hk <- scanone(Costus_prob, method="hk", pheno.col = "N_reabsorp", addcovar=covariates_Nreabs_select)
# Dropping 45 individuals with missing phenotypes.
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
#      lod
# 5%  4.14
# 10% 3.75
# 20% 3.30
Nreabs_10cutoff_lod =3.75
summary(Nreabs.scanone.hk[Nreabs.scanone.hk$lod > Nreabs_10cutoff_lod,])
write.table(summary(Nreabs.scanone.hk[Nreabs.scanone.hk$lod > Nreabs_10cutoff_lod,]), "N.reabsorption_sig.peaks_scanone-hk_hotone_covar_1kperm.txt", quote = F)

#plot
pdf("./QTL_plots/Nreabs-qtl.pdf",24,7)
par(mar = c(5, 5, 4, 2))
plot(Nreabs.scanone.hk,lodcolumn=1, 
     ylab="LOD score",
     main = "% leaf nitrogen reabsorption",
     cex.lab = 2,
     cex.main = 2.5)
abline(h=Nreabs_10cutoff_lod, col="red", lwd=3)
dev.off()


#### Percent chlorophyll reabsorbed during drought ####
Chlreabs.scanone.hk <- scanone(Costus_prob, method="hk", pheno.col = "chloro_reabsorp", addcovar=covariates_Chlreabs_select)
#  Dropping 12 individuals with missing phenotypes.
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
#      lod
# 5%  3.94
# 10% 3.58
# 20% 3.13
Chlreabs_10cutoff_lod = 3.58
summary(Chlreabs.scanone.hk[Chlreabs.scanone.hk$lod > Chlreabs_10cutoff_lod,])
write.table(summary(Chlreabs.scanone.hk[Chlreabs.scanone.hk$lod > Chlreabs_10cutoff_lod,]), "chloro.reabsorption_sig.peaks_scanone-hk_hotone_covar_1kperm.txt", quote = F)

#plot
pdf("./QTL_plots/Chlreabs-qtl.pdf",24,7)
par(mar = c(5, 5, 4, 2))
plot(Chlreabs.scanone.hk,lodcolumn=1, 
     ylab="LOD score",
     main = "% chlorophyll reabsorption",
     cex.lab = 2,
     cex.main = 2.5)
abline(h=Chlreabs_10cutoff_lod, col="red", lwd=3)
dev.off()

###### MULTI-QTL ANALYSIS ######
#### Seed dormancy ####
# create QTL object with both loci
# what="prob" will out the QTL genotype probabilities for use in HK regression
SD.qtl <- makeqtl(Costus_prob, chr=c("2","6"), pos=c(23.863645, 118.34823), what="prob")

# set covariates
covariates_SD_select <- pull.pheno(Costus_prob, c("hybrid_index","F21_259","F21_241"))

# fit the two locus additive model
# “drop one term at a time” table compares the fit of the two-QTL model to the reduced model where one QTL is omitted.
SD.out.fq <- fitqtl(Costus_prob, pheno.col="days_to_radicle", qtl=SD.qtl, method="hk", get.ests=TRUE, covar=covariates_SD_select, formula=y~Q1+Q2 + hybrid_index + F21_259 + F21_241)
summary(SD.out.fq)

# determine whether there is an interaction between the two QTL by fitting the model with the interaction
out.fqi <- fitqtl(Costus_prob, pheno.col="days_to_radicle", qtl=SD.qtl, method="hk", get.ests=TRUE, covar=covariates_SD_select, formula=y~Q1+Q2+Q1:Q2 + hybrid_index + F21_259 + F21_241)
summary(SD.out.fqi)


#### Growth rate ####
# create QTL object with both loci
# what="prob" will out the QTL genotype probabilities for use in HK regression
GR.qtl <- makeqtl(Costus_prob, chr="3", pos=11.007866, what="prob")

# set covariates
covariates_SD_select <- pull.pheno(Costus_prob, c("hybrid_index","F21_259","F21_241"))

# fit the two locus additive model
# “drop one term at a time” table compares the fit of the two-QTL model to the reduced model where one QTL is omitted.
GR.out.fq <- fitqtl(Costus_prob, pheno.col="GR_cm.dy", qtl=GR.qtl, method="hk", get.ests=TRUE, covar=covariates_GR_select, formula=y~Q1 + hybrid_index + coh5 + coh6 + coh7 + coh8 + coh9 + coh10 + coh11 + coh12 + coh13 + coh14 + coh15 + coh16 + coh18 + coh19 + coh20 + coh21 + coh22 + coh23 + coh26 + coh27 + coh29 + coh30 + coh31 + F21_259)
summary(GR.out.fq)

#### toughness ####
# create QTL object
# what="prob" will out the QTL genotype probabilities for use in HK regression
tough.qtl <- makeqtl(Costus_prob, chr="3", pos=11.845546, what="prob")

# set covariates
covariates_tough_select <- pull.pheno(Costus_prob, c("hybrid_index","coh5","coh6","coh9","coh10","coh11","coh13","coh16","coh18","coh21","coh23","coh26","coh31","F21_238","F21_259","F21_284"))

# fit qtl model
# “drop one term at a time” table compares the fit of the two-QTL model to the reduced model where one QTL is omitted.
tough.out.fq <- fitqtl(Costus_prob, pheno.col="ave_tough", qtl=tough.qtl, method="hk", get.ests=TRUE, covar=covariates_tough_select, formula=y~Q1 + hybrid_index + coh5 + coh6 + coh9 + coh10 + coh11 + coh13 + coh16 + coh18 + coh21 + coh23 + coh26 + coh31 + F21_238 + F21_259 + F21_284)
summary(tough.out.fq)

#### Percent N reabsorbed during drought ####
# create QTL object
# what="prob" will out the QTL genotype probabilities for use in HK regression
Nreabs.qtl <- makeqtl(Costus_prob, chr="1", pos=83.412886, what="prob")

# set covariates
covariates_Nreabs_select <- pull.pheno(Costus_prob, c("hybrid_index", "coh1","coh5","coh7","coh8","coh9","coh10","coh11","coh12","coh13","coh14","coh15","coh16","coh17","coh18","coh19","coh20","coh21","coh22","coh23","coh24","coh25","coh26","coh27","coh28","coh29","coh30","coh31","coh32","coh33","F21_238","F21_259","F21_241","F21_284"))

# fit the two locus additive model
# “drop one term at a time” table compares the fit of the two-QTL model to the reduced model where one QTL is omitted.
Nreabs.out.fq <- fitqtl(Costus_prob, pheno.col="N_reabsorp", qtl=Nreabs.qtl, method="hk", get.ests=TRUE, covar=covariates_Nreabs_select, formula=y~Q1 + hybrid_index + coh1 + coh5 + coh7 + coh8 + coh9 + coh10 + coh11 + coh12 + coh13 + coh14 + coh15 + coh16 + coh17 + coh18 + coh19 + coh20 + coh21 + coh22 + coh23 + coh24 + coh25 +  coh26 + coh27 + coh28 + coh29 + coh30 + coh31 + coh32 + coh33 +  F21_238 + F21_259 + F21_241 + F21_284)
summary(Nreabs.out.fq)

#### Percent chlorophyll reabsorbed during drought ####
Chlreabs.scanone.hk <- scanone(Costus_prob, method="hk", pheno.col = "chloro_reabsorp", addcovar=covariates_Chlreabs_select)

# create QTL object with both loci
# what="prob" will out the QTL genotype probabilities for use in HK regression
Chlreabs.qtl <- makeqtl(Costus_prob, chr="1", pos=71.772467, what="prob")

# set covariates
covariates_Chlreabs_select <- pull.pheno(Costus_prob, c("hybrid_index","coh5","coh12","coh13","coh14","coh15","coh16","coh19","coh20","coh21","coh22","coh24","coh25","coh27","coh28","coh30","coh31","coh32","F21_259","F21_241"))

# fit the two locus additive model
# “drop one term at a time” table compares the fit of the two-QTL model to the reduced model where one QTL is omitted.
Chlreabs.out.fq <- fitqtl(Costus_prob, pheno.col="N_reabsorp", qtl=Chlreabs.qtl, method="hk", get.ests=TRUE, covar=covariates_Chlreabs_select, formula=y~Q1 + hybrid_index + coh5 + coh12 + coh13 + coh14 + coh15 + coh16 + coh19 + coh20 + coh21 + coh22 + coh24 + coh25 + coh27 + coh28 + coh30 + coh31 + coh32 + F21_259 + F21_241)
summary(Chlreabs.out.fq)

# # refine the qtl positions
# # each QTL is moved to the position giving the highest likelihood,
# # and the entire process is repeated until no further improvement in likelihood can be obtained
# rqtl <- refineqtl(Costus_prob, pheno.col="days_to_radicle", qtl=qtl, method="hk")
# rqtl
# #       name chr     pos n.gen
# # Q1  2@23.8   2  23.809     3
# # Q2 6@118.2   6 118.165     3
# 
# # look for additional qtl, using the refined qtl positions
# out.aq_2qtl <- addqtl(Costus_prob, pheno.col="days_to_radicle", qtl=rqtl, method="hk", covar=covariates_SD_select, formula=y~Q1+Q2+Q1:Q2 + hybrid_index + F21_259 + F21_241)
# 
# plot(out.aq_2qtl)
# summary(out.aq_2qtl)
# out.aq_2qtl[out.aq_2qtl$lod > 3,] #terrible cutoff for me... need something way higher but need to figure out how... 
