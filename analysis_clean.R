source("ggcoxzph_fixed.R")

#define helper function#
values2deciles <- function(values){
  qua <- quantile(values, seq(0,0.9,0.1))
  res <- rep(0, length(values))
  for (cu in qua){
    idx <- values >= cu
    res[idx] <- res[idx] + 1
  }
  res <- res - 1
  return(res)
}

if (exists("adnimerge"))
  rm(adnimerge)

library(data.table)
library(ADNIMERGE)
library(survival)
library(survminer)
library(cutpointr)

### start options ###
#biomarker to use
csf_only <- F
pet_only <- F

#add flag about population structure
adjust_for_PC <- F

#use bellenguez instead of Kunkle et al
use_bellenguez <- F

#which PRS cutoff to use
use_prs <- "PRS_1e.08"
#use_prs <- "PRS_5e.08"
#use_prs <- "PRS_1e.04"
#use_prs <- "PRS_0.001"

### end options ###

###### PRECOMPUTED VALUES ##############################
#mean and sd of the tau meta ROI computed in CN subjects
tau_pet_cn_mean <- 1.196489
tau_pet_cn_sd   <- 0.1132466

csf_ptau_cn_mean <- 21.79055
csf_ptau_cn_sd   <- 9.171088

########################################################
## load tables ##

polygen <- read.csv("ADNI1_3_genetic_summary_CEU80_rel.csv")
if (use_bellenguez){
  #replce Kunkle PRS with Bellenguez PRS
  polygen <- polygen[,1:33]
  bellen  <- read.table("ADNI1-3_bellenguez_avgPRS_mod2.all.score", head=T)
  colnames(bellen) <- sub("X","PRS_",colnames(bellen))
  bellen <- bellen[,2:17]
  polygen <- merge(polygen, bellen, by.x="PTID",by.y="IID")
}

taupet  <- read.csv("UCBERKELEYAV1451_8mm_02_17_23_20Jul2023.csv")
amypet  <- read.csv("UCBERKELEYAV45_8mm_02_17_23_20Jul2023.csv")
amypet2  <- read.csv("UCBERKELEYFBB_8mm_02_17_23_20Jul2023.csv")

#load clean CSF biomarker values (prepared using 'prepare_csf_clean.R')
roche_csf <- read.csv("ROCHE_ELECSYS_biomk.csv")

#merge into adnimerge
adnimerge_csf <- merge(adnimerge, roche_csf, by.x=c("RID","VISCODE"), by.y=c("RID","VISCODE2"), suffixes = c("","_ROCHE"), all.x=TRUE)

########################################################

############# define amyloid and tau cutoffs ###########

### amyloid cutoffs ###
##hard defined cutoffs (csf)
Acut <- 880      #see cited publication

### tau cutoffs ###
#z-score thresholds
zcuts <- seq(1,3,0.25)

#tau-pet cutoffs based on mean and sd
tau_cuts <- zcuts * tau_pet_cn_sd + tau_pet_cn_mean

#define csf ptau cutoffs based on youden's index
csf_ocp2 <- c()
tmp_data <- merge(roche_csf, taupet, by=c("RID","VISCODE2"))
for (sss in tau_cuts){
  tmp_data[["PETpos"]] <- tmp_data$META_TEMPORAL_SUVR > sss
  ocp2 <- cutpointr(subset(tmp_data, !is.na(META_TEMPORAL_SUVR)), PTAU, PETpos, method=oc_youden_normal, na.rm=T)$optimal_cutpoint
  csf_ocp2 <- c(csf_ocp2, ocp2)
}
Tcuts <- round(csf_ocp2,2)

## use for sensitivity analysis (only CSF)
#Tcuts <- 22:31

#############################################

#pick a single cutoff (based on the index, used for screening)
#if no index is given, tau pet cutoff at z-score 2 is used

#tcut_idx <- NA
if (!exists("tcut_idx")){
  #z-score of 2
  tcut_idx <- 5
}

Tcut <- Tcuts[tcut_idx]
tau_cut <- tau_cuts[tcut_idx]

print(Tcut)
print(tau_cut)

################# establish Amyloid and Tau status #######################

### deal with abeta status ###
adnimerge_abeta  <- merge(adnimerge_csf, amypet,  by.x=c("RID","VISCODE"), by.y=c("RID","VISCODE2"), all.x=TRUE)
adnimerge_abeta2 <- merge(adnimerge_csf, amypet2, by.x=c("RID","VISCODE"), by.y=c("RID","VISCODE2"), all.x=TRUE)

#define amyloid positive
rm(aplus)
aplus   <- cbind(adnimerge_abeta$SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF > 0 , adnimerge_abeta2$SUMMARYSUVR_WHOLECEREBNORM_1.08CUTOFF > 0)
aplus   <- cbind(aplus, as.double(gsub(">","",adnimerge_abeta$ABETA42)) < Acut)

## switch to only CSF or only PET
if (csf_only){
  aplus[!is.na(aplus[,1]),1] <- NA
  aplus[!is.na(aplus[,2]),2] <- NA
}
if (pet_only){
  aplus[!is.na(aplus[,3]),3] <- NA
}

#define amyloid status for each visit (missing, + or -)
astatus <- apply(aplus,1,function(x){
  if (sum(!is.na(x))==0)
    return("Ana")
  if (sum(x, na.rm=T)>0){
    return("A+")
   } else {
    return("A-")
  }
})


### deal with tau status ###
adnimerge_tau <- merge(adnimerge_csf, taupet, by.x=c("RID","VISCODE"), by.y=c("RID","VISCODE2"), all.x=TRUE)
csf_tau <- as.double(gsub(">", "", gsub("<","",adnimerge_tau$PTAU_ROCHE)))

#define tau positive
rm(tplus)
tplus   <- (adnimerge_tau$META_TEMPORAL_SUVR / adnimerge_tau$INFERIORCEREBELLUM_SUVR) > tau_cut
tplus   <-  cbind(tplus, csf_tau > Tcut)

## switch to only CSF or only PET
if (csf_only){
  tplus[!is.na(tplus[,1]),1] <- NA
}
if (pet_only){
  tplus[!is.na(tplus[,2]),2] <- NA
}

tstatus <- apply(tplus,1,function(x){
  if (sum(!is.na(x))>0){
    if (sum(x, na.rm=T)>0){
      return("T+")
    } else {
      return("T-")
    }
  } else {
    return("Tna")
  }
})

#######################################
### add A and T status to adnimerge ###

colnames(adnimerge_abeta)[colnames(adnimerge_abeta)=="EXAMDATE.x"] <- "EXAMDATE"
adnimerge_bak <- adnimerge_abeta[,colnames(adnimerge)]

adnimerge_bak$A <- astatus
adnimerge_bak$T <- tstatus

ATN <- apply(cbind(astatus, tstatus),1,paste,collapse="")
ATN[ATN=="AnaTna"] <- NA
adnimerge_bak$ATN <- ATN


adnimerge_bak[adnimerge_bak$A == "Ana","A"] <- NA
adnimerge_bak[adnimerge_bak$T == "Tna","T"] <- NA

####################################################
#convert PRS to deciles ##
##get norm factors for PRS
bdat <- subset(adnimerge_bak, Years.bl==0)
bdat_gen <- merge(bdat, polygen, by="RID")
bdat_gen <- subset(bdat_gen, subset=(CEU>0.8) & (EXCLUDE==0))
decile <- values2deciles(bdat_gen[,use_prs])
bdat_gen$prs_decile <- decile
polygen2 <- merge(polygen, bdat_gen[,c("RID","prs_decile")], all.x=TRUE)
####################################################


### extract participants with longitudinal data and at least 2 visits with non-NA ATN status ###
longrid <- as.integer(names(which(table(subset(adnimerge_bak, !is.na(ATN))$RID) > 1)))

longATN <- subset(adnimerge_bak, is.element(RID, longrid))
