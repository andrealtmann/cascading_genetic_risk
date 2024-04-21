### models the A+T- to (A+)T+ conversion ###
# run analysis_clean.R before!
############################################

print_figures <- F
sensitivity <- F

##extract participants with A+T- profile
tmp <- unique(subset(longATN, ATN=="A+T-")$RID)
longApTm <- subset(longATN, is.element(RID, tmp))

#build data for survival
jjj <- sapply(unique(longApTm$RID), function(rid){
  tmp <- subset(longApTm, RID==rid)
  tmp <- tmp[order(tmp$Years.bl),]

  #earliest entry with A-T-
  idx1 <- min(which(tmp$ATN == "A+T-"))
  #latest entry with T
  if (sensitivity){
    idx2 <- max(which(!is.na(tmp$T)))
  } else {
    idx2 <- max(which(!is.na(tmp$A) & !is.na(tmp$T)))
  }
  #get 'conversion'
  if (sensitivity){
    idx3 <- grep("T\\+", tmp$ATN)
  } else {
    idx3 <- min(which(tmp$ATN == "A+T+"))
  }
  idx3 <- min(idx3[idx3 >= idx1])

  tc <- NA
  if (length(idx3)>0)
    tc <- tmp[idx3, "Years.bl"]

  res <- c(rid, tmp[ idx1, "Years.bl"], tmp[idx2, "Years.bl"], tc )

})

#make a data frame
long.base2 <- subset(longApTm, Years.bl==0)
dummy2.df <- merge(long.base2, t(jjj), by.x="RID", by.y="V1")
dummy2.df <- merge(dummy2.df, polygen2, by="RID")
dummy2.df <- subset(dummy2.df, EXCLUDE==0 & CEU>0.8)

#remove cases where end is before start
valid <- dummy2.df$V3
resp <- !is.na(dummy2.df$V4)
valid[resp] <- dummy2.df[resp, "V4"]
valid <- valid - dummy2.df$V2

dummy2.df <- subset(dummy2.df, valid > 0 & !is.na(rs7412_T))


zero <- rep(0, nrow(dummy2.df))
start <- dummy2.df$V2
stop <- dummy2.df$V3
resp <- !is.na(dummy2.df$V4)
print(table(resp))

stop[resp] <- dummy2.df[resp, "V4"]
stop <- stop - start
mysurv2 <- Surv(zero, stop, resp)

dummy2.df_short <- dummy2.df[,c("PTGENDER","rs7412_T","APOE4","PTEDUCAT","PC1","PC2","PC3","PC4","PC5",use_prs)]
dummy2.df_short[["AGE"]] <- dummy2.df$AGE + dummy2.df$V2
colnames(dummy2.df_short) <- c("Sex","APOEe2","APOEe4","Edu","PC1","PC2","PC3","PC4","PC5","PRS","Age")

#turn the selected PRS into a z-score
dummy2.df_short$PRS <- (dummy2.df_short$PRS - mean(dummy2.df_short$PRS))/sd(dummy2.df_short$PRS)

#dummy2.df_short$APOEe4 <- as.factor(dummy2.df_short$APOEe4)

if (adjust_for_PC){
  mm <- coxph(mysurv2 ~ Age + Sex + Edu + PC1 + PC2 + PC3 + PC4 + PC5 + APOEe2 + APOEe4 + PRS, data=dummy2.df_short)
} else {
  mm <- coxph(mysurv2 ~ Age + Sex + Edu + APOEe2 + APOEe4 + PRS, data=dummy2.df_short)
}

mm.ph <- cox.zph(mm)
print(summary(mm))

if (print_figures){
  pdf("~/work/paper/inprogress/Stage_Varying_Risk/Figure2_Tau.pdf", width=4.5)
  ggforest(mm, data=dummy2.df_short)
  dev.off()

  pdf("~/work/paper/revision/Stage_Varying_Risk/R1/FigureS3_Tau_ph_20240328.pdf", width=9)
  ggcoxzph(mm.ph)
  dev.off()
}


##concordance
# 81.20 % concordance between CSF and PET for tau

### computing C-index ###
if (0){
  library(rms)
  dummy2.df_short$APOEe4 <- as.integer(dummy2.df_short$APOEe4)

  mymod  <- cph(mysurv2 ~ Age + Sex + Edu + APOEe2 + APOEe4 + PRS, data=dummy2.df_short, x=T, y=T)
  suppressWarnings(myboot <- validate(mymod, B=1000))
  Cind <- myboot["Dxy","index.corrected"] / 2 + 0.5
  print( c("Full:", Cind))

  mymod  <- cph(mysurv2 ~ Age + Sex + Edu + PRS, data=dummy2.df_short, x=T, y=T)
  suppressWarnings(myboot <- validate(mymod, B=1000))
  Cind <- myboot["Dxy","index.corrected"] / 2 + 0.5
  print( c("No-APOE:", Cind))

  mymod  <- cph(mysurv2 ~ Age + Sex + Edu + APOEe2 + APOEe4, data=dummy2.df_short, x=T, y=T)
  suppressWarnings(myboot <- validate(mymod, B=1000))
  Cind <- myboot["Dxy","index.corrected"] / 2 + 0.5
  print( c("No-PRS:", Cind))

  mymod  <- cph(mysurv2 ~ Age + Sex + Edu, data=dummy2.df_short, x=T, y=T)
  suppressWarnings(myboot <- validate(mymod, B=1000))
  Cind <- myboot["Dxy","index.corrected"] / 2 + 0.5
  print( c("No Genetics:", Cind))


}
