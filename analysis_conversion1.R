
### models the A-T- to A+(T-) conversion ###
# run analysis_clean.R before!
############################################

print_figures <- F
sensitivity <- F

##extract participants with A-T- profile
tmp <- unique(subset(longATN, ATN=="A-T-")$RID)
longAmTm <- subset(longATN, is.element(RID, tmp))

#build data for survival
hhh <- sapply(unique(longAmTm$RID), function(rid){
  tmp <- subset(longAmTm, RID==rid)
  tmp <- tmp[order(tmp$Years.bl),]

  #earliest entry with A-T-
  idx1 <- min(which(tmp$ATN == "A-T-"))
  #latest entry with A
  if (sensitivity){
    idx2 <- max(which(!is.na(tmp$A)))
  } else {
    idx2 <- max(which(!is.na(tmp$A) & !is.na(tmp$T)))
  }
  #get 'conversion'
  if (sensitivity){
    idx3 <- grep("A\\+", tmp$ATN)
  } else {
    idx3 <- min(which(tmp$ATN == "A+T-"))
  }
  idx3 <- min(idx3[idx3 >= idx1])

  tc <- NA
  if (length(idx3)>0){
    #tc <- tmp[idx3[1], "Years.bl"]
    tc <- tmp[idx3, "Years.bl"]
  }

  res <- c(rid, tmp[ idx1, "Years.bl"], tmp[idx2, "Years.bl"], tc )
})


###################################################################
#make a data frame
long.base <- subset(longAmTm, Years.bl==0)
dummy.df <- merge(long.base, t(hhh), by.x="RID", by.y="V1")

dummy.df <- merge(dummy.df, polygen2, by="RID")
dummy.df <- subset(dummy.df, EXCLUDE==0 & CEU>0.8)

#remove cases where end is before start
#latest observation time
valid <- dummy.df$V3
#for converters set time to observation of conversion
resp <- !is.na(dummy.df$V4)
valid[resp] <- dummy.df[resp, "V4"]
#time is valid if conversion happens after start
valid <- valid - dummy.df$V2

dummy.df <- subset(dummy.df, valid > 0 & !is.na(rs7412_T))

zero <- rep(0, nrow(dummy.df))
start <- dummy.df$V2
stop <- dummy.df$V3
resp <- !is.na(dummy.df$V4)
stop[resp] <- dummy.df[resp, "V4"]
stop <- stop - start
print(table(resp))

mysurv <- Surv(zero, stop, resp)
dummy.df_short <- dummy.df[,c("PTGENDER","rs7412_T","APOE4","PTEDUCAT","PC1","PC2","PC3","PC4","PC5",use_prs)]
dummy.df_short[["AGE"]] <- dummy.df$AGE + dummy.df$V2
colnames(dummy.df_short) <- c("Sex","APOEe2","APOEe4","Edu","PC1","PC2","PC3","PC4","PC5","PRS","Age")
#turn the selected PRS into a z-score
dummy.df_short$PRS <- (dummy.df_short$PRS - mean(dummy.df_short$PRS))/sd(dummy.df_short$PRS)

#dummy.df_short$APOEe4 <- as.factor(dummy.df_short$APOEe4)

if (adjust_for_PC){
  nn <- coxph(mysurv ~ Age + Sex + Edu + PC1 + PC2 + PC3 + PC4 + PC5 + APOEe2 + APOEe4 + PRS, data=dummy.df_short)
} else {
  nn <- coxph(mysurv ~ Age + Sex + Edu + APOEe2 + APOEe4 + PRS, data=dummy.df_short)
}

nn.ph <- cox.zph(nn)
print(summary(nn))

if (print_figures){
  #png("~/work/paper/inprogress/Stage_Varying_Risk/Figure1_Amyloid.png", width=400, height=700, res=300)
  pdf("~/work/paper/inprogress/Stage_Varying_Risk/Figure1_Amyloid.pdf", width=4.5)
  ggforest(nn, data=dummy.df_short)
  dev.off()

  pdf("~/work/paper/revision/Stage_Varying_Risk/R1/FigureS1_Amyloid_ph_20240328.pdf", width=9)
  ggcoxzph_fixed(nn.ph)
  dev.off()
}


##concordance:
# 81.12 % AV41 and CSF
# 78.07 % FBB and CSF

### computing C-index ###
if (0){
  library(rms)
  dummy.df_short$APOEe4 <- as.integer(dummy.df_short$APOEe4)

  mymod  <- cph(mysurv ~ Age + Sex + Edu + APOEe2 + APOEe4 + PRS, data=dummy.df_short, x=T, y=T)
  suppressWarnings(myboot <- validate(mymod, B=1000))
  Cind <- myboot["Dxy","index.corrected"] / 2 + 0.5
  print( c("Full:", Cind))

  mymod  <- cph(mysurv ~ Age + Sex + Edu + PRS, data=dummy.df_short, x=T, y=T)
  suppressWarnings(myboot <- validate(mymod, B=1000))
  Cind <- myboot["Dxy","index.corrected"] / 2 + 0.5
  print( c("No-APOE:", Cind))

  mymod  <- cph(mysurv ~ Age + Sex + Edu + APOEe2 + APOEe4, data=dummy.df_short, x=T, y=T)
  suppressWarnings(myboot <- validate(mymod, B=1000))
  Cind <- myboot["Dxy","index.corrected"] / 2 + 0.5
  print( c("No-PRS:", Cind))

  mymod  <- cph(mysurv ~ Age + Sex + Edu, data=dummy.df_short, x=T, y=T)
  suppressWarnings(myboot <- validate(mymod, B=1000))
  Cind <- myboot["Dxy","index.corrected"] / 2 + 0.5
  print( c("No Genetics:", Cind))

}
