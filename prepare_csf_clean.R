

##sorting out the biomarker mess
upennbio <- read.csv("./UPENNBIOMK_MASTER_FINAL_21Nov2023.csv")
roche_elecsys <- subset(upennbio, RUNDATE>=as.Date("2016-11-17"))

#go through each RID and compute median of values at duplicate time points
roche_elecsys$ABETA42  <- as.double(gsub("<","", gsub(">","", roche_elecsys$ABETA42)))
roche_elecsys$TAU  <- as.double(gsub("<","", gsub(">","", roche_elecsys$TAU)))
roche_elecsys$PTAU <- as.double(gsub("<","", gsub(">","", roche_elecsys$PTAU)))

tag <- apply(roche_elecsys,1, function(x){ paste(as.integer(x["RID"]), x["VISCODE2"], sep="_") })

roche_elecsys_clean <- data.frame(t(sapply(unique(tag), function(tt){
  tmp <- roche_elecsys[tag == tt,]
  out_a <- tmp[1, c("PHASE","PTID","RID","VISCODE2","EXAMDATE")]
  out_b <- apply(tmp[,c("ABETA40","ABETA42","TAU","PTAU")],2, median, na.rm=T)
  return(unlist(c(out_a, out_b)))
})))

for ( bmk in c("RID","ABETA40","ABETA42","TAU","PTAU")){
  roche_elecsys_clean[[bmk]] <- as.double(roche_elecsys_clean[[bmk]])
}

write.csv(roche_elecsys_clean, "ROCHE_ELECSYS_biomk.csv", row.names=F)
