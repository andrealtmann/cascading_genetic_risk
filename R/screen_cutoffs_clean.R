
zero_to_A <- c()
A_to_T <- c()

for (tcut_idx in 1:10){
  message("Running analysis with cutoff ", tcut_idx)
  source("analysis_clean.R")

  source("analysis_conversion1.R")
  #get info from nn model:
  tmp <- c(Tcut, tau_cut, nn$n, nn$nevent, summary(nn)$coef[5:7,5])
  zero_to_A <- rbind(zero_to_A, tmp)

  source("analysis_conversion2.R")
  tmp2 <- c(Tcut, tau_cut, mm$n, mm$nevent, summary(mm)$coef[5:7,5])
  A_to_T <- rbind(A_to_T, tmp2)

  print(tmp)
  print(tmp2)

}
