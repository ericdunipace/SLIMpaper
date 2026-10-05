rm(list=ls())
library(plyr)
library(survival)
#code by Steffen
####################################################################################################
############### load DFCI GBM data-base ############################################################
{
  data.1       = "../Data/GBM/GBM_DataBase_07_09_18.csv"
  data.2       = "../Data/GBM/CNSMolecularDatabase_deidentified.csv"

  na.strings   = c( "Technical Failure", "NA",  "Unknown, N/A", "Equivical Methylation", "900 - not assessed", "", "Unknown - N/A")
  RWE.DB.1     = read.csv(data.1, header = T, na.strings=na.strings); dim(RWE.DB.1)
  RWE.DB.2     = read.csv(data.2, header = T, na.strings=na.strings); dim(RWE.DB.2)
  n.1          = nrow(RWE.DB.1)
  n.2          = nrow(RWE.DB.2)
}

####################################################################################################
###### treatment variable & standard-of-care #######################################################
{
  TMZ.1          = RWE.DB.1$TMZ.     == "Yes"
  RT.1           = RWE.DB.1$RT.      == "Yes"
  noBV.1         = RWE.DB.1$Avastin. == "No"
  NoDrug.1       = is.na(RWE.DB.1$Trial.Drug)
  is.SOC.1       = TMZ.1 & RT.1 & noBV.1 & NoDrug.1
  id.ndGBM.1     = RWE.DB.1$Timng.of.presentation.for.radiation == "New"
  id.ndGBM.SOC.1 = which(is.SOC.1 & id.ndGBM.1)

  TMZ.2          = RWE.DB.2$TMZ.     == "Yes"
  RT.2           = RWE.DB.2$RT.      == "Yes"
  noBV.2         = RWE.DB.2$Avastin. == "No"
  NoDrug.2       = is.na(RWE.DB.2$Trial.Drug)
  is.SOC.2       = TMZ.2 & RT.2 & noBV.2 & NoDrug.2
  id.ndGBM.2     = RWE.DB.2$Timng.of.presentation.for.radiation == "New"
  id.ndGBM.SOC.2 = which(is.SOC.2 & id.ndGBM.2)
}

####################################################################################################
###### Pre-treatment variables #####################################################################
{
  ################## Resection ##################
  RWE.DB.1$Resection = RWE.DB.1$Extent.Of.Resection...Surgical.Note
  RWE.DB.2$Resection = RWE.DB.2$Extent.Of.Resection...Surgical.Note

  ################## MGMT #######################
  MGMT.1     = rep(NA, n.1)
  id         = RWE.DB.1$MGMT.Result == "Unmethylated";
  MGMT.1[id] = 0
  id         = RWE.DB.1$MGMT.Result == "Methylated";
  MGMT.1[id] = 1
  id         = RWE.DB.1$MGMT.Result == "Partial Methylation"
  MGMT.1[id] = 1

  MGMT.2     = rep(NA, n.2)
  id         = RWE.DB.2$MGMT.Result == "Unmethylated";
  MGMT.2[id] = 0
  id         = RWE.DB.2$MGMT.Result == "Methylated";
  MGMT.2[id] = 1
  id         = RWE.DB.2$MGMT.Result == "Partial Methylation"
  MGMT.2[id] = 1

  RWE.DB.1$MGMT = MGMT.1
  RWE.DB.2$MGMT = MGMT.2


  ################## KPS - score
  KPS.s           = levels(RWE.DB.1$KPS.at.presentation.for.radiation)[c(2:9,1)]
  KPS.S           = c(seq(20, 100, 10), NA)

  id              = sapply(RWE.DB.1$KPS.at.presentation.for.radiation, function(x){ J=x == KPS.s; ifelse(any(J), which(J), 9) })
  KPS             = KPS.S[id]
  KPS.C           = 1*(KPS>=90)
  RWE.DB.1$KPS    = KPS
  RWE.DB.1$KPS.C  = KPS.C

  id              = sapply(RWE.DB.2$KPS.at.presentation.for.radiation, function(x){ J=x == KPS.s; ifelse(any(J), which(J), 9) })
  KPS             = KPS.S[id]
  KPS.C           = 1*(KPS>=90)
  RWE.DB.2$KPS    = KPS
  RWE.DB.2$KPS.C  = KPS.C

  ################## IDH
  RWE.DB.1$IDH   = 1*(RWE.DB.1$IDH.status == "Positive")
  RWE.DB.2$IDH   = 1*(RWE.DB.2$IDH.status == "Positive")

  ################## RPA
  RWE.DB.1$RPA     = NA
  RWE.DB.2$RPA     = NA

}

####################################################################################################
###### Outcome distribution ########################################################################

{
  RWE.DB.1$T.OS.D = RWE.DB.1$os
  RWE.DB.1$C.OS.D = 1*(RWE.DB.1$Last.Follow.Up.Status=="Deceased")

  RWE.DB.2$T.OS.D = RWE.DB.2$OS
  RWE.DB.2$T.OS.T = RWE.DB.2$OS.from.TR
  RWE.DB.2$C.OS.D = 1*(RWE.DB.2$Last.Follow.Up.Status=="Deceased")
}

####################################################################################################
# save in Trial format
####################################################################################################

x.names  = c("Age", "Gender",   "KPS", "KPS.C", "MGMT", "RPA", "Resection", "IDH1") # "RPA"
x.list   = c("agedx", "Gender", "KPS", "KPS.C", "MGMT", "RPA", "Resection", "IDH")   # "RPA"
y.names  = c("T.OS.T", "C.OS.D") # "T.OS.D

DB.1     = RWE.DB.1[id.ndGBM.SOC.1,  c(x.list,     "RT.Start.Date", "T.OS.D", "C.OS.D")]
DB.2     = RWE.DB.2[id.ndGBM.SOC.2,  c(x.list[-1], "RT.Start.Date", "T.OS.D", "C.OS.D", "T.OS.T", "Diagnosis.Date")]

DB       = data.frame(DB.1, DB.2)
DB       = DB[-which(is.na(DB$RT.Start.Date)),]
DB       = DB[DB$T.OS.D >= DB$T.OS.T,]
DB       = DB[-is.na(DB$MGMT),]

RWE.X    = DB[,x.list]
RWE.OS   = DB[,y.names]

colnames(RWE.X)  = x.names
colnames(RWE.OS) = y.names

RWE.X            = transform(RWE.X,
                             Gender    = mapvalues(x=Gender,    from= c("Female", "Male"),     to=c(1,0)),
                             Resection = mapvalues(x=Resection, from= levels(RWE.X$Resection), to=c("Biopsy", "Sub Total", "Gross Total")))

path.1 = "../Data/GBM/DFCI_GBM_DB_Sub_X_OS.RData"
# path.2 = "../Data/GBM/DFCI_GBM_DB_Sub_X_OS.RData"
save(RWE.X, RWE.OS, file=path.1)
# save(RWE.X, RWE.OS, file=path.2)


KM = survfit(Surv(RWE.OS$T.OS.T, RWE.OS$C.OS.D)~1)

plot(KM)
abline(v=365)
abline(h=.7)











