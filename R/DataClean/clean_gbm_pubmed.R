rm(list=ls())
library(plyr)

wd <- getwd()
#original code by Steffen

####################################################################################################
#####   Individual trials
####################################################################################################


{
  x.names   = c("Age", "Gender", "KPS",  "MGMT", "RPA",  "Resection", "IDH1")
  y.names   = c("T.OS", "C.OS")
  path      = "../Data/GBM/"
  setwd(path)

  # ####################################################################################################
  # #####   NCT01013285 -- UCLA Non-randomized trial
  # ####################################################################################################
  # Trial              = read.csv("RCT_2_NCT01013285/NCT01013285_RawData.csv");  head(Trial)  # read data
  # x.list             = c("age", "Gender", "KPS",  "MGMT", "RPA",  "Resection", "IDH1")
  # y.list             = c( "tts", "TTS.censor")
  # Trial.OS           = Trial[, y.list]
  # Trial.X            = Trial[, x.list]
  # colnames(Trial.X)  = x.names
  # colnames(Trial.OS) = y.names
  # Trial.X            = transform(Trial.X,
  #                                Gender    = mapvalues(x=Gender,    from= c("FEMALE", "MALE"),       to=c(1,0)),
  #                                RPA       = mapvalues(x=RPA,       from= c("III", "IV", "V", "VI"), to=3:6),
  #                                Resection = mapvalues(x=Resection, from= c("Biopsy", "Gross Total", "Sub Total", "NOT SPECIFIED"),
  #                                                      to  = c("Biopsy", "Gross Total", "Sub Total", NA)))
  # ########## NCT01013285 UCLA control arm (with NAs) - TMZ+RT
  # id                   = Trial$Control == 1
  # Trial.X.c            = Trial.X[id, ];       dim(Trial.X.c)
  # Trial.OS.c           = Trial.OS[id, ]
  # #                       save(Trial.X.c, Trial.OS.c,  file="Trial.c.NCT01013285.RData")
  #
  # ##########  NCT01013285 UCLA control arm 2 (with NAs) - TMZ+RT
  # id                   = Trial$Control == 2
  # Trial.X.c            = Trial.X[id, ];       dim(Trial.X.c)
  # Trial.OS.c           = Trial.OS[id, ]
  # #                       save(Trial.X.c, Trial.OS.c, file="Trial.c.NCT01013285.2.RData")
  #
  # ########## UCLA experimental arm (with NAs)  - TMZ+RT+BV
  # id                   = Trial$Control == 0
  # Trial.X.e            = Trial.X[id, ]
  # Trial.OS.e           = Trial.OS[id, ]
  # #                       save(Trial.X.e, Trial.OS.e,  file="Trial.e.NCT01013285.RData")
  #

  ####################################################################################################
  #####   PubMed-22120301
  ####################################################################################################
  Trial              = read.csv("PubMed22120301.csv");  head(Trial)  # read data
  x.list             = c("Age", "Sex", "Preop.KPS", "MGMT", "RPA", "OP.operation.")
  y.list             = c("OS", "OS.C")
  Trial.OS           = Trial[, y.list]
  Trial.OS[,1]       = Trial.OS[,1] * 30.5
  Trial.X            = Trial[, x.list]
  colnames(Trial.X)  = x.names[-7]
  colnames(Trial.OS) = y.names
  Trial.X            = transform(Trial.X,
                                 Gender    = mapvalues(x=Gender, from=c("Female", "Male"), to=c(1,0)),
                                 Resection = mapvalues(x=Resection, from= c("Total", "Subtotal"), to  = c("Gross Total", "Sub Total")) )


  ##########  control arm (without NAs)
  id                 = Trial$Arm == "Controll"
  Trial.X.c          = Trial.X[id, ]
  Trial.OS.c         = Trial.OS[id, ]
  save(Trial.X.c, Trial.OS.c, file="Trial.c.PubMed22120301.RData")

  ##########  experimental arm (without NAs)
  id                 = Trial$Arm == "Exp"
  Trial.X.e          = Trial.X[id, ]
  Trial.OS.e         = Trial.OS[id, ]
  save(Trial.X.e, Trial.OS.e, file="Trial.e.PubMed22120301.RData")

  # ####################################################################################################
  # #####   NCT00943826
  # ####################################################################################################
  # Trial1             = read.csv("RCT_1_NCT00943826/NCT00943826.csv");   head(Trial1)     # for outcome
  # Trial2             = read.csv("RCT_1_NCT00943826/NCT00943826_2.csv"); head(Trial2)     # for covariates
  # x.list             = c( "AGE",  "SEX", "KPS_BL", "MGMT",   "RNDRPA", "SURGTYP2")       # missing KPS, RPA, IDH1 (not included)
  # y.list             = c("TTDIED",  "CSDIED")
  # Trial.OS           = Trial1[, y.list]
  # Trial.X            = Trial2[, x.list]
  # colnames(Trial.X)  = x.names[-7]
  # colnames(Trial.OS) = y.names
  # Trial.X            = transform(Trial.X,
  #                                Gender    = mapvalues(x=Gender,    from = c("FEMALE", "MALE"),      to = c(1,0)),
  #                                KPS       = mapvalues(x=KPS,       from = c("50-80", "90-100", ""), to = c("50-80", "90-100", NA)),
  #                                RPA       = mapvalues(x=RPA,       from = levels(Trial.X$RPA),      to =  c(NA, 3:5, 3:5)),
  #                                MGMT      = mapvalues(x=MGMT,      from = c("METHYLATED", "NON-METHYLATED",  "MISSING"), to=c(1,0,NA)),
  #                                Resection = mapvalues(x=Resection, from = levels(Trial.X$Resection),to =c("Biopsy", "Gross Total", "Sub Total")) )
  #
  # ##########  control arm (with NAs) - RT & TMZ
  # id                 = which(Trial2[,"TRT1"] == "PLACEBO")
  # Trial.X.c          = Trial.X[id, ]
  # Trial.OS.c         = Trial.OS[id, ]
  # save(Trial.X.c, Trial.OS.c,  file="Trial.c.NCT00943826.RData")
  #
  # ##########  experimental arm (with NAs) - RT & TMZ
  # id                 = which(Trial2[,"TRT1"] == "BEVACIZUMAB")
  # Trial.X.e          = Trial.X[id, ]
  # Trial.OS.e         = Trial.OS[id, ]
  # save(Trial.X.e, Trial.OS.e,  file="Trial.e.NCT00943826.RData")
  #
  #
  # ####################################################################################################
  # #####   NCT00441142
  # ####################################################################################################
  # Trial1             = read.csv("RCT_5_NCT00441142/Subset_1.csv") # read data
  # Trial2             = read.csv("RCT_5_NCT00441142/Subset_2.csv") # read data
  # list.1             = c("Case.ID",  "Off.Tx.Code",  "Arm", "Age", "Gender", "KPS", "Type.of.Surg", "OS..days.", "Death.Status")
  # list.2             = c("Case.ID",  "Arm",  "MGMT.Status", "IDH.Status", "total.follow.up.time..days.", "Death.Status")
  # x.list             = c("Arm.x", "Age", "Gender", "KPS",  "MGMT", "Resection", "IDH1")
  # y.list             = c("OSd2", "OSc2")
  #
  # Trial              = merge(x=Trial1[, list.1], y=Trial2[, list.2], by="Case.ID", all.x = T, all.y = T)
  # colnames(Trial)    = c("ID", "Off.Tx", "Arm.x", "Age", "Gender", "KPS", "Resection", "OSd1", "OSc1", "Arm.y", "MGMT", "IDH1", "OSd2", "OSc2")
  # Trial              = Trial[ which(Trial$Off.Tx!="A"),]
  # Trial.X            = Trial[,x.list]
  # Trial.OS           = Trial[,y.list]
  # colnames(Trial.OS) = y.names
  #
  # Trial.X            = transform(Trial.X,
  #                                Gender    = mapvalues(x=Gender,    from= c("F", "M"),    to=c(1,0)),
  #                                IDH1      = mapvalues(x=IDH1,      from= c(0,1,-1,-2),   to=c(1,0, NA, NA)),
  #                                MGMT      = mapvalues(x=MGMT,      from= c(0,1,2,-1,-2), to=c(1,0,1,NA,NA)),
  #                                Resection = mapvalues(x=Resection, from= c(1,2,3),       to= c("Biopsy", "Sub Total", "Gross Total")) )
  #
  # ##########  control arm (with NAs) - RT & TMZ
  # id                 = which(Trial.X$Arm.x == "A")
  # Trial.X.c          = Trial.X[id, ]
  # Trial.OS.c         = Trial.OS[id, ]
  # save(Trial.X.c, Trial.OS.c,  file="Trial.c.NCT00441142.RData")
  #
  # ##########  experimental arm (with NAs)
  # id                 = which(Trial.X$Arm.x == "B")
  # Trial.X.e          = Trial.X[id, ]
  # Trial.OS.e         = Trial.OS[id, ]
  # save(Trial.X.e, Trial.OS.e,  file="Trial.e.NCT00441142.RData")
}
setwd(wd)
####################################################################################################
#####   Merge trials
####################################################################################################
{
  rm(list=ls())
  wd <- getwd()
  path      = "../Data/GBM/"
  setwd(path)
  cbind(dir())

  # ################################################ active control ####################################
  # ########## 1)  NCT01013285 - 1st UCLA control arm (with NAs) - TMZ+RT
  # load(file="Trial.c.NCT01013285.RData");
  # X = cbind(Trial.X.c,  Control="rwe.c1"); dim(X)
  # Y = cbind(Trial.OS.c, Control="rwe.c1"); dim(Y)
  #
  # ########## 2)  NCT01013285 - 2nd UCLA control arm (with NAs) - TMZ+RT
  # load(file="Trial.c.NCT01013285.2.RData");
  # X = rbind(X, cbind(Trial.X.c,  Control="rwe.c2")); dim(X)
  # Y = rbind(Y, cbind(Trial.OS.c, Control="rwe.c2")); dim(Y)
  #
  # ########## 3)  NCT00943826 - control arm (with NAs) - RT & TMZ
  # load(file="Trial.c.NCT00943826.RData");
  # X = rbind(X, cbind(Trial.X.c, "IDH1"=NA, Control="c1")); dim(X)
  # Y = rbind(Y, cbind(Trial.OS.c,           Control="c1")); dim(Y)

  ##########  4)  PubMed-22120301 control arm (without NAs) - RT & TMZ
  load(file="Trial.c.PubMed22120301.RData");
  # X = rbind(X, cbind(Trial.X.c, "IDH1"=NA, Control="c2")); dim(X)
  # Y = rbind(Y, cbind(Trial.OS.c,           Control="c2")); dim(Y)
  X = cbind(Trial.X.c, "IDH1"=NA, Control="c2"); dim(X)
  Y = cbind(Trial.OS.c,           Control="c2"); dim(Y)

  # ##########  5)  NCT00441142 control arm (with NAs) - RT & TMZ
  # load(file="Trial.c.NCT00441142.RData");
  # id  = c(2:5,8,6:7,9)
  # X = rbind(X, cbind(Trial.X.c,  "RPA"=NA, Control="c3")[,id]); dim(X)
  # Y = rbind(Y, cbind(Trial.OS.c,           Control="c3"));      dim(Y)


  # ################################################  experimental #####################################
  # ########## 6)  NCT01013285 -  UCLA experimental arm (with NAs)  - TMZ+RT+BV
  # load(file="Trial.e.NCT01013285.RData");
  # X = rbind(X, cbind(Trial.X.e,   Control="rwe.e1")); dim(X)
  # Y = rbind(Y, cbind(Trial.OS.e,  Control="rwe.e1")); dim(Y)
  #
  # ########## 7)  NCT00943826 - UCLA experimental arm (with NAs) - RT & TMZ
  # load(file="Trial.e.NCT00943826.RData");  names(Trial.X.e); dim(Trial.X.e) #
  # X = rbind(X, cbind(Trial.X.e, "IDH1"=NA, Control="e1")); dim(X)
  # Y = rbind(Y, cbind(Trial.OS.e,           Control="e1")); dim(Y)

  ########## 8)  PubMed-22120301 - china trial  experimental arm (without NAs)
  load(file="Trial.e.PubMed22120301.RData"); names(Trial.X.e); dim(Trial.X.e) # summary(Trial.X.e)
  X = rbind(X, cbind(Trial.X.e,  "IDH1"=NA,  Control="e2")); dim(X)
  Y = rbind(Y, cbind(Trial.OS.e,             Control="e2")); dim(Y)
  # X = cbind(Trial.X.e,  "IDH1"=NA,  Control="e2"); dim(X)
  # Y = cbind(Trial.OS.e,             Control="e2"); dim(Y)

  # ########## 9)  PubMed-NCT00441142 - DFCI trial experimental arm (without NAs)
  # load(file="Trial.e.NCT00441142.RData"); names(Trial.X.e); dim(Trial.X.e) # summary(Trial.X.e)
  # X = rbind(X, cbind(Trial.X.e,    "RPA"=NA, Control="e3")[,id]); dim(X)
  # Y = rbind(Y, cbind(Trial.OS.e,             Control="e3")); dim(Y)
}


####################################################################################################
#####  final changes
####################################################################################################
{
  X$KPS        = factor(x=X$KPS,  levels=c("60", "70", "80", "90", "100", "90-100", "50-80", NA), exclude = NA)
  X$MGMT       = factor(x=X$MGMT, levels=c(0,1, NA), exclude = NA)
  X$RPA        = factor(x=X$RPA,  levels=c(3:6, NA), exclude = NA)
  X$Resection  = factor(x=X$Resection,  levels=c(levels(X$Resection), NA), exclude = NA)
  X$IDH1       = factor(x=X$IDH1, levels=c(0,1, NA), exclude = NA)

  KPS.C                             = rep(NA, length(Y))
  KPS.C[which(X$KPS  == "90")]      = 1
  KPS.C[which(X$KPS  == "100")]     = 1
  KPS.C[which(X$KPS  == "90-100")]  = 1
  KPS.C[which(X$KPS  == "60")]      = 0
  KPS.C[which(X$KPS  == "70")]      = 0
  KPS.C[which(X$KPS  == "80")]      = 0
  KPS.C[which(X$KPS  == "50-80")]   = 0
  X$KPS.C                           = KPS.C


  path.1 = "Merged.ndGBM.RCTS.RData"
  # path.2 = "../Data/GBM/Merged.ndGBM.RCTS.RData"


  save(X, Y, file =path.1)
  # save(X, Y, file =path.2)
}
setwd(wd)

####################################################################################################
#####   Merge trials and DFCI DB
####################################################################################################
{
  rm(list=ls())
  path     = "../Data/GBM/"
  wd <- getwd()
  setwd(path)

  load("Merged.ndGBM.RCTS.RData");
  load("DFCI_GBM_DB_Sub_X_OS.RData");

  ## change col-names OS
  colnames(RWE.OS) = c("T.OS", "C.OS")
  RWE.OS.lab       = data.frame(RWE.OS,     Control="rwe.c0")
  RWE.X.lab        = data.frame(RWE.X[,-3], Control="rwe.c0")
  X.RCT.RWE        = rbind(RWE.X.lab,  X[,c(1:2,9,4:8)])
  Y.RCT.RWE        = rbind(RWE.OS.lab, Y)

  save(Y.RCT.RWE, X.RCT.RWE, file="Merged.ndGBM.RCTS.RWE.RData")
}
  setwd(wd)
####################################################################################################
####   Code for OS-12 Outcome
####################################################################################################
{
  rm(list=ls())
  wd <- getwd()
  path     = "../Data/GBM"
  setwd(path)

  load(file="Merged.ndGBM.RCTS.RWE.RData")

  t.e                      = 12 * 30.5    # time pt to evaluate OS (transform to binary)
  OS.12                    = rep(NA, nrow(Y.RCT.RWE))
  OS                       = Y.RCT.RWE$T.OS
  C.OS                     = Y.RCT.RWE$C.OS
  OS.12[OS>=t.e]           = 1
  OS.12[OS<t.e & C.OS==1 ] = 0
  Y.RCT.RWE$OS.12          = OS.12

  table(OS.12, Y.RCT.RWE$Control, useNA = "always")
  save(Y.RCT.RWE, X.RCT.RWE, file="Merged.ndGBM.RCTS.RWE.RData")


  ##################### no na
  id.vec    = c(paste0("rwe.c", 0:2), paste0("c", 1:3))
  x.var     = c("Age", "Gender", "KPS.C", "MGMT", "Resection", "Control" )

  id.sub    = sapply(X.RCT.RWE$Control, function(s) any(s==id.vec))
  X         = X.RCT.RWE[id.sub,x.var]
  Y         = Y.RCT.RWE[id.sub,]
  table(X$Control)

  id.sub    = (rowSums(is.na(X)) + rowSums(is.na(Y)) ) == 0
  X         = X[id.sub,]
  Y         = Y[id.sub,]
  table(X$Control)

  save(X, Y, file="Merged.ndGBM.RCTS.RWE.noNA.RData")


  ##################### select studies
  id.vec    = c(paste0("rwe.c", c(0,2)), paste0("c", 1:3))
  x.var     = c("Age", "Gender", "KPS.C", "MGMT", "Resection", "Control" )

  id.sub    = sapply(X.RCT.RWE$Control, function(s) any(s==id.vec))
  X         = X.RCT.RWE[id.sub,x.var]
  Y         = Y.RCT.RWE[id.sub,]
  table(X$Control)

  id.sub    = (rowSums(is.na(X)) + rowSums(is.na(Y)) ) == 0
  X         = X[id.sub,]
  Y         = Y[id.sub,]
  table(X$Control)



  #####################
  GLM        = glm(Y$OS.12 ~ Age +Gender +KPS.C+MGMT + KPS.C+Resection, data = X, family=binomial())
  x.new      = data.frame(Age=rep(median(X$Age),2), Gender=factor(c(0,1)), KPS.C=0, MGMT=factor(c(1,1)), Resection="Biopsy")

  LP         = predict(object = GLM, newdata=x.new)[1]
  D.HA       = .15
  TE.HA      = qlogis(D.HA+plogis(LP))-LP

  Y$P.OS.12.H0 = GLM$fitted.values
  Y$P.OS.12.HA = plogis(GLM$linear.predictors+TE.HA)

  mean(Y$P.OS.12.H0)
  mean(Y$P.OS.12.HA)

  save(X, Y, file="Merged.ndGBM.RCTS.RWE.Subset.RData")


  ##################### select studies for OS
  id.vec    = c(paste0("rwe.c", c(0,2)), paste0("c", 1:3))
  x.var     = c("Age", "Gender", "KPS.C", "MGMT", "Resection", "Control" )

  id.sub    = sapply(X.RCT.RWE$Control, function(s) any(s==id.vec))
  X         = X.RCT.RWE[id.sub,x.var]
  Y         = Y.RCT.RWE[id.sub,]
  table(X$Control)

  id.sub    = rowSums(is.na(X)) == 0
  X         = X[id.sub,]
  Y         = Y[id.sub,]
  table(X$Control)

  save(X, Y, file="Merged.ndGBM.RCTS.RWE.Subset.OS.RData")

}

setwd(wd)
