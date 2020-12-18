if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install()
}
# source("https://bioconductor.org/biocLite.R")
if(!requireNamespace("curatedOvarianData", quietly = TRUE) ) {
  BiocManager::install(c("curatedOvarianData"))
}

#### Data Cleaning/Examine ####
# source("data-raw/DataClean/downloadOvarian.R")
# data(package="curatedOvarianData")
require(curatedOvarianData)
require(SLIMpaper)

#### Load Data ####
data( TCGA_eset ) # has most observations
pheno <- phenoData(TCGA_eset)@data
express <- t(exprs(TCGA_eset))

#### remove cols pheno ####
no.var.pheno <- which(apply(pheno, 2, function(x) length(unique(x[!is.na(x)])) ) == 1)
if ( length(no.var.pheno) > 0) pheno <- pheno[,-no.var.pheno]

pheno$alt_sample_name <- NULL
pheno$flag_notes <- NULL
pheno$uncurated_author_metadata <- NULL


#### remove cols exprs ####
no.var.exprs <- which(apply(express, 2, function(x) length(unique(x[!is.na(x)])) ) == 1)
if ( length(no.var.exprs) > 0) express <- express[,-no.var.exprs]

#### Normalize Continuous Data ####
normalize_vars <- c("percent_normal_cells","percent_stromal_cells","percent_tumor_cells")
pheno[,normalize_vars] <- scale(pheno[,normalize_vars])
# express <- scale(express)

#### Set stage to factors ####

#### combine pheno/exprs ####
ovar <- cbind(pheno, express)

#### Check survival ####
sfit <- survival::survfit(formula = survival::Surv(days_to_death,
                                                   vital_status=="deceased") ~ -1,
                          data=ovar)
plot(sfit)

#### non-x vars ####
remove_all <- !(colnames(ovar) %in% c("unique_patient_ID", "days_to_tumor_recurrence","vital_status", "days_to_death",
                                      "os_binary", "relapse_binary","substage", "debulking")) &
  (colnames(ovar) %in% c("age_at_initial_pathologic_diagnosis", colnames(express)))
remove_rec <- !(colnames(ovar) %in% c(remove_all, "recurrence_status","site_of_tumor_first_recurrence")) & (colnames(ovar) %in% c("age_at_initial_pathologic_diagnosis","tumorstage", "summarygrade", colnames(express)))

#### fill in time for censored observations ####
event.recurr <- ovar$recurrence_status
time.recurr <- ovar$days_to_tumor_recurrence

event.death <- ovar$vital_status
time.death <- ovar$days_to_death

any(is.na(event.recurr))
any(is.na(event.death))
max.time <- max(c(time.recurr, time.death))
pmax.time <- pmax(time.recurr, time.death, na.rm=TRUE)
time.recurr[is.na(time.recurr)] <- pmax.time[is.na(time.recurr)]
time.recurr[is.na(time.recurr)] <- max.time

#### complete cases ####
complete.recurr <- complete.cases(cbind(time.recurr, event.recurr, ovar[, remove_rec]))
complete.death <- complete.cases(cbind(time.death, event.death, ovar[, remove_all]))


#### Truncate time and event variables ####
time.recurr <- time.recurr[complete.recurr]
event.recurr <- as.integer(event.recurr[complete.recurr] == "recurrence")

time.death <- time.death[complete.death]
event.death <- as.integer(event.death[complete.death] == "deceased")


#### Generate X variables ####
mm.recur <- model.matrix(formula("~."), data=ovar[complete.recurr, remove_rec])[,-1]
mm.death <- model.matrix(formula("~."), data=ovar[complete.death, remove_all])[,-1]

# mm.recur[,-1] <- scale(log(mm.recur[,-1]))
# mm.death[,-1] <- scale(log(mm.death[,-1]))

#### Reduce data size using cox regression ####
# set.seed(0193281)
set.seed(137371725) #from random.org
cl <- parallel::makeCluster(parallel::detectCores()-1)
doParallel::registerDoParallel(cl)
varsel.recur <- glmnet::cv.glmnet(mm.recur, cbind(time=time.recurr, status=event.recurr),
                                  family="cox",
                                  parallel=TRUE,
                                  nfolds = nrow(mm.recur))
varsel.death <- glmnet::cv.glmnet(mm.death, cbind(time=time.death, status=event.death),
                                  family="cox",
                                  parallel=TRUE,
                                  nfolds = nrow(mm.death))
parallel::stopCluster(cl)

# plot(varsel.recur)
mincvm <- which.min(varsel.recur$cvm)
covar.id <- which( as.numeric(coef(varsel.recur, varsel.recur$lambda[mincvm]))!=0 )
length(covar.id)

# sel <- max( which ( varsel.recur$cvm < (varsel.recur$cvm[mincvm] + varsel.recur$cvsd[mincvm]) ) )
# covar.id <- which( as.numeric(coef(varsel.recur, varsel.recur$lambda[sel]))!=0 )

# sel <- max( which ( varsel.recur$cvm < (varsel.recur$cvm[mincvm] + varsel.recur$cvup[mincvm]) ) )
# covar.id <- which( as.numeric(coef(varsel.recur, varsel.recur$lambda[sel]))!=0 )
# covar.id <- which( as.numeric(coef(varsel.recur, varsel.recur$lambda[mincvm]))!=0 )
# covar.small <- which( as.numeric(coef(varsel, varsel$lambda.min))!=0 )


mm.recur <- mm.recur[,covar.id]

mincvm <- which.min(varsel.death$cvm)
# sel <- max( which ( varsel.death$cvm < (varsel.death$cvm[mincvm] + varsel.death$cvup[mincvm]) ) )
# covar.id <- which( as.numeric(coef(varsel.death, varsel.death$lambda[sel]))!=0 )
covar.id <- which( as.numeric(coef(varsel.death, varsel.death$lambda[mincvm]))!=0 )
mm.death <- mm.death[,covar.id]

#currently not using death data so don't save it

#### create training, validation, and test data ####
shuffle.idx <- sample.int(nrow(mm.recur), nrow(mm.recur), replace = FALSE)
test.idx <- shuffle.idx[1]
validation.idx <- shuffle.idx[2:floor(nrow(mm.recur) * .1)]
train.idx <- shuffle.idx[!(shuffle.idx %in% c(test.idx, validation.idx))]

#### Save ovarian data as rds file ####
# output <- list(recurr = list(X     = mm.recur,
#                              event = event.recurr,
#                              time  = time.recurr),
#                death =  list(X     = mm.death,
#                              event = event.death,
#                              time  = time.death)
# )

pos <- which(apply(mm.recur, 2, function(x) all(x>0)))
mm.recur[,pos] <- log(mm.recur[,pos])

X <- mm.recur[train.idx,,drop=FALSE]
X <- scale(X)
center <- attributes(X)$`scaled:center`
scale <- attributes(X)$`scaled:scale`

ovar <- list(train = list(X = X,
                         event = event.recurr[train.idx],
                         time = time.recurr[train.idx]),
             val = list(X = scale(mm.recur[validation.idx,,drop = FALSE],
                                  center = center,
                                  scale = scale),
                        event = event.recurr[validation.idx],
                        time = time.recurr[validation.idx]),
             test = list(X = scale(mm.recur[test.idx,,drop = FALSE],
                                   center = center,
                                   scale = scale),
                         event = event.recurr[test.idx],
                         time = time.recurr[test.idx]))
usethis::use_data(ovar, compress="xz", version = 2)
# usethis::use_data(ovar, compress="xz", overwrite=TRUE, version = 2)
