if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install()
}
# source("https://bioconductor.org/biocLite.R")
if(!requireNamespace("curatedOvarianData", quietly = TRUE) ) {
  BiocManager::install(c("curatedOvarianData"))
}
