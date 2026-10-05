set.seed(09123085)
narray <- 100
seeds <- sample.int(.Machine$integer.max, narray)

saveRDS(seeds, "ovar_seeds.rds")
