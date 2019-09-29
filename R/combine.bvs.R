# combine.bvs <- function(sels = NULL, file.dir = NULL, file.names=NULL) {
#   if(!is.null(sels)){
#     all_probs <- sapply(sels, function(s) s$all_probs)
#     list_vis_covs <- lapply(sels, function(s) s$vis_mod)
#   } else if( ! is.null(file.names)){
#     n.run <- length(file.names)
#     runs.prob <- runs.mod <- vector("list", n.run)
#     for(idx in seq_along(file.names)) {
#       f <- file.names[idx]
#       sel <- readRDS(f)
#       runs.prob[[i]] <- sel$all_probs
#       runs.mod[[i]] <- sel$vis_mod
#     }
#     p <- ncol(sel$des_mat)
#     all_probs <- unlist(runs.prob)
#     list_vis_covs <- unlist(runs.mod, recursive=FALSE)
#   } else if (! is.null(file.dir)) {
#     file.names <- dir(file.dir)
#     n.run <- length(file.names)
#     runs.prob <- runs.mod <- vector("list", n.run)
#     for(idx in seq_along(file.names)) {
#       f <- file.names[idx]
#       sel <- readRDS(f)
#       runs.prob[[i]] <- sel$all_probs
#       runs.mod[[i]] <- sel$vis_mod
#     }
#     p <- ncol(sel$des_mat)
#     all_probs <- unlist(runs.prob)
#     list_vis_covs <- unlist(runs.mod, recursive=FALSE)
#   }
#   HashKey <- lapply(1:length(all_probs), function(i) rep(0, p))
#   for(i in seq_along(list_vis_covs)){
#     idx <- list_vis_covs[[i]]
#     HashKey[[i]][idx] <- 1
#     HashKey[[i]] <- c(HashKey[[i]])
#   }
#   HashKey <- unlist(HashKey)
#   uniques <- which(!duplicated(HashKey))
#   all_probs <- all_probs[uniques]
#   list_vis_covs <- list_vis_covs[uniques]
#   sout <- sort(all_probs, decreasing = T, index.return = T)
#   MaxMargs <- sout$x[1:outnum]
#   minds <- sout$ix[1:outnum]
#   max_marg <- MaxMargs[1]
#   indmax <- minds[1]
#   sel_model <- list_vis_covs[[indmax]] + 1
#   inc_probs <- BVSNLP:::inc_prob_calc(all_probs, list_vis_covs,
#                                       p)
#   median_model <- which(inc_probs >= 0.5)
#   return(list(inc_probs = inc_probs, MPM=median_model,
#              HPM = sel_model))
# }
#
# calc_key <- function(b){
#   a <- log(sort(b) + 1)
#   s <- length(a)
#   val <- 2^a + pi * a
#   out <- sum(val)
#   return(out)
# }
