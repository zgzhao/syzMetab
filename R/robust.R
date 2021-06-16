#' Calculate mutational robustness of metabolic network
#'
#' NOTE: you must set substrates and products before doing this!
#' @title calculate mutational robustness
#' @param object mgraph with substrates and products set
#' @param GMR numeric (default 0.01), gene mutation rate/probability between 0 and 1.
#' @param prob TRUE/FALSE (default). Return functional probability of the network if TRUE.
#' @param FUN function translates failure probability to robustness value. Default: FUN=function(x) {-log(x)}
#' @return a number for robustness
#' @author ZG Zhao
#' @export
setGeneric("robust", function(object, GMR=0.01, prob=FALSE, FUN=function(x){-log(x)}) standardGeneric("robust"))
setMethod("robust", "NULL", function(object, GMR, prob, FUN) return(-100))
setMethod("robust", "bgraph", function(object, GMR, prob, FUN){
    gsets <- attr(object, "stcuts")@genes
    vfac <- names(gsets)
    ## grouping
    gfac <- list()
    for(ff in vfac){
        if(ff %in% unlist(gfac)) next
        ffs <- vs_accessed_by(object,ff)
        ffs <- intersect(vfac, ffs)
        gfac <- c(gfac, list(ffs))
    }
    ## Joint probabilities
    pps <- prob_stcuts(gsets, gfac, GMR)

    ## NOTE:
    ## Mostly: prod(1- pps) = 1 - sum(pps)
    ## However, if pps is too small, prod(1-pps) results in 0!
    FP <- sum(pps)
    if(prob) 1 - FP
    else FUN(FP)
})

setMethod("robust", "stcuts", function(object, GMR, prob, FUN){
    if(is.empty(object)) return(-100)
    g <- make_bgraph(object)
    robust(g, GMR=GMR, prob=prob, FUN=FUN)
})

