#' Calculate mutational robustness of metabolic network
#'
#' NOTE: you must set substrates and products before doing this!
#' @title calculate mutational robustness
#' @param object mgraph with substrates and products set
#' @param GMR numeric (default 0.01), gene mutation rate/probability between 0 and 1.
#' @param prob TRUE/FALSE (default). Return functional probability of the network if TRUE.
#' @return a number for robustness
#' @author ZG Zhao
#' @export
setGeneric("robust", function(object, GMR=0.01, prob=FALSE) standardGeneric("robust"))
setMethod("robust", "bgraph", function(object, GMR, prob){
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
    rbs <- prod(1 - pps)
    if(! prob) rbs <- -log(1 - rbs)
    rbs
})

setMethod("robust", "stcuts", function(object, GMR, prob){
    g <- make_bgraph(object)
    robust(g, GMR=GMR, prob=prob)
})

