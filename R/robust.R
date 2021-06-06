#' Calculate mutational robustness of metabolic network
#'
#' NOTE: you must set substrates and products before doing this!
#' @title calculate mutational robustness
#' @param object mgraph with substrates and products set
#' @param s node names: primary substractes
#' @param p node names: end products
#' @param GMR numeric (default 0.01), gene mutation rate/probability between 0 and 1.
#' @param prob TRUE/FALSE (default). Return functional probability of the network if TRUE.
#' @return a number for robustness
#' @author ZG Zhao
#' @export
setGeneric("robust", function(object, s, p, GMR=0.01, prob=FALSE) standardGeneric("robust"))
setMethod("robust", "mgraph", function(object, s, p, GMR, prob){
    cat("Not implemented yet!\n")
    return(invisualble(NULL))

    ## empty network
    if(is.empty(object)) return(-1)
    if(GMR > 1 || GMR < 0) stop("Mutational probability of genes should be < 1 and > 0!")
    cutlists <- all_stcuts(object, s, p)
    if(! is.stcuts(cutlists)) return(0)

    ##======================================
    ## TODO: using PGM algorithms!
    gsets <- .xgsets(cutlists$gene.cuts)
    FPS <- sapply(gsets, FUN=function(x) {
        if(any("auto" %in% x)) return(1)
        GMR ^ length(x)
    })
    FPS <- sum(FPS)
    if(prob) rbs <- 1 - FPS
    else rbs <- -log(FPS)
    ##======================================
    rbs
})

## turn gene.cuts to independent sets.
.xgsets <- function(gcuts) {
    xgenes <- unique(unlist(gcuts))
    ss <- sapply(gcuts, FUN=function(x) xgenes %in% x)
    rownames(ss) <- xgenes
    sn <- apply(ss, 1, FUN = function(x) paste(as.integer(x), collapse=""))
    rex <- list()
    for(sx in unique(sn)) {
        ndx <- which(sn %in% sx)
        rex <- c(rex, list(xgenes[ndx]))
    }
    rex
}
