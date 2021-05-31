setIs("EntryList", "KEntityList")
setIs("ReactionList", "KEntityList")
setIs("RelationList", "KEntityList")
setIs("GraphicList", "KEntityList")
setAs("KEntityList", "list", function(from) {
    class(from) <- NULL
    return(from)
})
setAs("mgraph", "igraph", function(from) {
    class(from) <- "igraph"
    return(from)
})

#' @title object test
#' @description helper functions for object test
#' @name helpers_is
#' @aliases is.kdata is.rset is.igraph is.mgraph is.xgraph is.chemset
#' @details Most function names are self-explanatory:
#' - is.mpath: keggPATH test
#' - is.kdata: KEGGtest test
#' - is.rset: ReactionSet test
#' - is.igraph(x)
#' - is.mgraph(x)
#' - is.xgraph: igraph or mgraph
#' - is.chemset: test if primary substrates and end products were set.
#' - is.empty: test whether all objects are "empty". "Empty" objects are those of NULL, all NA, zero length and string containing spaces only. Graph objects without any edge are also regarded as "empty". S4 objects are not "empty" whatever data they have.
#' @return TRUE/FALSE
#' @author ZG Zhao
NULL

#' @export
is.mpath <- function(x) {
    inherits(x, "keggPATH")
}

#' @export
is.kdata<- function(x) {
    inherits(x, "KEntityList")
}

#' @export
is.rset <- function(x) {
    inherits(x, "ReactionSet")
}

#' @export
is.mgraph <- function(x) {
    inherits(x, "mgraph")
}

#' @export
is.xgraph <- function(x){
    is.igraph(x) || is.mgraph(x)
}

is.chemset <- function(x){
    ss <- getSubstrate(x)
    pp <- getProduct(x)
    if(is.empty(c(ss, pp))) return(FALSE)
    return(TRUE)
}

#' @export
is.empty <- function(...){
    ss <- lapply(list(...), FUN=function(x){
        if(is.igraph(x)) return(vcount(x) < 1)
        if(isS4(x)) return(FALSE)
        x <- unlist(x, recursive=TRUE)
        length(x) < 1 || all(is.na(x)) || is.null(x) || all(grepl("^\\s*$", x))
    })
    all(unlist(ss))
}

## return a graph without nodes and edges
.emptyGraph <- function(g){
    deleteNodes(g, nodes(g))
}
