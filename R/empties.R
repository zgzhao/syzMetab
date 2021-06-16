#' Uniform function for empty object test
#'
#' "Empty" refers to zero items or empty strings if character vector is given.
#' @title test empty
#' @param object vector, matrix, data.frame, list, ReactionSet, graph ...
#' @param ... multiple objects
#' @return TRUE/FALSE
#' @author ZG Zhao
#' @export
setGeneric("is.empty", function(object, ...) standardGeneric("is.empty"))
setMethod("is.empty", "NULL", function(object) TRUE)
setMethod("is.empty", "vector", function(object, ...) {
    xx <- unlist(c(object, ...))
    if (length(xx) < 1) return(TRUE)
    all(grepl("^\\s*$", xx))
})
setMethod("is.empty", "matrix", function(object, ...)
    is.empty(as.vector(c(object, ...))))

setMethod("is.empty", "list", function(object, ...)
    is.empty(unlist(c(object, ...))))
setMethod("is.empty", "sp.list", function(object) {
    class(object) <- "list"
    is.empty(object)
})

setMethod("is.empty", "data.frame", function(object) nrow(object) < 1)
setMethod("is.empty", "stcuts", function(object) is.empty(object@genes))
setMethod("is.empty", "KLists", function(object) {
    class(object) <- "list"
    is.empty(object)
})
setMethod("is.empty", "keggPATH", function(object) is.empty(object@reactions))
setMethod("is.empty", "xgraph", function(object) {
    vcount(object) < 1 || ecount(object) < 1
})
setMethod("is.empty", "igraph", function(object) vcount(object) < 1)

## return a graph without nodes and edges
#' @export
emptyGraph <- function(g){
    delete.vertices(g, vnames(g))
}
