setAs("KLists", "list", function(from) {
    class(from) <- NULL
    return(from)
})
setAs("xgraph", "igraph", function(from) {
    class(from) <- "igraph"
    return(from)
})
#' @title object test
#' @description functions for object test
#' @aliases is.klist is.rset is.mgraph is.rgraph is.ggraph is.xgraph is.chemset
#' @details Most function names are self-explanatory:
#' - is.kdset: KDataSet test
#' - is.klist: KList test
#' - is.rset: ReactionSet test
#' - is.mgraph(x)
#' - is.rgraph(x)
#' - is.ggraph(x)
#' - is.xgraph: mgraph, rgraph or mgraph
#' - is.chemset: test if primary substrates and end products were set.
#' - is.empty: test whether all objects are "empty". "Empty" objects are those of NULL, all NA, zero length and string containing spaces only. Graph objects without any edge are also regarded as "empty". S4 objects are not "empty" whatever data they have.
#' @param x object to be tested
#' @return TRUE/FALSE
#' @author ZG Zhao
#' @export
is.kdset <- function(x) {
    inherits(x, "KDataSet")
}

#' @export
is.klist <- function(x) {
    inherits(x, "KLists")
}

#' @export
is.rset <- function(x) {
    inherits(x, "ReactionSet")
}

#' @export
is.bgraph <- function(x) {
    inherits(x, "bgraph")
}
#' @export
is.mgraph <- function(x) {
    inherits(x, "mgraph")
}
#' @export
is.rgraph <- function(x) {
    inherits(x, "rgraph")
}
#' @export
is.ggraph <- function(x) {
    inherits(x, "ggraph")
}

#' @export
is.xgraph <- function(x){
    any(c("mgraph", "ggraph", "rgraph") %in% class(x))
}

is.chemset <- function(x){
    ss <- Substrates(x)
    pp <- Products(x)
    ! (is.empty(ss) || is.empty(pp))
}

is.stcuts <- function(x) {
    inherits(x, "stcuts")
}
