#' @export
print.KEGGdata <- function(x) {
    str(x)
}

#' @export
print.mgraph <- function(g) {
    cat("KEGG metabolic network represented by igraph object:\n")
    cat("\tCompounds/Nodes: ", vcount(g), "\n")
    cat("\tReactions: ", length(getReactions(g)), "\n")
    cat("\tEdges: ", ecount(g), "\n")
}

names.ReactionSet <- function(x){
    names(getReactions(x))
}
length.ReactionSet <- function(x){
    length(getReactions(x))
}
