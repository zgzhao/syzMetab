#' @export
print.KLists <- function(x) {
    str(x)
}

#' @export
print.mgraph <- function(g) {
    cat("KEGG metabolic network represented by igraph object:\n")
    cat("\tCompounds/Nodes: ", vcount(g), "\n")
    cat("\tReactions: ", length(Reactions(g)), "\n")
    cat("\tEdges: ", ecount(g), "\n")
}

names.ReactionSet <- function(x){
    names(Reactions(x))
}
length.ReactionSet <- function(x){
    length(Reactions(x))
}
