setMethod("show", "ReactionSet", function(object) {
    cat("Reactions: S4 object holding KEGG reaction data\n")
    cat("\tOrganism: ", object@organism,"\n")
    cat("\tReactions: ", length(Reactions(object)), "\n")
    cat("\tCompounds: ", length(Compounds(object)), "\n")
})

setMethod("show", "keggPATH", function(object) {
    cat("keggPATH object: S4 class holding KEGG pathway meta data\n")
    cat("\tName: ", object@pathInfo$name, "\n")
    cat("\tTitle: ", object@pathInfo$title, "\n")
    cat("\tOrganism: ", object@pathInfo$org,"\n")
    cat("\tReactions: ", length(object@reactions), "\n")
    cat("\tCompounds: ", length(Compounds(object)), "\n")
    cat("\tGenes/Orthologs: ", length(Genes(object)), "\n")
})

setMethod("show", "stcuts", function(object) {
    cat("S4 object: s-t cuts of metabolic graph\n")
    cat("Edge cuts (@edges):", length(object@edges), "\n")
    cat(str(object@edges), "\n")
    cat("Gene cuts (@genes):", length(object@genes), "\n")
    cat("\t", str(object@genes), "\n")
})
