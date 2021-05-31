setMethod("show", "ReactionSet", function(object) {
    cat("Reactions: S4 object holding KEGG reaction data\n")
    cat("\tOrganism: ", object@organism,"\n")
    cat("\tReactions: ", length(getReactions(object)), "\n")
    cat("\tCompounds: ", length(getCPDs(object)), "\n")
})

setMethod("show", "keggPATH", function(object) {
    cat("keggPATH object: S4 class holding KEGG pathway meta data\n")
    cat("\tName: ", object@pathInfo$name, "\n")
    cat("\tTitle: ", object@pathInfo$title, "\n")
    cat("\tOrganism: ", object@pathInfo$org,"\n")
    cat("\tReactions: ", length(object@reactions), "\n")
    cat("\tCompounds: ", length(getCPDs(object)), "\n")
    cat("\tGenes/Orthologs: ", length(getGenes(object)), "\n")
})
