
#' Make metabolic pathway (keggPATH object)
#'
#' This function reads from local folder (d.path) a KEGG xml file, which will be downloaded if not exists. Make sure internet is accessible.
#' Refer to \code{\link{keggPATH}} for class information.
#' TODO: make pathway from local KGML file regardless of d.path
#' @title make metabolic pathway
#' @param ko character, KEGG pathway ortholog (ko) index string such as ko00010, 00010 or ath00010
#' @param d.path file path for \code{\link{KEGG_get}}.
#' @return keggPATH object
#' @author ZG Zhao
#' @export
make_mpath <- function(ko, d.path="KEGG"){
    ko <- tolower(ko[1])
    new("keggPATH", ko, d.path)
}


#' Convert KOG names of a keggPATH object to real gene names
#'
#' In keggPATH objects, genes involved in reactions are generally represented by KEGG entry names (k genes or KOG).
#' @title General to organism-specific keggPATH conversion
#' @param kobj a keggPATH object
#' @param org abbreviation of an organism such as "ath"
#' @param d.path file path for \code{\link{KEGG_get}}
#' @return keggPATH object
#' @author ZG Zhao
#' @export
mpath_x_org <- function(kobj, org=NA, d.path="KEGG"){
    if(Organism(kobj) != "ko") stop("No a general KEGG pathway (eg. ko00010)!")
    if( is.empty(org) ) return(kobj)
    org <- tolower(org)
    gmap <- kogs_list(org, d.path)
    kogs <- names(gmap)
    ## filter reactions
    rns <- Reactions(kobj)
    rns <- lapply(rns, FUN=function(rr){
        gg <- intersect(rr$gene, kogs)
        gg <- if(length(gg) < 1) NA else unlist(gmap[gg])
        names(gg) <- NULL
        rr$gene <- sort(gg)
        rr
    })
    ss <- sapply(rns, FUN=function(x) ! is.empty(x$gene))
    rns <- rns[ss]
    class(rns) <- "ReactionList"
    kobj@reactions <- rns
    ## if(reaction.only) return(kobj)

    ## filter entries
    ents <- kobj@entries
    ents <- lapply(ents, FUN=function(rr){
        if(! rr$type %in% "ortholog") return(rr)
        oo <- intersect(rr$name, kogs)
        oo <- if(length(oo) < 1) NA else unlist(gmap[oo])
        names(oo) <- NULL
        ox <- grepl("^K[0-9]{5}$", oo)
        rr$name <- sort(oo[!ox])
        rr$type <- "gene"
        rr
    })
    ss <- sapply(ents, FUN=function(x) ! is.empty(x$name))
    ents <- ents[ss]
    class(ents) <- "EntryList"
    kobj@entries <- ents
    kobj@pathInfo <- lapply(kobj@pathInfo, FUN=function(x) {
        gsub("ko([0-9]+)", paste0(org, "\\1"), x)
    })
    kobj@pathInfo$org <- org
    kobj
}
