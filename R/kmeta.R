
#' Refer to \code{\link{KEGGmeta}} for class information.
#'
#' This function reads from local folder (d.path) a KEGG xml file, which will be downloaded if not exists. Make sure internet is accessible.
#' @title Generate KEGGmeta object
#' @param ko character, KEGG pathway ortholog (ko) index string such as ko00010, 00010 or ath00010
#' @param d.path file path for \code{\link{KEGG_get}}.
#' @return KEGGmeta object
#' @author ZG Zhao
#' @export
kmeta_from_ko <- function(ko, d.path="KEGG"){
    ko <- tolower(ko[1])
    new("KEGGmeta", ko, d.path)
}


#' Convert KOG names of a KEGGmeta object to real gene names
#'
#' In KEGGmeta objects, genes involved in reactions are generally represented by KEGG entry names (k genes or KOG).
#' @title General to organism-specific KEGGmeta conversion
#' @param kinfo a KEGGmeta object
#' @param org abbreviation of an organism such as "ath"
#' @param d.path file path for \code{\link{KEGG_get}}
#' @return KEGGmeta object
#' @author ZG Zhao
#' @export
kmeta_x_org <- function(kinfo, org=NA, d.path="KEGG"){
    if(kinfo@pathInfo$org != "ko") stop("No a general KEGG pathway (eg. ko00010)!")
    if( is.empty(org) ) return(kinfo)
    org <- tolower(org)
    gmap <- kogs_list(org, d.path)
    kogs <- names(gmap)
    ## filter reactions
    rns <- Reactions(kinfo)
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
    kinfo@reactions <- rns
    ## if(reaction.only) return(kinfo)

    ## filter entries
    ents <- kinfo@entries
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
    kinfo@entries <- ents
    kinfo@pathInfo <- lapply(kinfo@pathInfo, FUN=function(x) {
        gsub("ko([0-9]+)", paste0(org, "\\1"), x)
    })
    kinfo@pathInfo$org <- org
    kinfo
}
