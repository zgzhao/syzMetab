#' @name make_mgraph
#' @title make graph from KEGGmeta object
#' @aliases class_mgraph
#' @description Make mgraph object from various type of data: KEGGmeta, ReactionSet, ReactionList, list or KOs (KEGG pathway identifiers).
#' @details "mgraph" is a virtual S3 class defined in "gmetab" package.
#' - A mgraph object is actually a igraph object with a key attribute: "reactions", which is vital for downstream metwork analysis.
#' - Edges of a mgraph object have a "reactions" attribute holding the names of reactions involved, helpful for mapping edges back to reactions and genes.
#' - Users can set other graph, node or edge attributes
#' @usage make_mgraph(object, org=NULL)
#' @param object KEGGmeta, ReactionSet ReactionList, list or even igrph object.
#' @param org character, organism/species indentifier, parameter for "object" without organism info: ReactionList and list
#' @param d.path character, refer to \code{\link{KEGG_get}} for detail. Used when "object" is a KO vecter (KOs).
#' @author zhao
NULL

#' @export
setGeneric("make_mgraph", function(object, ...) standardGeneric("make_mgraph"))
setMethod("make_mgraph", "ReactionSet", function(object){
    ## basic method called by other methods!
    vv <- getCPDs(object)
    g <- make_empty_graph(n=length(vv))
    V(g)$name <- vv
    rlist <- getReactions(object)
    for(ndx in names(rlist)) {
        rx <- rlist[[ndx]]
        genes <- rx[["gene"]]
        ss <- rx[["substrate"]]
        pp <- rx[["product"]]
        ess <- t(expand.grid(ss, pp))
        g <- g %>% add.edges(ess, reaction=ndx)
    }
    attr(g, "reactions") <- object
    class(g) <- c("mgraph", class(g))
    g
})
setMethod("make_mgraph", "KEGGmeta", function(object){
    rsobj <- as_rset(object)
    make_mgraph(rsobj)
})
setMethod("make_mgraph", "character", function(object, d.path="KEGG"){
    org <- unique(gsub("[0-9]+", "", object))
    org <- setdiff(org, c("", "ko"))
    if(length(org) > 1) stop("Only one organism is allowed.")
    if(length(org) < 1) org <- "ko"
    kndx <- unique(gsub("[a-z]+", org, object))
    rtns <- list()
    for(kx in kndx) {
        xinfo <- kmeta_from_ko(kx, d.path)
        rtns <- c(rtns, getReactions(xinfo))
    }
    rsobj <- as_rset(rtns, org)
    make_mgraph(rsobj)
})
setMethod("make_mgraph", "ReactionList", function(object, org){
    rsobj <- as_rset(object, org)
    make_mgraph(rsobj)
})
setMethod("make_mgraph", "list", function(object, org){
    rvalid <- sapply(object, FUN=function(x){
        all(c("substrate", "product", "gene", "reversible") %in% names(x))
    })
    if(! all(rvalid)) stop("Invalid reaction list.")
    rsobj <- as_rset(rtns, org)
    make_mgraph(rsobj)
})

#' Normalize KEGG generic pathway to species specific pathway.
#'
#' Adjust and filter reactions, and rebuild graph object from reactions list.
#' @title general to organism-specific mgraph conversion
#' @param g mgraph object
#' @param org character
#' @param d.path file path for \code{\link{KEGG_get}}.
#' @return mgraph with organism set
#' @author ZG Zhao
#' @export
mgraph_x_org <- function(g, org, d.path = "KEGG"){
    if(! is.mgraph(g)) stop("Not a metabolic graph.")
    if(Species(g) != "ko") stop("Not a generic metabolic graph!")

    rtns <- getReactions(g)
    org <- org[1]
    gmap <- kogs_list(org, d.path)
    kogs <- names(gmap)
    rtns <- lapply(rtns, FUN=function(rr){
        gg <- intersect(rr$gene, kogs)
        gg <- if(length(gg) < 1) NA else unlist(gmap[gg])
        names(gg) <- NULL
        rr$gene <- gg
        rr$reversible <- FALSE
        rr
    })
    ss <- sapply(rtns, FUN=function(x) ! is.empty(x$gene))
    rtns <- rtns[ss]
    robj <- as_rset(rtns, org)
    g <- mgraph_from_rset(robj)
    g
}
