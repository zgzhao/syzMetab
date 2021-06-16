#' Get st_cuts info: edges and genes
#'
#' Get st-cuts of edges (reactions). The edge cut-sets are translated into gene cut-sets also. Results can apply to \code{\link{make_bgraph}}.
#' @title list all gene cut sets
#' @param object xgraph object
#' @param s source node
#' @param t target node
#' @return list
#' - min.genes
#' - min.edges
#' - edge.cuts
#' - gene.cuts
#' @author ZG Zhao
#' @export
#' @examples
#' ## NOT RUN
#' library(gmetab)
#' gm <- make_mgraph(c("ko00010", "ko00020"))
#' plot(gm)
#' chem1 <- c("C00031", "C00221", "C00267", "C01172", "C01451", "C06186")
#' chem2 <- "C00022"
#' gb <- make_bgraph(gm, s=chem1, t=chem2)
#' plot(gb)
#' xcuts <- make_stcuts(gm, s=chem1, t=chem2)
#' str(xcuts)
#' gx <- make_bgraph(xcuts)
#' plot(gx)
setGeneric("make_stcuts", function(object, s, t) standardGeneric("make_stcuts"))
setMethod("make_stcuts", "NULL", function(object, s, t) return(NULL))
setMethod("make_stcuts", "mgraph", function(object, s, t){
    object <- subset_graph(object, s, t, clean=FALSE)
    if(is.empty(object)) return(NULL)
    ## TIPS: st_cuts allows one s and one p only.
    ## Merge substrates and products, respectively.
    nns <- length(s)
    nnp <- length(t)
    if(nns > 0) {
        object <- add.edges(object, expand.grid("VCHEM1", s))
        s <- "VCHEM1"
    }
    if(nnp > 0) {
        object <- add.edges(object, expand.grid(t, "VCHEM2"))
        t <- "VCHEM2"
    }
    ## finding st-cuts
    xcuts <- st_cuts(object, s, t)
    ecuts <- lapply(xcuts[[1]], as_ids)
    vcuts <- lapply(xcuts[[2]], as_ids)

    ## map edge cuts to gene cuts
    r.names <- rnames(object)
    e.names <- enames(object)
    rlist <- Reactions(object)
    gcuts <- lapply(ecuts, FUN=function(enn){
        if(grepl("VCHEM", paste(enn, collapse=" "))) return(NULL)
        ss <- r.names[e.names %in% enn]
        gns <- lapply(rlist[ss], FUN=function(rr) rr[["gene"]])
        gns <- unique(unlist(gns))
        ## cut set with autonomous reaction CANNOT be interrupted!
        if(any(grepl("^auto", gns))) return(NULL)
        names(gns) <- NULL
        gns
    })
    ss <- sapply(gcuts, FUN=function(x) !is.empty(x))
    if(sum(unlist(ss)) < 1) return(NULL)
    ecuts <- ecuts[ss]
    gcuts <- gcuts[ss]
    names(gcuts) <- NULL
    new("stcuts", edges=ecuts, genes=gcuts)
})

#' Convert gene names and filter empty sets
#'
#' For time saving if you deal with a bunch of of organisms.
#' @title general (species="ko") stcuts to specices-specific stcuts
#' @param object stcuts object
#' @param org character
#' @param d.path file path for \code{\link{KEGG_get}}.
#' @return stcuts
#' @author ZG Zhao
#' @export
setorg_stcuts <- function(object, org, d.path="KEGG") {
    if(is.empty(object)) return(object)

    glist <- object@genes
    org <- tolower(org[1])
    gmap <- kogs_list(org, d.path)
    kogs <- names(gmap)
    glist <- lapply(glist, FUN=function(x){
        genes <- intersect(x, kogs)
        genes <- if(is.empty(genes)) NA else unlist(gmap[genes])
        .setGeneNames(genes)
    })
    ss <- sapply(glist, FUN=function(x) ! is.empty(x))
    object@edges <- object@edges[ss]
    object@genes <- glist
    object
}

