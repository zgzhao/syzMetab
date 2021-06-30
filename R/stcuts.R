#' Get st_cuts info: edges and genes
#'
#' Get st-cuts of edges (reactions). The edge cut-sets are translated into gene cut-sets also. Results can apply to \code{\link{make_bgraph}}.
#' @title list all gene cut sets
#' @param object xgraph object
#' @param s source node
#' @param t target node
#' @param minimal TRUE/FALSE. If TRUE (default), get st_min_cuts instead of all st_cuts.
#' @return list
#' - edges: edge cut sets
#' - genes: gent cut sets
#' @author ZG Zhao
#' @export
#' @examples
#' ## NOT RUN
#' library(syzMetab)
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
setGeneric("make_stcuts", function(object, s, t, minimal=TRUE) standardGeneric("make_stcuts"))
setMethod("make_stcuts", "NULL", function(object, s, t, minimal) return(NULL))
setMethod("make_stcuts", "mgraph", function(object, s, t, minimal){
    object <- subset_graph(object, s, t, clean=FALSE)
    if(is.empty(object)) return(NULL)
    s <- intersect(vnames(object), s)
    t <- intersect(vnames(object), t)
    if(length(s) > 1) {
        object <- merge_chems(object, s, "S")
        s <- "S"
    }
    if(length(t) > 1) {
        object <- merge_chems(object, t, "P")
        t <- "P"
    }

    ## finding st-cuts
    r.names <- rnames(object)
    e.names <- enames(object)
    rlist <- Reactions(object)
    if(minimal) {
        cps <- sapply(r.names, FUN=function(x){
            genes <- rlist[[x]][["gene"]]
            if(any(grepl("auto|HLINK", genes)))
                return(1000)
            length(genes)
        })
        xcuts <- st_min_cuts(object, s, t, cps)
        ecuts <- lapply(xcuts$cuts, as_ids)
    } else {
        xcuts <- st_cuts(object, s, t)
        ecuts <- lapply(xcuts$cuts, as_ids)
        ## vcuts <- lapply(xcuts$partition1s, as_ids)
    }

    ## map edge cuts to gene cuts
    gcuts <- sapply(ecuts, FUN=function(enn){
        ss <- r.names[e.names %in% enn]
        gns <- lapply(rlist[ss], FUN=function(rr) rr[["gene"]])
        gns <- sort(unique(unlist(gns)))
        ## cut set with autonomous reaction CANNOT be interrupted!
        if(any(grepl("^auto", gns))) return("")
        paste0(gns, collapse=" ")
    })
    ## length of ecuts may not equal to length of gcuts
    ecuts <- ecuts[gcuts != ""]
    gcuts <- gcuts[gcuts != "" & !duplicated(gcuts)]
    if(is.empty(gcuts)) gcuts <- list()
    else {
        gcuts <- lapply(gcuts, FUN=function(x) strsplit(x, " +")[[1]])
        names(gcuts) <- NULL
    }
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
    glist <- sapply(glist, FUN=function(x){
        genes <- intersect(x, kogs)
        if(is.empty(genes)) return("")
        genes <- unlist(gmap[genes])
        genes <- sort(.setGeneNames(genes))
        paste0(genes, collapse=" ")
    })
    glist <- glist[glist != ""]
    if(is.empty(glist)) glist <- list()
    else {
        glist <- glist[!duplicated(glist)]
        glist <- lapply(glist, FUN=function(x) strsplit(x, " +")[[1]])
        names(glist) <- NULL
        object@genes <- glist
    }
    object
}

