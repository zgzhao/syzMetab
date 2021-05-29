## NOTE:
## 1. DON'T TRY to put all essential infos of a KEGG reaction to nodes or edges of igraph object.
## 2. Reaction info are important for downstream calculations, they must be retained in graph object.
## 3. Save these info as attributes.

#' Reconstruction of KEGG metabolic network represented by an mgraph object.
#'
#' After conversion, compounds are represented by nodes and reactions by edges in the igraph object. Other information of KEGG metabolic networks are hold as attributes.
#' @title make graph from KEGGmeta object
#' @param kinfo \code{\link{KEGGmeta}} object
#' @return igraph object
#' @author zhao
#' @export
mgraph_from_kmeta <- function(kinfo) {
    if(! is.kmeta(kinfo)) stop("Not a KEGGmeta object!")
    rtns <- as_mreacts(kinfo)
    g <- mgraph_from_mreacts(rtns)
    g
}

mgraph_from_mreacts <- function(robj) {
    if(! is.mreacts(robj)) stop("Require MReactions object.")
    vv <- getCPDs(robj)
    g <- make_empty_graph(n=length(vv))
    V(g)$name <- vv
    rlist <- getReactions(robj)
    for(ndx in names(rlist)) {
        rx <- rlist[[ndx]]
        genes <- rx[["gene"]]
        ss <- rx[["substrate"]]
        pp <- rx[["product"]]
        ess <- t(expand.grid(ss, pp))
        g <- g %>% add.edges(ess, reaction=ndx)
    }
    attr(g, "reactions") <- robj
    class(g) <- c("mgraph", class(g))
    g
}

#' Reconstruct a KEGG pathway as mgraph network directly by providing the pathway identifiers (kos).
#'
#' KGML file will be downloaded online if you do not supply.
#' @title make graph from kos
#' @param kos character vector, KEGG pathway id(s). Mixing pathway ids of different organisms is not allowed.
#' @param d.path file path for \code{\link{KEGG_get}}.
#' @return mgraph object.
#' @author ZG Zhao
#' @export
mgraph_from_kos <- function(kos, d.path = "KEGG") {
    org <- unique(gsub("[0-9]+", "", kos))
    org <- setdiff(org, c("", "ko"))
    if(length(org) > 1) stop("Only one organism is allowed.")
    if(length(org) < 1) org <- "ko"
    kndx <- unique(gsub("[a-z]+", org, kos))
    rtns <- list()
    for(kx in kndx) {
        xinfo <- kmeta_from_ko(kx, d.path)
        rtns <- c(rtns, getReactions(xinfo))
    }
    robj <- as_mreacts(rtns, org)
    g <- mgraph_from_mreacts(robj)
    g
}

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
    robj <- as_mreacts(rtns, org)
    g <- mgraph_from_mreacts(robj)
    g
}

#' Append additional reactions to KEGG pathway.
#'
#' Reactions without KEGG index cannot be extrated from KGML files. This function help you add these reactions manually.
#' @title append reactions to mgraph manually
#' @param g graphMET object
#' @param df a data frame with at least four columns which supply "from", "to", "reversible", "genes" infos for the reactions. File path to a csv file containing such infos is also acceptable.
#' @return
#' @author ZG Zhao
#' @export
mgraph_append <- function(g, df) {
    if(! is.data.frame(df) && file.exists(df)) df <- read.csv(df)
    if(ncol(df) < 4) stop("Data should have at least four columns: from, to, reversible and genes.")
    if(all(c("from", "to", "reversible", "genes") %in% colnames(df)))
        df <- df[, c("from", "to", "reversible", "genes")]
    else
        warning("Column names do not match. The first four columns are used for from, to, reversible and genes!")

    robj <- getReactions(g, list.only=FALSE)
    for(i in 1:nrow(df)) {
        ss <- strsplit(df[i, 1], "[ ;,]+")[[1]]
        pp <- strsplit(df[i, 2], "[ ;,]+")[[1]]
        genes <- strsplit(df[i, 4], "[ ;,]+")[[1]]
        ne1 <- length(robj)
        robj <- mreacts_append(robj, ss, pp, genes)
        ne2 <- length(robj)
        ## add edge if only new reaction added
        if(ne2 > ne1) g <- g %>% add.edges(t(expand.grid(ss, pp)), reaction=paste0("RX", ne2))
        if(df[i, 3]) {
            robj <- mreacts_append(robj, pp, ss, genes)
            ne3 <- length(robj)
            if(ne3 > ne2) g <- g %>% add.edges(t(expand.grid(pp, ss)), reaction=paste0("RX", ne3))
        }
    }
    attr(g, "reactions") <- robj
    g
}

#' Get "clean" metabolic network with given substrates and products.
#'
#' Given primary substrates and end products, a "clean" metabolic network is a network in which all chemicals (vertices) and reactions (edges) are in the reaction chains from primary substrates to end products.
#' @title clean metabolic network
#' @param object mgraph object
#' @param s character vector, primary substrate names
#' @param p character vecter, end product names
#' @return mgraph object
#' @author ZG Zhao
#' @export
setGeneric("mgraph_clean", function(object, s, p) standardGeneric("mgraph_clean"))
setMethod("mgraph_clean", "mgraph", function(object, s, p){
    robj <- getReactions(object, list.only=FALSE)
    robj <- mreacts_merge_chems(robj, s, XCHEM1)
    robj <- mreacts_merge_chems(robj, p, XCHEM2)
    ## build graph for filter
    g <- mgraph_from_mreacts(robj)
    exx <- all_spaths_edges(g, XCHEM1, XCHEM2)
    ## build graph again
    g <- mgraph_from_mreacts(robj, e.names=exx)
    g
})


#' This function will: (1) merge any serial or parallel reactions; (2) update edge EP data accordingly.
#'
#' A simplified network contains no serial or parallel reactions.
#' @title Simplify the structure of a metabolic network
#' @param g metgraph object
#' @return metgraph object
#' @author ZG Zhao
#' @export
mgraph_simplify <- function(g) {
    ## NOTE: parallel robust
    if(! is.mgraph(g)) return(g)
    if(is.chemSet(g)) g <- mgraph_clean(g)
    g <- mgraph_update(g)

    while(TRUE) {
        dd <- graph::degree(g)
        ss <- dd$inDegree == 1 & dd$outDegree == 1
        if(sum(ss) < 1) break
        vx <- names(which(ss)[1])
        v1 <- unlist(inEdges(vx, g))
        v2 <- unlist(graph::adj(g, vx))
        p1 <- edgesEP(g, v1, vx)
        p2 <- edgesEP(g, vx, v2)
        px <- 1 - (1 - p1) * (1 - p2)
        g <- removeNode(vx, g)
        if(is.adjacentVs(g, v1, v2, "out")) px <- px * edgesEP(g, v1, v2)
        else g <- addEdge(v1, v2, g)
        names(px) <- NULL
        edgeData(g, v1, v2, "EP") <- px
    }
    ## update simple path info
    if(is.chemSet(g)) g@metaData$paths <- allReactionChains(g, VCHEM1, VCHEM2)
    g
}
