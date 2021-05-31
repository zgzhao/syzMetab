## NOTE:
## 1. DON'T TRY to put all essential infos of a KEGG reaction to nodes or edges of igraph object.
## 2. Reaction info are important for downstream calculations, they must be retained in graph object.
## 3. Save these info as attributes.


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

    robj <- Reactions(g, list.only=FALSE)
    for(i in 1:nrow(df)) {
        ss <- strsplit(df[i, 1], "[ ;,]+")[[1]]
        pp <- strsplit(df[i, 2], "[ ;,]+")[[1]]
        genes <- strsplit(df[i, 4], "[ ;,]+")[[1]]
        ne1 <- length(robj)
        robj <- rset_append(robj, ss, pp, genes)
        ne2 <- length(robj)
        ## add edge if only new reaction added
        if(ne2 > ne1) g <- g %>% add.edges(t(expand.grid(ss, pp)), reaction=paste0("RX", ne2))
        if(df[i, 3]) {
            robj <- rset_append(robj, pp, ss, genes)
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
mgraph_clean <- function(g, s, p){
    spp <- all_spaths_list(g, s, p)
    vss <- all_spaths_nodes(spp)
    ess <- all_spaths_edges(spp)
    vxx <- setdiff(vnames(g), vss)
    g <- vsdelete(g, vxx)
    exx <- setdiff(enames(g), ess)
    g <- esdelete(g, exx)
    ## save time-consuming results
    attr(g, "spaths") <- spp
    g
}

