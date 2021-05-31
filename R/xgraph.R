#' Get "clean" metabolic network with given substrates and products.
#'
#' Given primary substrates and end products, a "clean" metabolic network is a network in which all chemicals (vertices) and reactions (edges) are in the reaction chains from primary substrates to end products.
#' @title clean metabolic network
#' @param g xgraph object
#' @param s character vector, primary substrate names
#' @param p character vecter, end product names
#' @return mgraph object
#' @author ZG Zhao
#' @export
xgraph_clean <- function(g, s, p, mc.cores=detectCores() - 1){
    spp <- all_spaths_list(g, s, p, mc.cores)
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
