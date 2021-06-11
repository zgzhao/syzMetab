## game_mgraph <- function(nv, ne, nr, ng, seed=NULL){
##     if(!is.null(seed)) set.seed(seed)
##     g <- make_empty_graph(nv)
##     v.names <- .ordnames(1:nv, "C")
##     V(g)$name <- v.names
##     eii <- 0
##     for(va in v.names) {
##         for(vb in setdiff(v.names, va)) {
##             if(eii > ne) break
##             g <- add.edges(g, c(va, vb))
##         }
##     }
##     r.names <- .ordnames(1:nr, "R")
##     if(ne > nr) r.names <- rep(sample(r.names), length.out=ne)
##     E(g)$reaction <- r.names
## }
