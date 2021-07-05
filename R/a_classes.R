#' @name classes
#' @title package class designs
#' @description Virtual classes designed in mgraph package.
#' @details
#' Virtual classes:
#' - EntryList:: list
#' - ReactionList:: list
#' - RelationList:: list
#' - GraphicList:: list
#' - KLists:: EntryList, ReactionList, RelationList or GraphicList
#' - igraph: igraph
#' - mgraph: metabolite graph
#' - bgraph: bipartite graph of s-t cuts of genes
#' - rgraph: reaction graph
#' - ggraph: gene graph
#' - xgraph: virtual class of mgraph, rgraph, ggraph and bgraph
#'
#' S4 Classes:
#' - KDataSet: complete set of pathway data parsed from a KEGG xml file
#' - ReactionSet: reactions container used in xgraph object
#'
#' Short names or aliases:
#' - kdset: for `KDataSet`
#' - klist: for `KLists`
#' - rset: for `ReactionSet`
#' - rdata/rnames/rcount: reaction data/names/count
#' - edata/enames/ecount: edge data/names/count
#' - vdata/vnames/vcount: vertex data/names/count
#' @author zhao
#' @examples
#' library(syzMetab)
#' showClass("KDataSet")
#' showClass("ReactionSet")
#' showClass("mgraph")
#' showClass("xgraph")
NULL

## virtual classes
setClassUnion("EntryList", "list")
setClassUnion("ReactionList", "list")
setClassUnion("RelationList", "list")
setClassUnion("GraphicList", "list")
setClassUnion("KLists", c("EntryList", "ReactionList", "RelationList", "GraphicList"))
setClassUnion("igraph", "list")
setClassUnion("bgraph", "list")
setClassUnion("mgraph", "list")
setClassUnion("ggraph", "list")
setClassUnion("rgraph", "list")
setClassUnion("xgraph", c("mgraph", "ggraph", "rgraph", "bgraph"))
setClassUnion("sp.list", "list")

setClass("ReactionSet",
         slots=c(reaction="list",
                 organism="character"
                 ))
setMethod("initialize", "ReactionSet", function(.Object, reactions, organism=""){
    .Object@reaction <- reactions
    .Object@organism <- organism
    return(.Object)
})

setClass("KDataSet",
         slots=c(pathInfo="list",
                 entries="EntryList",
                 reactions="ReactionList",
                 relations="RelationList",
                 graphics="GraphicList",
                 chemicals="data.frame"
                 ))

setMethod("initialize", "KDataSet", function(.Object, ko, d.path) {
    if(missing(ko) || is.empty(ko)) {
        return(.Object)
    }
    ko <- tolower(ko)
    f.path <- KEGG_get(ko, d.path, f.type="kgml")
    kdata <- read_xml(f.path)
    kinfo <- xml_children(kdata)
    if(length(kinfo) > 0) {
        entries <- .parseEntryList(kinfo)
        .Object@entries   <- entries
        .Object@relations <- .parseRelationList(kinfo)
        .Object@reactions <- .parseReactionList(kinfo, entries)
        .Object@graphics  <- .parseGraphicsList(kinfo)
    }
    .Object@pathInfo <- as.list(xml_attrs(kdata))
    .Object@chemicals <- .parseCPD(ko, d.path)
    return(.Object)
})

setClass("stcuts", slots=c(edges="list", genes="list"))
