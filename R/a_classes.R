#' @name class_virtual
#' @title virtual classes
#' @aliases igraph mgraph KEntityList
#' @description Virtual classes designed in mgraph package.
#' @details There following classes are set as virutal:
#' - igraph
#' - mgraph
#' - KEntityList: actually list
#' @author zhao
NULL

#' @name class_ReactionSet
#' @title S4 class: ReactionSet
#' @description ReactionSet: S4 class holding reaction infos.
#' @details Three slots: reaction, organism and alias.
#' - reaction: list of reaction info (substrate, product and genes)
#' - organism: KEGG organism identifier
#' - alias: named list of chemical aliases. Set when merging chemicals.
#' @author ZG Zhao
NULL

#' @name class_keggPATH
#' @title S4 class: keggPATH
#' @aliases keggPATH
#' @seealso \code{\link{KEntityList}}, \code{\link{ReactionSet}}
#' @description "keggPATH" is S4 class designed for holding KEGG pathway information: entries, reactions, relations, graphics and other general information (pathInfo).
#' @details Typically, a keggPATH object can be obtained by \code{\link{make_mpath}} function, which will download and parse the KEGG xml file.
#' @author zhao
NULL

## virtual classes
setClassUnion("EntryList", "list")
setClassUnion("ReactionList", "list")
setClassUnion("RelationList", "list")
setClassUnion("GraphicList", "list")
setClassUnion("KEntityList", "list")
setClassUnion("igraph", "list")
setClassUnion("mgraph", "list")
setClassUnion("ggraph", "list")
setClassUnion("rgraph", "list")
setClassUnion("xgraph", "list")

setClass("ReactionSet",
         slots=c(reaction="ReactionList",
                 organism="character",
                 alias="list"))
setMethod("initialize", "ReactionSet", function(.Object, reactions, organism=""){
    .Object@reaction <- reactions
    .Object@organism <- organism
    .Object@alias <- list()
    return(.Object)
})

setClass("keggPATH",
         slots=c(pathInfo="list",
                 entries="EntryList",
                 reactions="ReactionList",
                 relations="RelationList",
                 graphics="GraphicList"))

setMethod("initialize", "keggPATH", function(.Object, ko, d.path) {
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
    return(.Object)
})

