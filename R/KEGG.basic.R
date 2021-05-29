#' API for downloding KEGG info with KEGG offical rest API.
#'
#' Download XML, image (png), AA sequence or NT sequence file from KEGG website. 
#' An alternative function is keggGet implemented in "KEGGREST" package.
#' @title KEGG download API
#' @param qid KEGG query id such as KEGG ortholog id like "ko00010", gene name like "ath:AT1G00010" or organism name like "ath" (for KEGG brite file only).
#' @param d.path Destination path (folder) to save or load downloaded file.
#' @param f.type Character string, file type of KEGG entry: "kgml", "aaseq", "ntseq" or "image". 
#' @param force.download TRUE/FALSE (default). Force update if TRUE and file exists.
#' @return path to the saved file.
#' @author ZG Zhao
#' @export
KEGG_get <- function(qid, d.path = "KEGG", f.type = c("kgml", "image", "htext"),
                     force.download = FALSE) {
    qid <- tolower(qid[1])
    if(grepl("^[a-z]+$", qid)) {
        ## for kegg brite file
        d.path <- .setKEGGLocalPath(d.path, "brite")
        f.path <- file.path(d.path, paste0(qid, "00001.txt"))
        xurl <- paste0("http://rest.kegg.jp/get/br:", qid, "00001")
    } else {
        ## others
        f.type <- intersect(c("kgml", "image", "htext"), tolower(f.type[1]))
        d.path <- .setKEGGLocalPath(d.path, f.type)
        if(length(f.type) < 1) return(invisible(NULL))
        
        ent <- tolower(qid)
        ext <- c(".xml", ".png", ".txt")
        names(ext) <- c("kgml", "image", "htext")
        ext <- ext[f.type]
        f.path <- file.path(d.path, paste0(ent, ext))
        
        kapi <- "http://rest.kegg.jp/get/"
        xtt <- if(f.type %in% c("kgml", "image")) paste0("/", f.type, "/") else ""
        xurl <- paste0(kapi, ent, xtt)
    }

    ## commit download
    if(! file.exists(f.path) || force.download) {
        if(! dir.exists(d.path)) dir.create(d.path, recursive = TRUE)
        try(download.file(xurl, f.path, method = "auto"), silent = TRUE, outFile = "kegg.log")
    }
    ## check after download
    if(! file.exists(f.path) || file.size(f.path) < 10) {
        warning("File does not exist or downloaded failed.")
        invisible(NULL)
    }
    return(f.path)
}


## parse KEGG XML tree infos, for internal use only.
.parseEntryList <- function(kdata) {
    itype <- xml_name(kdata)
    items <- kdata[itype == "entry"]
    if(length(items) < 1) return(NA)
    rids <- xml_attr(items, "id")
    rlist <- xml_attrs(items)
    rlist <- lapply(rlist, FUN=function(x) {
        rx <- x[setdiff(names(x), c("id", "link"))]
        rx <- gsub("\\b[a-z]+:", "", rx, ignore.case=TRUE)
        rx <- sapply(rx, FUN=function(z) strsplit(z, " +")[[1]])
        as.list(rx)
    })
    names(rlist) <- rids
    class(rlist) <- c("KEGGdata", class(rlist))
    rlist
}

.parseReactionList <- function(kdata, entries) {
    itype <- xml_name(kdata)
    items <- kdata[itype == "reaction"]
    if(length(items) < 1) return(NA)
    rlist <- lapply(items, FUN=function(x){
        xid <- xml_attr(x, "id")
        rtype <- xml_attr(x, "type") == "reversible"
        chem1 <- xml_attr(xml_find_all(x, ".//substrate"), "id")
        chem2 <- xml_attr(xml_find_all(x, ".//product"), "id")
        chem1 <- sapply(chem1, FUN=function(aa) entries[[aa]][["name"]])
        chem2 <- sapply(chem2, FUN=function(aa) entries[[aa]][["name"]])
        names(chem1) <- NULL
        names(chem2) <- NULL
        genes <- entries[[xid]][["name"]]
        list(substrate=chem1, product=chem2, reversible=rtype, gene=genes)
    })
    names(rlist) <- NULL
    class(rlist) <- c("KEGGdata", class(rlist))
    rlist
}

.parseGraphicsList <- function(kdata) {
    items <- xml_find_all(kdata, ".//graphics")
    ids <- xml_attr(xml_parent(items), "id")
    rlist <- xml_attrs(items)
    names(rlist) <- ids
    class(rlist) <- c("KEGGdata", class(rlist))
    rlist
}

.parseRelationList <- function(kdata) {
    itype <- xml_name(kdata)
    items <- kdata[itype == "relation"]
    if(length(items) < 1) return(NA)
    rlist <- lapply(items, FUN=function(x){
        rx1 <- xml_attrs(x)
        rx2 <- unlist(xml_attrs(xml_children(x)))
        c(rx1, rx2)
    })
    class(rlist) <- c("KEGGdata", class(rlist))
    rlist
}

## download and parse CPD info
.parseCPD <- function(KOPs, d.path="KEGG") {
    results <- NULL
    for (kpx in KOPs) {
        f.path <- KEGG_get(kpx, d.path=d.path, f.type="htext")
        dt <- readLines(f.path)
        ndx <- grep("^(COMPOUND|REFERENCE)", dt)
        ggs <- dt[ndx[1]:(ndx[2] - 1)]
        ggs <- gsub("^COMPOUND", "", ggs)
        ccs <- gsub("^ *([^ ]+).*$", "\\1", ggs)
        des <- gsub("^ *[^ ]+ +(.*)$", "\\1", ggs)
        dt <- data.frame(CPD = ccs, name = des)
        results <- rbind(results, dt)
    }
    results
}

## used by KEGGmeta initializer.
.transReacts <- function(reactions, entries) {
    rownames(entries) <- entries$id
    reactions$from <- sapply(reactions$substrate, FUN=function(xid){
        rex <- entries[xid, "name"]
        unique(unlist(rex))
    })
    reactions$to <- sapply(reactions$product, FUN=function(xid){
        rex <- entries[xid, "name"]
        unique(unlist(rex))
    })
    reactions$KOG <- sapply(reactions$id, FUN=function(xid){
        rex <- entries[xid, "name"]
        unique(unlist(rex))
    })
    cns <- colnames(reactions)
    colnames(reactions)[cns == "type"] <- "rev"
    reactions$rev <- reactions$rev == "reversible"
    reactions
}

.setKEGGLocalPath <- function(d.path, type){
    if(grepl("KEGG$", d.path, ignore.case = TRUE))
        return(file.path(d.path, type))
    ## left for user
    return(d.path)
}

.filterEList <- function(entries, by, vals) {
    by <- unlist(by)[1]
    ss <- sapply(entries, FUN=function(x) {
        if(by %in% names(x)) x[by] %in% vals
        else FALSE
    })
    entries[ss]
}
