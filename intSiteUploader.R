## check for presence of R packages
rPackages <- c("stats", "RMySQL", "GenomicRanges", "BiocGenerics", "parallel", "IRanges", "GenomeInfoDb")
rPackagesPresent <- is.element(rPackages, installed.packages()[,1])
if(any(!rPackagesPresent)){
    stop(paste(rPackages[!rPackagesPresent]), " is not available")
}
stopifnot(sapply(rPackages, require, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE))

options(stringsAsFactors=F)

## check if file exist and permission .my.cnf 
stopifnot(file.exists("~/.my.cnf"))
stopifnot(file.info("~/.my.cnf")$mode == as.octmode("600"))


##check for presence of command line stuff
commandLinePrograms <- c("mysql")
programsPresent <- !sapply(sprintf("which %s > /dev/null 2>&1", commandLinePrograms), system)
if(any(!programsPresent)){
  stop(paste(commandLinePrograms[!programsPresent]), " is not available")
}

## working directory (i.e. primary analysis directory) is passed in via command line
args <- commandArgs(trailingOnly=TRUE)
workingDir <- args[1]
if( interactive() | is.na(workingDir) ) workingDir <- "."
stopifnot(!is.na(workingDir))
workingDir <- normalizePath(workingDir, mustWork=TRUE)
setwd(workingDir)
message("Changed to directory: ", workingDir)

## get miseqid
miseqid <- unique(readLines("miseqid.txt"))
stopifnot(length(miseqid)==1)
message("miseqid: ", miseqid)

## get sample information
stopifnot(all(file.exists("sampleInfo.tsv", "completeMetadata.RData")))
metadata <- get(load('completeMetadata.RData'))
metadata <- subset(metadata, select=c("alias", "gender", "refGenome"))
names(metadata) <- c("sampleName", "gender", "refGenome")

## initialize connection to database
## ~/.my.cnf must be present
junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn <- dbConnect(MySQL(), group="intSitesDev237") 
stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)

## stop if any sample is already loaded
allSampleName <- dbGetQuery(dbConn, "SELECT DISTINCT sampleName FROM samples")
is.loaded <- metadata$sampleName %in% allSampleName$sampleName
if(any(is.loaded)) message(
    paste0("Sets already in the database: ",
           paste(metadata$sampleName[is.loaded], collapse="\n")))

if( any(grepl("^GTSP", metadata$sampleName[is.loaded], ignore.case=TRUE)) ) stop("GTSP sample already loaded, delete from the database or leave them alone")

metadata <- subset(metadata, !is.loaded)

## Get max sampleID, and start from max+1
stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)
currentMaxSampleID <- as.integer(suppressWarnings(dbGetQuery(dbConn, "SELECT MAX(sampleID) AS sampleID FROM samples;")))
if(is.na(currentMaxSampleID)) {
    nrows <- as.integer(dbGetQuery(dbConn, "SELECT count(*) FROM samples;"))
    if(nrows==0) currentMaxSampleID<-0
    if(nrows!=0) stop("Failed to get currentMaxSampleID")
}

## load table samples
metadata$sampleID <- seq(nrow(metadata))+currentMaxSampleID
metadata$miseqid <- miseqid
stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)
stopifnot( dbWriteTable(dbConn, "samples", metadata, append=T, row.names=F) )
## check wether load was successful
sample.tab <- suppressWarnings(dbReadTable(dbConn, "samples"))
merged.tab <- merge(metadata, sample.tab, by="sampleName", all.x=TRUE)
if( !all(merged.tab$sampleID.x==merged.tab$sampleID.y) ) {
    message("Sample ID error, check the following table")
    print(merged.tab)
}

## Get max siteID, and start from max+1
stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)
currentMaxSiteID <- as.integer(suppressWarnings(dbGetQuery(dbConn, "SELECT MAX(siteID) AS siteID FROM sites;")))
if(is.na(currentMaxSiteID)) { 
    nrows <- as.integer(dbGetQuery(dbConn, "SELECT count(*) FROM sites;"))
    if(nrows==0) currentMaxSiteID<-0
    if(nrows!=0) stop("Failed to get currentMaxSiteID")
}

## Get max MultihitID, and start from max+1
stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)
currentMaxMultihitID <- as.integer(suppressWarnings(dbGetQuery(dbConn, "SELECT MAX(multihitID) AS multihitID FROM multihitpositions;")))
if(is.na(currentMaxMultihitID)) {
    nrows <- as.integer(dbGetQuery(dbConn, "SELECT count(*) FROM multihitpositions;"))
    if(nrows==0) currentMaxMultihitID<-0
    if(nrows!=0) stop("Failed to get currentMaxSiteID")
}

## process by sample
for(i in seq(nrow(metadata))){
    file <- metadata[i,"sampleName"]
    message("\nProcessing: ", file)
    if(all(file.exists(paste0(file, "/sites.final.RData"), paste0(file, "/allSites.RData")))){
        load(paste0(file, "/sites.final.RData"))
        load(paste0(file, "/allSites.RData"))
        if(length(sites.final)>0){
            ##sites.final won't exist if there aren't sites, thus no need to check if sites.final has sites in it
            sites <- data.frame(
                "siteID"=seq(length(sites.final))+currentMaxSiteID,
                "sampleID"=metadata[i,"sampleID"],
                "position"=start(flank(sites.final, -1, start=T)),
                "chr"=as.character(seqnames(sites.final)),
                "strand"=as.character(strand(sites.final)) )
            
            ## change to the right class as database
            sites$siteID <- as(sites$siteID, "integer")
            sites$sampleID <- as(sites$sampleID, "integer")
            sites$position <- as(sites$position, "integer")
            sites$chr <- as(sites$chr, "character")
            sites$strand <- as(sites$strand, "character")
            
            ##Newer versions of intSiteCaller return allSites in the order dictated by
            ##sites.final.  This line allows import of 'legacy' output
            allSites <- allSites[unlist(sites.final$revmap)]
            
            ##could do the next three statements with aggregate, but this method is emperically 2x faster
            pcrBreakpoints <- sort(paste0(as.integer(Rle(sites$siteID, sapply(sites.final$revmap, length))),
                                          ".",
                                          start(flank(allSites, -1, start=F))))
            
            condensedPCRBreakpoints <- strsplit(unique(pcrBreakpoints), "\\.")
            
            pcrBreakpoints <- data.frame("siteID"=sapply(condensedPCRBreakpoints, "[[", 1),
                                         "breakpoint"=sapply(condensedPCRBreakpoints, "[[", 2),
                                         "count"=runLength(Rle(match(pcrBreakpoints, unique(pcrBreakpoints)))))
            
            ## change to the right class as database
            pcrBreakpoints$siteID <- as(pcrBreakpoints$siteID, "integer")
            pcrBreakpoints$breakpoint <- as(pcrBreakpoints$breakpoint, "integer")
            pcrBreakpoints$count <- as(pcrBreakpoints$count, "integer")
            
            ## load table sites            
            stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)
            message("Loading sites: ", nrow(sites), " entries")
            stopifnot( dbWriteTable(dbConn, "sites", sites, append=T, row.names=F) )
            ## check loaded sites
            sql <- sprintf("SELECT * FROM sites WHERE siteID>=%s AND siteID<=%s",
                           range(sites$siteID)[1],
                           range(sites$siteID)[2])
            sites.from.db <- suppressWarnings( dbGetQuery(dbConn, sql) )
            sites.from.db$siteID <- as(sites.from.db$siteID, "integer")
            sites.from.db$sampleID <- as(sites.from.db$sampleID, "integer")
            sites.from.db$position <- as(sites.from.db$position, "integer")
            sites.from.db$chr <- as(sites.from.db$chr, "character")
            sites.from.db$strand <- as(sites.from.db$strand, "character")
            
            sites <- plyr::arrange(sites, siteID, sampleID, position, chr, strand)
            sites.from.db <- plyr::arrange(sites.from.db, siteID, sampleID, position, chr, strand)
            if(!identical(sites, sites.from.db)) {
                save.image("debug.rdata")
                stop("sites, sites.from.db not identical")
            }
            
            ## load table pcrbreakpoints
            stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)
            message("Loading pcrbreakpoints: ", nrow(pcrBreakpoints), " entries")
            stopifnot( dbWriteTable(dbConn, "pcrbreakpoints", pcrBreakpoints, append=T, row.names=F) )
            ## check loaded pcrbreakpoints
            sql <- sprintf("SELECT * FROM pcrbreakpoints WHERE siteID>=%s AND siteID<=%s",
                           range(sites$siteID)[1],
                           range(sites$siteID)[2])
            pcrbreakpoints.from.db <- suppressWarnings( dbGetQuery(dbConn, sql) )
            pcrbreakpoints.from.db$siteID <- as(pcrbreakpoints.from.db$siteID, "integer")
            pcrbreakpoints.from.db$breakpoint <- as(pcrbreakpoints.from.db$breakpoint, "integer")
            pcrbreakpoints.from.db$count <- as(pcrbreakpoints.from.db$count, "integer")
            
            pcrBreakpoints <- plyr::arrange(pcrBreakpoints, siteID, breakpoint, count)
            pcrbreakpoints.from.db <- plyr::arrange(pcrbreakpoints.from.db, siteID, breakpoint, count)
            if(!identical(pcrBreakpoints, pcrbreakpoints.from.db)) {
                save.image("debug.rdata")
                stop("pcrBreakpoints, pcrbreakpoints.from.db not identical")
            }
            
            newMaxSiteID <- as.integer(suppressWarnings(dbGetQuery(dbConn, "SELECT MAX(siteID) AS siteID FROM sites;")))
            stopifnot(newMaxSiteID == currentMaxSiteID + nrow(sites))
            currentMaxSiteID <- newMaxSiteID
        }
    }
    
    if(file.exists(paste0(file, "/multihitData.RData"))){
        load(paste0(file, "/multihitData.RData"))
        
        if(length(multihitData[[1]])>0){
            multihitPositions <- multihitData[[2]]
            multihitLengths <- multihitData[[3]]
            stopifnot(length(multihitPositions)==length(multihitLengths))
            multihitPositions <- data.frame(
                "multihitID"=rep(seq(length(multihitPositions))+currentMaxMultihitID,
                    sapply(multihitPositions, length)),
                "sampleID"=metadata[i,"sampleID"],
                "position"=start(flank(unlist(multihitPositions), width=-1, start=TRUE, both=FALSE)),
                "chr"=as.character(seqnames(unlist(multihitPositions))),
                "strand"=as.character(strand(unlist(multihitPositions))) )
            
            ## change to the right class as database
            multihitPositions$multihitID <- as(multihitPositions$multihitID, "integer")
            multihitPositions$sampleID <- as(multihitPositions$sampleID, "integer")
            multihitPositions$position <- as(multihitPositions$position, "integer")
            multihitPositions$chr <- as(multihitPositions$chr, "character")
            multihitPositions$strand <- as(multihitPositions$strand, "character")
            
            multihitLengths <- data.frame(
                "multihitID"=rep(seq(length(multihitLengths))+currentMaxMultihitID,
                    sapply(multihitLengths, nrow)),
                "length"=as.integer(as.character(do.call(rbind, multihitLengths)$Var1)),
                "count"=do.call(rbind, multihitLengths)$Freq )
            
            ## change to the right class as database
            multihitLengths$multihitID <- as(multihitLengths$multihitID, "integer")
            multihitLengths$length <- as(multihitLengths$length, "integer")
            multihitLengths$count <- as(multihitLengths$count, "integer")
            
            ## load table multihitpositions
            stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)
            message("Loading multihitpositions:", nrow(multihitPositions), " entries")
            stopifnot( dbWriteTable(dbConn, "multihitpositions", multihitPositions, append=T, row.names=F) )
            ## check loaded multihitpositions
            sql <- sprintf("SELECT * FROM multihitpositions WHERE multihitID>=%s AND multihitID<=%s",
                           range(multihitPositions$multihitID)[1],
                           range(multihitPositions$multihitID)[2])
            multihitpositions.from.db <- suppressWarnings( dbGetQuery(dbConn, sql) )
            multihitpositions.from.db$multihitID <- as(multihitpositions.from.db$multihitID, "integer")
            multihitpositions.from.db$sampleID <- as(multihitpositions.from.db$sampleID, "integer")
            multihitpositions.from.db$position <- as(multihitpositions.from.db$position, "integer")
            multihitpositions.from.db$chr <- as(multihitpositions.from.db$chr, "character")
            multihitpositions.from.db$strand <- as(multihitpositions.from.db$strand, "character")
            
            multihitpositions.from.db <- plyr::arrange(multihitpositions.from.db, multihitID,sampleID, position,chr,strand)
            multihitPositions <- plyr::arrange(multihitPositions, multihitID, sampleID,position,chr,strand)
            if(!identical(multihitPositions, multihitpositions.from.db)) {
                save.image("debug.rdata")
                stop("multihitPositions, multihitpositions.from.db not identical")
            }
            
            ## load table multihitlengths
            stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)
            message("Loading multihitlengths: ", nrow(multihitLengths), " entries")
            stopifnot( dbWriteTable(dbConn, "multihitlengths", multihitLengths, append=T, row.names=F) )
            ## check loaded multihitlengths
            sql <- sprintf("SELECT * FROM multihitlengths WHERE multihitID>=%s AND multihitID<=%s",
                           range(multihitPositions$multihitID)[1],
                           range(multihitPositions$multihitID)[2])
            multihitlengths.from.db <- suppressWarnings( dbGetQuery(dbConn, sql) )
            multihitlengths.from.db$multihitID <- as(multihitlengths.from.db$multihitID, "integer")
            multihitlengths.from.db$length <- as(multihitlengths.from.db$length, "integer")
            multihitlengths.from.db$count <- as(multihitlengths.from.db$count, "integer")
            
            multihitlengths.from.db <- plyr::arrange(multihitlengths.from.db, multihitID, length, count)
            multihitLengths <- plyr::arrange(multihitLengths, multihitID, length, count)
            if(!identical(multihitLengths, multihitlengths.from.db)) {
                save.image("debug.rdata")
                stop("multihitLengths, multihitlengths.from.db not identical")
            }
            
            newMaxMultihitID <- as.integer(suppressWarnings(dbGetQuery(dbConn, "SELECT MAX(multihitID) AS multihitID FROM multihitpositions;")))
            stopifnot(newMaxMultihitID == currentMaxMultihitID + length(unique(multihitPositions$multihitID)))
            currentMaxMultihitID <- newMaxMultihitID  
        }
    }
}

dbDiscon <- dbDisconnect(dbConn)

