library(RMySQL, quietly=TRUE, verbose=FALSE)
library(dplyr, quietly=TRUE, verbose=FALSE)

codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))

source(file.path(codeDir, "utils.R"))
source(file.path(codeDir, "load_tables.R"))

check_presence_packages()
options(stringsAsFactors=F)

## check if file exist and permission .my.cnf 
stopifnot(file.exists("~/.my.cnf"))
stopifnot(file.info("~/.my.cnf")$mode == as.octmode("600"))

check_presence_command_line_tools()

## working directory (i.e. primary analysis directory) is passed in via command line
# and group for MySQL server
args <- commandArgs(trailingOnly=TRUE)
workingDir <- args[1]
mysql_group <- "intsites_miseq"
if ( ! is.na(args[2])) {
    mysql_group <- args[2]
} 
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
stopifnot(file.exists("completeMetadata.RData"))
metadata <- get(load('completeMetadata.RData'))
metadata <- subset(metadata, select=c("alias", "gender", "refGenome"))
names(metadata) <- c("sampleName", "gender", "refGenome")

## initialize connection to database
## ~/.my.cnf must be present
dbConn <- dbConnect(MySQL(), group=mysql_group) 

## stop if any sample is already loaded
read_conn <- create_src_mysql(dbConn)
is.loaded <- setNameExists(select(metadata, sampleName, refGenome), read_conn)
if(any(is.loaded)) message(
    paste0("Sets already in the database: ",
           paste(metadata[is.loaded, ], collapse="\n")))

if( any(grepl("^GTSP", metadata$sampleName[is.loaded], ignore.case=TRUE)) ) {
    stop("GTSP sample already loaded, delete from the database or leave them alone")
}

metadata <- subset(metadata, ! is.loaded)
if (nrow(metadata) == 0) {
    message("All samples are already in the DB. Nothing is pushed into DB.")
    q()
}
metadata$miseqid <- miseqid

dbGetQuery(dbConn, "START TRANSACTION;")
metadata <- write_table_samples(dbConn, metadata)

## Get max siteID, and start from max+1
currentMaxSiteID <- as.integer(dbGetQuery(dbConn, "SELECT MAX(siteID) AS siteID FROM sites;"))
if(is.na(currentMaxSiteID)) { 
    currentMaxSiteID<-100000000
}

## Get max MultihitID, and start from max+1
currentMaxMultihitID <- as.integer(dbGetQuery(dbConn, "SELECT MAX(multihitID) AS multihitID FROM multihitpositions;"))
if(is.na(currentMaxMultihitID)) {
    currentMaxMultihitID<-100000000
}

## process by sample and upload to sites, pcrbreakpoints, multihitpositions, multihitlengths
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
            message("Loading sites: ", nrow(sites), " entries")
            stopifnot( dbWriteTable(dbConn, "sites", sites, append=T, row.names=F) )
            
            ## load table pcrbreakpoints
            message("Loading pcrbreakpoints: ", nrow(pcrBreakpoints), " entries")
            stopifnot( dbWriteTable(dbConn, "pcrbreakpoints", pcrBreakpoints, append=T, row.names=F) )
            
            newMaxSiteID = currentMaxSiteID + nrow(sites)
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
            message("Loading multihitpositions:", nrow(multihitPositions), " entries")
            stopifnot( dbWriteTable(dbConn, "multihitpositions", multihitPositions, append=T, row.names=F) )
            
            ## load table multihitlengths
            message("Loading multihitlengths: ", nrow(multihitLengths), " entries")
            stopifnot( dbWriteTable(dbConn, "multihitlengths", multihitLengths, append=T, row.names=F) )
            newMaxMultihitID = currentMaxMultihitID + length(unique(multihitPositions$multihitID))
            currentMaxMultihitID <- newMaxMultihitID  
        }
    }
}

message("\nCOMMIT\n")
dbGetQuery(dbConn, "COMMIT;")

check_write_table_samples(dbConn, metadata)

dbDiscon <- dbDisconnect(dbConn)

