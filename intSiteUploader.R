## check for presence of R packages
rPackages <- c("stats", "RMySQL", "GenomicRanges", "BiocGenerics", "parallel")
rPackagesPresent <- is.element(rPackages, installed.packages()[,1])
if(any(!rPackagesPresent)){
    stop(paste(rPackages[!rPackagesPresent]), " is not available")
}
stopifnot(sapply(rPackages, require, character.only=TRUE))

library("stats")
library("RMySQL")
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

## get sample information
stopifnot(all(file.exists("sampleInfo.tsv", "processingParams.tsv")))
metadata <- read.table("sampleInfo.tsv", header=TRUE, stringsAsFactors=F)
processingParams <- read.table("processingParams.tsv", header=TRUE, stringsAsFactors=F)
junk <- merge(metadata, processingParams)
stopifnot(nrow(metadata) == nrow(junk))
metadata <- junk
metadata <- metadata[c("alias", "gender", "refGenome")]
names(metadata) <- c("sampleName", "gender", "refGenome")

## initialize connection to database
## ~/.my.cnf must be present
junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn <- dbConnect(MySQL(), group="intSitesDev237") 

## stop if any sample is already loaded
allSampleName <- dbGetQuery(dbConn, "SELECT DISTINCT sampleName FROM samples")
is.loaded <- metadata$sampleName %in% allSampleName$sampleName
if(any(is.loaded)) stop(
    paste0("Sets already in the database: ",
           paste(metadata$sampleName[is.loaded], collapse="\n")))


## Get max sampleID, and start from max+1
currentMaxSampleID <- as.integer(suppressWarnings(dbGetQuery(dbConn, "SELECT MAX(sampleID) AS sampleID FROM samples;")))
if(is.na(currentMaxSampleID)) currentMaxSampleID<-0

metadata$sampleID <- seq(nrow(metadata))+currentMaxSampleID

## load table samples
stopifnot( dbWriteTable(dbConn, "samples", metadata, append=T, row.names=F) )

## Get max siteID, and start from max+1
currentMaxSiteID <- as.integer(suppressWarnings(dbGetQuery(dbConn, "SELECT MAX(siteID) AS siteID FROM sites;")))
if(is.na(currentMaxSiteID)) currentMaxSiteID<-0

## Get max MultihitID, and start from max+1
currentMaxMultihitID <- as.integer(suppressWarnings(dbGetQuery(dbConn, "SELECT MAX(multihitID) AS multihitID FROM multihitpositions;")))
if(is.na(currentMaxMultihitID)) currentMaxMultihitID<-0

## process by sample
for(i in seq(nrow(metadata))){
    file <- metadata[i,"sampleName"]
    message("Processing: ", file)
    if(file.exists(paste0(file, "/sites.final.RData"), paste0(file, "/allSites.RData"))){
        load(paste0(file, "/sites.final.RData"))
        load(paste0(file, "/allSites.RData"))
        if(length(sites.final)>0){
            ##sites.final won't exist if there aren't sites, thus no need to check if sites.final has sites in it
            sites <- data.frame("sampleID"=metadata[i,"sampleID"],
                                "siteID"=seq(length(sites.final))+currentMaxSiteID,
                                "position"=start(flank(sites.final, -1, start=T)),
                                "chr"=as.character(seqnames(sites.final)),
                                "strand"=as.character(strand(sites.final)))
            
            currentMaxSiteID <- currentMaxSiteID + nrow(sites)
            
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
            
            ## load table samples            
            stopifnot( dbWriteTable(dbConn, "sites", sites, append=T, row.names=F) )
            ## load table pcrbreakpoints
            stopifnot( dbWriteTable(dbConn, "pcrbreakpoints", pcrBreakpoints, append=T, row.names=F) )
        }
    }
    
    if(file.exists(paste0(file, "/multihitData.RData"))){
        load(paste0(file, "/multihitData.RData"))
        
        if(length(multihitData[[1]])>0){
            multihitPositions <- multihitData[[2]]
            multihitLengths <- multihitData[[3]]
            stopifnot(length(multihitPositions)==length(multihitLengths))
            multihitPositions <- data.frame("sampleID"=metadata[i,"sampleID"],
                                            "multihitID"=rep(seq(length(multihitPositions))+currentMaxMultihitID,
                                                sapply(multihitPositions, length)),
                                            "position"=start(flank(unlist(multihitPositions), width=-1, start=TRUE, both=FALSE)),
                                            "chr"=as.character(seqnames(unlist(multihitPositions))),
                                            "strand"=as.character(strand(unlist(multihitPositions))))
            
            multihitLengths <- data.frame("multihitID"=rep(seq(length(multihitLengths))+currentMaxMultihitID,
                                              sapply(multihitLengths, nrow)),
                                          "length"=do.call(rbind, multihitLengths)$Var1,
                                          "count"=do.call(rbind, multihitLengths)$Freq)
            
            currentMaxMultihitID <- currentMaxMultihitID + length(multihitPositions)  
            
            ## load table multihitpositions
            stopifnot( dbWriteTable(dbConn, "multihitpositions", multihitPositions, append=T, row.names=F) )
            ## load table multihitlengths
            stopifnot( dbWriteTable(dbConn, "multihitlengths", multihitLengths, append=T, row.names=F) )
        }
    }
}

dbDiscon <- dbDisconnect(dbConn)

