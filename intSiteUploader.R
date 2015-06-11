#check for presence of R packages
rPackages <- c("stats", "RMySQL")
rPackagesPresent <- is.element(rPackages, installed.packages()[,1])
if(any(!rPackagesPresent)){
  stop(paste(rPackages[!rPackagesPresent]), " is not available")
}

library("stats")
library("RMySQL")
options(stringsAsFactors=F)

#check for presence of command line stuff
commandLinePrograms <- c("mysql")
programsPresent <- !sapply(sprintf("which %s > /dev/null 2>&1", commandLinePrograms), system)
if(any(!programsPresent)){
  stop(paste(commandLinePrograms[!programsPresent]), " is not available")
}

#working directory (i.e. primary analysis directory) is passed in via command line

args <- commandArgs(trailingOnly=TRUE)
workingDir <- args[1]
stopifnot(!is.na(workingDir))
workingDir <- normalizePath(".", mustWork=TRUE)
setwd(workingDir)

#can't rely on completeMetadata.RData as it may have been rm'd in cleanup
stopifnot(all(file.exists("sampleInfo.tsv", "processingParams.tsv")))

metadata <- read.table("sampleInfo.tsv", header=TRUE, stringsAsFactors=F)
processingParams <- read.table("processingParams.tsv", header=TRUE, stringsAsFactors=F)
junk <- merge(metadata, processingParams)
stopifnot(nrow(metadata) == nrow(junk))
metadata <- junk

metadata <- metadata[c("alias", "gender", "refGenome")]
names(metadata) <- c("sampleName", "gender", "refGenome")

all_cons <- dbListConnections(MySQL())
for (con in all_cons) {
  discCon <- dbDisconnect(con)
}
stopifnot(file.exists("~/.my.cnf"))
dbConn <- dbConnect(MySQL(), group="intSitesDev237") #~/.my.cnf must be present

#this isn't pretty, but it does the job, especially since we're not going to be loading in tons of samples at once
alreadyLoaded <- dbGetQuery(dbConn, paste0("SELECT DISTINCT sampleName
                                           FROM samples
                                           WHERE sampleName REGEXP ", dbQuoteString(dbConn, paste0(paste0("^", metadata$sampleName, "$"), collapse="|")), ";"))

if(nrow(alreadyLoaded) > 0){
  stop(paste0("sets already exist in the database: ", paste(alreadyLoaded$sampleName, collapse=", ")))
}

#assumes at least one sample is loaded into the DB
currentMaxSampleID <- as.integer(suppressWarnings(dbGetQuery(dbConn, "SELECT MAX(sampleID) AS sampleID FROM samples;")))
if(is.na(currentMaxSampleID)) currentMaxSampleID<-0

metadata$sampleID <- seq(nrow(metadata))+currentMaxSampleID

stopifnot( dbWriteTable(dbConn, "samples", metadata, append=T, row.names=F) )

#assumes at least one sample is already loaded into the DB
currentMaxSiteID <- as.integer(suppressWarnings(dbGetQuery(dbConn, "SELECT MAX(siteID) AS siteID FROM sites;")))
if(is.na(currentMaxSiteID)) currentMaxSiteID<-0

currentMaxMultihitID <- as.integer(suppressWarnings(dbGetQuery(dbConn, "SELECT MAX(multihitID) AS multihitID FROM multihitpositions;")))
if(is.na(currentMaxMultihitID)) currentMaxMultihitID<-0

for(i in seq(nrow(metadata))){
    file <- metadata[i,"sampleName"]
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
            
            stopifnot( dbWriteTable(dbConn, "sites", sites, append=T, row.names=F) )
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
            
            stopifnot( dbWriteTable(dbConn, "multihitpositions", multihitPositions, append=T, row.names=F) )
            stopifnot( dbWriteTable(dbConn, "multihitlengths", multihitLengths, append=T, row.names=F) )
        }
    }
}

dbDiscon <- dbDisconnect(dbConn)

