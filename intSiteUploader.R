library(RMySQL, quietly=TRUE, verbose=FALSE)
library(dplyr, quietly=TRUE, verbose=FALSE)
library(DBI, quietly=TRUE, verbose=FALSE)
library(yaml, quietly=TRUE, verbose=FALSE)
options(stringsAsFactors=F, useFancyQuotes=F)

codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))

# Primary analysis directory is passed in via command line.
args <- commandArgs(trailingOnly=TRUE)
workingDir <- args[1]


# Determine analysis directory and change into it.
if( is.na(workingDir) ) workingDir <- "."
workingDir <- normalizePath(workingDir, mustWork=TRUE)
setwd(workingDir)
message("Changed to directory: ", workingDir)

# Load configuration file
config <<- yaml.load_file(paste0(workingDir, '/INSPIIRED.yml'))

# Get the sample information.
stopifnot(file.exists("completeMetadata.RData"))
metadata <- get(load('completeMetadata.RData'))
metadata <- subset(metadata, select=c("alias", "gender", "refGenome"))
names(metadata) <- c("sampleName", "gender", "refGenome")


# Connect to my database
if (config$UseMySQL){
   stopifnot(file.exists("~/.my.cnf"))
   stopifnot(file.info("~/.my.cnf")$mode == as.octmode("600"))
   dbConn <- dbConnect(MySQL(), group=config$MySQLconnectionGroup)
}else{
   dbConn <- dbConnect(RSQLite::SQLite(), dbname=config$SQLiteIntSiteCallerDB)
}

samples <- dbGetQuery(dbConn, "select * from samples")

dupRunsTest <- function(x, samples=samples){
  t <- samples[samples$sampleName==x[1] & samples$refGenome==x[3],]
  if (nrow(t)>0) x[1]
}

# Test is a character vector with names of samples already in the sample database table
test <- apply(metadata, 1, dupRunsTest, samples=samples)

if ((nrow(metadata)) == (length(test))){
    message('All of the samples have already been uploaded to the data database')
    q()
}

if ( length(test) > 0 ){
  message('The following samples are already in the database:')
  for (i in names(test)) {
     message(i)
  } 
  q()
}


# Add a miseq column to the metadata table
stopifnot(!is.null(config$runId))
metadata$miseqid <- config$runId


insertSQL <- function (dbConn, x, colNames, table){
   fields <- paste(colNames, collapse=", ")
   values <- paste(sQuote(x), collapse=", ")
   q <- paste0('insert into ', table, ' (', fields, ') values (',  values, ')')
   r <- dbSendQuery(dbConn, q)
}


write_table_samples <- function(dbConn, metadata) {
   
    currentMaxSampleID <- as.integer(dbGetQuery(dbConn, "SELECT MAX(sampleID) AS sampleID FROM samples;"))

    if(is.na(currentMaxSampleID)) currentMaxSampleID <- 10000  # empty database condition

    # Create a column of sequence ids and append itto the metadata table.
    metadata$sampleID <- seq(nrow(metadata))+currentMaxSampleID

    null <- apply(metadata, 1, insertSQL, dbConn=dbConn, colNames=colnames(metadata), table='samples')
    metadata
}

check_write_table_samples <- function(dbConn, metadata) {
    sample.tab <- suppressWarnings(dbReadTable(dbConn, "samples"))
    merged.tab <- merge(metadata, sample.tab, by=c("sampleName", "refGenome"), all.x=TRUE)

    if( !all(merged.tab$sampleID.x==merged.tab$sampleID.y) ) {
        message("Sample ID error, check the following table")
        print(merged.tab)
    }
}


# SQLite is autocommit by deafult, no need to work with transaction sessions.
if (config$UseMySQL) dbGetQuery(dbConn, "START TRANSACTION;")

metadata <- write_table_samples(dbConn, metadata)


# Get max siteID, and start from max+1
currentMaxSiteID <- as.integer(dbGetQuery(dbConn, "SELECT MAX(siteID) AS siteID FROM sites;"))
if(is.na(currentMaxSiteID)) { 
    currentMaxSiteID<-100000000
}

# Get max MultihitID, and start from max+1
currentMaxMultihitID <- as.integer(dbGetQuery(dbConn, "SELECT MAX(multihitID) AS multihitID FROM multihitpositions;"))
if(is.na(currentMaxMultihitID)) {
    currentMaxMultihitID<-100000000
}

# process by sample and upload to sites, pcrbreakpoints, multihitpositions, multihitlengths
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
            
            # change to the right class as database
            pcrBreakpoints$siteID <- as(pcrBreakpoints$siteID, "integer")
            pcrBreakpoints$breakpoint <- as(pcrBreakpoints$breakpoint, "integer")
            pcrBreakpoints$count <- as(pcrBreakpoints$count, "integer")
            
            # load table sites            
            message("Loading sites: ", nrow(sites), " entries")
            ### stopifnot( dbWriteTable(dbConn, "sites", sites, append=T, row.names=F) )
            null <- apply(sites, 1, insertSQL, dbConn=dbConn, colNames=colnames(sites), table='sites')            

            # load table pcrbreakpoints
            message("Loading pcrbreakpoints: ", nrow(pcrBreakpoints), " entries")
            ### stopifnot( dbWriteTable(dbConn, "pcrbreakpoints", pcrBreakpoints, append=T, row.names=F) )
            null <- apply(pcrBreakpoints, 1, insertSQL, dbConn=dbConn, colNames=colnames(pcrBreakpoints), table='pcrbreakpoints')            

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
            
            # change to the right class as database
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
            
            # change to the right class as database
            multihitLengths$multihitID <- as(multihitLengths$multihitID, "integer")
            multihitLengths$length <- as(multihitLengths$length, "integer")
            multihitLengths$count <- as(multihitLengths$count, "integer")
            
            # load table multihitpositions
            message("Loading multihitpositions:", nrow(multihitPositions), " entries")
            ### stopifnot( dbWriteTable(dbConn, "multihitpositions", multihitPositions, append=T, row.names=F) )
            null <- apply(multihitPositions, 1, insertSQL, dbConn=dbConn, colNames=colnames(multihitPositions), table='multihitPositions')            

            # load table multihitlengths
            message("Loading multihitlengths: ", nrow(multihitLengths), " entries")
            ### stopifnot( dbWriteTable(dbConn, "multihitlengths", multihitLengths, append=T, row.names=F) )
            null <- apply(multihitLengths, 1, insertSQL, dbConn=dbConn, colNames=colnames(multihitLengths), table='multihitLengths')

            newMaxMultihitID = currentMaxMultihitID + length(unique(multihitPositions$multihitID))
            currentMaxMultihitID <- newMaxMultihitID  
        }
    }
}

if (config$UseMySQL) dbGetQuery(dbConn, "COMMIT;")

check_write_table_samples(dbConn, metadata)
dbDiscon <- dbDisconnect(dbConn)
