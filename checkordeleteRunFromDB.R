## check for presence of R packages
rPackages <- c("stats", "RMySQL")
rPackagesPresent <- is.element(rPackages, installed.packages()[,1])
if(any(!rPackagesPresent)){
    stop(paste(rPackages[!rPackagesPresent]), " is not available")
}
stopifnot(all(sapply(rPackages, require, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)))

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

## initialize connection to database
## ~/.my.cnf must be present
junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn <- dbConnect(MySQL(), group="intSitesDev237") 
stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)

allSampleName <- suppressWarnings( dbGetQuery(dbConn, "SELECT * FROM samples") )
##names(allSampleName) <- tolower(names(allSampleName))

## get argument
args <- commandArgs(trailingOnly=TRUE)
runid <- args[1]
if( is.na(runid) ) {
    message("Usage:\n\tRscript deleteRunFromDB.R miseqrunid") 
    message("Available runs:\n", paste(unique(allSampleName$miseqid), collapse="\n"))
    stop()
}

if(! runid %in% allSampleName$miseqid) {
    message(runid, " is not in database")
    message("Available miseqids are:\n", paste(unique(allSampleName$miseqid), collapse="\n"))
    stop()
}

message("Checking run: ", runid)


#' increase output width to console width
wideScreen <- function(howWide=as.numeric(strsplit(system('stty size', intern=T), ' ')[[1]])[2]) {
   options(width=as.integer(howWide))
}
wideScreen()

#' get siteID and multihitID from sampleID
#' @param sampleID integer vector
#' @param verbose TRUE/FALSE [FALSE]
#' @return list(siteID=siteID, multihitID=multihitID) list of integer vectors
#' @example get_siteID_multihitID_from_sampleID(5)
get_siteID_multihitID_from_sampleID <- function(sampleID, verbose=FALSE) {
    ##sampleID=1080
    ## get sites.siteID and pcrbreakpoints.siteID
    stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)
    sql <- sprintf("SELECT siteID FROM sites WHERE sampleID='%s'", sampleID)
    if(verbose) message(sql)
    res <- suppressWarnings( dbGetQuery(dbConn, sql) )
    siteID <- res$siteID
    
    ## get multihitID for multihitpositions and multihitlengths
    sql <- sprintf("SELECT multihitID FROM multihitpositions WHERE sampleID='%s'", sampleID)
    if(verbose) message(sql)
    res <- suppressWarnings( dbGetQuery(dbConn, sql) )
    multihitID <- unique(res$multihitID)
    return(list(siteID=siteID, multihitID=multihitID))
}

#' check pcr break points and multihitlength from site ids
#' @param sitelist list(siteID=siteID, multihitID=multihitID) list of integer vectors
#' @param verbose TRUE/FALSE [FALSE]
#' @param delete  TRUE/FALSE [FALSE]
#' @return list(nsites, npcr, nmultievents, nmultisites, nmultilength)
#' @example check_pcr_multihitlength_from_siteids(get_siteID_multihitID_from_sampleID(5))
#' 
check_pcr_multihitlength_from_siteids <- function(sitelist, verbose=FALSE, delete=FALSE) {
    nsites <- npcr <- nmultievents <- nmultisites <- nmultilength <- NA
    if(verbose) message("===========")
    stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)
    if( length(sitelist$multihitID)>0 ) {
        nmultievents  <- length(unique(sitelist$multihitID))
        sitesin <- paste("(", paste(unique(sitelist$multihitID), collapse=","), ")")
        sql <- paste("SELECT * FROM multihitlengths WHERE multihitID in", sitesin)
        if(verbose) message(sql)
        res <- suppressWarnings( dbGetQuery(dbConn, sql) )
        nmultilength <- nrow(res)
        if(verbose) message(nmultievents, "\t", nrow(res))
        
        sql <- paste("SELECT * FROM multihitpositions WHERE multihitID in", sitesin)
        if(verbose) message(sql)
        res <- suppressWarnings( dbGetQuery(dbConn, sql) )
        nmultisites <- nrow(res)
        if(verbose) message(nmultievents, "\t", nrow(res))
        
        if(delete) {
            sqld <- paste("DELETE FROM multihitlengths WHERE multihitID in", sitesin)
            if(verbose) message(sqld)
            dbSendQuery(dbConn, sqld)
            sql <- paste("SELECT * FROM multihitlengths WHERE multihitID in", sitesin)
            res <- suppressWarnings( dbGetQuery(dbConn, sql) )
            if( nrow(res)>0 ) stop(sqld, "\nNot successful")
            
            sqld <- paste("DELETE FROM multihitpositions WHERE multihitID in", sitesin)
            if(verbose) message(sqld)
            dbSendQuery(dbConn, sqld)
            sql <- paste("SELECT * FROM multihitpositions WHERE multihitID in", sitesin)
            res <- suppressWarnings( dbGetQuery(dbConn, sql) )
            if( nrow(res)>0 ) stop(sqld, "\nNot successful")
        }
    }
    
    if( length(sitelist$siteID)>0 ) {
        sitesin <- paste("(", paste(sitelist$siteID, collapse=","), ")")
        sql <- paste("SELECT * FROM pcrbreakpoints WHERE siteID in", sitesin)
        if(verbose) message(sql)
        res <- suppressWarnings( dbGetQuery(dbConn, sql) )
        npcr <- nrow(res)
        if(verbose) message(length(sitelist$siteID), "\t", nrow(res))
        
        sql <- paste("SELECT * FROM sites WHERE siteID in", sitesin)
        if(verbose) message(sql)
        res <- suppressWarnings( dbGetQuery(dbConn, sql) )
        nsites <- nrow(res)
        if(verbose) message(length(sitelist$siteID), "\t", nrow(res))
        
        if(delete) {
            sqld <- paste("DELETE FROM pcrbreakpoints WHERE siteID in", sitesin)
            message(sqld)
            dbSendQuery(dbConn, sqld)
            sql <- paste("SELECT * FROM pcrbreakpoints WHERE siteID in", sitesin)
            res <- suppressWarnings( dbGetQuery(dbConn, sql) )
            if( nrow(res)>0 ) stop(sqld, "\nNot successful")
            
            sqld <- paste("DELETE FROM sites WHERE siteID in", sitesin)
            message(sqld)
            dbSendQuery(dbConn, sqld)
            sql <- paste("SELECT * FROM sites WHERE siteID in", sitesin)
            res <- suppressWarnings( dbGetQuery(dbConn, sql) )
            if( nrow(res)>0 ) stop(sqld, "\nNot successful")
        }
    }
    
    return(list(nsites=nsites,
                npcr=npcr,
                nmultievents=nmultievents,
                nmultisites=nmultisites,
                nmultilength=nmultilength))
}
##check_pcr_multihitlength_from_siteids(get_siteID_multihitID_from_sampleID(1080))

junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
dbConn <- dbConnect(MySQL(), group="intSitesDev237") 
stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)

toDelete <- subset(allSampleName, miseqid==runid)

siteID_multihitID <- lapply(setNames(toDelete$sampleID, toDelete$sampleID),
                            function(sid)
                                check_pcr_multihitlength_from_siteids(get_siteID_multihitID_from_sampleID(as.integer(sid))))

siteID_multihitID.n <- as.data.frame(do.call(rbind, siteID_multihitID))
siteID_multihitID.n$sampleID <- as.integer(rownames(siteID_multihitID.n))

toDelete <- merge(subset(allSampleName, miseqid==runid),
                  siteID_multihitID.n, by="sampleID",
                  all.x=TRUE)
message("The following sets will be deleted, siteID and multihitID are counts of sites.")
print(toDelete, row.names=FALSE)
##write.table(toDelete, "", sep = "\t", na="NA", row.names=FALSE, quote=FALSE)

## double check before proceed to delete
cat("enter yes to continue: ")
yes <- readLines(con="stdin", 1)
if( yes!="yes" ) q()
cat("enter yes to continue, this is the LAST warning: ")
yes <- readLines(con="stdin", 1)
if( yes!="yes" ) q()

## delete rows from sites, pcr, multi, multi
null <- lapply(setNames(toDelete$sampleID, toDelete$sampleID),
               function(sid)
                   check_pcr_multihitlength_from_siteids(get_siteID_multihitID_from_sampleID(as.integer(sid)), delete=TRUE))

## delete rows from samples
if( length(toDelete$sampleID)>0 ) {
    junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
    dbConn <- dbConnect(MySQL(), group="intSitesDev237") 
    stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)
    
    sitesin <- paste("(", paste(toDelete$sampleID, collapse=","), ")")
    sqld <- paste("DELETE FROM samples WHERE sampleID in", sitesin)
    message(sqld)
    dbSendQuery(dbConn, sqld)
    sql <- paste("SELECT * FROM samples WHERE sampleID in", sitesin)
    res <- suppressWarnings( dbGetQuery(dbConn, sql) )
    if( nrow(res)>0 ) stop(sqld, "\nNot successful")
}


