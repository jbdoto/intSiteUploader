library("stats")

#working directory (i.e. primary analysis directory) is passed in via command line

args <- commandArgs(trailingOnly = TRUE)
workingDir = args[1]
stopifnot(!is.na(workingDir))
setwd(workingDir)

indexPath = get(load("indexPath.RData")) #already in the working directory
genomeName = strsplit(indexPath, ".2bit")[[1]][1]
genomeName = strsplit(genomeName, "/")[[1]][length(strsplit(genomeName, "/")[[1]])]

library("RMySQL") #also loads DBI
all_cons <- dbListConnections(MySQL())
for (con in all_cons) {
  discCon <- dbDisconnect(con)
}
dbConn = dbConnect(MySQL(), group="intSitesDEV-dev") #~/.my.cnf must be present 

filesToLoad = system("ls */sites.final.RData", intern=T)
filesToLoad = sapply(strsplit(filesToLoad, "/"), "[[" ,1)

#this isn't pretty, but it does the job, especially since we're not going to be loading in tons of samples at once
alreadyLoaded = dbGetQuery(dbConn, paste0("SELECT DISTINCT sampleName FROM samples WHERE sampleName IN (", paste0(dbQuoteString(dbConn, filesToLoad), collapse=","), ")"))

if(nrow(alreadyLoaded) > 0){
  stop(paste0("sets already exist in the database: ", paste(alreadyLoaded$sampleName, collapse=", ")))
}

#assumes at least one sample is loaded into the DB
currentMaxSampleID = as.integer(suppressWarnings(dbGetQuery(dbConn, "SELECT MAX(sampleID) AS sampleID FROM samples;")))

filesToLoad = data.frame("sampleID"=seq(length(filesToLoad))+currentMaxSampleID, "sampleName"=filesToLoad, "refGenome"=genomeName)

dbWriteTable(dbConn, "samples", filesToLoad, append=T, row.names=F)

#assumes at least one sample is loaded into the DB
currentMaxSiteID = as.integer(suppressWarnings(dbGetQuery(dbConn, "SELECT MAX(siteID) AS siteID FROM sites;")))

for(i in c(1:nrow(filesToLoad))){
  file = filesToLoad[i,"sampleName"]
  load(paste0(file, "/sites.final.RData"))
  load(paste0(file, "/allSites.RData"))
  load(paste0(file, "/multihitData.RData"))
  
  #sites.final won't exist if there aren't sites, thus no need to check if sites.final has sites in it
  sites = data.frame("sampleID"=filesToLoad[i,"sampleID"],
                     "siteID" = seq(length(sites.final))+currentMaxSiteID,
                     "position"=sites.final$intLoc,
                     "chr"=as.character(seqnames(sites.final)),
                     "strand"=as.character(strand(sites.final)),
                     "multihitID"=NA)
  
  currentMaxSiteID = currentMaxSiteID + nrow(sites)
  
  allSites = allSites[unlist(sites.final$revmap)]

  #could do the next three statements with aggregate, but this method is emperically 2x faster
  pcrBreakpoints = sort(paste0(as.integer(Rle(sites$siteID, sapply(sites.final$revmap, length))), ".", width(allSites)))
    
  condensedPCRBreakpoints = strsplit(unique(pcrBreakpoints), "\\.")
      
  pcrBreakpoints = data.frame("siteID" = sapply(condensedPCRBreakpoints, "[[", 1),
                              "distToBreakpoint" = sapply(condensedPCRBreakpoints, "[[", 2),
                              "counts" = runLength(Rle(match(pcrBreakpoints, unique(pcrBreakpoints)))))
    
  multihits = unlist(multihitData[[2]], use.names=F)
  
  if(length(multihits)>0){ #multihits could be empty though...
    multihits = data.frame("sampleID"=filesToLoad[i,"sampleID"],
                           "siteID"=seq(length(multihits))+currentMaxSiteID,
                           "position"=start(flank(multihits, width=-1, start=TRUE, both=FALSE)),
                           "chr"=as.character(seqnames(multihits)),
                           "strand"=as.character(strand(multihits)),
                           "multihitID"=multihits$ID)
  }else{
    multihits = data.frame()
  }
  
  currentMaxSiteID = currentMaxSiteID + nrow(multihits)  
  
  #these are processed as transactions by RMySQL under the hood
  dbWriteTable(dbConn, "sites", rbind(sites, multihits), append=T, row.names=F) 
  dbWriteTable(dbConn, "pcrbreakpoints", pcrBreakpoints, append=T, row.names=F)
  #is this^^^ safe or do I need to chunk it? Might be 1000's of rows for some samples
  
}

dbDiscon = dbDisconnect(dbConn)
  