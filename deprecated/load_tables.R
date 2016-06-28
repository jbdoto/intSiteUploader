write_table_samples <- function(dbConn, metadata) {
    ## Get max sampleID, and start from max+1
    currentMaxSampleID <- as.integer(dbGetQuery(dbConn, "SELECT MAX(sampleID) AS sampleID FROM samples;"))
    if(is.na(currentMaxSampleID)) { # empty DB
        currentMaxSampleID<-10000
    }

    ## load table samples
    metadata$sampleID <- seq(nrow(metadata))+currentMaxSampleID
    stopifnot( dbWriteTable(dbConn, "samples", metadata, append=T, row.names=F) )
    metadata
}

check_write_table_samples <- function(dbConn, metadata) {
    ## check wether load was successful
    sample.tab <- suppressWarnings(dbReadTable(dbConn, "samples"))
    merged.tab <- merge(metadata, sample.tab, by=c("sampleName", "refGenome"), all.x=TRUE)
    if( !all(merged.tab$sampleID.x==merged.tab$sampleID.y) ) {
        message("Sample ID error, check the following table")
        print(merged.tab)
    }
}
