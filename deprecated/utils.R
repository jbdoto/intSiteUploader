#' create srq_mysql without specifying params in ctor
create_src_mysql <- function(con) {
    info <- dbGetInfo(con)
    src_sql("mysql", con, info = info)
}

check_presence_packages <- function() {
    rPackages <- c("intSiteRetriever", "stats", "RMySQL", "GenomicRanges", 
        "BiocGenerics", "parallel", "IRanges", "GenomeInfoDb")
    rPackagesPresent <- is.element(rPackages, installed.packages()[,1])
    if(any(!rPackagesPresent)){
        stop(paste(rPackages[!rPackagesPresent]), " is not available")
    }
    stopifnot(sapply(rPackages, require, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE))
}

check_presence_command_line_tools <- function() {
    commandLinePrograms <- c("mysql")
    programsPresent <- !sapply(sprintf("which %s > /dev/null 2>&1", commandLinePrograms), system)
    if(any(!programsPresent)){
      stop(paste(commandLinePrograms[!programsPresent]), " is not available")
    }
}
