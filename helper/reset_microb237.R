#' This script reset the test database on microb237 according to the intsitesdev.sql file
#' This script can be executed from any directory
#' Rscript path/to/reset_microb237.R 

#### get codeDir ####
args <- commandArgs(trailingOnly = FALSE)
codeDir <- dirname(sub("--file=", "", grep("--file=", args, value=T)))
if ( length(codeDir) == 0 ) codeDir <- "/home/yinghua/intSiteUploader/helper"
codeDir <- normalizePath(file.path(codeDir, ".."))
stopifnot(length(list.files(path=codeDir, pattern="intSiteUploader.R$"))==1)

#### create intsitesdevtest on microb237 using intsitesdev.sql ####
test_db_cnf <- "~/.my.cnf.test"
test_db_cnf <- normalizePath(test_db_cnf, mustWork=TRUE)

if( any(grepl("microb98", readLines(test_db_cnf))) ) stop(
    "File ", test_db_cnf, " should not contain microb98 for safety measure" )

cmd <- sprintf("mysql --defaults-file=%s -e 'CREATE DATABASE IF NOT EXISTS intsitesdevtest' ", test_db_cnf)
message(cmd)
stopifnot( system(cmd)==0 )

cmd <- sprintf("mysql --defaults-file=%s intsitesdevtest < %s/intsitesdev.sql", test_db_cnf, codeDir)
message(cmd)
stopifnot( system(cmd)==0 )

