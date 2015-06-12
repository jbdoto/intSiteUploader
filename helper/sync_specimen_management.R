#' This script sync specimen_management from microb98 to microb237
#' @param database database to sync
#' @param microb98_db_cnf source database config file
#' @param test_db_cnf destination database config file
#' @example Rscript sync_specimen_management.R

## database to sync
database <- "specimen_management"

## source database server
microb98_db_cnf <- "~/.my.cnf.microb98"
microb98_db_cnf <- normalizePath(microb98_db_cnf, mustWork=TRUE)
stopifnot(file.info(microb98_db_cnf)$mode == as.octmode("600"))

## destination database server
test_db_cnf <- "~/.my.cnf.test"
test_db_cnf <- normalizePath(test_db_cnf, mustWork=TRUE)
stopifnot(file.info(test_db_cnf)$mode == as.octmode("600"))
if( any(grepl("microb98", readLines(test_db_cnf))) ) stop(
    "File ", test_db_cnf, " should not contain microb98 for safety measure" )


## check existance of source database
cmd <- sprintf("mysql --defaults-file=%s -Be \"SHOW DATABASES LIKE '%s'\"", microb98_db_cnf, database)
message(cmd)
if( !all(grepl(database, system(cmd, intern=TRUE))) ) stop(
                                         "Database ", database, " not exist")

## dump source database
cmd <- sprintf("mysqldump --defaults-file=%s --single-transaction %s > %s.sql", microb98_db_cnf, database, database)
message(cmd)
stopifnot( system(cmd)==0 )


## create database in destination if not exist
cmd <- sprintf("mysql --defaults-file=%s -e 'CREATE DATABASE IF NOT EXISTS %s' ", test_db_cnf, database)
message(cmd)
stopifnot( system(cmd)==0 )

## import into destination
cmd <- sprintf("mysql --defaults-file=%s %s < %s.sql", test_db_cnf, database, database)
message(cmd)
stopifnot( system(cmd)==0 )

