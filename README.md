# intSiteUploader

***


## Introduction
This code is designed to upload integration sites, PCR breakpoints, and multihits identified by `intSiteCaller` and upload them to the intsitesdev database.  The database is currently located at `microbxxx.med.upenn.edu:3306` and is described by the included schema, `insitesdev.sql`.


## Inputs

`intSiteUploader` expects to be passed the `primaryAnalysisDirectory` as defined when `intSiteCaller` was run.  This directory must have the following minimal file structure:

```
primaryAnalysisDirectory
├── sample1
│   ├── sites.final.RData
│   ├── multihitData.RData
│   └── allSites.RData
├─miseqid.txt
└─completeMetadata.Rdata
```

There can be as few or as many samples as the user desires in the
`primaryAnalysisDirectory`, so long as each sample is represented in `completeMetadata.Rdata`
.  See [intSiteCaller's
Documentation](http://www.github.com/esherm/intSiteCaller) for a description of
the values contained in these two metadata files. 

## Usage
Code example:
```
cd run20150505                                # a recent processed run folder
Rscript path/to/intSiteUploader.R .
Rscript intSiteUploader.R <primaryAnalysisDir> [mysql_group]
```

at present default group for `mysql_group` is `intSitesDev237`

Note:
* Run intSiteUploader.R only after running intSiteCaller,
* Only run one instance at a time,
* Code to initialize the database is located in the `helper` folder.

## Outputs

* `intSiteUploader.R` does not directly output any commands or data, as all MySQL commands are generated and executed internally.
* A script to check that the data is correctly uploaded is to be coded later.

## Dependencies

This code is dependent on the `RMySQL` and `stats` R packages as well as the `mysql` executable.

Additionally, the database connection information is expected to be stored in `~/.my.cnf` in the following format:

```
[intSitesDev]
user=foo
password=bar
host=microbxxx.med.upenn.edu
port=3306
database=intsitesdev
```

`user` and `password` coorespond to the account you wish to use to log into the MySQL database.

# Run identification

To distinguish different runs we read
run id from `myseqid.txt` that should be present in primary analysis folder.



## Testing 

sqllite version of db can be created by:

```
sqlite3 db.sqlite3 < integration_site_schema.sql
```
