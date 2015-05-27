# intSiteUploader

***


## Introduction
This code is designed to upload integration sites, PCR breakpoints, and multihits identified by `intSiteCaller` and upload them to the intsitesdev database.  The database is currently located at `microb98.med.upenn.edu:3309` and is described by the included schema, `insitesdev.sql`.


## Inputs

`intSiteUploader` expects to be passed the `primaryAnalysisDirectory` as defined when `intSiteCaller` was run.  This directory must have the following minimal file structure:

```
primaryAnalysisDirectory
├── sample1
│   ├── sites.final.RData
│   ├── multihitData.RData
│   └── allSites.RData
├── processingParams.csv
└── sampleInfo.csv
```

There can be as few or as many samples as the user desires in the `primaryAnalysisDirectory`, so long as each sample is represented in both `processingParams.csv` and `sampleInfo.csv`.  See [intSiteCaller's Documentation](http://www.github.com/esherm/intSiteCaller) for a description of the values contained in these two metadata files. 

## Usage

Uploading is started by running ```Rscript intSiteUploader.R <primaryAnalysisDir>```.

**At this time, it is highly recommmended to only run one instance of `intSiteUploader` at a time in order to avoid primary key collisions.**  Uploading an entire MiSeq run should take less than 2 minutes, so this restriction should be functionally inconsequential.

**DO NOT run `intSiteUploader` after `intSiteCaller` by using the `&&` command line operator.** `intSiteCaller.R` is only a 'kickoff' script and will exit with a status of 1 before the full analysis has completed.  It is highly recommended that the user manually confirm that the analysis succeeded and that vaguely reasonable results have been returned before running `intSiteUploader`.  There is currently no way to automatically purge a run from the database (not to mention this is poor practice in a production environment), so make sure you're happy with it before uploading!

## Outputs

`intSiteUploader.R` does not directly output any commands or data, as all MySQL commands are generated and executed internally.


## Dependencies

This code is dependent on the `RMySQL` and `stats` R packages as well as the `mysql` executable.

Additionally, the database connection information is expected to be stored in `~/.my.cnf` in the following format:

```
[intSitesDev]
user=foo
password=bar
host=microb98.med.upenn.edu
port=3309
database=intsitesdev
```

`user` and `password` coorespond to the account you wish to use to log into the MySQL database.

