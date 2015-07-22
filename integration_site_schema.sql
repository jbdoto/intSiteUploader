DROP TABLE IF EXISTS samples;
CREATE TABLE samples (
    sampleID int NOT NULL,
    sampleName varchar(255) NOT NULL,
    refGenome varchar(10) NOT NULL,
    gender char(1) NOT NULL,
    runId varchar(255),
    PRIMARY KEY (sampleID),
    CONSTRAINT uniq_samples UNIQUE (sampleName, refGenome)
);

-- unique hits
DROP TABLE IF EXISTS sites;
CREATE TABLE sites (
    siteID int NOT NULL,
    sampleID int NOT NULL,
    position int NOT NULL,
    chr varchar(255) NOT NULL,
    strand char(1) NOT NULL,
    PRIMARY KEY (siteID),
    FOREIGN KEY (sampleID) REFERENCES samples(sampleID)
);

DROP TABLE IF EXISTS pcrbreakpoints;
CREATE TABLE `pcrbreakpoints` (
    siteID int NOT NULL,
    breakpoint int NOT NULL,
    count int NOT NULL,
    PRIMARY KEY (siteID, breakpoint),
    FOREIGN KEY (siteID) REFERENCES sites(siteID) 
);

-- multihit schema 
DROP TABLE IF EXISTS multihitpositions;
CREATE TABLE multihitpositions (
    multihitID int NOT NULL,
    sampleID int NOT NULL,
    position int NOT NULL,
    chr varchar(255) NOT NULL,
    strand char(1) NOT NULL,
    PRIMARY KEY (multihitID, position, chr, strand),
    FOREIGN KEY (sampleID) REFERENCES samples(sampleID)
);

DROP TABLE IF EXISTS multihitlengths;
CREATE TABLE multihitlengths (
    multihitID int NOT NULL,
    length int NOT NULL,
    count int NOT NULL,
    PRIMARY KEY (multihitID, length),
    FOREIGN KEY (multihitID) REFERENCES multihitpositions(multihitID) 
);
