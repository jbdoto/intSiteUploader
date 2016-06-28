-- drop tables with FK constraints
DROP TABLE IF EXISTS pcrbreakpoints;
DROP TABLE IF EXISTS sites;
DROP TABLE IF EXISTS multihitlengths;
DROP TABLE IF EXISTS multihitpositions;
DROP TABLE IF EXISTS samples;

CREATE TABLE samples (
    sampleID int NOT NULL,
    sampleName varchar(255) NOT NULL,
    refGenome varchar(10) NOT NULL,
    gender char(1) NOT NULL,
    miseqid varchar(255),
    PRIMARY KEY (sampleID),
    CONSTRAINT uniq_samples UNIQUE (sampleName, refGenome),
    CONSTRAINT non_empty_ref CHECK (LENGTH(refGenome) >= 1),
    CONSTRAINT non_empty_sample CHECK (LENGTH(sampleName) >= 1),
    CONSTRAINT f_for_female_m_for_male CHECK(gender IN ('f', 'm'))
);

-- unique hits
CREATE TABLE sites (
    siteID int NOT NULL,
    sampleID int NOT NULL,
    position int NOT NULL,
    chr varchar(255) NOT NULL,
    strand char(1) NOT NULL,
    PRIMARY KEY (siteID),
    FOREIGN KEY (sampleID) REFERENCES samples(sampleID),
    CONSTRAINT non_empty_chr CHECK (LENGTH(chr) >= 1)
);

CREATE TABLE pcrbreakpoints (
    siteID int NOT NULL,
    breakpoint int NOT NULL,
    count int NOT NULL,
    PRIMARY KEY (siteID, breakpoint),
    FOREIGN KEY (siteID) REFERENCES sites(siteID) 
);

-- multihit schema 
CREATE TABLE multihitpositions (
    multihitID int NOT NULL,
    sampleID int NOT NULL,
    position int NOT NULL,
    chr varchar(255) NOT NULL,
    strand char(1) NOT NULL,
    PRIMARY KEY (multihitID, position, chr, strand),
    FOREIGN KEY (sampleID) REFERENCES samples(sampleID),
    CONSTRAINT non_empty_chr CHECK (LENGTH(chr) >= 1)
);

CREATE TABLE multihitlengths (
    multihitID int NOT NULL,
    length int NOT NULL,
    count int NOT NULL,
    PRIMARY KEY (multihitID, length),
    FOREIGN KEY (multihitID) REFERENCES multihitpositions(multihitID) 
);
