mysql -hmicrob98.med.upenn.edu -uxxxxxx -p -e '
USE intsites_miseq_dev;
SET FOREIGN_KEY_CHECKS = 0;
TRUNCATE samples;
TRUNCATE sites;
TRUNCATE pcrbreakpoints;
TRUNCATE multihitpositions;
TRUNCATE multihitlengths;
SET FOREIGN_KEY_CHECKS = 1;
EXIT
'
