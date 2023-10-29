-- This SQL script contains the database commands needed to update the table. It ensures the data gets added correctly without any integrity issues.
-- Authors: Kosar Hooshmand

ALTER TABLE NCIG_ids ADD COLUMN Md5Sum TEXT;
ALTER TABLE NCIG_ids ADD COLUMN Flowcell_ID TEXT;
ALTER TABLE NCIG_ids ADD COLUMN Lane INTEGER;
ALTER TABLE NCIG_ids ADD COLUMN DateSequenced TEXT;
ALTER TABLE NCIG_ids ADD COLUMN Sample_ID TEXT;
ALTER TABLE NCIG_ids ADD COLUMN Species TEXT;
ALTER TABLE NCIG_ids ADD COLUMN Barcode TEXT;
ALTER TABLE NCIG_ids ADD COLUMN RunIdentifier TEXT;
ALTER TABLE NCIG_ids ADD COLUMN DateCollected TEXT;
ALTER TABLE NCIG_ids ADD COLUMN Technician TEXT;
ALTER TABLE NCIG_ids ADD COLUMN Batch TEXT;
ALTER TABLE NCIG_ids ADD COLUMN Run TEXT;
ALTER TABLE NCIG_ids ADD COLUMN LibraryStrategy TEXT;
ALTER TABLE NCIG_ids ADD COLUMN Platform TEXT;
ALTER TABLE NCIG_ids ADD COLUMN InstrumentModel TEXT;
ALTER TABLE NCIG_ids ADD COLUMN SoftwareVersion TEXT;
ALTER TABLE NCIG_ids ADD COLUMN Sex TEXT
