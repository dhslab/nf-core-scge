#!/usr/bin/env bash

PWD=$(realpath .)

sed -e 's!REPLACEPATH!'$PWD'!g' demux_test_mastersheet.csv.orig > demux_test_mastersheet.csv
sed -e 's!REPLACEPATH!'$PWD'!g' analysis_test_mastersheet.csv.orig > analysis_test_mastersheet.csv
sed -e 's!REPLACEPATH!'$PWD'!g' cram_test_mastersheet.csv.orig > cram_test_mastersheet.csv
sed -e 's!REPLACEPATH!'$PWD'!g' bam_test_mastersheet.csv.orig > bam_test_mastersheet.csv
sed -e 's!REPLACEPATH!'$PWD'!g' fastq_test_mastersheet.csv.orig > fastq_test_mastersheet.csv
sed -e 's!REPLACEPATH!'$PWD'!g' gatherfastq_test_mastersheet.csv.orig > gatherfastq_test_mastersheet.csv
sed -e 's!REPLACEPATH!'$PWD'!g' demux_fastq/fastq_list.csv.orig > demux_fastq/fastq_list.csv
sed -e 's!REPLACEPATH!'$PWD'!g' demux_fastq/Reports/fastq_list.csv.orig > demux_fastq/Reports/fastq_list.csv

