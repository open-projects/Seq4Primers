#!/usr/bin/env python3

import os
import re
import sys
import argparse
import sqlite3
import csv


chr_sql = "CREATE TABLE `hgChrom` ( \
  `chrom` text, \
  `seq` text, \
  PRIMARY KEY (`chrom`) \
)"

loc_sql_1 = "CREATE TABLE `SNP_location` ( \
  `snp_id` text, \
  `chrom` text, \
  `pos_from` integer, \
  `len` integer, \
  `reference` text, \
  `is_common` integer \
)"
loc_sql_2 = "CREATE INDEX snp_idx ON SNP_location (snp_id)"
loc_sql_3 = "CREATE INDEX chrom_idx ON SNP_location (chrom)"
loc_sql_4 = "UPDATE SNP_location SET is_common = NULL WHERE is_common = '\\N'"

all_sql_1 = "CREATE TABLE `SNP_allele` ( \
  `snp_id` text, \
  `allele` text, \
  `freq` real DEFAULT NULL, \
  PRIMARY KEY (`snp_id`,`allele`) \
)"
all_sql_2 = "UPDATE SNP_allele SET freq = NULL WHERE freq = '\\N'"

def main():
    input_parser = argparse.ArgumentParser(description='MakeDB: the program for creating of SQLite data file for Seq4Primers program.')
    input_parser.add_argument('-s', metavar='chrom.csv', help='the file with chromosome sequences', required=True)
    input_parser.add_argument('-l', metavar='location.csv', help='the file with SNP coordinates', required=True)
    input_parser.add_argument('-a', metavar='allele.csv', help='the file with SNP alleles', required=True)
    input_parser.add_argument('-o', metavar='data.sqlite', default='data.sqlite', help='the output file', required=False)

    args = input_parser.parse_args()
    chr_file = args.s
    loc_file = args.l
    all_file = args.a
    out_file = args.o

    conn = None
    try:
        if os.path.isfile(out_file):
            os.remove(out_file)
        conn = sqlite3.connect(out_file)
    except Error as e:
        print(e)
        exit()

    cursor = conn.cursor()
    cursor.execute(chr_sql)
    cursor.execute(loc_sql_1)
    cursor.execute(all_sql_1)
    conn.commit()

    csv.field_size_limit(sys.maxsize)
    with open(chr_file, 'r') as fchr:
        csv_reader = csv.reader(fchr, delimiter='\t')
        for row in csv_reader:
            cursor.execute("INSERT INTO hgChrom (chrom, seq) VALUES (?, ?);", row)
        conn.commit()

    with open(loc_file, 'r') as floc:
        csv_reader = csv.reader(floc, delimiter='\t')
        for row in csv_reader:
            cursor.execute("INSERT INTO SNP_location (snp_id, chrom, pos_from, len, reference, is_common) VALUES (?, ?, ?, ?, ?, ?);", row)
        conn.commit()

    with open(all_file, 'r') as floc:
        csv_reader = csv.reader(floc, delimiter='\t')
        for row in csv_reader:
            cursor.execute("INSERT INTO SNP_allele (snp_id, allele, freq) VALUES (?, ?, ?);", row)
        conn.commit()

    cursor.execute(loc_sql_2)
    cursor.execute(loc_sql_3)
    cursor.execute(loc_sql_4)
    cursor.execute(all_sql_2)
    conn.commit()

    print("...done")
# end of main()

if __name__ == '__main__':
    main()
# end of script
