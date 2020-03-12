#!/usr/bin/env python3

import re
import argparse
import sqlite3
import datetime

select_seq = " \
    SELECT \
    chrom, pos_from, pos_from + len - 1, \
        LOWER(SUBSTR((SELECT seq FROM hgChrom AS t2 WHERE t2.chrom = t1.chrom), CASE WHEN pos_from - ? > 0 THEN pos_from - ? ELSE 1 END, ?)) || \
        LOWER(SUBSTR((SELECT seq FROM hgChrom AS t2 WHERE t2.chrom = t1.chrom), pos_from, len)) || \
        LOWER(SUBSTR((SELECT seq FROM hgChrom AS t2 WHERE t2.chrom = t1.chrom), pos_from + len, ?)) \
    FROM SNP_location AS t1 WHERE snp_id = ? \
"

select_snp = " \
    SELECT snp_id, pos_from, pos_from + len - 1, reference, GROUP_CONCAT(allele), GROUP_CONCAT(freq) \
    FROM SNP_location INNER JOIN SNP_allele USING(snp_id) WHERE \
    chrom = ? \
    AND \
    pos_from BETWEEN ? AND ? \
    AND \
    pos_from + len - 1 BETWEEN ? AND ? \
    GROUP BY snp_id \
    ORDER BY pos_from \
"

def main():
    input_parser = argparse.ArgumentParser(description='seq4primers: the program to generate GenBank files with subsequences which flank the SNPs of interests.')
    input_parser.add_argument('-s', metavar='SNP_ID', nargs='+', help='set of SNPs', required=True)
    input_parser.add_argument('-f', metavar='flank_size', type=int, default=500, help='flank size', required=False)
    input_parser.add_argument('-db', metavar='SQLite_database_file', default='db/data.sqlite', help='SQLite database file', required=False)
    input_parser.add_argument('-strict', action='store_true', help='use GenBank strict rules', required=False)

    args = input_parser.parse_args()
    snp_set = args.s
    flank_size = args.f
    db_file = args.db
    strict_rules = args.strict

    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)
        exit()

    now = datetime.datetime.now()
    with open("seq4primers.txt", 'w') as seq_file:
        cursor1 = conn.cursor()
        for snp_id in snp_set:
            with open("{}.gb".format(snp_id), 'w') as gb_file:
                header0 = 'LOCUS  ';
                header1 = '{:>24}'.format('Locus')
                header2 = ' MiHA ';
                header3 = '000000';
                header4 = ' bp    DNA     linear   UNA {}-{}-{}\n'.format(now.strftime("%d"), now.strftime("%b"), now.year)
                header5 = 'DEFINITION           \n'
                header6 = 'FEATURES             Location/Qualifiers\n'

                cursor1.execute(select_seq, [flank_size, flank_size, flank_size, flank_size, snp_id])
                for chrom, beg, end, seq in cursor1.fetchall():
                    header1 = '{:>24}'.format(chrom)
                    start_shift = beg - flank_size - 1
                    if start_shift < 0:
                        start_shift = 0
                        header5 = 'DEFINITION           {} {}:{}..{}\n'.format(snp_id, chrom, 1, len(seq))
                        print("WARNING: the left flank was truncated!")
                    else:
                        header5 = 'DEFINITION           {} {}:{}..{}\n'.format(snp_id, chrom, beg - flank_size, beg - flank_size + len(seq) - 1)
                    seq = seq.lower()
                    header3 = '{: > 9}'.format(len(seq))
                    gb_file.write(header0 + header1 + header2 + header3 + header4 + header5 + header6)

                    current_pos = 0
                    flanked_seq = ''

                    cursor2 = conn.cursor()
                    cursor2.execute(select_snp, [chrom, beg - flank_size, end + flank_size, beg - flank_size, end + flank_size])
                    for curr_snp_id, snp_beg, snp_end, ref, obs, freq in cursor2.fetchall():
                        if ref == '-':
                            snp_end, snp_beg = snp_beg, snp_end
                        if not freq:
                            print("NO FREQ: {}".format(curr_snp_id))
                            continue
                        if not ref:
                            exit("ERROR - no SNP record: {}\n".format(curr_snp_id))
                        if not ref:
                            exit("ERROR - wrong SNP allele format: {} => {}\n".format(curr_snp_id, obs))

                        ref = ref.lower();
                        snp_beg = snp_beg - start_shift
                        snp_end = snp_end - start_shift
                        snp_seq = '-'
                        if snp_end - snp_beg + 1:
                            snp_seq = seq[snp_beg - 1 : snp_end]
                        if snp_seq.upper() != ref.upper():
                            exit("ERROR - WRONG SNP POSITION: {} => {} <-> {}\n".format(curr_snp_id, ref, snp_seq))
                        if not check_allele(ref, obs):
                            exit("ERROR - no reference allele: {} => {} \n".format(curr_snp_id, obs))

                        code = obs + '/' + freq
                        if snp_id == curr_snp_id:
                            code = '[' + curr_snp_id + ':' + code + ']'
                        else:
                            code = '(' + curr_snp_id + ':' + code + ')'

                        current_len = snp_beg - current_pos - 1
                        flanked_seq = flanked_seq + seq[current_pos : current_pos + current_len]
                        flanked_seq = flanked_seq + code
                        current_pos = snp_end

                        obs_set = re.split(r', *', obs)
                        freq_set = re.split(r', *', freq)
                        alleles = list()
                        for i in range(0, len(obs_set)):
                            if obs_set[i] == '-':
                                alleles.append('del')
                            else:
                                alleles.append(obs_set[i] + '-' + freq_set[i])
                        seq = seq[0 : snp_beg - 1] + seq[snp_beg - 1 : snp_end].upper() + seq[snp_end:]

                        if strict_rules:
                            gb_file.write('    variation        {}..{}\n'.format(snp_beg, snp_end))
                        else:
                            gb_file.write('    Polymorphism     {}..{}\n'.format(snp_beg, snp_end))
                        if curr_snp_id == snp_id:
                            gb_file.write('                     /MiHA="{}"\n'.format(curr_snp_id))
                        gb_file.write('                     /label="{}"\n'.format(curr_snp_id))
                        gb_file.write('                     /frequency="{}"\n'.format('/'.join(alleles)))
                    flanked_seq = flanked_seq + seq[current_pos:]
                    seq_file.write(snp_id + "\t" + flanked_seq + "\n")
                    gb_file.write("ORIGIN\n")

                    i = 0
                    while i < len(seq):
                        gb_file.write('{:>9}'.format(i + 1))
                        sstr = seq[i : i + 60]
                        sstr = re.sub(r'(.{1,10})', r' \1', sstr)
                        gb_file.write(sstr + "\n")
                        i += 60
                    gb_file.write("//\n")
            print("the file for {} is created".format(snp_id))
        print("...done")
# end of main()

def check_allele(reference, observed):
    alleles = re.split(r', *', observed)
    for allele in alleles:
        if allele.upper() == reference.upper():
            return True
    return False
# end of check_allele()

if __name__ == '__main__':
    main()
# end of script
