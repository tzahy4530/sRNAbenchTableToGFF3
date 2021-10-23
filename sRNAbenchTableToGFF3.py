#!/usr/bin/python

import sys
import pandas as pd


# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

def run(input, output, additional=None):
    """
    :param additional: additonal sRNAbench output prediction file.
    :param input: sRNAbench output prediction files 'novel.txt', 'novel454.txt'
    :param output: output path of the GFF formatted file.
    :return:
    """
    version = "##gff-version 3\n"
    gff3_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gff3 = pd.DataFrame(columns=gff3_columns)
    table = pd.read_csv(input, sep='\t')

    if additional:
        table_to_add = pd.read_csv(additional, sep='\t')
        table = table.append(table_to_add)

    for index, row in table.iterrows():
        name = row['name']
        seqId = row['seqName']
        name5p = row['5pname']
        seq5p = row['5pseq']
        name3p = row['3pname']
        seq3p = row['3pseq']
        strand = row['strand']
        hairpin = row['hairpinSeq']
        start = row['start']
        end = row['end']

        gff_row = [[seqId, '.', 'pre_miRNA', start, end, '.', strand, '.', f'ID={name}']]

        if strand == '+':
            try:
                offset5p = len(hairpin.split(seq5p)[0])
                start5p = start + offset5p
                end5p = start + offset5p + len(seq5p) - 1
                gff_row.append([seqId, '.', 'miRNA', start5p, end5p, '.', strand, '.', f'ID={name5p}'])
            except:
                pass

            try:
                offset3p = len(hairpin.split(seq3p)[0])
                start3p = start + offset3p
                end3p = start + offset3p + len(seq3p) - 1
                gff_row.append([seqId, '.', 'miRNA', start3p, end3p, '.', strand, '.', f'ID={name3p}'])
            except:
                pass

        else:
            try:
                offset5p = len(hairpin.split(seq5p)[0])
                end5p = end - offset5p
                start5p = end - offset5p - len(seq5p) + 1
                gff_row.append([seqId, '.', 'miRNA', start5p, end5p, '.', strand, '.', f'ID={name5p}'])
            except:
                pass

            try:
                offset3p = len(hairpin.split(seq3p)[0])
                end3p = end - offset3p
                start3p = end - offset3p - len(seq3p) + 1
                gff_row.append([seqId, '.', 'miRNA', start3p, end3p, '.', strand, '.', f'ID={name3p}'])
            except:
                pass

        miRNAs = pd.DataFrame(gff_row, columns=gff3_columns)

        gff3 = gff3.append(miRNAs)

    with open(output, 'w') as file:
        file.write(version)

    gff3.to_csv(output, index=False, header=False, mode="a", sep='\t')


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    input = None
    output = None
    add = None
    args = []
    for i in range(1, len(sys.argv), 2):
        arg = sys.argv[i]
        if arg == '-i':
            input = sys.argv[i + 1]
        if arg == '-o':
            output = sys.argv[i + 1]
        if arg == '-a':
            add = sys.argv[i + 1]
        if arg == '--help' or arg == '-h':
            print(f'Manual:\n'
                  f' -i : sRNAbench prediction output, like novel.txt/novel451.txt.\n'
                  f' -o : output path.\n'
                  f' -a : additional input file.\n')
            sys.exit()

    if not input:
        raise ('Input path is required (-i <path>)')
    if not output:
        raise ('Output path is required (-o <path>)')
    run(input, output, add)
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
