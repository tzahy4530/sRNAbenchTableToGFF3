#!/usr/bin/python

import sys
import pandas as pd
ids_dic = {}

def handleGivenName(name, df, column):
    """
    This Function handle the given name of miRNA,
    name of miRNA used as a key, most be unique.
    :param name: sRNAbench given name.
    :param df: dataframe which including all the records.
    :param column: column name.
    :return: unique name.
    """
    if len(df[df[column] == name]) > 1:
        if name not in ids_dic:
            ids_dic[name] = 0
        ids_dic[name] += 1
        name = f'{name}_{ids_dic[name]}'

    return name

def run(input, output, additional=None, fasta_path=None):
    """
    This Function will create GFF3 file from the sRNAbench output.
    :param fasta_path: a path to create fasta file from the gff3 table.
    :param additional: additonal sRNAbench output prediction file.
    :param input: sRNAbench output prediction files 'novel.txt', 'novel454.txt'
    :param output: output path of the GFF formatted file.
    :return:
    """
    version = "##gff-version 3\n"
    gff3_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gff3 = pd.DataFrame(columns=gff3_columns)
    table = pd.read_csv(input, sep='\t')
    if fasta_path is not None:
        fasta_file = ''
        open(fasta_path, 'w').close()


    if additional:
        table_to_add = pd.read_csv(additional, sep='\t')
        table = table.append(table_to_add)

    for index, row in table.iterrows():
        name = handleGivenName(row['name'], table, 'name')
        seqId = row['seqName']
        name5p = handleGivenName(row['5pname'], table, '5pname')
        seq5p = row['5pseq']
        name3p = handleGivenName(row['3pname'], table, '3pname')
        seq3p = row['3pseq']
        strand = row['strand']
        hairpin = row['hairpinSeq']
        start = row['start']
        end = row['end']

        if row['5pRC'] >= row['3pRC']:
            name5p += '-m'
            name3p += '-s'
        else:
            name5p += '-s'
            name3p += '-m'

        seq5p_freq = len(table[(table['5pseq'] == seq5p) | (table['3pseq'] == seq5p)])
        seq3p_freq = len(table[(table['5pseq'] == seq3p) | (table['3pseq'] == seq3p)])

        name5p += f'-{seq5p_freq}'
        name3p += f'-{seq3p_freq}'
        
        if fasta_path is not None:
            if not pd.isnull(seq5p):
                fasta_file += f'>{name5p}\n{seq5p}\n'
            if not pd.isnull(seq3p):
                fasta_file += f'>{name3p}\n{seq3p}\n'

            if len(fasta_file) > 100000:
                with open(fasta_path, 'a+') as f:
                    f.write(fasta_file)
                fasta_file = ''

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

    if fasta_path is not None:
        with open(fasta_path, 'a+') as f:
            f.write(fasta_file)

    gff3.to_csv(output, index=False, header=False, mode="a", sep='\t')


if __name__ == '__main__':
    input = None
    output = None
    add = None
    fasta_path = None
    args = []
    for i in range(1, len(sys.argv), 2):
        arg = sys.argv[i]
        if arg == '-i':
            input = sys.argv[i + 1]
        if arg == '-o':
            output = sys.argv[i + 1]
        if arg == '-a':
            add = sys.argv[i + 1]
        if arg == '--create-fasta':
            fasta_path = sys.argv[i + 1]
        if arg == '--help' or arg == '-h':
            print(f'Manual:\n'
                  f' -i <path>: sRNAbench prediction output, like novel.txt/novel451.txt.\n'
                  f' -a <path>: additional input file.\n'
                  f' -o <path>: output path.\n'
                  f' --create-fasta <path>: create fasta file from the gff3 table.\n'
                  )

            sys.exit()

    if not input:
        raise ('Input path is required (-i <path>)')
    if not output:
        raise ('Output path is required (-o <path>)')
    run(input, output, add, fasta_path)
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
