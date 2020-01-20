import twobitreader
import argparse
import pandas as pd
import numpy as np

# GENE_TSV_HEADER = ['name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'proteinId', 'alignID']
GENE_TSV_HEADER = ['chrom', 'startTranscription', 'endTranscription', 'name', 'unimportant', 'strand', 'startTranslation', 'endTranslation', 'unimportant2', 'exonCount', 'exonSize', 'exonStart']

REVERSE_COMPLEMENT_MAP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def getSeq(gene_dataframe, twobitreader):
    for i, gene in gene_dataframe.iterrows():
        seq = twobitreader[gene['chrom']][int(gene['startTranscription']):int(gene['endTranscription'])]
        seq = np.array(list(seq))
        if gene['strand'] == '-':
            for k, v in REVERSE_COMPLEMENT_MAP.items():
                seq[seq == k] = v
        exon_starts = np.array(gene['exonStart'].split(',')[:-1], dtype=np.int64)
        exon_sizes = np.array(gene['exonSize'].split(',')[:-1], dtype=np.int64)
        intron_exon_lbl = np.zeros(len(seq))
        for start, size in zip(exon_starts, exon_sizes):
            intron_exon_lbl[start: start + size] = 1
        yield seq, intron_exon_lbl


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('2bit_path', help='2bit file we use to get sequences)')
    parser.add_argument('gene_tsv_path', help='gene tsv path')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    genes_df = pd.read_csv(args.gene_tsv_path, sep='\t', names=GENE_TSV_HEADER)
    genome_reader = twobitreader.TwoBitFile('./hg19.2bit')

    gene_gen = getSeq(genes_df, genome_reader)
    for i, (seq, lbl) in enumerate(gene_gen):
        print(i)
        print(seq)
        print(lbl)
        break

    print("end")
    # print(genome['chr9'][35684800:35684869])