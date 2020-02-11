#!/usr/bin/python
#Author: liuhui
#EMail: liuhui@bjfu.edu.cn / hui.liu@umu.se
#Description: generate 0 and 4fold sites

import sys
import math
from collections import OrderedDict


USAGE = "\nusage: python %s <gff3> <fasta> <0-fold output> <4-fold output>\n" % sys.argv[0]

if len(sys.argv) != 5:
    print USAGE
    sys.exit()


standardCodonTable = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

#https://github.com/russcd/vcf2MK/blob/master/src/fold_counts.h
def get_degeneracy(codon_table):
    nucs = ['A', 'T', 'G', 'C']
    codon_degen = OrderedDict()
    for codon in codon_table:
        degen_list = []
        for i in range(len(codon)):
            degen = 0
            for j in nucs:
                if i == 0:
                    mod_codon = j + codon[1:]
                elif i == 1:
                    mod_codon = codon[0] + j + codon[2]
                elif i == 2:
                    mod_codon = codon[0:2] + j

                if codon_table[codon] == codon_table[mod_codon]:
                    degen += 1

            degen_list.append(degen)

            for site in range(len(degen_list)):
                if degen_list[site] == 1:
                    degen_list[site] = 0
        codon_degen[codon] = degen_list

    return codon_degen

# https://www.slideserve.com/adelle/ejercicios-de-alineamiento-de-secuencias-clustalw-insertar-secuencias-de-fasta
# https://ppt-online.org/57106
codonDegenDict = {"AAA" : [0, 0, 2], "AAC" : [0, 0, 2], "AAG" : [0, 0, 2], "AAT" : [0, 0, 2],
                  "ACA" : [0, 0, 4], "ACC" : [0, 0, 4], "ACG" : [0, 0, 4], "ACT" : [0, 0, 4],
                  "AGA" : [2, 0, 2], "AGC" : [0, 0, 2], "AGG" : [2, 0, 2], "AGT" : [0, 0, 2],
                  "ATA" : [0, 0, 3], "ATC" : [0, 0, 3], "ATG" : [0, 0, 0], "ATT" : [0, 0, 3],
                  "CAA" : [0, 0, 2], "CAC" : [0, 0, 2], "CAG" : [0, 0, 2], "CAT" : [0, 0, 2],
                  "CCA" : [0, 0, 4], "CCC" : [0, 0, 4], "CCG" : [0, 0, 4], "CCT" : [0, 0, 4],
                  "CGA" : [2, 0, 4], "CGC" : [0, 0, 4], "CGG" : [2, 0, 4], "CGT" : [0, 0, 4],
                  "CTA" : [2, 0, 4], "CTC" : [0, 0, 4], "CTG" : [2, 0, 4], "CTT" : [0, 0, 4],
                  "GAA" : [0, 0, 2], "GAC" : [0, 0, 2], "GAG" : [0, 0, 2], "GAT" : [0, 0, 2],
                  "GCA" : [0, 0, 4], "GCC" : [0, 0, 4], "GCG" : [0, 0, 4], "GCT" : [0, 0, 4],
                  "GGA" : [0, 0, 4], "GGC" : [0, 0, 4], "GGG" : [0, 0, 4], "GGT" : [0, 0, 4],
                  "GTA" : [0, 0, 4], "GTC" : [0, 0, 4], "GTG" : [0, 0, 4], "GTT" : [0, 0, 4],
                  "TAA" : [0, 2, 2], "TAC" : [0, 0, 2], "TAG" : [0, 0, 2], "TAT" : [0, 0, 2],
                  "TCA" : [0, 0, 4], "TCC" : [0, 0, 4], "TCG" : [0, 0, 4], "TCT" : [0, 0, 4],
                  "TGA" : [0, 2, 0], "TGC" : [0, 0, 2], "TGG" : [0, 0, 0], "TGT" : [0, 0, 2],
                  "TTA" : [2, 0, 2], "TTC" : [0, 0, 2], "TTG" : [2, 0, 2], "TTT" : [0, 0, 2]}

def gff2dict(filename):
    """
    The phase of a CDS feature depends on the associated upstream CDS feature.
    If there is the length/3 of the previous CDS feature leaves a remainder of 1,
    your CDS feature requires a phase of 2 (the last base of the previous and the
    first two bases of this form a codon triplett). Standalone features and the
    first features always have a 0 phase value.
    """
    ids = {}
    cds_info = OrderedDict()
    with open(filename, 'r') as f:
        for line in f:
            if line[0] == '#': continue
            lsp = line.rstrip().split("\t")
            chr, source, type, start, end, score, strand, phase, attibutes = lsp
            d_attr = dict([v.split('=') for v in attibutes.strip(';').split(';')])
            if type == 'mRNA':
                id = d_attr['ID']
                ids.setdefault(chr, []).append([id, strand])
            elif type == 'CDS':
                id = d_attr['Parent']
                cds_info.setdefault(id, []).append([chr, int(start), int(end), strand])
    return ids, cds_info

def parseFasta(filename):
    fas = OrderedDict()
    id = None
    with open(filename, 'r') as fh:
        for line in fh:
            if line[0] == '>':
                header = line[1:].rstrip()
                id = header.split()[0]
                fas[id] = []
            else:
                fas[id].append(line.rstrip().upper())
        for id, seq in fas.items():
            fas[id] = ''.join(seq)
    return fas


def reverse_comp(sequence):
    comp_dict = {
        'A': 'T',
        'B': 'V',
        'C': 'G',
        'D': 'H',
        'G': 'C',
        'H': 'D',
        'M': 'K',
        'N': 'N',
        'R': 'Y',
        'S': 'S',
        'T': 'A',
        'U': 'A',
        'V': 'B',
        'W': 'W',
        'X': 'X',
        'Y': 'R'}
    #comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', '-': '-', 'N': 'N'}
    sequence = sequence.upper()
    sequence_rev = ''
    for i in range(1, len(sequence)+1):
        sequence_rev += comp_dict[sequence[-i]]
    return sequence_rev


def getCDS(cds_records, reference):
    seq = ''
    for record in cds_records:
        chr, start, end, strand = record
        seq += reference[(start-1):end]
    if strand == "-":
        seq = reverse_comp(seq)
        if len(seq) % 3 == 1:
            seq =seq[:-1]
        elif len(seq) % 3 == 2:
            seq =seq[:-2]
    return seq

def CDSsum(lst):
    # [['Chr15', 130096, 130344, '+'], ['Chr15', 131129, 131242, '+']]
    l = 0
    for i in lst:
        start, end = i[1:3]
        l += end - start + 1
    return l

def getcodon(cds_seq, n_bases, strand):
    number = int(math.ceil(n_bases/float(3)))
    if strand == "+":
        codon = cds_seq[(number-1)*3: number*3]
    else:
        cds_seq = cds_seq[::-1]
        codon = cds_seq[(number-1)*3: number*3]
        codon = codon[::-1]

    m = n_bases % 3
    if m == 0:
        codon_pos = 2
    elif m == 1:
        codon_pos = 0
    else:
        codon_pos = 1
    if strand == "-":
        codon_pos = 2 - codon_pos
    return codon, codon_pos

# input
#codonDegenDict = get_degeneracy(standardCodonTable)
# https://github.com/parkingvarsson/Degeneracy/blob/master/degeneracy.py
"""
codonDegenDict = {"ACA" : [0, 0, 4], "ACC" : [0, 0, 4], "ACG" : [0, 0, 4], "ACT" : [0, 0, 4],
                  "CCA" : [0, 0, 4], "CCC" : [0, 0, 4], "CCG" : [0, 0, 4], "CCT" : [0, 0, 4],
                  "CGA" : [2, 0, 4], "CGC" : [5, 0, 4], "CGG" : [2, 0, 4], "CGT" : [5, 0, 4],
                  "CTA" : [2, 0, 4], "CTC" : [5, 0, 4], "CTG" : [2, 0, 4], "CTT" : [5, 0, 4],
                  "GCA" : [0, 0, 4], "GCC" : [0, 0, 4], "GCG" : [5, 5, 4], "GCT" : [0, 0, 4],
                  "GGA" : [0, 0, 4], "GGC" : [0, 0, 4], "GGG" : [0, 0, 4], "GGT" : [0, 0, 4],
                  "GTA" : [0, 0, 4], "GTC" : [0, 0, 4], "GTG" : [0, 0, 4], "GTT" : [0, 0, 4],
                  "TCA" : [5, 5, 4], "TCC" : [5, 5, 4], "TCG" : [5, 5, 4], "TCT" : [5, 5, 4],
                  "TGA" : [0, 5, 5], "TAG" : [0, 5, 5], "TAA" : [0 ,5, 5],
                  "AGG" : [5, 0, 5], "AGA" : [5, 0, 5], "TTG" : [5, 0, 5], "TTA" : [5, 0, 5],
                  "GAG" : [0, 0, 5], "GAA" : [0, 0, 5], "GAT" : [0, 0, 5], "GAC" : [0, 0, 5],
                  "AAG" : [0, 0, 5], "AAA" : [0, 0, 5], "AAC" : [0, 0, 5], "AAT" : [0, 0, 5],
                  "ATA" : [0, 0, 5], "ATC" : [0, 0, 5], "ATT" : [0, 0, 5], "CAG" : [0, 0, 5],
                  "CAA" : [0, 0, 5], "CAT" : [0, 0, 5], "CAC" : [0, 0, 5], "TGC" : [0, 0, 5],
                  "TGT" : [0, 0, 5], "TAC" : [0, 0, 5], "TAT" : [0, 0, 5], "TTC" : [0, 0, 5],
                  "TTT" : [0, 0, 5],
                  "ATG" : [0, 0, 0], "TGG" : [0, 0, 0]
}
"""

ids_dict, cds_info_dict = gff2dict(sys.argv[1])
genome = parseFasta(sys.argv[2])
chroms = [chr for chr in genome]

# output
zero_fold_out = open(sys.argv[3], 'w')
four_fold_out = open(sys.argv[4], 'w')

zero_fold = OrderedDict()
four_fold = OrderedDict()

# loop for each chomosome
for chr in chroms:
    ref_seq = genome[chr]
    # scan each mRNA in a chorosome
    if chr not in ids_dict: continue
    for item in ids_dict[chr]:
        id, strand = item
        if id not in cds_info_dict: continue
        full_cds = cds_info_dict[id]
        length = CDSsum(full_cds)
        mod = length % 3
        # extract the whole CDS region
        cds_seq = getCDS(full_cds, ref_seq)
        n = 0
        if mod != 0 and strand == "-":
            first_CDS = full_cds[0][2] - full_cds[0][1] + 1
            if mod == first_CDS:
                full_cds = full_cds[1:]
            # mod == 2 and first_CDS == 1
            elif mod > first_CDS:
                full_cds[0][1] = full_cds[0][1] + 1
        # scan each part of the CDS
        for r in full_cds:
            start, end = r[1:3]
            for pos in range(start, end+1):
                n += 1
                codon, codon_pos = getcodon(cds_seq, n, strand)
                if codon not in codonDegenDict: continue
                if codon in ['TAA', 'TGA', 'TAG']: continue
                i_fold = codonDegenDict[codon][codon_pos]
                r = chr + "_" + str(pos-1)
                #print id, pos, codon, codon_pos, i_fold, r, cds_seq
                if i_fold == 0:
                    zero_fold[r] = [chr, str(pos-1), str(pos), id]
                elif i_fold == 4:
                    four_fold[r] = [chr, str(pos-1), str(pos), id]

for i in zero_fold:
    if i in four_fold: continue
    zero_fold_out.write("\t".join(zero_fold[i]) + "\n")
zero_fold_out.close()

for i in four_fold:
    if i in zero_fold: continue
    four_fold_out.write("\t".join(four_fold[i]) + "\n")
four_fold_out.close()
