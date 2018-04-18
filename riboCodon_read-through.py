import sys
import argparse
import re
from Bio import SeqIO
from collections import defaultdict

gencode = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}


def load_sequences (fastafile):
    fasta_sequences = SeqIO.parse(open(fastafile),'fasta')
    hSeq = {}
    hMap = defaultdict(list)
    hRT = defaultdict(list)
    for fasta in fasta_sequences:
        name, seq = fasta.id, str(fasta.seq)
        b = name.split('_')
        if (b[2] == '' or b[0].startswith('NR_')):
            continue
        ref = b[0] + "_" + b[1]
        st = int(b[2])
        en = int(b[3])
        if en<st or en==len(seq):
            continue
        for k in range(0, en - st, 3):
            #print seq[st + k:st + k + 3], gencode[seq[st + k:st + k + 3]]
            hMap[name].append(0)

        ####read-through ORF
        #print name
        #print seq
        k = en-st
        while(-1):
            if st+k+3>=len(seq):
                break
            c = seq[st + k:st + k + 3]
            a = gencode[c]
            #print c, a
            if (a == '_'):
                break
            else:
                hRT[name].append(0)

            k = k+3

        hSeq[name] = seq

    return hSeq, hMap, hRT

def load_refToSym (reffile):
    with open(reffile, 'rt') as infile:
        hRef = {}
        for l in infile:
            l = re.sub('\s+$', '', l)
            a = l.split('\t')
            hRef[a[1]] = a[0]

    return hRef

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fastafile", help="mRNA sequence with ID showing start and end of CDS with _", type=str)
    parser.add_argument("refToSym", help="file mapping gene symbols to RefSeq IDs", type=str)
    parser.add_argument("bedfile", help="ribosome footprint", type=str)
    parser.add_argument("--pos15", help="map A site (15nt shift)", dest='pos15', action='store_true')
    args = parser.parse_args()

    shift=0
    if (args.pos15):
        shift=3

    hSeq, hMap, hRT = load_sequences(args.fastafile)
    hRef = load_refToSym(args.refToSym)

    #NM_001190452_951_1026
    cnt = 0
    hRTcnt = {}
    hMcnt = {}
    with open(args.bedfile, 'rt') as infile:
        for l in infile:
            l = re.sub('\s+$', '', l)
            a = l.split('\t')
            name = a[0]
            b = name.split('_')
            if (b[2] == '' or b[0].startswith('NR_')):
                continue
            ref = b[0] + "_" + b[1]
            st = int(b[2])
            en = int(b[3])
            a[1] = int(a[1])
            a[2] = int(a[2])
            length = a[2]-a[1]
            if (length>33 or length<29):
                continue

            dist_from_start = a[1] - st
            dist_from_end = a[1] - en

            offset = 12+shift
            if (length>32):
                offset = 14+shift
            elif (length>30):
                offset = 13+shift

            index = int ((a[1]+offset - st)/3)
            cdslen = int((en-st)/3)

            if (index<0 or index>=len(hMap[name])+len(hRT[name])):
                continue
            if index<len(hMap[name]):
                hMap[name][index] += 1
                if name in hMcnt:
                    hMcnt[name] += 1
                else:
                    hMcnt[name] = 1
            elif index<len(hMap[name])+len(hRT[name]):
                if name in hRTcnt:
                    hRTcnt[name] += 1
                else:
                    hRTcnt[name] = 1

            if (cnt%100000==0):
                sys.stdout.write('.')
                sys.stdout.flush()
            cnt+=1

    infile.close()

    if (args.pos15):
        outfile4 = re.sub('.bed$','.rt-count.pos15.txt', args.bedfile)
    else:
        outfile4 = re.sub('.bed$','.rt-count.pos12.txt', args.bedfile)
        
    #with open(outfile1, 'wt') as out1, open(outfile2, 'wt') as out2, open(outfile3, 'wt') as out3, open(outfile4, 'wt') as out4:
    with open(outfile4, 'wt') as out4:
        for name in hMcnt:
            b = name.split('_')
            ref = b[0] + "_" + b[1]
            st = int(b[2])
            en = int(b[3])
            gene = hRef[ref]
            if name in hRTcnt:
                #print name, hMap[name]
                out4.write(name+"\t"+str(hMcnt[name])+"\t"+str(hRTcnt[name])+"\n")
            else:
                out4.write(name + "\t" + str(hMcnt[name]) + "\t0\n")

    out4.close()
