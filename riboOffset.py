import sys
import argparse
import re
from collections import defaultdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bedfile", help="ribosome footprint", type=str)
    parser.add_argument("offsetparam", help="offset parameters", type=str)
    parser.add_argument("--pos15", help="map A site (15nt shift)", dest='pos15', action='store_true')
    args = parser.parse_args()

    offsets = {}
    with open(args.offsetparam, 'rt') as offsetfile:
        for l in offsetfile:
            l = re.sub('\s+$', '', l)
            a = l.split('\t')
            offsets[a[0]] = int(a[1])
        
    shift=0
    if (args.pos15):
        shift=3

    #NM_001190452_951_1026
    cnt = 0
    outbed = re.sub('.bed$','.offset.bed', args.bedfile)
    with open(args.bedfile, 'rt') as infile, open(outbed, 'wt') as outfile:
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
            if (not (str(length) in offsets)):
                continue

            dist_from_start = a[1] - st
            dist_from_end = a[1] - en

            offset = offsets[str(length)]+shift

            index = int ((a[1]+offset - st)/3)
            outfile.write('%s\t%d\t%d\t%s\t%d\t+\n' % (name,a[1]+offset,a[1]+offset+length,a[3],index))
            if (cnt%100000==0):
                sys.stdout.write('.')
                sys.stdout.flush()
            cnt+=1

    infile.close()
