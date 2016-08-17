#different functions that get repeatedly used in dealing with .psl files from blat
import sys, subprocess, os
from Bio.Seq import Seq
from Bio import SeqIO
import parse_fa_chr

#reads in a PSL file and provides block sizes, query start positions, target start positions, mapping between query and target coordinates, and strand of hits
def readPSLs(psl_file,flip):
    f = open(psl_file,'r')
    num_blocks = 0
    for line in f:
        list = line.split('\t')
        chr_size = int(list[-7])
        gen_strand = list[8][1]
        block_sizes = map(int, list[-3].strip().strip(',').split(','))
        q_starts = map(int, list[-2].strip().strip(',').split(','))
        t_starts = map(int, list[-1].strip().strip(',').split(','))
        # fix coordinates depending on strand information
        if (gen_strand == '-') & (flip == 'y'):
            for i in range(len(t_starts)): t_starts[i] = chr_size - t_starts[i]
        qCoord2tCoord = {}
        for i in range(len(q_starts)):
            qstart = q_starts[i]
            tstart = t_starts[i]
            qCoord2tCoord[qstart] = tstart
            qCoord2tCoord[tstart] = qstart
    f.close()
    return block_sizes, q_starts, t_starts, qCoord2tCoord, gen_strand

# returns a list of a dictionaries of all relevant information in each line from PSL file. Can sort hits by positions using min_max_sort_option, can return info on only best match
def readPSLmultiLine(psl_file, min_max_sort_option,best_match_only,flip):
    psl_lines = []
    f = open(psl_file,'r')
    num_blocks = 0
    for line in f:
        list = line.strip().split('\t')
        chr_size = int(list[-7])
        gen_strand = list[8]
        block_sizes = map(int, list[-3].strip().strip(',').split(','))
        q_starts = map(int, list[-2].strip().strip(',').split(','))
        t_starts = map(int, list[-1].strip().strip(',').split(','))
        q_aln_bases = int(list[0])
        # fix coordinates depending on strand information
        if (gen_strand == '-') & (flip == 'y'):
            for i in range(len(t_starts)): t_starts[i] = chr_size - t_starts[i]

                #sort lists if option given
        if min_max_sort_option == 'y':
            q_starts = sorted(q_starts)
            t_starts = sorted(t_starts)
            block_sizes = [block_size for t_start,block_size in sorted(zip(t_starts,block_sizes))]
        qCoord2tCoord = {}
        for i in range(len(q_starts)):
            qstart = q_starts[i]
            tstart = t_starts[i]
            qCoord2tCoord[qstart] = tstart
            qCoord2tCoord[tstart] = qstart
        psl_lines.append({'block_sizes':block_sizes,'q_starts':q_starts,'t_starts':t_starts,'qCoord2tCoord':qCoord2tCoord,'gen_strand':gen_strand,'q_aln_bases':q_aln_bases})
    f.close()

    # if only want best match 
    if best_match_only == 'y':
        best_match_bases = 0
        best_match_line = {}
        for line in psl_lines:
            if line['q_aln_bases'] > best_match_bases:
                best_match_line = line
                best_match_bases = line['q_aln_bases']
        return best_match_line
    
    return psl_lines
