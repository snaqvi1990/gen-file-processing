#utilities for reading, handling, and writing MAF format multiple alignment files. Currently only deals with pairwise mafs

import sys,subprocess, os
import re

class pairMaf:
    def __init__(self,blocks_seqs={},seq_sizes=()):
        self.blocks = blocks_seqs # dict of tuples of blocks
        self.seqSizes = seq_sizes

    def subSeq(self,coord1,coord2,gap_strip='n',refidx=0):
        for block_pair in self.blocks:
            #print str(block_pair[refidx].start) + '\t' + str(block_pair[refidx].start+block_pair[refidx].length)
            if block_pair[refidx].contains(coord1):
                if block_pair[refidx].contains(coord1,coord2):
                    #print block_pair[refidx].start
                    return block_pair[refidx].subSeq(coord1,coord2,gap_strip)
                else:
                    #print block_pair[refidx].coords
                    #print block_pair[refidx].start
                    #print block_pair[refidx].start+block_pair[refidx].length
                    #print block_pair[refidx].seq
                    return block_pair[refidx].subSeq(coord1,(block_pair[refidx].start+block_pair[refidx].length),gap_strip)
            

class block:
    def __init__(self,start,length,src,strand,pairBlock='',seq=''):
        p = re.compile('[A,T,C,G,a,t,c,g,N,n]') # regexp for telling non-gap characters
        self.start = start
        self.length = length
        self.src = src
        self.strand = strand
        length2 = length
        self.seq = seq
        self.pair = pairBlock
        self.coords = []
        coord = start
        for i in range(length):
            if p.match(self.seq[i]): 
                coord += 1
                self.coords.append(coord)
            elif self.seq[i] == '-':
                self.coords.append('gap')
                length2 -= 1
            else:
                print 'something else'
        self.length = length2 -1
        self.start = start + 1
            

    def contains(self,q_coord1,q_coord2=-1):
            if (q_coord1 in range(self.start,self.start+self.length+1)):
                if (q_coord2 in range(self.start,self.start+self.length+1)) | (q_coord2 == -1):
                    return True
                else:
                    return False
    def subSeq(self,coord1,coord2,gap_strip='n'):
        p = re.compile('[A,T,C,G,a,t,c,g,N,n]') # regexp for telling non-gap characters
        #check if contains coordinates
        if self.contains(coord1,coord2) == False:
            raise ValueError('These coordinates are not part of this alignment block...')
        #first get out the ungapped subsequence from block and its paired alignment
        subseq = self.seq[self.coords.index(coord1):self.coords.index(coord2)+1]
        subseq_pair = self.pair.seq[self.coords.index(coord1):self.coords.index(coord2)+1]

        if gap_strip =='n':
            return subseq,subseq_pair
        elif gap_strip == 'y':
            subseq_nogap = []
            subseq_pair_nogap =[]
            for i in range(len(subseq)):
                if p.match(subseq[i]):
                    subseq_nogap.append(subseq[i])
                    subseq_pair_nogap.append(subseq_pair[i])
            return ''.join(subseq_nogap),''.join(subseq_pair_nogap)
        
        

#read a file and initialize a pairMaf
def readPairMAF(maf_file):
    blocks_seqs = {}
    maf = pairMaf()
    f = open(maf_file,'r')
    for line in f:
        if '#' in line: continue
        # if you're at the start of a block
        if line[0] == 'a':
            #print next(f).strip().split()
            s, src1, start1, size1, strand1, srcSize1, text1 = next(f).strip().split()
            s, src2, start2, size2, strand2, srcSize2, text2 = next(f).strip().split()
            #minus strand coordinates are relative to rev-com: fix
            start1,start2,size1,size2,srcSize1,srcSize2 = map(int,[start1,start2,size1,size2,srcSize1,srcSize2])
            
            if strand1 == '-':
                start1 = srcSize1 - start1
            if strand2 == '-':
                start2 = srcSize2 - start2

            block1 = block(start1,size1,src1,strand1,pairBlock='',seq=text1)
            block2 = block(start2,size2,src2,strand2,pairBlock=block1,seq=text2)
            block1.pair = block2
            blocks_seqs[(block1,block2)] = 1
    maf.blocks = blocks_seqs
    maf.seqSizes = (srcSize1,srcSize2)
    return maf
