import os,sys,re
import scipy
import numpy as np
import pysam 
from collections import Counter
from intervaltree import Interval, IntervalTree
import gzip as gz

def readIterator(filenames,chrom, start, end):
  for bamfile in filenames:
    if os.path.exists(bamfile) and (os.path.exists(bamfile.replace(".bam",".bai")) or os.path.exists(bamfile+".bai")):
      input_file = pysam.Samfile( bamfile, "rb" )
      for read in input_file.fetch(chrom, start, end):
        yield read
      input_file.close()



def define_large_bins(chromsizes, resolutions):
    bins = {}
    valid_chroms = {}
    lines = chromsizes.readlines()
    for resolution in resolutions:
        bins[resolution] = {}
    for resolution in resolutions:
        hindex = 0
        for line in lines:
            chromname, length = line.split()
            valid_chroms[chromname] = True
            for i in range(0,int(length),resolution):
                bins[resolution][(chromname, i)] = hindex
                hindex += 1
    return bins, valid_chroms


def calculatePerBase(filenames, tss, valid_chroms, label):
    '''the goal here is to take a list of BED-file formatted features & 
    extract all ZMWs (from aligned CCS) that match the bed entries.'''
    chc = False
    zmws = []
    for line in tss:
        split = line.split()
        #check if there is an optional label (for certain flat files)
        chrom, start, end = split[0], int(split[1]), int(split[2])
        if chrom not in valid_chroms: continue
        for read in readIterator(filenames, chrom, start, end):
            rname = read.qname
            rstart = read.reference_start
            rend = read.reference_end
            zmw = rname.split('/')[1]
            zmws.append("%s\t%s\t%s\t%s\t%s" % (zmw, chrom, rstart, rend, label))
    return zmws

def main():
    valid = {}
    tfile = gz.open(sys.argv[1])
    valid_chroms = open(sys.argv[2])
    for line in valid_chroms:
        split = line.split()
        valid[split[0]] = True
    tss = tfile.readlines()
    filename = sys.argv[3]
    label = sys.argv[4]
    filenames = [filename]
    zmws = calculatePerBase(filenames, tss, valid, label)
    for zmw in zmws:
        print(zmw)

if __name__ == "__main__":
    main()
