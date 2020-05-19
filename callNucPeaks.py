from __future__ import print_function, division
import sys, os
import numpy as np
import pbcore.io as pb
from Bio import Seq, SeqIO
from tqdm import tqdm
from scipy.signal import find_peaks
import pandas as pd

def callNucPeaks(npyf, refFile, sampname, heightcutoff=-10, windowsize=133):
    tipdms = np.load(npyf)
    
    for ir, record in enumerate(SeqIO.parse(refFile, 'fasta')):
        if ir > 0:
            raise InputError('Reference sequence has multiple entries')
        refseq = record.seq

    window = windowsize
    windowhalf = int((window-1)/2)
    #windmids = range(windowhalf, len(refseq)-windowhalf)

    peakdic = {'zmw':[], 'pos':[], 'height':[], 'prominence':[]}
    for izm in range(tipdms.shape[0]):
        rollingT = [np.nanmean(tipdms[izm, max(x-windowhalf,0):min(x+windowhalf+1,len(refseq))]) for x in range(0,len(refseq))]
        rollingT = np.asarray(rollingT)
        peaks, peakinf = find_peaks(-1 * rollingT, distance=147, height=heightcutoff, prominence=0.01)
        for ip, p in enumerate(peaks):
            peakdic['zmw'].append(izm)
            peakdic['pos'].append(peaks[ip])
            peakdic['height'].append(peakinf['peak_heights'][ip])
            peakdic['prominence'].append(peakinf['prominences'][ip])
            
    peakdf = pd.DataFrame(peakdic)
    peakdf.to_feather(sampname + '_peaks.feather')


def main():
    usenamepath = sys.argv[3] + '_' + sys.argv[4]
    print(usenamepath)
    heightparam = float(sys.argv[5]) if len(sys.argv) > 5 else -10
    windowparam = int(sys.argv[6]) if len(sys.argv) > 6 else 133
    callNucPeaks(sys.argv[1], sys.argv[2], usenamepath, heightparam)
    
    
if __name__ == "__main__":
    # Usage: ./callNucPeaks.py [onlyTipd.npy] [reference.fasta] [cell] [sampleName]
    main()
