'''
extractIPD.py
Colin McNally
2020/03/26

For a given sample to process, extract the IPD value for each base, and save the T IPDs and the binarized versions
'''


from __future__ import print_function, division
import sys, os
import pandas as pd
import numpy as np
import re
import argparse
import pbcore.io as pb
from pbcore.sequence import reverseComplement
from Bio import Seq, SeqIO
from tqdm import tqdm
import edlib
import pickle
from sklearn.mixture import GaussianMixture
from sklearn.exceptions import ConvergenceWarning
import multiprocessing as mp
import glob
import warnings
import Queue
import time

parser = argparse.ArgumentParser(description="Extract IPD values from an input sample")
parser.add_argument('referenceFile', nargs='?', default='/avicenna/vramani/analyses/pacbio/pbrun3456_SampleReference.csv',
                   help='The sample reference file from which to []')
parser.add_argument('sample', 
                    help='Either an integer index into the reference file, or a string identifier in the format [cell].[samplename]')
parser.add_argument('-o', '--outputlocation', help='Location to save the outputs')
parser.add_argument('-j', '--threads', type=int, help='Number of threads to use. Defaults to 1')
parser.add_argument('-q', '--quiet', action='store_true', help='Add to supress the progress bar')

usepercentiles = range(10,41)
baserc = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
p = re.compile(r'\d+[=XDI]')


# Extract IPD for each molecule from an array sample, return onlyT results in a numpy array
def extractIPDonlyTarray(cbamFile, bamFile, refFile):
    bam = pb.IndexedBamReader(bamFile)
    cbam = pb.IndexedBamReader(cbamFile)
    
    for ir, record in enumerate(SeqIO.parse(refFile, 'fasta')):
        if ir > 0:
            raise InputError('Reference fasta has multiple entries')
        refseq = record.seq

    #Find bases that are A or T
    refisa = [b == 'A' for b in refseq]
    refist = [b == 'T' for b in refseq]
    
    #find all ZMW with ccs that pass filter
    validzmw = []
    for cc in tqdm(cbam, position=0, leave=True, desc='Filtering CCS'):
        if (cc.readLength - len(refseq)) >= -2 and (cc.readLength - len(refseq)) <= 2:
            subrs = bam.readsByHoleNumber(cc.HoleNumber)
            if len(subrs) < 2:
                continue
            subror = np.array([x.isForwardStrand for x in subrs])
            if np.sum(subror == True) >= 3 and np.sum(subror == False) >= 3:
                validzmw.append(cc.HoleNumber)
    
    tipdms = np.empty((len(validzmw), len(refseq)), dtype='float32') # tipdms = T IPD means
    tipdms.fill(np.nan)
    
    for zind, zmw in enumerate(tqdm(validzmw, position=0, leave=True, desc='Getting IPDs')):
        subrs = bam.readsByHoleNumber(zmw)
        allipds = np.empty((len(subrs), len(refseq)))
        allipds.fill(np.nan)
        allpws = np.empty((len(subrs), len(refseq)))
        allpws.fill(np.nan)
        subrOrient = np.empty(len(subrs))

        for index, sr in enumerate(subrs):
            subrOrient[index] = sr.isForwardStrand
            cigarv = sr.unrolledCigar(orientation="genomic")
            #using subread base calls that are matches or mismatches
            #may want to consider ignoring mismatches in the future, especially for non-synethetic source DNA
            #insertions = 1, gaps = 2, leaving both out of downstream analysis
            usebases =  cigarv == 7 #np.logical_or(cigarv == 7, cigarv == 8)
            ipds = sr.baseFeature('Ipd',aligned=True, orientation="genomic")
            pws = sr.baseFeature('PulseWidth',aligned=True, orientation="genomic")
            refpos = sr.referencePositions(aligned=True, orientation="genomic")
            #use the reference positions to index the per base IPD values into their position in the reference sequence
            allipds[index, refpos[usebases]] = ipds[usebases]
            allpws[index, refpos[usebases]] = pws[usebases]

        allipds = allipds / np.mean(np.percentile(allipds[~np.isnan(allipds)],usepercentiles))
        
        # I expect to see RuntimeWarnings in this block
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            formean = np.nanmean(allipds[subrOrient == True,], axis=0)
            revmean = np.nanmean(allipds[subrOrient == False,], axis=0)
            
        #Take IPD values from the reverse strand if the reference is an A
        tipdms[zind, refisa] = revmean[refisa]
        #Take IPD values from the forward strand if the reference is a T
        tipdms[zind, refist] = formean[refist]
        
    return(tipdms)

# Fit a gaussian mixture model to the ipd of a single molecule, return the binarized values
def fitGaussian(ipds, initmean = False, zmw=-1):
    with np.errstate(invalid='ignore'):
        ridx = np.where(ipds > 0)[0]

    lrnn = np.log(ipds[ridx]).reshape(-1,1)
    
    gmm1 = GaussianMixture(1, covariance_type='spherical')
    if initmean:
        gmm2 = GaussianMixture(2, covariance_type='spherical', means_init=np.array([1.1, 3]).reshape((2,1)),
                              weights_init=np.array([.85, .15]), tol=1e-6)
    else:
        gmm2 = GaussianMixture(2, covariance_type='spherical')
    
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=ConvergenceWarning)
        gmm1.fit(lrnn)
        if not gmm1.converged_:
            print('zmw #%d did not converge on gmm1' % (zmw))
        gmm2.fit(lrnn)
        convround = 0
        if not gmm2.converged_:
            gmm2 = GaussianMixture(2, covariance_type='spherical', means_init=np.array([1.1, 3]).reshape((2,1)),
                                  n_init=10, init_params='random', tol=1e-5, max_iter=200)
            gmm2.fit(lrnn)
            convround = 1
            if not gmm2.converged_:
                gmm2 = GaussianMixture(2, covariance_type='spherical', means_init=np.array([1.1, 3]).reshape((2,1)),
                                      n_init=20, init_params='random', tol=1e-4, max_iter=400)
                gmm2.fit(lrnn)
                convround = 2
                if not gmm2.converged_:
                    gmm2 = GaussianMixture(2, covariance_type='spherical', means_init=np.array([1.1, 3]).reshape((2,1)),
                                      n_init=40, init_params='random', tol=1e-3, max_iter=600)
                    gmm2.fit(lrnn)
                    convround = 3
                    if not gmm2.converged_:
                        convround = 9
                        print('zmw #%d did not converge on gmm2 even after extensions' % (zmw))
    aicdif = gmm1.aic(lrnn) - gmm2.aic(lrnn)
    mixmeans = gmm2.means_.flatten()
    mixweights = gmm2.weights_.flatten()
    elevstate = np.argmax(mixmeans)
    resp = gmm2.predict_proba(np.log(ipds[ridx].reshape(-1,1)))
    respfull = np.empty(len(ipds), dtype='float32')
    respfull.fill(np.nan)
    respfull[ridx] = resp[:,elevstate]
    convergInf = str(convround) + '.' + str(gmm2.n_iter_)
    
    return (respfull, np.array([mixmeans[1-elevstate], mixmeans[elevstate]]), 
            np.array([mixweights[1-elevstate], mixweights[elevstate]]), aicdif, convergInf)

# binarize each molecule in an array sample, return as an numpy array
def binarizeIPD(tipdms):
    binargmm = np.empty(tipdms.shape, dtype='float32')
    binargmm.fill(np.nan)
    
    for izm in tqdm(range(tipdms.shape[0]), position=0, leave=True, desc='Binarizing'):
        respf, gmmeans, gmweights, aic, convinf = fitGaussian(tipdms[izm,:], initmean=True)
        binargmm[izm,:] = respf
    
    return binargmm

# Extract the IPD from each zmw passed in the zmwqueue, return a dictionary containing onlyT IPDs, the binarized versions,
# and other informatin to include in the zmwinfo file
def extractIPDonlyTgenomic(cbamfile, alncbamfile, sbamfile, zmwqueue, outqueue): #filebase, nit):
    with pb.IndexedBamReader(alncbamfile) as alncbam, pb.IndexedBamReader(sbamfile) as sbam,\
         pb.IndexedBamReader(cbamfile) as cbam:
        zmw = zmwqueue.get()
        while zmw is not None:

            cc = cbam.readsByHoleNumber(zmw)[0]
            ccread = cc.read(aligned=False, orientation='native')
            subrs = sbam.readsByHoleNumber(zmw)
            ccalns = alncbam.readsByHoleNumber(zmw)
            
            zmres = {}
            
            if len(ccalns) > 0:
                alnlen = np.array([ccal.readLength for ccal in ccalns])
                usealn = np.where(alnlen == np.max(alnlen))[0][0]
            elif len(ccalns) == 1:
                usealn = 0
            elif len(ccalns) == 0:
                usealn = None

            allipds = np.empty((len(subrs), len(ccread)), dtype='float32')
            allipds.fill(np.nan)
            subrOrient = np.empty(len(subrs), dtype='bool')

            for index, sr in enumerate(subrs):
                # Test if this subread aligns to the forward or reverse of the CCS
                forwardSread = sr.read(aligned=False, orientation='native')
                reverseSread = reverseComplement(forwardSread)
                faln = edlib.align(forwardSread, ccread, mode='NW', task='path')
                raln = edlib.align(reverseSread, ccread, mode='NW', task='path')
                if faln['editDistance'] < raln['editDistance']:
                    subrOrient[index] = True
                    alndir = faln
                    useread = forwardSread
                else:
                    subrOrient[index] = False
                    alndir = raln
                    useread = reverseSread

                # Use the alignment information to extract IPD at each base that aligns to the CCS
                origb = np.empty(len(useread), dtype=np.int16 )
                origb.fill(np.nan)
                ccb = np.empty(len(useread), dtype=np.int16)
                ccb.fill(np.nan)
                subI = 0
                ccI = 0
                for m in p.finditer(alndir['cigar']):
                    lg = int(m.group()[-len(m.group()):-1])
                    mtype = m.group()[-1]
                    if mtype == '=':
                        origb[subI:(subI + lg)] = range(subI, subI + lg)
                        ccb[subI:(subI + lg)] = range(ccI, ccI + lg)
                        subI += lg
                        ccI += lg
                    elif mtype == 'X':
                        subI += lg
                        ccI += lg
                    elif mtype == 'I':
                        subI += lg
                    elif mtype == 'D':
                        ccI += lg

                ccb = ccb[~np.isnan(ccb)]
                origb = origb[~np.isnan(origb)]
                if not subrOrient[index]:
                    for i in range(len(origb)):
                        origb[i] = -1 - origb[i]

                ipds = sr.baseFeature('Ipd',aligned=False, orientation="native")
                allipds[index, ccb] = ipds[origb]

            # Normalize the IPD values from each subread
            allipds = allipds / np.mean(np.percentile(allipds[~np.isnan(allipds)],usepercentiles))

            readisb = {refb:np.where([b == refb for b in ccread])[0] for refb in ['A','C','G','T']}

            # Take the mean IPD at each position
            with warnings.catch_warnings(): # ignoring warnings from taking the mean of columns that are all NaN
                warnings.simplefilter("ignore", category=RuntimeWarning)
                forwardMean = np.nanmean(allipds[subrOrient == True,:], axis=0)
                reverseMean = np.nanmean(allipds[subrOrient == False,:], axis=0)
            # get the mean at just Ts
            tonlyMean = np.empty(len(ccread), dtype='float32')
            tonlyMean.fill(np.nan)
            tonlyMean[readisb['T']] = forwardMean[readisb['T']]
            tonlyMean[readisb['A']] = reverseMean[readisb['A']]

            # Save useful information about this molecule
            zmres['zmw'] = zmw
            zmres['cclen'] = len(ccread)
            zmres['nsubr'] = len(subrs)
            zmres['naln'] = len(ccalns)
            if usealn is not None:
                zmres['chr'] = ccalns[usealn].referenceName
                zmres['refStart'] = ccalns[usealn].referenceStart
                zmres['refEnd'] = ccalns[usealn].referenceEnd
                zmres['alnStart'] = ccalns[usealn].aStart
                zmres['alnEnd'] = ccalns[usealn].aEnd
            else:
                zmres['chr'] = "noAlignment"
                zmres['refStart'] = -1
                zmres['refEnd'] = -1
                zmres['alnStart'] = -1
                zmres['alnEnd'] = -1
            with warnings.catch_warnings(): # ignoring warnings from taking the mean of all NaN
                warnings.simplefilter("ignore", category=RuntimeWarning)
                # Something is wrong with the read if these throw a warning, but still no need to print on command line
                # The stored value will be NaN
                zmres['basemeanA'] = np.nanmean(np.concatenate([forwardMean[readisb['A']], reverseMean[readisb['T']]],axis=None))
                zmres['basemeanC'] = np.nanmean(np.concatenate([forwardMean[readisb['C']], reverseMean[readisb['G']]],axis=None))
                zmres['basemeanG'] = np.nanmean(np.concatenate([forwardMean[readisb['G']], reverseMean[readisb['C']]],axis=None))
                zmres['basemeanT'] = np.nanmean(np.concatenate([forwardMean[readisb['T']], reverseMean[readisb['A']]],axis=None))

            zmres['onlyt'] = tonlyMean
            
            # binarize the onlyT ipds
            respf, gmmeans, gmweights, aic, convinf = fitGaussian(tonlyMean, initmean=True, zmw=zmw)
            zmres['bingmm'] = respf
            zmres['gmmlowmean'] = gmmeans[0]
            zmres['gmmlowweight'] = gmweights[0]
            zmres['gmmhighmean'] = gmmeans[1]
            zmres['gmmhighweight'] = gmweights[1]
            zmres['aicDiff'] = aic
            zmres['convergeInf'] = convinf
            
            outqueue.put(zmres) # put the results in the output queue
            zmwqueue.task_done()
            zmw = zmwqueue.get()
        zmwqueue.task_done()
        outqueue.put(None)


# Collect genomic results for each molecule, aggregate them in dictionaries and a data frame and save the results periodically
def listenerSaver(zmwqueue, outqueue, nthread, outbase, sampcn): #, totalccs):   
    MipdDic = {}
    BingDic = {}
    ccdicdat = {}
    for key in ['zmw', 'chr', 'refStart', 'refEnd', 'alnStart', 'alnEnd', 'naln', 'cclen', 'nsubr', 'basemeanA', 'basemeanC', 
                'basemeanG', 'basemeanT', 'gmmlowmean', 'gmmlowweight', 'gmmhighmean', 'gmmhighweight', 'aicDiff', 'convergeInf']:
        ccdicdat[key] = []
    seenEnd = 0
    savepart = 0
    lastqs = zmwqueue.qsize()
    #time.sleep(10)
    pbar = False
    while seenEnd < nthread:
        try:
            nextres = outqueue.get()
        except Queue.Empty:
            time.sleep(0.1)
            pass
        else:
            if nextres is None:
                seenEnd += 1
            else:
                if not pbar:
                    pbar = tqdm(desc='Processing zmws', total=(lastqs - nthread), smoothing=0.01)
                pbar.update(1)
                MipdDic[nextres['zmw']] = nextres['onlyt']
                BingDic[nextres['zmw']] = nextres['bingmm']
                
                for key in ccdicdat.keys():
                    ccdicdat[key].append(nextres[key])
                    
                if len(MipdDic.keys()) > 10000:
                    ccdat = pd.DataFrame(ccdicdat)
                    #reordering zmw info column names
                    ccdat = ccdat[['zmw', 'cclen', 'nsubr', 'chr', 'refStart', 'refEnd', 'alnStart', 'alnEnd', 'naln', 'basemeanA',
                                   'basemeanC', 'basemeanG', 'basemeanT', 'gmmlowmean', 'gmmlowweight', 'gmmhighmean', 'gmmhighweight',
                                   'aicDiff', 'convergeInf']] 
                    #need to convert strings to unicode or file can't be loaded in R
                    ccdat.chr = ccdat.chr.astype('unicode') 

                    #save files, one with zmw info, one with a dictionary of onlyT IPD along CCS
                    ccdat.to_pickle(os.path.join(outbase, 'processed', 'onlyT', 'tmp.' + sampcn + '_part' + str(savepart) +
                                                 '_onlyT_zmwinfo.pickle'))
                    with open(os.path.join(outbase, 'processed', 'onlyT', 'tmp.' + sampcn + '_part' + str(savepart) + 
                                           '_onlyT.pickle'), 'wb') as fout:
                        pickle.dump(MipdDic, fout, pickle.HIGHEST_PROTOCOL)
                    with open(os.path.join(outbase, 'processed', 'binarized', 'tmp.' + sampcn + '_part' + str(savepart) +
                                           '_bingmm.pickle'), 'wb') as fout:
                        pickle.dump(BingDic, fout, pickle.HIGHEST_PROTOCOL)
                    savepart += 1
                    MipdDic = {}
                    BingDic = {}
                    ccdicdat = {}
                    for key in ['zmw', 'chr', 'refStart', 'refEnd', 'alnStart', 'alnEnd', 'naln', 'cclen', 'nsubr', 'basemeanA',
                                'basemeanC', 'basemeanG', 'basemeanT', 'gmmlowmean', 'gmmlowweight', 'gmmhighmean', 'gmmhighweight',
                                'aicDiff', 'convergeInf']:
                        ccdicdat[key] = []
    pbar.close()
    ccdat = pd.DataFrame(ccdicdat)
    #reordering zmw info column names
    ccdat = ccdat[['zmw', 'cclen', 'nsubr', 'chr', 'refStart', 'refEnd', 'alnStart', 'alnEnd', 'naln', 'basemeanA', 'basemeanC',
                   'basemeanG', 'basemeanT', 'gmmlowmean', 'gmmlowweight', 'gmmhighmean', 'gmmhighweight', 'aicDiff']] 
    #need to convert strings to unicode or file can't be loaded in R
    ccdat.chr = ccdat.chr.astype('unicode') 
    
    #save files, one with zmw info, one with a dictionary of onlyT IPD along CCS
    nit = 0
    ccdat.to_pickle(os.path.join(outbase, 'processed', 'onlyT', 'tmp.' + sampcn + '_part' + str(savepart) + '_onlyT_zmwinfo.pickle'))
    with open(os.path.join(outbase, 'processed', 'onlyT', 'tmp.' + sampcn + '_part' + str(savepart) + '_onlyT.pickle'), 'wb') as fout:
        pickle.dump(MipdDic, fout, pickle.HIGHEST_PROTOCOL)
    with open(os.path.join(outbase, 'processed', 'binarized', 'tmp.' + sampcn + '_part' + str(savepart) + '_bingmm.pickle'), 'wb') as fout:
        pickle.dump(BingDic, fout, pickle.HIGHEST_PROTOCOL)
    
    
def main():
    args = parser.parse_args()
    sampleRef = pd.read_csv(args.referenceFile, index_col = 'index')
    strre = re.match('(.+)\.(.+)', args.sample)
    if re.match('[\d]+', args.sample) is not None:
        sampleInfo = sampleRef.loc[int(args.sample)]
    elif strre is not None:
        sampCell = strre.groups(0)[0]
        sampName = strre.groups(0)[1]
        sampleInfo = sampleRef.query('cell==@sampCell and sampleName==@sampName').squeeze()
    else:
        raise ValueError("Sample specification must be integer or '[cell].[samplename]'")
    if sampleInfo.shape[0] == 0:
        raise ValueError("No sample in the reference was matched by your sample specification")
    
    # Get output location, set up folders
    if args.outputlocation is not None:
        outBase = args.outputlocation
    else:
        outBase = os.path.join('/avicenna/vramani/analyses/pacbio', sampleInfo.cell)
        
    # make output folders if they don't exist
    if not os.path.exists(os.path.join(outBase, 'processed')):
        os.makedirs(os.path.join(outBase, 'processed'))
    if not os.path.exists(os.path.join(outBase, 'processed','onlyT')):
        os.makedirs(os.path.join(outBase, 'processed','onlyT'))
    if not os.path.exists(os.path.join(outBase, 'processed','binarized')):
        os.makedirs(os.path.join(outBase, 'processed','binarized'))
        
    if args.threads is not None:
        nthreads = args.threads
    else:
        nthreads = 1
    
    sampCN = sampleInfo.cell + '_' + sampleInfo.sampleName
    print("Sample: %s" % (sampCN))
    print("Process ID: %d" % (os.getpid()))
    
    if re.search('linker', sampleInfo.reference) is not None:
        if nthreads > 1:
            print('Multiple threads not implemented for array data')
        #sample is array / monosequence
        tipds = extractIPDonlyTarray(sampleInfo.ccsFile, sampleInfo.alignedSubreadsFile, sampleInfo.reference)
        np.save(file=os.path.join(outBase, 'processed', 'onlyT', sampCN + '_onlyT'), arr=tipds)
        bingmmar = binarizeIPD(tipds)
        np.save(file=os.path.join(outBase, 'processed', 'binarized', sampCN + '_bingmm'),
                arr=bingmmar)
        
    else:
        # sample is genomic
        alttotal = None #4000 # for testing, only do a subset of molecules
        perchunk = 5000
        # Create a list of all zmw hole numbers
        
        validzmw = []
        
        with pb.IndexedBamReader(sampleInfo.ccsFile) as cbam:
            nccs = len(cbam)
            for ic, cc in enumerate(cbam):
                if alttotal and ic >= alttotal:
                    break
                validzmw.append(cc.HoleNumber)
        
        chunksize = perchunk * nthreads
        
        onchunk = 0
        validzmwc = []
        validzmwc.append([])
        for zmw in validzmw:
            validzmwc[onchunk].append(zmw)
            if len(validzmwc[onchunk]) == chunksize:
                onchunk += 1
                validzmwc.append([])
        del validzmw
        nchunks = len(validzmwc)
        
        validzmwQ = mp.JoinableQueue()
        for chunk in range(nchunks):
            for i in range(len(validzmwc[chunk])):
                validzmwQ.put(validzmwc[chunk][i])
            for thr in range(nthreads):
                validzmwQ.put(None)
        del validzmwc
        
        # Start a process to read output and save to disk
        procq = mp.Queue()
        listenp = mp.Process(target=listenerSaver, args=(validzmwQ, procq, nthreads*nchunks, outBase, sampCN))
        listenp.daemon = True
        listenp.start()
        
        
        # start processes to extract IPD from each molecule
        
        for chunk in range(nchunks):
            workers = []
            for i in range(nthreads):
                x = mp.Process(target=extractIPDonlyTgenomic, args=(sampleInfo.ccsFile, sampleInfo.alignedCcsFile,
                                                                 sampleInfo.unalignedSubreadsFile, validzmwQ, procq))
                x.daemon = True
                workers.append(x)
                x.start()
            for thr in workers:
                thr.join()
        
        # wait until all extractIPDonlyTgenomic processes are finished
        listenp.join()

        # Results are now split between many tmp files. Identify them and join them together
        inffiles = [os.path.basename(x) for x in glob.glob(os.path.join(outBase, 'processed', 'onlyT', 'tmp.' +  sampCN +
                                                                        '_part*_onlyT_zmwinfo.pickle'))]
        onlytfiles = [os.path.basename(x) for x in glob.glob(os.path.join(outBase, 'processed', 'onlyT', 'tmp.' + sampCN +
                                                                          '_part*_onlyT.pickle'))]
        binfiles = [os.path.basename(x) for x in glob.glob(os.path.join(outBase, 'processed', 'binarized', 'tmp.' + sampCN +
                                                                        '_part*_bingmm.pickle'))]

        zminfoC = pd.DataFrame()
        zmtipdC = {}
        zmbinC = {}
        
        # combine zmw information files
        for inff in inffiles:
            zminfoC = pd.concat([zminfoC, pd.read_pickle(os.path.join(outBase,'processed','onlyT',inff))], sort=False)
        # reset indices for combined zmwinfo dataframe
        zminfoC = zminfoC.sort_values('zmw')
        zminfoC.reset_index(drop=True, inplace=True)
        
        # combine onlyT and binarized files
        for onlytf in onlytfiles:
            with open(os.path.join(outBase,'processed','onlyT',onlytf),'rb') as fopen:
                ottemp = pickle.load(fopen)
            zmtipdC.update(ottemp)
            
        for binf in binfiles:
            with open(os.path.join(outBase,'processed','binarized',binf),'rb') as fopen:
                ottemp = pickle.load(fopen)
            zmbinC.update(ottemp)
            
        # write combined files to pickle
        zminfoC.to_pickle(os.path.join(outBase, 'processed', 'onlyT', sampCN + '_onlyT_zmwinfo.pickle'))
        with open(os.path.join(outBase, 'processed', 'onlyT', sampCN + '_onlyT.pickle'), 'wb') as fout:
            pickle.dump(zmtipdC, fout, pickle.HIGHEST_PROTOCOL)
        with open(os.path.join(outBase, 'processed', 'binarized', sampCN + '_bingmm.pickle'), 'wb') as fout:
            pickle.dump(zmbinC, fout, pickle.HIGHEST_PROTOCOL)

        # delete the temporary files
        os.system('rm ' + os.path.join(outBase, 'processed', 'onlyT','tmp.' + sampCN + '_part*'))
        os.system('rm ' + os.path.join(outBase, 'processed', 'binarized','tmp.' + sampCN + '_part*'))
    
if __name__ == "__main__":
    # Usage: extractIPD.py [referencefile] sampleIdentifier
    main()