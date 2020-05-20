Processing SAMOSA data
=====================

The following README documents code & analyses performed in (and are largely reproduced from the methods section therein):

Abdulhay NJ, McNally CP, Hsieh LJ, Kasinathan S, Keith A, Estes LL, Karimzadeh M, Underwood JG, Goodarzi H, Narlikar GJ, Ramani V. "Massively multiplex single-molecule oligonucleosome footprinting." bioRxiv (2020)

*All scripts and binaries are provided as is, without any warranty and for use at your own risk. This is not the release of a software package. We are only providing this information and code in addition to a description of methods for making it easier to reproduce our analyses. We are __not__ providing any support for these scripts.* 

The directions below describe the process for starting from BAM files output from a PacBio Sequel I or II machine and eventually calculating a posterior probability of each adenine being methylated. To recapitulate our analysis, the intermediate files can be downloaded at https://doi.org/10.5281/zenodo.3834706. These files will be sufficient input for all of the of analyses shown in SAMOSA_analyses.ipynb

Two use cases
-------------
We apply our method to two use cases in the paper, and they differ in the computational workflow to analyze them. The first is for sequencing samples where every DNA molecule should have the same sequence, which is the case for our *in vitro* validation experiments presented in Figure 1. The second use case is for samples from cells containing varied sequences of DNA molecules. We will refer to the first as homogeneous samples, and the second as genomic samples. The workflow for genomic samples will be presented first in each sections, and the deviations for homogeneous samples detailed at the end.

Sequencing read processing
--------------------------
Sequencing reads were processed using software from Pacific Biosciences. The following describes the workflow for genomic samples.

1.  **Demultiplex reads**
Reads were demultiplexed using [`lima`](https://github.com/PacificBiosciences/barcoding). The flag `--same` was passed as libraries were generated with the same barcode on both ends. This produces a BAM file for the subreads of each sample.
2. **Generate Circular Consensus Sequences (CCS)**
CCS were generated for each sample using [`ccs`](https://github.com/PacificBiosciences/ccs). Default parameters were used other than setting the number of threads with `-j`. This produces a BAM file of CCS.
3. **Align CCS to the reference genome**
Alignment was done using [`pbmm2`](https://github.com/PacificBiosciences/pbmm2), and run on each CCS file, resulting in BAM files containing the CCS and alignment information.
4. **Generate missing indices**
Our analysis code requires pacbio index files (.pbi) for each BAM file. `pbmm2` does not generate index files, so missing indices were generated using `pbindex`.

For homogeneous samples, replace step 3 with this alternate step 3

3. **Align subreads to the reference genome**
[`pbmm2`](https://github.com/PacificBiosciences/pbmm2) was run on each subreads BAM file (the output of step 1) to align subreads to the reference sequence, producing a BAM file of aligned subreads.


Sample Reference Preparation
------------------------
Our script for analyzing samples relies on a CSV file input that contains information about each sample, including the locations of the relevant BAM files and a path to the reference genome. The CSV needs a header with the following columns:

* **index**: Integer indices for each sample. We write the table using `pandas` `.to_csv` function, with parameters `index=True, index_label='index'`
* **cell**: A unique name for the SMRT cell on which the sample was sequenced
* **sampleName**: The name of the sample
* **unalignedSubreadsFile**: This will be the file produced by step 1 above. This should be an absolute path to the file.
* **ccsFile**: This is the file produced by step 2 above
* **alignedSubreadsFile**: This is the file produced by the alternate step 3 above. It is required for homogeneous samples but can be left blank for genomic samples.
* **alignedCcsFile**: This is the file produced by step 3 above. It is required for genomic samples but can be left blank for homogeneous samples.
* **reference**: The file of the reference genome or reference sequence for the sample

Extracting IPD measurements and calling methylation
----------
The script `extractIPD.py` accesses the BAM files, reads the IPD values at each base and uses a gaussian mixture model to generate posterior probabilities of each adenine being methylated.

#### Requirements
`extractIPD` requires the `pbcore` python package, which only runs on Python 2. Thus a Python 2.7 environment is required. The following packages are required:
* `numpy`
* `pandas`
* `pbcore`
* `biopython`
* `tqdm`
* `edlib`
* `scikit-learn`

#### Input
`extractIPD` takes two positional arguments. The first is a path to the above sample reference CSV file. The second is a specification for which sample to run on. This can be either an integer index value, in which case `extractIPD` will run on the corresponding row. Alternatively it can be a string containing the cell and sampleName, separated by a period. Either way `extractIPD` will run on the specified sample using the paths to the BAM files contained within the CSV. The following arguments can be passed:

* **-o**, **--outputlocation**: Pass a path to the location to where the output should be saved. Within this location a 'processed' folder will be created, and then 'onlyT' and 'binarized' folders within that one, where the outputs will be placed. This argument is required unless you change the default output path within the script to one appropriate for your file system.
* **-j**, **--threads**: Pass an integer number of worker threads to use. An additional thread is used to compile the results produced by the workers, but this listener thread uses far less than a full CPU. This argument is ignored when processing homogeneous samples, in which case only a single thread is used.
* **-q**, **--quiet**: Suppress the progress bar.

#### Output
`extractIPD` produces the following three output files when run on genomic samples.

* **`processed/onlyT/{cell}_{sampleName}_onlyT_zmwinfo.pickle`**
This file is a `pandas` dataframe stored as a pickle, and can be read with the `pandas.read_pickle` function. This dataframe contains various information 
* **`processed/onlyT/{cell}_{sampleName}_onlyT.pickle`**
This file contains the normalized IPD value at every thymine. The data is stored as a dictionary object. The keys are the ZMW hole numbers (stored in the column 'zmw' in the zmwinfo dataframe), and the values are numpy arrays. The arrays are 1D with length equal  to the length of the CCS for that molecule. At bases that are A/T, there will be a normalized IPD value. Each G/C base and a few A/T bases for which an IPD value couldn't be measured will contain NaN.
* **`processed/binarized/{cell}_{sampleName}_bingmm.pickle`**
This file contains the posterior probability of each adenine being methylated. The data format is identical to the `_onlyT.pickle` file above, except the numpy array contains values between 0 and 1, where the higher values indicate a higher confidence that the adenine is methylated.

The `_onlyT.pickle` and `_bingmm.pickle` files can be read using the `load` function of the `pickle` package. When loading them in Python 3 use the parameter `encoding="latin1"`.

When run on homogeneous samples the following output files are alternately produced:

* **`processed/onlyT/{cell}_{sampleName}_onlyT.npy`**
This numpy array has a column for every base in the reference sequence, and a row for each DNA molecule that passes the filtering threshold. A normalized IPD value is stored for each adenine that could be measured at A/T bases, other bases are NaN.
* **`processed/binarized/{cell}_{sampleName}_bingmm.npy`**
This numpy array is the same shape as the `_onlyT.npy` file above. The values are posterior probabilities for an adenine being methylated, ranging from 0 to 1.

#### Tips
* As currently implemented, `extractIPD` will interpret a sample to be homogeneous if 'linker' is in the name of the reference, and genomic otherwise.
* When run on genomic data `extractIPD` will not properly terminate in response to CTRL-C. This is a result of the multi-processing implementation. The process can instead be terminated using unix's `kill`. `extractIPD` prints out its own process ID to facilitate this.

Analysis of methylation calls
-----------------
The code for further analyses of our results is contained in the jupyter notebook file SAMOSA_analyses.ipynb. This includes code to generate all results presented in the paper, and to plot each figure. The code is annotated within the notebook.

Extracting ZMWs near predicted transcription factor binding sites
-----------------
We used `zmw_selector.py` to extract all aligned ZMWs where a portion of the alignment falls within 1000 bp of a predicted trancription factor binding site. The script takes 4 command line arguments,
a list of sites formatted as 'chrid \t site \t strand', a list valid chromosomes / chromosome sizes, an aligned, sorted, and indexed BAM file (the aligned CCS reads in this case), and an integer for the window size. The output is a tab-separated file containing the ZMW hole numbers, the alignment start, end, and strand, and the feature strand. zmw_selector.py is run as follows: 

```python ../../zmw_selector.py ${flat} ${chrom_sizes} ${bam} 1000 > ${output}```

Extracting ZMWS that fall within specific BED intervals
-----------------
We used `zmw_selector_bed.py` to extract all aligned ZMWs where a portion of the alignment falls within a BED interval. The script takes 4 command line arguments: a GZIPPED BED File, a list of valid chromosomes and sizes, the aligned, sorted, and indexed CCS BAM file, and label. The output is a tab-separated file containing the ZMW hole numbers, the alignment start and end, and the provided label.

```python ../../zmw_selector_bed.py ${bed.gz} ${chrom_sizes} ${bam} ${bed}  > ${output}```




