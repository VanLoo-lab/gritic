# GRITIC
A tool for timing complex copy number gains. Provides gain timing estimates for segments with a total copy number of up to 9. Only copy number segments with 10 or more SNVs will be timed.
## Installation
First clone this repository and load anaconda onto your system. Then create a new conda environment
```
conda create -n gritic_env python=3 numba scipy pandas numpy networkx matplotlib
```
Then activate the environment
```
conda activate gritic_env
```
Note that this environment can be deactivated in the future using
```
conda activate
```

## Running
The easiest way to run GRITIC on a single sample is to use the following command
```
python /PATH_TO_GRITIC_DIR/rungritic_cmd.py -ARGS
```
GRITIC currently has six required arguments and one optional argument. These  are
- ```--mutation_table``` A path to the SNV table for the sample.
- ```--copy_number_table``` A path to the copy number table for the sample.
- ```--purity``` The estimated cellular purity for the sample.
- ```--wgd_status``` The whole genome duplication status for the sample (either True or False).
- ```--sample_id``` Sample ID, for output purposes only.
- ```--output``` The directory for the output, if it doesn't exist it will be created.
- ```--subclone_table``` (Optional) A path to the subclone table for the sample. If not provided it is assumed every SNV is clonal.

An example command to get GRITIC to run from within its directory is:

```
python rungritic.py --mutation_table examples/input_data/snv_table_example.tsv --subclone_table examples/input_data/subclone_table_example.tsv --copy_number_table examples/input_data/cn_table_example.tsv  --purity 0.39 --wgd_status T --output examples/output --sample_id TEST_ID
```


## Input Table Formats
The three input tables that GRITIC requires should be tab separated. Examples using simulated data are available in the example directory. We currently filter out any non-autosomal chromosomes.
### Mutation Table 
All SNVs for the sample. Requires the column names Chromosome, Position, Tumor_Ref_Count & Tumor_Alt_Count. 
### Copy Number Table 
The rounded allele-specific copy number profile for the sample. Requires the column names Chromosome, Segment_Start, Segment_End, Major_CN & Minor_CN. 

GRITIC will merge adjacent segments with the same allele specific copy number.
### Subclone Table (*Optional*)
The identified subclonal peaks and their relative sizes for the sample. All peaks with a cancer cell fraction of less than 1 should be included. GRITIC assumes a peak with a cancer cell fraction of 1.

This table requires the column names Subclone_CCF (the cancer cell fraction of the subclone) & Subclone_Fraction (the fraction of total SNVs present in the subclone).

This table is optional, if it is not included GRITIC will assume every SNV is clonal. This will likely bias the gains to be measured earlier.
## Output
### _gain_timing_best.tsv (Main Output)
This table gives estimates of the timing of clonal gains in mutation time with 90% confidence intervals derived from a bootstrapping procedure. For most complex gains there will be multiple rows per segment corresponding to the set of gains required to reach the final state.  The higher the timing estimate the later the gain. A timing of 1 indicates that the gain occurred very close to the emergence of the tumour's most recent common ancestor and 0 that the gain occured close to conception 

The column Pre_WGD_Probability gives the probability a gain occurred before a WGD. First_Independent_Gain_Probability gives the probability the gain was the first non-WGD gain for the segment.
### _gain_timing_full.tsv
A fuller description of the gain timing. This table gives timing for gains in all possible routes for each segment. The column Best_Route indicates which route is the best fit,  this is the only one present in* \*_gain_timing_best.tsv.*

The column BIC gives the Bayesian Information Criterion for each route. N_Events is the number of events the route takes to reach the final copy number state. Min_N_Events is the minimum number of events 
### tree_plots
Binary tree plots for the gain timings for each route.  Each plot has two binary trees corresponding to each parental allele. Blue nodes represent independent gains, yellow nodes WGD gains and red nodes the final alleles that were present at time of sampling.

