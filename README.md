# GRITIC
A tool for timing complex copy number gains in human cancers. Provides gain timing estimates for segments with a total copy number of up to 8+1. 

Each gain timing is measured in mutation time, a scale that ranges from 0 to 1. A timing of 0 indicates that the gain occured close to conception and 1 that the gain occurred very close to the emergence of the tumour's most recent common ancestor.

GRITIC is agnostic to reference genome and works with both human and mouse data.
## Installation
GRITIC can be installed using pip
```
pip install gritic
```
Python 3.X is required.
## Running
The easiest way to run GRITIC on a single sample is with the ```gritic``` command installed by pip.
```
gritic --help
```
This command has five required arguments.
- ```--mutation_table``` A path to the SNV table for the sample.
- ```--copy_number_table``` A path to the copy number table for the sample.
- ```--purity``` The estimated cellular purity for the sample.
- ```--sample_id``` Sample ID, for output purposes only.
- ```--output``` The directory for the output, GRITIC will then store all of the output at --output/SAMPLE_ID


There are also a number of optional arguments.
- ```--subclone_table``` A path to the subclone table for the sample. If not provided it is assumed every SNV is clonal.
- ```--wgd_status``` GRITIC will identify the WGD status of the sample. This can be overidden through this argument.
- ```--non_parsimony_penalty``` Apply a penalty on non-parsimonious routes in WGD tumours. Defalt is False. See [publication](https://aacrjournals.org/cancerdiscovery/article/14/10/1810/748591/The-History-of-Chromosomal-Instability-in-Genome) for more details.
- ```--plot_trees``` Plot the route trees for each segment. Default is True.
- ```--sample_sex``` Override the sample sex with ```XX``` or ```XY```. If omitted, GRITIC imputes ```XY``` when the copy number table contains a Y-chromosome segment and ```XX``` otherwise.
- ```--drop_unmatched_snvs``` In supplied-ID mode, drop mutation rows whose ```Segment_ID``` is absent from the copy number table and issue one warning reporting the number dropped. By default, unmatched IDs raise an error.


A command to run GRITIC using the example data is:

```
gritic --mutation_table examples/snv_table_example.tsv --subclone_table examples/subclone_table_example.tsv --copy_number_table examples/cn_table_example.tsv --purity 0.5 --wgd_status T --output examples/output --sample_id TEST_ID
```

GRITIC can also be run programmatically:

```python
import pandas as pd

from gritic import dataloader, gritictimer, sampletools

copy_number_table, mutation_table = dataloader.load_input_tables(
    'examples/cn_table_example.tsv',
    'examples/snv_table_example.tsv',
)
subclone_table = pd.read_csv(
    'examples/subclone_table_example.tsv',
    sep='\t',
    dtype={
        'Cluster': str,
        'Subclone_CCF': float,
        'Subclone_Fraction': float,
    },
)

sample = sampletools.Sample(
    mutation_table,
    copy_number_table,
    subclone_table,
    sample_id='TEST_ID',
    purity=0.5,
)
gritictimer.process_sample(
    sample,
    output_dir='examples/output',
    plot_trees=True,
    wgd_override=True,
)
```

## Input Table Formats
The three input tables that GRITIC requires should be tab separated. Examples using simulated data are available in the example directory. Data from any allele specific copy number caller, SNV caller and subclone caller can be used as long as the tables are formatted correctly.
### Mutation Table 
All SNVs for the sample. The columns ```Chromosome```, ```Tumor_Ref_Count``` and ```Tumor_Alt_Count``` are always required.

GRITIC supports two segment-assignment modes:

- If both the mutation and copy number tables contain ```Segment_ID```, GRITIC uses those IDs to associate mutations with the input copy number segments. ```Position``` is optional in this mode.
- Otherwise, the mutation table must contain ```Position```, which GRITIC uses to associate mutations with copy number segments.

In supplied-ID mode, every copy number row must have a ```Segment_ID``` that is not null, empty, or whitespace-only, and copy number ```Segment_ID``` values must be globally unique. Every mutation row must also have a ```Segment_ID``` that is not null, empty, or whitespace-only. By default, every mutation ```Segment_ID``` must be present in the copy number table, and ```Chromosome``` must match the corresponding copy number row. If ```--drop_unmatched_snvs``` is supplied, mutation rows whose ```Segment_ID``` is absent from the copy number table are instead dropped with a single warning reporting the number of unmatched SNVs. Missing or blank IDs and chromosome mismatches remain errors.

Programmatic callers can enable the same behavior with ```drop_unmatched_snvs=True``` in ```dataloader.load_input_tables```, or in ```sampletools.Sample``` when constructing a sample directly from DataFrames.

### Copy Number Table 
The rounded allele-specific copy number profile for the sample. Requires the column names ```Chromosome```, ```Segment_Start```, ```Segment_End```, ```Major_CN``` & ```Minor_CN```. ```Major_CN``` must be larger than or equal to ```Minor_CN```. ```Segment_ID``` is optional and it is only used when it is also present in the mutation table.

By default, in both assignment modes, GRITIC merges adjacent segments with the same allele-specific copy number and generates the final segment IDs from ```Chromosome```, ```Segment_Start``` and ```Segment_End``` (for example, ```1-0-199```).

When sample sex is not supplied, GRITIC imputes ```XY``` if the copy number table contains any Y-chromosome segment and ```XX``` otherwise. An explicitly supplied sex takes precedence. This inference assumes that Y-chromosome segments have not been removed from the copy number input.
### Subclone Table (*Optional*)
The identified subclonal peaks and their relative sizes for the sample. All peaks with a cancer cell fraction of less than 0.9 should be included. Any peaks higher than this are automatically filtered. GRITIC assumes a peak with a cancer cell fraction of 1. 

This table requires the column names ```Subclone_CCF``` (the cancer cell fraction of the subclone) & ```Subclone_Fraction``` (the fraction of total SNVs present in the subclone). An ```Cluster``` column is also required as an index for the subclones.

Only subclones with more than 10% of the subclonal SNVs are included. If there are more than two subclones, GRITIC will reformat the sample to have two subclones. The subclone with the largest CCF is unmodified and the remaining clones grouped together. This speeds convergence while maintaining an accurate estimation of clonal vs subclonal SNV composition.

This table is optional, if it is not included GRITIC will assume every SNV is clonal. This will likely bias the gains to be measured earlier.


## Output
GRITIC produces a number of outputs to describe the timing of copy number gains in a given tumour. We recommend only considering gained segments with 10 or more SNVs. 

### _posterior_timing_table_summary.tsv (Main Output)
Gives a summary of the gain timing infotmation for each gained segment. For each gained segment, the sequential gains are labelled with ```Gain_Index```. The median and the 95% percentiles for the gain timing for each sequential gain are reported.

In WGD tumours, the number of gains that arise independently of the WGD will vary depending on the route. Only gains that arise independently of the WGD in at least 80% of posterior samples are reported, this is recorded in the ```Proportion``` column. The timing of the WGD and the probability that each gain arose before the WGD is also recorded.

Two tables are produced, with and without a penalty on non-parsimony. Without a penalty is recommended by default, see our preprint for more details.

### _posterior_timing_table.tsv
This table gives 250 direct samples from the timing posterior for each gained segment. Both the gain and the WGD timing are given in the case of a WGD tumour.

Each individual sample from the posterior for a given segment is labelled with ```Posterior_Sample_Index```. The route for each posterior sample is given by the route column.

Again, as the number of gains that arise independently of the WGD will vary depending on the route, the number of gains per ```Posterior_Sample_Index``` can vary in WGD tumours.

Two tables are produced, with and without a penalty on non-parsimony. Without a penalty is recommended by default, see our preprint for more details.

### _gain_timing_table.tsv
This table gives more details about each possible route for each gained segment. As well as giving the timing for each gains in the route, it also records the route probability that is used when sampling the posterior. The total number of events implied for each route is also given, as well as the density of the timing samples. 

Please see the supplementary materials of the accompanying preprint for more details on the route densities.

### _wgd_calling_info.tsv
A table with a single row that gives WGD calling information for the sample. 

If the most common major copy number by total base pair length (```Major_CN_Mode```) is 1 then the sample is non-WGD. If the ```Major_CN_Mode``` is > 2 then the sample has likely undergone multiple WGDs, but this history is not currently supported by GRITIC.

If the ```Major_CN_Mode``` is 2 then a secondary test is conducted to test the WGD status. All major copy number two segments are timed and the point in mutation time that overlaps with the greatest proportion of segments by width is recorded. If this is greater than 60% then the sample is identified as WGD. The best overlap timing and the proportion are recorded in the columns ```Best_Overlap_Timing``` and ```Overlap_Proportion``` respectively. The timing of the WGD is then jointly estimated from the overlapping segments and recorded in the ```Timing```,```CI_High``` and ```CI_Low``` columns.

### _gain_timing_table_wgd_segments.tsv
This table is produced for WGD tumours. It gives the timing of the major copy number two segments which are assumed to primarily arise through a WGD. The ```Best_Overlap_Timing``` column records the point in mutation time that overlaps with the 89% confidence intervals of the timing of the highest proportion of basepair length of major copy number two segments. The proportion is recorded in ```Proportion``` and whether the given segment's timing confidence intervals intersect with the timing by ```Intersecting```.


### _mutation_table.tsv
The mutation table processed by GRITIC to give additional SNV multiplicity information. SNV multiplicity probabilities are given by the ```Prob_Mult_``` columns. The expected number of mutations with the same coverage as the given SNV that would be expected to have less than three reads supporting the alternate allele is given by the ```Three_Reads_Correction_Mult``` columns. This is approximately the number of SNVs that would be missed by the variant caller at this coverage level.

### _tree_plots
Binary tree plots for the gain timings for each route for a given segment. Each plot has two binary trees corresponding to each parental allele. Blue nodes represent independent gains, yellow nodes WGD gains and red nodes the final alleles that were present at time of sampling.

### _timing_dicts
Python dictionary objects containing the stored gain timing and multiplicity posterior samples for each route for each gained segment. Not necessary for most use cases. Can be read using the ```pickle``` and ```bz2``` modules in python. For example:
```
import pickle
import bz2

with bz2.BZ2File('SAMPLE_ID_timing_dicts/1-0-199_timing_dict.bz2', 'rb') as f:
    timing_dict = pickle.load(f)

```

The keys of each dict correspond to the routes for the sample. Within each route there are ```Timing``` and ```Mult``` keys. The ```Timing``` entry gives the timing samples of the independent gains indexed by the corresponding node in the tree. A WGD timing entry is also given if applicable.

The ```Mult``` entry gives the multiplicity proportions corresponding to each timing sample. It is a N_SamplesxN_Multiplicities numpy array. Across the columns, the multiplicities are orderred from 1 to the major copy number of the segment, followed by the subclonal multiplicity probabilities.
