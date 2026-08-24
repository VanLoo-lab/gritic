# GRITIC
A tool for timing complex copy number gains in cancer genomes. Provides gain timing estimates for segments with a total copy number of up to 8+1.

Each gain timing is measured in mutation time, a scale that ranges from 0 to 1. A timing of 0 indicates that the gain occured close to conception and 1 that the gain occurred very close to the emergence of the tumour's most recent common ancestor.

GRITIC is agnostic to reference genome. The number of numbered autosomes is configurable, and both X/Y and Z/W sex-chromosome systems are supported.
## Installation
GRITIC can be installed using pip
```
pip install gritic
```
Python 3.10 or newer is required.
## Running
The easiest way to run GRITIC on a single sample is with the ```gritic``` command installed by pip.
```
gritic --help
```
This command has five required arguments.
- ```--mutation-table``` A path to the SNV table for the sample.
- ```--copy-number-table``` A path to the copy number table for the sample.
- ```--purity``` The estimated cellular purity for the sample. It must be finite, greater than 0, and no greater than 1.
- ```--sample-id``` Sample ID used as an output-directory component and filename prefix. It must be a cross-platform-safe filename component.
- ```--output``` The directory for the output, GRITIC will then store all of the output at --output/SAMPLE_ID


There are also a number of optional arguments.
- ```--subclone-table``` A path to the subclone table for the sample. If not provided it is assumed every SNV is clonal.
- ```--wgd-count {0,1}``` Override GRITIC's inferred whole-genome-duplication count. Counts below 0 or above 1 are rejected because GRITIC currently supports at most one WGD. If omitted, GRITIC infers the count.
- ```--plot-trees``` Plot the route trees for each segment. This is an opt-in switch and is disabled by default.
- ```--sample-sex``` Override the sample sex with ```XX```, ```XY```, ```ZZ```, or ```ZW```. If omitted, GRITIC infers the sex-chromosome system and karyotype from exact X, Y, Z, and W copy-number chromosome labels.
- ```--autosome-count``` The number of numbered autosomes in the organism. This defines the accepted numbered chromosome labels and the chromosomes eligible for WGD inference. The default is 22.
- ```--drop-unmatched-chromosomes``` Drop copy-number and mutation rows whose chromosome is not one of the configured autosomes or present sex chromosomes, with warnings reporting the number of rows dropped. By default, any such chromosome is an error.
- ```--drop-unmatched-snvs``` Drop mutation rows that cannot be associated with a copy-number segment by either supplied ```Segment_ID``` or genomic ```Position```, with one warning reporting the number dropped. By default, unmatched mutations raise an error.
- ```--merge-adjacent-segments``` Merge coordinate-adjacent copy-number segments having identical ```Major_CN``` and ```Minor_CN```. This is disabled by default.
- ```--min-mutation-alt-count``` Minimum ```Tumor_Alt_Count``` needed to retain a mutation. The default is 3.
- ```--min-mutation-coverage``` Minimum ```Tumor_Ref_Count + Tumor_Alt_Count``` needed to retain a mutation. The default is 10.
- ```--max-subclone-ccf``` Maximum ```Subclone_CCF``` retained as a subclone. The default is 0.9.
- ```--min-subclone-fraction``` A subclone's normalized share of the surviving subclonal mutation fractions must be strictly greater than this threshold. The default is 0.1.

Long option names use hyphens, not underscores. ```Sample_ID``` values are validated rather than sanitized: they cannot be empty, path components such as ```.``` or ```..```, contain path separators, Windows-forbidden filename characters, Unicode control/format/surrogate characters, end in a dot or space, use a reserved Windows device name, or make a derived GRITIC filename exceed the usual 255-unit component limits.


A command to run GRITIC using the example data is:

```
gritic --mutation-table examples/snv_table_example.tsv --subclone-table examples/subclone_table_example.tsv --copy-number-table examples/cn_table_example.tsv --purity 0.5 --wgd-count 1 --output examples/output --sample-id TEST_ID
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
    merge_cn=False,
    min_mutation_alt_count=3,
    min_mutation_coverage=10,
    max_subclone_ccf=0.9,
    min_subclone_fraction=0.1,
    autosome_count=22,
    drop_unmatched_chromosomes=False,
)
gritictimer.process_sample(
    sample,
    output_dir='examples/output',
    plot_trees=True,
    wgd_count=1,
)
```

## Input Table Formats
The required copy-number and mutation tables, and the optional subclone table, should be tab separated. Examples using simulated data are available in the example directory. Data from any allele specific copy number caller, SNV caller and subclone caller can be used as long as the tables are formatted correctly.
### Mutation Table 
All SNVs for the sample. The columns ```Chromosome```, ```Tumor_Ref_Count``` and ```Tumor_Alt_Count``` are always required. A chromosome must be a canonical integer from 1 through ```--autosome-count``` or a sex chromosome present under the selected/inferred karyotype; one leading ```chr``` prefix is accepted and removed. Thus XX accepts X, XY accepts X/Y, ZZ accepts Z, and ZW accepts Z/W. Both read counts must be finite non-negative integers within the signed 64-bit range; their sum must be greater than zero and must also fit within that range. Integer-equivalent spellings are canonicalized to signed-64-bit integers. Every mutation table must also contain either ```Mutation_ID``` or ```Position```. ```Mutation_ID``` is loaded as literal text, so values such as ```000123```, ```NA```, and ```NULL``` are preserved. When present, every ```Position``` must instead be a finite non-negative integer within the signed 64-bit range; integer-equivalent spellings such as ```000123``` or ```123.0``` are canonicalized to ```123```.

Input mutation columns named ```Segment_Start```, ```Segment_End```, ```Major_CN```, or ```Minor_CN``` are ignored because GRITIC always annotates those values from the copy-number table.

GRITIC supports two segment-assignment modes:

- If both the mutation and copy number tables contain ```Segment_ID```, GRITIC uses those IDs to associate mutations with the input copy number segments. ```Position``` is not needed for assignment.
- Otherwise, the mutation table must contain ```Position```, which GRITIC uses to associate mutations with copy-number segments under the zero-based, half-open interval rule ```Segment_Start <= Position < Segment_End```.

In supplied-ID mode, every copy number row must have a ```Segment_ID``` that is not null, empty, or whitespace-only, and copy number ```Segment_ID``` values must be unique within the copy number table. Every mutation row must also have a ```Segment_ID``` that is not null, empty, or whitespace-only. By default, every mutation ```Segment_ID``` must be present in the copy number table, and ```Chromosome``` must match the corresponding copy number row. In position-assignment mode, every mutation must fall inside a copy-number segment on the same chromosome. If ```--drop-unmatched-snvs``` is supplied, either kind of unmatched mutation is instead dropped with one count warning. Missing or blank supplied IDs and chromosome mismatches remain errors.

The selected ```Mutation_ID``` or canonical integer ```Position``` value must be unique within its source segment. An explicit ```Mutation_ID``` column takes precedence for every row; GRITIC does not fall back to ```Position``` for blank values in that column. The same selected value may be reused in different source segments, including source segments that GRITIC subsequently merges.

Programmatic callers enable the same behavior with ```drop_unmatched_snvs=True``` in ```sampletools.Sample```. When using ```dataloader.load_input_tables``` with supplied segment IDs, pass the same value there as well so unmatched IDs are permitted and removed during loading; the command-line interface forwards it to both stages. The separate ```drop_unmatched_chromosomes=True``` API argument corresponds to ```--drop-unmatched-chromosomes``` and filters invalid chromosome rows from both input tables before their numeric fields and segment relationships are validated.

### Copy Number Table 
The rounded allele-specific copy number profile for the sample. Requires the column names ```Chromosome```, ```Segment_Start```, ```Segment_End```, ```Major_CN``` & ```Minor_CN```. Chromosome labels follow the same rules as the mutation table. Segment coordinates are zero-based, half-open, non-negative signed-64-bit integers: ```Segment_Start``` is included, ```Segment_End``` is excluded, and ```Segment_End``` must be strictly greater than ```Segment_Start```. A one-base segment is therefore ```[a, a + 1)``` and has width ```Segment_End - Segment_Start```. Both allele-specific copy numbers must also be finite signed-64-bit integers and non-negative, and ```Major_CN``` must be greater than or equal to ```Minor_CN```. Violations raise an error. ```Segment_ID``` is optional and it is only used when it is also present in the mutation table. Mixed sex-chromosome systems make sex inference ambiguous, and numbered chromosomes above ```--autosome-count``` are unmatched.

GRITIC orders copy-number rows by natural chromosome order, ```Segment_Start```, and ```Segment_End``` immediately after validating the input coordinates. Half-open intervals on the same chromosome must not overlap. By default GRITIC preserves separate copy-number segments and generates final segment IDs from ```Chromosome```, ```Segment_Start```, and ```Segment_End```. Supplying ```--merge-adjacent-segments``` (or ```merge_cn=True``` through the API) merges equal-copy-number intervals only when the following ```Segment_Start``` equals the preceding ```Segment_End```. Gapped intervals are not merged, and the resulting merged interval is revalidated against the same coordinate bounds.

When sample sex is not supplied, an exact Y chromosome label implies ```XY``` and W implies ```ZW```; if neither is present, X implies ```XX``` and Z implies ```ZZ```. If no sex chromosome is represented, GRITIC defaults to ```XX```, so callers using another system should supply ```--sample-sex```. Inputs that mix the X/Y and Z/W systems cannot be inferred and fail. An explicit karyotype takes precedence. For ```XX```, X is present and Y is unmatched; for ```XY```, both X and Y are present. For ```ZZ```, Z is present and W is unmatched; for ```ZW```, both Z and W are present. Unmatched rows fail unless chromosome dropping is enabled. Normal-cell X/Y copy numbers are 2/0 for XX and 1/1 for XY; the corresponding Z/W values are 2/0 for ZZ and 1/1 for ZW. The ```autosome_count``` API argument and ```--autosome-count``` CLI option define the numbered autosomes and default to 22.
### Subclone Table (*Optional*)
The identified subclonal peaks and their relative sizes for the sample. GRITIC filters candidate peaks using the configurable maximum CCF and minimum normalized-fraction thresholds described below.

This table requires the column names ```Subclone_CCF``` (the cancer cell fraction of the subclone) & ```Subclone_Fraction``` (the fraction of total SNVs present in the subclone). A ```Cluster``` column is also required as an index for the subclones. CCF and fraction values must be finite numeric (not boolean) values between 0 and 1 inclusive, and the fractions must sum to no more than 1; invalid prior values fail before filtering.

By default, only subclones with ```Subclone_CCF <= 0.9``` and more than 10% of the surviving subclonal mutation fraction are included; the thresholds are configurable with ```--max-subclone-ccf``` and ```--min-subclone-fraction```. If no subclones remain, GRITIC uses its clonal-only model as if no subclone table was supplied. If there are more than two subclones, GRITIC reformats the sample to have two: the subclone with the largest CCF is unmodified and the remaining clones are grouped together.

GRITIC always derives ```N_SNVs``` from the retained mutation count and each retained/combined ```Subclone_Fraction```, overwriting an input ```N_SNVs``` column. ```Subclone_Fraction``` contributes to the clonal/subclonal mixture prior used during route sampling, while ```Subclone_CCF``` determines the expected VAF of each subclone state.

This table is optional, if it is not included GRITIC will assume every SNV is clonal. This will likely bias the gains to be measured earlier.


## Output
GRITIC produces a number of outputs to describe the timing of copy number gains in a given tumour. We recommend only considering gained segments with 10 or more SNVs. 

Each run requires its ```OUTPUT/SAMPLE_ID``` directory to be absent or empty. GRITIC rejects a nonempty sample directory rather than appending to, or mixing output with, an earlier run.

### _posterior_timing_table_summary_penalty_&lt;True|False&gt;.tsv (Main Output)
Gives a summary of the gain timing information for each gained segment. For each gained segment, the sequential gains are labelled with ```Gain_Index```. The median and the 95% interval for the gain timing of each sequential gain are reported.

In WGD tumours, the number of gains that arise independently of the WGD will vary depending on the route. Only gains that arise independently of the WGD in at least 80% of posterior samples are reported, this is recorded in the ```Proportion``` column. The timing of the WGD and the probability that each gain arose before the WGD is also recorded.

Two summary tables are produced for every run that produces timing output, including non-WGD runs. The filename ending in ```_penalty_False.tsv``` summarizes draws from the ordinary route probabilities. The filename ending in ```_penalty_True.tsv``` summarizes a second set of draws after each route probability has been reweighted once by the post-hoc non-parsimony penalty and renormalized. No command-line option is required to produce either result.

### _posterior_timing_table_penalty_&lt;True|False&gt;.tsv
This table gives 100 direct samples from the timing posterior for each gained segment. Both the gain and the WGD timing are given in the case of a WGD tumour.

Each individual sample from the posterior for a given segment is labelled with ```Posterior_Sample_Index```. The route for each posterior sample is given by the route column.

Again, as the number of gains that arise independently of the WGD will vary depending on the route, the number of gains per ```Posterior_Sample_Index``` can vary in WGD tumours.

Two posterior-sample tables are produced for every run that produces timing output, with the same ```_penalty_False.tsv``` and ```_penalty_True.tsv``` mode tags as their corresponding summary tables. The unpenalized and penalized draws are generated separately. The penalty is applied once, post-hoc, when constructing the ```_penalty_True.tsv``` draw set; it is not incorporated into route fitting or applied again during summarization.

### _gain_timing_table.tsv
This table gives more details about each possible route for each gained segment. As well as giving the timing for each gain in the route, it records the ordinary, unpenalized route probability. The total number of events implied for each route is also given, as well as the density of the timing samples.

The ```Probability``` values in this table are used directly for the ```_penalty_False.tsv``` posterior draws. For each segment, GRITIC derives the route probabilities for the ```_penalty_True.tsv``` draws by multiplying each route's ordinary probability by ```exp(-2.7 * Average_N_Events)``` exactly once and renormalizing across routes. This post-hoc calculation does not modify the gain timing table. Tree output and downstream mutation timing therefore continue to use the unpenalized route probabilities. See [the publication](https://aacrjournals.org/cancerdiscovery/article/14/10/1810/748591/The-History-of-Chromosomal-Instability-in-Genome) for details of the penalty.

Please see the supplementary materials of the accompanying preprint for more details on the route densities.

### _wgd_calling_info.tsv
A table with a single row that gives WGD calling information for the sample. 

In automatic mode, GRITIC first calculates the modal ```Major_CN``` weighted by half-open segment width (```Segment_End - Segment_Start```) across the configured numbered autosomes; sex chromosomes do not contribute. If no configured autosomal segment is present, inference fails. If this ```Major_CN_Mode``` is 1, the inferred WGD count is 0. Modes other than 1 or 2 are rejected because the corresponding histories are not currently supported.

If the ```Major_CN_Mode``` is 2, GRITIC times major-copy-number-two segments with at least 10 retained mutations on the configured numbered autosomes. For each successfully timed segment it forms a 90% posterior credible interval and finds the mutation-time point covered by the greatest total segment width. ```Overlap_Proportion``` is that covered width divided by the total width of all eligible segments having finite timing intervals. An overlap of at least 60% gives an inferred WGD count of 1; a lower overlap gives 0. Thus 0.6 means that intervals for at least 60% of the eligible genomic span share one time point—it is not a posterior probability or a fraction of mutations or segments. The best overlap timing and proportion are recorded in ```Best_Overlap_Timing``` and ```Overlap_Proportion```. For a WGD call, the overlapping segments are pooled and refit by minor-copy-number class, and their timing densities are combined into the WGD timing distribution recorded in ```Timing```, ```CI_High```, and ```CI_Low```. If no segment is eligible or none yields a finite interval, GRITIC stops with an explicit error.

Supplying ```--wgd-count 0``` bypasses WGD timing and forces non-WGD treatment. Supplying ```--wgd-count 1``` forces WGD treatment but still requires GRITIC to estimate a WGD timing from eligible major-copy-number-two segments. Omitting the option uses the inference above. GRITIC warns when the supplied count conflicts with modal major copy number.

### _gain_timing_table_wgd_segments.tsv
This table is produced while evaluating WGD timing. It gives the timing of eligible major copy number two segments. The ```Best_Overlap_Timing``` column records the point in mutation time that overlaps with the 90% posterior credible intervals of the greatest total width of major-copy-number-two segments. The overlap proportion is recorded in ```Overlap_Proportion```, and whether each segment's timing interval intersects the selected point is recorded in ```Intersecting```.


### _mutation_table.tsv
The mutation table processed by GRITIC to give additional SNV multiplicity information. ```Segment_ID``` is GRITIC's final coordinate-derived segment ID, after copy-number merging when explicitly enabled.

Every output row contains these mutation identity, mapping, and provenance columns:

- ```Source_Segment_ID``` is the input segment ID when both input tables supplied matching ```Segment_ID``` columns. For position-based copy-number assignment it is the final assigned segment ID.
- ```Mutation_ID``` contains the literal input ```Mutation_ID``` when supplied and is blank otherwise.
- ```Position``` contains the canonical non-negative integer input position when supplied and is blank otherwise. The column is always present.
- ```GRITIC_Mutation_ID``` is the canonical sample-unique identifier derived from the source segment plus ```Mutation_ID``` when supplied, otherwise from the source segment plus ```Position```. Its two components are URL-escaped and separated by ```:``` (for example, source segment ```001``` and selected value ```0007``` produce ```001:0007```). Consumers should treat this value as opaque.
- ```Segment_Mutation_Index``` is a zero-based, consecutive index within the final ```Segment_ID```. GRITIC assigns it by sorting ```GRITIC_Mutation_ID``` lexicographically within each segment, so it is independent of input row order. MUTIC uses this index as the row number in the segment's posterior timing matrix.

```Mutation_ID``` and ```Position``` are retained only as provenance after ```GRITIC_Mutation_ID``` is generated. Neither is a downstream mutation key; downstream mutation referencing always uses ```GRITIC_Mutation_ID```.

SNV multiplicity probabilities are given by the ```Prob_Mult_``` columns. ```Alt_Count_Correction_Mult_``` and ```Alt_Count_Correction_Subclone_``` columns estimate the probability that a mutation at the modelled segment coverage and multiplicity would satisfy the configured minimum alternate-read threshold. GRITIC uses these probabilities to correct for alternate-read ascertainment; the separate minimum-total-coverage filter is not included in this correction.

### _tree_plots
Binary tree plots for the gain timings for each route for a given segment. Each plot has two binary trees corresponding to each parental allele. Blue nodes represent independent gains, yellow nodes WGD gains and red nodes the final alleles that were present at time of sampling.

### _timing_dicts
Python dictionary objects containing the stored gain timing and multiplicity posterior samples for each route for each gained segment. Not necessary for most use cases. Can be read using the ```pickle``` and ```bz2``` modules in python. For example:
```
import pickle
import bz2

with bz2.BZ2File('SAMPLE_ID_timing_dicts/1-0-200_timing_dict.bz2', 'rb') as f:
    timing_dict = pickle.load(f)

```

The keys of each dict correspond to the routes for the sample. Within each route there are ```Timing``` and ```Mult``` keys. The ```Timing``` entry gives the timing samples of the independent gains indexed by the corresponding node in the tree. A WGD timing entry is also given if applicable.

The ```Mult``` entry gives the multiplicity proportions corresponding to each timing sample. It is a N_SamplesxN_Multiplicities numpy array. Across the columns, the multiplicities are orderred from 1 to the major copy number of the segment, followed by the subclonal multiplicity probabilities.
