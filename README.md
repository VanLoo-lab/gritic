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

### Required run arguments

- ```--mutation-table``` A path to the SNV table for the sample.
- ```--copy-number-table``` A path to the copy number table for the sample.
- ```--purity``` The estimated cellular purity for the sample. It must be finite, greater than 0, and no greater than 1.
- ```--sample-id``` Sample ID used as an output-directory component and filename prefix. It must be a cross-platform-safe filename component.
- ```--output``` The directory for the output, GRITIC will then store all of the output at --output/SAMPLE_ID

There are also a number of optional arguments.

All probability-, proportion-, quantile-, and interval-width inputs use the
unit interval from 0 to 1; parameters whose semantics exclude zero state that
restriction explicitly.

### Genome and input handling

- ```--autosome-count``` The number of numbered autosomes in the organism. This defines the accepted numbered chromosome labels and the chromosomes eligible for WGD inference. The default is 22.
- ```--sample-sex``` Override the sample sex with ```XX```, ```XY```, ```ZZ```, or ```ZW```. If omitted, GRITIC infers the sex-chromosome system and karyotype from exact X, Y, Z, and W copy-number chromosome labels.
- ```--drop-unmatched-chromosomes``` Drop copy-number and mutation rows whose chromosome is not one of the configured autosomes or present sex chromosomes, with warnings reporting the number of rows dropped. By default, any such chromosome is an error.
- ```--drop-unmatched-snvs``` Drop mutation rows that cannot be associated with a copy-number segment by either supplied ```Segment_ID``` or genomic ```Position```, with one warning reporting the number dropped. By default, unmatched mutations raise an error.
- ```--drop-unrecognized-phasing``` Drop mutation rows whose non-missing ```Phasing``` value is not ```major``` or ```minor```, with one warning reporting the number dropped. By default, unrecognized phasing labels raise an error.
- ```--no-merge-adjacent-segments``` Preserve adjacent copy-number segments as separate intervals. By default, GRITIC merges coordinate-adjacent segments having identical ```Major_CN``` and ```Minor_CN```.

### Mutation filtering and detection correction

- ```--min-mutation-alt-count``` Minimum ```Tumor_Alt_Count``` needed to retain a mutation. The default is 3.
- ```--min-mutation-coverage``` Minimum ```Tumor_Ref_Count + Tumor_Alt_Count``` needed to retain a mutation. The default is 10.
- ```--coverage-vaf-quantile``` Observed-SNV VAF quantile used to select mutations for the mean-coverage estimate in the exact Poisson-thinning detection correction. It is expressed on the unit interval and defaults to 0.9, the unit-interval form of the historical 90th-percentile code value. This replaces ```--coverage-vaf-percentile```; the former option and percentage-scaled values are not accepted. This high-VAF coverage heuristic differs from the publication/thesis method, which averages detection power over every SNV's observed depth in the segment. It affects the likelihood correction in both subclone-fraction-prior modes and the prior adjustment only in ```adjusted``` mode.

### Subclone handling

- ```--subclone-table``` A path to the subclone table for the sample. If not provided it is assumed every SNV is clonal.
- ```--clip-subclone-ccf``` Clip finite ```Subclone_CCF``` values outside the interval from 0 to 1 to the nearest boundary before validation and filtering. This is disabled by default; non-numeric and non-finite values remain errors. With the default CCF filters, values clipped to either boundary are subsequently excluded.
- ```--min-subclone-ccf``` Minimum ```Subclone_CCF``` retained as a subclone, inclusive. It must be greater than 0 and no greater than ```--max-subclone-ccf```. The default is 0.01.
- ```--max-subclone-ccf``` Maximum ```Subclone_CCF``` retained as a subclone, inclusive. It must be no lower than ```--min-subclone-ccf```. The default is 0.9.
- ```--min-subclone-fraction``` A subclone's normalized share of the surviving subclonal mutation fractions must be strictly greater than this threshold. The default is 0.1.
- ```--subclone-fraction-prior {adjusted,supplied}``` Select how the supplied fractions of called/input SNVs assigned to each subclone are used in the mutation-share prior. ```adjusted``` (the default) divides the supplied clonal and subclonal mutation shares by their estimated detection probabilities separately for each copy-number segment and renormalizes them. ```supplied``` uses the fractions directly and reproduces the publication/thesis sample-wide prior. This option does not modify ```Subclone_CCF```, and the likelihood detection correction remains active per segment in both modes.

With the ```adjusted``` prior, clone-share draws use a separately corrected prior for each copy-number segment. With the ```supplied``` publication/thesis prior, they use the supplied sample-wide fractions directly. The likelihood's detection correction is calculated per segment in both modes. Route geometry is independent of either prior and is shared in memory across bounded batches of segments with the same copy-number and WGD state.

Posterior intervals use the shortest contiguous empirical HPD interval by default. Each interval family has a ```--...-interval-width``` option expressed as a proportion greater than 0 and at most 1, and a matching ```--...-interval-method {hpd,equal-tailed}``` option. Percentage-scaled widths such as 90 or 95 are rejected rather than converted. Here HPD means the narrowest contiguous interval containing at least the requested fraction of the empirical draws. The output schemas contain one low/high pair, so multimodal posteriors are represented by one contiguous interval rather than a disjoint highest-density set.

### Inference model and WGD calling

- ```--wgd-count {0,1}``` Override GRITIC's inferred whole-genome-duplication count. Counts below 0 or above 1 are rejected because GRITIC currently supports at most one WGD. If omitted, GRITIC infers the count.
- ```--random-seed``` Seed both NumPy and Numba stochastic inference with an integer from 0 through ```2**32 - 1```. If omitted, GRITIC does not impose a seed.
- ```--wgd-overlap-interval-width```, ```--wgd-overlap-interval-method {hpd,equal-tailed}``` control the hidden candidate-segment bounds used for WGD overlap and can therefore change WGD inference (default width: 0.9, or 90%).
- ```--unordered-balanced-route-prior``` Reproduce the former uniform prior over unordered allele-route pairs while retaining the current ordered route model. In a balanced ```N+N``` segment, an ordered route with different Major and Minor component histories receives half the prior weight of a route with identical component histories. Reciprocal orientations remain separate, so together they carry one prior unit. This is disabled by default.

### Reported timing intervals

- ```--route-gain-interval-width```, ```--route-gain-interval-method {hpd,equal-tailed}``` control the route-conditional gain bounds in both gain timing tables (default width: 0.95, or 95%).
- ```--wgd-timing-interval-width```, ```--wgd-timing-interval-method {hpd,equal-tailed}``` control the final sample-level WGD bounds reused by the JSON, route table, and yellow tree nodes (default width: 0.9, or 90%).
- ```--posterior-summary-interval-width```, ```--posterior-summary-interval-method {hpd,equal-tailed}``` control both gain and gain-conditioned WGD bounds in posterior summaries (default width: 0.95, or 95%).

### Tree plots

- ```--plot-trees``` Plot the route trees for each segment. This is an opt-in switch and is disabled by default.
- ```--tree-gain-interval-width```, ```--tree-gain-interval-method {hpd,equal-tailed}``` control blue gain-node labels in tree PDFs (default width: 0.9, or 90%).

Long option names use hyphens, not underscores. ```Sample_ID``` values are validated rather than sanitized: they cannot be empty, path components such as ```.``` or ```..```, contain path separators, Windows-forbidden filename characters, Unicode control/format/surrogate characters, end in a dot or space, use a reserved Windows device name, or make a derived GRITIC filename exceed the usual 255-unit component limits.


A command to run GRITIC using the example data is:

```
gritic --mutation-table examples/snv_table_example.tsv --copy-number-table examples/cn_table_example.tsv --purity 0.5 --sample-id TEST_ID --output examples/output --subclone-table examples/subclone_table_example.tsv --wgd-count 1
```

GRITIC can also be run programmatically:

```python
import pandas as pd

from gritic import dataloader, gritictimer, intervaltools, sampletools

copy_number_table, mutation_table = dataloader.load_input_tables(
    'examples/cn_table_example.tsv',
    'examples/snv_table_example.tsv',
    drop_unrecognized_phasing=False,
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
    merge_cn=True,
    min_mutation_alt_count=3,
    min_mutation_coverage=10,
    coverage_vaf_quantile=0.9,
    min_subclone_ccf=0.01,
    max_subclone_ccf=0.9,
    min_subclone_fraction=0.1,
    clip_subclone_ccf=False,
    autosome_count=22,
    drop_unmatched_chromosomes=False,
    drop_unrecognized_phasing=False,
)
gritictimer.process_sample(
    sample,
    output_dir='examples/output',
    plot_trees=True,
    min_wgd_overlap=0.6,
    wgd_count=1,
    subclone_fraction_prior='adjusted',
    unordered_balanced_route_prior=False,
    random_seed=20260828,
    interval_config=intervaltools.TimingIntervalConfig(
        route_gain=intervaltools.IntervalSpec(0.95),
        tree_gain=intervaltools.IntervalSpec(0.9),
        wgd_overlap=intervaltools.IntervalSpec(0.9),
        sample_wgd=intervaltools.IntervalSpec(0.9),
        posterior_summary=intervaltools.IntervalSpec(0.95),
    ),
)
```

## Input Table Formats
The required copy-number and mutation tables, and the optional subclone table, should be tab separated. Examples using simulated data are available in the example directory. Data from any allele specific copy number caller, SNV caller and subclone caller can be used as long as the tables are formatted correctly.

Column names expected or generated by GRITIC use the exact, case-sensitive ```Pascal_Snake_Case``` spellings documented below. Additional caller-supplied columns are accepted without a naming-style check and are not renamed by GRITIC.

### Mutation Table 
All SNVs for the sample. The columns ```Chromosome```, ```Tumor_Ref_Count``` and ```Tumor_Alt_Count``` are always required. A chromosome must be a canonical integer from 1 through ```--autosome-count``` or a sex chromosome present under the selected/inferred karyotype; one leading ```chr``` prefix is accepted and removed. Thus XX accepts X, XY accepts X/Y, ZZ accepts Z, and ZW accepts Z/W. Both read counts must be finite non-negative integers within the signed 64-bit range; their sum must be greater than zero and must also fit within that range. Integer-equivalent spellings are canonicalized to signed-64-bit integers. Every mutation table must also contain either ```Mutation_ID``` or ```Position```. ```Mutation_ID``` is loaded as literal text, so values such as ```000123```, ```NA```, and ```NULL``` are preserved. When present, every ```Position``` must instead be a finite non-negative integer within the signed 64-bit range; integer-equivalent spellings such as ```000123``` or ```123.0``` are canonicalized to ```123```.

```Phasing``` is optional. Non-missing phasing labels are stripped of surrounding whitespace and converted to lowercase, so values such as ```Major ```, ```MAJOR```, and ``` major ``` all become ```major```. The only recognized labels are ```major``` and ```minor```; missing values remain unphased. By default, any other non-missing label raises an error. Supplying ```--drop-unrecognized-phasing``` drops the affected mutation rows instead and emits one warning with the number dropped and the unrecognized values.

After copy-number assignment, every mutation with ```Phasing=minor``` must have an assigned ```Minor_CN``` greater than zero. GRITIC raises an error rather than dropping the mutation or coercing its phasing, whether copy number was assigned through supplied ```Segment_ID``` values or by genomic position.

An original balanced ```2+2``` segment has one WGD-specific exception. Only while GRITIC evaluates a possible WGD and builds pooled timing for a called WGD, it retains each mutation's ```Phasing``` label but evaluates every SNV against one representative duplicated allele because the two homologs' routes and timings are constrained to be identical. If WGD is rejected, GRITIC uses the ordinary independent-route ```2+2``` model and phasing remains active.

The ```Phasing=minor``` invariant still accepts this WGD-specific case: the representative allele's internal geometry is ```2+0```, but the mutation table's emitted, original ```Minor_CN``` remains ```2```.

Input mutation columns named ```Segment_Start```, ```Segment_End```, ```Major_CN```, or ```Minor_CN``` are ignored because GRITIC always annotates those values from the copy-number table.

GRITIC supports two segment-assignment modes:

- If both the mutation and copy number tables contain ```Segment_ID```, GRITIC uses those IDs to associate mutations with the input copy number segments. ```Position``` is not needed for assignment.
- Otherwise, the mutation table must contain ```Position```, which GRITIC uses to associate mutations with copy-number segments under the zero-based, half-open interval rule ```Segment_Start <= Position < Segment_End```.

In supplied-ID mode, every copy number row must have a ```Segment_ID``` that is not null, empty, or whitespace-only, and copy number ```Segment_ID``` values must be unique within the copy number table. Every mutation row must also have a ```Segment_ID``` that is not null, empty, or whitespace-only. By default, every mutation ```Segment_ID``` must be present in the copy number table, and ```Chromosome``` must match the corresponding copy number row. In position-assignment mode, every mutation must fall inside a copy-number segment on the same chromosome. If ```--drop-unmatched-snvs``` is supplied, either kind of unmatched mutation is instead dropped with one count warning. Missing or blank supplied IDs and chromosome mismatches remain errors.

GRITIC does not model mutations in zero-copy ```0+0``` segments. In supplied-ID
mode, mutations assigned to such otherwise valid input segments are dropped with
one aggregated warning reporting the affected mutation count, segment count, and
source ```Segment_ID``` values. This is independent of
```--drop-unmatched-snvs``` because those IDs are matched rather than unknown.
If no mutations remain after this removal, GRITIC raises an error.

The selected ```Mutation_ID``` or canonical integer ```Position``` value must be unique within its source segment. An explicit ```Mutation_ID``` column takes precedence for every row; GRITIC does not fall back to ```Position``` for blank values in that column. The same selected value may be reused in different source segments, including source segments that GRITIC subsequently merges.

Programmatic callers enable the same behavior with ```drop_unmatched_snvs=True``` in ```sampletools.Sample```. When using ```dataloader.load_input_tables``` with supplied segment IDs, pass the same value there as well so unmatched IDs are permitted and removed during loading; the command-line interface forwards it to both stages. Likewise, ```drop_unrecognized_phasing=True``` corresponds to ```--drop-unrecognized-phasing``` in both APIs. The separate ```drop_unmatched_chromosomes=True``` API argument corresponds to ```--drop-unmatched-chromosomes``` and filters invalid chromosome rows from both input tables before their numeric fields and segment relationships are validated.

### Copy Number Table 
The rounded allele-specific copy number profile for the sample. Requires the column names ```Chromosome```, ```Segment_Start```, ```Segment_End```, ```Major_CN``` & ```Minor_CN```. Chromosome labels follow the same rules as the mutation table. Segment coordinates are zero-based, half-open, non-negative signed-64-bit integers: ```Segment_Start``` is included, ```Segment_End``` is excluded, and ```Segment_End``` must be strictly greater than ```Segment_Start```. A one-base segment is therefore ```[a, a + 1)``` and has width ```Segment_End - Segment_Start```. Both allele-specific copy numbers must also be finite signed-64-bit integers and non-negative, and ```Major_CN``` must be greater than or equal to ```Minor_CN```. Violations raise an error. ```Segment_ID``` is optional and it is only used when it is also present in the mutation table. Mixed sex-chromosome systems make sex inference ambiguous, and numbered chromosomes above ```--autosome-count``` are unmatched.

GRITIC orders copy-number rows by natural chromosome order, ```Segment_Start```, and ```Segment_End``` immediately after validating the input coordinates. Half-open intervals on the same chromosome must not overlap. By default GRITIC merges equal-copy-number intervals only when the following ```Segment_Start``` equals the preceding ```Segment_End```. Gapped intervals are not merged, and the resulting merged interval is revalidated against the same coordinate bounds. Supplying ```--no-merge-adjacent-segments``` (or ```merge_cn=False``` through the API) instead preserves each input copy-number segment and generates final segment IDs from ```Chromosome```, ```Segment_Start```, and ```Segment_End```.

When sample sex is not supplied, an exact Y chromosome label implies ```XY``` and W implies ```ZW```; if neither is present, X implies ```XX``` and Z implies ```ZZ```. If no sex chromosome is represented, GRITIC defaults to ```XX```, so callers using another system should supply ```--sample-sex```. Inputs that mix the X/Y and Z/W systems cannot be inferred and fail. An explicit karyotype takes precedence. For ```XX```, X is present and Y is unmatched; for ```XY```, both X and Y are present. For ```ZZ```, Z is present and W is unmatched; for ```ZW```, both Z and W are present. Unmatched rows fail unless chromosome dropping is enabled. Normal-cell X/Y copy numbers are 2/0 for XX and 1/1 for XY; the corresponding Z/W values are 2/0 for ZZ and 1/1 for ZW. The ```autosome_count``` API argument and ```--autosome-count``` CLI option define the numbered autosomes and default to 22.
### Subclone Table (*Optional*)
The identified subclonal peaks and their assigned mutation fractions for the sample. GRITIC filters candidate peaks using configurable lower and upper CCF bounds and the minimum normalized-fraction threshold described below.

This table requires the column names ```Subclone_CCF``` (the inferred fraction of cancer cells containing the subclone) & ```Subclone_Fraction``` (the fraction of called/input SNVs assigned to the subclone). ```Subclone_Fraction``` is not a cellular fraction and does not estimate mutations absent from the input call set. A ```Cluster``` column is also required as an index for the subclones. CCF and fraction values must be finite numeric (not boolean) values between 0 and 1 inclusive, and the fractions must sum to no more than 1. By default, invalid prior values fail before filtering. Supplying ```--clip-subclone-ccf``` (or ```clip_subclone_ccf=True``` through the API) instead clips finite out-of-range CCF values to 0 or 1 before the configured CCF filters run; it does not relax validation of non-finite CCFs or any ```Subclone_Fraction``` value.

By default, only subclones with ```0.01 <= Subclone_CCF <= 0.9``` and more than 10% of the surviving subclonal mutation fraction are included. Both CCF boundaries are inclusive. The thresholds are configurable with ```--min-subclone-ccf```, ```--max-subclone-ccf```, and ```--min-subclone-fraction```; the minimum CCF must be greater than 0 and cannot exceed the maximum. If no subclones remain, GRITIC uses its clonal-only model as if no subclone table was supplied. If there are more than two subclones, GRITIC reformats the sample to have two: the subclone with the largest CCF is unmodified and the remaining clones are grouped together.

GRITIC always derives ```N_SNVs``` from the retained mutation count and each retained/combined ```Subclone_Fraction```, overwriting an input ```N_SNVs``` column. ```Subclone_Fraction``` contributes to the clonal/subclonal mixture prior used during route sampling, while ```Subclone_CCF``` determines the expected VAF of each subclone state.

This table is optional, if it is not included GRITIC will assume every SNV is clonal. This will likely bias the gains to be measured earlier.


## Output
GRITIC produces a number of outputs to describe the timing of copy number gains in a given tumour. We recommend only considering gained segments with 10 or more SNVs. 

Each run requires its ```OUTPUT/SAMPLE_ID``` directory to be absent or empty. GRITIC rejects a nonempty sample directory rather than appending to, or mixing output with, an earlier run.

### _posterior_timing_table_summary_penalty_&lt;True|False&gt;.tsv (Main Output)
Gives a summary of the gain timing information for each gained segment. For each gained segment, the sequential gains are labelled with the 1-based ```Gain_Index```, where 1 is the earliest sampled independent gain. The median and the configured posterior-summary interval are reported; this is a 95% HPD interval by default.

In WGD tumours, the number of gains that arise independently of the WGD will vary depending on the route. Only gains that arise independently of the WGD in at least 80% of posterior route draws are reported; this is recorded in the ```Proportion``` column. The denominator includes every internal route draw, including draws of routes with no independent gains, rather than being inferred from gain rows alone. Gain and WGD intervals on a row are both calculated from the subset of route-marginalized draws in which that kth gain exists. The timing of the WGD and the probability that each gain arose before the WGD are also recorded.

Two summary tables are produced for every run that produces timing output, including non-WGD runs. The filename ending in ```_penalty_False.tsv``` summarizes draws using ```Probability``` from the route table. The filename ending in ```_penalty_True.tsv``` summarizes a second set of draws using ```Penalized_Probability```. No command-line option is required to produce either result.

The gain-draw and route-ledger data frames used to calculate these summaries are internal. GRITIC no longer writes the legacy ```_posterior_timing_table_penalty_<True|False>.tsv``` or ```_posterior_route_draw_table_penalty_<True|False>.tsv``` files. Raw route-conditional timing samples remain available in the timing stores.

### _route_table.tsv
This table contains exactly one row for each possible route of each timed segment. Its first three columns and key are ```(Sample_ID, Segment_ID, Route)```. It stores both the ordinary ```Probability``` and post-hoc ```Penalized_Probability```, average event and loss counts, timing-geometry proposal density and runtime, and the segment and WGD metadata. The density diagnostic covers the hit-and-run timing coordinates; independently sampled Dirichlet clone shares are not chain dimensions. Routes with no independently timeable gains are retained here. Its ```WGD_Timing``` interval fields exactly repeat the configured sample-level WGD summary written to the calling-info JSON.

```Route``` is an opaque, order-sensitive identity of the complete allele route: its internal identity is formed from the ```Major``` component followed by the ```Minor``` component, including their WGD annotations. An absent minor allele has one explicit empty-component identity. GRITIC does not recover allele identity from component size and does not sort the two component hashes. Consequently, all route identifiers deliberately differ from identifiers produced by the former unordered-hash implementation. A route table, gain timing table, timing store, and tree plots from an older run must not be mixed with current output; regenerate the complete sample output together.

For an allele-balanced ```N+N``` state, GRITIC enumerates the full ordered Cartesian product of the single-allele histories. If histories ```A``` and ```B``` differ, both ```Major=A, Minor=B``` and ```Major=B, Minor=A``` are separate routes with separate identifiers; a diagonal ```A/A``` route occurs once. By default, each ordered route receives one equal prior unit before its marginal likelihood is normalized, equivalent to independently drawing the two allele histories from a uniform single-allele route prior while conditioning on the shared copy-number and WGD setting.

With ```--unordered-balanced-route-prior```, a diagonal route retains prior weight 1 and each orientation of an off-diagonal pair receives weight 0.5 before evidence normalization. Thus ```A/B``` and ```B/A``` together carry the same one prior unit as ```A/A```, exactly reproducing the former uniform prior over unordered route pairs after expansion into the ordered state space. This option changes only route-level prior weights: it does not merge reciprocal routes, alter route-conditional timing samples, or restore the former one-sided handling of phased balanced routes. It is a no-op for unbalanced states and balanced states having only an identical-component route. In-memory geometry sharing is unchanged because route geometry is independent of route probabilities.

When a balanced segment contains no phased SNVs, the two orientations of an off-diagonal pair have exactly the same likelihood. GRITIC fits one orientation and constructs the other by explicitly exchanging its Major and Minor timing and multiplicity arrays; both ordered routes remain in the output and receive either one prior unit each by default or half a unit each under ```--unordered-balanced-route-prior```. With any phased SNV, GRITIC evaluates the likelihood and likelihood-weighted posterior separately for both orientations before applying the selected route prior. The orientations share one likelihood-independent in-memory geometry sample, with the reciprocal route produced by explicitly transforming those proposal arrays; segment-specific clone shares are then applied in memory and no temporary gzip proposal cache is written. WGD ```2+2``` timing is the deliberate exception: all SNVs are first treated as one representative allele and the model uses the single pseudo-```2+0``` route because the two homolog histories and timings are constrained to be identical.

```Penalized_Probability``` is calculated by multiplying each route's ordinary probability by ```exp(-2.7 * Average_N_Events)``` exactly once and renormalizing across routes within the segment. The original ```Probability``` remains unchanged, and each probability column sums to one within a valid segment. The number of independent gains does not alter either route-selection probability. Tree output and downstream mutation timing continue to use the ordinary probabilities. See [the publication](https://aacrjournals.org/cancerdiscovery/article/14/10/1810/748591/The-History-of-Chromosomal-Instability-in-Genome) for details of the penalty.

### _gain_timing_table.tsv
This table contains exactly one row per independently timeable gain node per route, keyed by ```(Sample_ID, Segment_ID, Route, Node)```. ```Node``` is written as an integer structural tree identifier; it is not a chronological gain index and need not be consecutive. The table stores ```Node_Phasing``` and the median and configured route-gain interval of that node's 1,000 likelihood-weighted samples conditional on the route. The default is a 95% HPD interval. Route probabilities and other route-level fields are not duplicated here. A route without an independently timeable gain has no row in this table; no artificial ```NA``` node is written.

```Node_Phasing``` is a structural route-component label and always uses the exact, case-sensitive value ```Major``` or ```Minor```. Every node in a route with only one extant allele component is labelled ```Major```. For a route with two components, each node retains the explicit allele role assigned when that ordered route was constructed; GRITIC never recovers this role from component size. In a balanced state with different component histories, the reciprocal role assignment is represented by a separate route. The two names are biologically exchangeable in the absence of phased evidence, in which case the reciprocal routes receive equal probabilities. These labels are separate from mutation-table ```Phasing```, whose canonical values remain lowercase ```major``` and ```minor``` (or missing for an unphased mutation).

When constructing posterior-summary draws, GRITIC first selects a route from the unique rows in ```_route_table.tsv```. It then selects one posterior-array index for that route and uses the same index for the WGD and every independent-gain timing, preserving their joint dependence. The gain timing table supplies node identity and phasing metadata; the full aligned posterior arrays remain in ```_timing_dicts```.

Please see the supplementary materials of the accompanying preprint for more details on the route densities.

### _wgd_calling_info.json
A JSON object that gives WGD calling information for the sample. Missing or nonfinite numeric values are represented by JSON ```null```. Its keys are ```WGD_Timing```, ```WGD_Timing_CI_Low```, ```WGD_Timing_CI_High```, ```Major_CN_Mode```, ```Overlap_Proportion```, ```WGD_Status```, and ```Best_Overlap_Timing```.

In automatic mode, GRITIC first calculates the modal ```Major_CN``` weighted by half-open segment width (```Segment_End - Segment_Start```) across the configured numbered autosomes; sex chromosomes do not contribute. If no configured autosomal segment is present, inference fails. If this ```Major_CN_Mode``` is 1, the inferred WGD count is 0. Modes other than 1 or 2 are rejected because the corresponding histories are not currently supported.

If the ```Major_CN_Mode``` is 2, GRITIC times major-copy-number-two segments with at least 10 retained mutations on the configured numbered autosomes. For each successfully timed segment it forms the configured hidden WGD-overlap interval and finds the mutation-time point covered by the greatest total segment width. The default is a 90% HPD interval. ```Overlap_Proportion``` is that covered width divided by the total width of all eligible segments having finite timing intervals. An overlap of at least 60% gives an inferred WGD count of 1; a lower overlap gives 0. Thus 0.6 means that intervals for at least 60% of the eligible genomic span share one time point—it is not a posterior probability or a fraction of mutations or segments. The best overlap timing and proportion are recorded in ```Best_Overlap_Timing``` and ```Overlap_Proportion```. For a WGD call, the overlapping segments are pooled and refit by minor-copy-number class, and their timing densities are combined into 500 sample-level WGD draws. ```WGD_Timing``` and its interval fields summarize those draws using the configured WGD-timing interval, a 90% HPD interval by default. If no segment is eligible or none yields a finite interval, GRITIC stops with an explicit error.

Supplying ```--wgd-count 0``` bypasses WGD timing and forces non-WGD treatment. Supplying ```--wgd-count 1``` forces WGD treatment but still requires GRITIC to estimate a WGD timing from eligible major-copy-number-two segments. Omitting the option uses the inference above. GRITIC warns when the supplied count conflicts with modal major copy number.

### _gain_timing_table_wgd_segments.tsv
This table is produced while evaluating WGD timing. It gives the preliminary non-WGD node posterior for every eligible major-copy-number-two segment. Its displayed ```Timing_CI_Low``` and ```Timing_CI_High``` use the same configured route-gain interval as the main gain timing table, 95% HPD by default. ```Intersecting```, ```Best_Overlap_Timing```, and ```Overlap_Proportion``` are instead calculated from the separate, unrounded hidden WGD-overlap intervals, 90% HPD by default. The visible 95% bounds therefore need not be the bounds that determine ```Intersecting```.


### _mutation_table.tsv
The mutation table processed by GRITIC to give additional SNV multiplicity information. ```Segment_ID``` is GRITIC's final coordinate-derived segment ID, after copy-number merging when explicitly enabled.

```Phasing``` contains the canonical lowercase values ```major``` or ```minor```, or is blank for an unphased mutation.

Every output row contains these mutation identity, mapping, and provenance columns:

- ```Source_Segment_ID``` is the input segment ID when both input tables supplied matching ```Segment_ID``` columns. For position-based copy-number assignment it is the final assigned segment ID.
- ```Mutation_ID``` contains the literal input ```Mutation_ID``` when supplied and is blank otherwise.
- ```Position``` contains the canonical non-negative integer input position when supplied and is blank otherwise. The column is always present.
- ```GRITIC_Mutation_ID``` is the canonical sample-unique identifier derived from the source segment plus ```Mutation_ID``` when supplied, otherwise from the source segment plus ```Position```. Its two components are URL-escaped and separated by ```:``` (for example, source segment ```001``` and selected value ```0007``` produce ```001:0007```). Consumers should treat this value as opaque.
- ```Segment_Mutation_Index``` is a zero-based, consecutive index within the final ```Segment_ID```. GRITIC assigns it by sorting ```GRITIC_Mutation_ID``` lexicographically within each segment, so it is independent of input row order. MUTIC uses this index as the row number in the segment's posterior timing matrix.

```Mutation_ID``` and ```Position``` are retained only as provenance after ```GRITIC_Mutation_ID``` is generated. Neither is a downstream mutation key; downstream mutation referencing always uses ```GRITIC_Mutation_ID```.

SNV multiplicity probabilities are given by the ```Prob_Mult_``` columns. ```Alt_Count_Correction_Mult_``` and ```Alt_Count_Correction_Subclone_``` columns contain the exact probability that a Poisson alternate-read count, whose mean is the selected segment coverage multiplied by the modelled VAF, satisfies the configured minimum alternate-read threshold. GRITIC uses these probabilities to correct for alternate-read ascertainment; the separate minimum-total-coverage filter is not included in this correction.

### _subclone_table.tsv
The retained and normalized subclone inputs with the fixed columns ```Cluster```, ```Subclone_CCF```, ```Subclone_Fraction```, and ```N_SNVs```. This file is always written. If no subclone table was supplied, or no candidate survived filtering, it contains the header and zero rows.

### _tree_plots
Binary tree plots for the gain timings for each route for a given segment. Each plot has two binary trees corresponding to each parental allele. Blue nodes represent independent gains and show the configured tree-gain interval from the same 1,000 route-conditional samples, 90% HPD by default. Yellow WGD nodes exactly repeat the configured final sample-level WGD interval used in the JSON and route table. Red nodes are the final alleles present at sampling.

### _timing_dicts
Compressed timing stores containing the gain-timing and multiplicity posterior samples for every route of each gained segment. These stores are not necessary for most use cases. Each logical store is a required pair:

- ```SEGMENT_ID_timing_dict.npz``` contains the numeric tables in linear order as ```table_000000```, ```table_000001```, and so on. It is a compressed NumPy archive and never contains pickled or object-dtype arrays.
- ```SEGMENT_ID_timing_dict.manifest.json``` maps the original nested hierarchy onto those table indexes. It also records the format version, archive filename and SHA-256, and each table's dtype and shape.

Pooled WGD stores use the same pair with a ```WGD_minor_cn_N``` identifier in place of ```SEGMENT_ID```. Both members of a pair are required. Load a store through GRITIC so the pair and manifest are validated and the original dictionary hierarchy is reconstructed:

```
from gritic.timingio import load_timing_archive

timing_dict = load_timing_archive(
    'SAMPLE_ID_timing_dicts/1-0-200_timing_dict.npz'
)

```

The reconstructed dictionary keys correspond to the routes for the sample. Within each route there are ```Timing```, ```Mult```, and ```Raw_Samples``` keys. The ```Timing``` entry gives raw stored timing distributions for independent gains indexed by the corresponding structural node in the tree; it does not store a chronological ```Gain_Index```. A WGD timing entry is also given if applicable. No credible interval is applied to timing stores: consumers choose their own interval width and method from the stored samples.

The ```Mult``` entry gives the multiplicity proportions corresponding to each timing sample. It is a N_SamplesxN_Multiplicities numpy array. Across the columns, the multiplicities are orderred from 1 to the major copy number of the segment, followed by the subclonal multiplicity probabilities.

```Raw_Samples``` contains the aligned proposal-sampling ```Timing```, ```Mult```, ```WGD_Timing```, and ```LL``` arrays used before likelihood resampling. The likelihood-resampled route output does not contain an ```LL``` array.
