# GRITIC

GRITIC estimates the timing and route histories of clonal copy number gains from single-sample bulk whole-genome sequencing of cancer genomes. It supports complex allele-specific copy-number states and up to one whole-genome duplication (WGD); see [timing eligibility](#timing-eligibility).

Gain timing is measured in mutation time, from 0 near conception to 1 at the emergence of the tumour's most recent common ancestor.

GRITIC is agnostic to reference genome. The number of numbered autosomes is configurable, and both X/Y and Z/W sex-chromosome systems are supported.

The method is described in [Baker et al. (2024), *The History of Chromosomal Instability in Genome-Doubled Tumors*][publication].

## Installation

Python 3.12 or newer is required. Install the version documented here from this repository checkout:

```bash
python -m pip install .
```

## Quick start

Run the examples from the repository root. Each run writes directly to `--sample-dir`, which must not already exist, even if empty. Use `--overwrite` to reuse it. Both examples below use `examples/output/TEST_ID`.

### Command line

```bash
gritic \
    --mutation-table examples/snv_table_example.tsv \
    --copy-number-table examples/cn_table_example.tsv \
    --subclone-table examples/subclone_table_example.tsv \
    --purity 0.5 \
    --sample-id TEST_ID \
    --sample-dir examples/output/TEST_ID \
    --wgd-count 1 \
    --random-seed 20260828 \
    --plot-trees
```

Use `gritic --help` for command-line help or see [run options](#run-options).

### Python API

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
    dtype={'Cluster': str},
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
    sample_dir='examples/output/TEST_ID',
    plot_trees=True,
    wgd_count=1,
    random_seed=20260828,
)
```

Python callers configure [timing intervals](#timing-intervals) through the `interval_config` argument to `process_sample`, using `distributiontools.TimingIntervalConfig` and `distributiontools.IntervalSpec`.

For unmatched-SNV dropping with supplied segment IDs, `dataloader.load_input_tables(..., drop_unmatched_snvs=True)` removes unmatched rows during loading. For position-based assignment, pass `drop_unmatched_snvs=True` to `sampletools.Sample`.

## Input tables

The required copy-number and mutation tables, and the optional subclone table, should be tab separated. Example tables are available in [`examples/`](examples/). Inputs can come from any allele-specific copy-number caller, SNV caller, or subclone caller.

Column names are case-sensitive. Additional caller-supplied columns are accepted.

### Mutation table

All SNVs for the sample. The columns `Chromosome`, `Tumor_Ref_Count`, and `Tumor_Alt_Count` are required. Chromosome labels follow the [chromosome-handling rules](#chromosome-handling). Both read counts must be non-negative integers with a positive sum. Every mutation table must also contain either `Mutation_ID` or `Position`. `Mutation_ID` is loaded as literal text, preserving values such as `000123`, `NA`, and `NULL`. When present, `Position` must be a non-negative integer.

`Phasing` is optional and accepts `major` and `minor` case-insensitively, ignoring surrounding whitespace. Missing values remain unphased. Supplying `--drop-unrecognized-phasing` drops affected mutation rows and emits one warning with the number dropped and the unrecognized values.

After copy-number assignment, every mutation with `Phasing=minor` must have an assigned `Minor_CN` greater than zero.

Input mutation columns named `Segment_Start`, `Segment_End`, `Major_CN`, or `Minor_CN` are ignored because GRITIC always annotates those values from the copy-number table.

GRITIC supports two segment-assignment modes:

- If both the mutation and copy number tables contain `Segment_ID`, GRITIC uses those IDs to associate mutations with the input copy number segments. `Position` is not needed for assignment.
- Otherwise, the mutation table must contain `Position`, which assigns mutations to segments on the same chromosome using the [coordinate convention](#copy-number-table).

In supplied-ID mode, both tables require nonblank `Segment_ID` values, and those IDs must be unique in the copy-number table. Every mutation's `Segment_ID` must match a copy-number row, including its `Chromosome`. If `--drop-unmatched-snvs` is supplied, mutations unmatched by either assignment mode are dropped with one count warning. Missing or blank supplied IDs and chromosome mismatches remain errors.

GRITIC does not model mutations in zero-copy `0+0` segments. In supplied-ID
mode, mutations assigned to such otherwise valid input segments are dropped with
one aggregated warning reporting the affected mutation count, segment count, and
source `Segment_ID` values. This is independent of
`--drop-unmatched-snvs` because those IDs are matched rather than unknown.

The selected `Mutation_ID` or canonical integer `Position` value must be unique within its source segment. An explicit `Mutation_ID` column takes precedence for every row; GRITIC does not fall back to `Position` for blank values in that column. The same selected value may be reused in different source segments. With supplied segment IDs, this also applies to source segments that GRITIC subsequently merges.

### Copy number table

The rounded allele-specific copy-number profile for the sample requires `Chromosome`, `Segment_Start`, `Segment_End`, `Major_CN`, and `Minor_CN`. Chromosome labels follow the [rules below](#chromosome-handling). Segments use nonempty, zero-based, half-open intervals with non-negative integer coordinates. Intervals on the same chromosome must not overlap. Allele-specific copy numbers are non-negative integers with `Major_CN >= Minor_CN`. Supplied `Segment_ID` values enable the assignment mode described under [Mutation table](#mutation-table).

#### Timing eligibility

The supported allele-specific copy-number states are:

| `Major_CN` | Permitted `Minor_CN` |
| --- | --- |
| 1 | 0–1 |
| 2 | 0–2 |
| 3–5 | 0 through `Major_CN` |
| 6 | 0–4 |
| 7 | 0–3 |
| 8 | 0–1 |

The `1+0` and `1+1` states use uniform no-gain timing. In WGD runs, major-copy-number-two segments are excluded from the ordinary per-segment fit. [WGD timing estimation](#wgd-timing-estimation) uses those on configured autosomes with at least 10 retained SNVs. Balanced `2+2` segments use a [WGD-specific allele model](#balanced-wgd-segments).

#### Segment merging

By default, GRITIC merges consecutive segments on the same chromosome having identical `Major_CN` and `Minor_CN`. The merged interval runs from the first segment's start through the last segment's end, including any intervening uncovered bases. Use `--max-merge-gap` to limit gaps or `--no-merge-adjacent-segments` to preserve input segments; the Python arguments are `max_merge_gap` and `merge_cn=False`. Final segment IDs are generated from `Chromosome`, `Segment_Start`, and `Segment_End` after merging.

#### Chromosome handling

Both input tables accept numbered chromosomes from 1 through `--autosome-count` and the sex chromosomes present under the selected karyotype. Chromosome labels are case-sensitive; one leading `chr` prefix is removed.

| Karyotype | Present sex chromosomes | Normal X or Z copies | Normal Y or W copies |
| --- | --- | --- | --- |
| XX | X | 2 | 0 |
| XY | X, Y | 1 | 1 |
| ZZ | Z | 2 | 0 |
| ZW | Z, W | 1 | 1 |

When sample sex is not supplied, GRITIC infers it from the copy-number table: Y implies `XY` and W implies `ZW`; otherwise X implies `XX` and Z implies `ZZ`. If no sex chromosome is represented, GRITIC defaults to `XX`, so callers using another system should supply `--sample-sex`. Inputs mixing the X/Y and Z/W systems cannot be inferred. Chromosomes outside the configured set follow `--drop-unmatched-chromosomes` behavior.

### Subclone table

The optional subclone table gives the identified subclonal peaks and their assigned mutation fractions for the sample.

Required columns are `Cluster` (the subclone identifier), `Subclone_CCF` (cancer cell fraction), and `Subclone_Fraction` (the fraction of input SNVs assigned to the subclone). CCF determines the expected VAF of a subclone state; mutation shares determine its mixture prior. `Subclone_Fraction` is not a cellular fraction and does not estimate mutations absent from the input call set. Fractions must sum to no more than 1. Values are validated before filtering.

Candidates are filtered using the [subclone-handling options](#subclone-handling). If no subclones remain, GRITIC uses its clonal-only model. If there are more than two subclones, GRITIC groups them into two: the subclone with the largest CCF is unmodified and the remaining clones are combined by summing their fractions and taking their fraction-weighted mean CCF.

Retained `Subclone_Fraction` values keep their input scale; the remaining share is treated as clonal.

GRITIC derives `N_SNVs` from the retained mutation count and each retained/combined `Subclone_Fraction`, overwriting an input `N_SNVs` column.

## Run options

### Required run arguments

- `--mutation-table` A path to the [mutation table](#mutation-table) for the sample.
- `--copy-number-table` A path to the [copy-number table](#copy-number-table) for the sample.
- `--purity` The estimated cellular purity for the sample; must be greater than 0.
- `--sample-id` Sample ID used as an output filename prefix. It must be a cross-platform-safe filename component.
- `--sample-dir` Directory for this sample. GRITIC writes outputs directly here without appending the sample ID. Pass the same path to MUTIC and SIGTIC.

Probability, proportion, quantile, and interval-width inputs use `[0, 1]`; parameters that exclude zero state this explicitly.

### Output handling

- `--overwrite` Reuse the existing `--sample-dir` directory and its subdirectories. Files with the same names as new outputs are overwritten; all other files are preserved. Parent directories may already exist without this switch. The Python API accepts `overwrite=True` for the same behavior.

### Genome and input handling

- `--autosome-count` The number of numbered autosomes in the organism. This defines the accepted numbered chromosome labels and the chromosomes eligible for WGD inference. The default is 22.
- `--sample-sex` Override the inferred karyotype with `XX`, `XY`, `ZZ`, or `ZW`. See [chromosome handling](#chromosome-handling) for inference rules and normal copy numbers.
- `--drop-unmatched-chromosomes` Drop copy-number and mutation rows whose chromosome is not one of the configured autosomes or present sex chromosomes, with warnings reporting the number of rows dropped. By default, any such chromosome is an error.
- `--drop-unmatched-snvs` Drop mutation rows that cannot be associated with a copy-number segment by either supplied `Segment_ID` or genomic `Position`, with one warning reporting the number dropped. By default, unmatched mutations raise an error.
- `--drop-unrecognized-phasing` Drop mutation rows whose non-missing `Phasing` value is not `major` or `minor`, with one warning reporting the number dropped. By default, unrecognized phasing labels raise an error.
- `--no-merge-adjacent-segments` Preserve input copy-number segments separately. See [segment merging](#segment-merging) for the default behavior.
- `--max-merge-gap N` Merge consecutive equal-copy-number segments only when the gap is at most `N` bases. If omitted, there is no maximum; use `0` to merge only intervals that touch.

### Mutation filtering and detection correction

- `--min-mutation-alt-count` Minimum `Tumor_Alt_Count` needed to retain a mutation. The default is 3.
- `--min-mutation-coverage` Minimum `Tumor_Ref_Count + Tumor_Alt_Count` needed to retain a mutation. The default is 10.
- `--coverage-vaf-quantile` Observed-SNV VAF quantile used to select mutations for the mean-coverage estimate in the detection correction. The default is 0.9. See [detection correction and mutation-share priors](#detection-correction-and-mutation-share-priors).

### Subclone handling

- `--subclone-table` A path to the [subclone table](#subclone-table) for the sample. If omitted, GRITIC assumes every SNV is clonal, which can bias gain timings earlier.
- `--clip-subclone-ccf` Clip out-of-range `Subclone_CCF` values before validation and filtering. This is disabled by default. With the default CCF filters, values clipped to either boundary are subsequently excluded.
- `--min-subclone-ccf` Minimum `Subclone_CCF` retained as a subclone, inclusive. The default is 0.01.
- `--max-subclone-ccf` Maximum `Subclone_CCF` retained as a subclone, inclusive. The default is 0.9.
- `--min-subclone-fraction` A subclone's normalized share of the subclonal mutation fractions after CCF filtering must be strictly greater than this threshold. The default is 0.1.
- `--subclone-fraction-prior {adjusted,supplied}` Use detection-adjusted mutation fractions (`adjusted`, the default) or the supplied fractions (`supplied`) in the mutation-share prior. See [detection correction and mutation-share priors](#detection-correction-and-mutation-share-priors).

The CCF bounds must satisfy `0 < min_subclone_ccf <= max_subclone_ccf`.

### Inference model and WGD calling

Every run requires configured autosomal segments with a segment-width-weighted modal `Major_CN` of 1 or 2, including runs with a supplied WGD count.

- `--wgd-count {0,1}` Override GRITIC's inferred WGD count. A count of 0 bypasses WGD timing; a count of 1 still requires a timing estimate from eligible major-copy-number-two segments. If omitted, GRITIC infers the count. See [WGD timing estimation](#wgd-timing-estimation).
- `--random-seed` Seed stochastic inference with an unsigned 64-bit integer (`0` through `2**64 - 1`).
- `--unordered-balanced-route-prior` Use a uniform prior over unordered allele-route pairs. See [balanced route priors](#balanced-route-priors) for the weighting of ordered routes. This is disabled by default.

### Timing intervals

Posterior intervals default to contiguous empirical highest posterior density (HPD) intervals. Widths are probability mass greater than 0 and at most 1; method options accept `hpd` or `equal-tailed`. MUTIC and SIGTIC use the same `--posterior-summary-interval-width` and `--posterior-summary-interval-method` options, also defaulting to 95% HPD summary intervals.

Sample intervals are computed with the NumPy array interface of [ArviZ Stats](https://python.arviz.org/projects/stats/en/stable/array_stats_only.html): `hdi(method="nearest")` for HPD and `eti` for equal-tailed intervals. ArviZ's nearest HDI uses sorted endpoints `floor(width * number_of_draws)` indices apart; small-sample bounds can therefore differ from the previous implementation. A width of `1` uses the full observed range for either method.

| Interval family | Options | Default width | Controls |
| --- | --- | --- | --- |
| Route gain | `--route-gain-interval-width`, `--route-gain-interval-method` | 0.95 | Route-conditional gain bounds in both gain timing tables. |
| WGD overlap | `--wgd-overlap-interval-width`, `--wgd-overlap-interval-method` | 0.9 | Internal candidate-segment bounds that determine WGD overlap and can change WGD inference. |
| Sample WGD | `--wgd-timing-interval-width`, `--wgd-timing-interval-method` | 0.9 | Final sample-level WGD bounds shared by the calling-info JSON, route table, and yellow tree nodes. |
| Posterior summary | `--posterior-summary-interval-width`, `--posterior-summary-interval-method` | 0.95 | Gain and gain-conditioned WGD bounds in posterior summaries. |
| Tree gain | `--tree-gain-interval-width`, `--tree-gain-interval-method` | 0.9 | Blue gain-node labels in tree PDFs. |

### Tree plots

- `--plot-trees` Enable route-tree plots for each segment.

## Outputs

Outputs are written directly under `--sample-dir`. The file and directory names below are suffixes prefixed by `SAMPLE_ID`. We recommend only considering gained segments with 10 or more SNVs.

| File or directory suffix | Contents |
| --- | --- |
| [`_posterior_timing_table_summary_penalty_*.tsv`](#posterior-timing-summaries) | Main gain-timing summaries, with ordinary and penalized route weights. |
| [`_route_table.tsv`](#_route_tabletsv) | Route probabilities and segment metadata. |
| [`_gain_timing_table.tsv`](#_gain_timing_tabletsv) | Gain-node timing conditional on each route. |
| [`_wgd_calling_info.json`](#_wgd_calling_infojson) | Sample WGD call and timing. |
| [`_gain_timing_table_wgd_segments.tsv`](#_gain_timing_table_wgd_segmentstsv) | Preliminary segment timings used to evaluate WGD. |
| [`_tree_plots`](#_tree_plots) | Optional route-tree PDFs. |
| [`_mutation_table.tsv`](#_mutation_tabletsv) | Processed mutations and downstream identifiers. |
| [`_count_group_table.tsv`](#_count_group_tabletsv) | Distinct read-count pairs. |
| [`_phase_group_table.tsv`](#_phase_group_tabletsv) | Read-count groups subdivided by phasing. |
| [`_likelihood_context_table.tsv`](#_likelihood_context_tabletsv) | Distinct copy-number observation contexts. |
| [`_segment_context_table.tsv`](#_segment_context_tabletsv) | Segment-to-context mapping. |
| [`_count_group_likelihood_table.tsv`](#_count_group_likelihood_tabletsv) | Multiplicity likelihoods by count group and context. |
| [`_segment_group_table.tsv`](#_segment_group_tabletsv) | Mutation counts by segment and phase group. |
| [`_subclone_table.tsv`](#_subclone_tabletsv) | Processed subclone inputs. |
| [`_timing_dicts`](#timing-archives) | Posterior archives for downstream mutation timing. |

### Posterior timing summaries

`_posterior_timing_table_summary_penalty_<True|False>.tsv` summarizes independent-gain timing across routes, keyed by (`Sample_ID`, `Segment_ID`, `Gain_Index`). `Gain_Index` is a 1-based chronological index of independent gains.

| Result columns | Interpretation |
| --- | --- |
| `Timing_Median`, `Timing_Low_CI`, `Timing_High_CI` | Gain-timing median and configured posterior-summary interval. |
| `Proportion` | Fraction of posterior route draws in which this independent gain exists. |
| `WGD_Timing_Median`, `WGD_Timing_Low_CI`, `WGD_Timing_High_CI` | WGD-timing median and interval, conditional on this gain existing. |
| `Pre_WGD_Probability`, `Post_WGD_Probability` | Probabilities of the gain preceding or following WGD, conditional on this gain existing. |

Only gains with `Proportion >= 0.8` are reported. Its denominator includes all route draws for the segment, including routes with no independent gains.

Gain and WGD statistics on a row use the same subset of draws in which that gain exists. WGD-related fields are blank in non-WGD runs. Segment coordinates, copy numbers, mutation count, mutation rate, and WGD status accompany these results.

Two summary tables are produced for every run that produces timing output, including non-WGD runs. The filename ending in `_penalty_False.tsv` summarizes draws using `Probability` from the route table. The filename ending in `_penalty_True.tsv` summarizes a second set of draws using `Penalized_Probability`. See [posterior sampling](#posterior-sampling) for how these draws are constructed.

### _route_table.tsv

This table contains one row for each possible route of each timed segment, keyed by `(Sample_ID, Segment_ID, Route)`. Routes with no independently timeable gains are retained here.

`Timing_Representation` is `Route_Particles` for ordinary gained routes and `Uniform_No_Gain` when every extant copy spans the complete `[0,1]` interval. The table stores ordinary `Probability` and post-hoc `Penalized_Probability`, average event and loss counts, a timing-space sampling-density diagnostic (`Density`), runtime (`Time`), and segment and WGD metadata.

`Route` is an opaque, order-sensitive identifier of the complete allele route. See [posterior sampling](#posterior-sampling) for the density diagnostic.

`Penalized_Probability` is calculated by multiplying each route's ordinary probability by `exp(-2.7 * Average_N_Events)` and renormalizing across routes within the segment. Tree output uses ordinary probabilities, as does downstream mutation timing by default. See [Baker et al. (2024)][publication] for details of the penalty.

### _gain_timing_table.tsv

This table contains one row per independently timeable gain node per route and is keyed by four identifier columns: (`Sample`, `Segment`, `Route`, `Node`). A route without an independently timeable gain has no row in this table.

The remaining columns are `Node_Phasing`, `Timing` (the median), and `Timing_CI_Low` and `Timing_CI_High` (the configured route-gain interval). These timing statistics summarize the node's 1,000 likelihood-weighted samples conditional on the route.

`Node_Phasing` labels each node by its route component, `Major` or `Minor`. In a route with one extant allele component, every node is labelled `Major`; routes with two components use their assigned allele roles. Mutation-table `Phasing` uses lowercase `major` and `minor`, with missing values for unphased mutations.

### _wgd_calling_info.json

A JSON object that gives WGD calling information for the sample. Missing or nonfinite numeric values are represented by JSON `null`. Its keys are `WGD_Timing`, `WGD_Timing_CI_Low`, `WGD_Timing_CI_High`, `Major_CN_Mode`, `Overlap_Proportion`, `WGD_Status`, and `Best_Overlap_Timing`.

See [WGD timing estimation](#wgd-timing-estimation) for eligibility, the overlap decision, and the final timing estimate.

### _gain_timing_table_wgd_segments.tsv

This table is produced while evaluating WGD timing. It gives the preliminary non-WGD node posterior for every eligible major-copy-number-two segment. Displayed `Timing_CI_Low` and `Timing_CI_High` use the route-gain interval. `Intersecting`, `Best_Overlap_Timing`, and `Overlap_Proportion` use the separate, unrounded internal WGD-overlap intervals. Displayed bounds therefore need not determine `Intersecting`; see [timing intervals](#timing-intervals) for the settings.

### _tree_plots

Binary tree plots for the gain timings of each route in a segment. Each plot has one or two allele trees, according to the route. Blue nodes show independent gains with the tree-gain interval from the 1,000 route-conditional samples. Yellow nodes show WGD timing; red nodes are the extant copies at sampling.

### Mutation and subclone data

#### _mutation_table.tsv

The processed mutation table, including identifiers used by downstream mutation timing. `Segment_ID` is GRITIC's final coordinate-derived segment ID after the selected merge behavior.

Every output row contains these mutation identity, mapping, and provenance columns:

- `Source_Segment_ID` is the input segment ID when both input tables supplied matching `Segment_ID` columns. For position-based copy-number assignment it is the final assigned segment ID.
- `Mutation_ID` contains the literal input `Mutation_ID` when supplied and is blank otherwise.
- `Position` contains the canonical input position when supplied and is blank otherwise.
- `GRITIC_Mutation_ID` is the canonical sample-unique identifier derived from the source segment plus `Mutation_ID` when supplied, otherwise from the source segment plus `Position`. Its two components are URL-escaped and separated by `:`. Consumers should treat this value as opaque.
- `Segment_Mutation_Index` is a zero-based, consecutive index within the final `Segment_ID`. GRITIC assigns it by sorting `GRITIC_Mutation_ID` lexicographically within each segment.
- `Phase_Group_ID` is a sample-wide, zero-based identifier that links the mutation to one row of `_phase_group_table.tsv`. Together with `Segment_ID`, it also links to the mutation count in `_segment_group_table.tsv`.

Use `GRITIC_Mutation_ID` for downstream mutation joins.

Follow `Phase_Group_ID` through the phase- and count-group dictionaries to recover phasing and read counts. Follow `Segment_ID` through `_segment_context_table.tsv` and `_likelihood_context_table.tsv` to recover `Major_CN`, `Minor_CN`, and total copy number; `Gain_Type` is derived from the major/minor combination. The alternate-read correction is an internal segment/state quantity and is not emitted.

#### _count_group_table.tsv

This sample-wide dictionary stores one row for each distinct `(Tumor_Ref_Count, Tumor_Alt_Count)` pair in the retained sample. Its columns are `Sample_ID`, `Count_Group_ID`, `Tumor_Ref_Count`, and `Tumor_Alt_Count`. `Count_Group_ID` is zero-based and consecutive across the sample.

#### _phase_group_table.tsv

This sample-wide dictionary subdivides count groups by mutation phasing. Its columns are `Sample_ID`, `Phase_Group_ID`, `Count_Group_ID`, and `Phasing`, whose values are `non_phased`, `major`, or `minor`. `Phase_Group_ID` is zero-based and consecutive across the sample; `Count_Group_ID` references `_count_group_table.tsv`.

#### _likelihood_context_table.tsv

This dictionary contains one row for each distinct read-count observation context, with columns `Sample_ID`, `Likelihood_Context_ID`, `Major_CN`, `Minor_CN`, and `Normal_Total_CN`. Context IDs are zero-based and consecutive across the sample. `Normal_Total_CN` distinguishes observation models such as a haploid normal sex chromosome from a diploid normal autosome.

#### _segment_context_table.tsv

This table maps each final segment with retained mutations to one observation context. Its columns are `Sample_ID`, `Segment_ID`, and `Likelihood_Context_ID`.

#### _count_group_likelihood_table.tsv

This sparse table stores one likelihood vector for each `(Likelihood_Context_ID, Count_Group_ID)` pair used by at least one segment. Its fixed leading columns are `Sample_ID`, `Likelihood_Context_ID`, and `Count_Group_ID`, followed by `Prob_Mult_1`, `Prob_Mult_2`, and so on through the largest sample-wide major copy number, then `Prob_Subclone_0`, `Prob_Subclone_1`, and so on. Multiplicity columns above the row context's `Major_CN` are blank. Applicable entries are finite, nonnegative, and sum to one per row.

#### _segment_group_table.tsv

This sparse association table stores `Sample_ID`, `Segment_ID`, `Phase_Group_ID`, and `N_Mutations`. Its compound key is `(Segment_ID, Phase_Group_ID)`, and `N_Mutations` counts retained mutations in the segment belonging to that phase group.

#### _subclone_table.tsv

The retained and processed subclone inputs with the fixed columns `Cluster`, `Subclone_CCF`, `Subclone_Fraction`, and `N_SNVs`. This file is always written. If no subclone table was supplied, or no candidate survived filtering, it contains the header and zero rows.

## Model and archive details

### Balanced route priors

For an allele-balanced `N+N` state, GRITIC defaults to a uniform prior over the ordered Cartesian product of single-allele histories, conditional on the copy-number and WGD setting.

With `--unordered-balanced-route-prior`, an `A/A` route retains weight 1 and each orientation of an `A/B` pair receives weight 0.5, so the pair carries the same total prior weight as `A/A`. This changes only route-level prior weights and has no effect for unbalanced states or balanced states with only an identical-component route.

Without phased SNVs, reciprocal routes have equal likelihoods and probabilities. With phased SNVs, GRITIC evaluates their likelihoods and posteriors separately before applying the selected prior.

### Balanced WGD segments

When evaluating a possible WGD or pooling timing for a called WGD, GRITIC constrains the two homolog histories and timings of an original balanced `2+2` segment to be identical. It retains mutation `Phasing` labels but evaluates every SNV against one representative duplicated allele. The likelihood context retains the original `2+2` copy number, while timing uses a single pseudo-`2+0` route. If WGD is rejected, the ordinary independent-route `2+2` model uses mutation phasing.

Pooled [timing archives](#timing-archives) record target `2+2` WGD metadata and model `2+0` non-WGD metadata; their phasing masks and state slices address the representative major-allele geometry.

### WGD timing estimation

GRITIC calculates the modal `Major_CN`, weighted by segment width, across the configured autosomes. With automatic WGD inference, mode 1 gives a WGD count of 0. Mode 2 triggers timing of autosomal major-copy-number-two segments with at least 10 retained SNVs; at least one must produce a finite timing interval.

GRITIC finds the mutation-time point covered by the greatest total segment width across the candidate segments' internal WGD-overlap intervals. `Overlap_Proportion` is that covered width divided by the total width of eligible segments with finite intervals. It measures shared genomic span. An overlap of at least 60% gives an inferred WGD count of 1; a lower overlap gives 0 and emits a warning. Python callers can set this threshold with `process_sample(..., min_wgd_overlap=...)`.

The best point and overlap proportion are recorded in `Best_Overlap_Timing` and `Overlap_Proportion`. Overlapping segments are pooled and refit by minor-copy-number class, and their timing densities are combined into 500 sample-level draws. For a WGD call, `WGD_Timing` and its interval fields summarize those draws using the configured sample-WGD interval.

The [WGD-count option](#inference-model-and-wgd-calling) can override the inferred count. GRITIC warns when the supplied count conflicts with modal major copy number.

### Detection correction and mutation-share priors

GRITIC uses each segment's read counts to fit clonal and subclonal mutation shares, with the supplied subclone CCFs determining the expected VAFs of the subclonal states. The retained input fractions inform a weak Dirichlet prior; the clonal input fraction is one minus their sum.

With `--subclone-fraction-prior adjusted`, GRITIC divides the clonal and subclonal input fractions by their estimated detection probabilities separately for each segment, then renormalizes them. With `supplied`, it uses the fractions directly. Both modes set the Dirichlet parameters to one plus the resulting fractions and apply the segment-specific detection correction in the likelihood.

The correction uses Poisson thinning to calculate the probability of meeting `--min-mutation-alt-count`, given each state's expected VAF and the segment's estimated mean coverage. Coverage is averaged over SNVs with observed VAF greater than the `--coverage-vaf-quantile` quantile minus 0.01. This quantile affects the likelihood correction in both prior modes and the prior adjustment only in `adjusted` mode.

The `supplied` prior follows [Baker et al. (2024)][publication], Supplementary Methods §8.3.4, equation 22. The publication's detection correction averages detection power over individual SNV depths (§8.3.5, equation 23); the current implementation uses the mean-coverage approach above.

### Posterior sampling

GRITIC retains 1,000 likelihood-resampled particles for each sampled route. The `Density` diagnostic reports the fraction of tested timing points with another sampled point within the sampling neighborhood. It covers hit-and-run timing coordinates; independently sampled Dirichlet clone shares are not chain dimensions. Analytic uniform routes require no chain and report density 1.

Posterior-summary draws select a route, then use one posterior-array index for its WGD and all independent-gain timings, preserving their joint dependence. The [gain timing table](#_gain_timing_tabletsv) supplies node identity and phasing metadata; the aligned arrays and authoritative route weights are in the [timing archives](#timing-archives).

The gain-draw and route-ledger data frames used to calculate these summaries are internal.

### Timing archives

The `_timing_dicts` directory contains compressed posterior archives used by MUTIC for mutation timing. Gain-summary tables can be inspected directly.

`Route_Particles` segments store gain timing and multiplicity particles. A `Uniform_No_Gain` segment writes no store when it has no subclones; with subclones it writes only the fitted clone-share posterior. Each logical store that exists is a required pair:

- `SEGMENT_ID_timing_dict.npz` contains the numeric tables in linear order as `table_000000`, `table_000001`, and so on. It is a compressed NumPy archive and never contains pickled or object-dtype arrays.
- `SEGMENT_ID_timing_dict.manifest.json` maps the original nested hierarchy onto those table indexes. It also records the format version, archive filename and SHA-256, and each table's dtype and shape.

Pooled WGD stores use the same pair with a `WGD_minor_cn_N` identifier in place of `SEGMENT_ID`. From the sample output directory, load a store through GRITIC to validate the pair and reconstruct the dictionary hierarchy:

```python
from gritic.timingio import load_timing_archive

timing_dict = load_timing_archive(
    'SAMPLE_ID_timing_dicts/1-0-200_timing_dict.npz'
)
```

The reconstructed dictionary keys are opaque route identifiers. Each `Route_Particles` entry contains the following numeric arrays:

- `Probability` and `Penalized_Probability` are one-element float64 arrays containing the authoritative route weights used for mutation timing.
- `Timing` is a float64 `N_Particles x N_Timing_Nodes` matrix. `Timing_Node_ID` is the aligned int64 identifier for each column, corresponding to `Node` in the [gain timing table](#_gain_timing_tabletsv).
- `WGD_Timing` is the aligned float64 WGD vector and `Mult` is the aligned float64 multiplicity matrix. Non-WGD models use an all-NaN `WGD_Timing` vector; WGD models use finite values.
- `Interval_Start_Source`, `Interval_End_Source`, `Interval_Multiplicity`, and `Interval_Phasing` describe every mutation-bearing route interval without requiring the consumer to reconstruct a route tree. Source 0 is constant zero, source 1 is constant one, source 2 is `WGD_Timing`, and source `3+j` is column `j` of `Timing`. Phasing is a bit mask: 1 is `non_phased`, 2 is `major`, and 4 is `minor`.
- `State_Column_Offsets` and `State_Columns` map phasings onto `Mult` columns. The four offsets delimit three slices in `non_phased`, `major`, `minor` order. `Mult` contains the corresponding unphased, major-phased, and minor-phased multiplicity blocks followed by shared subclone columns, which occur at the tail of every applicable slice.
- `Target_Major_CN`, `Target_Minor_CN`, and `Target_WGD_Status` describe the genomic segment to which inference applies. `Model_Major_CN`, `Model_Minor_CN`, and `Model_WGD_Status` describe the fitted particle geometry. `N_Subclones` gives the shared number of subclone columns.
- `Archive_Kind` is 0 for an ordinary segment store and 1 for a pooled-WGD store. The [balanced WGD model](#balanced-wgd-segments) explains the target/model metadata and representative-allele encoding for pooled `2+2` segments.

A subclonal `Uniform_No_Gain` route contains only `Clone_Share`, whose columns contain fitted mutation shares for the clonal cluster followed by each subclone.

The 1,000 rows of `Timing`, `WGD_Timing`, and `Mult` are aligned; use the same row index for one joint posterior draw.

[publication]: https://aacrjournals.org/cancerdiscovery/article/14/10/1810/748591/The-History-of-Chromosomal-Instability-in-Genome
