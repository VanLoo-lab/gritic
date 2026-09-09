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

Run this example from the repository root. It writes directly to `examples/output/TEST_ID`, which must not already exist, even if empty. To reuse it, supply `--overwrite`.

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

Use `gritic --help` for command-line help or see [run options](#run-options). For programmatic use, see [Python API](#python-api).

## Input tables

Supply tab-separated copy-number and mutation tables, plus an optional subclone table. Inputs from any caller are accepted if they follow these schemas; see [`examples/`](examples/).

Column names are case-sensitive. Additional caller-supplied columns are accepted.

### Mutation table

All SNVs for the sample, with either `Mutation_ID` or `Position` required alongside the three mandatory columns below.

| Column | Requirement | Values |
| --- | --- | --- |
| `Chromosome` | Required | Labels following the [chromosome-handling rules](#chromosome-handling). |
| `Tumor_Ref_Count` | Required | Number of reads supporting the reference allele. |
| `Tumor_Alt_Count` | Required | Number of reads supporting the alternate allele. |
| `Mutation_ID` | Required unless `Position` is supplied | Literal text, preserving values such as `000123`, `NA`, and `NULL`. |
| `Position` | Required unless `Mutation_ID` is supplied | Zero-based genomic position. |
| `Phasing` | Optional | `major` or `minor`, case-insensitive, ignoring surrounding whitespace. Missing values remain unphased. |

Supplying `--drop-unrecognized-phasing` drops mutations with unrecognized phasing and reports the number dropped and the unrecognized values.

Input mutation columns named `Segment_Start`, `Segment_End`, `Major_CN`, or `Minor_CN` are ignored because GRITIC always annotates those values from the copy-number table.

GRITIC supports two segment-assignment modes:

- If both the mutation and copy number tables contain `Segment_ID`, GRITIC uses those IDs to associate mutations with the input copy number segments. `Position` is not needed for assignment.
- Otherwise, the mutation table must contain `Position`, which assigns mutations to segments on the same chromosome using the [coordinate convention](#copy-number-table).

In supplied-ID mode, both tables require nonblank `Segment_ID` values, and those IDs must be unique in the copy-number table. Every mutation's `Segment_ID` must match a copy-number row, including its `Chromosome`. If `--drop-unmatched-snvs` is supplied, GRITIC drops mutations unmatched by either assignment mode and reports the number dropped. Missing or blank supplied IDs and chromosome mismatches remain errors.

GRITIC does not model mutations in zero-copy `0+0` segments. In supplied-ID
mode, GRITIC drops mutations assigned to these segments and reports the mutation
count, segment count, and source `Segment_ID` values. This is independent of
`--drop-unmatched-snvs` because those IDs are matched rather than unknown.

The selected `Mutation_ID` or canonical integer `Position` value must be unique within its source segment. An explicit `Mutation_ID` column takes precedence for every row; GRITIC does not fall back to `Position` for blank values in that column. The same selected value may be reused in different source segments. With supplied segment IDs, this also applies to source segments that GRITIC subsequently merges.

### Copy number table

The rounded allele-specific copy-number profile for the sample requires `Chromosome`, `Segment_Start`, `Segment_End`, `Major_CN`, and `Minor_CN`. Chromosome labels follow the [rules below](#chromosome-handling). Segments use zero-based, half-open intervals. Intervals on the same chromosome must not overlap. Supplied `Segment_ID` values enable the assignment mode described under [Mutation table](#mutation-table).

#### Timing eligibility

The supported allele-specific copy-number states are:

| `Major_CN` | Permitted `Minor_CN` |
| --- | --- |
| 1–5 | 0 through `Major_CN` |
| 6 | 0–4 |
| 7 | 0–3 |
| 8 | 0–1 |

The `1+0` and `1+1` states use uniform no-gain timing. In WGD runs, major-copy-number-two segments are excluded from the ordinary per-segment fit. [WGD timing estimation](#wgd-timing-estimation) uses those on configured autosomes with at least 10 retained SNVs. Balanced `2+2` segments use a [WGD-specific allele model](#balanced-wgd-segments).

#### Segment merging

By default, GRITIC merges consecutive segments on the same chromosome having identical `Major_CN` and `Minor_CN`. The merged interval runs from the first segment's start through the last segment's end, including any intervening uncovered bases. Use `--max-merge-gap` to limit gaps or `--keep-adjacent-segments` to preserve input segments. Final segment IDs are generated from `Chromosome`, `Segment_Start`, and `Segment_End` after merging.

#### Chromosome handling

Both input tables accept numbered chromosomes from 1 through `--autosome-count` and the sex chromosomes present under the selected karyotype. Chromosome labels are case-sensitive; `chr` prefixes are removed.

| Karyotype | Present sex chromosomes | Normal X or Z copies | Normal Y or W copies |
| --- | --- | --- | --- |
| XX | X | 2 | 0 |
| XY | X, Y | 1 | 1 |
| ZZ | Z | 2 | 0 |
| ZW | Z, W | 1 | 1 |

When sample sex is not supplied, GRITIC infers it from the copy-number table: Y implies `XY` and W implies `ZW`; otherwise X alone implies XX and Z alone implies ZZ. If no sex chromosome is represented, GRITIC defaults to `XX`, so callers using another system should supply `--sample-sex`. Chromosomes outside the configured set follow `--drop-unmatched-chromosomes` behavior.

### Subclone table

The optional subclone table gives the identified subclonal peaks and their assigned mutation fractions for the sample.

Required columns are:

- `Cluster`: the subclone identifier.
- `Subclone_CCF`: cancer cell fraction. CCF determines the expected VAF of a subclone state.
- `Subclone_Fraction`: the fraction of input SNVs assigned to the subclone. Mutation shares determine its mixture prior.

Candidates are filtered using the [subclone-handling options](#subclone-handling). If no subclones remain, GRITIC uses its clonal-only model. If there are more than two subclones, GRITIC groups them into two: the subclone with the largest CCF is unmodified and the remaining clones are combined by summing their fractions and taking their fraction-weighted mean CCF.

Retained `Subclone_Fraction` values keep their input scale; the remaining share is treated as clonal.

## Run options

### Required run arguments

- `--mutation-table` A path to the [mutation table](#mutation-table) for the sample.
- `--copy-number-table` A path to the [copy-number table](#copy-number-table) for the sample.
- `--purity` The estimated cellular purity for the sample.
- `--sample-id` Sample ID used as an output filename prefix. It must be a cross-platform-safe filename component.
- `--sample-dir` Sample output directory.

### Output handling

- `--overwrite` Reuse the existing `--sample-dir` directory and its subdirectories. Files with the same names as new outputs are overwritten; all other files are preserved.
- `--plot-trees` Enable route-tree plots for each segment.

### Genome and input handling

- `--autosome-count` Number of numbered autosomes in the organism (default: 22). This defines the accepted numbered chromosome labels and the chromosomes eligible for WGD inference.
- `--sample-sex` Override the inferred karyotype with `XX`, `XY`, `ZZ`, or `ZW`. See [chromosome handling](#chromosome-handling) for inference rules and normal copy numbers.
- `--drop-unmatched-chromosomes` Drop copy-number and mutation rows whose chromosome is not one of the configured autosomes or present sex chromosomes and report the number of rows dropped. By default, any such chromosome is an error.
- `--drop-unmatched-snvs` Drop mutation rows that cannot be associated with a copy-number segment by either supplied `Segment_ID` or genomic `Position` and report the number dropped. By default, unmatched mutations raise an error.
- `--drop-unrecognized-phasing` Drop mutation rows whose non-missing `Phasing` value is not `major` or `minor` and report the number dropped. By default, unrecognized phasing labels raise an error.
- `--keep-adjacent-segments` Preserve input copy-number segments separately. See [segment merging](#segment-merging) for the default behavior.
- `--max-merge-gap N` Merge consecutive equal-copy-number segments only when the gap is at most `N` bases. If omitted, there is no maximum; use `0` to merge only intervals that touch.

### Mutation filtering and detection correction

- `--min-mutation-alt-count` Minimum `Tumor_Alt_Count` needed to retain a mutation (default: 3).
- `--min-mutation-coverage` Minimum `Tumor_Ref_Count + Tumor_Alt_Count` needed to retain a mutation (default: 10).
- `--coverage-vaf-quantile` Observed-SNV VAF quantile used to select mutations for the mean-coverage estimate in the detection correction (default: 0.9). See [detection correction and mutation-share priors](#detection-correction-and-mutation-share-priors).

### Subclone handling

- `--subclone-table` A path to the [subclone table](#subclone-table) for the sample. If omitted, GRITIC assumes every SNV is clonal, which can bias gain timings earlier.
- `--clip-subclone-ccf` Clip out-of-range `Subclone_CCF` values before filtering. With the default CCF filters, values clipped to either boundary are subsequently excluded.
- `--min-subclone-ccf` Minimum `Subclone_CCF` retained as a subclone, inclusive (default: 0.01).
- `--max-subclone-ccf` Maximum `Subclone_CCF` retained as a subclone, inclusive (default: 0.9).
- `--min-subclone-fraction` A subclone's normalized share of the subclonal mutation fractions after CCF filtering must be strictly greater than this threshold (default: 0.1).
- `--subclone-fraction-prior {adjusted,supplied}` Use detection-adjusted mutation fractions or the supplied fractions in the mutation-share prior (default: `adjusted`). See [detection correction and mutation-share priors](#detection-correction-and-mutation-share-priors).

The value of `--min-subclone-ccf` must be greater than 0.

### Inference model and WGD calling

Every run requires configured autosomal segments with a segment-width-weighted modal `Major_CN` of 1 or 2, including runs with a supplied WGD count.

- `--wgd-count {0,1}` Override GRITIC's inferred WGD count. If omitted, GRITIC infers the count. GRITIC warns when the supplied count conflicts with modal major copy number. See [WGD timing estimation](#wgd-timing-estimation).
- `--random-seed` Seed stochastic inference with an unsigned 64-bit integer. If omitted, GRITIC generates a random seed using sources provided by the operating system.
- `--unordered-balanced-route-prior` Use a uniform prior over unordered allele-route pairs. See [balanced route priors](#balanced-route-priors) for the weighting of ordered routes.

### Timing intervals

Intervals default to highest posterior density (HPD, `hpd`); `equal-tailed` is also available. Interval widths specify probability mass.

| Interval family | Options | Default width | Controls |
| --- | --- | --- | --- |
| Route gain | `--route-gain-interval-width`, `--route-gain-interval-method` | 0.95 | Route-conditional gain bounds in both gain timing tables. |
| WGD overlap | `--wgd-overlap-interval-width`, `--wgd-overlap-interval-method` | 0.9 | Internal candidate-segment bounds that determine WGD overlap and can change WGD inference. |
| Sample WGD | `--wgd-timing-interval-width`, `--wgd-timing-interval-method` | 0.9 | Final sample-level WGD bounds shared by the calling-info JSON, route table, and yellow tree nodes. |
| Posterior summary | `--posterior-summary-interval-width`, `--posterior-summary-interval-method` | 0.95 | Gain and gain-conditioned WGD bounds in posterior summaries. |
| Tree gain | `--tree-gain-interval-width`, `--tree-gain-interval-method` | 0.9 | Blue gain-node labels in tree PDFs. |

## Outputs

Append each suffix below to `SAMPLE_ID`. We recommend only considering gained segments with 10 or more SNVs.

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
| `WGD_Timing_Median`, `WGD_Timing_Low_CI`, `WGD_Timing_High_CI` | WGD-timing median and interval. |
| `Pre_WGD_Probability`, `Post_WGD_Probability` | Probabilities of the gain preceding or following WGD. |

Only gains with `Proportion >= 0.8` are reported. Its denominator includes all route draws for the segment, including routes with no independent gains.

Gain and WGD statistics on a row use the same subset of draws in which that gain exists. WGD-related fields are blank in non-WGD runs. Segment coordinates, copy numbers, mutation count, mutation rate, and WGD status accompany these results.

`_penalty_False.tsv` uses `Probability`; `_penalty_True.tsv` uses a separate set of draws with `Penalized_Probability`. See [posterior sampling](#posterior-sampling) for how these draws are constructed.

### _route_table.tsv

This table contains one row for each possible route of each timed segment, keyed by `(Sample_ID, Segment_ID, Route)`. Routes with no independently timeable gains are retained here.

`Timing_Representation` is `Route_Particles` for ordinary gained routes and `Uniform_No_Gain` when every extant copy spans the `[0,1]` interval. The table stores ordinary `Probability` and post-hoc `Penalized_Probability`, average event and loss counts, a timing-space sampling-density diagnostic (`Density`), runtime (`Time`), and segment and WGD metadata.

`Route` is an order-sensitive identifier of the complete allele route. See [posterior sampling](#posterior-sampling) for the density diagnostic.

`Penalized_Probability` is calculated by multiplying each route's ordinary probability by `exp(-2.7 * Average_N_Events)` and renormalizing across routes within the segment. Tree output uses ordinary probabilities. See [Baker et al. (2024)][publication] for details of the penalty.

### _gain_timing_table.tsv

This table contains one row per independently timeable gain node per route and is keyed by four identifier columns: (`Sample`, `Segment`, `Route`, `Node`).

The remaining columns are `Node_Phasing`, `Timing` (the median), and `Timing_CI_Low` and `Timing_CI_High` (the configured route-gain interval). Timing statistics are conditional on the route.

`Node_Phasing` is `Major` or `Minor`; a route with one extant allele component labels every node `Major`.

### _wgd_calling_info.json

Keys are `WGD_Timing`, `WGD_Timing_CI_Low`, `WGD_Timing_CI_High`, `Major_CN_Mode`, `Overlap_Proportion`, `WGD_Status`, and `Best_Overlap_Timing`; missing or nonfinite numeric values are `null`.

See [WGD timing estimation](#wgd-timing-estimation) for eligibility, the overlap decision, and the final timing estimate.

### _gain_timing_table_wgd_segments.tsv

This table is produced while evaluating WGD timing. It gives the preliminary non-WGD node posterior for every eligible major-copy-number-two segment. Displayed `Timing_CI_Low` and `Timing_CI_High` use the route-gain interval. `Intersecting`, `Best_Overlap_Timing`, and `Overlap_Proportion` use the separate, unrounded internal WGD-overlap intervals. Displayed bounds therefore need not determine `Intersecting`; see [timing intervals](#timing-intervals) for the settings.

### _tree_plots

Binary tree plots for the gain timings of each route in a segment. Each plot has one or two allele trees, according to the route. Blue nodes show independent gains with the tree-gain interval. Yellow nodes show WGD timing; red nodes are the extant copies at sampling.

### Processed inputs and downstream tables

#### _mutation_table.tsv

The processed mutation table, including identifiers used by downstream mutation timing. `Segment_ID` is GRITIC's final coordinate-derived segment ID after the selected merge behavior.

Every output row contains these mutation identity, mapping, and provenance columns:

- `Source_Segment_ID` is the input segment ID when both input tables supplied matching `Segment_ID` columns. For position-based copy-number assignment it is the final assigned segment ID.
- `Mutation_ID` contains the literal input `Mutation_ID` when supplied and is blank otherwise.
- `Position` contains the canonical input position when supplied and is blank otherwise.
- `GRITIC_Mutation_ID` is the canonical sample-unique identifier derived from the source segment plus `Mutation_ID` when supplied, otherwise from the source segment plus `Position`.
- `Segment_Mutation_Index` is a zero-based, consecutive index within the final `Segment_ID`.
- `Phase_Group_ID` links the mutation to one row of [`_phase_group_table.tsv`](#_phase_group_tabletsv). Together with `Segment_ID`, it also links to the mutation count in [`_segment_group_table.tsv`](#_segment_group_tabletsv).

Use `GRITIC_Mutation_ID` for downstream mutation joins.

Follow `Phase_Group_ID` through the [phase](#_phase_group_tabletsv) and [count](#_count_group_tabletsv) dictionaries to recover phasing and read counts. Follow `Segment_ID` through [`_segment_context_table.tsv`](#_segment_context_tabletsv) and [`_likelihood_context_table.tsv`](#_likelihood_context_tabletsv) to recover `Major_CN`, `Minor_CN`, and total copy number.

#### _subclone_table.tsv

The retained and processed subclone inputs with the fixed columns `Cluster`, `Subclone_CCF`, `Subclone_Fraction`, and `N_SNVs`. This file is always written. If no subclone table was supplied, or no candidate survived filtering, it contains the header and zero rows.

GRITIC derives `N_SNVs` from the retained mutation count and each retained/combined `Subclone_Fraction`, overwriting an input `N_SNVs` column.

## Model details

### Balanced route priors

For an allele-balanced `N+N` state, GRITIC defaults to a uniform prior over the ordered Cartesian product of single-allele histories, conditional on the copy-number and WGD setting.

With `--unordered-balanced-route-prior`, an `A/A` route retains weight 1 and each orientation of an `A/B` pair receives weight 0.5, so the pair carries the same total prior weight as `A/A`.

Without phased SNVs, reciprocal routes have equal likelihoods and probabilities. With phased SNVs, GRITIC evaluates their likelihoods and posteriors separately before applying the selected prior.

### Balanced WGD segments

When evaluating a possible WGD or pooling timing for a called WGD, GRITIC constrains the two homolog histories and timings of an original balanced `2+2` segment to be identical. It retains mutation `Phasing` labels but evaluates every SNV against one representative duplicated allele. The likelihood context retains the original `2+2` copy number, while timing uses a single pseudo-`2+0` route. If WGD is rejected, the ordinary independent-route `2+2` model uses mutation phasing.

See [timing archives](#timing-archives) for the pooled-WGD encoding.

### WGD timing estimation

GRITIC calculates the modal `Major_CN`, weighted by segment width, across the configured autosomes. With automatic WGD inference, mode 1 gives a WGD count of 0. Mode 2 triggers timing of autosomal major-copy-number-two segments with at least 10 retained SNVs; at least one must produce a finite timing interval.

`Best_Overlap_Timing` is the mutation-time point covered by the greatest total segment width across the candidate segments' internal WGD-overlap intervals. `Overlap_Proportion` is that covered width divided by the total width of eligible segments with finite intervals. An overlap of at least 60% gives an inferred WGD count of 1; a lower overlap gives 0 and emits a warning.

Overlapping segments are pooled and refit by minor-copy-number class, and their timing densities are combined into 500 sample-level draws. For a WGD call, `WGD_Timing` and its interval fields summarize those draws using the configured sample-WGD interval.

### Detection correction and mutation-share priors

GRITIC uses each segment's read counts to fit clonal and subclonal mutation shares, with the supplied subclone CCFs determining the expected VAFs of the subclonal states. The clonal input fraction is one minus the sum of the retained subclonal fractions. GRITIC constructs a weak Dirichlet prior with parameters equal to one plus the supplied or adjusted fractions.

With `--subclone-fraction-prior adjusted`, GRITIC divides the clonal and subclonal input fractions by their estimated detection probabilities separately for each segment, then renormalizes them. With `supplied`, it uses the fractions directly. Both modes apply the segment-specific detection correction in the likelihood.

The correction uses Poisson thinning to calculate the probability of meeting `--min-mutation-alt-count`, given each state's expected VAF and the segment's estimated mean coverage. Coverage is averaged over SNVs with observed VAF greater than the `--coverage-vaf-quantile` quantile minus 0.01. This quantile affects the likelihood correction in both prior modes and the prior adjustment only in `adjusted` mode.

The `supplied` prior follows [Baker et al. (2024)][publication], Supplementary Methods §8.3.4, equation 22. The publication's detection correction averages detection power over individual SNV depths (§8.3.5, equation 23); the current implementation uses the mean-coverage approach above.

### Posterior sampling

GRITIC retains 1,000 likelihood-resampled particles for each sampled route. The `Density` diagnostic reports the fraction of tested timing points with another sampled point within the sampling neighborhood. It covers hit-and-run timing coordinates. Analytic uniform routes require no chain and report density 1.

Each posterior-summary draw samples a route and its gain and WGD timings jointly, preserving dependence between events. The [gain timing table](#_gain_timing_tabletsv) supplies node identity and phasing metadata; the [timing archives](#timing-archives) contain the joint timing draws and route weights used for downstream inference.

## Downstream formats

`Count_Group_ID`, `Phase_Group_ID`, and `Likelihood_Context_ID` are sample-wide, zero-based, consecutive identifiers.

### _count_group_table.tsv

One row per distinct `(Tumor_Ref_Count, Tumor_Alt_Count)` pair in the retained sample, with columns `Sample_ID`, `Count_Group_ID`, `Tumor_Ref_Count`, and `Tumor_Alt_Count`.

### _phase_group_table.tsv

Count groups subdivided by mutation phasing, with columns `Sample_ID`, `Phase_Group_ID`, `Count_Group_ID`, and `Phasing`. Phasing values are `non_phased`, `major`, or `minor`; `Count_Group_ID` references [`_count_group_table.tsv`](#_count_group_tabletsv).

### _likelihood_context_table.tsv

One row per distinct read-count observation context, with columns `Sample_ID`, `Likelihood_Context_ID`, `Major_CN`, `Minor_CN`, and `Normal_Total_CN` (normal total copy number).

### _segment_context_table.tsv

Each final segment with retained mutations mapped to one observation context, with columns `Sample_ID`, `Segment_ID`, and `Likelihood_Context_ID`.

### _count_group_likelihood_table.tsv

One likelihood vector per `(Likelihood_Context_ID, Count_Group_ID)` pair used by at least one segment. The fixed leading columns are `Sample_ID`, `Likelihood_Context_ID`, and `Count_Group_ID`, followed by `Prob_Mult_1`, `Prob_Mult_2`, and so on through the largest sample-wide major copy number, then `Prob_Subclone_0`, `Prob_Subclone_1`, and so on. Multiplicity columns above the row context's `Major_CN` are blank. Applicable entries sum to one per row.

### _segment_group_table.tsv

Columns `Sample_ID`, `Segment_ID`, `Phase_Group_ID`, and `N_Mutations`, keyed by (`Segment_ID`, `Phase_Group_ID`). `N_Mutations` counts retained mutations in the segment belonging to that phase group.

### Timing archives

The `_timing_dicts` directory contains compressed posterior archives for mutation timing.

`Route_Particles` segments store gain timing and multiplicity particles. A `Uniform_No_Gain` segment writes no store when it has no subclones; with subclones it writes only the fitted clone-share posterior. Each archive consists of two files:

- `SEGMENT_ID_timing_dict.npz`
- `SEGMENT_ID_timing_dict.manifest.json`

Pooled WGD stores use the same pair with a `WGD_minor_cn_N` identifier in place of `SEGMENT_ID`. See [Reading timing archives](#reading-timing-archives) for a Python example.

The rows of `Timing`, `WGD_Timing`, and `Mult` are aligned; use the same row index for one joint posterior draw.

A subclonal `Uniform_No_Gain` route contains only `Clone_Share`, whose columns contain fitted mutation shares for the clonal cluster followed by each subclone.

<details>
<summary>Archive schema</summary>

Each `Route_Particles` entry contains the following numeric arrays:

- `Probability` and `Penalized_Probability` are one-element float64 arrays containing the route weights used for downstream inference.
- `Timing` is a float64 `N_Particles x N_Timing_Nodes` matrix. `Timing_Node_ID` is the aligned int64 identifier for each column, corresponding to `Node` in the [gain timing table](#_gain_timing_tabletsv).
- `WGD_Timing` is the aligned float64 WGD vector and `Mult` is the aligned float64 multiplicity matrix. Non-WGD models use an all-NaN `WGD_Timing` vector; WGD models use finite values.
- `Interval_Start_Source`, `Interval_End_Source`, `Interval_Multiplicity`, and `Interval_Phasing` describe every mutation-bearing route interval without requiring the consumer to reconstruct a route tree. Source 0 is constant zero, source 1 is constant one, source 2 is `WGD_Timing`, and source `3+j` is column `j` of `Timing`. Phasing is a bit mask: 1 is `non_phased`, 2 is `major`, and 4 is `minor`.
- `State_Column_Offsets` and `State_Columns` map phasings onto `Mult` columns. The four offsets delimit three slices in `non_phased`, `major`, `minor` order. `Mult` contains the corresponding unphased, major-phased, and minor-phased multiplicity blocks followed by shared subclone columns, which occur at the tail of every applicable slice.
- `Target_Major_CN`, `Target_Minor_CN`, and `Target_WGD_Status` describe the genomic segment to which inference applies. `Model_Major_CN`, `Model_Minor_CN`, and `Model_WGD_Status` describe the fitted particle geometry. `N_Subclones` gives the shared number of subclone columns.
- `Archive_Kind` is 0 for an ordinary segment store and 1 for a pooled-WGD store. For pooled `2+2` segments, target metadata describe `2+2` with WGD and model metadata describe `2+0` without WGD. Phasing masks and state slices address the representative major allele; see the [balanced WGD model](#balanced-wgd-segments).

</details>

## Python API

### Running GRITIC

Run this example from the repository root. It writes directly to `examples/output/TEST_ID`. Pass `overwrite=True` to `gritictimer.process_sample` to reuse an existing sample directory; files with matching names are overwritten and other files are preserved.

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

### Input handling

For unmatched-SNV dropping with supplied segment IDs, `dataloader.load_input_tables(..., drop_unmatched_snvs=True)` removes unmatched rows during loading. For position-based assignment, pass `drop_unmatched_snvs=True` to `sampletools.Sample`.

For [segment merging](#segment-merging), `sampletools.Sample` accepts `max_merge_gap` to limit gaps and `merge_cn=False` to preserve input segments.

`min_subclone_ccf` must be greater than 0.

### Timing intervals and WGD inference

Python callers configure intervals through the `interval_config` argument to `process_sample`, using `distributiontools.TimingIntervalConfig` and `distributiontools.IntervalSpec`.

Set the [WGD overlap threshold](#wgd-timing-estimation) with `process_sample(..., min_wgd_overlap=...)`; the default is 0.6.

### Reading timing archives

From the sample output directory, load a store through GRITIC to validate the pair and reconstruct the dictionary hierarchy:

```python
from gritic.timingio import load_timing_archive

timing_dict = load_timing_archive(
    'SAMPLE_ID_timing_dicts/1-0-200_timing_dict.npz'
)
```

The reconstructed dictionary keys are route identifiers. See [Timing archives](#timing-archives) for the array definitions and alignment rules.

[publication]: https://aacrjournals.org/cancerdiscovery/article/14/10/1810/748591/The-History-of-Chromosomal-Instability-in-Genome
