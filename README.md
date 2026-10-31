# 🧬 Protein List Builder for GO Analysis

A **Python + Streamlit application** that adds UniProt-derived Gene Ontology (GO) annotations to a protein half-life list and helps explore half-life distributions, statistical differences, and overlap between GO terms.

The application reads the `T-half` and `UniProtGO` sheets from a single Excel workbook, matches gene names, and generates an annotated protein list.

> **Scope:** This application uses previously calculated half-lives for annotation and distribution comparisons. It does not estimate half-lives, query UniProt online, propagate annotations through the GO hierarchy, or perform conventional GO enrichment analyses such as over-representation analysis.
>
> This README describes the supplied `app.py`. The application code has not been modified.

---

## 1. Features

| Feature | Description |
|---|---|
| Gene-name matching | Matches the first column of `T-half` against `UniProtGO.From` |
| Annotation-row selection | Prioritizes reviewed records, then selects the longest candidate |
| GO and Entry transfer | Copies CC, BP, and MF annotations and the selected Entry |
| Term summaries | Splits semicolon-delimited annotations and lists unique terms with occurrence counts |
| Boxplots | Displays half-life distributions on a logarithmic axis and highlights terms with the highest and lowest medians |
| Per-term tests | Compares annotated and non-annotated rows using Mann–Whitney U tests, Cliff's delta, and BH adjustment |
| Term overlap | Calculates Jaccard and overlap coefficients |
| Conditional comparisons | Excludes rows carrying a specified GO annotation from the analysis population before comparison |
| Export | Downloads an annotated workbook, category-specific PNG plots, and category-specific statistical workbooks |

The three categories are **cellular component (CC)**, **biological process (BP)**, and **molecular function (MF)**.

---

## 2. Installation and Startup

### Required packages

```bash
python -m pip install streamlit pandas openpyxl matplotlib
```

The supplied code does not specify a minimum Python version or pinned dependency versions. For reproducibility, record the versions used in your validated environment.

### Launch the application

Save the supplied file as `app.py`, open a terminal in the same directory, and run:

```bash
streamlit run app.py
```

If you rename the file to `protein_list_go_app.py`, use:

```bash
streamlit run protein_list_go_app.py
```

Use the browser address shown in the startup output. Plotting uses Matplotlib's non-interactive `Agg` backend.

### Minimal project structure

```text
project/
├── app.py
├── README.md
└── README-en.md
```

Upload the input workbook through the application interface. The code does not require external API keys or database connection settings.

---

## 3. Input Workbook Requirements

### 3.1 File types and sheet names

- Accepted upload formats: **`.xlsx` and `.xlsm`**.
- Default sheet names: **`T-half` and `UniProtGO`**.
- Different sheet names can be specified in the sidebar.
- The first row of each sheet is read as the header.
- Output is a new `.xlsx` workbook. The application does not preserve input macros, formatting, or additional sheets.

### 3.2 T-half sheet

The application reads the **first two columns by position**. Their input headers may differ, but the first column must contain gene names and the second must contain half-lives.

| Position | Suggested header | Content |
|---|---|---|
| A | `Gene Name` | String matched against `UniProtGO.From` |
| B | `HalfLife` | Half-life; plots and statistical outputs interpret the values as days |

Example:

| Gene Name | HalfLife |
|---|---:|
| GENE-A | 12.5 |
| GENE-B | 48.0 |
| GENE-C | 7.2 |

- Columns after the first two are not carried into the completed T-half table.
- Input row order and duplicate rows are retained. Duplicate genes are not merged or averaged.
- No unit conversion is performed. Convert values to days before uploading if necessary.

### 3.3 UniProtGO sheet

The following **seven columns are referenced by header name**. Their order is not fixed, but their names must match exactly.

| Required header | Purpose |
|---|---|
| `From` | Gene-name matching key |
| `Entry` | UniProt accession of the selected record |
| `Reviewed` | Record-selection priority |
| `Length` | Selection of the longest eligible candidate |
| `Gene Ontology (cellular component)` | CC annotations |
| `Gene Ontology (biological process)` | BP annotations |
| `Gene Ontology (molecular function)` | MF annotations |

Comments in the source describe expected column letters A, B, C, H, I, J, and L. However, the actual UniProtGO processing uses **header names**, not those fixed positions.

Example GO cell:

```text
example component [GO:1234567]; another component [GO:7654321]
```

The names and identifiers above illustrate the format only. Use appropriate annotations from your actual data source.

- Use a standard semicolon `;` as the separator.
- Leading and trailing whitespace is removed from each split term; empty terms are ignored.
- Term summaries and ordinary per-term tests use the **complete term string** as the key. Different labels with the same GO ID are treated as different terms.
- GO ID extraction and label hiding recognize the format `[GO:digits]`.

---

## 4. Matching and Record-Selection Rules

### 4.1 Gene-name matching

The first column of `T-half` and `UniProtGO.From` are converted to strings and matched by **exact equality**.

- Gene-name case and surrounding whitespace are not normalized automatically.
- Aliases, obsolete names, and strings containing multiple genes are not resolved or split.
- Remove blank gene identifiers before uploading. Missing values may match one another after conversion to the string `nan`.
- Unmatched input rows remain in the output, with empty GO and Entry fields.

### 4.2 Multiple candidates for the same From value

1. Trim whitespace from `Reviewed` and convert it to lowercase.
2. Prefer candidates whose normalized value is exactly **`reviewed`**.
3. If reviewed candidates exist, select among them; otherwise, use all candidates.
4. Convert `Length` to numeric and select the candidate with the largest value.
5. Missing or non-numeric lengths are treated as `-1`.
6. If maximum lengths are tied, select the first candidate in input order.

Example:

| Reviewed | Length | Selected? |
|---|---:|---|
| unreviewed | 900 | No |
| reviewed | 450 | No |
| reviewed | 600 | **Yes** |

**A single record supplies all three GO categories.** CC, BP, and MF are not selected independently. The selected Entry is copied into columns F, G, and H. Missing GO annotations in that record are not filled from other candidates.

---

## 5. Workflow

1. Confirm the input sheet names in the sidebar.
2. Adjust frequency filters and plot settings if required.
3. Upload the Excel workbook.
4. Check the reported row counts for both sheets.
5. Click **“▶ 実行 (Protein list を作成)”** — Run / Create protein list.
6. Review matching counts, unique-term counts, and the previews of the first 50 rows.
7. Download the completed protein-list workbook.
8. Review the CC, BP, and MF boxplots and save the required PNG files.
9. Select a category in the supplementary analysis section and inspect per-term tests, overlap, and conditional comparisons.
10. Download a statistical workbook for each category you want to retain.

The summary metrics labeled as total and matched genes are **input-row counts**, not counts of unique gene identifiers. A row is counted as matched when a candidate exists, even if the selected record contains no GO annotations.

---

## 6. Completed T-half Column Layout

| Column | Output header | Content |
|---|---|---|
| A | `Gene Name` | Original first column |
| B | `HalfLife` | Original second column |
| C | `GO cellular component` | CC from the selected record |
| D | `GO biological process` | BP from the selected record |
| E | `GO molecular function` | MF from the selected record |
| F | `Entry (cellular component)` | Selected Entry |
| G | `Entry (biological process)` | Selected Entry |
| H | `Entry (molecular function)` | Selected Entry |
| I | `Unique GO cellular component` | Unique CC terms |
| J | `Unique GO biological process` | Unique BP terms |
| K | `Unique GO molecular function` | Unique MF terms |
| L | `Count (cellular component)` | Occurrences of the term in column I of the same row |
| M | `Count (biological process)` | Occurrences of the term in column J of the same row |
| N | `Count (molecular function)` | Occurrences of the term in column K of the same row |

### Important: row alignment

- **Columns A–H describe individual input rows.**
- **Columns I–N are independent category-level summaries.**
- A term in column I does not necessarily belong to the gene in column A of that row.
- CC, BP, and MF term lists are also independent of one another.
- Unique terms are listed in **first-occurrence order**, not frequency order.
- The output row count is the maximum of the input row count and the three unique-term counts. Shorter columns are padded with blanks.

Term lists and occurrence counts include annotations from all input rows, including rows with invalid half-lives. In contrast, plots and statistical analyses generally use rows with `HalfLife > 0`. Consequently, counts in columns L–N may differ from plot sample sizes.

---

## 7. Boxplots and Filters

### 7.1 Default settings

| Setting | Default | Behavior |
|---|---:|---|
| Minimum protein count | 20 | Lower bound on the number of eligible half-life observations for a term |
| Maximum fraction | 90% | Upper bound relative to all eligible rows |
| Display limit: Top N | 0 | Zero means all terms; a positive value keeps the N most frequent terms |
| Highlight highest/lowest median terms | 10 | Highest medians are red; lowest medians are blue |
| Plot only highest/lowest median terms | Off | Shows at most twice the highlight count, in descending median order |
| Hide GO IDs | Off | Removes GO IDs from plot labels |
| Show all data points | Off | Overlays observations with reproducible horizontal jitter |

### 7.2 Filter population and order

Let $N$ be the total number of eligible rows, $n_t$ the number of half-life observations for a term, $n_{min}$ the minimum count, and $f_{max}$ the maximum fraction. A term passes when:

$$
n_{min}\leq n_t\leq f_{max}N
$$

- $N$ includes **all rows whose numeric half-life is greater than zero**, including unmatched rows and rows without annotations in the selected category.
- Frequency conditions are applied first, followed by frequency sorting and the Top N limit. Median-based selection is performed afterward.
- Standard plots use frequency order. Highest/lowest-only plots use descending median order for the selected terms.
- If the highest and lowest groups overlap, the red highest-median group takes priority. The blue group may therefore contain fewer terms than requested.
- When the highlight count is zero, enabling highest/lowest-only mode does not restrict the plotted terms.

### 7.3 Plot specifications

- Y-axis: `T1/2 (days)`, on a **logarithmic scale**.
- Red: highest-median terms; blue: lowest-median terms; light purple: other terms.
- Boxes: first to third quartile; center line: median.
- Whiskers: Matplotlib's standard 1.5 IQR rule.
- Outliers are shown even when all-point display is disabled.
- PNG export resolution: **150 dpi**.
- All-point jitter uses `random.Random(0)` for reproducible positions.

Only the display axis is logarithmic. The code does not log-transform half-lives before statistical testing.

---

## 8. Supplementary Analysis: Per-Term Tests

For each term passing the frequency filters and Top N limit in the selected category, the application compares:

- **Annotated group:** eligible rows carrying the term.
- **Non-annotated group:** all other eligible rows.

The non-annotated group includes unmatched rows and rows without any annotation in that category. Absence of an annotation does not establish absence of the biological function.

### 8.1 Test implementation

The application implements a **two-sided Mann–Whitney U test** with:

- Normal approximation of the rank statistic.
- Tie correction.
- No continuity correction.
- No exact-test calculation.
- Missing P values when either group is empty or the corrected variance is non-positive.

This is a rank-based comparison of distributions. Do not interpret the test solely as a comparison of the displayed medians.

### 8.2 Effect size: Cliff's delta

With the annotated group denoted by $x$ and the non-annotated group by $y$, the implementation calculates:

$$
\delta=\frac{2U_x}{n_xn_y}-1
$$

- Positive: annotated observations tend to have longer half-lives.
- Negative: annotated observations tend to have shorter half-lives.
- Near zero: little net rank-based dominance.

The results table is sorted by descending delta, not by ascending P value or adjusted P value.

### 8.3 Multiple-testing adjustment

P values are adjusted using **Benjamini–Hochberg (BH)** and reported in the column `q (BH)`.

**The adjustment covers the terms tested in the currently selected category under the current frequency filters and Top N setting.** It does not combine all three categories into a single testing family. Changing Top N may change both the tested terms and their adjusted values.

Restricting the plot to highest/lowest-median terms does not restrict testing to that visual subset. Tests still cover all terms passing frequency filtering and Top N.

The interface counts terms with `q < 0.05` and positive or negative delta. Red and blue boxplot highlights indicate median ranking, **not statistical significance**.

---

## 9. Overlap Between Terms

For highest- and lowest-median terms, the application forms **sets of gene names** carrying each annotation and calculates:

$$
\mathrm{Jaccard}(A,B)=\frac{|A\cap B|}{|A\cup B|}
$$

$$
\mathrm{Overlap}(A,B)=\frac{|A\cap B|}{\min(|A|,|B|)}
$$

For example, if $|A|=20$, $|B|=40$, and the intersection contains 10 genes, Jaccard is 0.20 and the overlap coefficient is 0.50.

- The workbook includes a matrix and a pairwise table.
- The interface displays the first 30 pairs in descending Jaccard order.
- If the highlight count is zero, overlap analysis still uses the **10 highest- and 10 lowest-median terms**.
- Set operations remove duplicate gene names. This differs from the input-row basis used for boxplots and tests, so duplicate genes can produce different counts across analyses.

**Known limitation:** If only one term is selected, pairwise-table construction fails because the resulting table has no `Jaccard` column. Adjust filters so that at least two terms remain when using supplementary analysis.

---

## 10. Conditional Comparisons

Within the selected category, the application removes all rows carrying an exclusion GO ID **from the entire analysis population**, then compares rows with and without the target GO ID.

Defaults:

| Role | GO ID | Description in the code |
|---|---|---|
| Target term | `GO:0043209` | myelin sheath |
| Exclusion term | `GO:0005739` | mitochondrion |

Exclusion applies to both the target-annotated group and its comparison group. This is a subset-based sensitivity analysis, not covariate adjustment in a regression model.

Outputs include the target count before exclusion, overlap count and fraction, target-annotated count after exclusion, both group medians, delta, and P value.

- The target and exclusion terms need not pass the ordinary frequency or Top N filters.
- IDs are matched within the selected category. Reconsider the input IDs when changing categories.
- If an annotation lacks a bracketed GO ID, its entire term string is used as the matching key instead.
- BH adjustment is not applied to the conditional-comparison P value.
- A comparison is not valid when exclusion leaves an empty group, including when target and exclusion IDs are identical.

---

## 11. Downloaded Files

### 11.1 Annotated protein list

Filename: `ProteinList_GO.xlsx`

| Sheet | Content |
|---|---|
| `T-half` | Completed columns A–N |
| `UniProtGO` | Table data read from the source UniProtGO sheet |

Output sheet names remain fixed even if custom input sheet names are used. This is a new workbook reconstructed from table data, not an in-place workbook edit preserving other sheets, original formatting, formula structures, or macros.

### 11.2 Boxplots

| Category | Filename |
|---|---|
| CC | `boxplot_cellular_component.png` |
| BP | `boxplot_biological_process.png` |
| MF | `boxplot_molecular_function.png` |

### 11.3 Statistical results

The filename includes the category, for example: `GO_term_statistics_cellular_component.xlsx`.

| Sheet | Content |
|---|---|
| `Term tests` | Per-term sample counts, medians, delta, P values, and BH-adjusted values |
| `Jaccard matrix` | Jaccard matrix for selected terms |
| `Jaccard pairs` | Pairwise set sizes, intersection sizes, Jaccard, and overlap coefficients |
| `Conditional test` | Target and exclusion IDs and comparison results |
| `Settings` | Category, frequency thresholds, Top N, overlap-selection count, test method, and population size |

Although the download button displays a generic filename, the actual filename includes the category. Each export covers the currently selected category; switch categories and save separately to retain all three.

---

## 12. Data Quality and Implementation Notes

### Duplicate terms within a cell

Remove duplicates such as `term A; term A` before uploading.

- The unique-term list contains the term once, but occurrence counts and boxplot half-life observations count it twice.
- Per-term tests use set-based membership and count that row once.
- Consequently, frequency filtering and plot counts may disagree with statistical sample counts.

### Invalid half-life values

Missing, non-numeric, and non-positive values are excluded from plots and statistics but remain in the completed workbook. Positive infinity is not explicitly filtered out, so prepare **finite, strictly positive numeric values** for analysis.

### Empty T-half sheet

A sheet with two headers but no data rows may fail during processing. Provide at least one valid input row.

### Preview limitation

The I–N preview retains only rows where **column I, the CC term list, is non-empty**. Additional BP or MF terms beyond the length of the CC list may therefore be absent from the preview. Check the exported workbook for the complete lists.

### Session state and reruns

- Matching results are retained in the session. **Click Run again** after changing sheet names or source data.
- File-change detection uses filename and size. A replacement with the same filename and size may not automatically clear previous results.
- Plot filters and display settings apply to the saved matching result. PNGs are reused or regenerated according to their settings.
- If a new processing attempt fails while an earlier result is stored, resolve the error and rerun. Do not mistake an earlier result for the failed new run.

---

## 13. Troubleshooting

| Problem | Check or action |
|---|---|
| Sheet not found | Check the displayed sheet list and enter exact names in the sidebar |
| Required columns missing | Check all seven UniProtGO headers, including whitespace, punctuation, and spelling |
| Low matching rate | Check gene-name whitespace, case, aliases, and concatenated identifiers |
| Empty GO annotations | Check whether the selected reviewed/longest record actually contains annotations |
| No terms displayed | Lower the minimum count, increase the maximum fraction, and check the number of positive half-lives |
| Jaccard error in supplementary analysis | Ensure at least two terms remain after filtering |
| Zero rows in conditional comparison | Check the selected category, IDs, and bracketed GO ID format |
| Missing P values or adjusted values | Check for empty groups or zero rank variance, such as identical values throughout |
| Plot too wide | Use Top N or highest/lowest-median-only display |
| Japanese X-axis title missing | The application omits that title when no supported Japanese font is detected; English GO labels use the default font |
| Macros absent from output | Output is a new `.xlsx`; retain the original `.xlsm` separately |

---

## 14. Information to Retain When Sharing Results

For reproducibility and interpretation, save the following alongside the output files:

- Input workbook and the source, retrieval date, and organism of the UniProtGO annotations.
- Input sheet names, duplicate-gene and missing-value handling, and half-life units.
- Minimum count, maximum fraction, Top N, and analysis category.
- Median-highlight count and plot-display options.
- The term family included in BH adjustment.
- Target and exclusion GO IDs for conditional comparisons.
- Python and dependency versions.

GO terms frequently share the same proteins. Differences observed for several terms are not necessarily independent evidence. Interpret effect sizes, distributions, overlap, and annotation coverage together.
