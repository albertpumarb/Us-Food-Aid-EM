# US Food Aid and Civil Conflict

Problem Set 3 replication of Nunn and Qian (2014), American Economic Review
104(6), 1630–1666. We reproduce Table 2, column 6, and the recipient-group
comparison in Figures 3 and 4. The extension adds country-specific linear
trends to both stages of the baseline IV model.

## Reproduce the results

Run from the repository root with R 4.2.2 and its bundled recommended package
`foreign` 0.8-83. No additional packages are required.

    Rscript --vanilla replicate.R

`replicate.R` reads `data/FAid_Final.dta`, constructs the analysis sample and
controls, estimates both specifications, and checks the baseline against Table 2.
It creates the following report inputs and supporting files under `results/`.

| Output | Contents |
| --- | --- |
| `main.tex` | Published and replicated estimates for Table 2, column 6 |
| `extension.tex` | Baseline and country-trend estimates and first-stage diagnostics |
| `conflict_and_wheat.pdf` | Two-panel reconstruction of Figures 3 and 4 |
| `estimates.csv`, `figure_data.csv`, `session.txt` | Numeric results, plotted data, and R session details |

## Data source

`data/FAid_Final.dta` is the unmodified analysis input from the
[authors' replication package](https://www.aeaweb.org/articles?id=10.1257/aer.104.6.1630),
available under **Additional Materials**. Dataset copyright 2014 American
Economic Association, licensed under
[Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/).
The full archive, original Stata code, and original output are available from
the publisher.
