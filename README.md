# Rethinking Parameter Invariance in Cognitive Modeling

Code for the paper

> Yu, K., & Robinson, M. M. (2026). *Rethinking parameter invariance in cognitive modeling* [Manuscript submitted for publication]. KU Leuven.

DOI: [10.5281/zenodo.23733805](https://doi.org/10.5281/zenodo.23733805)

The repository contains the code for the formal examples and the simulations in
the paper, the scripts that draw each figure of the paper and the supplement, and
the saved results of the simulations. The supplementary material itself is in
`supplement.pdf`. `code/README.md` lists which script produces which figure and
result, and `code/DATA_DICTIONARY.md` describes the result files.

All commands are run from the repository root.

## Reproducing the paper

The analyses were run with R 4.5.2. The first two commands install the recorded
package versions into a project library; the third rebuilds the figures.

```sh
Rscript -e 'install.packages("renv")'
Rscript --vanilla -e 'renv::restore()'
Rscript --vanilla code/run_all.R
```

The last command takes about ten seconds. It reruns Simulation 1 and the
supplementary checks, draws the other figures from the saved simulation results,
and compares every saved result with the checksums in `code/results/manifest.csv`.

To run Simulations 2A and 2B and the finite-validation analysis again, which takes
about eight minutes:

```sh
Rscript --vanilla code/run_all.R --resimulate-all
```

## License

MIT, see `LICENSE`.
