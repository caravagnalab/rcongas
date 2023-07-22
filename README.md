
# Rcongas <a href='caravagnalab.github.io/rcongas'><img src='man/figures/logo.png' align="right" height="120" /></a>

<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

This is the `(R)CONGAS+` R package, an interface to run the
probabilistic methods implemented in the Python
[CONGAS+](https://github.com/caravagnalab/CONGASp) package using Pyro. These
methods implement several statistical models to genotype Copy Number
Alterations from single-cell RNA and ATAC sequencing, integrating at the
same time bulk DNA sequencing.

This package implements S3 objects to preprocess and visualize input
single-cell data, create and visualize model fits. The current package
is an extension of the original, single-molecule version; it retains the
same name to avoid confusion.

#### Citation

[![](https://img.shields.io/badge/doi-10.1101/2021.02.02.429335-red.svg)](https://doi.org/10.1101/2023.04.01.535197)

If you use `Rcongas`, please cite these two papers:
- *A Bayesian method to cluster single-cell RNA sequencing data using copy number alterations*. Salvatore Milite, Riccardo Bergamin, Lucrezia Patruno, Nicola Calonaci, Giulio Caravagna. [Bioinformatics](https://doi.org/10.1093/bioinformatics/btac143) 2022.

- *A Bayesian method to infer copy number clones from single-cell RNA and ATAC sequencing.* Lucrezia Patruno1, Salvatore Milite, Riccardo Bergamin Nicola Calonaci, Alberto D’Onofrio, Fabio Anselmi, Marco Antoniotti, Alex Graudenzi, Giulio Caravagna.
    [biorXiv
    preprint](https://www.biorxiv.org/content/10.1101/2023.04.01.535197v1),
    2023



#### Help and support

[![](https://img.shields.io/badge/GitHub%20Pages-https://caravagnalab.github.io/rcongas/-steelblue.svg)](https://caravagnalab.github.io/rcongas)

### Installation

You can install the released version of `Rcongas` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("caravagnalab/Rcongas")
```

------------------------------------------------------------------------

#### Copyright and contacts

Salvatore Milite, MSc, Riccardo Bergamin, PhD, Nicola Calonaci, PhD and Giulio Caravagna,
PhD. *University of Trieste, Trieste, Italy*.
Lucrezia Patruno, PhD. *UCL Cancer Institute, London, UK*
Alex Graudenzi, PhD, Prof. Marco Antoniotti *University of Milano-Bicocca, Milan, Italy*

[![](https://img.shields.io/badge/Email-gcaravagn@gmail.com-steelblue.svg)](mailto:gcaravagn@gmail.com)
[![](https://img.shields.io/badge/CDS%20Lab%20Github-caravagnalab-seagreen.svg)](https://github.com/caravagnalab)
[![](https://img.shields.io/badge/CDS%20Lab%20webpage-https://www.caravagnalab.org/-red.svg)](https://www.caravagnalab.org/)
