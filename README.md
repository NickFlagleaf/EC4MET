# EC4MET <a href="https://nickflagleaf.github.io/EC4MET/index.html"><img src="man/figures/logo.png" align="right" height="200"/></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/NickFlagleaf/EC4MET/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/NickFlagleaf/EC4MET/actions/workflows/R-CMD-check.yaml)

[![DOI](https://zenodo.org/badge/938421607.svg)](https://doi.org/10.5281/zenodo.15400244)

<!-- badges: end -->

## Derive ECs for METs

EC4MET derives environmental data and defines covariates (ECs) for crop trial environments in Australia to model environmental effects and genotype by environment interactions.

See details published in [Prediction of Australian wheat genotype by environment interactions and mega-environments](https://doi.org/10.1007/s00122-025-05023-6).


An example workflow is detailed in the [vignette](https://nickflagleaf.github.io/EC4MET/articles/EC4MET-workflow-example.html).

------------------------------------------------------------------------

### Installation

To install from Github using devtools:

```         
devtools::install_github("NickFlagleaf/EC4MET",build_vignettes = T)
```

------------------------------------------------------------------------


### Data licenses and use

EC4MET uses freely available weather and soil data from multiple sources including [SILO](https://www.longpaddock.qld.gov.au/silo/),
[BARRA-R2](https://opus.nci.org.au/spaces/NDP/pages/264241166/BOM+BARRA2+ob53), [SLGA](https://esoil.io/TERNLandscapes/Public/Pages/SLGA/index.html),
and [CMIP6 QDC](https://doi.org/10.25919/03by-9y62) which are all licensed under Creative Commons Attribution 4.0 International ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)).