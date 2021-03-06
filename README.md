
<!-- README.md is generated from README.Rmd. Please edit that file -->

# poisrrr

The package poisrrr implements the Poisson Reduced Rank method
introduced in [C. Jentsch, E. R. Lee and E. Mammen (2020) Time-dependent
Poisson reduced rank models for political text data analysis.
Computational Statistics and Data Analysis, 142,
106813](https://doi.org/10.1016/j.csda.2019.106813). See also [C.
Jentsch, E. Mammen and E. R. Lee (2021) Poisson reduced rank models with
an application to political text data. Biometrika, 108, 2, 455 -
468](https://doi.org/10.1093/biomet/asaa063)

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("chroetz/poisrrr")
```

## Example

Create a Term-Document-Matrix form quanteda’s inaugural address corpus.

``` r
library(magrittr)
quanteda::data_corpus_inaugural %>% 
  quanteda::tokens(remove_punct = TRUE, remove_symbols = TRUE, remove_numbers = TRUE) %>% 
  quanteda::dfm(verbose = FALSE) %>% 
  as.matrix() %>% 
  t() ->
  tdm
tdm <- tdm[rowSums(tdm) > 5, ] # remove rare words
```

Apply the method for *K* = 2 dimensions and plot the resulting plane
with document positions.

``` r
library(poisrrr)
K <- 2
theta <- estim(tdm, K, verbose=FALSE)
lst <- theta2plist(theta, K)
v <- lst$v
plot(NA, xlim=c(-0.2, 0.22), ylim=c(-0.3, 0.25), xlab="Dimension 1", ylab="Dimension 2")
points(v)
lines(v)
labels <- rownames(v)
labels[-c(1,4,7,19,20,22,29,32,33,37,38,39,40,46,55,56,58)] <- NA
text(v, labels=labels, cex=0.8, pos=3, offset=0.2)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />
