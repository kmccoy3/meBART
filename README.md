
# Bayesian Additive Regression Trees with Measurement Error (meBART)


$y_i = f(X_i) + \varepsilon_i, ~~~~~~~~ \varepsilon_i \sim \mathcal{N}(0, \sigma^2)$

$X_i^* = X_i + e_i$

# Installation Instructions

```{bash}
install.packages("remotes")
remotes::install_github("kmccoy3/meBART")
```

# Notes

To add manual:
`devtools::build_manual()`

To re-do Roxygen documentation:
`devtools::document()`

Add documentation for whole package:
`use_package_doc()`