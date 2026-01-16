
# measurement error Bayesian additive regression trees (meBART)

We propose measurement error Bayesian additive regression trees (meBART), a novel extension to the BART model ([Chipman et al., 2010](https://doi.org/10.1214/09-AOAS285)) that directly incorporates measurement error in the independent variable(s). BART assumes that data $\{(X_i, y_i); i=1,\ldots, n\}$ is generated according to some unknown function $f$ that can be modeled as a sum of decision trees:

$$y_i = f(X_i) + \varepsilon_i, ~~~~~~~~ \varepsilon_i \sim \mathcal{N}(0, \sigma^2)$$

$$f(X_i) = \sum\limits_{h=1}^m g(X_i; T_h, M_h)$$

meBART, however, represents the true unobserved feature values $X_i$ as latent variables in a Bayesian hierarchical model, and uses these latent variables as inputs to the prediction model. Our observed data is thus $\{(X_i^*, y_i); i=1,\ldots, n\}$, where $X_i^*$ is a noisy realization of $X_i$. We assume that the measurement error is additive and normal:

$$X_i^* = X_i + e_i ~~~~~~~~ e_i \sim \mathcal{N}(0, \sigma_e^2)$$


This code repository accompanies the preprint titled "Tree-Based Prediction Models for Noisy Input Data" (McCoy et al., 2026). The manuscript is currently being prepared for submission to a suitable journal. This README will be updated as the article gets closer to publication. That being said, the source code contained in this repository is complete, and the code used to generate the plots contained in the manuscript can be found in `/vignettes/`.



# Installation Instructions

In order to download and use meBART, run the R code below:

```{r}
install.packages("remotes")
remotes::install_github("kmccoy3/meBART")
```

# Relevant Works

The original BART article:
> Chipman, H. A., George, E. I., & McCulloch, R. E. (2010). BART: Bayesian additive regression trees. *The Annals of Applied Statistics*, 4(1), 266 - 298

The R BART package, on which the meBART package is built:
> Sparapani, R., Spanbauer, C., & McCulloch, R. (2021). Nonparametric machine learning and efficient computation with Bayesian additive regression trees: The BART R package. *Journal of Statistical Software*, 97, 1-66.

# Citation

To cite this repository, simply cite our manuscript:

```{TeX}
@article{mccoy2025tree,
  title={Tree-Based Prediction Models for Noisy Input Data},
  author={McCoy, Kevin and Wooten, Zachary and Peterson, Christine B.},
  journal={TBD},
  volume={TBD},
  pages={TBD},
  year={2026},
  publisher={TBD}
}
```

