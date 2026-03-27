# fddml

**Fuzzy Difference-in-Differences with Machine Learning**

R package implementing the DML-Wald and DML-TC estimators from Mena (2026), "Double Debiased Machine Learning for Difference-in-Differences under Imperfect Compliance."

## Installation

```r
# install.packages("remotes")
remotes::install_github("andres-mena/fddml")
```

## Quick Start

```r
library(fddml)

# Bundled INPRES application (Duflo 2001, CD'H 2018)
data(duflo)
fit <- fddml(duflo$Y, duflo$D, duflo$G, duflo$Ti, duflo$X,
             estimand = "both", method = "lasso", K = 5,
             cluster = duflo$cluster)
print(fit)
summary(fit)
plot(fit)
```

## Features

- **Two estimands**: DML-Wald (Assumptions 1-6) and DML-TC (Assumptions 1-3, 4')
- **Flexible ML**: Lasso, Random Forest, neural networks, OLS, or custom functions
- **Cross-fitting**: K-fold sample splitting for valid post-selection inference
- **Optimal trimming**: Data-driven propensity score trimming (Theorem 2)
- **Inference**: Analytical, cluster-robust, or multiplier bootstrap SEs
- **Bundled data**: INPRES school construction application

## Modular API

```r
# Step-by-step for advanced users
nuis <- fddml_nuisance(Y, D, G, Ti, X, method = "rf", K = 5)
trim <- fddml_trim(nuis$pG_raw, method = "auto")
wald <- fddml_wald(W, nuis)
inf  <- fddml_inference(wald, cluster = district_id, se_type = "cluster")
```

## References

- Mena, A. (2026). Double Debiased Machine Learning for DID under Imperfect Compliance. Brown University.
- De Chaisemartin, C. and D'Haultfoeuille, X. (2018). Fuzzy Differences-in-Differences. *Review of Economic Studies*, 85(2).
- Chernozhukov, V. et al. (2018). Double/Debiased Machine Learning. *Econometrics Journal*, 21(1).

## License

MIT
