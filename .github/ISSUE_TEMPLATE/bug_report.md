---
name: Bug report
about: Something does not work as documented
labels: bug
---

**What happened**

A clear description of the problem and of what you expected instead.

**Reproducible example**

```r
library(lagprotest)
sim <- simulate_lagged_series(n_times = 12, true_lag = 2, seed = 1)
# code that shows the problem
```

**Session info**

<details>

```r
sessionInfo()
```

</details>
