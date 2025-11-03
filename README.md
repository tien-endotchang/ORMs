# Understanding and Using the Relative Importance Measures Based on Orthogonalization and Reallocation

**R Version:** 4.4.1 (2024-06-14 ucrt)

This repository contains the R scripts used for our comprehensive Monte Carlo simulation studies on ORMs (Orthogonalization-Reallocation Measures).  
All scripts and their purposes are summarized below.

---

## Scripts Overview

1. `sim.R`
- Performs comprehensive Monte Carlo simulations across all **four orthogonalization** and **four reallocation** methods.  
- **Output:** Automatically creates a `res/` folder containing the simulation results.

### `plot.R`
- Generates **Figures 4 – 10** and **Table 1** from the paper.  
- **Output:** Plots and tables are saved to the `fig/` directory.

### `example.R`
- Reproduces **Tables 3 and 4**.  
- **Output:** Tables are saved to the `fig/` directory.

### `utils_sim.R`
- Helper functions for the simulation routines (`sim.R`).

### `utils_plot.R`
- Helper functions used for generating and formatting plots (`plot.R`).

---

## Running the Simulation

To fully reproduce the results from our study, run the scripts in the following order:

```r
source("sim.R")
source("plot.R")
source("example.R")
```

## Using Saved Data
If you prefer not to rerun the full simulation:
1. Use the aggregated data already provided in the `res/` folder.
2. Run only:

```r
source("plot.R")
```

This will generate all figures and tables directly from the pre-computed data.

---
## Path and Encoding Notice
Please avoid using excessively long paths or non-ASCII characters (e.g., Mandarin) in your working directory names.
Such paths may cause I/O errors (e.g., `Invalid argument`) when reading or writing `.rds` files on Windows systems.
