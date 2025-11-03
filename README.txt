R Version: 4.4.1 (2024-06-14 ucrt)

This repository contains the R scripts used for our simulation study. The scripts and their functions are detailed below:

**Scripts:**

*   `sim.R`
    *   Performs comprehensive Monte Carlo simulations for all four orthogonalization methods and all four reallocation methods.
    *   **Output:** Creates a `res` folder containing simulation results.

*   `plot.R`
    *   Generates Figures 4-10 and Table 1.
    *   **Output:** Plots and table saved to the working directory.

*  `example.R`
     *  Generates Tables 3 and 4.
     *  **Output:** Tables saved to the working directory.

*   `utils_sim.R`
    *   Contains functions used within the simulation scripts and related tasks.

*   `utils_plot.R`
    *   Contains functions for generating the figures.

**Running the Simulation:**

To fully replicate our simulation study, run the scripts in the following order:

1.  `sim.R`
2.  `plot.R`
3. `example.R`

**Using Pre-Aggregated Data:**

The aggregated data used in our study is located in the `res` folder. 
To reproduce the figures, begin by running the `plot.R` script, using the data contained within that folder.

**Path**
Please avoid too long path name (or Mandarin characters in the path)