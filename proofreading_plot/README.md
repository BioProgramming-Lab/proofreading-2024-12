# proofreading_plot
This folder is for running simulations for figures, excluding those running on HPC

## IMPORTANT!!!
This part of the project uses branches for different plots.
- **main**: simulation and plot with standard extracellular proofreading circuit.
- **detailed-formula-for-figS5**: addition of a detailed formulation of eta (fidelity) in classic kinetic proofreading, in response to Document S2\

# For figure3b&3e&4d: Retrieving Simulation Results and Analyzing Gradients
. This code retrieves results from simulations and computes gradients for different parameter combinations.  
. It enables exploration of gradients for specific results, aiding in understanding their trends and behavior.  
. The outputs are not the same as the figures in the referenced paper, but the observed trends remain consistent.

Below is the conversion table for the naming of chemical species in code and the article.
| Naming in Code    | Naming in Article |
| --------          | -------           |
| A                 | R                 |
| B                 | W                 |
| C                 | E                 |
| R                 | Receptor          |
| AC                | ER                |
| BC                | EW                |
| AR                | not applicable    |
| BR                | not applicable    |
| _ac               | self-activation   |
| _rp               | mutual inhibition |
| bassline          | without proofreading|
