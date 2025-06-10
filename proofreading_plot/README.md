# proofreading_plot
This folder is for running simulations for figures exculding those running on HPC

## IMPORTANT!!!
This part of the project uses branches for different plot.
- **main**: simulation and plot with standard extracellular proofreading circuit.
- **detailed-formula-for-figS5**: addition in the formulation of eta (fidelity) in classic kinetic proofreading, in respond to Document S2\

# For figure3b&3e&4d: Retrieving Simulation Results and Analyzing Gradients
. This code retrieves results from simulations and computes gradients for different parameter combinations.  
. It enables exploration of gradients for specific results, aiding in understanding their trends and behavior.  
. The output are not the same as figures in the referenced paper, but the observed trends remain consistent.

Here below is the conversion table for the naming of chemical species in code and article.
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
| bassline          | withtout proofreading|
