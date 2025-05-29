# proofreading_hpc_simulation_and_plot
this folder is for running simulations on Westlake HPC with multithreading and LHS enabled.

Here below is the conversion table for the naming of chemical species in code and article.
| Naming in Code    | Naming in Article |
| --------          | -------           |
| A                 | R                 |
| B                 | W                 |
| C                 | E                 |
| R                 | Recptor           |
| AC                | ER                |
| BC                | EW                |
| AR                | not applicable    |
| BR                | not applicable    |

## file list
- **create_parameters.py**: create parameters used for simulation, output into a folder of csvs, each csv contains a chunk of parameters for multiple runs in single thread , which speed up simulation.
- **run_simulation.py**: actual simulation in here, read a folder of parameters, output into **parameters.csv**, **with_feedback.csv**, **without_feedback.csv**, use filelock to ensure parameters and result is synchronized.
- **workflow.sh**: wrapper for **run_simulation.py**, for submitting to HPC with following settings.
```shell
#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH -q huge
#SBATCH -J proofreading_simulation
#SBATCH -c 1
#SBATCH -a 1-800
#SBATCH -o log/array-job.%A.%a.log
```
- **simulation_plot.ipynb**: analyzing simulation result and ploting.
- **rd_solver.py**: adapted from Zhuo & Yuanqi with modification.

## Licensing and Attribution

The **sankey_mod.py** file is modified from project [**mpl_sankey**](https://github.com/toobaz/mpl_sankey), originally licensed under the LGPL-3.0. The modified file is also released under the LGPL-3.0.

### Original File

- **File:** __init__.py
- **Author:** Pietro Battiston <me@pietrobattiston.it>
- **License:** LGPL-3.0

### Modified File

- **File:** sankey_mod.py
- **Modifier:** Yuxuan Yan <yanyuxuan@westlake.edu.cn>
- **Date:** 2024/12/27
- **Modifications:** modified display style
- **License:** LGPL-3.0

The **bounded_1d_kde.py** file is modified from project [**pesummary**](https://github.com/pesummary/pesummary), originally licensed under the MIT. The modified file is also released under the MIT.

### Original File

- **File:** bounded_1d_kde.py
- **Author:** Charlie Hoy <charlie.hoy@ligo.org> & Michael Puerrer <michael.puerrer@ligo.org>
- **License:** MIT

### Modified File

- **File:** sankey_mod.py
- **Modifier:** Yuxuan Yan <yanyuxuan@westlake.edu.cn>
- **Date:** 2025/03/21
- **Modifications:** make the file can be used individually and do not depend on the entire pesummary package.
- **License:** MIT