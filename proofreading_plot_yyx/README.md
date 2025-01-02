# proofreading_plot_yyx

## file list
- **create_parameters.py**: create parameters used for simulation, output into a folder of csvs, each csv contains a chunk of parameters for multiple runs in single thread , which speed up simulation.
- **figure_2_lhs.py**: actual simulation in here, read a folder of parameters, output into **parameters.csv**, **with_feedback.csv**, **without_feedback.csv**, use filelock to ensure parameters and result is synchronized.
- **workflow.sh**: wrapper for **figure_2_lhs.py**, for submitting in HPC.
```shell
#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH -q huge
#SBATCH -J proofreading_simulation
#SBATCH -c 1
#SBATCH -a 1-800
#SBATCH -o log/array-job.%A.%a.log
```
- **figure_2_plot.ipynb**: analyzing simulation result and ploting.
- **rd_solver.py**: adapted from Zhuo & Yuanqi

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