import time
import numpy as np
from rd_solver import *
import os
import zarr
from shared import *


diff_coeffs = DIFFUSION()

# read pre-generated parameters
task_id = int(os.getenv("SLURM_ARRAY_TASK_ID"))
if task_id is None:
    raise Exception(
        "Unable to find environment variable SLURM_ARRAY_TASK_ID"
    )

# read parameter
store = "parameters_20250225"
z = zarr.open(
    store=store,
    mode="r",
)
parameters = z[
    ((task_id - 1) * chunk_size):(task_id * chunk_size)
]

# open result zarr
store = "result_20250225"
z = zarr.open(
    store=store,
    mode="a",
)

   

results = {}
for i, parameter in enumerate(parameters.reshape((-1, num_dimensions))):
    # randomize parameters
    D_0 = parameter[0]
    k_kinase = parameter[1]
    k_phosphotase = parameter[2]
    c_0 = parameter[3]

    # run with different beta (receiver region j_A0 factor)
    for j, meta_parameter in enumerate(meta_parameters):
        c_0_tuple = (
            # c_A
            np.full(n_gridpoints, c_0),
            # c_B
            np.full(n_gridpoints, c_0),
            # c_C
            np.full(n_gridpoints, c_0),
            # c_Ap
            np.zeros(n_gridpoints),
            # c_Bp
            np.zeros(n_gridpoints),
            # c_CAp
            np.zeros(n_gridpoints),
            # c_CBp
            np.zeros(n_gridpoints),
        )

        # * run simulation
        # diffusion rate of each molecule
        diff_coeffs.D_A = D_0  # um^2/s
        diff_coeffs.D_B = D_0  # um^2/s
        diff_coeffs.D_C = D_0  # um^2/s
        diff_coeffs.D_Ap = D_0  # um^2/s
        diff_coeffs.D_Bp = D_0  # um^2/s
        diff_coeffs.D_complex = D_0  # um^2/s

        rxn_params = RXN_params_yuanqi(
            k_kin=k_kinase,
            k_p=k_phosphotase,
        )

        derivs_0 = np.array([
            -rxn_params.k_kin,
            -rxn_params.k_kin,
            0,
            rxn_params.k_kin,
            rxn_params.k_kin,
            0,
            0
        ])

        res = np.array(RD_solve(
            c_0_tuple, t, L=L, derivs_0=derivs_0, derivs_L=0,
            diff_coeff_fun=Diff_fun, diff_coeff_params=(diff_coeffs,),
            rxn_fun=RD_rxn, rxn_params=(rxn_params,),
            rtol=1e-6, atol=1e-6
        ))

        results[(i, j)] = res


# storage configuration
# root
# ├─ meta_parameter
# |  ├─ time
# |  |  ├─ chemical species
# |  |  |  n_sim X n_grid (chunk: chunk_size X n_grid)
# |  |  |
# ...


start = time.time()
for j, meta_parameter in enumerate(meta_parameters):
    group_meta_parameter = z[meta_parameter.as_string()]
    # group_meta_parameter.attrs.update(meta_parameter)
    for _t in t_record:
        array_s_sim_l = group_meta_parameter[str(_t)]

        chunk_s = []

        for s in Species:
            chunk = []
            for i, parameter in enumerate(parameters.reshape((-1, num_dimensions))):
                chunk.append(results[(i, j)][s.value, list(t).index(_t), :])

            chunk_s.append(chunk)
        array_s_sim_l[
            :, ((task_id - 1) * chunk_size):(task_id * chunk_size), :
        ] = chunk_s

print(time.time() - start)
