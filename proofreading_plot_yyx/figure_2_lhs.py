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
store = "parameters_20250218_jR"
z = zarr.open(
    store=store,
    mode="r",
)
parameters = z[
    ((task_id - 1) * chunk_size):(task_id * chunk_size)
]

# open result zarr
store = "result_20250320_sami_error_rate"
z = zarr.open(
    store=store,
    mode="a",
)


results = {}
for i, parameter in enumerate(parameters.reshape((-1, num_dimensions))):
    # randomize parameters
    D_0 = parameter[0]
    j_A0 = parameter[1]
    j_R0 = parameter[2]
    koff_AR0 = parameter[3]
    gamma0 = parameter[4]

    # run with different beta (receiver region j_A0 factor)
    for j, meta_parameter in enumerate(meta_parameters):
        c_0_tuple = (
            # c_A
            np.zeros(meta_parameter["n_gridpoints"]),
            # c_B
            np.zeros(meta_parameter["n_gridpoints"]),
            # c_C
            np.zeros(meta_parameter["n_gridpoints"]),
            # c_R
            np.zeros(meta_parameter["n_gridpoints"]),
            # c_AC
            np.zeros(meta_parameter["n_gridpoints"]),
            # c_BC
            np.zeros(meta_parameter["n_gridpoints"]),
            # c_AR
            np.zeros(meta_parameter["n_gridpoints"]),
            # c_BR
            np.zeros(meta_parameter["n_gridpoints"]),
        )

        # * run simulation
        # diffusion rate of each molecule
        diff_coeffs.D_A = D_0
        diff_coeffs.D_B = D_0
        diff_coeffs.D_C = D_0
        diff_coeffs.D_complex = D_0

        # j_A: sender_region secretion of free A
        j_A = np.zeros(meta_parameter["n_gridpoints"])
        j_A[0:meta_parameter["sender_region"]] = j_A0 * \
            meta_parameter["sender_ratio"]
        # j_B: sender_region secretion of free B
        j_B = np.zeros(meta_parameter["n_gridpoints"])
        j_B[0:meta_parameter["sender_region"]] = j_A0 * \
            meta_parameter["sender_ratio"]
        # j_C: sender_region secretion of free C
        j_C = np.zeros(meta_parameter["n_gridpoints"])
        j_C[0:meta_parameter["sender_region"]] = j_A0 * 2 * \
            meta_parameter["sender_ratio"]
        # # j_R: the secretion rate of free receptor
        j_R = np.zeros(meta_parameter["n_gridpoints"])
        j_R[meta_parameter["n_gridpoints"] -
            meta_parameter["receiver_region"]:] = j_R0

        # self_activation production of free A
        j_self_activation_ar_on_a = np.zeros(meta_parameter["n_gridpoints"])
        j_self_activation_ar_on_a[
            meta_parameter["n_gridpoints"] - meta_parameter["receiver_region"]:
        ] = j_A0 * 0
        # self_activation production of free B
        j_self_activation_br_on_b = np.zeros(meta_parameter["n_gridpoints"])
        j_self_activation_br_on_b[
            meta_parameter["n_gridpoints"] - meta_parameter["receiver_region"]:
        ] = j_A0 * 0

        # self_activation production rate of A+C
        j_self_activation_ac_on_ac = np.zeros(meta_parameter["n_gridpoints"])
        j_self_activation_ac_on_ac[
            meta_parameter["n_gridpoints"] - meta_parameter["receiver_region"]:
        ] = j_A0 * meta_parameter["b_ac"]
        # self_activation production rate of B+C
        j_self_activation_bc_on_bc = np.zeros(meta_parameter["n_gridpoints"])
        j_self_activation_bc_on_bc[
            meta_parameter["n_gridpoints"] - meta_parameter["receiver_region"]:
        ] = j_A0 * meta_parameter["b_ac"]

        # mutual_inhibition AC on BC
        j_mutual_inhibition_ac_on_bc = np.zeros(meta_parameter["n_gridpoints"])
        j_mutual_inhibition_ac_on_bc[
            meta_parameter["n_gridpoints"] - meta_parameter["receiver_region"]:
        ] = j_A0 * meta_parameter["b_rp"]
        # mutual_inhibition BC on AC
        j_mutual_inhibition_bc_on_ac = np.zeros(meta_parameter["n_gridpoints"])
        j_mutual_inhibition_bc_on_ac[
            meta_parameter["n_gridpoints"] - meta_parameter["receiver_region"]:
        ] = j_A0 * meta_parameter["b_rp"]

        # mutual_inhibition AC on BC
        j_ac_rp = np.zeros(meta_parameter["n_gridpoints"])
        j_ac_rp[
            meta_parameter["n_gridpoints"] - meta_parameter["receiver_region"]:
        ] = j_A0 * meta_parameter["b_ac_rp"]
        # mutual_inhibition BC on AC
        j_ac_rp = np.zeros(meta_parameter["n_gridpoints"])
        j_ac_rp[
            meta_parameter["n_gridpoints"] - meta_parameter["receiver_region"]:
        ] = j_A0 * meta_parameter["b_ac_rp"]

        production_rate = (
            j_A, j_B, j_C,
            j_self_activation_ar_on_a, j_self_activation_br_on_b,
            j_self_activation_ac_on_ac, j_self_activation_bc_on_bc,
            j_mutual_inhibition_ac_on_bc, j_mutual_inhibition_bc_on_ac,
            j_ac_rp,
            j_R,
        )
        rxn_params = RXN_params_yuanqi(
            r_AR=koff_AR0, r_BR=koff_AR0, gamma=gamma0,
            n_ac=meta_parameter["n_ac"],
            k_ac=meta_parameter["k_ac"],
            n_rp=meta_parameter["n_rp"],
            k_rp=meta_parameter["k_rp"],
            error_rate=meta_parameter["error_rate"],
            sami_error_rate=meta_parameter["sami_error_rate"]
        )

        res = np.array(RD_solve(
            c_0_tuple, t, L=L, derivs_0=0, derivs_L=0,
            diff_coeff_fun=Diff_fun, diff_coeff_params=(diff_coeffs,),
            rxn_fun=RD_rxn, rxn_params=(rxn_params, production_rate),
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
