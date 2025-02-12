from contextlib import ExitStack
import numpy as np
from rd_solver import *
import os
import fcntl
import csv
from itertools import product
import zarr
from enum import Enum

from shared import *



diff_coeffs = DIFFUSION()

# read pre-generated parameters
task_id = int(os.getenv("SLURM_ARRAY_TASK_ID"))
if task_id is None:
    raise Exception(
        "Unable to find environment variable SLURM_ARRAY_TASK_ID"
    )

parameter_dir = "parameters_20250212"
z = zarr.open(
    store=parameter_dir,
    mode="r",
)
parameters = z[
    ((task_id - 1) * chunk_size):(task_id * chunk_size)
]

output_dir = "result_20250212"
z = zarr.open_group(
    store=output_dir,
    mode="a",
)




# setup hill function parameters
# (beta, k, n)
meta_parameters = []

# ------------ 20250115
# meta_parameters += [
#     # just to take a note on the previously default parameter
#     # DEFAULT_META_PARAMETERS = {
#     #     "b": 0,
#     #     "k": 1,
#     #     "n": 1,
#     #     "sender_region": 200,
#     #     "receiver_region": 251,
#     #     "sender_ratio": 1,
#     #     "n_gridpoints": 451,
#     #     "proofreading_basal": 0,
#     #     "mutual_inhibition": 0,
#     #     "receptor_preequilibium": 0
#     # }


#     # base line (no feedback)
#     MetaParameter(),

#     # hill function parameters
#     MetaParameter(n=1, k=1, b=1),
#     MetaParameter(n=3, k=1, b=1),
#     MetaParameter(n=1, k=3, b=1),
#     MetaParameter(n=1, k=6, b=1),
#     MetaParameter(n=1, k=1, b=0.5),
#     MetaParameter(n=1, k=1, b=2),

#     # increase sender secretion rate
#     MetaParameter(sender_ratio=3),

#     # extend sender region to the right
#     MetaParameter(sender_region=451),

#     # extend sender region X3 to the left
#     MetaParameter(sender_region=600, n_gridpoints=851),

#     # [0, 0, 0],  # Case 1: no feedback
#     # same as baseline
#     # [1, 0, 0], # Case 2:  basal secretion in proofreading region
#     MetaParameter(proofreading_basal=1),
#     # [0, 1, 0], # Case 3:  self-activation
#     # same as n=1, k=1, b=1
#     # [1, 0, 1], # Case 4:  mutual inhibition
#     MetaParameter(proofreading_basal=1, mutual_inhibition=1),
#     # [1, 1, 1] # Case 5:  self-activation + mutual inhibition
#     MetaParameter(n=1, k=1, b=1, proofreading_basal=1, mutual_inhibition=1),
#     # pre equilibrium receptor concentration
#     MetaParameter(receptor_preequilibium=1),
# ]


# --------------- 20250116
# just to take a note on the previously default parameter
# DEFAULT_META_PARAMETERS = {
#     "b_ac_rp": 0,
#     "b_ac": 0,
#     "k_ac": 1,
#     "n_ac": 1,
#     "b_rp": 0,
#     "n_rp": 1,
#     "k_rp": 1,
#     "sender_region": 200,
#     "receiver_region": 251,
#     "sender_ratio": 1,
#     "n_gridpoints": 451,
#     "receptor_preequilibium": 0
# }

# meta_parameters.append(MetaParameter())
# for b_ac, b_rp, k_ac, k_rp, n_ac, n_rp in product(
#     [0, 1, 2],
#     [0, 1, 2],
#     [1, 10, 100],
#     [1, 10, 100, 1000],
#     [1, 3],
#     [1, 3],
# ):
#     if b_ac == 0 and b_rp == 0:
#         continue
#     meta_parameters.append(
#         MetaParameter(
#             b_ac=b_ac, b_rp=b_rp, k_ac=k_ac, k_rp=k_rp, n_ac=n_ac, n_rp=n_rp,
#         )
#     )
# for b_ac_rp, k_ac, k_rp, n_ac, n_rp in product(
#     [1, 2],
#     [1, 10, 100],
#     [1, 10, 100, 1000],
#     [1, 3],
#     [1, 3],
# ):
#     meta_parameters.append(
#         MetaParameter(
#             b_ac_rp=b_ac_rp, k_ac=k_ac, k_rp=k_rp, n_ac=n_ac, n_rp=n_rp,
#         )
#     )

# ------------ 20250121
# meta_parameters += [
#     # just to take a note on the previously default parameter
#     # DEFAULT_META_PARAMETERS = {
#     #     # beta for activator, repressor co-affect one same output
#     #     # k and n is extracted from the standalone version
#     #     "b_ac_rp": 0,
#     #     # beta for activator affect one individual output
#     #     "b_ac": 0,
#     #     "k_ac": 1,
#     #     "n_ac": 1,
#     #     # beta for repressor affect one individual output
#     #     "b_rp": 0,
#     #     "n_rp": 1,
#     #     "k_rp": 1,
#     #     "sender_region": 200,
#     #     "receiver_region": 251,
#     #     "sender_ratio": 1,
#     #     "n_gridpoints": 451,
#     #     "receptor_preequilibium": 0
#     # }

#     # MetaParameter(
#     #     b_ac_rp=1,
#     # ),
#     # MetaParameter(
#     #     b_ac_rp=2,
#     # ),
#     # MetaParameter(
#     #     b_ac_rp=2, n_rp=3
#     # ),
#     # MetaParameter(
#     #     sender_ratio=3
#     # )
# ]

meta_parameters += [
    # base line (no feedback)
    MetaParameter(),

    # hill function parameters
    MetaParameter(n_ac=1, k_ac=1, b_ac=1),
    MetaParameter(n_ac=3, k_ac=1, b_ac=1),
    MetaParameter(n_ac=1, k_ac=3, b_ac=1),
    MetaParameter(n_ac=1, k_ac=6, b_ac=1),
    MetaParameter(n_ac=1, k_ac=1, b_ac=0.5),
    MetaParameter(n_ac=1, k_ac=1, b_ac=2),

    # increase sender secretion rate
    MetaParameter(sender_ratio=3),

    # extend sender region to the right
    MetaParameter(sender_region=451),

    MetaParameter(b_rp=1),
    MetaParameter(b_ac_rp=1),
]



for i, parameter in enumerate(parameters.reshape((-1, 5))):
    # randomize parameters
    D_0 = parameter[0]
    j_A0 = parameter[1]
    j_R0 = parameter[2]
    koff_AR0 = parameter[3]
    gamma0 = parameter[4]

    # run with different beta (receiver region j_A0 factor)
    for meta_parameter in meta_parameters:
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
        # j_R: the secretion rate of free receptor
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
        )

        res = np.array(RD_solve(
            c_0_tuple, t, L=L, derivs_0=0, derivs_L=0,
            diff_coeff_fun=Diff_fun, diff_coeff_params=(diff_coeffs,),
            rxn_fun=RD_rxn, rxn_params=(rxn_params, production_rate),
            rtol=1.49012e-8, atol=1.49012e-8
        ))

        # storage configuration
        # root
        # ├─ meta_parameter
        # |  ├─ time
        # |  |  ├─ chemical species
        # |  |  |  n_sim X n_grid (chunk: chunk_size X n_grid)
        # |  |  |
        # ...

        group_meta_parameter = z.require_group(meta_parameter.as_string())
        group_meta_parameter.attrs.update(meta_parameter)
        for _t in t_record:
            group_t = group_meta_parameter.require_group(str(_t))
            for s in Species:
                array_s = group_t.require_array(
                    name=s.name,
                    dtype=res.dtype,
                    shape=(chunk_size * num_thread, n_gridpoints),
                    chunks=(chunk_size, n_gridpoints),
                )

                array_s[
                    ((task_id - 1) * chunk_size) + i, :
                ] = res[s.value, list(t).index(_t), :]
