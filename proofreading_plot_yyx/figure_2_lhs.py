from contextlib import ExitStack
import numpy as np
from rd_solver import *
import os
import fcntl
import csv


class RXN_params_yuanqi(object):
    """
    Container for reaction parameters
    """

    def __init__(self, **kwargs):
        """
        k_AC: association rate for A and Receptor
        r_AC: dissociation rate for AA'
        k_AR: association rate for A and Receptor
        r_AR: dissociation rate for A-Receptor complex
        k_BC: association rate for B and Receptor
        r_BC: dissociation rate for BA'
        k_BR: association rate for B and Receptor
        r_BR: dissociation rate for B-Receptor complex
        RTotal: Total receptor number
        gamma: receptor turnover rate ([X-R] --> R)
        deg: universal degradation rate for free proteins
        """
        # RTotal: Total receptor number (nM)
        # self.RTotal = 2.7 # nM  0.0375 – 2.7 nM (corresponding to 18 –1300 molecules/um2 and a 800 um Matrigel layer on the cell surface)
        # k (nM-1*s-1) and r (s-1) of A+A' <-> AA'
        self.k_AC = 1e-4   # /nM/min
        self.r_AC = 1e-4  # /min
        # k (nM-1*s-1) and r (s-1) of B+A' <-> BA'
        self.k_BC = 1e-4  # nM/smin
        self.r_BC = 1e-3  # /min
        # k (nM-1*s-1) and r (s-1) of A+Receptor <-> A-Receptor complex
        self.k_AR = 4.5e-4  # nM/min
        self.r_AR = 1e-3  # /s
        # k (nM-1*s-1) and r (s-1) of B+Receptor <-> B-Receptor complex
        self.k_BR = 4.5e-4  # nM/min
        self.r_BR = 1e-3  # /min
        # Recycling rate (s-1) of A-Receptor and B-Receptor
        self.gamma = 4e-4

        self.ratio = self.r_BC / self.r_AC

        # Not considering basal degradation of species
        self.deg = 2e-5

        # Hill function
        self.k = 1
        self.n = 3

        # Put in params that were specified in input
        for entry in kwargs:
            setattr(self, entry, kwargs[entry])

# * initialize parameters


# Time points
t = np.arange(3600 * 1.0, 3600 * 48, 3600)

# number of grid points
n_gridpoints = 451  # refer to the BMP model,because when their parameters are simulated, when L=3000, the value of citrine is close to 0
grid_spacing = 10  # um refer to 293T cell length
sender_region = 200

# Physical length of system
L = grid_spacing * (n_gridpoints - 1)  # um


diff_coeffs = DIFFUSION()
c_0_tuple = (
    np.zeros(n_gridpoints), np.zeros(n_gridpoints),
    np.zeros(n_gridpoints), np.zeros(n_gridpoints),
    np.zeros(n_gridpoints), np.zeros(n_gridpoints),
    np.zeros(n_gridpoints), np.zeros(n_gridpoints),
)

# read pre-generated parameters
task_id = os.getenv("SLURM_ARRAY_TASK_ID")
if task_id is None:
    raise Exception(
        "Unable to find environment variable SLURM_ARRAY_TASK_ID"
    )
parameter_folder = "parameters_20241223"
parameters = np.genfromtxt(
    "{}/{}.csv".format(parameter_folder, int(task_id) - 1), delimiter=","
)
output_folder = "result_20241223"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    print("The new directory is created!")

# setup betas
betas = [0.5, 1, 2]

for parameter in parameters:
    # randomize parameters
    D_0 = parameter[0]
    j_A0 = parameter[1]
    j_R0 = parameter[2]
    koff_AR0 = parameter[3]
    gamma0 = parameter[4]

    # use a dictionary to save all filenames and result
    result_dict = {
        "parameters.csv": parameter,
    }

    # run with different beta (receiver region j_A0 factor)
    for beta in betas:
        # * run simulation
        with_positive_feedbacks = [True, False]
        for with_positive_feedback in with_positive_feedbacks:
            diff_coeffs.D_A = D_0
            diff_coeffs.D_B = D_0
            diff_coeffs.D_C = D_0
            diff_coeffs.D_complex = D_0

            # Initialize
            j_A = np.zeros(n_gridpoints)
            j_A[0:sender_region] = j_A0
            j_B = np.zeros(n_gridpoints)
            j_B[0:sender_region] = j_A0
            j_C = np.zeros(n_gridpoints)
            j_C[0:sender_region] = j_A0 * 2
            j_R = np.zeros(n_gridpoints)
            j_R = j_R0

            j_a = np.zeros(n_gridpoints)
            j_a[sender_region:] = j_A0 * 0
            j_b = np.zeros(n_gridpoints)
            j_b[sender_region:] = j_A0 * 0

            if with_positive_feedback:
                j_ac = np.zeros(n_gridpoints)

                j_ac[sender_region:] = j_A0 * beta
                j_bc = np.zeros(n_gridpoints)
                j_bc[sender_region:] = j_A0 * beta

                production_rate = (j_A, j_B, j_C, j_a, j_b, j_ac, j_bc, j_R)
                rxn_params = RXN_params_yuanqi(
                    r_AR=koff_AR0, r_BR=koff_AR0, gamma=gamma0
                )

                result_dict["with_feedback_{}.csv".format(beta)] = np.array(RD_solve(
                    c_0_tuple, t, L=L, derivs_0=0, derivs_L=0,
                    diff_coeff_fun=Diff_fun, diff_coeff_params=(diff_coeffs,),
                    rxn_fun=RD_rxn, rxn_params=(rxn_params, production_rate),
                    rtol=1.49012e-8, atol=1.49012e-8
                ))[:, -1, :]

            else:
                j_ac = np.zeros(n_gridpoints)
                j_ac[sender_region:] = j_A0 * 0
                j_bc = np.zeros(n_gridpoints)
                j_bc[sender_region:] = j_A0 * 0

                production_rate = (j_A, j_B, j_C, j_a, j_b, j_ac, j_bc, j_R)
                rxn_params = RXN_params_yuanqi(
                    r_AR=koff_AR0, r_BR=koff_AR0, gamma=gamma0
                )

                result_dict["without_feedback_{}.csv".format(beta)] = np.array(RD_solve(
                    c_0_tuple, t, L=L, derivs_0=0, derivs_L=0,
                    diff_coeff_fun=Diff_fun, diff_coeff_params=(diff_coeffs,),
                    rxn_fun=RD_rxn, rxn_params=(rxn_params, production_rate),
                    rtol=1.49012e-8, atol=1.49012e-8
                ))[:, -1, :]

    # use fcntl save all data at same time
    with ExitStack() as stack:
        files = [
            stack.enter_context(
                open(os.path.join(output_folder, fname), "a")
            ) for fname in result_dict
        ]

        for f in files:
            fcntl.flock(f, fcntl.LOCK_EX)

        for f, k in zip(files, result_dict):
            csv_writer = csv.writer(f)
            if k == "parameters.csv":
                csv_writer.writerow(result_dict[k])
            else:
                csv_writer.writerows(result_dict[k])

        for f in files:
            fcntl.flock(f, fcntl.LOCK_UN)
