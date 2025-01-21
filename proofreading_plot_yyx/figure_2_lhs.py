from contextlib import ExitStack
import numpy as np
from rd_solver import *
import os
import fcntl
import csv
from itertools import product


class RXN_params_yuanqi(object):
    """
    Container for reaction parameters
    """
    def __init__(self, **kwargs):
        """
        k_kin: kinase rate on membrane boundary(x=0)
        k_p: phosphotase rate in cytosol 
        k_CAp: association rate for Ap and A'
        r_CAp: dissociation rate for A'-Ap complex
        k_CBp: association rate for Ap and A'
        r_CBp: dissociation rate for A'-Bp complex
        """
        self.k_kin = 0.2 #nM/s
        self.k_p = 5 #/s
        self.k_CAp = 0.1 # /nM/s = 10^5 /M/s
        self.r_CAp = 0.1 # /s
        self.k_CBp = 0.1 # /nM/s = 10^5 /M/s
        self.r_CBp = 1 # /s

        self.ratio = self.r_CBp / self.r_CAp
        
        # Put in params that were specified in input
        for entry in kwargs:
            setattr(self, entry, kwargs[entry])

# * initialize parameters


# number of grid points 
n_gridpoints = 101
grid_spacing = 0.1 # µm
# time 60s 
t = np.arange(1, 61)

# Physical length of system
L = grid_spacing * (n_gridpoints - 1) # µm = 10µm


diff_coeffs = DIFFUSION()

# read pre-generated parameters
task_id = os.getenv("SLURM_ARRAY_TASK_ID")
if task_id is None:
    raise Exception(
        "Unable to find environment variable SLURM_ARRAY_TASK_ID"
    )
parameter_folder = "parameters_20250117_intracellular"
parameters = np.genfromtxt(
    "{}/{}.csv".format(parameter_folder, int(task_id) - 1), delimiter=","
)
output_folder = "result_20250117_intracellular"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    print("The new directory is created!")

   

for parameter in parameters.reshape((-1, 3)):
    # randomize parameters
    D_0 = parameter[0]
    K_kinase = parameter[1]
    K_phosphotase = parameter[2]

    # use a dictionary to save all filenames and result
    result_dict = {
        "parameters.csv": parameter,
    }

    # run with different beta (receiver region j_A0 factor)
    for c_AB, c_C in [(100, 50), (100, 100), (100, 200)]:
        c_0_tuple = (
            # c_A
            np.full(n_gridpoints, c_AB),
            # c_B
            np.full(n_gridpoints, c_AB),
            # c_C
            np.full(n_gridpoints, c_C),
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
            k_kin=K_kinase,
            k_p=K_phosphotase,
        )

        result_dict["result_C0{}.csv".format(C_0)] = np.array(RD_solve(
            c_0_tuple, t, L=L, derivs_0=0, derivs_L=0,
            diff_coeff_fun=Diff_fun, diff_coeff_params=(diff_coeffs,),
            rxn_fun=RD_rxn, rxn_params=(rxn_params,),
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

