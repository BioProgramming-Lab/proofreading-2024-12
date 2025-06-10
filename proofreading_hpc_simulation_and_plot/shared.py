from enum import Enum
import numpy as np

# * initialize parameters

# Time points
t = np.arange(0, 3600 * 49, 3600)
t_record = t[[-1]]

# number of grid points
n_gridpoints = 451  # refer to the BMP model,because when their parameters are simulated, when L=3000, the value of citrine is close to 0
grid_spacing = 10  # um refer to 293T cell length
sender_region = 200

# Physical length of system
L = grid_spacing * (n_gridpoints - 1)  # um

# parameter infomation
num_thread = 800
chunk_size = 50
num_samples = num_thread * chunk_size
num_dimensions = 5


class Species(Enum):
    A = 0
    B = 1
    C = 2
    R = 3
    AC = 4
    BC = 5
    AR = 6
    BR = 7


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
        self.k_AC = 1e-4   # /nM/s
        self.r_AC = 1e-4  # /s
        # k (nM-1*s-1) and r (s-1) of B+A' <-> BA'
        self.k_BC = 1e-4  # /nM/s
        self.r_BC = 1e-3  # /s
        # k (nM-1*s-1) and r (s-1) of A+Receptor <-> A-Receptor complex
        self.k_AR = 4.5e-4  # /nM/s
        self.r_AR = 1e-3  # /s
        # k (nM-1*s-1) and r (s-1) of B+Receptor <-> B-Receptor complex
        self.k_BR = 4.5e-4  # /nM/s
        self.r_BR = 1e-3  # /s
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


DEFAULT_META_PARAMETERS = {
    # beta for activator, repressor co-affect one same output
    # k and n is extracted from the standalone version
    "b_ac_rp": 0,
    # beta for activator affect one individual output
    "b_ac": 0,
    "k_ac": 1,
    "n_ac": 1,
    # beta for repressor affect one individual output
    "b_rp": 0,
    "n_rp": 1,
    "k_rp": 1,
    "sender_region": 200,
    "receiver_region": 251,
    "sender_ratio": 1,
    "n_gridpoints": 451,
    "receptor_preequilibium": 0,
    "error_rate": 0,
    "sami_error_rate": 0,
}


class MetaParameter(dict):
    def __init__(self, **kwargs):
        super().__init__(DEFAULT_META_PARAMETERS)
        self.update(kwargs)

    def as_string(self):
        keys = []
        for k in self:
            if self[k] != DEFAULT_META_PARAMETERS[k]:
                keys.append(k)
        if len(keys) == 0:
            return "baseline"
        else:
            return "_".join(map(lambda x: x + str(self[x]), keys))


# setup hill function parameters
# (beta, k, n)
meta_parameters = [
    # base line (no feedback)
    # MetaParameter(),

    # sa error rate
    # MetaParameter(n_ac=1, k_ac=1, b_ac=1),
    # MetaParameter(n_ac=1, k_ac=1, b_ac=1, error_rate=1e-4),
    # MetaParameter(n_ac=1, k_ac=1, b_ac=1, error_rate=1e-3),
    # MetaParameter(n_ac=1, k_ac=1, b_ac=1, error_rate=1e-2),
    # MetaParameter(n_ac=1, k_ac=1, b_ac=1, error_rate=1e-1),
    # MetaParameter(n_ac=1, k_ac=1, b_ac=1, error_rate=0.5),
    # MetaParameter(n_ac=1, k_ac=1, b_ac=1, error_rate=1),

    # sami error rate
    # MetaParameter(b_ac_rp=1),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-4),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-3),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-2),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-1),
    # MetaParameter(b_ac_rp=1, sami_error_rate=0.5),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1),

    # sami error rate hill function rescue
    MetaParameter(b_ac_rp=1, sami_error_rate=0),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-4),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-3),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-2),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-1),
    MetaParameter(b_ac_rp=1, sami_error_rate=0, k_ac = 3),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-4, k_ac = 3),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-3, k_ac = 3),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-2, k_ac = 3),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-1, k_ac = 3),
    MetaParameter(b_ac_rp=1, sami_error_rate=0, k_ac = 6),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-4, k_ac = 6),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-3, k_ac = 6),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-2, k_ac = 6),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-1, k_ac = 6),
    MetaParameter(b_ac_rp=1, sami_error_rate=0, n_ac = 3),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-4, n_ac = 3),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-3, n_ac = 3),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-2, n_ac = 3),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-1, n_ac = 3),
    MetaParameter(b_ac_rp=1, sami_error_rate=0, n_ac=6),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-4, n_ac = 6),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-3, n_ac = 6),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-2, n_ac = 6),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-1, n_ac = 6),
    MetaParameter(b_ac_rp=1, sami_error_rate=0, k_rp = 3),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-4, k_rp = 3),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-3, k_rp = 3),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-2, k_rp = 3),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-1, k_rp = 3),
    MetaParameter(b_ac_rp=1, sami_error_rate=0, k_rp = 6),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-4, k_rp = 6),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-3, k_rp = 6),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-2, k_rp = 6),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-1, k_rp = 6),
    MetaParameter(b_ac_rp=1, sami_error_rate=0, n_rp = 3),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-4, n_rp = 3),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-3, n_rp = 3),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-2, n_rp = 3),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-1, n_rp = 3),
    MetaParameter(b_ac_rp=1, sami_error_rate=0, n_rp=6),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-4, n_rp = 6),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-3, n_rp = 6),
    # MetaParameter(b_ac_rp=1, sami_error_rate=1e-2, n_rp = 6),
    MetaParameter(b_ac_rp=1, sami_error_rate=1e-1, n_rp = 6),
]


# Replace with your MinIO credentials and endpoint
minio_endpoint = "http://10.16.21.7:9000"
access_key = "minioadmin"
secret_key = "minioadmin"
parameter_bucket = "parameter"  # Replace with the actual bucket name
result_bucket = "result1"
parameter_file = "parameter_20250216.zarr"
result_file = "result_20250216.zarr"
