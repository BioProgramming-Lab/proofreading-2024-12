from enum import Enum
import numpy as np

# * initialize parameters

# Time points
t = np.arange(0, 61)
t_record = t[[-1]]

# number of grid points
n_gridpoints = 101  # refer to the BMP model,because when their parameters are simulated, when L=3000, the value of citrine is close to 0
grid_spacing = 0.1  # um refer to 293T cell length

# Physical length of system
L = grid_spacing * (n_gridpoints - 1)  # um

# parameter infomation
num_thread = 800
chunk_size = 50
num_samples = num_thread * chunk_size
num_dimensions = 4


class Species(Enum):
    A = 0
    B = 1
    C = 2
    Ap = 3
    Bp = 4
    CAp = 5
    CBp = 6


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
        self.k_kin = 0.2  # nM/s
        self.k_p = 5  # /s
        self.k_CAp = 0.1  # /nM/s = 10^8 /M/s
        self.r_CAp = 0.1  # /s
        self.k_CBp = 0.1  # /nM/s = 10^8 /M/s
        self.r_CBp = 1  # /s

        self.ratio = self.r_CBp / self.r_CAp

        # Put in params that were specified in input
        for entry in kwargs:
            setattr(self, entry, kwargs[entry])


DEFAULT_META_PARAMETERS = {}


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
    MetaParameter(),
]


# Replace with your MinIO credentials and endpoint
minio_endpoint = "http://10.16.21.7:9000"
access_key = "minioadmin"
secret_key = "minioadmin"
parameter_bucket = "parameter"  # Replace with the actual bucket name
result_bucket = "result1"
parameter_file = "parameter_20250216.zarr"
result_file = "result_20250216.zarr"
