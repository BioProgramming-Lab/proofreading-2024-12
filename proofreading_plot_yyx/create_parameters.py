from scipy.stats import qmc
import csv
import numpy as np
import zarr
from shared import *

# create sampler and sample parameters
sampler = qmc.LatinHypercube(d=num_dimensions)
lhs_samples = sampler.random(n=num_samples)

# scale to the range of each parameter
# D_0 = np.random.uniform(0, 50)
# j_A0 = 10 ** np.random.uniform(-5, -3)
# j_R0 = 10 ** np.random.uniform(-5, -3)
# koff_AR0 = 10 ** np.random.uniform(-5, -3)
# gamma0 = 2 * 10 ** np.random.uniform(-20, -3)

lhs_samples[:, 0] = lhs_samples[:, 0] * 50

log_min = -5
log_max = -2
lhs_samples[:, 1:4] = 10 ** (
    log_min + (log_max - log_min) * lhs_samples[:, 1:4]
)

log_min = np.log10(2.8e-4)
log_max = np.log10(2.8e-3)
lhs_samples[:, 4] = 10 ** (
    log_min + (log_max - log_min) * lhs_samples[:, 4]
)

z = zarr.create_array(
    store="parameters_20250212",
    data=lhs_samples,
    chunks=(chunk_size, num_dimensions),
)
