from scipy.stats import qmc
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

log_min = np.log10(0.0375)
log_max = np.log10(2.7)
lhs_samples[:, 2] = 10 ** (
    log_min + (log_max - log_min) * lhs_samples[:, 2]
)

log_min = -5
log_max = -2
lhs_samples[:, 4] = 10 ** (
    log_min + (log_max - log_min) * lhs_samples[:, 4]
)

store = "parameters_20250217"
z = zarr.open(
    store=store,
    mode="w",
    shape=lhs_samples.shape,
    chunks=(chunk_size, num_dimensions),
    dtype=lhs_samples.dtype,
)

z[:] = lhs_samples