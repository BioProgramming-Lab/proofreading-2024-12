from scipy.stats import qmc
import csv
import numpy as np

# define the number of samples and dimensions
num_samples = 40000
num_dimensions = 5

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
log_max = -3
lhs_samples[:, 1:4] = 10 ** (
    log_min + (log_max - log_min) * lhs_samples[:, 1:4]
)
log_min = -5
log_max = -2.5
lhs_samples[:, 4] = 10 ** (
    log_min + (log_max - log_min) * lhs_samples[:, 4]
)

# parameter_filename = "parameters.csv"
# with open(parameter_filename, "a") as f_p:
#     csv_writer = csv.writer(f_p)
#     csv_writer.writerows(lhs_samples)

# save individual chunks for processing
num_thread = 800
parameter_folder = "parameters_20241223"
for i, chunk in enumerate(np.vsplit(lhs_samples, num_thread)):
    with open("{}/{}.csv".format(parameter_folder, i), "w") as f_p:
        csv_writer = csv.writer(f_p)
        csv_writer.writerows(chunk)
