from shared import *
import time
import numpy as np
import zarr

# open result zarr
store = "result_20250408_error_rate_rescue"
z = zarr.open(
    store=store,
    mode="w",
)

start = time.time()
for j, meta_parameter in enumerate(meta_parameters):
    group_meta_parameter = z.require_group(meta_parameter.as_string())
    group_meta_parameter.attrs.update(meta_parameter)
    for _t in t_record:
        array_s_sim_l = group_meta_parameter.require_array(
            name=str(_t),
            dtype=np.float64,
            shape=(len(Species), chunk_size * num_thread, n_gridpoints),
            chunks=(len(Species), chunk_size, n_gridpoints),
        )

print(time.time() - start)
