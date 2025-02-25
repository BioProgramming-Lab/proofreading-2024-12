import numpy as np
from shared import *
import zarr

# this decorator look up the subset indices recorded in Simulations.subsets
# and combine those indices with np.bitwise_and
# ! have to specificly use subsets=[] and beta=[].


def subset_decorator(func):
    def wrapper(
        self,
        *args,
        subsets=None,
        not_subsets=None,
        **kwargs
    ):
        # subset is not pass to inner function, subsets only works for the
        # wrapper
        output = func(self, *args, **kwargs)
        indices = np.ones(self.n_sims, dtype=bool)

        if subsets is not None:
            for subset in subsets:
                if subset in self.subsets:
                    indices = np.bitwise_and(indices, self.subsets[subset])
                    continue
                # if cannot find subset in dict, then try to construct one from
                # beta
                subset = (
                    subset,
                    kwargs.get("beta"),
                )
                if subset in self.subsets:
                    indices = np.bitwise_and(indices, self.subsets[subset])
                    continue
                # still cannot find the subset
                print("no subset named {}, skipped".format(subset))

        if not_subsets is not None:
            for subset in not_subsets:
                if subset in self.subsets:
                    indices = np.bitwise_and(
                        indices, np.bitwise_not(self.subsets[subset])
                    )
                    continue
                # if cannot find subset in dict, then try to construct one from
                # beta
                subset = (
                    subset,
                    kwargs.get("beta"),
                )
                if subset in self.subsets:
                    indices = np.bitwise_and(
                        indices, np.bitwise_not(self.subsets[subset])
                    )
                    continue
                # still cannot find the subset
                print("no subset named {}, skipped".format(subset))

        if isinstance(output, np.ndarray):
            return output[indices]
        elif isinstance(output, tuple):
            for i in range(len(output)):
                output[i] = output[i][indices]
        else:
            raise RuntimeError("Unreachable")
        return output

    return wrapper


class Simulations:
    def __init__(
        self,
        result_path="./",
        parameter_path="",
        betas=None,
        t=None,
        l=-1,
        species=Species
    ):
        # betas need to be a list

        # initial dicts
        self.subsets = {}
        self.parameters = {}
        self.solve = {}
        self.meta_parameters = {}

        self.t = t
        self.species = species

        # load parameters
        z = zarr.open(parameter_path, "r")
        self.parameters[None] = z[:]
        self.n_sims = self.parameters[None].shape[0]
        # self.genrate_koff_gamma0()
        print("Parameters loaded.")
        print(self.get_parameters_shape())

        # load result.
        z = zarr.open(result_path, "r")
        for b in betas if betas is not None else z.keys():
            print("reading {}".format(b))
            group_b = z[b]
            self.meta_parameters[b] = group_b.attrs
            if self.t is None:
                self.t = 0
                t_str = "0"
                for _t in group_b.keys():
                    if float(_t) > self.t:
                        self.t = float(_t)
                        t_str = _t
            self.solve[b] = group_b[t_str][:, :, l].swapaxes(0, 1)

        self.check_validity()
        print("Simulations result loaded.")
        print(self.get_solve_shape())

    def get_keys(self):
        return self.solve.keys()

    def check_validity(self):
        # should call this as soon as data loaded
        for k in self.solve:
            subset = np.ones(self.n_sims, dtype=np.bool_)
            for s in self.species:
                subset = np.bitwise_and(
                    subset,
                    np.bitwise_and(
                        np.isfinite(self.solve[k][:, s.value]),
                        self.solve[k][:, s.value] > 0
                    )
                )

            self.add_subset(("valid", k), subset)

    # no longer in use
    def genrate_koff_gamma0(self):
        self.parameters[None] = np.append(
            self.get_parameters(),
            self.get_parameters()[:, [3]] +
            self.get_parameters()[:, [4]],
            axis=-1
        )

    def get_parameters_shape(self):
        shapes = {}
        for k in self.parameters:
            shapes[("parameters", k)] = self.get_parameters().shape
        return shapes

    def get_solve_shape(self):
        return {
            k: v.shape
            for k, v in self.solve.items()
        }

    @subset_decorator
    def get_parameters(self, subsets=None, not_subsets=None):
        return self.parameters[None]

    @subset_decorator
    def get_parameters_log(self, subsets=None, not_subsets=None):
        parameters = self.get_parameters()
        parameters[:, 1:] = np.log10(parameters[:, 1:])
        return parameters

    @subset_decorator
    def get_solve(self, beta=0, subsets=None, not_subsets=None):
        return self.solve[beta]

    @subset_decorator
    def get_solve_log(self, beta=0, subsets=None, not_subsets=None):
        return np.log10(self.get_solve(beta))

    @subset_decorator
    def get_AC(self, beta=0, subsets=None, not_subsets=None):
        return self.get_solve(beta)[:, self.species.CAp.value]

    @subset_decorator
    def get_AC_log(self, beta=0, subsets=None, not_subsets=None):
        return np.log10(self.get_AC(beta))

    @subset_decorator
    def get_fidelity(self, beta=0, subsets=None, not_subsets=None):
        # enhanced fidelity
        return self.get_AC(beta) / RXN_params_yuanqi.ratio / \
            self.get_solve(beta)[:, self.species.CBp.value]

    @subset_decorator
    def get_fidelity_log(self, beta=0, subsets=None, not_subsets=None):
        # enhanced fidelity
        return np.log10(self.get_fidelity(beta))
    
    @subset_decorator
    def get_r_total(self, beta=0, subsets=None, not_subsets=None):
        return self.get_solve(beta)[:, [self.species.AR.value, self.species.BR.value, self.species.R.value]].sum(axis=1)
    
    @subset_decorator
    def get_r_total_log(self, beta=0, subsets=None, not_subsets=None):
        return np.log10(self.get_r_total(beta))

    @subset_decorator
    def get_parameters_fidelity_concentration(
        self, beta=0, subsets=None, not_subsets=None
    ):
        return (
            self.get_parameters(),
            self.get_fidelity(beta),
            self.get_AC(beta)
        )

    @subset_decorator
    def get_parameters_fidelity_concentration_log(
        self, beta=0, subsets=None, not_subsets=None
    ):
        return (
            self.get_parameters_log(),
            self.get_fidelity_log(beta),
            self.get_AC_log(beta)
        )

    def add_subset(self, name, subset):
        if len(subset) != self.n_sims:
            print(
                "subset length not equals to simulation counts, {} and {}"
                .format(len(subset), self.n_sims)
            )
            return
        self.subsets[name] = subset

    def get_subset(self, subsets=None, not_subsets=None):
        indices = np.ones(self.n_sims, dtype=bool)

        if subsets is not None:
            for subset in subsets:
                if subset in self.subsets:
                    indices = np.bitwise_and(indices, self.subsets[subset])
                    continue
                # still cannot find the subset
                print("no subset named {}, skipped".format(subset))

        if not_subsets is not None:
            for subset in not_subsets:
                if subset in self.subsets:
                    indices = np.bitwise_and(
                        indices, np.bitwise_not(self.subsets[subset])
                    )
                    continue
                # still cannot find the subset
                print("no subset named {}, skipped".format(subset))
        return indices
