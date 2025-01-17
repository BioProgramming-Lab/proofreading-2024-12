import numpy as np
import rust_lib_proofreading
from RXN_params_yuanqi import RXN_params_yuanqi

# this decorator look up the subset indices recorded in Simulations.subsets
# and combine those indices with np.bitwise_and
# ! have to specificly use subses=[] and beta=[].


def subset_decorator(func):
    def wrapper(self, *args, subsets=None, not_subsets=None, **kwargs):
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
    def __init__(self, path="./", betas=[""], n_species=7):
        # betas need to be a list

        # initial dicts
        self.subsets = {}
        self.parameters = {}
        self.solve = {}

        # load parameters
        self.parameters[None] = np.array(
            rust_lib_proofreading.read_parameters(path + "/parameters.csv")
        )
        self.n_sims = self.parameters[None].shape[0]
        print("Parameters loaded.")
        print(self.get_parameters_shape())

        # load result
        for b in betas:
            self.solve[b] = np.array(rust_lib_proofreading.read_solve_at_end(
                path + "/result_{}.csv".format(b).replace("_.csv", ".csv"), n_species
            ))

        self.check_validity()
        print("Simulations result loaded.")
        print(self.get_solve_shape())

    def get_keys(self):
        return self.solve.keys()

    def check_validity(self):
        # should call this as soon as data loaded
        for k in self.solve:
            self.add_subset(
                ("valid", k),
                np.all(
                    np.bitwise_and(
                        np.isfinite(self.solve[k]),
                        self.solve[k] > 0
                    ),
                    axis=-1
                )
            )

    def get_parameters_shape(self):
        shapes = {}
        for k in self.parameters:
            shapes[("parameters", k)] = self.get_parameters().shape
        return shapes

    def get_solve_shape(self):
        # should manually check data size after data loaded
        shapes = {}
        for k in self.solve:
            shapes[k] = self.solve[k].shape
        return shapes

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
        return np.log10(
            self.get_solve(beta)
        )

    @subset_decorator
    def get_AC(self, beta=0, subsets=None, not_subsets=None):
        return self.get_solve(beta)[:, 5]

    @subset_decorator
    def get_AC_log(self, beta=0, subsets=None, not_subsets=None):
        return np.log10(self.get_AC(beta))

    @subset_decorator
    def get_fidelity(self, beta=0, subsets=None, not_subsets=None):
        # enhanced fidelity
        return self.get_AC(beta) / \
            self.get_solve(beta)[:, 6] / \
            RXN_params_yuanqi().ratio

    @subset_decorator
    def get_fidelity_log(self, beta=0, subsets=None, not_subsets=None):
        # enhanced fidelity
        return np.log10(self.get_fidelity(beta))

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
