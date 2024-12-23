# rust_lib_proofreading

## Requirement
- **pyenv-virtualenv**, if you use conda/mamba, you **DONOT** need to manually install this.
- **maturin**, to build PyO3 project
- **rustup**, **cargo**, for rust compiling

## Installation
please follow the official [PyO3 guide](https://pyo3.rs/v0.23.3/index.html) if customization is needed.
```shell
conda activate your_env
maturin develop
```

## Usage
```python
import rust_lib_proofreading

parameters = np.array(rust_lib_proofreading.read_parameter("path/to/parameters.csv"))
n_species = 8
sol = np.array(rust_lib_proofreading.read_solve_at_end("path/to/sol.csv"), n_species)

print(parameters.shape)
print(sol.shape)
```