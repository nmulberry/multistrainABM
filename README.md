# multistrainABM

An agent-based model for multistrain dynamics inspired by the pneumococcus. The model allows for an arbitrary number of circulating strains.


## Compiling the model
First, clone [ODE-event-solvers](https://github.com/nmulberry/ode-event-solvers) and update the path to this package in `Cargo.toml`.
Now you should be able to compile using cargo:
```
    cargo build --release
```
This should build the binary 
```
    target/release/multiabm_samplestrains
```
## Running analyses
Model analyses are performed in `R`. For a simple example, run `R/sample_run.R`.
The analyses from the paper can be reproduced by running `R/main.R`.
