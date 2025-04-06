# Molecular Dynamics Simulation of a Van der Waals Gas

Simple molecular dynamics program that deals with a Van der Waals gas of
particles. This program was developed as a part of the final assignment of the
Advanced Computer Tools subject.

The project includes both a serial version of the code and a version
parallelized using OpenMPI, divided in two folders (`serial/` and `parallel/`).

## Compiling and running the simulation

Each version of the project includes a Makefile which should be used in order to
compile its source code. In order to download and build this program one should
follow the next steps:

```
$ git clone https://github.com/Eines-Informatiques-Avancades/Project_2025.git
$ cd Project_2025
```

Once inside the root directory of the project one must move into one of the 2
folder containing either the parallelized version of the code or the serial one
and run:

```
$ make
```

Compiling this project will produce several independent binaries
which can be used to perform different tasks.

- `vdw_gas`: the main simulation program. Will perform the simulation and
    output the results into different text files. Results are outputted in
    reduced units.

- `binning`: performs data sampling of the results of the Van der Waals gas
    molecular dynamics simulation from a given text file outputted by `vdw_gas`.

Edit the `input_parameters.in` file before running the simulation to set the
system parameters. There's no need to recompile the program once this file is
edited, as it is read at runtime.

The user might want to execute the different parts of the serial program
sequentally. In order to ease and automatize this process, the `run.sh` shell
script is included. This will log the program's execution and move the output
files into a separate folder.

The Makefile inside the `parallel/` folder also includes rules to compile and
run the code in the cerqt2 computing cluster. In case one wishes to run the code
in said cluster, it should be compiled in it too, so that the produced binary is
optimized for its hardware.

`make cluster-compile` will send a job to compile the parallel version of the
code in cerqt2.

`make cluster-run-` followed by either 1, 2, 4, 8, 16, 32 or 40 will execute the
simulation using the number of cores corresponding to that number. `make
cluster-run-all` will send one job for each of those numbers so that the user
can compare the performance of the code with different numbers of processors.

## Plotting the results

After running the simulation, output files with the results will be produced. In
order to better review and analyse these results, a set of Python scripts are
included under the `plot/` folder. These scripts make use of several Python
libraries to produce different plots of the results (see
[Dependencies](#dependencies)).

The scripts can be run manually (one by one) or using the `plot.sh` shell
script, which will execute all of them sequentially.

## Dependencies

- A Fortran compiler (E.g.: [GNU Fortran](https://gcc.gnu.org/fortran/))
- `make`
- [OpenMPI](https://www.open-mpi.org/) (parallel code only)

In order to visualize the results of the simulation, several Python plotting
scripts are provided. These scripts have the following dependencies:

- [Python 3](https://www.python.org/)
  - [NumPy](https://numpy.org/)
  - [Matplotlib](https://matplotlib.org/)
  - [imageio](https://imageio.readthedocs.io/en/stable/)
  - [pandas](https://pandas.pydata.org/)
- [cm-super TeX font package](https://ctan.org/pkg/cm-super)

## Project structure

The main simulation program consists of a main file (`vdw_gas.f90`) and several
modules stored under the `include/` directory, which contain the necessary
subroutines needed for the simulation separated by topic or task.

As previously mentioned, system parameters are set under `input_parameters.in`,
which is read at runtime.

Another standalone program, `binning.f90`, is also included to perform the
statistical analysis of the results produced by the simulation.

The project tree for each of the two versions has the following structure:

```
.
|-- include/
|   |-- lj_forces.f90       <- Lennard-Jones forces computation.
|   |-- geometry.f90        <- PBC
|   |-- initial_conf.f90    <- Initial configuration generation.
|   |-- integrators.f90     <- Time-step integrators.
|   |-- io.f90              <- Read/write from/to external files.
|   |-- thermodynamics.f90  <- Computation of different measures.
|   `-- thermostat.f90      <- Implementation of the Andersen thermostat.
|-- plot/                   <- Python plotting scripts for the results.
|   `-- ...
|-- test/                   <- Files related to the program testing.
|   `-- ...
|-- Makefile
|-- binning.f90             <- Binning statistical analysis of the results.
|-- input_parameters.in     <- System and simulation parameter definition.
|-- run.sh                  <- Wrapper script for executing the full program.
`-- vdw_gas.f90             <- Main simulation program.
```

## Credits and contributors

- Ricard Rodríguez: system initialization, periodic boundary condition, project
    coordination
- Oriol Miró: integrators
- Alejandro Díaz: forces
- Joan Serrano: statistical analysis
- Huang Haoyu: post-trajectory analysis
- Itziar Rabal: final testing

The Matplotlib Python plotting scripts under the `plot/` subfolder make
use of the `science.mplstyle` Matplotlib style, which is part of the
[SciencePlots project](https://github.com/garrettj403/SciencePlots) by [John
Garrett](https://github.com/garrettj403).
