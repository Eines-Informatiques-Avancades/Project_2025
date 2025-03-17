# Molecular Dynamics Simulation of a Van der Waals Gas

Simple molecular dynamics program that deals with a Van der Waals gas of
particles. This program was developed as a part of the final assignment of the
Advanced Computer Tools subject.

## Compiling and running the simulation

This project includes a Makefile which should be used in order to compile its
source code. In order to download and build this program one should follow the
next steps:

```
$ git clone https://github.com/Eines-Informatiques-Avancades/Project_2025.git
$ cd Project_2025
$ make
```

Compiling this project will produce several independent binaries
which can be used to perform different tasks.

- `vdw_gas`: the main simulation program. Will perform the simulation and
    output the results into different text files.

- `binning`: performs data sampling of the results of the Van der Waals gas
    molecular dynamics simulation from a given text file outputted by `vdw_gas`.

- `jackknife`: performs data resampling of the results of the Van der Waals gas
    molecular dynamics simulation from a given text file outputted by `binning`.


The user might want to execute the different parts of the program sequentally.
In order to ease and automatize this process, the `run.sh` shell script is
included. This will log the program's execution and move the output files into a
separate folder.

Edit the `input_parameters.in` file before running the simulation to set the
system parameters. There's no need to recompile the program once this file is
edited, as it is read at runtime.

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

In order to visualize the results of the simulation, several Python plotting
scripts are provided. These scripts have the following dependencies:

- [Python 3](https://www.python.org/)
  - [NumPy](https://numpy.org/)
  - [Matplotlib](https://matplotlib.org/)
  - [imageio](https://imageio.readthedocs.io/en/stable/)
- [cm-super TeX font package](https://ctan.org/pkg/cm-super)

## Project structure

The project consists of a main file (`vdw_gas.f90`) and a subroutines file
(`subroutines.f90`) that sources the contents of the `include/` directory, which
stores all the different subroutines used in the code, separated into different
files by topic or task.

A single module called `subroutines` is used due to the fact that different
subroutines might interact with each other.

As previously mentioned, system parameters are set under `input_parameters.in`,
which is read at runtime.

The project tree has the following structure (excluding git-related files, such
as this README).

```
.
|-- include/
|   |-- forces.f90          <- Lennard-Jones forces computation.
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
|-- jackknife.f90           <- Jackknife statistical analysis of the results.
|-- run.sh                  <- Wrapper script for executing the full program.
|-- subroutines.f90         <- Subroutines module. Sources 'include/' contents.
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
