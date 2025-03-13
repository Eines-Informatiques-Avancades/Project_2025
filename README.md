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

This will produce the `vdw_gas` executable file, which can be ran by:

```
$ ./vdw_gas
```

Edit the `input_parameters.in` file before running the simulation to set the
system parameters. There's no need to recompile the program once this file is
edited, as it is read at runtime.

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
|   |-- forces.f90
|   |-- geometry.f90
|   |-- initial_conf.f90
|   |-- integrators.f90
|   |-- io.f90
|   |-- thermodynamics.f90
|   `-- thermostat.f90
|-- plot/
|   `-- ...
|-- test/
|   `-- ...
|-- Makefile
|-- input_parameters.in
|-- subroutines.f90
`-- vdw_gas.f90
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
