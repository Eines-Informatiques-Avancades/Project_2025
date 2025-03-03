# Molecular Dynamics Simulation of a Van der Waals Gas

Simple molecular dynamics program that deal with a Van der Waals gas of
particles. This program was developed as a part of the final assignment of the
Advanced Computer Tools subject.

## Compiling and running the simulation

The project can be cloned and compiled by running:

```
$ git clone https://github.com/Eines-Informatiques-Avancades/Project_2025.git
$ cd Project_2025
$ make
```

This will produce the `vdw_gas` executable file, which can be ran by:

```
$ ./vdw_gas
```

## Dependencies

- A Fortran compiler (E.g.: [GNU Fortran](https://gcc.gnu.org/fortran/))
- `make`

## Project structure

The project consists of a main file (`vdw_gas.f90`) and a subroutines file
(`subroutines.f90`) that sources the contents of the `include/` directory, which
stores all the different subroutines used in the code, separated into different
files by topic or task.

A single module called `subroutines` is used due to the fact that different
subroutines might interact with each other.

## Credits and contributors

- Ricard Rodríguez: system initialization, periodic boundary condition, project
    coordination
- Oriol Miró: integrators
- Alejandro Díaz: forces
- Joan Serrano: statistical analysis
- Huang Haoyu: post-trajectory analysis
- Itziar Rabal: final testing
