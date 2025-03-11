!
! vdw_gas
! Molecular Dynamics Simulation of a Van der Waals Gas
! Ricard Rodriguez
!
! Main program corresponding to the "Molecular Dynamics Simulation of a Van der
! Waals Gas" project for the Advanced Computer Tools subject. Please check the
! rest of the files containing the module(s) necessary for the program
! execution.
!

program vdw_gas
    use subroutines

    implicit none

    integer :: part_num, time
    real :: part_density, system_size, volume, cutoff
    real, allocatable :: positions(:, :), forces(:, :)
    character(6) :: lattice_type
    character(50) :: positions_file

    ! Generate initial lattice.
    lattice_type = 'SC'
    part_num = 125
    system_size = 500

    volume = system_size**(3.)  ! System is a cubic box.
    cutoff = 0.5*system_size    ! Cutoff radius for molecular interactions.

    call gen_initial_conf(lattice_type, system_size, part_num, part_density, positions)

    ! Center initial config at the origin of coordinates.
    call apply_pbc(positions, system_size)

    positions_file = 'positions.xyz'

    ! Create a new positions_file or replace the existing one.
    open(4, file = positions_file, status = 'replace')
    write(4, *) '# time, particle coordinates'
    close(4)

    do time = 1, 5
        call write_positions_xyz(part_num, time, positions, positions_file)
    end do
end program vdw_gas
