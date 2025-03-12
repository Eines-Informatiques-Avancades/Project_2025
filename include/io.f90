!
! io.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Ricard Rodriguez, Itziar Rabal
!
! Contains all the subroutines related to writing to or reading from external
! files.
!

! Ricard Rodriguez
!
! Write the coordinates of all particles of the system in the current point of
! time in an specified external file with an XYZ format.
! The external file must already exist and results will be appended to it.
subroutine write_positions_xyz(part_num, time, positions, output_file)
    implicit none

    integer, intent(in) :: part_num
    real, intent(in) :: time
    character(50), intent(in) :: output_file
    real, allocatable, intent(in) :: positions(:, :)

    integer :: i

    open(4, file = output_file, status = 'old', action = 'write', access = 'append')

    do i = 1, part_num
        write(4, *) time, positions(i, :)
    end do

    close(4)
end subroutine

! Itziar Rabal
!
! Reads external input files that follow the structure of test_input.dat
subroutine read_input(input_file, part_num, system_size, lattice_type, timestep, step_num)
    implicit none

    character(50), intent(in) :: input_file
    integer, intent(out) :: part_num, step_num
    real, intent(out) :: system_size, timestep
    character(6), intent(out) :: lattice_type

    integer :: ios, nlines, i
    character(70), allocatable :: line(:)

    nlines = 0

    open(10, file = input_file, status = 'old', iostat = ios)
    if (ios /= 0) then
        print *, 'Error reading input file ', input_file
        stop
    endif

    do
        read(10, *, iostat = ios)
        if (ios /= 0) exit
        nlines = nlines + 1
    enddo

    rewind(10)

    allocate(line(nlines))

    do i = 1, nlines
        read(10, '(A)') line(i)
    enddo

    read(line(3)(17:20), '(I3)') part_num
    read(line(4)(17:20), '(F12.6)') system_size
    read(line(5)(17:19), '(A)') lattice_type
    read(line(6)(17:23), '(F12.6)') timestep
    read(line(7)(17:23), '(I6)') step_num

    deallocate(line)

    close(10)
end subroutine read_input
