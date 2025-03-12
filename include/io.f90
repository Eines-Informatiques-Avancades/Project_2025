!
! io.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Ricard Rodriguez
!
! Contains all the subroutines related to writing to or reading from external
! files.
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
