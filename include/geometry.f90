!
! geometry.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Ricard Rodriguez
!
! Contains all the subroutines dependent or related to the system geometry.
!

module geometry
    use global_vars

    implicit none

    contains
        ! Apply PBC conditions to the positions array.
        ! This subroutine must be run everytime the positions of the particles are
        ! updated.
        subroutine apply_pbc(positions)
            implicit none

            real, allocatable, intent(inout) :: positions(:, :)

            integer :: i, j

            do i = 1, part_num
                do j = 1, 3
                    positions(i, j) = pbc(positions(i, j), system_size)
                end do
            end do
        end subroutine apply_pbc

        ! Apply PBC to a certain distance, given the box size.
        function pbc(distance, box_size)
            implicit none

            real :: distance, box_size, pbc

            if (distance > box_size/2) then
                distance = distance - box_size
            else if (distance < -box_size/2) then
                distance = distance + box_size
            end if

            pbc = distance
        end function pbc
end module geometry
