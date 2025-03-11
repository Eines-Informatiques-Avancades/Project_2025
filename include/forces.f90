!
! forces.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Alejandro Díaz, Ricard Rodríguez
!
! Compute forces derived from a Lennard-Jones potential.
!

subroutine compute_forces(part_num, positions, forces, system_size, cutoff)
    implicit none
    
    integer, intent(in) :: part_num
    real, intent(in) :: positions(3, part_num)
    real, intent(in) :: system_size, cutoff
    real, intent(out) :: forces(3, part_num)
    
    integer :: i, j, k
    real :: r, r2, f
    real :: r_vec(3)

    ! In reduced units: epsilon = sigma = 1.
    forces = 0.0

    do i = 1, part_num - 1
        do j = i + 1, part_num
            ! Compute distance between particles i and j.
            r_vec(1) = positions(1, i) - positions(1, j)
            r_vec(2) = positions(2, i) - positions(2, j)
            r_vec(3) = positions(3, i) - positions(3, j)
            
            ! Apply PBC
            do k = 1, 3
                r_vec(k) = pbc(r_vec(k), system_size)
            end do
            
            r2 = r_vec(1)**2 + r_vec(2)**2 + r_vec(3)**2
            r = sqrt(r2)

            if (r < cutoff .and. r > 0.0) then
                ! Force related to a Lennard-Jones potential.
                f = 48.0 / (r**14) - 24.0 / (r**8)
                                
                ! Update forces.
                forces(1, i) = forces(1, i) + f * r_vec(1)
                forces(2, i) = forces(2, i) + f * r_vec(2)
                forces(3, i) = forces(3, i) + f * r_vec(3)
                
                forces(1, j) = forces(1, j) - f * r_vec(1)
                forces(2, j) = forces(2, j) - f * r_vec(2)
                forces(3, j) = forces(3, j) - f * r_vec(3)
            end if
        end do
    end do
end subroutine compute_forces
