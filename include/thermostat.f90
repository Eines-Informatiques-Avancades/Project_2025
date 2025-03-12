!
! thermostat.f90
! Molecular Dynamics simulation of 3D Van der Waals gas
! Oriol Miro
!
! Andersen thermostat for Molecular Dynamics.
! A 3-dimensional cubic system is assumed.
!

! Andersen thermostat, updates the velocity magnitude and direction of a random number of particles to one that follows the Maxwell-Boltzmann distribution.

subroutine andersen_thermostat(part_num, temperature, collision_frequence, velocities)
    implicit none

    integer, intent(in) :: part_num
    real, allocatable, intent(inout) :: velocities
    real, intent(in) :: temperature, collision frequence
    real :: sigma
    integer :: i

    allocate(velocities(part_num, 3))
    sigma = sqrt(temperature)

    do i = 1, part_num
        if (random_number() < collision_frequence) then
                velocities(i,1) = sigma * random_gaussian()
                velocities(i,2) = sigma * random_gaussian()
                velocities(i,3) = sigma * random_gaussian()
        end if
    end do

end subroutine andersen_thermostat


! Function that returns a random number following the gaussian distribution.

function random_gaussian() result(r)
    implicit none
    real(8) :: r, u1, u2
    call random_number(u1)
    call random_number(u2)
    r = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * 3.141592653589793d0 * u2)
end function random_gaussian