!
! read_input.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Itziar Rabal
!
! Subroutine that reads external input files that follow the structure of test_input.dat
!

subroutine read_input(part_num, system_size, lattice_type, timestep, step_num)
    implicit none
    integer, intent(out) :: part_num, step_num
    real, intent(out) :: system_size, timestep
    character(6), intent(out) :: lattice_type

    character(70) :: filename
    integer :: io, nlines, i
    character(70), allocatable :: line(:)
    
    nlines = 0

    print*, "enter input file"
    read(*,*)filename

    open(10, file = filename, status = "old", iostat = io)
    if (io /= 0) then
        print*, "error reading input file"
        stop
    endif
    
    do 
        read(10,*, iostat = io)
        if (io /= 0) exit
        nlines = nlines + 1
    enddo
    
    rewind(10)

    allocate(line(nlines))

    do i =  1, nlines
        read(10,"(A)")line(i)
    enddo

    read(line(3)(17:20), "(I3)")part_num
    read(line(4)(17:20), "(I3)")system_size
    read(line(5)(17:19), "(A)")lattice_type
    read(line(6)(17:23), "(F12.6)")timestep
    read(line(7)(17:23), "(I6)")step_num
    
    deallocate(line)

    close(10)
end subroutine read_input
