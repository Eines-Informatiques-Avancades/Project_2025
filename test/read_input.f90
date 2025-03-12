
program main
    use, intrinsic :: iso_fortran_env, only : dp => real64

    implicit none
    integer :: part_num, system_size
    character(len = 2) :: lattice_type

    call read_input(part_num, system_size, lattice_type)

    print("(I3,x,I3,x,A)"), part_num, system_size, lattice_type
contains

subroutine read_input(part_num, system_size, lattice_type)
    implicit none
    integer, intent(out) :: part_num, system_size
    character(len = 2), intent(out) :: lattice_type

    integer :: io
    character(len = 14) :: filename
    character(len = 70) :: line

    filename = "test_input.dat"

    open(10, file = filename, iostat = io)
    if (io == 0) then
        print*, "testing file opened for reading"
    endif

    
    read(10,*)  !skip title line
    read(10,*)  !skip title line
    read(10,*)  !skip title line

    read(10,"(A)")line
    print*, line

    read(line(17:20), "(I3)")part_num
    print*, part_num

    read(10,"(A)")line
    print*, line

    read(line(17:20), "(I3)")system_size
    print*, system_size

    read(10,"(A)")line
    print*, line

    read(line(17:19), "(A)")lattice_type
    print*, lattice_type

    close(10)
end subroutine read_input

end program main
