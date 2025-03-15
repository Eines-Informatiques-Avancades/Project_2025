!
! jackknife.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Ricard Rodriguez, Joan Serrano
!
! Jackknife program corresponding to the final project for the Advanced
! Computer Tools subject.
!

program jackknife
    implicit none

    integer(8) :: i, j, bin_size, ios, line_num, header_lines, num_values
    integer(8), allocatable :: bin_size_values(:)
    real(8) :: mean_x2, mean_x_2, diff_x
    real(8) :: fm_x, fv_x, fe_x  ! Jackknife estimators
    real(8), allocatable :: x(:), x2(:), xj(:), x2j(:)  ! Data arrays
    character(80) :: input_datafile, output_file, line

    print *
    print *, repeat('-', 80)
    print *, 'Jackknife data analysis for Lennard-Jones Binning'
    print *, repeat('-', 80)
    print *

    ! Input file
    input_datafile = 'binning_lj.out'
    print *, 'Reading data from ', trim(adjustl(input_datafile))

    ! Number of header lines
    header_lines = 1

    ! Open file and count lines
    open(4, file = input_datafile, status = 'old')
    do i = 1, header_lines
        read(4, fmt = '(A)')  ! Skip header
    end do

    line_num = 0
    do
        read(4, fmt = '(A)', iostat = ios) line
        if (ios /= 0) exit
        line_num = line_num + 1
    end do

    allocate(bin_size_values(line_num), x(line_num), x2(line_num))

    ! Read data
    rewind(4)
    do i = 1, header_lines
        read(4, fmt = '(A)')  ! Skip header again
    end do

    do i = 1, line_num
        read(4, *) bin_size_values(i), x(i), x2(i)
    end do

    close(4)

    output_file = 'jackknife_lj.out'
    open(12, file = output_file, action = 'write', status = 'replace')
    write(12, '(A)') '# Bin size   |   Jackknife mean (X)     Variance (X)    Error bar (X)'

    ! Perform jackknife for the read data
    allocate(xj(line_num), x2j(line_num))

    do j = 1, line_num
        bin_size = bin_size_values(j)
        call datjack(line_num, x, xj)
        call datjack(line_num, x2, x2j)

        ! Calculate <X^2>
        mean_x2 = sum(x2j(:))/line_num
        ! Calculate <X>^2
        mean_x_2 = (sum(xj(:))/line_num)**2
        ! Difference <X^2> - <X>^2
        diff_x = mean_x2 - mean_x_2

        ! Jackknife estimators
        call stebj0(line_num, xj, fm_x, fv_x, fe_x)

        print *, repeat('-', 30)
        print *, 'm = ', bin_size
        print *, '<X^2> = ', mean_x2
        print *, '<X>^2 = ', mean_x_2
        print *, '<X^2> - <X>^2 = ', diff_x
        print *, 'Jackknife mean (X) = ', fm_x
        print *, 'Jackknife variance (X) = ', fv_x
        print *, 'Jackknife error bar (X) = ', fe_x

        ! Write results
        write(12, '(A)') repeat('-', 30)
        write(12, '(A, I10)') 'm = ', bin_size
        write(12, '(A, ES20.10)') '<X^2> = ', mean_x2
        write(12, '(A, ES20.10)') '<X>^2 = ', mean_x_2
        write(12, '(A, ES20.10)') '<X^2> - <X>^2 = ', diff_x
        write(12, '(A, ES20.10)') 'Jackknife mean (X) = ', fm_x
        write(12, '(A, ES20.10)') 'Jackknife variance (X) = ', fv_x
        write(12, '(A, ES20.10)') 'Jackknife error bar (X) = ', fe_x
    end do
    close(12)

    deallocate(bin_size_values, x, x2, xj, x2j)

contains

    subroutine datjack(n, x, xj)
        implicit none
        integer(8), intent(in) :: n
        real(8), intent(in) :: x(n)
        real(8), intent(out) :: xj(n)
        integer(8) :: i
        real(8) :: xsum, factor
        
        xsum = sum(x)
        factor = 1.0/real(n - 1, 8)
        do i = 1, n
            xj(i) = factor * (xsum - x(i))
        end do
    end subroutine datjack

    subroutine stebj0(n, fj, fm, fv, fe)
        implicit none
        integer(8), intent(in) :: n
        real(8), intent(in) :: fj(n)
        real(8), intent(out) :: fm, fv, fe
        integer(8) :: i
        
        fm = sum(fj) / n
        fv = sum((fj - fm)**2) * (n - 1) / n
        fe = sqrt(fv / n)
    end subroutine stebj0

end program jackknife
