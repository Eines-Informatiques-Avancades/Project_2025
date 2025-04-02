program test_main
    use mpi
    implicit none

    integer, parameter :: part_num = 125
    real(8), parameter :: system_size = 10.0

    integer :: ierr
    integer(kind=MPI_INTEGER_KIND) :: rank, nproc
    integer :: actual_steps
    real(8), allocatable :: x(:, :), y(:, :), z(:, :), time(:)
    character(50) :: positions_file, rdf_file, rmsd_file

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    positions_file = 'positionstest.xyz'
    rdf_file = 'rdf_output.dat'
    rmsd_file = 'rmsd_output.dat'

    call read_trajectory(positions_file, x, y, z, time, actual_steps)
    call compute_rdf(x, y, z, rdf_file, actual_steps)
    call compute_rmsd(x, y, z, time, rmsd_file)

    deallocate(x, y, z, time)
    call MPI_Finalize(ierr)

    if (rank == 0) then
        print *, '✅ All analysis completed successfully.'
    end if

contains

subroutine read_trajectory(positions_file, x, y, z, time, actual_steps)
    use mpi
    implicit none
    character(50), intent(in) :: positions_file
    real(8), allocatable, intent(out) :: x(:, :), y(:, :), z(:, :), time(:)
    integer, intent(out) :: actual_steps

    integer :: part, step, ios, n, ierr, count
    integer, parameter :: part_num = 125
    real(8), parameter :: system_size = 6.0
    integer(kind=MPI_INTEGER_KIND) :: rank, nproc
    integer :: i
    character(50) :: atom_type
    real :: t
    character(len=1000) :: line

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    if (rank == 0) then
        print *, '🔍 Scanning trajectory file to determine actual step count...'
        open(4, file = positions_file, status = 'old', action = 'read')
        actual_steps = 0
        do
            read(4, '(A)', iostat=ios) line
            if (ios /= 0) exit
            actual_steps = actual_steps + 1
            do count = 1, part_num + 1
                read(4, '(A)', iostat=ios)
                if (ios /= 0) exit
            end do
        end do
        close(4)
        print *, '📊 Detected steps:', actual_steps
    end if

    call MPI_Bcast(actual_steps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    if (rank == 0) then
        print *, '📦 Allocating arrays: part_num =', part_num, ', steps =', actual_steps
        print *, 'Estimated memory per array: ', 8.d0 * part_num * actual_steps / 1024 / 1024, 'MB'
    end if

    allocate(x(part_num, actual_steps), y(part_num, actual_steps), z(part_num, actual_steps), time(actual_steps))

    if (rank == 0) then
        open(4, file = positions_file, status = 'old', action = 'read')
        do step = 1, actual_steps
            read(4, *) n
            read(4, *) t
            do part = 1, part_num
                read(4, *) atom_type, x(part, step), y(part, step), z(part, step)
            end do
            time(step) = t
        end do
        close(4)
        print *, '✅ Finished reading trajectory.'
        print *, 'atom:', atom_type
    end if

    call MPI_Bcast(x, part_num * actual_steps, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(y, part_num * actual_steps, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(z, part_num * actual_steps, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(time, actual_steps, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if (rank /= 0) then
        print *, '--- Rank', rank, 'received data sample ---'
        print *, 'Time step 1 time value:', time(1)
        print *, 'First 3 atoms at step 1:'
        do i = 1, min(3, part_num)
            print *, 'Atom', i, ':', x(i,1), y(i,1), z(i,1)
        end do
        print *, '----------------------------------------'
    end if
end subroutine read_trajectory

! (compute_rdf 和 compute_rmsd 保持一致，只需确保 rank, nproc 定义为 MPI_INTEGER_KIND)



subroutine compute_rdf(x, y, z, rdf_file, step_num)
    use mpi
    implicit none
    real(8), intent(in) :: x(:, :), y(:, :), z(:, :)
    character(50), intent(in) :: rdf_file
    integer, intent(in) :: step_num

    integer :: i, j, k, time_index
    real(8), parameter :: dr = 0.01
    real(8), parameter :: system_size = 10.0
    integer, parameter :: part_num = 125
    logical, parameter :: debug = .true.   ! ✅ 调试开关

    real(8) :: maximum_radius, volume, density
    integer :: bins
    real(8), allocatable :: h(:), rdf(:), r_values(:), h_total(:)
    real(8) :: r, r_sq, dx, dy, dz, dv, r_lo, r_hi, const, nid
    integer :: bin_index, size, ierr, rank
    integer :: frame_chunk, start, fin

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    maximum_radius = system_size / 2
    bins = int(maximum_radius / dr)
    volume = system_size**3
    density = part_num / volume

    allocate(h(bins), rdf(bins), r_values(bins), h_total(bins))
    h = 0.0
    h_total = 0.0

    frame_chunk = step_num / size
    start = rank * frame_chunk + 1
    fin = start + frame_chunk - 1
    if (rank == size - 1) then
        fin = step_num
    endif
    if (debug) print *, '🔧 Rank', rank, 'handles steps:', start, 'to', fin

    if (start > fin) return  ! 🔒 防止空循环

    do time_index = start, fin
        do i = 1, part_num - 1
            do j = i + 1, part_num
                dx = x(j, time_index) - x(i, time_index)
                dy = y(j, time_index) - y(i, time_index)
                dz = z(j, time_index) - z(i, time_index)

                dx = dx - system_size * anint(dx / system_size)
                dy = dy - system_size * anint(dy / system_size)
                dz = dz - system_size * anint(dz / system_size)

                r_sq = dx**2 + dy**2 + dz**2
                r = sqrt(r_sq)
		
		

                if (debug .and. time_index == start .and. i == 1 .and. j <= 3) then
                    print *, '📏 Rank', rank, 'Step', time_index, 'Pair (', i, ',', j, ')'
                    print *, '    dx =', dx, 'dy =', dy, 'dz =', dz, 'r =', r
                end if

                if (r < maximum_radius) then
                    bin_index = floor(r / dr) + 1
                    if (debug .and. time_index == start .and. i == 1 .and. j <= 3) then
                        print *, '📊 bin_index =', bin_index, 'before h =', h(bin_index)
                    end if
                    h(bin_index) = h(bin_index) + 2
                    if (debug .and. time_index == start .and. i == 1 .and. j <= 3) then
                        print *, 'after h =', h(bin_index)
                    end if
                end if

            end do
        end do
    end do

    if (debug) print *, '🧮 Rank', rank, 'local h() total =', sum(h)

    call MPI_Reduce(h, h_total, bins, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (rank == 0) then
        if (debug) print *, '📦 Total h_total sum =', sum(h_total)
        const = 4.0 * 3.14159265358979 / 3.0
        do k = 1, bins
            r_lo = (k - 1) * dr
            r_hi = r_lo + dr
            dv = const * (r_hi**3 - r_lo**3)
            nid = dv * density
            rdf(k) = h_total(k) / (part_num * step_num * nid)
            r_values(k) = (k - 0.5) * dr

            if (debug .and. k <= 5) then
                print *, '📐 bin', k, ': r_lo =', r_lo, 'r_hi =', r_hi, 'dv =', dv, 'nid =', nid
                print *, '    rdf =', rdf(k), 'h_total =', h_total(k)
            end if
        end do

        open(12, file = rdf_file, status = 'replace')
        do k = 1, bins
            write(12, '(F10.5, F15.8)') r_values(k), rdf(k)
        end do
        close(12)
        print *, '✅ RDF calculation completed and saved to ', rdf_file
    end if

    deallocate(h, h_total, rdf, r_values)
end subroutine compute_rdf

subroutine compute_rmsd(x, y, z, time, rmsd_file)
    use mpi
    implicit none
    real(8), intent(in) :: x(:, :), y(:, :), z(:, :), time(:)
    character(50), intent(in) :: rmsd_file

    integer :: i, j, nproc, rank, ierr, partial_steps, frame_chunk, start, fin, actual_step
    integer :: step_num, part_num
    real(8), allocatable :: total_rmsd(:), partial_rmsd(:)
    real(8) :: dx, dy, dz, sum_sq
    integer, allocatable :: counts(:), displs(:)

    part_num = size(x, 1)
    step_num = size(x, 2)

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

    frame_chunk = step_num / nproc
    start = rank * frame_chunk + 1
    fin = start + frame_chunk - 1
    if (rank == nproc - 1) fin = step_num

    partial_steps = fin - start + 1
    allocate(partial_rmsd(partial_steps))

    if (rank == 0) then
        allocate(total_rmsd(step_num))
        allocate(counts(nproc), displs(nproc))
        do i = 0, nproc - 1
            if (i == nproc - 1) then
                counts(i + 1) = step_num - i * frame_chunk
            else
                counts(i + 1) = frame_chunk
            end if
        end do
        displs(1) = 0
        do i = 2, nproc
            displs(i) = displs(i - 1) + counts(i - 1)
        end do
    end if

    do j = 1, partial_steps
        actual_step = start + j - 1
        sum_sq = 0.0
        do i = 1, part_num
            dx = x(i, actual_step) - x(i, 1)
            dy = y(i, actual_step) - y(i, 1)
            dz = z(i, actual_step) - z(i, 1)
            sum_sq = sum_sq + (dx**2 + dy**2 + dz**2)
        end do
        partial_rmsd(j) = sqrt(sum_sq / part_num)
    end do

    call MPI_Gatherv(partial_rmsd, partial_steps, MPI_DOUBLE_PRECISION, &
                     total_rmsd, counts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if (rank == 0) then
        open(12, file = rmsd_file, status = 'replace')
        do j = 1, step_num
            write(12, *) time(j), total_rmsd(j)
        end do
        close(12)
        print *, '📦 RMSD calculation completed and saved to ', rmsd_file
        deallocate(total_rmsd, counts, displs)
    end if

    deallocate(partial_rmsd)
end subroutine compute_rmsd

end program

