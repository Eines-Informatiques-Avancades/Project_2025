!=======================================================================
! mpi_binning.f90
! Programa MPI para cálculo de binning con lectura optimizada y
! medición de tiempos en cada sección.
!
! Se espera que el archivo de entrada tenga el siguiente formato:
!   time, lj_potential, kinetic_energy, total_energy, temperature_inst
!   1, -76.3191071, 451.121307, 374.802185, 2.42538333
!   1, -78.2791214, 471.210419, 392.931305, 2.53338933
!   ...
!
! Compilar con: mpif90 binning_par_4.f90 -o binning_par_4 -lm
! Ejecutar con: mpirun -np 8 ./binning_par_4
!=======================================================================

program mpi_binning
    use mpi
    implicit none

    integer :: i, j, ierr, my_rank, nprocs
    integer :: n_rows, n_cols, unit_out
    character(len=256) :: input_file, output_file
    character(len=1024) :: header_line
    real(8) :: t_total_start, t_total_end
    real(8) :: t_read_start, t_read_end
    real(8) :: t_bcast_start, t_bcast_end
    real(8) :: t_binning_start, t_binning_end

    ! Parámetro: número máximo de columnas (se puede aumentar si es necesario)
    integer, parameter :: MAX_COLS = 200
    character(len=64), dimension(MAX_COLS) :: col_names

    real(8), allocatable :: data(:, :)

    ! Variables para lectura optimizada:
    character(len=:), allocatable :: file_buffer
    character(len=1024), allocatable :: file_lines(:)
    integer :: file_size, num_lines

    !-----------------------------------------------------------------------
    ! 1. Inicializar MPI y medir tiempo total
    !-----------------------------------------------------------------------
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    t_total_start = MPI_Wtime()

    !-----------------------------------------------------------------------
    ! 2. (Solo rank 0) Leer el archivo completo y separar en líneas
    !-----------------------------------------------------------------------
    if (my_rank == 0) then
        input_file  = 'thermodynamics.dat'
        output_file = 'binning_' // trim(input_file)
        print *, 'Procesando el archivo: ', trim(input_file)
        t_read_start = MPI_Wtime()
        call read_file_buffer(input_file, file_buffer, file_size, ierr)
        if (ierr /= 0) then
            print *, 'Error al leer el archivo.'
            call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)
        end if
        call split_lines(file_buffer, file_lines, num_lines)
        if (num_lines < 2) then
            print *, 'Error: el archivo debe tener al menos cabecera + 1 línea de datos.'
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! La primera línea es la cabecera
        header_line = file_lines(1)
        call parse_header(header_line, col_names, n_cols)
        if (n_cols < 1) then
            print *, 'Error: cabecera vacía o mal leída.'
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! Las líneas restantes son datos
        n_rows = num_lines - 1
        allocate(data(n_rows, n_cols))
        do i = 1, n_rows
            call parse_data_line(file_lines(i+1), data(i, :), n_cols)
        end do
        t_read_end = MPI_Wtime()
    end if

    !-----------------------------------------------------------------------
    ! 3. Broadcast de dimensiones, datos y cabecera, y medir tiempo de broadcast.
    !-----------------------------------------------------------------------
    call MPI_Bcast(n_rows, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(n_cols, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (my_rank /= 0) then
        allocate(data(n_rows, n_cols))
    end if
    t_bcast_start = MPI_Wtime()
    call MPI_Bcast(data, n_rows*n_cols, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(col_names, MAX_COLS, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    t_bcast_end = MPI_Wtime()

    !-----------------------------------------------------------------------
    ! 4. (Solo rank 0) Abrir archivo de salida.
    !-----------------------------------------------------------------------
    if (my_rank == 0) then
        unit_out = 20
        open(unit_out, file=output_file, status='replace', action='write')
        write(unit_out, '(A)') '# Resultados de binning para datos multicolumna'
        write(unit_out, '(A)') '# Formato: ColName, BinSize, Mean, Variance, StdDev'
    end if

    !-----------------------------------------------------------------------
    ! 5. Ejecutar sección de binning y medir tiempo de cómputo paralelo.
    !-----------------------------------------------------------------------
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    t_binning_start = MPI_Wtime()
    do j = 2, n_cols
        call do_binning_for_column(data(:, j), n_rows, col_names(j), unit_out, my_rank, nprocs)
    end do
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    t_binning_end = MPI_Wtime()

    !-----------------------------------------------------------------------
    ! 6. Medir tiempo total y mostrar todos los tiempos (solo rank 0).
    !-----------------------------------------------------------------------
    if (my_rank == 0) then
        t_total_end = MPI_Wtime()
        print *, "Tiempo total del programa: ", t_total_end - t_total_start, " seg"
        print *, "Tiempo de lectura del archivo: ", t_read_end - t_read_start, " seg"
        print *, "Tiempo de broadcast de datos: ", t_bcast_end - t_bcast_start, " seg"
        print *, "Tiempo de cómputo (binning): ", t_binning_end - t_binning_start, " seg"
        close(unit_out)
        print *, 'Cálculo de binning completado. Resultados en: ', trim(output_file)
    end if

    if (allocated(file_buffer)) deallocate(file_buffer)
    if (allocated(file_lines)) deallocate(file_lines)
    if (allocated(data)) deallocate(data)
    call MPI_Finalize(ierr)
    
contains

    subroutine read_file_buffer(fname, buffer, fsize, ierr)
        implicit none
        character(len=*), intent(in) :: fname
        character(len=:), allocatable, intent(out) :: buffer
        integer, intent(out) :: fsize, ierr
        integer :: unit, ios
        inquire(file=fname, size=fsize)
        allocate(character(len=fsize) :: buffer)
        open(newunit=unit, file=fname, status='old', access='stream', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = ios
            return
        end if
        read(unit) buffer
        close(unit)
        ierr = 0
    end subroutine read_file_buffer

    subroutine split_lines(buffer, lines, num_lines)
        implicit none
        character(len=*), intent(in) :: buffer
        character(len=1024), allocatable, dimension(:) :: lines
        integer, intent(out) :: num_lines
        integer :: i, start, len_buf

        len_buf = len_trim(buffer)
        num_lines = 0
        start = 1
        do i = 1, len_buf
            if (buffer(i:i) == char(10)) then
                num_lines = num_lines + 1
            end if
        end do
        if (buffer(len_buf:len_buf) /= char(10)) num_lines = num_lines + 1
        allocate(lines(num_lines))
        num_lines = 0
        start = 1
        do i = 1, len_buf
            if (buffer(i:i) == char(10)) then
                num_lines = num_lines + 1
                lines(num_lines) = buffer(start:i-1)
                start = i + 1
            end if
        end do
        if (start <= len_buf) then
            num_lines = num_lines + 1
            lines(num_lines) = buffer(start:len_buf)
        end if
    end subroutine split_lines

    subroutine parse_header(line_in, col_names, n_cols)
        implicit none
        character(len=*), intent(in) :: line_in
        character(len=64), dimension(:), intent(out) :: col_names
        integer, intent(out) :: n_cols
        character(len=1024) :: tmp_line
        integer :: pos, len_line, token_len, i
        character(len=64) :: token

        tmp_line = adjustl(line_in)
        if (len_trim(tmp_line) > 0 .and. tmp_line(1:1) == '#') then
            tmp_line = adjustl(tmp_line(2:))
        end if

        ! --- Remove commas from the header line ---
        do i = 1, len(tmp_line)
            if (tmp_line(i:i) == ',') then
                tmp_line(i:i) = ' '
            end if
        end do

        len_line = len_trim(tmp_line)
        n_cols = 0
        pos = 1
        do while (pos <= len_line)
            token_len = index(tmp_line(pos:), ' ')
            if (token_len == 0) then
                token = tmp_line(pos:)
                pos = len_line + 1
            else
                token = tmp_line(pos: pos+token_len-2)
                pos = pos + token_len
            end if
            if (len_trim(token) > 0) then
                n_cols = n_cols + 1
                col_names(n_cols) = token
            end if
        end do
    end subroutine parse_header

    subroutine parse_data_line(line_in, arr, n_cols)
        implicit none
        character(len=*), intent(in) :: line_in
        real(8), intent(out) :: arr(:)
        integer, intent(in) :: n_cols
        character(len=1024) :: tmp_line
        integer :: pos, len_line, token_len, token_count
        character(len=64) :: token
        token_count = 0
        tmp_line = adjustl(line_in)
        len_line = len_trim(tmp_line)
        pos = 1
        do while (pos <= len_line .and. token_count < n_cols)
            token_len = index(tmp_line(pos:), ' ')
            if (token_len == 0) then
                token = tmp_line(pos:)
                pos = len_line + 1
            else
                token = tmp_line(pos: pos+token_len-2)
                pos = pos + token_len
            end if
            if (len_trim(token) > 0) then
                token_count = token_count + 1
                read(token, *) arr(token_count)
            end if
        end do
        if (token_count /= n_cols) then
            print *, 'Error al parsear la línea de datos: ', trim(line_in)
            print *, 'Se esperaban ', n_cols, ' tokens, se encontraron ', token_count
        end if
    end subroutine parse_data_line

    subroutine do_binning_for_column(col_data, n_values, col_name, unit_out, my_rank, nprocs)
        use mpi
        implicit none
        real(8), intent(in) :: col_data(:)
        integer, intent(in) :: n_values
        character(len=*), intent(in) :: col_name
        integer, intent(in) :: unit_out, my_rank, nprocs
        integer(8) :: i, j, k
        integer(8) :: max_num_bins, num_bins, bin_size
        integer(8) :: local_bin_start, local_bin_end, local_count
        real(8) :: global_sum, global_sum2, global_mean, variance, stddev
        real(8) :: local_sum, local_sum2
        real(8), allocatable :: local_avg(:)
        real(8), dimension(2) :: local_vals, global_vals
        integer :: ierr

        max_num_bins = 1
        do while (n_values / (2**max_num_bins) > 10)
            max_num_bins = max_num_bins + 1
        end do

        do i = 0, max_num_bins - 1
            bin_size = 2**i
            num_bins = n_values / bin_size

            local_bin_start = int(my_rank * dble(num_bins) / dble(nprocs)) + 1
            local_bin_end   = int((my_rank+1) * dble(num_bins) / dble(nprocs))
            local_count = local_bin_end - local_bin_start + 1
            if (local_count < 1) cycle

            allocate(local_avg(local_count))
            do j = local_bin_start, local_bin_end
                local_avg(j - local_bin_start + 1) = 0.0d0
                do k = (j - 1)*bin_size + 1, j*bin_size
                    local_avg(j - local_bin_start + 1) = local_avg(j - local_bin_start + 1) + col_data(k)
                end do
                local_avg(j - local_bin_start + 1) = local_avg(j - local_bin_start + 1) / bin_size
            end do

            local_sum = sum(local_avg(1:local_count))
            local_sum2 = sum(local_avg(1:local_count)**2)
            deallocate(local_avg)

            local_vals(1) = local_sum
            local_vals(2) = local_sum2

            call MPI_Allreduce(local_vals, global_vals, 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
            global_sum = global_vals(1)
            global_sum2 = global_vals(2)
            global_mean = global_sum / dble(num_bins)
            variance = (global_sum2 / dble(num_bins)) - global_mean**2
            if (variance < 0.0d0) variance = 0.0d0
            stddev = sqrt(variance)

            if (my_rank == 0) then
                write(unit_out, '(A, 1X, I6, 1X, F15.8, 1X, F15.8, 1X, F15.8)') &
                     trim(col_name), int(bin_size), global_mean, variance, stddev
            end if
        end do
    end subroutine do_binning_for_column

end program mpi_binning
