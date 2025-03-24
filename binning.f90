!
! binning.f90
! Molecular Dynamics Simulation of a Van der Waals Gas
! Joan Serrano, Ricard Rodríguez
!
! Binning program corresponding to the final project for the Advanced Computer
! Tools subject.
!

program binning
    implicit none

    integer :: i, j
    integer :: n_rows, n_cols
    integer :: ios, unit_in, unit_out
    character(len=256) :: input_file, output_file
    character(len=1024) :: line, header_line

    ! Máximo posible de columnas que leeremos (ajústalo si necesitas más)
    integer, parameter :: MAX_COLS = 50

    ! Arreglo para almacenar los nombres de columna (hasta MAX_COLS)
    character(len=64), dimension(MAX_COLS) :: col_names

    ! Aquí guardaremos los datos en un arreglo 2D:
    ! data(row, col)
    real(8), allocatable :: data(:, :)

    !----------------------------------------------------------------------
    ! 1. Definir el archivo de entrada y de salida
    !----------------------------------------------------------------------
    input_file  = 'thermodynamics.dat'
    output_file = 'binning_' // trim(input_file) ! Archivo de salida con resultados

    print *, 'Binning for multi-column file: ', trim(input_file)
    print *, 'Resultados se guardarán en:    ', trim(output_file)

    !----------------------------------------------------------------------
    ! 2. Abrir el archivo de entrada para:
    !    - Leer la línea de encabezado
    !    - Determinar el número de columnas (n_cols)
    !    - Contar cuántas filas de datos hay (n_rows)
    !----------------------------------------------------------------------
    unit_in = 10
    open(unit_in, file=input_file, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Error al abrir el archivo: ', trim(input_file)
        stop
    end if

    ! Lee la primera línea como encabezado
    read(unit_in, '(A)', iostat=ios) header_line
    if (ios /= 0) then
        print *, 'Error leyendo la cabecera del archivo.'
        stop
    end if

    ! ---------------------------------------------------------------------
    ! 2.1 Parsear la línea de encabezado para encontrar cuántas columnas hay
    !     y guardar sus nombres en col_names(1..n_cols).
    !     Suponemos que el archivo está separado por comas o espacios.
    !     Se ha agregado una rutina para remover el carácter '#' si está presente.
    ! ---------------------------------------------------------------------
    n_cols = 0
    call parse_header(header_line, col_names, n_cols)

    if (n_cols < 1) then
        print *, 'No se encontraron columnas en la cabecera.'
        stop
    end if

    ! ---------------------------------------------------------------------
    ! 2.2 Contar cuántas filas de datos hay
    !     (Leemos línea a línea hasta el final)
    ! ---------------------------------------------------------------------
    n_rows = 0
    do
        read(unit_in, '(A)', iostat=ios) line
        if (ios /= 0) exit
        n_rows = n_rows + 1
    end do

    close(unit_in)

    if (n_rows == 0) then
        print *, 'No hay filas de datos en el archivo.'
        stop
    end if

    !----------------------------------------------------------------------
    ! 3. Reservar memoria para los datos: data(n_rows, n_cols)
    !----------------------------------------------------------------------
    allocate(data(n_rows, n_cols))

    !----------------------------------------------------------------------
    ! 4. Volver a abrir el archivo y leer nuevamente:
    !    - Saltamos la línea de encabezado
    !    - Leemos cada fila de datos numéricos
    !----------------------------------------------------------------------
    open(unit_in, file=input_file, status='old', action='read', iostat=ios)
    read(unit_in, '(A)', iostat=ios) header_line  ! Saltar cabecera

    do i = 1, n_rows
        ! Leemos la línea completa como texto
        read(unit_in, '(A)', iostat=ios) line
        if (ios /= 0) then
            print *, 'Error leyendo línea ', i
            exit
        end if
        ! Parseamos esta línea para extraer n_cols valores reales
        call parse_data_line(line, data(i,:), n_cols)
    end do
    close(unit_in)

    !----------------------------------------------------------------------
    ! 5. Abrir el archivo de salida y escribir encabezado
    !----------------------------------------------------------------------
    unit_out = 20
    open(unit_out, file=output_file, status='replace', action='write')
    write(unit_out, '(A)') '# Binning results for multi-column data'
    write(unit_out, '(A)') '# Formato: ColName, BinSize, Mean, Variance, StdDev'

    !----------------------------------------------------------------------
    ! 6. Para cada columna, aplicamos el binning
    !----------------------------------------------------------------------
    do j = 2, n_cols
        call do_binning_for_column(data(:, j), n_rows, col_names(j), unit_out)
    end do

    close(unit_out)
    deallocate(data)

    print *, 'Proceso de binning completado. Revisa ', trim(output_file)

contains
    !======================================================================
    ! SUBRUTINA: parse_header
    !   - Recibe la línea de cabecera (header_line)
    !   - Extrae nombres de columna y cuenta cuántos hay
    !   - Si la cabecera comienza con '#' se elimina.
    !======================================================================
    subroutine parse_header(header_line, col_names, n_cols)
        implicit none
        character(len=*), intent(in) :: header_line
        character(len=64), dimension(:), intent(out) :: col_names
        integer, intent(out) :: n_cols

        integer :: pos, start, len_line
        character(len=64) :: token
        character(len=1024) :: tmp_line

        tmp_line = header_line
        ! Remover espacios a la izquierda
        tmp_line = adjustl(tmp_line)
        ! Si la línea comienza con '#', eliminarla
        if (len_trim(tmp_line) > 0 .and. tmp_line(1:1) == '#') then
            tmp_line = adjustl(tmp_line(2:))
        end if

        len_line = len_trim(tmp_line)
        n_cols = 0
        pos = 1

        ! Se asume que el separador es espacio o coma.
        ! Un truco sencillo es reemplazar comas por espacios.
        call replace_char(tmp_line, ',', ' ')

        ! Bucle para extraer "tokens" (nombres de columna)
        do while (pos <= len_line)
            ! Saltar espacios iniciales
            do while (pos <= len_line .and. tmp_line(pos:pos) == ' ')
                pos = pos + 1
            end do
            if (pos > len_line) exit

            start = pos
            ! Avanzar hasta espacio o fin de línea
            do while (pos <= len_line .and. tmp_line(pos:pos) /= ' ')
                pos = pos + 1
            end do

            ! Extraer el token
            token = tmp_line(start:pos-1)
            n_cols = n_cols + 1
            if (n_cols <= size(col_names)) then
                col_names(n_cols) = token
            else
                print *, 'ADVERTENCIA: más columnas de las esperadas.'
                exit
            end if
        end do
    end subroutine parse_header

    !======================================================================
    ! SUBRUTINA: parse_data_line
    !   - Parsea una línea con n_cols valores reales
    !   - Asume que están separados por espacios (o comas, que convertimos)
    !======================================================================
    subroutine parse_data_line(line, arr, n_cols)
        implicit none
        character(len=*), intent(in) :: line
        real(8), intent(out) :: arr(:)
        integer, intent(in) :: n_cols

        character(len=1024) :: tmp_line
        integer :: ios

        ! Reemplazamos comas por espacios (por si el archivo es CSV real)
        tmp_line = line
        call replace_char(tmp_line, ',', ' ')

        ! Ahora leemos con list-directed I/O (asume n_cols valores reales)
        read(tmp_line, *, iostat=ios) (arr(j), j=1,n_cols)
        if (ios /= 0) then
            print *, 'Error parseando la línea: ', trim(line)
        end if
    end subroutine parse_data_line

    !======================================================================
    ! SUBRUTINA: replace_char
    !   - Reemplaza todas las ocurrencias de un caracter por otro
    !======================================================================
    subroutine replace_char(str, ch_old, ch_new)
        implicit none
        character(len=*), intent(inout) :: str
        character(len=1), intent(in) :: ch_old, ch_new
        integer :: i, l

        l = len_trim(str)
        do i = 1, l
            if (str(i:i) == ch_old) str(i:i) = ch_new
        end do
    end subroutine replace_char

    !======================================================================
    ! SUBRUTINA: do_binning_for_column
    !   - Recibe un arreglo 1D con los datos de la columna
    !   - Calcula binning y escribe resultados
    !======================================================================
    subroutine do_binning_for_column(col_data, n_values, col_name, unit_out)
        implicit none
        real(8), intent(in) :: col_data(:)
        integer, intent(in) :: n_values
        character(len=*), intent(in) :: col_name
        integer, intent(in) :: unit_out

        integer(8) :: i, j, k
        integer(8) :: max_num_bins, num_bins, bin_size, bin_start, bin_end
        real(8)    :: sum_bin_variable, var_variable
        real(8)    :: sum_var_variable, mean_variable
        real(8), allocatable :: avg_bin(:)

        ! 1. Determinar el número máximo de niveles de binning
        max_num_bins = 1
        do while (n_values / (2**max_num_bins) > 100)
            max_num_bins = max_num_bins + 1
        end do

        ! 2. Para cada nivel de binning
        do i = 0, max_num_bins - 1
            bin_size = 2**i
            num_bins = n_values / bin_size

            allocate(avg_bin(num_bins))

            ! Calcular la media de cada bin
            do j = 1, num_bins
                bin_start = (j - 1)*bin_size + 1
                bin_end   = j*bin_size
                sum_bin_variable = 0.0d0

                do k = bin_start, bin_end
                    sum_bin_variable = sum_bin_variable + col_data(k)
                end do

                avg_bin(j) = sum_bin_variable / bin_size
            end do

            ! Calcular la media global de los valores binned
            mean_variable = sum(avg_bin) / num_bins

            ! Calcular la varianza
            sum_var_variable = 0.0d0
            do j = 1, num_bins
                sum_var_variable = sum_var_variable + (avg_bin(j) - mean_variable)**2
            end do

            ! Esta es la fórmula original del código
            var_variable = sum_var_variable / real((num_bins - 1)*num_bins, 8)

            ! Escribimos la línea en el archivo de salida
            write(unit_out, '(A, 1X, I6, 1X, ES15.8, 1X, ES15.8, 1X, ES15.8)') &
                trim(col_name), bin_size, mean_variable, var_variable, sqrt(var_variable)

            deallocate(avg_bin)
        end do
    end subroutine do_binning_for_column
end program binning
