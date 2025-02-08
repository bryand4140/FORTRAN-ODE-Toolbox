module IO_Toolbox
    use MOD_Select_Kind, only: pv

    implicit none

    private

    public :: write_matrix, read_matrix


contains

subroutine write_matrix(matrix, filename, path, scientific)
    ! This subroutine writes a matrix to a CSV file.

    ! Input:
    !   matrix:     the matrix to write to the csv file
    !   filename:   the name of the csv file to write to. Must have .csv extension. 
    !   path:       optional argument specifying the directory path for the output file
    !   scientific: optional argument to specify whether to use scientific notation
    !               for the output values. 

    !Note on Paths: For windows, use double backslashes in the path name. Do not include the
    !               filename in the path. Example: 'C:\\Users\\username\\Desktop'

    implicit none

    ! Declare the input arguments
    real(pv), intent(in) :: matrix(:,:)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in), optional :: path
    logical, intent(in), optional :: scientific

    ! Declare local variables
    integer :: i, j
    integer :: status
    character(len=30) :: value
    logical :: use_scientific
    character(len=256) :: full_path

    ! Determine if scientific notation should be used
    if (present(scientific)) then
        use_scientific = scientific
    else
        use_scientific = .false.
    end if

    ! Construct the full path
    if (present(path)) then
        ! Check if path ends with separator, if not add it
        if (path(len_trim(path):len_trim(path)) == '/' .or. &
            path(len_trim(path):len_trim(path)) == '\') then
            full_path = trim(path) // trim(filename)
        else
            full_path = trim(path) // '/' // trim(filename)
        end if
    else
        full_path = trim(filename)
    end if

    ! Open the file for writing
    open(unit=10, file=full_path, status='replace', action='write', &
         form='formatted', iostat=status)
    if (status /= 0) then
        print*, 'Error opening file ', full_path
        stop
    end if

    ! Write the matrix to the file
    print*, " " ! Spacer
    write(*,"(A,A)") "Writing matrix to file ", trim(full_path)
    do i = 1, size(matrix, 1)
        do j = 1, size(matrix, 2)
            if (use_scientific) then
                write(value, '(E15.6)') matrix(i,j)
            else
                write(value, '(F15.6)') matrix(i,j)
            end if
            if (j < size(matrix, 2)) then
                write(10, '(A)', advance='no') trim(value) // ','
            else
                write(10, '(A)', advance='no') trim(value)
            end if
        end do
        write(10, *)
    end do

    ! Close the file
    close(10)
    write(*,"(A)") "File write Complete."
    print*, " " ! Spacer

end subroutine write_matrix


subroutine read_matrix(matrix, filename)
    ! This subroutine reads a CSV file into a matrix
    implicit none

    ! Declare the input arguments
    character(len=*), intent(in) :: filename
    real(pv), allocatable, intent(out) :: matrix(:,:)

    ! Declare local variables
    integer :: i, j, ios, nrows, ncols, pos
    character(len=1000) :: line
    character(len=15) :: token
    real, allocatable :: temp_matrix(:,:)

    ! Open the file for reading
    open(unit=10, file=filename, status='old', action='read', &
         form='formatted', iostat=ios)
    if (ios /= 0) then
        print*, 'Error opening file ', filename
        stop
    end if

    ! Determine the number of rows and columns
    nrows = 0
    ncols = 0
    do
        read(10, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (nrows == 0) then
            ncols = 1
            do pos = 1, len_trim(line)
                if (line(pos:pos) == ',') ncols = ncols + 1
            end do
        end if
        nrows = nrows + 1
    end do

    ! Allocate the matrix
    allocate(temp_matrix(nrows, ncols))

    ! Rewind the file to the beginning
    rewind(10)

    ! Read the data into the matrix
    i = 1
    do
        read(10, '(A)', iostat=ios) line
        if (ios /= 0) exit
        pos = 1
        do j = 1, ncols
            call get_token(line, pos, token)
            read(token, *) temp_matrix(i, j)
        end do
        i = i + 1
    end do

    ! Close the file
    close(10)

    ! Return the matrix
    matrix = temp_matrix

    contains
        subroutine get_token(input_line, input_pos, output_token)
            ! This subroutine extracts a token from a line starting at position pos
            character(len=*), intent(in) :: input_line
            integer, intent(inout) :: input_pos
            character(len=15), intent(out) :: output_token
            integer :: start, end

            start = input_pos
            end = index(input_line(start:), ',') - 1
            if (end == -1) then
                end = len_trim(input_line) - start + 1
            end if
            output_token = input_line(start:start+end-1)
            input_pos = start + end + 1
        end subroutine get_token

end subroutine read_matrix



end module IO_Toolbox