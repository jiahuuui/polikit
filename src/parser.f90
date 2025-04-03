! This module read data file with different formats and analysis parameters.
! When performing analysis, this module is called within specific functional parts.
MODULE parser
    use precision

    implicit none
    save
    character(len=80) :: path, file_name   ! path that contain files, or file names
    character(len=20) :: coption, doption
    character(len=20), allocatable :: fnames(:)
    integer :: fnumber, frame_interval
    character(len=30) :: cutoff_str, pbc_str
    logical :: static

    real(dp), allocatable :: cutoffs(:)
    integer :: pbcs(3)

contains

SUBROUTINE get_input_options()
    IMPLICIT NONE
!
    character(len=80) :: args
    integer(inp) :: n, i

    n = iargc()

    if(n == 0) then
        call help_msg()
    end if

    do i = 1, n
    call get_command_argument(i, args)
        select case (args)
        case ("-p")         !check if pbc is applied
            call get_command_argument(i+1, pbc_str)
            call get_pbc(pbc_str)
            print '(a,i2,i2,i2)', ' Periodic boundary conditions are: ', pbcs
!             read (args, *) pbc
!             if (pbc == 1) print *, "Periodic boundary conditions ... True"
        case ("-r")         !read cutoff distance
            call get_command_argument(i+1, cutoff_str)
            call get_cutoff(cutoff_str)
            print *, 'Cutoff values are:', cutoffs
!             read (args, *) cutoff_str
        case ("-f")         !data is a single file
            static = .true.
            print *, "Performing static analysis ... True"
            path=''
            call get_command_argument(i+1, args)
            file_name = trim(args)
            print *, "File name is ", trim(file_name)
        case ("-d")         !data is a directory
            static = .false.
            call get_command_argument(i+1, args)
            path = trim(args)
            call read_file_names()
            call get_command_argument(i+2, args)
            read (args, *) frame_interval
        case ("-c")        !computing coption
            call get_command_argument(i+1, coption)
            if(len_trim(coption)==0) then
                print *, 'Error: no argument specified for computing!'
                stop
            end if
        case ("-o")
            call get_command_argument(i+1, doption)   ! dump option
            if(len_trim(doption)==0) then
                print *, 'Error: no argument specified for output!'
            end if
        case ("--help", "-h")     !help document
            call help_msg()
        case ("--version", "-v")
            CALL version_msg()
        end select
    end do
    return
END SUBROUTINE

! read file names in the directory and sort the files.
SUBROUTINE read_file_names()
    IMPLICIT NONE
    integer :: i,j,k,m, tmp
    integer, allocatable, dimension(:) :: t
    character(len=40) :: tmp_name

    call system("ls "//trim(path)//"> fname.tmp")

    fnumber = get_n_lines("fname.tmp")

    allocate(fnames(fnumber))
    allocate(t(fnumber))

    open(unit=30, file="fname.tmp", status='old', iostat=ierr, iomsg=emsg)

    do i = 1, fnumber
        read (30,*, iostat=ierr) fnames(i)
    end do
    close(unit=30)
    call system("rm fname.tmp")

    ! This part is to sort the files.
    do i = 1, fnumber
        call get_digits(fnames(i), t(i))
        do j = 1, i-1
            if (t(i) < t(j)) then
                tmp=t(j)
                t(j) = t(i)
                t(i) = tmp

                tmp_name = fnames(j)
                fnames(j) = fnames(i)
                fnames(i) = tmp_name
            endif
        enddo
    end do

!     ! Pring the sorted files.
!     do i = 1, fnumber
!         print *, fnames(i)
!     end do

    deallocate(t)
END SUBROUTINE

! get the number of lines(frame number) from a file
FUNCTION get_n_lines(filename) RESULT(nlines)
    use precision
    implicit none
    integer :: nlines
    character(len=*) :: filename

    open(unit=21, file= filename, status='old', iostat=ierr, iomsg=emsg)
    nlines = 0
    do
        read(21, *, iostat=ierr, iomsg=emsg)
        if(ierr/=0) exit    ! print emsg if need to debug here
        nlines = nlines+1
    enddo
    close(unit=21)
    print '("Found ",i0," files in directory.")', nlines

END FUNCTION

! Get the number part of the file names.
subroutine get_digits(filename, num)
    implicit none
    ! IN:
    CHARACTER(LEN=*), INTENT(IN) :: filename
    ! OUT:
    INTEGER, INTENT(OUT) :: num
    ! Private:
    INTEGER :: i, start, end
    CHARACTER(LEN=256) :: num_str

    start = 0
    end = 0
    num_str = ''

    ! Find first digit in the filename
    DO i = 1, LEN_TRIM(filename)
        IF (filename(i:i) >= '0' .AND. filename(i:i) <= '9') THEN
            start = i
            EXIT
        END IF
    END DO

    ! Find last digit
    IF (start > 0) THEN
        DO i = start, LEN_TRIM(filename)
            IF (.NOT. (filename(i:i) >= '0' .AND. filename(i:i) <= '9')) THEN
                end = i - 1
                EXIT
            END IF
        END DO
        IF (end == 0) end = LEN_TRIM(filename)  ! If only numbers at the end
        num_str = filename(start:end)
    END IF

    ! Convert string to number
    READ (num_str, *, IOSTAT=ierr) num
    IF (ierr /= 0) print *, 'No number found!'  ! Default if no number found

end subroutine get_digits

! Read pair-wise cutoffs.
subroutine get_cutoff(str_in)
    implicit none
    ! IN:
    character(len=*), intent(in) :: str_in
    ! PRIV:
    integer :: p, k, i, n

    if (index(str_in, ",") == 0) then
        if (.not. allocated(cutoffs)) allocate(cutoffs(1))
        read(str_in, *) cutoffs(1)
    else
        do i = 1, len_trim(str_in)
            if (str_in(i:i) == ',') then
                n=n + 1
            end if
        end do
        if (.not. allocated(cutoffs)) allocate(cutoffs(n+1))
        p = 1
        k = 0
        i = 1
        do while(index(str_in(p:), ",") /= 0)
            k = index(str_in(p:), ",") + k
            read(str_in(p:k-1),*) cutoffs(i)
            p = k+1
            i = i+1
        end do
        read(str_in(p:),*) cutoffs(i)
    end if
end subroutine get_cutoff

! Read PBCs if more than 1 are given.
subroutine get_pbc(str_in)
    implicit none
    ! IN:
    character(len=*), intent(in) :: str_in
    ! PRIV:
    integer :: p, k, i

    if (index(str_in, ",") == 0) then
!         if (.not. allocated(pbcs)) allocate(pbcs(3))
        read(str_in, *) pbcs(1)
        pbcs = pbcs(1)
    else
!         if (.not. allocated(pbcs)) allocate(pbcs(3))
        p = 1
        k = 0
        i = 1
        do while(index(str_in(p:), ",") /= 0)
            k = index(str_in(p:), ",") + k
            read(str_in(p:k-1),*) pbcs(i)
            p = k+1
            i = i+1
        end do
        read(str_in(p:),*) pbcs(i)
    end if
end subroutine get_pbc

SUBROUTINE help_msg()

print *, " Example usage: './src/polikit -f ../test/ga2o3_test.xyz -p 1 -r 2.32 -c p'   (polyhedral analysis)"
print *, "                './src/polikit -f ../test/ga2o3_test.xyz -p 1 -r 2.32 -c b'   (bond angle analysis)"
print *, "                './src/polikit -f ../test/ga2o3_test.xyz -p 1 -r 10 -c g'     (radial distribution)"
print *, "                './src/polikit -f ../test/ga2o3_test.xyz -p 1 -r 5 -c w'      (Wendt-Abraham parameter)"
print *, "                './src/polikit -f ../test/ga2o3_test.xyz -p 1 -r 2.3 -c h'    (Honeycutt-Anderson parameters)"
print *, "                './src/polikit -f ../test/ga2o3_test.xyz -p 1 -r 2.3 -c r'    (ring statistics analysis)"
print *, "                './src/polikit -d ../test/test_dir/ 3 -p 1 -r 2.3 -c p'       (dynamic neighbor change)"
print *, " Variables:   "
print *, "   -f [string]        File name, supports .xyz, .dump, .data formats, incompatible with '-d' option."
print *, "                      Also be careful to match the data in each colume."
print *, "   -d [string] [int]  Directory name and interval in dynamic analysis. Interval 0 will just perform  "
print *, "                      static analysis without comparing. Incompatible with '-f' option."
print *, "   -r [float]         The cutoff value."
print *, "   -p [1 .or. 0]      Whether periodic boundary condition is activated."
print *, "   -c [string]        Analyzing options. 'p' - polyhedral analysis; 't' - tct analysis; 'r' - ring"
print *, "                      statistics analysis; 'b' - bond angle distribution; 'g' - radial distribution function."
print *, "                      'w' - WA parameter; 'h' - HA parameter."
! print *, "   -o [string]        Atomic output options, won't dump atomic file if not set. 'n' - atomic coordination;"
! print *, "                        't' - tct results; 'p' - poly. neighbor; 'l' - linked state."

END SUBROUTINE

SUBROUTINE version_msg()

print *,    "Polikit - Atomistic Simulation Analysis Tool"
print *,    "    V0.3"
print *,    "Bug report: zjh239@foxmail.com"

END SUBROUTINE

END MODULE parser
