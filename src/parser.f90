! This module read data file with different formats and analysis parameters.
! When performing analysis, this module is called within specific functional parts.
MODULE parser
    use precision
    use flags

    implicit none
    save
    character(len=80) :: path, file_name   ! path that contain files, or file names
    character(len=20) :: coption, doption
    character(len=20), allocatable :: fnames(:)
    integer :: fnumber, frame_interval
    character(len=30) :: cutoff_str, pbc_str
    logical :: static

    real(dp), allocatable :: cutoffs(:)

    real(dp) :: rdf_r, d2min_r
    integer :: pbcs(3)

contains

SUBROUTINE get_input_options()
    IMPLICIT NONE
!
    character(len=80) :: args
    integer(inp) :: n, i, k

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
        case ("-f")
            !data is a single file
            static = .true.
            print *, "Input data is a single file ... True"
            path=''
            call get_command_argument(i+1, args)
            file_name = trim(args)
            print *, "File name is ", trim(file_name)
        case ("-d")
            !data is a directory
            static = .false.
            print *, "Input data is a directory ... True"
            call get_command_argument(i+1, args)
            path = trim(args)
            call read_file_names()
        case('-os') ! offset
            if (static .eqv. .true.) stop '--- Should not have offset in static analysis mode.'
            call get_command_argument(i+1, args)
            read (args, *, iostat=k) frame_interval
            print *, 'Frame interval is ', frame_interval
        case ("-c")        !computing coption
                stop '--- Note: computing option is deprecated!'

        case('-rdf')
            flag_rdf = .true.
            call get_command_argument(i+1, args)
            read (args, *) rdf_r
            print *, 'RDF cutoff value is:', rdf_r
        case('-d2min')
            flag_nfd = .true.
            flag_d2min = .true.
            call get_command_argument(i+1, args)
            read (args, *) d2min_r
            print *, 'D2min cutoff value is:', d2min_r
        case('-bad')
            flag_nfd = .true.
            flag_bad = .true.
            call get_command_argument(i+1, args)
            if (verify('-', args) /= 0) then
                cutoffs = get_cutoff(args)
                print *, 'Cutoff values are:', cutoffs
            end if
        case('-ring')
            flag_nf = .true.
            flag_rstat = .true.
            call get_command_argument(i+1, args)
            if (verify('-', args) /= 0) then
                cutoffs = get_cutoff(args)
                print *, 'Cutoff values are:', cutoffs
            end if
        case('-tct')
            flag_nfd = .true.
            flag_tct = .true.
            call get_command_argument(i+1, args)
            if (verify('-', args) /= 0) then
                cutoffs = get_cutoff(args)
                print *, 'Cutoff values are:', cutoffs
            end if
        case('-poly')
            flag_nf = .true.
            flag_poly = .true.
            call get_command_argument(i+1, args)
            if (verify('-', args) /= 0) then
                cutoffs = get_cutoff(args)
                print *, 'Cutoff values are:', cutoffs
            end if
        case('-wa')
            flag_wa = .true.
            call get_command_argument(i+1, args)
            read (args, *) rdf_r
            print *, 'RDF cutoff value is:', rdf_r
        case('-ha')
            flag_nf = .true.
            flag_ha = .true.
            call get_command_argument(i+1, args)
            if (verify('-', args) /= 0) then
                cutoffs = get_cutoff(args)
                print *, 'Cutoff values are:', cutoffs
            end if
        case('-cluster')
            flag_cluster = .true.
            call get_command_argument(i+1, args)
            if (verify('-', args) /= 0) then
                cutoffs = get_cutoff(args)
                print *, 'Cutoff values are:', cutoffs
            end if
        case('-nc')
            flag_nc=.true.
            flag_nf = .true.
            call get_command_argument(i+1, args)
            if (verify('-', args) /= 0) then
                cutoffs = get_cutoff(args)
                print *, 'Cutoff values are:', cutoffs
            end if

        case ("--help", "-h")     !help document
            call help_msg()
        case ("--version", "-v")
            CALL version_msg()
        case default
            if(verify('-', args) == 0) then
                print *, '!ERROR: Input contains unknown variable: ', args
                stop
            end if
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
    print '(" Found ",i0," files in directory.")', nlines

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
function get_cutoff(str_in) result(r_list)
    implicit none
    ! IN:
    character(len=*), intent(in) :: str_in
    ! PRIV:
    integer :: p, k, i, n
    ! out:
    real(dp), allocatable :: r_list(:)

    if (allocated(r_list)) return

    if (index(str_in, ",") == 0) then
        allocate(r_list(1))
        read(str_in, *) r_list(1)
    else
        do i = 1, len_trim(str_in)
            if (str_in(i:i) == ',') then
                n=n + 1
            end if
        end do
        allocate(r_list(n+1))
        p = 1
        k = 0
        i = 1
        do while(index(str_in(p:), ",") /= 0)
            k = index(str_in(p:), ",") + k
            read(str_in(p:k-1),*) r_list(i)
            p = k+1
            i = i+1
        end do
        read(str_in(p:),*) r_list(i)
    end if
end function get_cutoff

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

print *, " Example usage: './polikit -f ../test/ga2o3_test.xyz -p 1 -poly 2.3'   (polyhedral analysis)"
print *, "                './polikit -f ../test/ga2o3_test.xyz -p 1 -bad 2.3'   (bond angle analysis)"
print *, "                './polikit -f ../test/ga2o3_test.xyz -p 1 -rdf 10'     (radial distribution)"
print *, "                './polikit -f ../test/ga2o3_test.xyz -p 1 -wa 5'      (Wendt-Abraham parameter)"
print *, "                './polikit -f ../test/ga2o3_test.xyz -p 1 -ha 2.3'    (Honeycutt-Anderson parameters)"
print *, "                './polikit -f ../test/ga2o3_test.xyz -p 1 -ring 2.3'    (ring statistics analysis)"
print *, "                './polikit -d ../test/test_dir/ 3 -p 1 -r 2.3 -c p'       (dynamic neighbor change)"
print *, " Variables:   "
print *, "   -f [string]     File name, supports .xyz, .dump, .data formats, incompatible with '-d' option."
print *, "                   Also be careful to match the data in each colume."
print *, "   -d [string]     Directory name and interval in dynamic analysis. Interval 0 will just perform  "
print *, "                   static analysis without comparing. Incompatible with '-f' option."
print *, "   -p [1 .or. 0]   Whether periodic boundary condition is activated."
print *, '   -os [int]       Offset for dynamic analysis, should be used with -d.'
print *, "   -[key] [float]  Analyzing options. 'poly' - polyhedral analysis; 'd2min' - non-affine displacement"
print *, "                 analysis; 'ring' - ring statistics analysis; 'bad' - bond angle distribution; 'rdf' - "
print *, "                 radial distribution function; 'wa' - WA parameter; 'ha' - HA parameter."
! print *, "   -o [string]        Atomic output options, won't dump atomic file if not set. 'n' - atomic coordination;"
! print *, "                        't' - tct results; 'p' - poly. neighbor; 'l' - linked state."

END SUBROUTINE

SUBROUTINE version_msg()

print *,    "Polikit - Atomistic Simulation Analysis Tool"
print *,    "    V0.4"
print *,    "Bug report: zjh239@foxmail.com"

END SUBROUTINE

END MODULE parser
