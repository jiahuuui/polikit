! This module read data file with different formats and analysis parameters.
! When performing analysis, this module is called within specific functional parts.
MODULE parser
    use precision

    implicit none
    save
    character(len=80) :: path, file_name   ! path that contain files, or file names
    character(len=20) :: coption, doption
    character(len=20), allocatable :: fnames(:)
    integer :: fnumber, pbc, frame_interval
    real(dp) :: cutoff
    logical :: static

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
            call get_command_argument(i+1, args)
            read (args, *) pbc
            if (pbc == 1) print *, "Periodic boundary conditions ... True"
        case ("-r")         !read cutoff distance
            call get_command_argument(i+1, args)
            read (args, *) cutoff
        case ("-f")         !data is a single file
            static = .true.
            print *, "Performing static analysis ... True"
            path=''
            call get_command_argument(i+1, args)
            file_name = trim(args)
            print *, "File name is ", file_name
        case ("-d")         !data is a directory
            static = .false.
            call get_command_argument(i+1, args)
            path = trim(args)
            call read_file_names
            call get_command_argument(i+2, args)
            read (args, *) frame_interval
        case ("-c")        !computing coption
            call get_command_argument(i+1, coption)
            if(len(trim(coption))==0) then
                print *, 'Error: no argument specified for computing!'
            end if
        case ("-o")
            call get_command_argument(i+1, doption)   ! dump option
            if(len(trim(doption))==0) then
                print *, 'Error: no argument specified for output!'
            end if
        case ("--help")     !help document
            call help_msg()
        CASE ("-h")
            CALL help_msg()
        case ("--version")
            CALL version_msg()
        CASE ("-v")
            CALL version_msg()
        end select
    end do
    return
END SUBROUTINE

! read file names and tell it is a .xyz file or .dump file.
SUBROUTINE read_file_names()
    use precision
    IMPLICIT NONE
    !
    integer :: i,j,k,m, tmp
    integer, allocatable :: t(:)
    character(len=20) :: tmp_name
    call system("ls "//trim(path)//"> fname.tmp")

    fnumber = get_n_lines("fname.tmp")

    allocate(fnames(fnumber))
    allocate(t(fnumber))

    open(unit=30, file="fname.tmp", status='old', iostat=ierr, iomsg=emsg)

    do i = 1, fnumber
        read (30,*, iostat=ierr) fnames(i)
        print *, fnames(i)
    end do
    close(unit=30)
    call system("rm fname.tmp")

! this part is to distinguish .xyz and .dump
    m = len(trim(fnames(1)))

!     if(fnames(1)(m-3:m)=='.xyz') then
    if (verify('xyz', fnames(1)) == 0) then
        k = 4
    else if (verify('dump', fnames(1)) == 0) then
        k = 5
    end if

    do i = 1, fnumber
        m = len(trim(fnames(i)))
        read (fnames(i)(1:m-k),*, iostat = ierr) t(i)

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

    deallocate(t)
    return
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
    print '("Got ",i2," files from directory.")', nlines

END FUNCTION

SUBROUTINE help_msg

print *,    "Example usage: './polikit -f abc.xyz -p 1 -r 2.32 -c pt -o ntpl'       (for single file analysis)"
print *,    "               './polikit -d xyzfiles 10 -p 1 -r 2.32 -c pt -o ntpl'   (for multi file analysis)"
print *,    "   -f [string]        File name, which the data will be loaded from, incompatible with '-d' option."
print *,    "   -d [string] [int]  Directory name and interval in dynamic analysis. Files in the directory will "
print *,    "                      be handled one by one, incompatible with '-f' option."
print *,    "   -r [float]         The cutoff value."
print *,    "   -p [1 .or. 0]      Whether periodic boundary condition is considered."
print *,    "   -c [string]        Analyzing options. 'p' - polyhedral analysis; 't' - tct analysis; 'r' - ring"
print *,    "                      statistics analysis; 'b' - bond angle distribution."
print *,    "   -o [string]        Atomic output options, won't dump atomic file if not set. 'n' - atomic coordination;"
print *,    "                        't' - tct results; 'p' - poly. neighbor; 'l' - linked state."

END SUBROUTINE

SUBROUTINE version_msg

print *,    "Polikit - Atomistic Simulation Analysis Tool"
print *,    "    V0.3"
print *,    "Bug report: zjh263@126.com"

END SUBROUTINE

END MODULE parser
