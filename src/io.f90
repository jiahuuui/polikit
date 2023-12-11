! This module read data file with different formats and analysis parameters.
! When performing analysis, this module is called within specific functional parts.
MODULE io
    use precision

    implicit none
    save
    character(len=80) :: dpath, fname
    character(len=20) :: coption, doption
    character(len=20), allocatable :: fnames(:)
    integer(4) :: natom, pbc, o_type, interval, fnumber
    integer(4), allocatable :: ptype(:)
    real(dp), allocatable :: xyz(:,:)
    real(dp) :: lx, ly, lz, cutoff
    logical :: static

contains

SUBROUTINE helpmessage

print *,    "Example usage: 'x.exe -f abc.xyz -p 1 -r 2.32 -c pt -o ntpl'       (for single file analysis)"
print *,    "               'x.exe -d xyzfiles 10 -p 1 -r 2.32 -c pt -o ntpl'   (for multi file analysis)"
print *,    "   -f [filename]       Load data which is a single file, incompatible with '-d' option."
print *,    "   -d [dir] [inv]      Load files in the directory, incompatible with '-f' option."
print *,    "   -r [float]          The cutoff value."
print *,    "   -p [1 .or. 0]       Apply periodic boundary condition or not."
print *,    "   -c [pt]             Compute option string. 'p' - polyhedral analysis; 't' - tct analysis."
print *,    "   -o [string]         Atomic output option string, won't dump atomic file if not set. 'n' - atomic coordination;"
print *,    "                       't' - tct results; 'p' - poly. neighbor; 'l' - linked state."

END SUBROUTINE

SUBROUTINE initialize
    IMPLICIT NONE
    character(len=80) :: args
    integer(4) :: n, i

    n = iargc()

    if(n == 0) then
        call helpmessage()
    end if

    do i = 1, n
        call getarg(i, args)
            select case (args)
                case ("-p")         !check if pbc is applied
                    call getarg(i+1, args)
                    read (args, *) pbc
                case ("-r")         !read cutoff distance
                    call getarg(i+1, args)
                    read (args, *) cutoff
                case ("-f")         !data is a single file
                    static = .true.
                    dpath=''
                    call getarg(i+1, args)
                    fname = trim(args)
                case ("-d")         !data is a directory
                    static = .false.
                    call getarg(i+1, args)
                    dpath = trim(args)
                    call readfnames
                    call getarg(i+2, args)
                    read (args, *) interval
                case ("-c")        !computing coption
                    call getarg(i+1, coption)
                    if(len(trim(coption))==0) then
                        print *, 'Error: no argument specified for computing!'
                    end if
                case ("-o")
                    call getarg(i+1, doption)   ! dump option
                    if(len(trim(doption))==0) then
                        print *, 'Error: no argument specified for output!'
                    end if
                case ("--help")     !help document
                    call helpmessage()
            end select
    end do
    return
END SUBROUTINE

! read arguments include filename, pbc, cutoff, atomtype
SUBROUTINE get_xyz
    IMPLICIT NONE
    integer :: slength
    print *, 'Entering get xyz subroutine, read fname is ', trim(dpath)//fname
    slength = len(trim(fname))
    ! check the format
    if (fname(slength-2:slength)=='xyz') then
        call xyz_file
    else if (fname(slength-3:slength)=='dump') then
        call dump_file
    end if

    print *, 'Leaving get xyz subroutine ...'
END SUBROUTINE get_xyz

!this sub would be called if data is in .xyz format.
SUBROUTINE xyz_file
    IMPLICIT NONE
    integer(4) :: i
    character(len=10), allocatable :: typechar(:)

    open (unit=20, file=trim(dpath)//fname, status='old', iostat=ierr, iomsg=emsg)
    if (ierr /=0) print *, emsg
    read (20,*, iostat=ierr) natom
    print *, natom," atoms read from ",fname
    read (20,*, iostat=ierr)

    allocate(xyz(natom, 3), STAT=ierr, ERRMSG=emsg)
    allocate(typechar(natom), STAT=ierr, ERRMSG=emsg)

    do i=1, natom
        read (20,*, iostat=ierr) typechar(i), xyz(i,1), xyz(i,2), xyz(i,3)
    end do

    call type_convert(typechar)

    lx = MAXVAL(xyz(:,1)) - MINVAL(xyz(:,1))
    ly = MAXVAL(xyz(:,2)) - MINVAL(xyz(:,2))
    lz = MAXVAL(xyz(:,3)) - MINVAL(xyz(:,3))

    deallocate(typechar)

END SUBROUTINE

! this sub would be called if data is in .dump format(used for LAMMPS).
SUBROUTINE dump_file
    IMPLICIT NONE
    integer :: i

    open(unit=20, file=trim(dpath)//fname, status='old', iostat=ierr, iomsg=emsg)

    read (20,*, iostat=ierr) ! ITEM: BOX BOUNDS xx yy zz
    read (20,*, iostat=ierr) ! xlo xhi
    read (20,*, iostat=ierr) ! ylo yhi
    read (20,*, iostat=ierr) ! zlo zhi


END SUBROUTINE


! if type names are string in data file, we need to convert it to integer. Oxygen atom is
! automatically distinguished by assuming oxygen takes the largest part.
SUBROUTINE type_convert(charin)
    IMPLICIT NONE

    character(len=10), intent(in) :: charin(natom)
    character(len=10) :: typename(10)
    integer :: ntype, i, j

    allocate(ptype(natom), STAT=ierr, ERRMSG=emsg)
    ptype = 0
    ntype = 1
    typename(1) = charin(1)

    do i = 1, natom     !compare char type with each existing type, add it to typename if not exist
        do j = 1, ntype
            if(charin(i) == typename(j)) then
                ptype(i) = j
                exit
            else if(j==ntype) then
                ntype = ntype + 1
                typename(ntype) = charin(i)
                ptype(i) = ntype
            end if
        end do

    end do
    print *, "typename    typeid    number"

    o_type = 1
    do i = 1, ntype
        print *, typename(i), i, count(ptype == i)  !atomintype(i)
        if (count(ptype==i) > count(ptype==o_type)) then
            o_type = i
        end if
    end do

END SUBROUTINE

SUBROUTINE readfnames
    IMPLICIT NONE
    integer :: i, k, j
    integer, allocatable :: t(:)

    call system("ls "//trim(dpath)//" > fname.tmp")
    fnumber = nlines("fname.tmp")

    allocate(fnames(fnumber))
    allocate(t(fnumber))

    open(unit=30, file="fname.tmp", status='old', iostat=ierr, iomsg=emsg)

    do i = 1, fnumber
        read (30,*, iostat=ierr) fnames(i)
    end do
    close(unit=30)
    call system("rm fname.tmp")

! this part is to distinguish .xyz and .dump
    k = len(trim(fnames(1)))
!     read(fnames(1)(1:k-4),*, iostat=ierr) t(1)  ! if read integer success, the file is .xyz format.
    if(fnames(1)(k-4:k)=='.xyz') then
        k = 4
    else if (fnames(1)(k-5:k)=='.dump') then
        k = 5
    end if

    do i = 1, fnumber
        j = len(trim(fnames(i)))
        read (fnames(i)(1:j-k),*, iostat = ierr) t(i)
    end do

    call bubble_sort(fnumber, t)

    deallocate(fnames)
    allocate(fnames(fnumber))

    if(k==4) then
        do i = 1, fnumber
            write (fnames(i),*) t(i),".xyz"
            print *, fnames(i)
        end do
    else if(k==5) then
        do i = 1, fnumber
            write (fnames(i),*) t(i),".dump"
            print *, fnames(i)
        end do
    end if

    return
END SUBROUTINE

! get the number of lines(atom number) from a file
FUNCTION nlines(filename)
    implicit none
    integer :: nlines
    character(len=*), intent(in) :: filename

    open(unit=21, file= filename, status='old', iostat=ierr, iomsg=emsg)
    nlines = 0
    do
        read(21, *, iostat=ierr, iomsg=emsg)
        if(ierr/=0) exit    ! print emsg if need to debug here
        nlines = nlines+1
    enddo
    close(unit=21)
    print *, 'Got',nlines,' files from directory.'

END FUNCTION

END MODULE io
