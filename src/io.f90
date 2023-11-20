! This module read data file with different formats and analysis parameters.
! When performing analysis, this module is called within specific functional parts.
MODULE io
    use precision
    implicit none
    save
    character(len=80) :: dpath
    character(len=20) :: fname, demand
    character(len=20), allocatable :: fnames(:)
    integer(4) :: natom, pbc, o_type, x_type, interval, fnumber
    integer(4), allocatable :: ptype(:)
    real(dp), allocatable :: xyz(:,:)
    real(dp) :: lx, ly, lz, cutoff
    logical :: static

contains

SUBROUTINE helpmessage

print *,    "Example usage: 'x.exe -f abc.xyz -p 1 -r 2.32'"
print *,    "   -f [filename]       Load data which is a single file, incompatible with '-d' option."
print *,    "   -d [dir] [inv]      Load files in the directory, incompatible with '-f' option."
print *,    "   -r [float]          The cutoff value."
print *,    "   -p [1 .or. 0]       Apply periodic boundary condition or not."
print *,    "   -o [pt]             String that controls what to compute, 'p' means polyhedral analysis, 't' means tct analysis."

END SUBROUTINE

SUBROUTINE initialize
    IMPLICIT NONE
    character*20 :: args
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
                    call getarg(i+1, args)
                    fname = trim(args)
                case ("-d")         !data is a directory
                    static = .false.
                    call getarg(i+1, args)
                    dpath = trim(args)
                    call readfnames
                    call getarg(i+2, args)
                    read (args, *) interval
                case ("-o")        !computing demand
                    call getarg(i+1, demand)
                case ("--help")     !help document
                    call helpmessage()
            end select
    end do

END SUBROUTINE

! read arguments include filename, pbc, cutoff, atomtype
SUBROUTINE get_xyz
    IMPLICIT NONE
    character(len=20) :: f
    integer :: slength
    call initialize

    f = trim(fname)
    slength = len(trim(fname))
    ! check the format
    if (fname(slength-2:slength)=='xyz') then
        call xyz_file
    else if (args(slength-3:slength)=='dump') then
        call dump_file
    end if

    return
END SUBROUTINE get_xyz

!this sub would be called if data is in .xyz format.
SUBROUTINE xyz_file
    IMPLICIT NONE
    integer(4) :: i
    character*10, allocatable :: typechar(:)

    open (unit=20, file=fname, status='old', iostat=ierr, iomsg=emsg)
    read (20,*, iostat=ierr) natom
    print *, natom," atoms read from file;"
    read (20,*, iostat=ierr)

    allocate(xyz(1:natom, 3), STAT=ierr, ERRMSG=emsg)
    allocate(typechar(1:natom), STAT=ierr, ERRMSG=emsg)

    do i=1, natom
        read (20,*, iostat=ierr) typechar(i), xyz(i,1), xyz(i,2), xyz(i,3)
    end do

    call type_convert(typechar)

    lx = MAXVAL(xyz(:,1)) - MINVAL(xyz(:,1))
    ly = MAXVAL(xyz(:,2)) - MINVAL(xyz(:,2))
    lz = MAXVAL(xyz(:,3)) - MINVAL(xyz(:,3))

    return

END SUBROUTINE

! this sub would be called if data is in .dump format(used for LAMMPS).
SUBROUTINE dump_file
    IMPLICIT NONE
    integer :: i

    open(unit=20, file=fname, status='old', iostat=ierr, iomsg=emsg)

    read (20,*, iostat=ierr) ! ITEM: BOX BOUNDS xx yy zz
    read (20,*, iostat=ierr) ! xlo xhi
    read (20,*, iostat=ierr) ! ylo yhi
    read (20,*, iostat=ierr) ! zlo zhi


END SUBROUTINE


! if type names are string in data file, we need to convert it to integer.
SUBROUTINE type_convert(charin)
    IMPLICIT NONE

    character*10, intent(in) :: charin(natom)
    character*10 :: typename(10)
    integer :: ntype, i, j

    allocate(ptype(1:natom), STAT=ierr, ERRMSG=emsg)
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
    do i = 1, ntype
        if(typename(i)=='O') then
            o_type = i
        else
            x_type = i
        end if
        print *, typename(i), i, count(ptype == i)  !atomintype(i)
    enddo

    return
END SUBROUTINE

SUBROUTINE readfnames
    IMPLICIT NONE
    integer :: i, k, j
    integer, allocatable :: t(:)


    call system("ls "//dpath//" > fname.tmp")
    fnumber = nlines("fname.tmp")
    allocate(fnames(fnumber))
    allocate(t(fnumber))

    open(unit=30, file="fname.tmp", status='old', iostat=ierr, iomsg=emsg)
    ! use linked list to store fnames

    do i = 1, fnumber
        read (30,*, iostat=ierr) fnames(i)
    end do
    close(unit=30)
    call system("rm fname.tmp")

    k = len(trim(fnames(1)))
    read(fnames(1)(1:k-4),*, iostat=ierr) t(1)  ! if read integer success, the file is .xyz format.
    if(ierr==0) then
        k = 4
    else
        k = 5
    end if

    do i = 1, fnumber
        j = trim(fnames(i))
        read (fnames(i)(1:j-4)), *, iostat = ierr) t(i)
    end do

    call bubble_sort(fnumber, t)

    deallocate(fnames)
    allocate(fnames)

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

FUNCTION nlines(filename)
    implicit none
    integer :: nlines
    character, intent(in) :: filename

    open(unit=21, file= filename, status='old', iostat=ierr, iomsg=emsg)
    nlines = 0
    do
        read(20, *, iostat=ierr)
        if(ierr/=0) exit
        nlines = nlines+1
    enddo
    close(unit=21)
    return

END FUNCTION

END MODULE io
