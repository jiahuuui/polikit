!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a test main program for input/out put and settings

program main
    use io
    use nf
    use topo_analysis
    use tct
    implicit none
! computing options
    logical :: flag_poly = .false.  ! if poly_analysis is performed.
    logical :: flag_tct = .false.   ! if tct_analysis is performed.
    logical :: flag_dump = .false.  ! if to dump atomic properties.
    logical :: flag_nc = .false.    ! if neighbor change is performed.
! dump options
    logical :: dflag_cn = .false.
    logical :: dflag_poly = .false.
    logical :: dflag_pn = .false.
    logical :: dflag_tct = .false.

    integer :: clength, olength, i, fcounter = 0

    type intll(na)
        integer(4), len :: na
        integer(4), dimension(na) :: list
        type(intll(na)), pointer :: next
    end type intll

    type reall(na)
        integer(4), len :: na
        real(dp), dimension(na) :: list
        type(reall(na)), pointer :: next
    end type reall

!     logical ::
!     real(dp), dimension(natom) :: diff
    call initialize
!     get flags
    clength = len(trim(coption))
    do i = 1, clength
        select case (coption(i:i+1))
            case("p")
                flag_poly = .true.
            case("t")
                flag_tct = .true.
            case("d")
                flag_dump = .true.
            case('n')
                flag_nc = .true.
        end select
    end do

    olength = len(trim(doption))
    if(olength /= 0) flag_dump = .true.
    do i = 1, olength-1
        select case (coption(i:i+1))
            case("n")
                dflag_cn = .true.
            case("p")
                dflag_poly = .true.
            case("t")
                dflag_pn = .true.
            case("l")
                dflag_tct = .true.
        end select
    end do

! call certain subs according to flags

    if (static .eqv. .true.) then
        call static_analysis
    else       ! that static .eqv. .false.
        call dynamic_analysis
!         print *,
    end if
contains

subroutine static_analysis
    implicit none

    if (flag_poly .eqv. .true.) then
        call poly_neighbor      !how does it know if we need static or dynamic results?
    end if
    if (flag_tct .eqv. .true.) then
        print *, "Error: TCT analysis can not be static."
        stop
    end if
    ! block below has condition of: static == true, coption has 'p' included.
!     if(flag_poly .eqv. .true.) then

!         do i=1,natom
!             print *, xyz(i,1), xyz(i,2), xyz(i,3), ln(i)
!         end do
!     end if

    if(flag_dump .eqv. .true.) then
        ! how to dump
    end if
end subroutine

! this subroutine is to check if the interval meet, and whether the results should be compared and exported.
subroutine dynamic_analysis
    implicit none
    ! data must be module variables so that they won't lost
    integer(4) :: ntemp, lnc(4)
    integer(4), allocatable :: diff(:)
    
    type(intll(natom)), pointer :: head
    type(intll(natom)), pointer :: tail
    type(intll(natom)), pointer :: tmp

    nullify(head)
    nullify(tail)
    nullify(tmp)

    do fcounter=1, fnumber
        fname = trim(fnames(fcounter))
        print *, 'Performing analysis on ',fname
        call tct_calculate

        ! fisrt use linked list to store the computed results.

        if(.not. associated(head)) then
            allocate(head, stat=ierr)
            tail => head
            nullify(tail%next)
            tail%list = constrain(:,1)
            ntemp = natom
            allocate(diff(natom))
        else
            allocate(tail%next, stat=ierr)
            tail => tail%next
            nullify(tail%next)
            tail%list = constrain(:,1)
            if (size(tail%list) /= ntemp) then
                print *, '!Error: Unmatched atom number between frames!'
                stop
            end if
            ntemp = natom
        end if

        ! perform comparison if applicable, and deallocate the results after store them in linked list.
        if (fcounter > interval) then
            print *, 'At frame',fcounter,'/',fnumber,', comparison performed.'
!             ! compute difference between first and end of the linked list
            ! lnc = ln_change(natom, )
            diff = tail%list - head%list ! this line need update, the diff is a integer but the intersection of the list need to be compared actually.
!             ! release the first element of the linked list.
            tmp => head%next
            deallocate(head) ! if neccesary to do so?
            allocate(head)
            head%list = tmp%list
            head%next => tmp%next
            print *, count(diff==1), "   ",count(diff==2),"    ",count(diff==3)
!           call memcln
!             ! export results
        else
            print *, 'At frame',fcounter,'/',fnumber,', comparison is not performed.'
        end if

        call memcln

        print *, ' '
        print *, '========================================================================='
        print *, ' '
    end do

!     ref = ln
!     cframe
end subroutine

subroutine compare_nc
    implicit none


end subroutine

subroutine memcln
    implicit none

    if (allocated(poly))        deallocate(poly)
    if (allocated(ln))          deallocate(ln)

    if (allocated(n_neighbor))  deallocate(n_neighbor)
    if (allocated(neighbor))    deallocate(neighbor)

    if (allocated(delta))       deallocate(delta)
    if (allocated(constrain))   deallocate(constrain)

    if (allocated(ptype))       deallocate(ptype)
    if (allocated(xyz))         deallocate(xyz)

end subroutine

END PROGRAM main

