!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a test main program for input/out put and settings

program main
    use precision
!     use control
    use io
    use topo_analysis
!     use tct

    implicit none
    integer :: delength, i, fcounter = 0
    logical :: flag_poly = .false. , flag_tct = .false.
!     real(dp), dimension(natom) :: diff


    call initialize

!     read demand
    delength = len(trim(demand))
    do i = 1, delength
        select case (demand(i:i+1))
            case("p")
                flag_poly = .true.
            case("t")
                flag_tct = .true.
        end select
    end do

!     select case (trim(demand))
!         case("p")
!             flag_poly = .true.
!         case("t")
!             flag_tct = .true.
!         case("pt")
!
!     end select

    if (static .eqv. .true.) then
        if (flag_poly .eqv. .true.) then
            call poly_neighbor      !how does it know if we need static or dynamic results?
            call export
        end if
        if (flag_tct .eqv. .true.) then
!             call tct      !how does it know if we need static or dynamic results?
!             call export
                print *, "hi"
        end if
!         call export
    else       ! that static .eqv. .false.
        do i=1, fnumber
            fname = trim(fnames(i))
            call poly_neighbor
            call compare        !how does it know which set of data to compare? add atomic attributes.
        end do
!         print *,


    end if

contains

subroutine export_static
    implicit none
    ! block below has condition of: static == true, demand has 'p' included.
    if(static .eqv. .true. .and. flag_poly .eqv. .true.) then
        print *, "      CS       ES       FS"
        print *, count(ln==1), "   ",count(ln==2),"    ",count(ln==3)
        do i=1,natom
            print *, xyz(i,1), xyz(i,2), xyz(i,3), ln(i)
        end do
    end if


end subroutine

! this subroutine is to check if the interval meet, and whether the results should be compared and exported.
subroutine compare_export
    implicit none
    ! data must be in the module so that they won't lost
!     integer(4):: na
!     na = natom
    real(dp) :: diff(natom)
    type linkedl(na)
        integer(4), len :: na
        real(dp), dimension(na) :: list
        type(linkedl(*)), pointer :: next
    end type linkedl

    type(linkedl(natom)), pointer :: head, tail, tmp

    fcounter = fcounter + 1

    ! fisrt use linked list to store the computed results.
    if(.not. associated(head)) then
        allocate(head, stat=ierr)
        tail => head
        nullify(tail%next)
        tail%list = ln
    else
        allocate(tail%next, stat=ierr)
        tail => tail%next
        nullify(tail%next)
        tail%list = ln
    end if

    ! perform comparison if applicable.
    if (fcounter > interval) then
        ! compute difference between first and end of the linked list
        diff = tail%list - head%list

        ! release the first element of the linked list.
        tmp => head%next
        deallocate(head)
        allocate(head)
        head%list = tmp%list
        head%next => tmp%next

        ! export results
        call export
    end if

    if(static .eqv. .false. .and. flag_poly .eqv. .true.) then
!         print *, count(diff==1), "   ",count(diff==2),"    ",count(diff==3)

    end if
!     ref = ln
!     cframe
end subroutine

END PROGRAM main

