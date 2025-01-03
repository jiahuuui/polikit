!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! main program

program main
    use parser
    use neighbor_finder
    use poly_analysis
    use rdf
    use tct
    use bad
    use dynamic_data
    use rings_simple
    use data_input
    implicit none
! computing options
    LOGICAL :: flag_nf = .false.    ! neighbor finder flag.
    LOGICAL :: flag_nfd = .false.
    logical :: flag_poly = .false.  ! poly_analysis flag.
    logical :: flag_tct = .false.   ! tct_analysis flag.
    logical :: flag_dump = .false.  ! if to dump atomic properties.
    logical :: flag_nc = .false.    ! if neighbor change is performed.
    logical :: flag_bad = .false.   ! if bond angle distribution is analyzed.
    logical :: flag_rstat = .false. ! ring statistics analysis.
    logical :: flag_rdf = .false.   ! radial distribution function.

! dump options
    logical :: dflag_cn = .false.
    logical :: dflag_poly = .false.
    logical :: dflag_pn = .false.
    logical :: dflag_tct = .false.

    integer :: clength, olength
    character(len=1) :: c_slice, o_slice

    real :: start, finish
! we have integer type linked list and real type linked list, they are used in tct or
! momentary d2min analysis when analysis results in a certain range need to be stored
! and the head and tail elements are the ones compared.


    call cpu_time(start)
        ! put code to test here

    call get_input_options() ! get the command line input strings
!     get flags

    IF (verify('p',coption) == 0) THEN
        print *, "Perform polyhedral analysis ... True"
        flag_poly = .true.
        flag_nf = .true.
    END IF
    IF (verify('t',coption) == 0) THEN
        print *, "Perform topological contraint analysis ... True"
        flag_nfd = .true.
        flag_tct = .true.
!     ELSE IF (verify(coption,'d') == 0) THEN
!         flag_dump = .true.
    END IF
    IF (verify('n',coption) == 0) THEN
        flag_nf = .true.
        flag_nc = .true.
    END IF
    IF (verify('b',coption) == 0) THEN
        print *, "Perform bond angle analysis ... True"
        flag_nfd = .true.
        flag_bad = .true.
    ENDIF
    IF (verify('r', coption) == 0) THEN
        print *, "Perform ring statistics analysis ... True"
        flag_nf = .true.
        flag_rstat = .true.
    ENDIF
    IF (verify('g', coption) == 0) THEN
        print *, "Calculate radial distribution function ... True"
        flag_rdf = .true.
    ENDIF

    olength = len(trim(doption))
    if(olength /= 0) then
        flag_dump = .true.
        IF (verify('n', doption) == 0) THEN
            dflag_cn = .true.
        END IF
        IF (verify('p', doption) == 0) THEN
            dflag_poly = .true.
        END IF
        IF (verify('t', doption) == 0) THEN
            dflag_pn = .true.
        END IF
        IF (verify('l', doption) == 0) THEN
            dflag_tct = .true.
        END IF
    end if

    if (static .eqv. .true.) then
        print *, "Starting static analysis ..."
        call static_analysis()
    else
        print *, "Starting dynamic analysis ..."
        call dynamic()
!         print *, 'dynamic analysis report error'
!         print *, (interval < 0), (flag_nc .eqv. .true.)
!         stop
    end if
    call cpu_time(finish)
    print '("Wall time = ",f7.2," seconds.")', finish-start
contains

subroutine static_analysis()
    implicit none

    CALL get_data_from_file(file_name, path)

    if (flag_nf)  call find_neighbors() ! if (flag_nf)  call neighbor_finder_old

    if (flag_nfd) call find_neighbors_d()

    if (flag_rdf) call get_rdf()

    if (flag_poly) call poly_neighbor()

    if (flag_bad) call bond_angle_distribution()

    if (flag_rstat) call rsa_simple()
!         call poly_neighbor      !how does it know if we need static or dynamic results?
!     end if
    if (flag_tct) then
        print *, "Error: TCT analysis can not be static."
        stop
    end if

    ! block below has condition of: static == true, coption has 'p' included.

!     if(flag_dump) stop
        ! how to dump
!     end if
!
!     if (flag_rings) call rings
!     if (flag_cav) call cavity

end subroutine static_analysis

! this subroutine is for dynamic comparison, which means a constant interval between current
! frame and reference frame is given. It checks if the interval meet and the results should be compared and exported.

SUBROUTINE dynamic()
    IMPLICIT NONE
    integer :: fcounter

    print *, fnumber
    do fcounter = 1, fnumber

        file_name = trim(fnames(fcounter))
        print *, 'Performing dynamic analysis on ',file_name

        call static_analysis()

        call collect_data()

        call compare_data()

        call mem_clean()
    end do

END SUBROUTINE dynamic

subroutine collect_data()
    implicit none

    if (flag_nf)  call collect_neighbor() ! neighbor_finder
!     if (flag_nf)  call neighbor_finder_old

!     if (flag_nfd) call collect_neighbor_d()

!     if (flag_poly) call collect_poly()

!     if (flag_bad) call collect_tct()

end subroutine collect_data

subroutine collect_neighbor()
    implicit none

    if (.not. allocated(neigh_list_bf)) then
        allocate(neigh_list_bf(frame_interval,natom,11))
        neigh_list_bf = 0
        print *, 'Initializaing neighbor list data container ...'
    end if

    neigh_list_bf(:-2,:,:) = neigh_list_bf(2:,:,:)
    neigh_list_bf(-1,:,1) = neigh_list%n_neighbor
    neigh_list_bf(-1,:,2:) = neigh_list%neighbors

!     d_neighbors(:-2) = d_neighbors(2:)

!     if (.not. allocated(d_neighbors)) then
!         allocate(neighbor_list(atom_number = natom, capacity = 10) :: d_neighbors(frame_interval))
!         print *, 'Initializaing neighbor list data container ...'
!     end if
!
!     d_neighbors(:-2) = d_neighbors(2:)
!     print *, neigh_list%n_neighbor(100)
!     print *, d_neighbors(3)%n_neighbor(100)! = 0
!     d_neighbors(frame_interval)%neighbors = 0
!     d_neighbors(frame_interval) = neigh_list

end subroutine collect_neighbor

subroutine compare_data()
    implicit none

    if (flag_nf)  call compare_neighbor() ! neighbor_finder
!     if (flag_nf)  call neighbor_finder_old

!     if (flag_nfd) call compare_neighbor_d()

!     if (flag_poly) call compare_poly()

!     if (flag_bad) call compare_tct()

end subroutine compare_data

subroutine compare_neighbor()
    implicit none
    integer :: i, j
    integer :: cn_increase, cn_decrease, cn_changed
    cn_increase = 0
    cn_decrease = 0
    cn_changed = 0

    associate(this_frame => neigh_list_bf(-1,:,:), last_frame => neigh_list_bf(-2,:,:))
    do i = 1, natom
!         print *, this_frame(i,1), last_frame(i,1)
        if (this_frame(i,1) > last_frame(i,1)) then
            cn_increase = cn_increase + 1
        elseif (this_frame(i,1) < last_frame(i,1)) then
            cn_decrease = cn_decrease + 1
        elseif (this_frame(i,1) == last_frame(i,1)) then
            do j = 2, this_frame(i,1)
                if (this_frame(i,j) /= last_frame(i,j)) then
                    cn_changed = cn_changed + 1
                    exit
                end if
            end do
        end if
    end do
    print *, 'cn_increase cn_decrease cn_changed'
!     print *, cn_increase, cn_decrease, cn_changed
    end associate

    print *, cn_increase, cn_decrease, cn_changed
end subroutine compare_neighbor

subroutine mem_clean
    implicit none

!     if (allocated(poly))        deallocate(poly)
!     if (allocated(ln))          deallocate(ln)

!     if (allocated(n_neighbor))  deallocate(n_neighbor)
!     if (allocated(neighbor))    deallocate(neighbor)
!     if (allocated(delta))       deallocate(delta)

!     if (allocated(constrain))   deallocate(constrain)

    if (flag_nf .or. flag_nfd)    call clean_neighbor
    if (flag_rdf)   call clean_rdf

    call clean_xyz_data

end subroutine

END PROGRAM main

