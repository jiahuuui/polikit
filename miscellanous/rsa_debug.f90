MODULE rings_simple
    USE precision
    USE data_input
    USE neighbor_finder, only: neigh_list
    USE stdlib_array
    IMPLICIT NONE
    integer, parameter :: maxlvl = 7   ! Max length, decided by ring size limit
    integer, parameter :: max_ringcap = 20
!     TYPE ring
!         integer :: l
!         integer :: element(20)
!         integer :: sorted(20)
!     END TYPE

    integer :: ring_cap, ring_size  ! Main ring list capacity and current size.
    real(dp) :: tstart, tcheck, tneilist, tfindring, tcheckrepi, tcheckpr, taddring
    integer :: crude_ring_num

!     type(ring) :: emptyring

CONTAINS

! Find the pos of the ring in the main ring list.
subroutine if_ring_exist(a_ring, mainringlist, pos, exist)
    implicit none
    ! IN:
    type(ring), intent(in) :: a_ring
    type(ring), allocatable, intent(in) :: mainringlist(:)
    ! OUT:
    integer, intent(out) :: pos
    logical, intent(out) :: exist
    ! Priv:
    integer :: n_rings, n_elem, i, j, elements(20), lastcol

    pos = -1
    exist = .true.

    n_rings = ring_size
    n_elem = a_ring%l
    elements(:n_elem) = a_ring%sorted(:n_elem)

    j = 1
    i = 1

    do while (i <= n_rings)
    if (j>1) then
                if (elements(j-1) /= mainringlist(i)%sorted(j-1)) then
                    ! print *, 'Different branch: ', elements(j-1), mainringlist(i)%sorted(j-1)
                    pos = i
                    exist = .false.
                    return
                end if
            end if
!             print *, i, j
        if (elements(j) < mainringlist(i)%sorted(j)) then
            ! print *, 'Smaller index: ', elements(j), mainringlist(i)%sorted(j)
            pos = i
            exist = .false.
            return
        elseif (elements(j) == mainringlist(i)%sorted(j)) then
            if (j == n_elem) then
                ! print *, 'Found the ring: '
                exist = .true.
                pos = i
                return
            end if
            lastcol = elements(j)
            j = j+1
!                 i = i-1
        else !(elements(j) > mainringlist(i)%element(j))
            if (i == n_rings) then
                ! print *, 'At end, cant find ring.'
                pos = i+1
                exist = .false.
                return
            end if
            i = i+1
        end if
    end do

    if (n_rings == 0) then
        exist = .false.
        pos = 1
        return
    end if

    print *, 'quiting ring exist checking', ring_size

end subroutine if_ring_exist

subroutine new_check_rp_d(ar, mainringlist, pos, goal_found)
    implicit none
    ! IN:
    type(ring), intent(in) :: ar
    type(ring), allocatable, intent(in) :: mainringlist(:)
    ! OUT:
    integer, intent(inout) :: pos
    logical, intent(inout) :: goal_found
    ! Private:
    integer :: low, high, middle, level, goal, row, n_elem, last_low, last_high

    n_elem = ar%l
    low = 1
    high = ring_size
    level = 1
    pos = 1

    goal_found = .true.

    do while(level <= n_elem)
        goal = ar%sorted(level)

        last_low = low
        last_high = high

        DO WHILE(low <= high)! .AND. pos == -1)
            print *, low, high, goal, level
            ! If item out of range, return
            if (goal < mainringlist(low)%sorted(level)) then
                pos = low
                goal_found = .false.
                return
            else if (goal > mainringlist(high)%sorted(level)) then
                pos = high+1
                goal_found = .false.
                return
            end if

            ! Now searching element middle.
            middle = (low + high)/2
            IF (goal == mainringlist(middle)%sorted(level)) THEN
                pos = middle
                exit
            ELSE IF (goal < mainringlist(middle)%sorted(level)) THEN
                high = middle-1
            ELSE
                low = middle+1
            END IF
        END DO
        ! Get the new range of the list.
        low = middle
        high = middle
        do row = pos, last_low, -1
            if(mainringlist(row)%sorted(level) == mainringlist(pos)%sorted(level)) then
                low = row
            else
                exit
            end if
        end do
        do row = pos, last_high, 1
            if(mainringlist(row)%sorted(level) == mainringlist(pos)%sorted(level)) then
                high = row
            ELSE
                exit
            end if
        end do

        level = level+1
    end do

    if (ring_size == 0) then
        goal_found = .false.
        pos = 1
        return
    end if
end subroutine new_check_rp_d

! Find the pos of the ring in the main ring list.
subroutine if_ring_exist_d(a_ring, mainringlist, pos, exist)
    implicit none
    ! IN:
    type(ring), intent(in) :: a_ring
    type(ring), allocatable, intent(in) :: mainringlist(:)
    ! OUT:
    integer, intent(out) :: pos
    logical, intent(out) :: exist
    ! Priv:
    integer :: n_rings, n_elem, i, j, elements(20), lastcol

    pos = -1
    exist = .true.

    n_rings = ring_size
    n_elem = a_ring%l
    elements(:n_elem) = a_ring%sorted(:n_elem)

    j = 1
    i = 1

    do while (i <= n_rings)
        print *, i, j
        if (j>1) then
            if (elements(j-1) /= mainringlist(i)%sorted(j-1)) then
                ! print *, 'Different branch: ', elements(j-1), mainringlist(i)%sorted(j-1)
                pos = i
                exist = .false.
                return
            end if
        end if

        if (elements(j) < mainringlist(i)%sorted(j)) then
            ! print *, 'Smaller index: ', elements(j), mainringlist(i)%sorted(j)
            pos = i
            exist = .false.
            return
        elseif (elements(j) == mainringlist(i)%sorted(j)) then
            if (j == n_elem) then
                ! print *, 'Found the ring: '
                exist = .true.
                pos = i
                return
            end if
            lastcol = elements(j)
            j = j+1
!                 i = i-1
        else !(elements(j) > mainringlist(i)%element(j))
            if (i == n_rings) then
                ! print *, 'At end, cant find ring.'
                pos = i+1
                exist = .false.
                return
            end if
            i = i+1
        end if
    end do

    if (n_rings == 0) then
        exist = .false.
        pos = 1
        return
    end if

    print *, 'quiting ring exist checking', ring_size

end subroutine if_ring_exist_d

FUNCTION checkShortCut_d(rr) RESULT(ifpr)
    IMPLICIT NONE
    ! IN:
    type(ring), intent(in) :: rr    ! Input ring
    ! OUT:
    logical :: ifpr
    ! PRIVATE:
    integer :: src1, src2, src3   ! check-node and its mid-node(even ring) or mid-nodes(odd ring).
    logical :: isodd
    integer, allocatable, dimension(:) :: elem  ! Stores two repeated ring list elements for a constant offset.
    integer, allocatable, dimension(:) :: head1, head2  ! Heads list of the wave.
    integer, allocatable, dimension(:) :: last1, last2  ! Heads list of last level.
    integer, allocatable, dimension(:) :: scndlast1, scndlast2
    integer, allocatable, dimension(:) :: tmp
    integer :: lvl, j, n, m, l, distance
    integer :: brlen, clen, mxlvl  ! branch length, current length

    call cpu_time(tstart)

    isodd = .false.
    if (mod(rr%l,2)/=0) isodd = .true.

    ifpr = .true.
    allocate(elem(rr%l*2))
    elem = [rr%element(:rr%l), rr%element(:rr%l)]   ! Avoid seg. fault.

    brlen = ceiling((rr%l+1)/2.)
    mxlvl = ceiling(brlen/2.)

    associate(n_neighbor => neigh_list%n_neighbor, neighbor => neigh_list%neighbors)
    do m = 1, brlen !-1

        distance = 0

        src1 = elem(m)
        src2 = elem(m+brlen-1)
        print *, 'The pair is', src1, src2
        if (isodd) src3 = elem(m+brlen)

        allocate(last1(1), source = 0)
        allocate(last2(1), source = 0)

        allocate(head1(1), source = src1)
        if (isodd) then
            allocate(head2(2), source = [src2, src3])
        else
            allocate(head2(1), source = src2)
        endif

        do lvl = 2, mxlvl
            ! update head1.
            call move_alloc(last1, scndlast1)
            call move_alloc(head1, last1)
            allocate(head1(1), source=0)

            n = 1
            do while (n <= size(last1))
                do j = 1, n_neighbor(last1(n))
                    ! if this atom is already in the wave, cycle
                    if (any(scndlast1==neighbor(last1(n), j))) cycle
                    if (any(last1==neighbor(last1(n), j))) cycle
                    if (any(head1==neighbor(last1(n), j))) cycle

                    ! if not cycled, push it to the list
                    if (head1(1)==0) then
                        head1(1) = neighbor(last1(n), j)
                    else
                        allocate(tmp(size(head1)+1))
                        tmp(:size(head1)) = head1
                        tmp(size(head1)+1) = neighbor(last1(n), j)
!                         deallocate(head1)
                        call move_alloc(tmp, head1)
                    end if
                end do

!                 if (n >= size(last1)) exit
                n = n+1
            end do
            clen = lvl*2 - 1
            print 122, clen, head1
            ! check if the branches meet.
            do n = 1, size(head1)
                if (any(head2 == head1(n))) then
                    ! print *, 'Short cut found, length is', 2*i-2, ', meet at ', tmp1(n)
                    ifpr = .false.
                    return
                end if
            end do

            if (clen == brlen) then
                if (m == brlen) then
                    ifpr = .true.
                    return
                else
                    cycle
                end if
            end if
            ! update head2.
            call move_alloc(last2, scndlast2)
            call move_alloc(head2, last2)
            allocate(head2(1), source=0)
            n = 1
            do while(n <= size(last2))
                do j = 1, n_neighbor(last2(n))
                    ! if this atom is already in the wave, cycle
                    if (any(scndlast2==neighbor(last2(n), j))) cycle
                    if (any(last2==neighbor(last2(n), j))) cycle
                    if (any(head2==neighbor(last2(n), j))) cycle

                    ! if not cycled, push it to the list
                    if (head2(1)==0) then
                        head2(1) = neighbor(last2(n), j)
                    else
                        allocate(tmp(size(head2)+1))
                        tmp(:size(head2)) = head2
                        tmp(size(head2)+1) = neighbor(last2(n), j)
!                         deallocate(head2)
                        call move_alloc(tmp, head2)
                    end if
                end do

!                 if (n >= size(last2)) exit
                n = n+1
            end do
            clen = lvl*2
            print 122, clen, head2
            ! check if the branches meet.
            do n = 1, size(head2)
                if (any(head1 == head2(n))) then
                    ifpr = .false.
                    return
                end if
            end do
            ! deallocate second last lists.
            if (allocated(scndlast1)) deallocate(scndlast1)
            if (allocated(scndlast2)) deallocate(scndlast2)
        end do

    if (allocated(head1)) deallocate(head1)
    if (allocated(head2)) deallocate(head2)

    if (allocated(last1)) deallocate(last1)
    if (allocated(last2)) deallocate(last2)

    end do

    end associate

    deallocate(elem)

    call cpu_time(tcheck)
    tcheckpr = tcheckpr + (tcheck - tstart)

122 format ( "("I0")", *(2x,I0))
END FUNCTION checkShortCut_d

END MODULE rings_simple
