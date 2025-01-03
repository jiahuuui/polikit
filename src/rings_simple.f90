MODULE rings_simple
    USE precision
    USE data_input
    USE neighbor_finder, only: neigh_list
    USE stdlib_array
    IMPLICIT NONE
    TYPE ring
        integer :: l
        integer :: element(20)
    END TYPE

CONTAINS

SUBROUTINE rsa_simple()     ! Ring statistics analysis simple
    IMPLICIT NONE
    ! PRIVATE:
    integer :: atom     ! Center node index
    integer :: maxlvl   ! Max length, decided by ring size limit
    type(ring), allocatable, dimension(:) :: ringList
    integer(inp), allocatable, dimension(:,:) :: pathArray

    maxlvl = 4
    print *, maxlvl
    allocate(ringList(1))
    ringList(1)%l = 0
    ringList(1)%element = 0

    DO atom = 1, natom
        print *, atom
        CALL create_path_list(atom, maxlvl, pathArray)

        CALL find_rings(pathArray, ringList)

        ! CALL rm_not_pr(ringList)
    END DO

    call print_ringno(ringList%l)

END SUBROUTINE rsa_simple

! This subroutine creates the shortest paths list of a given center node(atom).
SUBROUTINE create_path_list(id_in, lvlim, pathArray)
    IMPLICIT NONE
    ! IN:
    integer(inp), intent(in) :: id_in, lvlim
    ! OUT:
    integer(inp), allocatable, dimension(:,:), intent(out) :: pathArray ! pathlist
    ! PRIVATE:
    integer(inp), allocatable, dimension(:,:) :: tmp
    integer(inp) :: i, row, k !
    integer(inp) :: lvl  !

    allocate(pathArray(1,1), source=id_in) !     lvl = 1

    associate(n_neighbor => neigh_list%n_neighbor, neighbor => neigh_list%neighbors)
    do lvl= 1, lvlim
        if (n_neighbor(id_in)==0) exit
        ! level is the current distance to the center node.
        ! Expand the path length for one step
        allocate(tmp(size(pathArray(:,1)), lvl+1), source=0)
        tmp(:,:lvl) = pathArray
        deallocate(pathArray)
        call move_alloc(tmp, pathArray)

        row = 1   ! loop over the second last column
        do  ! j = 1, size(pathArray(:,lvl))
            ! end point of odd ring
            if (lvl > 1) then
                if (any(pathArray(:, lvl-1) == pathArray(row,lvl))) then
                    go to 102
                else if (pathArray(row,lvl) == 0) then
                    go to 102
                end if
            end if
            k = 1
            do ! k = 1, n_neighbor(pathArray(row,lvl))
                ! skip if the atom already exist in former level.
                if (lvl > 1) then
                    if (any(pathArray(:, lvl-1) == neighbor(pathArray(row,lvl), k))) then
                        go to 101
                    end if
                end if

                if (pathArray(row,lvl+1)==0) then
                    pathArray(row,lvl+1) = neighbor(pathArray(row,lvl), k)
                else
                    ! Expand the path list for new path
                    allocate(tmp(size(pathArray(:,1))+1, size(pathArray(1,:))))
                    tmp(:size(pathArray(:,1)),:) = pathArray
                    deallocate(pathArray)  ! deallocated?
                    call move_alloc(tmp, pathArray)

                    row = row+1
                    pathArray(row:,:) = eoshift(pathArray(row:,:), shift=-1)
                    pathArray(row,:) = pathArray(row-1,:)
                    pathArray(row,lvl+1) = neighbor(pathArray(row,lvl), k)
                    ! because moved row row-1 to row, so here the neighbor list is still the same
                end if
101             if (k >= n_neighbor(pathArray(row,lvl))) exit
                k = k+1
            end do

102         if (row==size(pathArray(:,lvl))) exit
            row = row+1

        end do
!     call printa(pathArray)
    end do
    end associate
END SUBROUTINE create_path_list

! This subroutine finds all the possible rings around a center atom, given the
!   constructed shortest paths list. Push all the found rings to a data container.
SUBROUTINE find_rings(pathlist, mainringlist)
! natomring, ringatom, noprlist, numnopr, noprindex)
IMPLICIT NONE
    ! IN:
    integer(inp), allocatable, dimension(:,:), intent(in) :: pathlist
    ! inOUT:
    type(ring), allocatable, dimension(:), intent(inout) :: mainringlist
    ! PRIVATE:
    type(ring), allocatable, dimension(:) :: ringlist
    integer(inp) :: mxrow
    logical, allocatable, dimension(:,:) :: vis
    logical, allocatable, dimension(:) :: bpoint
    integer(inp), allocatable, dimension(:,:) :: vispl
    integer(inp), dimension(:), allocatable :: smlst
    integer(inp) :: rma, rmb, id, q, t(1), lvl, n, mxlvl, l, j, id2, id0, n1, n2
    integer(inp) :: row, row_2
    real(dp) :: start, end

    allocate(ringList(1))
    ringList(1)%l = 0
    ringList(1)%element = 0

    mxrow = size(pathlist(:,1))
    mxlvl = size(pathlist(1,:))
    allocate(bpoint(mxrow), source = .false.)

    id0 = pathlist(1,1)
    associate(n_neighbor => neigh_list%n_neighbor, neighbor => neigh_list%neighbors)
        allocate(smlst(n_neighbor(id0)*(n_neighbor(id0)-1)/2))
    end associate

    allocate(vis(mxrow,mxrow), source=.true.)
    forall (j = 1:mxrow)
        vis(j,j) = .false.
    end forall

    bpoint(1) = .true.
    rma = 1
    id  = pathlist(1,2)
    do row = 1, mxrow
        if (pathlist(row,2) /= id) then
            rmb = row-1
            vis(rma:rmb, rma:rmb) = .false.
            bpoint(row) = .true.

            id = pathlist(row,2)
            rma = row
        else if (row==mxrow) then
            vis(rma:row, rma:row) = .false.
        end if
    end do

    do lvl = 3, mxlvl
        id  = pathlist(1,lvl)
        do row = 1, mxrow
            if (pathlist(row, lvl) == 0) then
                vis(row,:) = .false.
                vis(:,row) = .false.
                cycle
            end if
            if (bpoint(row) .or. pathlist(row, lvl) /= id) then
                bpoint(row) = .true.
                id = pathlist(row,lvl)

                id2 = id
                do row_2= row, mxrow

                    if (bpoint(row_2) .or. pathlist(row_2, lvl) /= id2) then

                        if (vis(row, row_2)) then
                        id2 = pathlist(row_2,lvl)

                            ! check for odd ring
                            if   (pathlist(row_2, lvl-1)==pathlist(row, lvl)&
                            .and. pathlist(row_2, lvl)==pathlist(row, lvl-1)) then
            !                     print *, 'Odd ring, branch 1: ', pathlist(row,:lvl)
            !                     print *, '          branch 2: ', pathlist(row_2,:lvl)

                                call add_ring(pathlist(row,:lvl), pathlist(row_2,:lvl), ringlist)

                                vis(row_2,:) = .false.
                                vis(:,row_2) = .false.
                                vis(row,:) = .false.
                                vis(:,row) = .false.
                            end if

                            ! check for even rings
                            if (pathlist(row_2, lvl)==pathlist(row, lvl)) then
            !                     print *, 'Even ring, branch 1:', pathlist(row,:lvl)
            !                     print *, '           branch 2:', pathlist(row_2,:lvl)
                                call add_ring(pathlist(row,:lvl), pathlist(row_2,:lvl), ringlist)
                            end if
                        end if
                    end if
                end do
            end if
    !     call printl(vb)
        end do

        if (all(vis .eqv. .false.)) then
    !         print *, 'Visibility array full of FALSE at level ', lvl,', quit.'
            exit
        end if
    end do

    call rm_not_pr(ringList)

    call add_ringlist(mainringlist, ringList)

    deallocate(bpoint)
END SUBROUTINE find_rings

SUBROUTINE rm_not_pr(ringlist)
    IMPLICIT NONE
    ! INOUT:
    type(ring), intent(inout), allocatable, dimension(:) :: ringlist
    ! PRIVATE:
    integer(inp) :: rsize, nring
    logical, allocatable, dimension(:) :: pr
    type(ring), dimension(:), allocatable :: tmplist

    rsize = size(ringlist)
    allocate(pr(rsize))

    pr = .true.
    ! print *, lbound(pr), ubound(pr), rsize
    do nring = 1, rsize

        pr(nring) = checkShortCut(ringlist(nring))
!         print *, nring, pr(nring)
    end do

    if (any(pr .eqv. .false.)) then
        allocate(tmplist, source = ringlist(trueloc(pr)))
        deallocate(ringlist)
        deallocate(pr)
        call move_alloc(tmplist, ringlist)
    end if

END SUBROUTINE rm_not_pr

FUNCTION checkShortCut(rr) RESULT(ifpr)
    IMPLICIT NONE
    ! IN:
    type(ring), intent(in) :: rr    ! Input ring
    ! OUT:
    logical :: ifpr
    ! PRIVATE:
    integer :: src1, src2   ! Two source node.

    integer, allocatable, dimension(:) :: elem  ! Stores two repeated ring list elements for a constant offset.
    integer, allocatable, dimension(:) :: head1, head2  ! Heads list of the wave.
    integer, allocatable, dimension(:) :: last1, last2  ! Heads list of last level.
    integer, allocatable, dimension(:) :: tmp1, tmp2, wt
    integer :: i, j, k, n, m
    

    ifpr = .true.
    allocate(elem(rr%l*2))
    elem = [rr%element(:rr%l), rr%element(:rr%l)]

    k = rr%l/2 + 1

    associate(n_neighbor => neigh_list%n_neighbor, neighbor => neigh_list%neighbors)
ck: do m = 1, k-1

    src1 = elem(m)
    src2 = elem(m+k-1)

    if (allocated(head1)) deallocate(head1)
    if (allocated(head2)) deallocate(head2)

    allocate(head1(n_neighbor(src1)), source = neighbor(src1, :n_neighbor(src1)))
    allocate(head2(n_neighbor(src2)), source = neighbor(src2, :n_neighbor(src2)))

    if (allocated(last1)) deallocate(last1)
    if (allocated(last2)) deallocate(last2)

    allocate(last1(1), source = src1)
    allocate(last2(1), source = src2)

    do i  = 3, (k+1)/2
        ! print *, 'level is ', i, ', heads are ', elem(m), elem(m+k-1)
        allocate(tmp1(1), source = 0)
        allocate(tmp2(1), source = 0)

        n = 1
        do
            do j = 1, n_neighbor(head1(n))
                ! if this atom is already in the wave, cycle
                if (any(last1==neighbor(head1(n), j))) cycle
                if (any(head1==neighbor(head1(n), j))) cycle
                if (any(tmp1==neighbor(head1(n), j))) cycle

                ! if not cycled, push it to the list
                if (tmp1(1)==0) then
                    tmp1(1) = neighbor(head1(n), j)
                else
                    allocate(wt(size(tmp1)+1))
                    wt(:size(tmp1)) = tmp1
                    wt(size(tmp1)+1) = neighbor(head1(n), j)
                    call move_alloc(wt, tmp1)
                end if
            end do

            if (n >= size(head1)) exit
            n = n+1
        end do

        do n = 1, size(tmp1)
            if (any(head2 == tmp1(n)) .and. i<(k+1)/2) then
                ! print *, 'Short cut found, length is', 2*i-2, ', meet at ', tmp1(n)
                ifpr = .false.
                exit ck
            end if
        end do

        n = 1
        do
            do j = 1, n_neighbor(head2(n))
                ! if this atom is already in the wave, cycle
                if (any(last2==neighbor(head2(n), j))) cycle
                if (any(head2==neighbor(head2(n), j))) cycle
                if (any(tmp2==neighbor(head2(n), j))) cycle

                ! if not cycled, push it to the list
                if (tmp2(1)==0) then
                    tmp2(1) = neighbor(head2(n), j)
                else
                    allocate(wt(size(tmp2)+1))
                    wt(:size(tmp2)) = tmp2
                    wt(size(tmp2)+1) = neighbor(head2(n), j)
                    call move_alloc(wt, tmp2)
                end if
            end do
            if (n >= size(head2)) exit
            n = n+1
        end do

        do n = 1, size(tmp1)
            if (any(tmp2 == tmp1(n)) .and. i<(k+1)/2) then
                ! print *, 'Short cut found, length is', 2*i-1, ', meet at ', tmp1(n)
                exit ck
            end if
        end do

        deallocate(last1)
        deallocate(last2)

        call move_alloc(head1, last1)
        call move_alloc(head2, last2)

        call move_alloc(tmp1, head1)
        call move_alloc(tmp2, head2)
    end do
    end do ck

    end associate

    deallocate(elem)
END FUNCTION checkShortCut

! This subroutine adds the short ring list to the total ring list.
subroutine add_ringlist(mainringlist, tmpringlist)
    implicit none
    ! IN:
    type(ring), intent(in) :: tmpringlist(:)
    ! INOUT:
    type(ring), allocatable, intent(inout) :: mainringlist(:)
    ! private:
    type(ring), allocatable :: tmp(:)
    integer :: l1, l2

    l1 = size(mainringlist)
    l2 = size(tmpringlist)

    allocate(tmp(l1+l2))
    tmp(:l1) = mainringlist
    tmp(l1+1:l1+l2) = tmpringlist
    deallocate(mainringlist)
    call move_alloc(tmp, mainringlist)

end subroutine add_ringlist

! This subroutine adds a new ring type element to the ring list. The inputs are two
!   integer type lists, means the two branches of a ring.
SUBROUTINE add_ring(branch1, branch2, ringlist)
    IMPLICIT NONE
    ! IN:
    integer(inp), intent(in) :: branch1(:), branch2(:)
    ! INOUT:
    type(ring), allocatable, intent(inout) :: ringlist(:)
    ! PRIVATE:
    type(ring), allocatable :: tmp(:)
    type(ring) :: ar
    logical :: isodd
    integer(inp) :: k, i, t

    k = size(branch1)
    ar%l = 0
    ar%element = 0

    if (branch1(k) == branch2(k)) then
        isodd=.false.
    else
        isodd=.true.
    end if

    if (isodd) then
        ar%l = 2*k-3
        ar%element(:ar%l) = [branch1(:k-2), branch2(k:2:-1)]
    else
        ar%l = 2*k-2
        ar%element(:ar%l) = [branch1(:k-1), branch2(k:2:-1)]
    end if

    k = size(ringlist)
    if (k==1 .and. ringlist(1)%l==0) then
        ringlist(1) = ar

    else
        allocate(tmp(k+1))
        tmp(:k) = ringlist
        deallocate(ringlist)
        tmp(k+1) = ar
        call move_alloc(tmp, ringlist)

    end if

END SUBROUTINE add_ring

subroutine print_ringno(ring_l)
  implicit none
  ! IN:
  integer, dimension(:), intent(in) :: ring_l
  !
  integer :: maxn, i
  integer, allocatable :: rank(:), amount(:)

!     ringList%l
!   associate(ring_l => ringList%l)
    print *, size(ring_l), maxval(ring_l)
    maxn = maxval(ring_l)

    if (maxn > 20) maxn = 20

    allocate(rank(0:maxn), amount(0:maxn))

    rank = [(i, i=0, maxn, 2)]
    do i = 0, maxn, 2
      amount(i) = count(ring_l == i)
    end do

    print *, "Size (of cation) = ", rank
    print *, "          Amount = ", amount(0:maxn:2)

!   end associate
end subroutine print_ringno

pure function randomness(ringa)
    implicit none
    ! Input:
    type(ring), intent(in) :: ringa
    ! Output:
    real(dp) :: randomness

    integer :: i, id, t, r_count

    t = coord_data%ptype(ringa%element(1))

    do i = 2, ringa%l
        id = ringa%element(i)
        if (t /= coord_data%ptype(id)) r_count = r_count + 1
        t = coord_data%ptype(id)
    end do
    if (t /= coord_data%ptype(1)) r_count = r_count + 1

    randomness = r_count/ringa%l

end function

END MODULE rings_simple
