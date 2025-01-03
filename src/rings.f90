MODULE rings
    use precision
    use data_input
    use neighbor_finder, only: neigh_list
    use stdlib_array
    implicit none

    type intlinlis
        integer :: id
        type(intlinlis), pointer :: last
        type(intlinlis), pointer :: next
    end type intlinlis

    TYPE ring
        integer, len :: cap
        integer :: l
        type(intlinlis), dimension(cap) :: element
    END TYPE

    integer :: criteria

!     INTERFACE PrintArray
!     MODULE PROCEDURE :: printa, printl, printr
!     END INTERFACE

contains

! SUBROUTINE testdata
!     IMPLICIT NONE

!     deallocate(neighbor)
!     deallocate(n_neighbor)

!     allocate(neighbor(17,4))
!     allocate(n_neighbor(17))

!     natom = 18
!     n_neighbor = [4,2,4,1,1,2,3,3,2,2,3,4,3,3,4,2,2,2]
!     neighbor = transpose(reshape(&
!     [2,7,11,15, &
!      1,3,0,0, &
!      2,4,5,15, &
!      3,0,0,0, &
!      3,0,0,0, &
!      8,18,0,0, &
!      1,8,10,0,&
!      6,7,9,0,&
!      8,12,0,0,&
!      7,11,0,0,&
!      1,10,12,0,&
!      9,11,13,14,&
!      12,18,0,0,&
!      17,15,12,0,&
!      1,3,14,16,&
!      15,17,0,0,&
!      14,16,0,0,&
!      6,13,0,0], [4,18]))



! END SUBROUTINE testdata

SUBROUTINE rsa  ! ring statistics analysis
    IMPLICIT NONE
    integer(inp) :: n, m
    integer(inp) :: maxlv     ! The length limit of paths, maximum ring size = 2*maxlv-2
    integer(inp), allocatable, dimension(:,:) :: path_array, prindex, noprindex
    logical, allocatable, dimension(:,:) :: vis, vb
    type(ring), allocatable, dimension(:) :: prlist, noprlist
    real :: ctime, ltime
    integer(inp), allocatable, dimension(:) :: numpr, numnopr

    maxlv = 9
    criteria = 3
    ! 1 - Flawed primitive ring criteria (debug usage);
    ! 2 - Shortest path criteria;
    ! 3 - Corrected primitive ring criteria (this one should always be used, unless clearer criteria is found).

    allocate(prindex(natom, 800), source = 0)
    allocate(numpr(natom), source = 0)
    allocate(prlist(1))
    prlist(1)%l = 0
    prlist(1)%element = 0

    allocate(noprindex(natom, 800), source = 0)
    allocate(numnopr(natom), source = 0)
    allocate(noprlist(1))
    noprlist(1)%l = 0
    noprlist(1)%element = 0

!     call testdata
!     ringlist(1)%l = 6
!     ringlist(1)%element(:6) = [1,7,8,9,12,11]
    print *, 'Criteria is ', criteria
    call CPU_TIME(ctime)
!     do n = 1, 1
    do n = 16983,16983
!         print *, numpr(n), numnopr(n)
        if (n_neighbor(n) == 0) cycle
        if (mod(n,1000) == 0) then
            call CPU_TIME(ltime)
            print '("Last 1000 atoms takes ",f7.3," seconds.")', ltime-ctime
            ctime = ltime
        end if

        call get_pl(n, maxlv, path_array)
!         call printa(path_array)
        m = size(path_array(:,1))
        call load_rings(m, vis, vb, path_array, n, prlist, numpr, prindex, noprlist, numnopr, noprindex)
!         call printl(vis)
        call check_rings(m, path_array, vis, vb, prlist, numpr, prindex, noprlist, numnopr, noprindex)
        call check_shortcut(prlist(9))
    end do
    call printr(prlist)
    call printr(noprlist)

    print *, count(prlist%l== 4), count(prlist%l== 6), count(prlist%l== 8), count(prlist%l== 10),&
    count(prlist%l== 12), count(prlist%l== 14), count(prlist%l== 16), count(prlist%l== 18)

    deallocate(prlist)
    deallocate(prindex)
    deallocate(numpr)

    deallocate(noprlist)
    deallocate(noprindex)
    deallocate(numnopr)
!     print *, 'Working on it...'
END SUBROUTINE rsa

! Get the path list start from a node. The principle is, from each node on the path to the first node is a
! shortest path. Should be careful about nodes that appear already exist.
SUBROUTINE create_path_list(id_in, lvlim, pathArray)
    IMPLICIT NONE
    ! IN:
    integer, intent(in) :: id_in, lvlim
    ! OUT:
    integer, allocatable, dimension(:,:), intent(out) :: pathArray ! pathlist
    ! PRIVATE:
    integer, allocatable, dimension(:,:) :: tmp
    integer :: i, row, k !
    integer :: lvl  !

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

! Check for existing rings for a node, and modify the visibility array
! according to the existing rings.
SUBROUTINE load_rings(m, vis, vb, pathlist, node,&
                    ringlist, natomring, atomring,&
                    noprlist, numnopr, noprindex)
    IMPLICIT NONE
    integer(inp), intent(in) :: m, node
    integer(inp) :: j, k
    logical, allocatable, dimension(:,:), intent(out) :: vis, vb
    integer(inp), allocatable, dimension(:,:), intent(in) :: pathlist, atomring, noprindex
    type(ring), allocatable, intent(in) :: ringlist(:), noprlist(:)
    integer(inp), allocatable :: branch1(:), branch2(:)
    integer(inp), allocatable, dimension(:) :: natomring, numnopr

    allocate(vis(m,m), source=.true.)
    allocate(vb(m,m), source=.true.)
    forall (j = 1:m)
        vis(j,j) = .false.
    end forall

    do j = 1, natomring(node)
        k = atomring(node, j)
        call splitring(ringlist(k), node, branch1, branch2)
        if (criteria == 1) then
            call mod_pr(branch1, branch2, vis, pathlist)
        else if (criteria == 2) then
            call mod_sp(branch1, branch2, vis, pathlist)
        else if (criteria == 3) then
            call norep(branch1, branch2, vis, pathlist)
            call mod_pr(branch1, branch2, vb, pathlist)
        end if
    end do

    do j = 1, numnopr(node)
        k = noprindex(node, j)
        call splitring(noprlist(k), node, branch1, branch2)
        if (criteria == 1) then
            call mod_pr(branch1, branch2, vis, pathlist)
        else if (criteria == 2) then
            call mod_sp(branch1, branch2, vis, pathlist)
        else if (criteria == 3) then
            call norep(branch1, branch2, vis, pathlist)
            call mod_pr(branch1, branch2, vb, pathlist)
        end if
    end do
END SUBROUTINE load_rings

! This subroutine find all the rings based on the path lists and VA mask.
SUBROUTINE check_rings(m, pathlist, vis, vb,&
                    ringlist, natomring, ringatom,&
                    noprlist, numnopr, noprindex)
    IMPLICIT NONE
    logical, allocatable, dimension(:,:), intent(inout) :: vis, vb
    logical, allocatable, dimension(:) :: bpoint

    integer(inp), allocatable, dimension(:,:), intent(in) :: pathlist
    integer(inp), allocatable, dimension(:,:), intent(inout) :: ringatom, noprindex
    integer(inp), allocatable, dimension(:), intent(inout) :: natomring, numnopr
    integer(inp), intent(in) :: m
    integer(inp), allocatable, dimension(:,:) :: vispl
    integer(inp), dimension(:), allocatable :: smlst
    integer(inp) :: rma, rmb, p, id, q, t(1), lvl, n, mxlvl, l, j, o, id2, id0, n1, n2

    type(ring), allocatable, dimension(:), intent(inout) :: ringlist, noprlist

    mxlvl = size(pathlist(1,:))
    allocate(bpoint(m), source = .false.)

    id0 = pathlist(1,1)
    allocate(smlst(n_neighbor(id0)*(n_neighbor(id0)-1)/2))

    bpoint(1) = .true.
    rma = 1
    id  = pathlist(1,2)
    do p = 1, m
        if (pathlist(p,2) /= id) then
            rmb = p-1
            vis(rma:rmb, rma:rmb) = .false.
            bpoint(p) = .true.

            id = pathlist(p,2)
            rma = p
        else if (p==m) then
            vis(rma:p, rma:p) = .false.
        end if
    end do

    do lvl = 3, size(pathlist(1,:))
    id  = pathlist(1,lvl)
    do p = 1, m
        if (pathlist(p, lvl) == 0) then
            vis(p,:) = .false.
            vis(:,p) = .false.
            cycle
        end if
        if (bpoint(p) .or. pathlist(p, lvl) /= id) then
            bpoint(p) = .true.
            id = pathlist(p,lvl)

            id2 = id
            do o= p, m

            if (bpoint(o) .or. pathlist(o, lvl) /= id2) then

            if (vis(p, o)) then
            id2 = pathlist(o,lvl)

                ! check for odd ring
                if   (pathlist(o, lvl-1)==pathlist(p, lvl)&
                .and. pathlist(o, lvl)==pathlist(p, lvl-1)) then
!                     print *, 'Odd ring, branch 1: ', pathlist(p,:lvl)
!                     print *, '          branch 2: ', pathlist(o,:lvl)
                    if (criteria == 1) then
                        call mod_pr(pathlist(p,:lvl), pathlist(o,:lvl), vis, pathlist)
                        call addring(pathlist(p,:lvl), pathlist(o,:lvl), ringlist, natomring, ringatom)
                    else if (criteria == 2) then
                        call mod_sp(pathlist(p,:lvl), pathlist(o,:lvl), vis, pathlist)
                        call addring(pathlist(p,:lvl), pathlist(o,:lvl), ringlist, natomring, ringatom)
                    else if (criteria == 3) then
                        call norep(pathlist(p,:lvl), pathlist(o,:lvl), vis, pathlist)
                        if (vb(p,o)) then
                            call addring(pathlist(p,:lvl), pathlist(o,:lvl), ringlist, natomring, ringatom)
                        else
                            call addring(pathlist(p,:lvl), pathlist(o,:lvl), noprlist, numnopr, noprindex)
                        end if
                        call mod_pr(pathlist(p,:lvl), pathlist(o,:lvl), vb, pathlist)
                    end if

                    vis(o,:) = .false.
                    vis(:,o) = .false.
                    vis(p,:) = .false.
                    vis(:,p) = .false.
                end if

                ! check for even rings
                if (pathlist(o, lvl)==pathlist(p, lvl)) then
!                     print *, 'Even ring, branch 1:', pathlist(p,:lvl)
!                     print *, '           branch 2:', pathlist(o,:lvl)
                    if (criteria == 1) then
                        call mod_pr(pathlist(p,:lvl), pathlist(o,:lvl), vis, pathlist)
                        call addring(pathlist(p,:lvl), pathlist(o,:lvl), ringlist, natomring, ringatom)
                    else if (criteria == 2) then
                        call mod_sp(pathlist(p,:lvl), pathlist(o,:lvl), vis, pathlist)
                        call addring(pathlist(p,:lvl), pathlist(o,:lvl), ringlist, natomring, ringatom)
                    else if (criteria == 3) then
                        call norep(pathlist(p,:lvl), pathlist(o,:lvl), vis, pathlist)
                        if (vb(p,o)) then
                            call addring(pathlist(p,:lvl), pathlist(o,:lvl), ringlist, natomring, ringatom)
                        else
                            call addring(pathlist(p,:lvl), pathlist(o,:lvl), noprlist, numnopr, noprindex)
                        end if
                        call mod_pr(pathlist(p,:lvl), pathlist(o,:lvl), vb, pathlist)
                    end if
                end if

            end if
            end if
            end do

        end if
!     call printl(vb)
    end do

    do j = 1, n_neighbor(id0)
        if (branch1(2)==neighbor(id0, j)) n1 = j
    end do
    do j = 1, n_neighbor(id0)
        if (branch2(2)==neighbor(id0, j)) n2 = j
    end do
    if (n1>n2) then
        j = n2
        n2 = n1
        n1 = j
    end if
    j = (n1-1)*n+(-n1+1)*n1/2+(n2-n1)
    if (smlst(j) /= 0 .and. size < smlst(j)) smlst(j) = size

    if (all(vis .eqv. .false.)) then
!         print *, 'Visibility array full of FALSE at level ', lvl,', quit.'
        exit
    end if
    end do
    deallocate(bpoint)
END SUBROUTINE check_rings

! After finding a node in a ring, the two branches from the node to the prime mid-node are
! formed. They are used in modifying the visibility array.
SUBROUTINE splitring(rin, node, branch1, branch2)
    IMPLICIT NONE
    type(ring), intent(in) :: rin
    integer(inp), intent(in) :: node
    integer(inp) :: t(1), ab, k, j
    integer(inp), allocatable :: tmp(:)
    integer(inp), allocatable, intent(out) :: branch1(:), branch2(:)
    allocate(tmp(3*rin%l))
    tmp = [rin%element(:rin%l), rin%element(:rin%l), rin%element(:rin%l)]

    ab = rin%l-2*rin%l/2
    if (ab == 0) then
        k = 1+rin%l/2
    else if (ab == 1) then
        k = 2+rin%l/2
    end if

    allocate(branch1(k))
    allocate(branch2(k))

    t = findloc(rin%element, node)
    j = t(1)+rin%l

    branch1 = tmp(j:j+k-1)
    branch2 = tmp(j:j-k+1:-1)

END SUBROUTINE splitring

! Make a ring from two branches and push them to the ringlist.
SUBROUTINE addring(branch1, branch2, ringlist, natomring, atomring)
    IMPLICIT NONE
    integer(inp), intent(in) :: branch1(:), branch2(:)
    type(ring), allocatable, intent(inout) :: ringlist(:)
    type(ring), allocatable :: tmp(:)
    type(ring) :: ar
    logical :: isodd
    integer(inp) :: k, i, t
    integer(inp), allocatable, intent(inout) :: natomring(:), atomring(:,:)

    k = size(branch1)

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

        do i = 1, ar%l
        t = ar%element(i)
        natomring(t) = natomring(t)+1
        atomring(t, natomring(t)) = k
        end do
    else
        allocate(tmp(k+1))
        tmp(:k) = ringlist
        deallocate(ringlist)
        tmp(k+1) = ar
        call move_alloc(tmp, ringlist)

        do i = 1, ar%l
        t = ar%element(i)
        if (natomring(t)>=size(atomring(1,:))) then
            print *, 'Warning: ring list about to full.'
        end if
        natomring(t) = natomring(t)+1
        atomring(t, natomring(t)) = k+1
        end do
    end if

!     print *, 'Adding ring'
    ! Put the ring id to the ring list of each member.

END SUBROUTINE addring

! Avoid repitition by modify the visibility array.
SUBROUTINE norep(branch1, branch2, vis, pathlist)
    IMPLICIT NONE
    integer(inp), allocatable, intent(in) :: pathlist(:,:)
    integer(inp), intent(in) :: branch1(:), branch2(:)
    logical, intent(inout) :: vis(:,:)
    integer(inp) :: k, m, a, i, j
    logical :: isodd
    logical, allocatable :: mask1(:), mask2(:), tmp(:)

    k = size(branch1)
    m = size(vis(:,1))
    allocate(mask1(m), source = .true.)
    allocate(mask2(m), source = .true.)
    allocate(tmp(m), source = .true.)

    if (branch1(k) == branch2(k)) then
        isodd=.false.
    else
        isodd=.true.
    end if

    if (.not. isodd) then
!         do i = 2, k-1
!             mask1 = .true.
!             mask2 = .true.
!
!             do j = 2,i
!                 tmp = pathlist(:,j)==branch1(j)
!                 mask1 = mask1 .and. tmp
!             end do
!             a = k+1-i
!     !         branch2(:a)
!             do j = 2, a
!                 tmp = pathlist(:,j) == branch2(j)
!                 mask2 = mask2 .and. tmp
!             end do
!             vis(trueloc(mask1), trueloc(mask2)) = .false.
!             vis(trueloc(mask2), trueloc(mask1)) = .false.
!         end do
!     else !if (isodd .eqv. .false.) then
    ! even ring
        mask1 = .true.
        mask2 = .true.

        do i = 2, k
            tmp = pathlist(:,i) == branch1(i)
            mask1 = mask1 .and. tmp

            tmp = pathlist(:,i) == branch2(i)
            mask2 = mask2 .and. tmp
        end do

        vis(trueloc(mask1), trueloc(mask2)) = .false.
        vis(trueloc(mask2), trueloc(mask1)) = .false.
    end if

    deallocate(mask1)
    deallocate(mask2)
    deallocate(tmp)
END SUBROUTINE norep

! Modify the visibility array according to the shortest path definition.
SUBROUTINE mod_sp(branch1, branch2, vis, pathlist)
    IMPLICIT NONE
    integer(inp), allocatable, intent(in) :: pathlist(:,:)
    integer(inp), intent(in) :: branch1(:), branch2(:)
    logical, intent(inout) :: vis(:,:)
    integer(inp) :: k, m, a, i, j
    logical :: isodd
    logical, allocatable :: mask1(:), mask2(:), tmp(:)

    k = size(branch1)
    m = size(vis(:,1))
    allocate(mask1(m), source = .true.)
    allocate(mask2(m), source = .true.)
    allocate(tmp(m), source = .true.)

    if (branch1(k) == branch2(k)) then
        isodd=.false.
    else
        isodd=.true.
    end if

    mask1 = .true.
    mask2 = .true.

    tmp = pathlist(:,2)==branch1(2)
    mask1 = mask1 .and. tmp

    tmp = pathlist(:,2) == branch2(2)
    mask2 = mask2 .and. tmp

    vis(trueloc(mask1), trueloc(mask2)) = .false.
    vis(trueloc(mask2), trueloc(mask1)) = .false.

    deallocate(mask1)
    deallocate(mask2)
    deallocate(tmp)
END SUBROUTINE mod_sp

! Modify the visibility array according to the primitive ring definition.
SUBROUTINE mod_pr(branch1, branch2, vis, pathlist)
    IMPLICIT NONE
    integer(inp), allocatable, intent(in) :: pathlist(:,:)
    integer(inp), intent(in) :: branch1(:), branch2(:)
    logical, intent(inout) :: vis(:,:)
    integer(inp) :: k, m, a, i, j
    logical :: isodd
    logical, allocatable :: mask1(:), mask2(:), tmp(:)

    k = size(branch1)
    m = size(vis(:,1))
    allocate(mask1(m), source = .true.)
    allocate(mask2(m), source = .true.)
    allocate(tmp(m), source = .true.)

    if (branch1(k) == branch2(k)) then
        isodd=.false.
    else
        isodd=.true.
    end if

    if (isodd) then
        do i = 2, k-1
            mask1 = .true.
            mask2 = .true.

            do j = 2,i
                tmp = pathlist(:,j)==branch1(j)
                mask1 = mask1 .and. tmp
            end do
            a = k+1-i
    !         branch2(:a)
            do j = 2, a
                tmp = pathlist(:,j) == branch2(j)
                mask2 = mask2 .and. tmp
            end do
            vis(trueloc(mask1), trueloc(mask2)) = .false.
            vis(trueloc(mask2), trueloc(mask1)) = .false.
        end do
    else !if (isodd .eqv. .false.) then
    ! even ring
        do i = 2, k
            mask1 = .true.
            mask2 = .true.
    !         branch1(:i)
            do j = 2,i
                tmp = pathlist(:,j) == branch1(j)
                mask1 = mask1 .and. tmp
            end do
            a = k+2-i
    !         branch2(:a)
            do j = 2, a
                tmp = pathlist(:,j) == branch2(j)
                mask2 = mask2 .and. tmp
            end do
            vis(trueloc(mask1), trueloc(mask2)) = .false.
            vis(trueloc(mask2), trueloc(mask1)) = .false.
        end do
    end if

    deallocate(mask1)
    deallocate(mask2)
    deallocate(tmp)
END SUBROUTINE mod_pr

! Form a new ring from two existing rings, gives the two branches.
! SUBROUTINE twoto1(r1, r2, branch1, branch2)
!     IMPLICIT NONE
!     type(ring), intent(in) :: r1, r2
!     integer(inp), allocatable, dimension(:), intent(out) :: branch1, branch2
!
!
!
! END SUBROUTINE twoto1

!
! ! Print integer 2-D array.
! SUBROUTINE printa(array)
!     IMPLICIT NONE
!     integer(inp), allocatable, intent(in) :: array(:,:)
!     integer :: i
!     do i = 1, size(array(:,1))
!         print *, array(i,:)
!     end do
!     print *, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
! END SUBROUTINE printa
!
! ! Print logical 2-D array.
! SUBROUTINE printl(array)
!     IMPLICIT NONE
!     logical, allocatable, intent(in) :: array(:,:)
!     integer :: i
!     do i = 1, size(array(:,1))
!         print *, array(i,:)
!     end do
!     print *, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
! END SUBROUTINE printl
!
! ! Print rings.
! SUBROUTINE printr(ringlist)
!     IMPLICIT NONE
!     type(ring), allocatable, intent(in) :: ringlist(:)
!     integer :: i, l
!     do i = 1, size(ringlist(:))
!         l = ringlist(i)%l
!         print *, i, '(ID) ', ringlist(i)%l, '(RSIZE) ', ringlist(i)%element(:l)
!     end do
!     print *, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
! END SUBROUTINE printr

! Check specifically if a ring is primitive.
! Only for even rings.
FUNCTION check_shortcut(rr) RESULT(ifpr)
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
END FUNCTION check_shortcut

END MODULE rings
