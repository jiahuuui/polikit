! From a center atom, this function gives every atom in the system a level value
! that denotes the shortest distance from the atom to the center atom.
SUBROUTINE dijkstra(center, lvllim, n_neighbor, neighbors, lvldist, natom)
  IMPLICIT NONE

  ! Declare inputs and variables
  INTEGER, INTENT(IN) :: center, lvllim, natom
  INTEGER, INTENT(IN) :: n_neighbor(:), neighbors(:,:)
  INTEGER, INTENT(OUT) :: lvldist(:)

  ! Local variables
  INTEGER :: t, quebgn, quend, nodcrt
  INTEGER :: lnkscrt, nodprb, lvlprb
  INTEGER, ALLOCATABLE :: queue(:)
  INTEGER :: dat

  ! Initialize variables
  ALLOCATE(queue(natom))
  DO t = 1, natom
    lvldist(t) = lvllim + 2
  END DO
  lvldist(center) = 0
  queue(1) = center
  quebgn = 1
  quend = 1
  dat = 0

  ! Main loop for Dijkstra's algorithm
  DO WHILE (quebgn <= quend)
    nodcrt = queue(quebgn)
    lvlprb = lvldist(nodcrt) + 1
    quebgn = quebgn + 1

    ! Iterate through the neighbors of the current node
    DO lnkscrt = 1, n_neighbor(nodcrt)
      nodprb = neighbors(nodcrt, lnkscrt)

      IF (lvldist(nodprb) > lvlprb) THEN
        lvldist(nodprb) = lvlprb

        ! Add to the queue if level is below the requirement
        IF (lvlprb < lvllim) THEN
          quend = quend + 1
          queue(quend) = nodprb
          dat = dat + 1
        END IF
      END IF
    END DO
  END DO

  ! Deallocate queue
  DEALLOCATE(queue)

END SUBROUTINE dijkstra

! function pair_search(atom1, atom2, ref_dis) result(has_sc)
!     implicit none
!     integer, intent(in) :: atom1, atom2 ! A node and its mid-node.
!     integer, intent(in) :: ref_dis      ! Reference distance.
!     logical :: has_sc
!
!     if (allocated(head1)) deallocate(head1)
!     if (allocated(head2)) deallocate(head2)
!
!     allocate(head1(n_neighbor(atom1)), source = neighbor(atom1, :n_neighbor(atom1)))
!     allocate(head2(n_neighbor(atom2)), source = neighbor(atom2, :n_neighbor(atom2)))
!
!     if (allocated(last1)) deallocate(last1)
!     if (allocated(last2)) deallocate(last2)
!
!     allocate(last1(1), source = atom1)
!     allocate(last2(1), source = atom2)
!
!     do i  = 3, (k+1)/2
!         ! print *, 'level is ', i, ', heads are ', elem(m), elem(m+k-1)
!         allocate(tmp1(1), source = 0)
!         allocate(tmp2(1), source = 0)
!
!         n = 1
!         do
!             do j = 1, n_neighbor(head1(n))
!                 ! if this atom is already in the wave, cycle
!                 if (any(last1==neighbor(head1(n), j))) cycle
!                 if (any(head1==neighbor(head1(n), j))) cycle
!                 if (any(tmp1==neighbor(head1(n), j))) cycle
!
!                 ! if not cycled, push it to the list
!                 if (tmp1(1)==0) then
!                     tmp1(1) = neighbor(head1(n), j)
!                 else
!                     allocate(wt(size(tmp1)+1))
!                     wt(:size(tmp1)) = tmp1
!                     wt(size(tmp1)+1) = neighbor(head1(n), j)
!                     call move_alloc(wt, tmp1)
!                 end if
!             end do
!
!             if (n >= size(head1)) exit
!             n = n+1
!         end do
!
!         do n = 1, size(tmp1)
!             if (any(head2 == tmp1(n)) .and. i<(k+1)/2) then
!                 ! print *, 'Short cut found, length is', 2*i-2, ', meet at ', tmp1(n)
!                 ifpr = .false.
!                 exit ck
!             end if
!         end do
!
!         n = 1
!         do
!             do j = 1, n_neighbor(head2(n))
!                 ! if this atom is already in the wave, cycle
!                 if (any(last2==neighbor(head2(n), j))) cycle
!                 if (any(head2==neighbor(head2(n), j))) cycle
!                 if (any(tmp2==neighbor(head2(n), j))) cycle
!
!                 ! if not cycled, push it to the list
!                 if (tmp2(1)==0) then
!                     tmp2(1) = neighbor(head2(n), j)
!                 else
!                     allocate(wt(size(tmp2)+1))
!                     wt(:size(tmp2)) = tmp2
!                     wt(size(tmp2)+1) = neighbor(head2(n), j)
!                     call move_alloc(wt, tmp2)
!                 end if
!             end do
!             if (n >= size(head2)) exit
!             n = n+1
!         end do
!
!         do n = 1, size(tmp1)
!             if (any(tmp2 == tmp1(n)) .and. i<(k+1)/2) then
!                 ! print *, 'Short cut found, length is', 2*i-1, ', meet at ', tmp1(n)
!                 exit ck
!             end if
!         end do
!
!         deallocate(last1)
!         deallocate(last2)
!
!         call move_alloc(head1, last1)
!         call move_alloc(head2, last2)
!
!         call move_alloc(tmp1, head1)
!         call move_alloc(tmp2, head2)
!     end do
!
! end function pair_search
