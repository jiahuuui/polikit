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
