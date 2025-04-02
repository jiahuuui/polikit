module bad  ! bond angle distribution
    use precision
    use neighbor_finder, only: neigh_list, delta
    use data_input, only: coord_data, natom, o_type
    contains

subroutine bond_angle_distribution()
  implicit none
  integer(inp) :: id, j, k, c, p
  integer :: oxo_id, xox_id
  integer :: oxo_cap, xox_cap
  integer :: neigh_1, neigh_2

  real(dp) :: angle = 0.0
  real(dp), dimension(:), allocatable :: oxo_angle, xox_angle, tmp

  print *, 'Performing bond angle distribution analysis...'

  oxo_cap = 100
  xox_cap = 100

    allocate(oxo_angle(oxo_cap))
    allocate(xox_angle(xox_cap))

    oxo_id = 1
    xox_id = 1

  !   open(27, file = 'oxo.bin', access = 'direct', form = 'unformatted', recl = 4)
    associate(n_neighbor => neigh_list%n_neighbor, ptype => coord_data%ptype)

    do id = 1, natom
      if (n_neighbor(id) < 2) cycle

      do j = 1, n_neighbor(id)-1
        neigh_1 = neigh_list%neighbors(id, j)
        if (ptype(neigh_1) == ptype(id)) cycle

        do k = j+1, n_neighbor(id)
          neigh_2 = neigh_list%neighbors(id, k)
          if (ptype(neigh_2) == ptype(id)) cycle

          angle= r2d*acos((delta(id,j,1)*delta(id,k,1)+delta(id,j,2)*delta(id,k,2)&
          +delta(id,j,3)*delta(id,k,3))/(norm2(delta(id,j,:))*norm2(delta(id,k,:))))

          if (ptype(id) /= o_type) then
            oxo_angle(oxo_id) = angle
            oxo_id = oxo_id + 1

            if (oxo_id == oxo_cap) then
              oxo_cap = oxo_cap*2
              call move_alloc(oxo_angle, tmp)
              allocate(oxo_angle(oxo_cap))
              oxo_angle(:oxo_id) = tmp
              deallocate(tmp)
            end if
          else
            xox_angle(xox_id) = angle
            xox_id = xox_id + 1

            if (xox_id == xox_cap) then
              xox_cap = xox_cap*2
              call move_alloc(xox_angle, tmp)
              allocate(xox_angle(xox_cap))
              xox_angle(:xox_id) = tmp
              deallocate(tmp)
            end if
          end if

        enddo
      enddo
    enddo
  end associate

  print *, "b1| O-X-O bond angle histgram ===================="
    call hist_of_bad(oxo_angle, oxo_id-1)
  print *, "b2| X-O-X bond angle histgram ===================="
    call hist_of_bad(xox_angle, xox_id-1)
  print *, "=================================================="

!   close(27, status='delete')
end subroutine bond_angle_distribution

! Sort the bond angles and draw histogram
subroutine hist_of_bad(array, cap)
  implicit none
  ! INOUT:
  real(dp), intent(inout) :: array(:)
  integer, intent(in) :: cap
  !
  integer :: bin(180), i, k
  real(dp) :: binedge, bin_center

  binedge = 1.0
  bin_center = 0.5

  call quicksort_nr(array(:cap))

  bin = 0

  k = 1
  i = 1
  do while(i <= cap)

    if (array(i) < binedge) then
      bin(k) = bin(k)+1

    else
      k = k+1
      binedge = binedge + 1.0
      i = i-1
    end if
    i = i+1
  end do

  do i = 1,180
    print "(f8.2, '    ', i0)", bin_center, bin(i)
    bin_center = bin_center + 1.0
  end do

end subroutine hist_of_bad

end module bad
