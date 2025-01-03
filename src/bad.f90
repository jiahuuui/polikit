module bad  ! bond angle distribution
    use precision
    use neighbor_finder, only: neigh_list, delta
    use data_input, only: coord_data, natom, o_type
    contains

subroutine bond_angle_distribution()
  implicit none
  integer(inp) :: id, j, k, c, p
  integer :: oxo_number, xox_number, oxo_id, xox_id
  real(dp) :: angle = 0.0
  real(dp), allocatable :: oxo_angle(:), xox_angle(:)

  print *, 'Performing bond angle distribution analysis...'

  oxo_number = 0 ! O-X-O bond number
  xox_number = 0 ! X-O-X bond number

  associate(n_neighbor => neigh_list%n_neighbor, ptype => coord_data%ptype)

    do id = 1, natom
      if (n_neighbor(id) > 1) then
        if (ptype(id) /= o_type) then
          oxo_number = oxo_number + (n_neighbor(id)*(n_neighbor(id)-1))/2
        else
          xox_number = xox_number + (n_neighbor(id)*(n_neighbor(id)-1))/2
        end if
      end if
    end do

    print *, 'O-X-O bond angle number is ', oxo_number
    print *, 'X-O-X bond angle number is ', xox_number

    allocate(oxo_angle(oxo_number))
    allocate(xox_angle(xox_number))

    oxo_id = 1
    xox_id = 1

  !   open(27, file = 'oxo.bin', access = 'direct', form = 'unformatted', recl = 4)

    do id = 1, natom
        do j = 1, n_neighbor(id)-1
          do k = j+1, n_neighbor(id)
            angle= r2d*acos((delta(id,j,1)*delta(id,k,1)+delta(id,j,2)*delta(id,k,2)&
            &+delta(id,j,3)*delta(id,k,3))/(norm2(delta(id,j,:))*norm2(delta(id,k,:))))
            if (ptype(id) /= o_type) then
              oxo_angle(oxo_id) = angle
  !             write(27 , rec = c) angle
              oxo_id = oxo_id + 1
            else
              xox_angle(xox_id) = angle
              xox_id = xox_id + 1
            end if
          enddo
        enddo
    enddo
  end associate
  print *, "=========== O-X-O bond angle histgram ============"
    call hist_of_bad(oxo_angle)
  print *, "=========== O-X-O bond angle histgram ============"
    call hist_of_bad(xox_angle)
  print *, "=================================================="

!   close(27, status='delete')
end subroutine bond_angle_distribution

subroutine hist_of_bad(array)
  implicit none
  real(dp), intent(inout) :: array(:)
  integer(inp) :: bin(180), i, k
  real(dp) :: binedge, bin_center

  binedge = 1.0
  bin_center = 0.5

  call quicksort_nr(array)

  bin = 0

  k = 1

  i = 1
  do
    if (array(i) < binedge) then
      bin(k) = bin(k)+1
    else
      k = k+1
      binedge = binedge + 1.0
      i = i-1
    end if
    i = i+1
    if (i > size(array)) exit
  end do

  do i = 1,180
    print *, bin_center, bin(i)
    bin_center = bin_center + 1.0
  end do

end subroutine hist_of_bad

end module bad
