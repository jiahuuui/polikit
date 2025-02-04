module neighbor_finder
  use precision
  use data_types
  use parser, only: cutoff, pbc
  use data_input, only: coord_data, natom

  implicit none
!   save
  ! neighbor list
  type neighbor_list(atom_number, capacity)
      integer, len :: atom_number, capacity
      integer, dimension(atom_number) :: n_neighbor = 0
      integer, dimension(atom_number, capacity) :: neighbors = 0
  end type

  ! bin type define
  type bins(capacity)
      integer, len :: capacity   !capacity of the bin
      integer :: n = 0
      integer :: x_pbc = 0, y_pbc = 0, z_pbc = 0
      integer, dimension(capacity) :: ids = 0  !id of the atoms in bin
  end type

  real(dp), allocatable :: delta(:,:,:) ! neighbor(natom, 10), n_neighbor(natom)

  type(neighbor_list(atom_number = :, capacity = :)), allocatable :: neigh_list

  contains

! Divide the simulation box into bins, and put atoms into their corresponding bins.
SUBROUTINE create_bins(rCut, cellbins, xbin_max, ybin_max, zbin_max)
  IMPLICIT NONE
  ! in:
  real(dp), intent(in) :: rCut
  ! inOUT:
  type(bins(capacity = :)), intent(inout), allocatable, dimension(:,:,:) :: cellbins
  integer, intent(inout) :: xbin_max, ybin_max, zbin_max
  ! PRIVATE:
  real(dp), allocatable, dimension(:,:) :: realxyz
  REAL(dp) :: x_min, y_min, z_min
  INTEGER :: xbin, ybin, zbin, atom, bincap

  bincap = rCut**3 * 2
  print *, 'Deviding simulation box into bins ...'
  print *, "Capacity of bins is set to: ", bincap

  associate(xyz0 => coord_data%coord, lx => coord_data%lx, ly => coord_data%ly, lz => coord_data%lz)

    x_min = coord_data%xmin
    y_min = coord_data%ymin
    z_min = coord_data%zmin

    allocate(realxyz(natom, 3), STAT=ierr, ERRMSG=emsg)
    realxyz(:,1) = xyz0(:,1) - x_min
    realxyz(:,2) = xyz0(:,2) - y_min
    realxyz(:,3) = xyz0(:,3) - z_min

    xbin_max = CEILING(lx/rCut) - 1
    ybin_max = CEILING(ly/rCut) - 1
    zbin_max = CEILING(lz/rCut) - 1

    print *, "Number of bins on dimension X: ", xbin_max
    print *, "Number of bins on dimension Y: ", ybin_max
    print *, "Number of bins on dimension Z: ", zbin_max

  end associate

  allocate(bins(bincap) :: cellbins(0:xbin_max+1,0:ybin_max+1,0:zbin_max+1))
  cellbins%n = 0
  cellbins%x_pbc = 0
  cellbins%y_pbc = 0
  cellbins%z_pbc = 0

  do atom = 1, natom

    xbin = CEILING(realxyz(atom, 1)/rCut)
    ybin = CEILING(realxyz(atom, 2)/rCut)
    zbin = CEILING(realxyz(atom, 3)/rCut)

    if (xbin > xbin_max) xbin = xbin_max
    if (ybin > ybin_max) ybin = ybin_max
    if (zbin > zbin_max) zbin = zbin_max
    if (xbin == 0) xbin = 1
    if (ybin == 0) ybin = 1
    if (zbin == 0) zbin = 1

    cellbins(xbin,ybin,zbin)%n = cellbins(xbin,ybin,zbin)%n + 1
    cellbins(xbin,ybin,zbin)%ids(cellbins(xbin,ybin,zbin)%n) = atom

  end do

  if (pbc == 1) then
    do xbin = 0, xbin_max + 1
    do ybin = 0, ybin_max + 1
    do zbin = 0, zbin_max + 1

      associate(x_pbc => cellbins(xbin, ybin, zbin)%x_pbc, &
                y_pbc => cellbins(xbin, ybin, zbin)%y_pbc, &
                z_pbc => cellbins(xbin, ybin, zbin)%z_pbc)

        if (xbin == 0) x_pbc = 1
        if (xbin == xbin_max + 1) x_pbc = -1
        if (ybin == 0) y_pbc = 1
        if (ybin == ybin_max + 1) y_pbc = -1
        if (zbin == 0) z_pbc = 1
        if (zbin == zbin_max + 1) z_pbc = -1
        cellbins(xbin, ybin, zbin)%n = &
        cellbins(xbin + x_pbc*xbin_max, ybin + y_pbc*ybin_max, zbin + z_pbc*zbin_max)%n
        cellbins(xbin, ybin, zbin)%ids = &
        cellbins(xbin + x_pbc*xbin_max, ybin + y_pbc*ybin_max, zbin + z_pbc*zbin_max)%ids
      end associate
    end do
    end do
    end do
  endif

  print *, 'Leaving bins construction subroutine ...'

END SUBROUTINE create_bins

SUBROUTINE find_neighbors()
  IMPLICIT NONE
  !
  type(bins(capacity = :)), allocatable, dimension(:,:,:) :: cellbins
  integer :: xbin_max, ybin_max, zbin_max
  integer :: xbin, ybin, zbin, atom, atom2, i, p, q, o
  integer :: id, checkid
  real(dp) :: d, r

  r = cutoff**2

  call create_bins(cutoff, cellbins, xbin_max, ybin_max, zbin_max)
  print *, 'Constructing neighbor list from bins data ...'

  allocate(neighbor_list(atom_number = natom, capacity = 10) :: neigh_list)
  neigh_list%n_neighbor = 0
  neigh_list%neighbors = 0

  associate(xyz => coord_data%coord, n_neighbor => neigh_list%n_neighbor, neighbor => neigh_list%neighbors)

  do xbin = 1, xbin_max
  do ybin = 1, ybin_max
  do zbin = 1, zbin_max

      do atom = 1, cellbins(xbin, ybin, zbin)%n
        id = cellbins(xbin, ybin, zbin)%ids(atom)

        do p = -1, 1
        do q = -1, 1
        do o = -1, 1

        associate(checked_n => cellbins(xbin+p, ybin+q, zbin+o)%n,&
                  checked_ids => cellbins(xbin+p, ybin+q, zbin+o)%ids,&
                  x_pbc => cellbins(xbin+p, ybin+q, zbin+o)%x_pbc, &
                  y_pbc => cellbins(xbin+p, ybin+q, zbin+o)%y_pbc, &
                  z_pbc => cellbins(xbin+p, ybin+q, zbin+o)%z_pbc)

          do atom2 = 1, checked_n

            checkid = checked_ids(atom2)
            if (checkid < id) then   !to avoid repeat calculation

              d = (xyz(checkid,1) - x_pbc*coord_data%lx - xyz(id,1))**2& !x
                    +(xyz(checkid,2) - y_pbc*coord_data%ly - xyz(id,2))**2& !y
                    +(xyz(checkid,3) - z_pbc*coord_data%lz - xyz(id,3))**2  !z

              if (d < r) then
                n_neighbor(checkid) = n_neighbor(checkid) + 1
                n_neighbor(id) = n_neighbor(id) + 1
                neighbor(checkid, n_neighbor(checkid)) = id
                neighbor(id, n_neighbor(id)) = checkid
              endif
            endif
          end do

        end associate

        end do
        end do
        end do
      end do

  end do
  end do
  end do

  call print_cn

    do i = 1, natom
      call bubble_sort(n_neighbor(i), neighbor(i,:))
    enddo

  end associate

END SUBROUTINE find_neighbors

SUBROUTINE find_neighbors_d()
  IMPLICIT NONE
  !
  type(bins(capacity = :)), allocatable, dimension(:,:,:) :: cellbins
  integer :: xbin_max, ybin_max, zbin_max
  integer :: xbin, ybin, zbin, atom, atom2, i, p, q, o
  integer :: id, checkid
  real(dp) :: d, r, x_tmp, y_tmp, z_tmp

  r = cutoff**2

  call create_bins(cutoff, cellbins, xbin_max, ybin_max, zbin_max)
  print *, 'Constructing neighbor list from bins data ...'

  allocate(neighbor_list(atom_number = natom, capacity = 10) :: neigh_list)
  allocate(delta(natom, 10, 3), STAT=ierr, ERRMSG=emsg)

  neigh_list%n_neighbor = 0
  neigh_list%neighbors = 0
  delta = 0

  associate(xyz => coord_data%coord, n_neighbor => neigh_list%n_neighbor, neighbor => neigh_list%neighbors)

  do xbin = 1, xbin_max
  do ybin = 1, ybin_max
  do zbin = 1, zbin_max

      do atom = 1, cellbins(xbin, ybin, zbin)%n
        id = cellbins(xbin, ybin, zbin)%ids(atom)

        do p = -1, 1
        do q = -1, 1
        do o = -1, 1

        associate(checked_n => cellbins(xbin+p, ybin+q, zbin+o)%n,&
                  checked_ids => cellbins(xbin+p, ybin+q, zbin+o)%ids,&
                  x_pbc => cellbins(xbin+p, ybin+q, zbin+o)%x_pbc, &
                  y_pbc => cellbins(xbin+p, ybin+q, zbin+o)%y_pbc, &
                  z_pbc => cellbins(xbin+p, ybin+q, zbin+o)%z_pbc)

          do atom2 = 1, checked_n

            checkid = checked_ids(atom2)
            if (checkid < id) then   !to avoid repeat calculation
              x_tmp = xyz(checkid,1) - x_pbc*coord_data%lx
              y_tmp = xyz(checkid,2) - y_pbc*coord_data%ly
              z_tmp = xyz(checkid,3) - z_pbc*coord_data%lz

              d = (x_tmp - xyz(id,1))**2 + (y_tmp - xyz(id,2))**2 + (z_tmp - xyz(id,3))**2

              if (d < r) then
                n_neighbor(checkid) = n_neighbor(checkid) + 1
                n_neighbor(id) = n_neighbor(id) + 1
                neighbor(checkid, n_neighbor(checkid)) = id
                neighbor(id, n_neighbor(id)) = checkid

                delta(id,n_neighbor(id),1) = x_tmp - xyz(id,1)
                delta(id,n_neighbor(id),2) = y_tmp - xyz(id,2)
                delta(id,n_neighbor(id),3) = z_tmp - xyz(id,3)

                delta(checkid,n_neighbor(checkid),1) = xyz(id,1) - x_tmp
                delta(checkid,n_neighbor(checkid),2) = xyz(id,2) - y_tmp
                delta(checkid,n_neighbor(checkid),3) = xyz(id,3) - z_tmp
              endif
            endif
          end do

        end associate

        end do
        end do
        end do
      end do

  end do
  end do
  end do

  call print_cn

    do i = 1, natom
      call bubble_sort(n_neighbor(i), neighbor(i,:))
    enddo

  end associate

END SUBROUTINE find_neighbors_d

subroutine print_cn()
  implicit none
  integer :: maxn, i
  integer, allocatable :: rank(:), amount(:)

  associate(n_neighbor => neigh_list%n_neighbor)
    maxn = maxval(n_neighbor)

    if (maxn > 12) maxn = 12

    allocate(rank(0:maxn), amount(0:maxn))

    rank = [(i, i=0, maxn)]
    do i = 0, maxn
      amount(i) = count(n_neighbor == i)
    end do

    print *, ' ### Coordination Distribution'
    print *, '***************************'
    print 117, ' | CN    | ',rank
    print 117, ' | Count | ',amount
    print *, '***************************'
    117 format (a11,*(i6, ' | '))

  end associate
end subroutine print_cn

subroutine clean_neighbor
  implicit none

  if (allocated(neigh_list)) deallocate(neigh_list)
  if (allocated(delta)) deallocate(delta)

end subroutine clean_neighbor

end module neighbor_finder
