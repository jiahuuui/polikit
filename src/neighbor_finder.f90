module neighbor_finder
  use precision
  use data_types
  use parser, only: cutoffs, pbc
  use data_input, only: coord_data, natom, ntype

  implicit none
!   save
  ! neighbor list
  type neighbor_list(atom_number, capacity)
      integer, len :: atom_number, capacity
      integer, dimension(atom_number) :: n_neighbor = 0
      integer, dimension(atom_number, capacity) :: neighbors = 0
  end type

  ! bin type define, abolished due to stability issue.
  type bins(capacity)
      integer, len :: capacity   !capacity of the bin
      integer :: n = 0
      integer :: x_pbc = 0, y_pbc = 0, z_pbc = 0
      integer, dimension(capacity) :: ids = 0  !id of the atoms in bin
  end type

  real(dp), save, allocatable :: delta(:,:,:)

  type(neighbor_list(atom_number = :, capacity = :)), allocatable :: neigh_list

  integer, dimension(:,:,:), allocatable :: cells_n, cells_xpbc, cells_ypbc, cells_zpbc
  integer, dimension(:,:,:,:), allocatable :: cells_ids

  contains

! Divide the simulation box into bins, and put atoms into their corresponding bins.
SUBROUTINE create_bins(rCut, xbin_max, ybin_max, zbin_max)
  IMPLICIT NONE
  ! in:
  real(dp), intent(in) :: rCut
  ! inOUT:
  integer, intent(inout) :: xbin_max, ybin_max, zbin_max
  ! PRIVATE:
  real(dp), allocatable, dimension(:,:) :: realxyz
  REAL(dp) :: x_min, y_min, z_min
  INTEGER :: xbin, ybin, zbin, atom, bincap

199 format (a, i0)
  bincap = rCut**3

  print 199, " Capacity of bins is set to: ", bincap

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

    print 199, " Number of bins on dimension X: ", xbin_max
    print 199, " Number of bins on dimension Y: ", ybin_max
    print 199, " Number of bins on dimension Z: ", zbin_max

  end associate

  allocate(cells_n(0:xbin_max+1,0:ybin_max+1,0:zbin_max+1))
  allocate(cells_xpbc(0:xbin_max+1,0:ybin_max+1,0:zbin_max+1))
  allocate(cells_ypbc(0:xbin_max+1,0:ybin_max+1,0:zbin_max+1))
  allocate(cells_zpbc(0:xbin_max+1,0:ybin_max+1,0:zbin_max+1))
  allocate(cells_ids(0:xbin_max+1,0:ybin_max+1,0:zbin_max+1, bincap))

  print *, 'Deviding box into bins, memory cost of bins: ', sizeof(cells_ids)
  cells_n = 0
  cells_xpbc = 0
  cells_ypbc = 0
  cells_zpbc = 0
  cells_ids = 0

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

    cells_n(xbin,ybin,zbin) = cells_n(xbin,ybin,zbin) + 1
    if (cells_n(xbin,ybin,zbin) > bincap) then
        print *, " ! Bin capacity fulfilled due to unknown reason, stopping."
        stop
    end if
    cells_ids(xbin,ybin,zbin,cells_n(xbin,ybin,zbin)) = atom

  end do

  if (pbc == 1) then
    do xbin = 0, xbin_max + 1
    do ybin = 0, ybin_max + 1
    do zbin = 0, zbin_max + 1

      associate(x_pbc => cells_xpbc(xbin, ybin, zbin), &
                y_pbc => cells_ypbc(xbin, ybin, zbin), &
                z_pbc => cells_zpbc(xbin, ybin, zbin))

        if (xbin == 0) x_pbc = 1
        if (xbin == xbin_max + 1) x_pbc = -1
        if (ybin == 0) y_pbc = 1
        if (ybin == ybin_max + 1) y_pbc = -1
        if (zbin == 0) z_pbc = 1
        if (zbin == zbin_max + 1) z_pbc = -1
        cells_n(xbin, ybin, zbin) = &
        cells_n(xbin + x_pbc*xbin_max, ybin + y_pbc*ybin_max, zbin + z_pbc*zbin_max)
        cells_ids(xbin, ybin, zbin, :) = &
        cells_ids(xbin + x_pbc*xbin_max, ybin + y_pbc*ybin_max, zbin + z_pbc*zbin_max, :)
      end associate
    end do
    end do
    end do
  endif

  print *, 'Leaving bins construction subroutine ...'

END SUBROUTINE create_bins

SUBROUTINE find_neighbors(flag_delta)
  IMPLICIT NONE
  ! IN:
  logical, intent(in) :: flag_delta
  !
  integer :: xbin_max, ybin_max, zbin_max
  integer :: xbin, ybin, zbin, atom, atom2, i, p, q, o
  integer :: id, checkid, type1, type2
  real(dp) :: d, x_tmp, y_tmp, z_tmp
!   real(dp) :: r

  real(dp), allocatable :: r(:,:)

  allocate(r(ntype, ntype))

  if (size(cutoffs) == 1) then
    r = cutoffs(1)
  else
    i = 1
    do p = 1, ntype
      do q = p, ntype
        r(p,q) = cutoffs(i)
        r(q,p) = cutoffs(i)
        i = i+1
      end do
    end do
  end if

  r = r**2
  d = maxval(r(:,:))

  call create_bins(d, xbin_max, ybin_max, zbin_max)

  if (.not. allocated(neigh_list)) allocate(neighbor_list(atom_number = natom, capacity = 10) :: neigh_list)
  print *, 'Constructing neighbor list, memory cost: ', sizeof(neigh_list%neighbors)
  neigh_list%n_neighbor = 0
  neigh_list%neighbors = 0

  if ((.not. allocated(delta)) .and. flag_delta) then
    allocate(delta(natom, 10, 3), STAT=ierr, ERRMSG=emsg)
    delta = 0.
    print *, 'Constructing delta array, memory cost: ', sizeof(delta)
  end if

  associate(xyz => coord_data%coord, n_neighbor => neigh_list%n_neighbor, neighbor => neigh_list%neighbors)

  do xbin = 1, xbin_max
  do ybin = 1, ybin_max
  do zbin = 1, zbin_max

      do atom = 1, cells_n(xbin, ybin, zbin)
        id = cells_ids(xbin, ybin, zbin, atom)
        type1 = coord_data%ptype(id)

        do p = -1, 1
        do q = -1, 1
        do o = -1, 1

        associate(checked_n => cells_n(xbin+p, ybin+q, zbin+o),&
                  checked_ids => cells_ids(xbin+p, ybin+q, zbin+o, :),&
                  x_pbc => cells_xpbc(xbin+p, ybin+q, zbin+o), &
                  y_pbc => cells_ypbc(xbin+p, ybin+q, zbin+o), &
                  z_pbc => cells_zpbc(xbin+p, ybin+q, zbin+o))

          do atom2 = 1, checked_n
            checkid = checked_ids(atom2)
            type2 = coord_data%ptype(checkid)

            if (checkid < id) then   !to avoid repeat calculation

              x_tmp = xyz(checkid,1) - x_pbc*coord_data%lx
              y_tmp = xyz(checkid,2) - y_pbc*coord_data%ly
              z_tmp = xyz(checkid,3) - z_pbc*coord_data%lz

              d = (x_tmp - xyz(id,1))**2& !x
                    +(y_tmp - xyz(id,2))**2& !y
                    +(z_tmp - xyz(id,3))**2  !z

              if (d < r(type1, type2)) then
                n_neighbor(checkid) = n_neighbor(checkid) + 1
                n_neighbor(id) = n_neighbor(id) + 1
                neighbor(checkid, n_neighbor(checkid)) = id
                neighbor(id, n_neighbor(id)) = checkid

                if (flag_delta) then
                delta(id,n_neighbor(id),1) = x_tmp - xyz(id,1)
                delta(id,n_neighbor(id),2) = y_tmp - xyz(id,2)
                delta(id,n_neighbor(id),3) = z_tmp - xyz(id,3)

                delta(checkid,n_neighbor(checkid),1) = xyz(id,1) - x_tmp
                delta(checkid,n_neighbor(checkid),2) = xyz(id,2) - y_tmp
                delta(checkid,n_neighbor(checkid),3) = xyz(id,3) - z_tmp
                end if
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
    print 117, 'c| Count | ',amount
    print *, '***************************'
    117 format (a11,*(i6, ' | '))

  end associate
end subroutine print_cn

subroutine clean_bins()
  implicit none

  if (allocated(cells_n)) deallocate(cells_n)
  if (allocated(cells_xpbc)) deallocate(cells_xpbc)
  if (allocated(cells_ypbc)) deallocate(cells_ypbc)
  if (allocated(cells_zpbc)) deallocate(cells_zpbc)
  if (allocated(cells_ids)) deallocate(cells_ids)

end subroutine clean_bins

subroutine clean_neighbor
  implicit none

  call clean_bins()

  if (allocated(neigh_list)) then
    neigh_list%n_neighbor = 0
    neigh_list%neighbors = 0
  end if

  if (allocated(delta)) delta = 0.

end subroutine clean_neighbor

end module neighbor_finder
