
! This is the old neighbor finder subroutine. Need to be put into the neighbor
! finder module when using.
SUBROUTINE  neighbor_finder_old
  IMPLICIT NONE
  integer(inp) :: o, p, q, xbin_max, ybin_max, zbin_max
  integer(inp), allocatable :: xbin(:), ybin(:), zbin(:), id(:)
  real(dp), allocatable :: realxyz(:,:),xyz(:,:)
  LOGICAL :: isghost = .false.
  REAL(dp) :: r,d,x_min,y_min,z_min
  integer(inp) :: i,j,n,l,k,xpbc=0,ypbc=0,zpbc=0&
    &,xbin_t,ybin_t,zbin_t

  type(bins(capacity = :)), allocatable, dimension(:,:,:) :: cellbins

  print *, 'Entering neighbor finder subroutine ...'

  allocate(xbin(natom), STAT=ierr, ERRMSG=emsg)
  allocate(ybin(natom), STAT=ierr, ERRMSG=emsg)
  allocate(zbin(natom), STAT=ierr, ERRMSG=emsg)

  allocate(id(natom*2), STAT=ierr, ERRMSG=emsg)
  allocate(realxyz(natom, 3), STAT=ierr, ERRMSG=emsg)
  allocate(xyz(natom*2,3), STAT=ierr, ERRMSG=emsg)

  allocate(neighbor_list(atom_number = natom, capacity = 10) :: neigh_list)

  associate(xyz0 => coord_data%coord, n_neighbor => neigh_list%n_neighbor, neighbor => neigh_list%neighbors,&
            lx => coord_data%lx, ly => coord_data%ly, lz => coord_data%lz)

  x_min=minval(xyz0(:,1))
  y_min=minval(xyz0(:,2))
  z_min=minval(xyz0(:,3))

  realxyz(:,1)=xyz0(:,1)-x_min
  realxyz(:,2)=xyz0(:,2)-y_min
  realxyz(:,3)=xyz0(:,3)-z_min

  r = cutoff**2
  xyz(:natom,:) = realxyz(:,:)

  xbin_max = CEILING(lx/cutoff)-1
  ybin_max = CEILING(ly/cutoff)-1
  zbin_max = CEILING(lz/cutoff)-1

  allocate(bins(30) :: cellbins(0:xbin_max+1,0:ybin_max+1,0:zbin_max+1))

  j = 1
  do i = 1, natom
    id(i) = i
    ! print *, id(i)
    xbin(i) = CEILING(xyz(i,1)/cutoff)
    ybin(i) = CEILING(xyz(i,2)/cutoff)
    zbin(i) = CEILING(xyz(i,3)/cutoff)

    if (xbin(i) > xbin_max) xbin(i) = xbin_max
    if (ybin(i) > ybin_max) ybin(i) = ybin_max
    if (zbin(i) > zbin_max) zbin(i) = zbin_max
    if (xbin(i) == 0) xbin(i) = 1
    if (ybin(i) == 0) ybin(i) = 1
    if (zbin(i) == 0) zbin(i) = 1

    cellbins(xbin(i),ybin(i),zbin(i))%n = cellbins(xbin(i),ybin(i),zbin(i))%n+1
    cellbins(xbin(i),ybin(i),zbin(i))%ids(cellbins(xbin(i),ybin(i),zbin(i))%n) = i

    n=1
    l=1
    k=1
    if (pbc == 1) then
      if (xbin(i) == 1) then
        xpbc = 1
        isghost = .true.
      else if (xbin(i) == xbin_max) then
        xpbc = -1
        n=-1
        isghost = .true.
      endif
      if (ybin(i) == 1) then
        ypbc = 1
        isghost = .true.
      else if (ybin(i) == ybin_max) then
        ypbc = -1
        l=-1
        isghost = .true.
      endif
      if (zbin(i) == 1) then
        zpbc = 1
        isghost = .true.
      else if (zbin(i) == zbin_max) then
        zpbc = -1
        k=-1
        isghost = .true.
      endif
      if (isghost) then
        !so natom+j is the mask id of atom i
        do p=0,xpbc,n
          do q=0,ypbc,l
            do o=0,zpbc,k
              if (p==0.and.q==0.and.o==0) then
                cycle
              endif

              id(natom+j) = i
              ! print *, id(natom+j)
              xbin_t = xbin(i)+p*xbin_max
              ybin_t = ybin(i)+q*ybin_max
              zbin_t = zbin(i)+o*zbin_max

              xyz(natom+j,1) = xyz(i,1) + lx*p
              xyz(natom+j,2) = xyz(i,2) + ly*q
              xyz(natom+j,3) = xyz(i,3) + lz*o

              cellbins(xbin_t,ybin_t,zbin_t)%n= cellbins(xbin_t,ybin_t,zbin_t)%n+1
              cellbins(xbin_t,ybin_t,zbin_t)%ids(cellbins(xbin(i),ybin(i),zbin(i))%n) = natom+j
              j=j+1

            enddo
          enddo
        enddo
        xpbc = 0
        ypbc = 0
        zpbc = 0
        isghost=.false.
      endif
    endif

  end do

  print *, 'sort atoms to bins--done, total number:', natom+j-1

  n_neighbor = 0
  neighbor = 0

  !!!!!!!!!
  ! This part can be looped according to bins, to enable parallel computing.
  do i = 1, natom

      do p = -1, 1
      do q = -1, 1
      do o = -1, 1
        if (cellbins(xbin(i)+p,ybin(i)+q,zbin(i)+o)%n /= 0) then
          do k = 1, cellbins(xbin(i)+p,ybin(i)+q,zbin(i)+o)%n
            l = cellbins(xbin(i)+p,ybin(i)+q,zbin(i)+o)%ids(k)

            if (i < id(l)) then   !to avoid repeat calculation
              d = (xyz(i,1)-xyz(l,1))**2& !x
                    &+(xyz(i,2)-xyz(l,2))**2& !y
                    &+(xyz(i,3)-xyz(l,3))**2  !z

              if (d < r) then
                n_neighbor(id(i)) = n_neighbor(id(i)) + 1
                n_neighbor(id(l)) = n_neighbor(id(l)) + 1
                neighbor(id(i),n_neighbor(id(i))) = id(l)
                neighbor(id(l),n_neighbor(id(l))) = id(i)

!                 delta(id(i),n_neighbor(id(i)),1)=xyz(l,1)-xyz(i,1)
!                 delta(id(i),n_neighbor(id(i)),2)=xyz(l,2)-xyz(i,2)
!                 delta(id(i),n_neighbor(id(i)),3)=xyz(l,3)-xyz(i,3)
!
!                 delta(id(l),n_neighbor(id(l)),1)=xyz(i,1)-xyz(l,1)
!                 delta(id(l),n_neighbor(id(l)),2)=xyz(i,2)-xyz(l,2)
!                 delta(id(l),n_neighbor(id(l)),3)=xyz(i,3)-xyz(l,3)
              endif
            endif
          enddo
        endif
      end do
      end do
      end do
  end do
  !!!!!!!!!

  call print_cn

  do i = 1, natom
    call bubble_sort(n_neighbor(i), neighbor(i,:))
  enddo

  end associate

print *, 'Leaving neighbor finder subroutine, constructing neighbor list -- done'

END SUBROUTINE neighbor_finder_old
