module neigh_finder
  use precision
  use io, xyz0 => xyz, atom_number => natom
  implicit none
  save
  type bins(cap)
      integer(4), len :: cap   !capacity of the bin
      integer(4) :: n, n_list(cap)  !id of the bin
  end type
!   note: for gfortran 13.2.1, deallocate a array of bins type will result in error. Thus in
! SUBROUTINEs below, the deallocate command is comment out at this moment.
  INTEGER(4), allocatable :: neighbor(:,:), n_neighbor(:)
  real(dp), allocatable :: delta(:,:,:)
                          ! neighbor(natom, 10), n_neighbor(natom)

  contains

  SUBROUTINE  neighbor_finder
    IMPLICIT NONE
    INTEGER(4) :: o, p, q, xbin_max, ybin_max, zbin_max
    integer(4) :: bincap = 10
    integer(4), allocatable :: xbin(:), ybin(:), zbin(:), id(:)
    real(dp), allocatable :: realxyz(:,:),xyz(:,:)
    LOGICAL :: isghost = .false.
    REAL(dp) :: r,d,x_min,y_min,z_min
    integer(4) :: i,j,n,l,k,xpbc=0,ypbc=0,zpbc=0&
      &,xbin_t,ybin_t,zbin_t
    type(bins(bincap)), allocatable :: cellbins(:,:,:)

    call get_xyz

    allocate(xbin(atom_number), STAT=ierr, ERRMSG=emsg)
    allocate(ybin(atom_number), STAT=ierr, ERRMSG=emsg)
    allocate(zbin(atom_number), STAT=ierr, ERRMSG=emsg)

    allocate(id(atom_number*2), STAT=ierr, ERRMSG=emsg)
    allocate(realxyz(atom_number, 3), STAT=ierr, ERRMSG=emsg)
    allocate(xyz(atom_number*2,3), STAT=ierr, ERRMSG=emsg)

    allocate(neighbor(atom_number, 10), STAT=ierr, ERRMSG=emsg)
    allocate(n_neighbor(atom_number), STAT=ierr, ERRMSG=emsg)


    x_min=minval(xyz0(:,1))
    y_min=minval(xyz0(:,2))
    z_min=minval(xyz0(:,3))

    realxyz(:,1)=xyz0(:,1)-x_min
    realxyz(:,2)=xyz0(:,2)-y_min
    realxyz(:,3)=xyz0(:,3)-z_min

    r = cutoff**2
    xyz(:atom_number,:) = realxyz(:,:)

    xbin_max = CEILING(lx/cutoff)-1
    ybin_max = CEILING(ly/cutoff)-1
    zbin_max = CEILING(lz/cutoff)-1

    allocate(cellbins(0:xbin_max+1,0:ybin_max+1,0:zbin_max+1))

! initialize cell index
    do j = 0, xbin_max+1
      do i = 0,ybin_max+1
        do n=0,zbin_max+1
          cellbins(j,i,n)%n = 0
          cellbins(j,i,n)%n_list=0
        enddo
      enddo
    enddo

    j = 1
    do i = 1, atom_number
      id(i) = i
      ! print *, id(i)
      xbin(i) = CEILING(xyz(i,1)/cutoff)
      ybin(i) = CEILING(xyz(i,2)/cutoff)
      zbin(i) = CEILING(xyz(i,3)/cutoff)

      if (xbin(i) > xbin_max) then
        xbin(i) = xbin_max
      endif
      if (ybin(i) > ybin_max) then
        ybin(i) = ybin_max
      endif
      if (zbin(i) > zbin_max) then
        zbin(i) = zbin_max
      endif
      if (xbin(i) < 1) then
        xbin(i) = 1
      endif
      if (ybin(i) < 1) then
        ybin(i) = 1
      endif
      if (zbin(i) < 1) then
        zbin(i) = 1
      endif

      cellbins(xbin(i),ybin(i),zbin(i))%n = cellbins(xbin(i),ybin(i),zbin(i))%n+1
      cellbins(xbin(i),ybin(i),zbin(i))%n_list(cellbins(xbin(i)&
      &,ybin(i),zbin(i))%n) = i
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
        if (isghost .eqv. .true.) then
          !so atom_number+j is the mask id of atom i
          do p=0,xpbc,n
            do q=0,ypbc,l
              do o=0,zpbc,k
                if (p==0.and.q==0.and.o==0) then
                  cycle
                endif

                id(atom_number+j) = i
                ! print *, id(atom_number+j)
                xbin_t = xbin(i)+p*xbin_max
                ybin_t = ybin(i)+q*ybin_max
                zbin_t = zbin(i)+o*zbin_max

                xyz(atom_number+j,1) = xyz(i,1) + lx*p
                xyz(atom_number+j,2) = xyz(i,2) + ly*q
                xyz(atom_number+j,3) = xyz(i,3) + lz*o

                cellbins(xbin_t,ybin_t,zbin_t)%n= cellbins(xbin_t,ybin_t,zbin_t)%n+1
                cellbins(xbin_t,ybin_t,zbin_t)%n_list(cellbins(xbin(i),ybin(i),&
                &zbin(i))%n) = atom_number+j
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

    print *, 'sort atoms to bins--done, total number:', atom_number+j-1

    n_neighbor=0
    neighbor=0

    do i = 1, atom_number
        n = 1
        do p = -1, 1
          do q = -1, 1
            do o = -1, 1
              if (cellbins(xbin(i)+p,ybin(i)+q,zbin(i)+o)%n /= 0) then
                do k = 1, cellbins(xbin(i)+p,ybin(i)+q,zbin(i)+o)%n
                  l = cellbins(xbin(i)+p,ybin(i)+q,zbin(i)+o)%n_list(k)

                  if (i < id(l)) then   !to avoid repeat calculation
                    d = (xyz(i,1)-xyz(l,1))**2& !x
                          &+(xyz(i,2)-xyz(l,2))**2& !y
                          &+(xyz(i,3)-xyz(l,3))**2  !z
                        if (d < r) then
                          n_neighbor(id(i)) = n_neighbor(id(i)) + 1
                          n_neighbor(id(l)) = n_neighbor(id(l)) + 1
                          neighbor(id(i),n_neighbor(id(i))) = id(l)
                          neighbor(id(l),n_neighbor(id(l))) = id(i)
                        endif
                  endif
                enddo
              endif
            end do
          end do
        end do
    end do
!     print *, 'constructing neighbor list--done'

    do i = 1, atom_number
      call bubble_sort(n_neighbor(i), neighbor(i,:))
    enddo
!   deallocate(cellbins)
  deallocate(xbin)
  deallocate(ybin)
  deallocate(zbin)

  deallocate(id)
  deallocate(realxyz)
  deallocate(xyz)

    return

  END SUBROUTINE

  SUBROUTINE neighbor_finder_d
    IMPLICIT NONE
    integer(4) :: bincap = 10
    INTEGER(4) :: o, p, q, xbin_max, ybin_max, zbin_max
    integer(4), allocatable :: xbin(:), ybin(:), zbin(:), id(:)
    real(dp), allocatable :: realxyz(:,:),xyz(:,:)
    LOGICAL :: isghost = .false.
    REAL(dp) :: r,d,x_min,y_min,z_min
    integer(4) :: i,j,n,l,k,xpbc=0,ypbc=0,zpbc=0&
      &,xbin_t,ybin_t,zbin_t
    type(bins(bincap)), allocatable :: cellbins(:,:,:)

    call get_xyz

    allocate(xbin(atom_number), STAT=ierr, ERRMSG=emsg)
    allocate(ybin(atom_number), STAT=ierr, ERRMSG=emsg)
    allocate(zbin(atom_number), STAT=ierr, ERRMSG=emsg)

    allocate(id(atom_number*2), STAT=ierr, ERRMSG=emsg)
    allocate(realxyz(atom_number, 3), STAT=ierr, ERRMSG=emsg)
    allocate(xyz(atom_number*2,3), STAT=ierr, ERRMSG=emsg)

    allocate(neighbor(atom_number, 10), STAT=ierr, ERRMSG=emsg)
    allocate(n_neighbor(atom_number), STAT=ierr, ERRMSG=emsg)

    allocate(delta(atom_number,10,3), STAT=ierr, ERRMSG=emsg)

        x_min=minval(xyz0(:,1))
    y_min=minval(xyz0(:,2))
    z_min=minval(xyz0(:,3))

    realxyz(:,1)=xyz0(:,1)-x_min
    realxyz(:,2)=xyz0(:,2)-y_min
    realxyz(:,3)=xyz0(:,3)-z_min

    r = cutoff**2
    xyz(:atom_number,:) = realxyz(:,:)

    xbin_max = CEILING(lx/cutoff)-1
    ybin_max = CEILING(ly/cutoff)-1
    zbin_max = CEILING(lz/cutoff)-1

    allocate(cellbins(0:xbin_max+1,0:ybin_max+1,0:zbin_max+1), STAT = ierr)
    print *, ierr
! initialize cell index
    do j = 0, xbin_max+1
      do i = 0,ybin_max+1
        do n=0,zbin_max+1
          cellbins(j,i,n)%n = 0
          cellbins(j,i,n)%n_list=0
        enddo
      enddo
    enddo

    j = 1
    do i = 1, atom_number
      id(i) = i
      ! print *, id(i)
      xbin(i) = CEILING(xyz(i,1)/cutoff)
      ybin(i) = CEILING(xyz(i,2)/cutoff)
      zbin(i) = CEILING(xyz(i,3)/cutoff)

      if (xbin(i) > xbin_max) then
        xbin(i) = xbin_max
      endif
      if (ybin(i) > ybin_max) then
        ybin(i) = ybin_max
      endif
      if (zbin(i) > zbin_max) then
        zbin(i) = zbin_max
      endif
      if (xbin(i) < 1) then
        xbin(i) = 1
      endif
      if (ybin(i) < 1) then
        ybin(i) = 1
      endif
      if (zbin(i) < 1) then
        zbin(i) = 1
      endif

      cellbins(xbin(i),ybin(i),zbin(i))%n = cellbins(xbin(i),ybin(i),zbin(i))%n+1
      cellbins(xbin(i),ybin(i),zbin(i))%n_list(cellbins(xbin(i)&
      &,ybin(i),zbin(i))%n) = i
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
        if (isghost .eqv. .true.) then
          !so atom_number+j is the mask id of atom i
          do p=0,xpbc,n
            do q=0,ypbc,l
              do o=0,zpbc,k
                if (p==0.and.q==0.and.o==0) then
                  cycle
                endif

                id(atom_number+j) = i
                ! print *, id(atom_number+j)
                xbin_t = xbin(i)+p*xbin_max
                ybin_t = ybin(i)+q*ybin_max
                zbin_t = zbin(i)+o*zbin_max

                xyz(atom_number+j,1) = xyz(i,1) + lx*p
                xyz(atom_number+j,2) = xyz(i,2) + ly*q
                xyz(atom_number+j,3) = xyz(i,3) + lz*o

                cellbins(xbin_t,ybin_t,zbin_t)%n= cellbins(xbin_t,ybin_t,zbin_t)%n+1
                cellbins(xbin_t,ybin_t,zbin_t)%n_list(cellbins(xbin(i),ybin(i),&
                &zbin(i))%n) = atom_number+j
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

    print *, 'sort atoms to bins--done, total number:', atom_number+j-1

    n_neighbor=0
    neighbor=0

    do i = 1, atom_number
        n = 1
        do p = -1, 1
          do q = -1, 1
            do o = -1, 1
              if (cellbins(xbin(i)+p,ybin(i)+q,zbin(i)+o)%n /= 0) then
                do k = 1, cellbins(xbin(i)+p,ybin(i)+q,zbin(i)+o)%n
                  l = cellbins(xbin(i)+p,ybin(i)+q,zbin(i)+o)%n_list(k)

                  if (i < id(l)) then   !to avoid repeat calculation
                    d = (xyz(i,1)-xyz(l,1))**2& !x
                          &+(xyz(i,2)-xyz(l,2))**2& !y
                          &+(xyz(i,3)-xyz(l,3))**2  !z
                    if (d < r) then
                      n_neighbor(id(i)) = n_neighbor(id(i)) + 1
                      n_neighbor(id(l)) = n_neighbor(id(l)) + 1
                      neighbor(id(i),n_neighbor(id(i))) = id(l)
                      neighbor(id(l),n_neighbor(id(l))) = id(i)

                      delta(id(i),n_neighbor(id(i)),1)=xyz(l,1)-xyz(i,1)
                      delta(id(i),n_neighbor(id(i)),2)=xyz(l,2)-xyz(i,2)
                      delta(id(i),n_neighbor(id(i)),3)=xyz(l,3)-xyz(i,3)

                      delta(id(l),n_neighbor(id(l)),1)=xyz(i,1)-xyz(l,1)
                      delta(id(l),n_neighbor(id(l)),2)=xyz(i,2)-xyz(l,2)
                      delta(id(l),n_neighbor(id(l)),3)=xyz(i,3)-xyz(l,3)
                    endif
                  endif
                enddo
              endif
            end do
          end do
        end do
    end do
!     print *, 'constructing neighbor list--done'

    do i = 1, atom_number
      call bubble_sort(n_neighbor(i), neighbor(i,:))
    enddo

!   deallocate(cellbins)
  deallocate(xbin)
  deallocate(ybin)
  deallocate(zbin)

  deallocate(id)
  deallocate(realxyz)
  deallocate(xyz)

  return
  END SUBROUTINE

end module neigh_finder
