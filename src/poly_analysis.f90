! f2py3 -c [file_name] -m topology
! python3 [python_interface]

! coarse grain neighbors
module poly_analysis
  use precision
  use data_types
  use data_input, only: coord_data, o_type, natom
  use neighbor_finder, only: neigh_list
  implicit none

  ! polyhedron
  type polyhedron(atom_number, capacity)
      integer, len :: atom_number, capacity
      integer, dimension(atom_number) :: center_type = 0
      integer, dimension(atom_number) :: topo_type = 0
      integer, dimension(atom_number) :: poly_neigh_number = 0
      integer, dimension(atom_number, capacity) :: poly_list = 0
  end type

  type(polyhedron(atom_number = :, capacity = :)), allocatable :: polys

contains

SUBROUTINE poly_neighbor()
  !!!!!!!!!!!!!!!!!
  !
  ! This subroutine computes the neighbor number and neighbor list(poly_n, poly_list) of
  ! polyhedra from given atom neighbor number and neighbor list(ref_n, ref_list). type_x
  ! and type_o are either 1 or 0, given by io module so that this subroutine knows the
  ! atom is which type by comparing ptype with them.
  !
  !!!!!!!!!!!!!!!!!!
  IMPLICIT NONE
  integer(inp) :: m,n,i,j, k, same_neigh_num
  integer(inp), pointer :: ref_n(:), ref_list(:,:)
  print *, 'Polyhedral analysis subroutine ... Start'

  allocate(polyhedron(atom_number = natom, capacity = 20) :: polys, STAT=ierr, ERRMSG=emsg)

  associate(ref_n => neigh_list%n_neighbor, ref_list => neigh_list%neighbors, &
            ln => polys%topo_type, ptype => coord_data%ptype, &
            poly_n => polys%poly_neigh_number, poly_list => polys%poly_list)

  do i = 1, natom

    if (ptype(i)==o_type) cycle   !oxygen is not center atom
    do n = 1, ref_n(i)
      if (ptype(ref_list(i,n))==o_type) then
        do m = 1, ref_n(ref_list(i,n))
          if (ptype(ref_list(ref_list(i,n),m))/=o_type .and. ref_list(ref_list(i,n),m)/=i) then
            k = 1
            do
              if (poly_n(i)==0) then
                poly_n(i)=poly_n(i)+1
                poly_list(i,1)=ref_list(ref_list(i,n),m)
                exit
              endif
              if (ref_list(ref_list(i,n),m)==poly_list(i,k)) exit
              k = k + 1
              ! print *, k
              if (k > poly_n(i)) then
                poly_n(i) = poly_n(i)+1
                poly_list(i,k) = ref_list(ref_list(i,n),m)
                if (poly_n(i) == 20) then
                  print *, 'Warning: maximum polyhedral neighbor length reached.'
                endif
                exit
              endif
            enddo
          endif
        enddo
      endif
    enddo

    call bubble_sort(poly_n(i), poly_list(i,:))

  enddo

! Check the neighbor list and count number of shared neighbor atom.
  do i = 1,natom

    if (poly_n(i)/=0) then
      do j =1, poly_n(i)
        same_neigh_num = 0
        m=1
        n=1
      do while(m <= ref_n(i) .and. n <= ref_n(poly_list(i,j)) )
  !       print *, m, n
            if (ref_list(i,m) > ref_list(poly_list(i,j),n)) then
              n = n + 1
            else if (ref_list(i,m) < ref_list(poly_list(i,j),n)) then
              m = m + 1
            else if (ref_list(i,m) == ref_list(poly_list(i,j),n)) then
              m = m + 1
              n = n + 1
              same_neigh_num = same_neigh_num + 1
            endif
        enddo

        if (ln(i)< same_neigh_num) then
          ln(i)= same_neigh_num
        endif
      enddo
    endif
  enddo

    print *, "      CS       ES       FS"
    print *, count(ln==1), "   ",count(ln==2),"    ",count(ln==3)
    print *, 'Polyhedral analysis ... Done'
  end associate
END SUBROUTINE poly_neighbor


subroutine neighbor_change(natom,old_n,old_list,new_n,new_list,n_change)
  ! calculate the exclusive-or elements of old and new neighbor list, count them as the neighbor change.

    implicit none
!     integer, parameter :: dp = selected_real_kind(15, 307)
    integer(inp), intent(in) :: natom, old_n(natom), new_n(natom), &
          &old_list(natom,20), new_list(natom,20)
    integer(inp), intent(out) :: n_change(natom)
    integer(inp) :: i, m, n,j
    integer(inp), allocatable :: unshared(:)

do i = 1, natom
  ! print *, i
  ! call exclusiveor(old_n(i),old_list(i),new_n(i),new_list(i),n_change(i),unshared(i))
    allocate(unshared(old_n(i)+new_n(i)))
    m = 1
    n = 1
    j = 1 !exclusive-or count
    do while(m <= old_n(i) .and. n <= new_n(i) )
        if (old_list(i,m) > new_list(i,n)) then
            unshared(j) = new_list(i,n)
           n = n + 1
           j = j + 1
        else if (old_list(i,m) < new_list(i,n)) then
            unshared(j) = old_list(i,m)
           m = m + 1
           j = j + 1
        else if (old_list(i,m) == new_list(i,n)) then
           m = m + 1
           n = n + 1
        endif
    enddo
    do while(m <= old_n(i))
        unshared(j) = old_list(i,m)
       m = m + 1
       j = j + 1
    enddo
    do while(n <= new_n(i))
        unshared(j) = new_list(i,n)
       n = n + 1
       j = j + 1
    enddo
    j=j-1   !j is the length of changed neighbors
    m=1
    n=1
    n_change(i) = j
    deallocate(unshared)
enddo
    return
end subroutine neighbor_change

! get coarse grain neighbor change event number and IPR.
subroutine neighbor_change_ipr(natom,old_n,old_list,new_n,new_list,ipr_list,n_change,ipr)
    implicit none
    integer(inp), intent(in) :: natom, old_n(natom), new_n(natom), &
          &old_list(natom,20), new_list(natom,20)
    integer(inp), intent(out) :: n_change(natom)
    real(dp), intent(out) :: ipr(natom)
    integer(8),intent(inout) :: ipr_list(natom,200,2)
    integer(inp) :: i, m, n,j
    integer(8) :: nu,de
    integer(inp), allocatable :: unshared(:)
    ! calculate the exclusive-or elements of old and new neighbor list, count them as the neighbor change.
    do i = 1, natom
      ! print *, i
      ! call exclusiveor(old_n(i),old_list(i),new_n(i),new_list(i),n_change(i),unshared(i))
        allocate(unshared(old_n(i)+new_n(i)))
        m = 1
        n = 1
        j = 1 !exclusive-or count
        nu=0
        de=0
        do while(m <= old_n(i) .and. n <= new_n(i) )
            if (old_list(i,m) > new_list(i,n)) then
                unshared(j) = new_list(i,n)
              n = n + 1
              j = j + 1
            else if (old_list(i,m) < new_list(i,n)) then
                unshared(j) = old_list(i,m)
              m = m + 1
              j = j + 1
            else if (old_list(i,m) == new_list(i,n)) then
              m = m + 1
              n = n + 1
            endif
        enddo
        do while(m <= old_n(i))
            unshared(j) = old_list(i,m)
          m = m + 1
          j = j + 1
        enddo
        do while(n <= new_n(i))
            unshared(j) = new_list(i,n)
          n = n + 1
          j = j + 1
        enddo
        j=j-1   !j is the length of changed neighbors
        m=1
        n=1

        do while (m<=j)
          ! print *, unshared(m),ipr_list(i,n,1)
          if (unshared(m)>ipr_list(i,n,1)) then
            if (ipr_list(i,n,1)==0) then
              ipr_list(i,n,1)=unshared(m)
              ipr_list(i,n,2)=ipr_list(i,n,2)+1
              m=m+1
            endif
            n=n+1
          elseif (unshared(m)<ipr_list(i,n,1)) then
            ipr_list(i,n+1:200,:)=ipr_list(i,n:199,:) !move the list one step to right
            ipr_list(i,n,1)=unshared(m)
            ipr_list(i,n,2)=ipr_list(i,n,2)+1
            m=m+1
            n=n+1
          else
            ipr_list(i,n,2)=ipr_list(i,n,2)+1
            m=m+1
            n=n+1
          endif
        enddo
        n_change(i) = j
        do n=1,200
          ! print *, n
          nu=nu+(ipr_list(i,n,2))**4
          de=de+(ipr_list(i,n,2))**2
        enddo
        ! print *, nu, de
        if (de/=0) then
          ipr(i)=nu*1.0_dp/(de*1.0_dp)**2
        endif
        ! print *, ipr(i)
        deallocate(unshared)
    enddo

end subroutine neighbor_change_ipr

!!!!!!!!
! This subroutine compares ln at one frame and a reference frame and conclude
! whether the polyhedron changed its ln. (a bit ambiguous because unlike coordination change,
! no changing direction is specified.)
!!!!!!!!
subroutine ln_change(natom,new_ln,ref_ln,out)
  implicit none
  integer(inp),intent(in) :: natom,new_ln(natom),ref_ln(natom)
  integer(inp),intent(out) :: out(natom,4)
  integer(inp) :: i
  ! print *, nc(:100)
  out=0
  ! print *, ref_ln(:100), new_ln(:100)

  do i=1,natom
    if (ref_ln(i)==1 .and. new_ln(i)==1) then
      out(i,1)=out(i,1)+1
    elseif (ref_ln(i)==1 .and. new_ln(i)==2) then
      out(i,2)=out(i,2)+1
    elseif (ref_ln(i)==2 .and. new_ln(i)==1) then
      out(i,3)=out(i,3)+1
    elseif (ref_ln(i)==2 .and. new_ln(i)==2) then
      out(i,4)=out(i,4)+1
    endif

  enddo

  return
end subroutine ln_change


! If atom number is N, there is a list(N) to store the flexibility count.
! It is the number that a polyhedron has changed its type during deformation.
! Devide the count M with frame number Z, we have M/Z, which is the flexibility.
! subroutine flexible()
!   implicit none
!
!   if (.not. allocated) allocate(flex_list(atom_number = ..) :: flex)
!
!   do i = 1, atom_number
!     if (flex%ln /= ln(i)) flex%count = flex%count + 1
!   end do
!
! end subroutine flexible

end module poly_analysis
