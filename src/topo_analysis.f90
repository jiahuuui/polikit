! f2py3 -c [file_name] -m topology
! python3 [python_interface]

! coarse grain neighbors
module topo_analysis
  use precision
  use io
  use neigh_finder, ref_n => n_neighbor, ref_list => neighbor
  implicit none
  save
  type polyhedron
    integer(4) :: nei, neilist(20)
  end type

  type(polyhedron), allocatable :: poly(:)
  integer(4), allocatable :: poly_n(:), poly_list(:,:), ln(:)

contains

SUBROUTINE poly_neighbor
  !!!!!!!!!!!!!!!!!
  !
  ! This subroutine computes the neighbor number and neighbor list(poly_n, poly_list) of
  ! polyhedra from given atom neighbor number and neighbor list(ref_n, ref_list). type_x
  ! and type_o are either 1 or 0, given by io module so that this subroutine knows the
  ! atom is which type by comparing ptype with them.
  !
  !!!!!!!!!!!!!!!!!!
  IMPLICIT NONE
  integer(4) :: m,n,i,j,k

  call neighbor_finder   ! and get_xyz is called in neighbor_finder


  allocate(poly(natom), STAT=ierr, ERRMSG=emsg)
  allocate(ln(natom), STAT=ierr, ERRMSG=emsg)

  ln = 0
  forall (n=1: natom)
!   print *, n
    poly(n)%nei = 0
    poly(n)%neilist = 0
  end forall

  do i = 1,natom
!     print *, i
    if (ptype(i)==o_type) cycle   !oxygen is not center atom
    do n=1, ref_n(i)
      if (ptype(ref_list(i,n))==o_type) then
        do m=1, ref_n(ref_list(i,n))
          if (ptype(ref_list(ref_list(i,n),m))==x_type .and. ref_list(ref_list(i,n),m)/=i) then
            k=1
            do
              if (poly(i)%nei==0) then
                poly(i)%nei=poly(i)%nei+1
                poly(i)%neilist(1)=ref_list(ref_list(i,n),m)
                exit
              endif
              if (ref_list(ref_list(i,n),m)==poly(i)%neilist(k)) exit
              k=k+1
              ! print *, k
              if (k>poly(i)%nei) then
                poly(i)%nei=poly(i)%nei+1
                poly(i)%neilist(k)=ref_list(ref_list(i,n),m)
                if (poly(i)%nei==20) then
                  print *, 'warning: maximum polyhedral neighbor length reached.'
                endif
                exit
              endif
            enddo
          endif
        enddo
      endif
    enddo

    call bubble_sort(poly(i)%nei, poly(i)%neilist)

  enddo

! print *, "poly neighbor constructed."
ln = 0
do i = 1,natom

  if (poly(i)%nei/=0) then
    do j =1, poly(i)%nei
      k=0
      m=1
      n=1
     do while(m <= ref_n(i) .and. n <= ref_n(poly(i)%neilist(j)) )
!       print *, m, n
          if (ref_list(i,m) > ref_list(poly(i)%neilist(j),n)) then
             n = n + 1
          else if (ref_list(i,m) < ref_list(poly(i)%neilist(j),n)) then
             m = m + 1
          else if (ref_list(i,m) == ref_list(poly(i)%neilist(j),n)) then
             m = m + 1
             n = n + 1
             k=k+1
          endif
      enddo

      if (ln(i)<k) then
        ln(i)=k
      endif
    enddo
  endif
enddo

  print *, 'polyhedra neighbor analysis -- done'
  return
END SUBROUTINE poly_neighbor


subroutine neighbor_change(natom,old_n,old_list,new_n,new_list,n_change)
  ! calculate the exclusive-or elements of old and new neighbor list, count them as the neighbor change.

    implicit none
!     integer, parameter :: dp = selected_real_kind(15, 307)
    integer(4), intent(in) :: natom, old_n(natom), new_n(natom), &
          &old_list(natom,20), new_list(natom,20)
    integer(4), intent(out) :: n_change(natom)
    integer(4) :: i, m, n,j
    integer(4), allocatable :: unshared(:)

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
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer(4), intent(in) :: natom, old_n(natom), new_n(natom), &
          &old_list(natom,20), new_list(natom,20)
    integer(4), intent(out) :: n_change(natom)
    real(dp), intent(out) :: ipr(natom)
    integer(8),intent(inout) :: ipr_list(natom,200,2)
    integer(4) :: i, m, n,j
    integer(8) :: nu,de
    integer(4), allocatable :: unshared(:)
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
    return

end subroutine neighbor_change_ipr

subroutine ln_analysis(natom,new_ln,ref_ln,nc,out)
!!!!!!!!
! This subroutine compares ln at one frame and a reference frame and conclude
! whether the polyhedron changed its ln. (a bit ambiguous because unlike coordination change, no changing direction is specified.)
!!!!!!!!
implicit none
integer(4),intent(in) :: natom,new_ln(natom),ref_ln(natom),nc(natom)
integer(4),intent(out) :: out(natom,4)
integer(4) :: i
! print *, nc(:100)
out=0
! print *, ref_ln(:100), new_ln(:100)

do i=1,natom
  if (ref_ln(i)==1 .and. new_ln(i)==1) then
    out(i,1)=out(i,1)+nc(i)
  elseif (ref_ln(i)==1 .and. new_ln(i)==2) then
    out(i,2)=out(i,2)+nc(i)
  elseif (ref_ln(i)==2 .and. new_ln(i)==1) then
    out(i,3)=out(i,3)+nc(i)
  elseif (ref_ln(i)==2 .and. new_ln(i)==2) then
    out(i,4)=out(i,4)+nc(i)
  endif

enddo

return
end subroutine ln_analysis

subroutine classify(natom, out_length, old_n, new_n, old_list, new_list,&
    & n_increase, n_decrease, n_change, n_start, n_end, cn_out)
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer(4), intent(in) :: natom, out_length, old_n(natom), new_n(natom),&
            &old_list(natom,10), new_list(natom,10), n_start, n_end
    integer(4), intent(out) :: cn_out(out_length)
    ! integer(4), intent(in) :: old_list(_n), new_list(new_n)
    integer(4), intent(inout) :: n_increase, n_decrease, n_change
    integer(4) :: i, j
do j = 1, natom
    if (old_n(j) < new_n(j)) then
       n_increase = n_increase + 1
    else if (old_n(j) > new_n(j)) then
       n_decrease = n_decrease + 1
    else
        do i = 1, old_n(j)   !old_n == new_n at this moment
            if (old_list(j,i) /= new_list(j,i)) then
               n_change = n_change + 1
                exit
            endif
        enddo
    endif
    do i = n_start, n_end
      if (new_n(j) == i) then
        cn_out(i-n_start+1) = cn_out(i-n_start+1) + 1
        exit
      endif
    enddo
enddo
    return
end subroutine classify



end module topo_analysis
