module tct
  use precision
  use io, n_atom => natom
  use neigh_finder, vectors => delta
  implicit none
  save
  real(dp) :: bond_length(n_atom,10,20),bond_angle(n_atom,30,20)
  integer(4) :: length_header(n_atom,10),angle_header(n_atom,30,2), counter(n_atom,10,2)
  real(dp) :: constrain(n_atom,2)


  contains

subroutine tct_calculate
! Constraint is calculated as follow:
! 1. For each atom, 2 lists are created for bond lengths and bond angles correspond
!   to it. It is a n*20 list for each atom, 20 is the number of frames will be
!   used to calculate standard deviation.
! 2. When having enough data, which means a bond is stable and unbroken for enough
!   time, standard deviation of its corralated angles are calculated.
! 3. If standard deviation is lower than given threshold, it is considered as an
!   active constraint.

  implicit none

  integer(4) :: i,j,k,l,topa,topb,tm,tm_list(10)
  real(dp) :: tct,angle,length
  real(dp), parameter :: pi = 3.141592653589793 ,r2d = 180.0/pi

  call neigh_finder_d

  do i =1, n_atom
    bond_length(i,:,2:20)=bond_length(i,:,:19)  !move bond_length
    bond_angle(i,:,2:20)=bond_angle(i,:,:19)  !move bond angle
    do j=1,10
      if (counter(i,j,2)/=20) then
        counter(i,j,1)=0
      endif
      if (counter(i,j,2)>0) then
        counter(i,j,2)=counter(i,j,2)-1
      endif
    enddo
    j = 1
    do !clear space for bond length
      if (counter(i,j,1)==0 .and. length_header(i,j) > 0) then
        bond_length(i,j:-1,:)=bond_length(i,j+1:,:)
        bond_length(i,10,:)=0
        l=1
        do !l = 1,20
          do k = 1,2
            !clear space for bond angle
            if (angle_header(i,l,k)==length_header(i,j)) then
              bond_angle(i,l:-1,:)=bond_angle(i,l+1:,:)
              bond_angle(i,30,:)=0
              angle_header(i,l:-1,:)=angle_header(i,l+1:,:)
              angle_header(i,30,:)=0
              if (l<29) then
              l=l-1
              endif
              exit
            endif
          enddo
          l=l+1
          ! if (l==30) then
          !   print *, i, j, l
          ! endif
          if (l>30) exit
        enddo
        length_header(i,j:-1)=length_header(i,j+1:)
        length_header(i,10)=0
        counter(i,j:-1,:)=counter(i,j+1:,:)
        counter(i,10,:)=0
        j=j-1
      end if
      j=j+1
      if (j>10) exit
    enddo

    do j = 1,n_neighbor(i)
      length = norm2(vectors(i,j,:))
      do l = 1,10
        if (neighbor(i,j)==length_header(i,l)) then
          bond_length(i,l,1)=length
          counter(i,l,1)=counter(i,l,1)+1
          counter(i,l,2)=20
          exit  !exit if this bond length is stored
        elseif (neighbor(i,j)<length_header(i,l)) then
          length_header(i,l+1:)=length_header(i,l:-1)
          length_header(i,l)=neighbor(i,j)
          bond_length(i,l+1:,:)=bond_length(i,l:-1,:)
          bond_length(i,l,1)=length
          counter(i,l+1:,:)=counter(i,l:-1,:)
          counter(i,l,1)=1
          counter(i,l,2)=20
          exit
        elseif (length_header(i,l)==0) then
          length_header(i,l+1:)=length_header(i,l:-1)
          length_header(i,l)=neighbor(i,j)
          bond_length(i,l,1)=length
          counter(i,l,1)=counter(i,l,1)+1
          counter(i,l,2)=20
          exit
        endif
        if (l==10) then
          print *, 'Can"t find a place to insert bond length data.'
        endif
      enddo
      do l = j+1,n_neighbor(i)
        topa=neighbor(i,j)
        topb=neighbor(i,l)
        angle= r2d*acos((vectors(i,j,1)*vectors(i,l,1)+vectors(i,j,2)*&
        &vectors(i,l,2)+vectors(i,j,3)*vectors(i,l,3))/(norm2(vectors(i,j,:))*&
        &norm2(vectors(i,l,:))))

        do k = 1,30
          if (topa==angle_header(i,k,1) .and. topb==angle_header(i,k,2)) then
            bond_angle(i,k+1:,:)=bond_angle(i,k:-1,:)
            bond_angle(i,k,1)=angle
            exit
          elseif (topa<angle_header(i,k,1)) then
            angle_header(i,k+1:,:)=angle_header(i,k:-1,:)
            angle_header(i,k,1)=topa
            angle_header(i,k,2)=topb
            bond_angle(i,k+1:,:)=bond_angle(i,k:-1,:)
            bond_angle(i,k,1)=angle
            exit
          elseif (topa==angle_header(i,k,1) .and. topb<angle_header(i,k,2)) then
            angle_header(i,k+1:,:)=angle_header(i,k:-1,:)
            angle_header(i,k,1)=topa
            angle_header(i,k,2)=topb
            bond_angle(i,k+1:,:)=bond_angle(i,k:-1,:)
            bond_angle(i,k,1)=angle
            exit
          elseif (topa>angle_header(i,k,1) .and. k>1 .and. topa<=angle_header(i,k-1,1)) then
            angle_header(i,k+1:,:)=angle_header(i,k:-1,:)
            angle_header(i,k,1)=topa
            angle_header(i,k,2)=topb
            bond_angle(i,k+1:,:)=bond_angle(i,k:-1,:)
            bond_angle(i,k,1)=angle
            exit
          elseif (angle_header(i,k,1)==0) then
            angle_header(i,k+1:,:)=angle_header(i,k:-1,:)
            angle_header(i,k,1)=topa
            angle_header(i,k,2)=topb
            bond_angle(i,k+1:,:)=bond_angle(i,k:-1,:)
            bond_angle(i,k,1)=angle
            exit
          endif
          if (k==30) then
            print *, 'Can"t find a place to insert bond angle data.'
          endif
        enddo
      enddo
    enddo

    tm=0
    tm_list=0
    do j=1,10
      if (counter(i,j,1)>=20 .and. counter(i,j,2)==20) then
        tm=tm+1
        tm_list(tm)=length_header(i,j)
        !stddva of bond length
        tct=stddev(bond_length(i,j,:))
        if (tct<1.) then
          constrain(i,1) = constrain(i,1)+1.
        endif
        do l=j+1,10
          if (counter(i,l,1)>=20 .and. counter(i,l,2)==20) then
            do k=1,30
              if (angle_header(i,k,1)==length_header(i,j) .and. &
              &angle_header(i,k,2)==length_header(i,l)) then
                tct = stddev(bond_angle(i,k,:))
                if (tct<15.) then
                  constrain(i,2) = constrain(i,2)+1.
                endif
                exit
              endif
            enddo
          endif
        enddo
      endif
    enddo
    constrain(i,1)=0.5*constrain(i,1)
    if (constrain(i,2)>=1.) then
      constrain(i,2)=angle2ncbb(constrain(i,2))
      ! sqrt(8.*constrain(i,2)+1.)-2.
    endif
  enddo

  ! print *, bond_length(886315,2,:)
  ! print *, vectors(886315,2,:)

  return
  contains
  function stddev(array)
    ! a function to calculate standard deviation from a 1D array.
    implicit none
    real(dp), intent(in) :: array(:)
    real(dp) :: bar, sm, stddev
    integer(4) :: p, q
    sm=0.
    q=size(array)
    bar = sum(array)/q
    do p = 1,q
      sm = sm+(array(p)-bar)**2.
    enddo
    stddev=sqrt(sm/q)
    return
  end function stddev

  function angle2ncbb(angle_number)
    ! a step function to get ncbb from stable angle number, can be a little bit ambigious.
    implicit none
    real(dp), intent(in) :: angle_number
    real(dp) :: angle2ncbb
    if (angle_number < 1.) then
      angle2ncbb = 0.
    else if (angle_number < 3.) then
      angle2ncbb = 1.
    else if (angle_number < 6.) then
      angle2ncbb = 3.
    else if (angle_number < 10.) then
      angle2ncbb = 5.
    else if (angle_number < 15.) then
      angle2ncbb = 7.
    else if (angle_number < 21.) then
      angle2ncbb = 9.
    else if (angle_number < 28.) then
      angle2ncbb = 11.
    else
      angle2ncbb = 13.
    end if
    return
  end function angle2ncbb

end subroutine tct_calculate


subroutine poly_volume(n_atom,n_neigh,delta,ptype,ctype,n_ctype,hull_v,std_volume)
  implicit none
  integer(4), intent(in) :: n_atom,n_neigh(n_atom),ptype(n_atom),ctype,n_ctype
  real(dp), intent(in) :: delta(n_atom,10,3),hull_v(n_ctype)
  real(dp), intent(out) :: std_volume(n_ctype)
  real(dp) :: edge(n_atom,20), vfactor,mean
  integer(4) :: n,k,i,j,c
  vfactor=6.*sqrt(2.)
  c=1
  do n=1,n_atom
    if (ptype(n)==ctype .and. n_neigh(n)==4) then
      k=1
      do i=1,n_neigh(n)-1
        do j=i+1,n_neigh(n)
          edge(n,k)=norm2((delta(n,j,:)-delta(n,i,:)))
          k=k+1
        enddo
      enddo
      mean = sum(edge(n,:k-1))/(k-1)
      std_volume(c)=hull_v(c)/mean**3*vfactor
      c=c+1
    endif
  enddo
end subroutine poly_volume

end module tct
