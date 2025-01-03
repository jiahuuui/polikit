MODULE data_types
USE precision

IMPLICIT NONE

! linked list
type intll(na)
    integer(inp), len :: na
    integer(inp), dimension(na) :: list
    type(intll(na)), pointer :: next
end type intll

type reall(na)
    integer(inp), len :: na
    real(dp), dimension(na) :: list
    type(reall(na)), pointer :: next
end type reall

INTERFACE print_array
    MODULE PROCEDURE :: printa, printl!, printr
END INTERFACE

contains
! Print integer 2-D array.
SUBROUTINE printa(array)
    IMPLICIT NONE
    integer(inp), allocatable, intent(in) :: array(:,:)
    integer :: i
    do i = 1, size(array(:,1))
        print *, array(i,:)
    end do
    print *, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
END SUBROUTINE printa

! Print logical 2-D array.
SUBROUTINE printl(array)
    IMPLICIT NONE
    logical, allocatable, intent(in) :: array(:,:)
    integer :: i
    do i = 1, size(array(:,1))
        print *, array(i,:)
    end do
    print *, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
END SUBROUTINE printl

! ! Print rings.
! SUBROUTINE printr(ringlist)
!     IMPLICIT NONE
!     type(ring), allocatable, intent(in) :: ringlist(:)
!     integer :: i, l
!     do i = 1, size(ringlist(:))
!         l = ringlist(i)%l
!         print *, i, '(ID) ', ringlist(i)%l, '(RSIZE) ', ringlist(i)%element(:l)
!     end do
!     print *, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
! END SUBROUTINE printr


!     subroutine push_back_int()
!         implicit none
!
!         allocate(tail%next, stat=ierr)
!         tail => tail%next
!         nullify(tail%next)
!         tail%list = constrain(:,1)
!         if (size(tail%list) /= ntemp) then
!             print *, '!Error: Unmatched atom number between frames!'
!             stop
!         end if
!         ntemp = natom
!
!     end subroutine push_back_int
!
!     subroutine push_back_real()
!         implicit none
!
!         allocate(tail%next, stat=ierr)
!         tail => tail%next
!         nullify(tail%next)
!         tail%list = constrain(:,1)
!         if (size(tail%list) /= ntemp) then
!             print *, '!Error: Unmatched atom number between frames!'
!             stop
!         end if
!         ntemp = natom
!
!     end subroutine push_back_real

END MODULE data_types
