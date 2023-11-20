MODULE precision
    IMPLICIT NONE
    save
    integer, parameter :: dp = 4
!     integer, parameter :: dp = selected_real_kind(15, 307)
    integer :: ierr
    character(len=80) :: emsg

    contains

subroutine bubble_sort(last,array)
    implicit none
    ! last is the position of last element in array.
    integer(4), intent(inout) :: array(:) !length should match the main subroutine
    integer(4),intent(in) :: last
    integer(4) :: temp, i, j, k

    do i=last-1,1,-1
        do j=1,i
            if (array(j+1).lt.array(j)) then
                temp=array(j+1)
                array(j+1)=array(j)
                array(j)=temp
            endif
        enddo
    enddo
    return
end subroutine bubble_sort

END MODULE precision
