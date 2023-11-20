!Delaunay triangulation
MODULE delaunay
    use precision
    use io
    implicit none

    type del2d
        integer ::

    end type

contains

    SUBROUTINE delaunay2d
        implicit NONE

        call get_xyz

!         construct a big enough triangle
        x_min=minval(xy(:,1))
        y_min=minval(xy(:,2))

        x_max=maxval(xy(:,1))
        y_max=maxval(xy(:,2))

        do i=1,natom
!             insert points and divide triangle
        end do

!         deal with boundary condition

        return
    END SUBROUTINE

END MODULE delaunay
