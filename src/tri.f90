!Delaunay triangulation
MODULE tri
    use precision
    use io
    implicit none

    type alpha
        integer(inp) :: id
        integer(inp) :: vertex(4) ! id of its four vertex atoms
        real(dp) :: cenxyz(3)   ! coordination of circumscribed circle center
        real(dp) :: circumr     ! radius of circumscribed circle
        integer(inp)) :: face(4,2)    ! four faces, 2 means id (the atom id opposite to the face) and neighboring alpha id
    end type

contains

SUBROUTINE delaunay
    implicit NONE

    ! if lx ly lz > 10 angstrom, use 10 as extending distance, other wise use lx ly lz.
    rext =

    xyz(x<rext)
    xyz(y<rext)
    xyz(z<rext)
    xyz(x,y<rext)
    xyz(y,z<rext)
    xyz(x,z<rext)
    xyz(x,y,z<rext)

    ! build super-triangle

    a%vertex = (,,,)
    a%cenxyz =
    a%circumr =
    a%face = (,,,;,)

    do i = 1, next
        do j = 1, ntri
            if (xyz(i) -> a(j)%cenxyz) < a(j)%circumr then
                polylist.append(j)
            end if
        end do
        do j in polylist
            do k = 1,4
                if face(k,2) not in polylist then
                    ! remove old alphas
                    ! construct new alpha with face vertex and xyz(i)
                end if
            end do
        end do
    end do
    ! only store alphas that within or partly within the real box
    do i = 1, ntri

    end do

END SUBROUTINE delaunay

END MODULE tri
