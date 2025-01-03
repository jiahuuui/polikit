module dynamic_data
    use precision
    use neighbor_finder
!     use parser
    implicit none

!     type int_linklist(length)
!         integer, len :: length
!         integer, dimension(length) :: list
!         type(int_linklist(length)), pointer :: next
!     end type int_linklist
!
!     type real_linklist(length)
!         integer, len :: length
!         real(dp), dimension(length) :: list
!         type(real_linklist(length)), pointer :: next
!     end type real_linklist
!
!     type(int_linklist(natom)), pointer :: head_i
!     type(int_linklist(natom)), pointer :: tail_i
!     type(int_linklist(natom)), pointer :: tmp_i
!
!     type(int_linklist(natom)), pointer :: head_r
!     type(int_linklist(natom)), pointer :: tail_r
!     type(int_linklist(natom)), pointer :: tmp_r

    integer, allocatable, dimension(:,:,:) :: neigh_list_bf
    integer, allocatable, dimension(:,:,:) :: poly_list_bf

    type(neighbor_list(atom_number = :, capacity = :)), allocatable, dimension(:) :: d_neighbors
!     type(polyhedron), allocatable, dimension(:) :: d_polys
!     type(length_data), allocatable, dimension(:) :: d_bondlengths
!     type(angle_data), allocatable, dimension(:) :: d_bondangles
    contains




end module dynamic_data
