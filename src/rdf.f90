module rdf
    use precision
    use data_types
    use parser, only: cutoff, pbc
    use data_input, only: coord_data, natom
    use neighbor_finder

    type histogram(bin_count)
        integer, len :: bin_count
        real(dp), dimension(bin_count, 2) :: bound = 0.
        real(dp), dimension(bin_count) :: pos = 0.
        integer, dimension(bin_count) :: n = 0
        real(dp), dimension(bin_count) :: rdf = 0.
    end type

    type(histogram(bin_count = :)), allocatable :: rdf_data

    contains

subroutine get_rdf()
    implicit none
    ! in:
    ! cutoff, xyz
    real(dp) :: bin_size, atom_density, r, d, ideal_count
    integer :: bin_count = 400

    integer :: xbin_max, ybin_max, zbin_max, atom, atom2
    integer :: xbin, ybin, zbin, id, checkid, p, q, o
    type(bins(capacity=:)), allocatable, dimension(:,:,:) :: mesh

    print *, "Radial distribution function calculation ... Start'"

    r = cutoff**2
    bin_size = cutoff/bin_count
    call get_bin_pos(cutoff, bin_count, rdf_data)

    atom_density = natom / (coord_data%lx * coord_data%ly * coord_data%lz)

    call create_bins(cutoff, mesh, xbin_max, ybin_max, zbin_max)

    do xbin = 1, xbin_max
    do ybin = 1, ybin_max
    do zbin = 1, zbin_max
        do atom = 1, mesh(xbin, ybin, zbin)%n
            id = mesh(xbin, ybin, zbin)%ids(atom)

            do p = -1, 1
            do q = -1, 1
            do o = -1, 1

            associate(checked_n => mesh(xbin+p, ybin+q, zbin+o)%n, checked_ids => mesh(xbin+p, ybin+q, zbin+o)%ids,&
            x_pbc => mesh(xbin+p, ybin+q, zbin+o)%x_pbc, y_pbc => mesh(xbin+p, ybin+q, zbin+o)%y_pbc, &
            z_pbc => mesh(xbin+p, ybin+q, zbin+o)%z_pbc, xyz => coord_data%coord)

                do atom2 = 1, checked_n
                    checkid = checked_ids(atom2)
                    if (checkid < id) then   !to avoid repeat calculation
                        d = (xyz(checkid,1) - x_pbc*coord_data%lx - xyz(id,1))**2& !x
                        + (xyz(checkid,2) - y_pbc*coord_data%ly - xyz(id,2))**2& !y
                        + (xyz(checkid,3) - z_pbc*coord_data%lz - xyz(id,3))**2  !z

                        if (d < r) then
                            call push_to_histbin(rdf_data, sqrt(d))
                        endif
                    endif
                end do

            end associate

            end do
            end do
            end do
        end do

    end do
    end do
    end do
    do p = 1, bin_count
        ideal_count = 2*pi*rdf_data%pos(p)**2*bin_size*atom_density*natom
        rdf_data%rdf(p) = rdf_data%n(p)/ideal_count
    end do

    call print_rdf()

end subroutine get_rdf

subroutine get_bin_pos(cutoff, bin_count, hist)
    implicit none
    ! Input:
    real(dp), intent(in) :: cutoff
    integer, intent(in) :: bin_count
    ! InOutput:
    type(histogram(bin_count = :)), allocatable, intent(inout) :: hist
    !
    real(dp) :: bin_size, half_bin
    integer :: i

    print *, 'Calculating bin positions in RDF histogram ...'

    allocate(histogram(bin_count = bin_count) :: hist)
    hist%n = 0
    hist%rdf = 0.

    bin_size = cutoff/bin_count
    half_bin = bin_size/2.0

    associate(bound => hist%bound, pos => hist%pos)
        do i = 1, bin_count
            pos(i) = bin_size * (i - 1) + half_bin  ! center
            bound(i, 1) = bin_size * (i - 1)        ! lower bound
            bound(i, 2) = bin_size * i              ! upper bound
        end do
    end associate

end subroutine get_bin_pos

subroutine push_to_histbin(hist, a)
    implicit none
    ! inout:
    type(histogram(bin_count = :)), allocatable, intent(inout) :: hist
    ! in:
    real(dp), intent(in) :: a
    !
    integer :: i

    do i = 1, hist%bin_count
        if (a < hist%bound(i,2)) then
            hist%n(i) = hist%n(i) + 1
            exit
        else
            continue
        end if
    end do

end subroutine push_to_histbin

subroutine print_rdf()
    implicit none
    integer :: i

    print *, '================================================'
    print *, 'r     g(r)'
    do i = 1, rdf_data%bin_count
        print *, rdf_data%pos(i), '  ',rdf_data%rdf(i)
    end do
    print *, '================================================'

end subroutine print_rdf

subroutine clean_rdf()
    implicit none

    if (allocated(rdf_data)) deallocate(rdf_data)
end subroutine clean_rdf

end module rdf
