module rdf
    use precision
    use data_types
    use parser, only: cutoffs, pbc
    use data_input, only: coord_data, natom, ntype, atom_frac
    use neighbor_finder

    type histogram(bin_count)
        integer, len :: bin_count
        real(dp), dimension(bin_count, 2) :: bound = 0.
        real(dp), dimension(bin_count) :: pos = 0.
        integer, dimension(bin_count) :: n = 0
        real(dp), dimension(bin_count) :: rdf = 0.
    end type

    real(dp), dimension(:,:), allocatable :: bin_info
    ! [lower, center, upper][bin_count]
    integer, dimension(:,:,:), allocatable :: rdf_raw
    ! [type1][type2][bin_count]
    real(dp), dimension(:,:), allocatable :: rdf_data

    contains

subroutine get_rdf()
    implicit none
    ! in:
    ! cutoff, xyz
    real(dp) :: bin_size, atom_density, r, d, cutoff
    integer :: bin_count = 400

    integer :: xbin_max, ybin_max, zbin_max, atom, atom2
    integer :: xbin, ybin, zbin, p, q, o, k
    integer :: id, checkid, type1, type2

    real(dp), allocatable :: ideal_count(:)
    integer, allocatable :: sum_rdf(:)

    print *, "Radial distribution function calculation ... Start'"

    cutoff = cutoffs(1)
    r = cutoff**2
    bin_size = cutoff/bin_count
    call get_bin_pos(cutoff, bin_count) ! initialize bin_info, rdf_raw, rdf_data

    atom_density = natom / (coord_data%lx * coord_data%ly * coord_data%lz)

    call create_bins(cutoff, xbin_max, ybin_max, zbin_max)
!$omp parallel do
    do xbin = 1, xbin_max
!     PRINT*, "Hello from thread", xbin, OMP_GET_THREAD_NUM()
    do ybin = 1, ybin_max
    do zbin = 1, zbin_max
        do atom = 1, cells_n(xbin, ybin, zbin)
            id = cells_ids(xbin, ybin, zbin, atom)
            type1 = coord_data%ptype(id)

            do p = -1, 1
            do q = -1, 1
            do o = -1, 1

            associate(checked_n => cells_n(xbin+p, ybin+q, zbin+o), checked_ids => cells_ids(xbin+p, ybin+q, zbin+o, :),&
            x_pbc => cells_xpbc(xbin+p, ybin+q, zbin+o), y_pbc => cells_ypbc(xbin+p, ybin+q, zbin+o), &
            z_pbc => cells_zpbc(xbin+p, ybin+q, zbin+o), xyz => coord_data%coord)

                do atom2 = 1, checked_n
                    checkid = checked_ids(atom2)
                    type2 = coord_data%ptype(checkid)

                    if (checkid < id) then   !to avoid repeat calculation

                        d = (xyz(checkid,1) - x_pbc*coord_data%lx - xyz(id,1))**2& !x
                        + (xyz(checkid,2) - y_pbc*coord_data%ly - xyz(id,2))**2& !y
                        + (xyz(checkid,3) - z_pbc*coord_data%lz - xyz(id,3))**2  !z

                        if (d < r) then
                            call push_to_histbin(type1, type2, sqrt(d))
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
!$omp end parallel do

    k = 1+ntype*(ntype+1)/2

    allocate(rdf_data(k, bin_count))
    allocate(ideal_count(bin_count))
    allocate(sum_rdf(bin_count))

    ideal_count = 4*pi*bin_info(1,:)**2*bin_size*atom_density*natom
    sum_rdf = [(sum(rdf_raw(:,:,o)), o=1, bin_count)]
!     sum(rdf_raw, dim=3)
    rdf_data(1,:) = sum_rdf/ideal_count

    o = 2
    do p = 1, ntype
        do q = p, ntype
            ideal_count = 4*pi*bin_info(1,:)**2*bin_size*atom_density*natom*atom_frac(p)*atom_frac(q)
            rdf_data(o,:) = rdf_raw(p,q,:)/ideal_count
            o = o+1
        end do
    end do

end subroutine get_rdf

subroutine calculate_rdf()
    implicit none

    call get_rdf()
    call print_rdf()

end subroutine calculate_rdf

! Calculate the Wendt-Abraham parameter from RDF data.
subroutine wa_parameter()
    implicit none
    integer, dimension(1) :: p
    integer :: minp, maxp
    real(dp) :: gmin, gmax
    real(dp) :: r_wa    ! Wendt-Abraham parameter
    real(dp) :: r_mwa   ! Modified Wendt-Abraham parameter

    call get_rdf()
    p = maxloc(rdf_data(1,:))
    maxp = p(1)
    gmax = rdf_data(1,maxp)

    p = minloc(rdf_data(1,maxp:))
    minp = p(1)+maxp
    gmin = rdf_data(1,minp)

    r_wa = gmin/gmax
    r_mwa = r_wa**2

    print *, 'Minloc and Maxloc:', minp, maxp
    print *, 'Min and Max value:', gmin, gmax
    print *, 'R_wa (Wendt-Abraham parameter) is: ', r_wa
    print *, 'R_mwa (Modified WA parameter) is: ', r_mwa

end subroutine wa_parameter

! initialize histogram, rdf_raw, rdf_data.
subroutine get_bin_pos(cutoff, bin_count)
    implicit none
    ! Input:
    real(dp), intent(in) :: cutoff
    integer, intent(in) :: bin_count
    ! Private:
    real(dp) :: bin_size, half_bin
    integer :: i

    if (.not. allocated(bin_info)) allocate(bin_info(3, bin_count))
    if (.not. allocated(rdf_raw)) allocate(rdf_raw(ntype, ntype, bin_count))
!     if (.not. allocated(rdf_data)) allocate(rdf_data(ntype, ntype, bin_count))

    bin_info = 0.
    rdf_raw = 0
!     rdf_data = 0.

    bin_size = cutoff/bin_count
    half_bin = bin_size/2.0

    do i = 1, bin_count
        bin_info(2,i) = bin_size*(i-1)+half_bin ! center
        bin_info(1,i) = bin_size*(i-1)          ! lower bound
        bin_info(3,i) = bin_size*i              ! upper bound
    end do

end subroutine get_bin_pos

! Push raw data to the histogram by pair types.
subroutine push_to_histbin(type1, type2, a)
    implicit none
    ! in:
    integer, intent(in) :: type1, type2
    real(dp), intent(in) :: a
    ! private
    integer :: i, bc

    bc = size(bin_info(1,:))

    do i = 1, bc
        if (a < bin_info(3,i)) then
            rdf_raw(type1, type2, i) = rdf_raw(type1, type2, i)+1
            rdf_raw(type2, type1, i) = rdf_raw(type2, type1, i)+1
            exit
        else
            continue
        end if
    end do

end subroutine push_to_histbin

subroutine print_rdf()
    implicit none
    integer :: i, j, t
    character(len=8) :: str
    character(:), allocatable :: head
    real(dp), allocatable :: row(:)

    allocate(row(2+ntype*(ntype+1)/2))
    head = 'r     g(r)'

    print *, '================================================'
    do i = 1, ntype
        do j = i, ntype
            t = i*10+j
            write (str, '(I8)') t
            head  = head//str
        end do
    end do

    print *, head
    do i = 1, size(bin_info(1,:))
        row(1) = bin_info(1,i)
        row(2:) = rdf_data(:,i)
        print *, row
    end do
    print *, '================================================'

end subroutine print_rdf

subroutine clean_rdf()
    implicit none

    if (allocated(bin_info)) deallocate(bin_info)
    if (allocated(rdf_raw)) deallocate(rdf_raw)
    if (allocated(rdf_data)) deallocate(rdf_data)
end subroutine clean_rdf

end module rdf
