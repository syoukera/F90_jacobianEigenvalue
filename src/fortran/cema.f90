module cema
    use, intrinsic :: iso_c_binding
    use globals, only : nf, nrf, nrb, nrp
    implicit none

    character(len=10), allocatable :: species_names_cema(:)  ! species_names in c

    ! pyJac
    double precision, allocatable, target :: y(:)
    double precision, allocatable, target :: jac(:, :)
    double precision, allocatable, target :: conc(:)
    double precision, allocatable, target :: fwd_rxn_rates(:)
    double precision, allocatable, target :: rev_rxn_rates(:)
    double precision, allocatable, target :: pres_mod(:)
    double precision, allocatable :: rop(:) ! rate ob progress
    double precision, target :: y_N
    double precision, target :: mw_avg ! mass-averaged density
    double precision, target :: rho ! average molecular weight

    double precision, allocatable :: wr(:), wi(:)   ! real and imaginary component of eigenvalue
    double precision, allocatable :: vl(:, :)       ! left eigenvecto
    double precision, allocatable :: vr(:, :)       ! right eigenvector
    double precision, allocatable :: work(:)        ! work array
    
    double precision, allocatable :: a_exp(:), b_exp(:)
    double precision, allocatable :: EP(:), EI(:), PP(:), PI(:)
    double precision :: wr_max, EP_sum, PP_sum

    ! indices for reaction picked from spec_rates.c
    integer(4), allocatable :: stoich_coeffs(:, :)
    integer(4), allocatable :: list_i_rev_rates(:)
    integer(4), allocatable :: list_k_pres_mod(:)
    
    ! indices for sort
    integer(4), allocatable :: indices(:)
    double precision, parameter :: threshold_negative = 1d0
    integer(4), parameter :: k_negative = 7 ! n_element + 1
        
    ! LAPACK
    external dgeev
    character(1) :: jobvl, jobvr
    integer :: n, lda, ldvl, ldvr, lwork, info

    interface

        subroutine eval_jacob(t, pres, y, jac) bind(C)
            use, intrinsic :: iso_c_binding
            real(c_double), value :: t, pres
            type(c_ptr), value :: y
            type(c_ptr), value :: jac
        end subroutine eval_jacob

        subroutine eval_conc(T, pres, y, y_N, mw_avg, rho, conc) bind(C)
            use, intrinsic :: iso_c_binding
            real(c_double), value :: T, pres
            type(c_ptr), value :: y
            type(c_ptr), value :: y_N
            type(c_ptr), value :: mw_avg
            type(c_ptr), value :: rho
            type(c_ptr), value :: conc
        end subroutine eval_conc

        subroutine eval_rxn_rates(T, pres, conc, fwd_rxn_rates, rev_rxn_rates) bind(C)
            use, intrinsic :: iso_c_binding
            real(c_double), value :: T, pres
            type(c_ptr), value :: conc
            type(c_ptr), value :: fwd_rxn_rates
            type(c_ptr), value :: rev_rxn_rates
        end subroutine eval_rxn_rates
        
        subroutine get_rxn_pres_mod(T, pres, conc, pres_mod) bind(C)
            use, intrinsic :: iso_c_binding
            real(c_double), value :: T, pres
            type(c_ptr), value :: conc
            type(c_ptr), value :: pres_mod
        end subroutine get_rxn_pres_mod

    end interface

contains

    subroutine allocation_cema()
        implicit none
        
        ! pyJac
        allocate(y(nf))
        allocate(jac(nf, nf))
        allocate(conc(nf))
        allocate(fwd_rxn_rates(nrf))
        allocate(rev_rxn_rates(nrb))
        allocate(pres_mod(nrp))
        allocate(rop(nrf))

        allocate(wr(nf))
        allocate(wi(nf))
        allocate(vl(nf, nf))
        allocate(vr(nf, nf))
        allocate(work(4 * nf))
        allocate(a_exp(nf))
        allocate(b_exp(nf))
        allocate(EP(nf))
        allocate(EI(nf))
        allocate(PP(nrf))
        allocate(PI(nrf))
        allocate(species_names_cema(nf))

        ! indices for reaction picked from spec_rates.c
        allocate(stoich_coeffs(nf, nrf))
        allocate(list_i_rev_rates(nrf))
        allocate(list_k_pres_mod(nrf))

        allocate(indices(nf))

    end subroutine allocation_cema

    subroutine deallocation_cema()
        implicit none
        
        deallocate(y)
        deallocate(jac)
        deallocate(conc)
        deallocate(fwd_rxn_rates)
        deallocate(rev_rxn_rates)
        deallocate(pres_mod)
        deallocate(rop)

        deallocate(wr)
        deallocate(wi)
        deallocate(vl)
        deallocate(vr)
        deallocate(work)
        deallocate(a_exp)
        deallocate(b_exp)
        deallocate(EP)
        deallocate(EI)

        ! indices for reaction picked from spec_rates.c
        deallocate(stoich_coeffs)
        deallocate(list_i_rev_rates)
        deallocate(list_k_pres_mod)

        deallocate(indices)

    end subroutine deallocation_cema

    subroutine initialize_cema()
        implicit none
        
        y     = 0.0d0
        jac   = 0.0d0
        fwd_rxn_rates = 0.0d0
        rev_rxn_rates = 0.0d0
        pres_mod = 0.0d0
        rop = 0.0d0

        wr    = 0.0d0
        wi    = 0.0d0
        vl    = 0.0d0
        vr    = 0.0d0
        work  = 0.0d0
        a_exp = 0.0d0
        b_exp = 0.0d0

        EP = 0.0d0
        EI = 0.0d0
        PP = 0.0d0
        PI = 0.0d0
        
        ! ! indices for reaction picked from spec_rates.c
        ! stoich_coeffs = 0
        ! list_i_rev_rates = 0
        ! list_k_pres_mod = 0

        ! LAPACK
        jobvl = 'V'
        jobvr = 'V'
        n = nf
        lda = nf
        ldvl = nf
        ldvr = nf
        lwork = 4 * nf
        info = 1

    end subroutine initialize_cema

    subroutine read_species_names_cema()
        implicit none
        character(len=10) :: dummy_name  ! temporary character for name
        character(len=100) :: line  ! temporary character for line
        integer :: ios, idx
    
        open(unit=10, file="src/c/mechanism.h", status="old", action="read")
    
        ! skip header
        do
            read(10, '(A)', iostat=ios) line
            if (ios /= 0) exit  ! exit at end of file
    
            ! skip until "/* Species Indexes"]
            if (index(line, "/* Species Indexes") > 0) then
                exit
            end if
        end do
    
        ! read species names
        do
            read(10, '(A)', iostat=ios) line
            if (ios /= 0) exit  ! exit at end of file
            if (index(line, "*/") > 0) exit  ! exit at end of species
    
            ! print *, line
    
            read(line, '(I3, A)', iostat=ios) idx, dummy_name
            dummy_name = trim(adjustl(dummy_name))
            species_names_cema(idx + 1) = dummy_name
        end do
    
        close(10)

    end subroutine read_species_names_cema

    subroutine read_indices_cema()
        implicit none
        integer :: idf

        open(newunit=idf, file="ipynb/stoich_coeffs.bin", form="unformatted", access="stream", action="read")
        read(idf) stoich_coeffs
        close(idf)
        
        open(newunit=idf, file="ipynb/list_i_rev_rates.bin", form="unformatted", access="stream", action="read")
        read(idf) list_i_rev_rates
        close(idf)
        
        open(newunit=idf, file="ipynb/list_k_pres_mod.bin", form="unformatted", access="stream", action="read")
        read(idf) list_k_pres_mod
        close(idf)
    
    end subroutine read_indices_cema

    subroutine calc_cema(y_local,p_local,t_local,lambda_e,index_EI,index_PI)
        implicit none
        double precision, intent(in)::y_local(1:nf)
        double precision, intent(in)::p_local
        double precision, intent(in)::t_local
        double precision, intent(out)::lambda_e
        double precision, intent(out)::index_EI
        double precision, intent(out)::index_PI

        double precision, parameter :: t = 0.0d0 ! dummy parameter for time

        integer :: i, j, i_wr

        ! call allocation_cema()
        call initialize_cema()

        ! Tamaoki mechanism
        y(1)  = t_local     ! T
        y(2)  = y_local(3)  ! HE
        y(3)  = y_local(2)  ! AR
        y(4)  = y_local(8)  ! H2
        y(5)  = y_local(5)  ! O2
        y(6)  = y_local(4)  ! H
        y(7)  = y_local(6)  ! O
        y(8)  = y_local(7)  ! OH
        y(9)  = y_local(10) ! HO2
        y(10) = y_local(9)  ! H2O
        y(11) = y_local(11) ! H2O2
        y(12) = y_local(33) ! OHD-OH
        y(13) = y_local(18) ! N
        y(14) = y_local(12) ! NH3
        y(15) = y_local(13) ! NH2
        y(16) = y_local(14) ! NH
        y(17) = y_local(20) ! NNH
        y(18) = y_local(19) ! NO
        y(19) = y_local(23) ! N2O
        y(20) = y_local(15) ! HNO
        y(21) = y_local(17) ! HON
        y(22) = y_local(16) ! H2NO
        y(23) = y_local(25) ! HNOH
        y(24) = y_local(24) ! NH2OH
        y(25) = y_local(22) ! NO2
        y(26) = y_local(21) ! HONO
        y(27) = y_local(27) ! HNO2
        y(28) = y_local(28) ! NO3
        y(29) = y_local(26) ! HONO2
        y(30) = y_local(29) ! N2H2
        y(31) = y_local(30) ! H2NN
        y(32) = y_local(32) ! N2H4
        y(33) = y_local(31) ! N2H3
        ! y(34) = y_local(1)  ! N2

        
        ! calclate Jacobian using pyJac
        call eval_jacob(t, p_local, c_loc(y), c_loc(jac))

        ! calclate eigenvalues and eigenvectors using dgeev in LAPACK
        call dgeev(jobvl, jobvr, n, jac, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)

        if (info /= 0) then
            print *, "Error: dgeev failed with INFO =", info
            stop 1
        end if

        ! get maximum eigenvalue
        i_wr = maxloc(wr, 1)
        lambda_e = wr(i_wr)

        ! ! pick negative eigenvalue when maximum eigenvalue is small
        ! if (lambda_e < threshold_negative) then
            
        !     ! get sorted indices for wr
        !     call sort_indices_by_abs()

        !     ! update eigenvalue and index
        !     i_wr = indices(k_negative)
        !     lambda_e = wr(i_wr)

        ! end if
        
        ! calculate EP
        EP_sum = 0.0d0
        do j = 1, nf
            a_exp(j) = vr(j, i_wr)  ! 右固有ベクトル
            b_exp(j) = vl(j, i_wr)  ! 左固有ベクトル
            EP(j) = a_exp(j) * b_exp(j)
            EP_sum = EP_sum + abs(EP(j))
        end do

        ! calculate EI
        do j = 1, nf
            EI(j) = abs(EP(j)) / EP_sum
        end do

        ! get index of maximum EI
        index_EI = maxloc(EI, 1)
        ! test to call pyJac function
        call eval_conc(t_local, p_local, c_loc(y(2)), c_loc(y_N), c_loc(mw_avg), c_loc(rho), c_loc(conc))
        

        call eval_rxn_rates(t_local, p_local, c_loc(conc), c_loc(fwd_rxn_rates), c_loc(rev_rxn_rates))
        call get_rxn_pres_mod(t_local, p_local, c_loc(conc), c_loc(pres_mod))
        
        ! calculate rop
        do i = 1, nrf

            ! rop from forward/backword reaction
            if (list_i_rev_rates(i) /= -1) then        
                rop(i) = fwd_rxn_rates(i) - rev_rxn_rates(list_i_rev_rates(i)+1) ! adjust 0 strat index
            else 
                rop(i) = fwd_rxn_rates(i)
            end if

            ! add pressure module to rop
            if (list_k_pres_mod(i) /= -1) then
                rop(i) = rop(i)*pres_mod(list_k_pres_mod(i)+1) ! adjust 0 strat index
            end if

        end do

        ! calulate Participation Pointer (PP)
        PP = 0.0d0
        do i = 1, nrf
            ! sum product along species
            do j = 1, nf-1
                PP(i) = PP(i) + a_exp(j+1)*dble(stoich_coeffs(j, i))
            end do 

            ! multiply rop
            PP(i) = PP(i) * rop(i)

            ! sum of PP
            PP_sum = PP_sum + abs(PP(i))

        end do

        ! calculate PI
        do i = 1, nrf
            PI(i) = abs(PP(i)) / PP_sum
        end do

        ! get index of maximum PI
        index_PI = maxloc(PI, 1)

    end subroutine calc_cema

    ! sort indices of wr by abs
    subroutine sort_indices_by_abs()

        integer :: i, j, temp_idx
        
        ! initialize indices
        indices = [(i, i = 1, nf)]

        do i = 1, nf-1
            do j = i+1, nf
                if (abs(wr(indices(i))) < abs(wr(indices(j)))) then
                    temp_idx = indices(i)
                    indices(i) = indices(j)
                    indices(j) = temp_idx
                end if
            end do
        end do

    end subroutine sort_indices_by_abs

end module cema
