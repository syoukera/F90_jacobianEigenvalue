module cema
    use, intrinsic :: iso_c_binding
    use globals, only : nf
    implicit none

    real(8), parameter :: t = 0.0d0
    real(8), parameter :: pres = 101325.0d0 ! ambient pressure
    character(len=10), allocatable :: species_names_cema(:)  ! species_names in c

    ! pyJac
    real(8), allocatable, target :: y(:)
    real(8), allocatable, target :: jac(:, :)
    real(8), allocatable :: wr(:), wi(:)   ! real and imaginary component of eigenvalue
    real(8), allocatable :: vl(:, :)       ! left eigenvecto
    real(8), allocatable :: vr(:, :)       ! right eigenvector
    real(8), allocatable :: work(:)        ! work array
    
    real(8), allocatable :: a_exp(:), b_exp(:)
    real(8), allocatable :: EP(:), EI(:)
    real(8) :: wr_max, EP_sum
        
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
    end interface

contains

    subroutine allocation_cema()
        implicit none
        
        allocate(y(nf))
        allocate(jac(nf, nf))
        allocate(wr(nf))
        allocate(wi(nf))
        allocate(vl(nf, nf))
        allocate(vr(nf, nf))
        allocate(work(4 * nf))
        allocate(a_exp(nf))
        allocate(b_exp(nf))
        allocate(EP(nf))
        allocate(EI(nf))
        allocate(species_names_cema(nf))

    end subroutine allocation_cema

    subroutine deallocation_cema()
        implicit none
        
        deallocate(y)
        deallocate(jac)
        deallocate(wr)
        deallocate(wi)
        deallocate(vl)
        deallocate(vr)
        deallocate(work)
        deallocate(a_exp)
        deallocate(b_exp)
        deallocate(EP)
        deallocate(EI)

    end subroutine deallocation_cema

    subroutine initialize_cema()
        implicit none
        
        y     = 0.0d0
        jac   = 0.0d0
        wr    = 0.0d0
        wi    = 0.0d0
        vl    = 0.0d0
        vr    = 0.0d0
        work  = 0.0d0
        a_exp = 0.0d0
        b_exp = 0.0d0
        EP    = 0.0d0
        EI    = 0.0d0

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
        integer :: i, ios, idx
    
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
    
        ! print *, "Read", idx+1, "species:"
        ! do i = 1, idx+1
        !     print *, i, trim(species_names_cema(i))
        ! end do

    end subroutine read_species_names_cema

    subroutine calc_cema(y_local,temp,lambda_e,index_EI)
        implicit none
        double precision, intent(in)::y_local(1:nf),temp
        double precision, intent(out)::lambda_e
        double precision, intent(out)::index_EI
        integer :: i, j, i_wr

        ! call allocation_cema()
        call initialize_cema()

        ! y = (/1.00000000d+03, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
        !       1.96270854d-01, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
        !       0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
        !       0.00000000d+00, 6.54236179d-02, 0.00000000d+00, 0.00000000d+00, &
        !       0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
        !       0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
        !       0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
        !       0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
            !   0.00000000d+00 /)

        ! Tamaoki mechanism
        y(1)  = temp        ! T
        y(2)  = y_local(3)  ! HE
        y(3)  = y_local(2)  ! AR
        y(4)  = y_local(8)  ! H2
        y(5)  = y_local(5)  ! O2
        y(6)  = y_local(4)  ! H
        y(7)  = y_local(6)  ! O
        y(8)  = y_local(7)  ! OH
        y(9)  = y_local(10)	! HO2
        y(10) = y_local(9)  ! H2O
        y(11) = y_local(11)	! H2O2
        y(12) = y_local(33)	! OHD-OH
        y(13) = y_local(18)	! N
        y(14) = y_local(12)	! NH3
        y(15) = y_local(13)	! NH2
        y(16) = y_local(14)	! NH
        y(17) = y_local(20)	! NNH
        y(18) = y_local(19)	! NO
        y(19) = y_local(23)	! N2O
        y(20) = y_local(15)	! HNO
        y(21) = y_local(17)	! HON
        y(22) = y_local(16)	! H2NO
        y(23) = y_local(25)	! HNOH
        y(24) = y_local(24)	! NH2OH
        y(25) = y_local(22)	! NO2
        y(26) = y_local(21)	! HONO
        y(27) = y_local(27)	! HNO2
        y(28) = y_local(28)	! NO3
        y(29) = y_local(26)	! HONO2
        y(30) = y_local(29)	! N2H2
        y(31) = y_local(30)	! H2NN
        y(32) = y_local(32)	! N2H4
        y(33) = y_local(31)	! N2H3
        ! y(34) = y_local(1)  ! N2
        
        ! calclate Jacobian using pyJac
        call eval_jacob(t, pres, c_loc(y), c_loc(jac))

        ! calclate eigenvalues and eigenvectors using dgeev in LAPACK
        call dgeev(jobvl, jobvr, n, jac, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)

        if (info /= 0) then
            print *, "Error: dgeev failed with INFO =", info
            stop 1
        end if

        ! get maximum eigenvalue
        i_wr = maxloc(wr, 1)
        ! wr_max = 0.0d0
        ! do i = 1, nf
        !     ! print "(E13.6)", wr(i)
        !     if (wr(i) > wr_max) then
        !         wr_max = wr(i)
        !         i_wr = i
        !     end if
        ! end do

        ! print *, "Maximum Eigenvalue ", i_wr, ":", wr(i_wr)
        lambda_e = wr(i_wr)

        ! calculate EP
        ! print *, "EP:"
        EP_sum = 0.0d0
        do j = 1, nf
            a_exp(j) = vr(j, i_wr)  ! 右固有ベクトル
            b_exp(j) = vl(j, i_wr)  ! 左固有ベクトル
            EP(j) = a_exp(j) * b_exp(j)
            EP_sum = EP_sum + abs(EP(j))
            ! print "(E13.6)", EP(j)
        end do

        ! calculate EI
        ! print *, "EI:"
        do j = 1, nf
            EI(j) = abs(EP(j)) / EP_sum
            ! print "(E13.6)", EI(j)
        end do

        ! get index of maximum EI
        index_EI = maxloc(EI, 1)

    end subroutine calc_cema

end module cema
