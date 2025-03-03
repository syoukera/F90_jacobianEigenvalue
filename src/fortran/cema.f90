module cema
    use, intrinsic :: iso_c_binding
    implicit none

    ! 定数定義
    integer, parameter :: NSP = 33  ! NSP の値は適切に設定すること
    real(8), parameter :: t = 0.0d0
    real(8), parameter :: pres = 101325.0d0

    ! pyJac
    real(8), allocatable, target :: y(:)
    real(8), allocatable, target :: jac(:, :)
    real(8), allocatable :: wr(:), wi(:)   ! 固有値の実部と虚部
    real(8), allocatable :: vl(:, :)       ! 左固有ベクトル
    real(8), allocatable :: vr(:, :)       ! 右固有ベクトル
    real(8), allocatable :: work(:)      ! 作業配列
    
    real(8), allocatable :: a_exp(:), b_exp(:)
    real(8), allocatable :: EP(:), EI(:)
    real(8) :: wr_max, EP_sum
        
    ! LAPACK
    external dgeev
    character(1) :: jobvl = 'V'
    character(1) :: jobvr = 'V'
    integer, parameter :: n = NSP
    integer, parameter :: lda = NSP
    integer, parameter :: ldvl = NSP
    integer, parameter :: ldvr = NSP
    integer, parameter :: lwork = 4 * NSP
    integer :: info

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
        
        allocate(y(NSP))
        allocate(jac(NSP, NSP))
        allocate(wr(NSP))
        allocate(wi(NSP))
        allocate(vl(NSP, NSP))
        allocate(vr(NSP, NSP))
        allocate(work(4 * NSP))
        allocate(a_exp(NSP))
        allocate(b_exp(NSP))
        allocate(EP(NSP))
        allocate(EI(NSP))

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
        
        jobvl = 'V'
        jobvr = 'V'
        info = 1

    end subroutine initialize_cema

    subroutine calc_cema(y_local,temp,cema)
        use globals, only : nf
        implicit none
        double precision, intent(in)::y_local(1:nf),temp
        double precision, intent(out)::cema
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
        
        ! eval_jacob を呼び出し
        call eval_jacob(t, pres, c_loc(y), c_loc(jac))

        ! LAPACK の dgeev を使って固有値と固有ベクトルを計算
        call dgeev(jobvl, jobvr, n, jac, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)

        if (info /= 0) then
            print *, "Error: dgeev failed with INFO =", info
            stop 1
        end if

        ! 最大固有値を探す
        wr_max = 0.0d0
        i_wr = 2
        do i = 1, NSP
            ! print "(E13.6)", wr(i)
            if (wr(i) > wr_max) then
                wr_max = wr(i)
                i_wr = i
            end if
        end do

        ! print *, "Maximum Eigenvalue ", i_wr, ":", wr(i_wr)

        cema = wr(i_wr)
        ! cema = 1.0

        ! ! EP 計算
        ! print *, "EP:"
        ! EP_sum = 0.0d0
        ! do j = 1, NSP
        !     a_exp(j) = vr(j, i_wr)  ! 右固有ベクトル
        !     b_exp(j) = vl(j, i_wr)  ! 左固有ベクトル
        !     EP(j) = a_exp(j) * b_exp(j)
        !     EP_sum = EP_sum + abs(EP(j))
        !     ! print "(E13.6)", EP(j)
        ! end do

        ! ! EI 計算
        ! print *, "EI:"
        ! do j = 1, NSP
        !     EI(j) = abs(EP(j)) / EP_sum
        !     print "(E13.6)", EI(j)
        ! end do

        ! call deallocation_cema()

    end subroutine calc_cema

end module cema
