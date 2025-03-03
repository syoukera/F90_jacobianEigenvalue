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
    integer :: lwork, info
    real(8), allocatable :: a_exp(:), b_exp(:)
    real(8), allocatable :: EP(:), EI(:)
    real(8) :: wr_max, EP_sum
    integer :: i, j, i_wr

    ! LAPACK
    external dgeev
    character(1) :: jobvl, jobvr
    integer :: n, lda, ldvl, ldvr

    interface
        subroutine eval_jacob(t, pres, y, jac) bind(C)
            use, intrinsic :: iso_c_binding
            real(c_double), value :: t, pres
            type(c_ptr), value :: y
            type(c_ptr), value :: jac
        end subroutine eval_jacob
    end interface

contains

    subroutine init_cema()
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

        y = (/1.00000000d+03, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
              1.96270854d-01, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
              0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
              0.00000000d+00, 6.54236179d-02, 0.00000000d+00, 0.00000000d+00, &
              0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
              0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
              0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
              0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
              0.00000000d+00 /)

        ! y(:) = (/1000.d0, 0.d0, 0.04030952d0, 0.d0, 0.95969048d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)

    end subroutine init_cema

    subroutine calc_cema()
        implicit none

        ! 初期化
        jobvl = 'V'
        jobvr = 'V'
        n = NSP
        lda = NSP
        ldvl = NSP
        ldvr = NSP
        lwork = 4 * NSP

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
            if (wr(i) > wr_max) then
                wr_max = wr(i)
                i_wr = i
            end if
        end do

        print *, "Maximum Eigenvalue ", i_wr, ":", wr(i_wr)

        ! EP 計算
        print *, "EP:"
        EP_sum = 0.0d0
        do j = 1, NSP
            a_exp(j) = vr(j, i_wr)  ! 右固有ベクトル
            b_exp(j) = vl(j, i_wr)  ! 左固有ベクトル
            EP(j) = a_exp(j) * b_exp(j)
            EP_sum = EP_sum + abs(EP(j))
            print "(E13.6)", EP(j)
        end do

        ! EI 計算
        print *, "EI:"
        do j = 1, NSP
            EI(j) = abs(EP(j)) / EP_sum
            print "(E13.6)", EI(j)
        end do

    end subroutine calc_cema

end module cema
