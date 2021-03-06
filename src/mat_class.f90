module mat_class
    
    use constants, only : wp
    use vec4_class

    !not fully implmented matix class
    !minimum implmented for neural sdf type
    type :: mat
        real(kind=wp) :: vals(4, 4)
        contains

        generic   :: operator(/)     => mat_div_scal
        generic   :: operator(*)     => mat_mult_scal, scal_mult_mat, mat_mult_mat
        generic   :: operator(+)     => mat_add_scal, scal_add_mat
        generic   :: operator(-)     => mat_minus_scal

        procedure, pass(a), private :: mat_div_scal

        procedure, pass(a), private :: mat_mult_mat
        procedure, pass(a), private :: mat_mult_scal
        procedure, pass(b), private :: scal_mult_mat

        procedure, pass(a), private :: mat_add_scal
        procedure, pass(b), private :: scal_add_mat

        procedure, pass(a), private :: mat_minus_scal

    end type mat


    interface mat
        module procedure mat_init
    end interface mat

private
public :: mat
contains

    
    type(mat) function mat_init(array)

        implicit none

        real(kind=wp) :: array(16)
        integer :: i, cnt

        cnt = 1

        do i = 1, 4
            mat_init%vals(:, i) = array(cnt:cnt+3)
            cnt = cnt + 4
        end do 


    end function mat_init


    type(mat) function mat_add_scal(a, b)

        implicit none

        class(mat),    intent(IN) :: a
        real(kind=wp), intent(IN) :: b

        mat_add_scal%vals = a%vals + b

    end function mat_add_scal


    type(mat) function scal_add_mat(a, b)

        implicit none

        class(mat),    intent(IN) :: b
        real(kind=wp), intent(IN) :: a

        scal_add_mat%vals = b%vals + a

    end function scal_add_mat


    type(mat) function mat_minus_scal(a, b)

        implicit none

        class(mat),    intent(IN) :: a
        real(kind=wp), intent(IN) :: b

        mat_minus_scal%vals = a%vals - b

    end function mat_minus_scal

    type(mat) function mat_div_scal(a, b)

        implicit none

        class(mat),    intent(IN) :: a
        real(kind=wp), intent(IN) :: b

        mat_div_scal%vals = a%vals / b

    end function mat_div_scal

    type(mat) function mat_mult_scal(a, b)

        implicit none

        class(mat),    intent(IN) :: a
        real(kind=wp), intent(IN) :: b

        mat_mult_scal%vals = a%vals * b

    end function mat_mult_scal

    type(mat) function scal_mult_mat(a, b)

        implicit none

        class(mat),    intent(IN) :: b
        real(kind=wp), intent(IN) :: a

        scal_mult_mat%vals = b%vals * a

    end function scal_mult_mat

    type(vec4) function mat_mult_mat(a, b)

        use vec4_class

        implicit none

        class(mat), intent(IN) :: a
        type(vec4), intent(IN) :: b

        real(kind=wp) :: tmp(4)

        tmp = matmul(a%vals, [b%x, b%y, b%z, b%p])
        mat_mult_mat = vec4(tmp(1), tmp(2), tmp(3), tmp(4))

    end function mat_mult_mat



end module mat_class

! Program p
    
!     use mat_class
!     use vec4_class

!     implicit none
    
!     real(kind=wp) :: array(16)
!     type(mat) :: m
!     type(vec4) :: v4

!     v4 = vec4(1., 1., 1., 1.)

!     array = [1., 1., 1., 1., 2., 2., 2., 2., 3., 3., 3., 3., 4., 4., 4., 4.]
!     m = mat(array)
!     write(*,"(4f9.5)")m%vals

!     m = m + 1.
!     print*," "
!     write(*,"(4f9.5)")m%vals

!     m = 1. + m
!     print*," "
!     write(*,"(4f9.5)")m%vals


!     m = m - 2.
!     print*," "
!     write(*,"(4f9.5)")m%vals

!     m = m / 2.
!     print*," "
!     write(*,"(4f9.5)")m%vals


!     m = m * 2.
!     print*," "
!     write(*,"(4f9.5)")m%vals


!     ! m = 2. * m
!     ! print*," "
!     ! write(*,"(4f9.5)")m%vals


!     v4 = m * v4
!     print*," "
!     write(*,"(4f9.5)")v4
! end program p