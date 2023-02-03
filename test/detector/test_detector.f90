module testsDetectorMod

    use detector_mod
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use constants, only : wp

    implicit none

    private
    public :: detector_suite

    contains

    subroutine detector_suite(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Check_hit_circle", hit_circle) &
                ! ! new_unittest("Vector_subtract", vector_sub), &
                ! ! new_unittest("Vector_multiply", vector_mult), &
                ! ! new_unittest("Vector_dot", vector_dot) &
                ]

    end subroutine detector_suite

    subroutine hit_circle(error)

        use vector_class, only : vector

        type(error_type), allocatable, intent(out) :: error

        type(hit_t) :: hitpoint
        type(circle_dect) :: a
        type(vector) :: pos, dir
        integer :: layer, nbins
        real(kind=wp) :: radius, maxval, val, t
        logical :: flag

        pos = vector(0._wp, 0._wp, 0._wp)
        dir = vector(1._wp, 0._wp, 0._wp)
        layer = 1
        radius = 0.5
        nbins = 100
        maxval = 100._wp
        a = circle_dect(pos, dir, layer, radius, nbins, maxval, .false.)

        pos = vector(0._wp, 0._wp, 0._wp)
        dir = vector(1._wp, 0._wp, 0._wp)
        val = 1._wp
        layer = 1
        hitpoint = hit_t(pos, dir, val, layer)

        flag = a%check_hit(hitpoint, t)

        call check(error, flag, .true.)
        if(allocated(error))return

    end subroutine hit_circle

end module testsDetectorMod