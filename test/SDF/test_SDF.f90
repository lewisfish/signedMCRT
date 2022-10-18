module testsSDFMod

    use sdfs
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use constants, only : wp

    implicit none

    contains

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        ! testsuite = [ &
                ! new_unittest("Vector_add", vector_add) &
                ! ! new_unittest("Vector_subtract", vector_sub), &
                ! ! new_unittest("Vector_multiply", vector_mult), &
                ! ! new_unittest("Vector_dot", vector_dot) &
                ! ]

    end subroutine collect_suite1
end module testsSDFMod