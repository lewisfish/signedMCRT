module util

    use testdrive, only : testsuite_type

    implicit none

    private
    public :: grow_suite

    contains

    subroutine grow_suite(suite1, suite2)

        type(testsuite_type), allocatable :: suite1(:), suite2(:), tmp(:)

        allocate(tmp(size(suite1) + size(suite2)))
        tmp(1:size(suite1)) = suite1
        tmp(size(suite1)+1:) = suite2
        deallocate(suite1)
        call move_alloc(tmp, suite2)

    end subroutine grow_suite
end module util
program test_all
    use, intrinsic :: iso_fortran_env, only: error_unit
    
    use util, only : grow_suite
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type, context_t
    
    use testsDetectorMod
    use testScatterMod
    use testsVec4Mod
    use testsIOMod
    use testsMatrixMod
    use testsPhotonMod
    use testsSDFMod 
    use testsVecMod
    use testsFresnelMod
    use testsPiecewiseMod
    use testsrandom

    implicit none

    type(testsuite_type), allocatable :: testsuites(:), tmp(:)
    type(context_t) :: context
    integer :: i, stat
    character(len=*), parameter :: fmt='("#", *(1x, a))'

    testsuites = [new_testsuite("Suite: Detector Class", detector_suite, context), &
                  new_testsuite("Suite: End to End tests", End_to_End_suite, context) &
                 ]
    
    call random_suite(tmp, context)
    call grow_suite(tmp, testsuites)
             
    call Piecewise_suite(tmp, context)
    call grow_suite(tmp, testsuites)
    
    call Vector_suite(tmp, context)
    call grow_suite(tmp, testsuites)

    call vec4_suite(tmp, context)
    call grow_suite(tmp, testsuites)

    call Matrix_suite(tmp, context)
    call grow_suite(tmp, testsuites)

    call photon_suite(tmp, context)
    call grow_suite(tmp, testsuites)

    call Fresnel_suite(tmp, context)
    call grow_suite(tmp, testsuites)

    call SDF_suite(tmp, context)
    call grow_suite(tmp, testsuites)

    do i = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(i)%name
        call run_testsuite(testsuites(i)%collect, error_unit, stat, parallel=.true., context=context)
    end do
    call context%report()
    if(stat > 0)error stop 1
end program test_all