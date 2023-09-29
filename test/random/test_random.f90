module testsrandom

    use random
    use testdrive, only : new_unittest, unittest_type, error_type, check, new_testsuite, testsuite_type, context_t, test_failed
    use constants, only : wp

    implicit none

    private
    public :: random_suite

    contains

    subroutine random_suite(testsuites, context)

        type(testsuite_type), allocatable, intent(out) :: testsuites(:)
        type(context_t) :: context

        testsuites = [new_testsuite("Test Random distributions", collect_suite1, context),&
                      new_testsuite("Test Sequences", collect_suite2, context)&
                     ]

    end subroutine random_suite

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Initalise RNG", test_initRNG),&
                new_unittest("Gaussian Dist.", test_rang),&
                new_unittest("Uniform Dist.", test_uni),&
                new_unittest("Random number wrapper", test_ran2),&
                new_unittest("Random Integer.", test_int)&
                ]

    end subroutine collect_suite1

    subroutine collect_suite2(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Test Seq.", test_Seq)&
                ]

    end subroutine collect_suite2


    subroutine test_initRNG(error)

        use utils, only : str

        type(error_type), allocatable, intent(out) :: error
        integer, allocatable :: seed(:)
        integer              :: n, i

        call random_seed(size=n)
        allocate(seed(n))
        call init_rng(spread(123456789, 1, 8), .true.)
        call random_seed(get=seed)

        do i = 1, n
            if(seed(i) == 123456789)then
                call test_failed(error, "init_rng check failed!", "Expected seed to not be 123456789 "//str(seed(i)))
            end if
            if(allocated(error))return
        end do
        
        call init_rng(spread(123456789, 1, 8))
        call random_seed(get=seed)

        do i = 1, n
            call check(error, seed(i), 123456789)
            if(allocated(error))return
        end do

        call init_rng(fwd=.false.)
        call random_seed(get=seed)

        do i = 1, n
            call check(error, seed(i), 1234567)
            if(allocated(error))return
        end do

    end subroutine test_initRNG
    
    subroutine test_rang(error)

        use utils, only : str

        type(error_type), allocatable, intent(out) :: error

        integer :: i
        integer, parameter :: n = 100000
        real(kind=wp) :: x, y, mean, vals(n), sig

        mean = 0.0_wp
        vals = 0.0_wp

        do i = 1, n
            call rang(x, y, 0._wp, 1._wp)
            mean = mean + x
            vals(i) = x
        end do

        mean = sum(vals) / n
        sig = sqrt(sum(abs(vals - mean)**2) / (n))

        call check(error, mean, 0.0_wp, thr=0.001_wp)
        if(allocated(error))return
    
        call check(error, sig, 1.0_wp, thr=0.001_wp)
        if(allocated(error))return

    end subroutine test_rang

    subroutine test_uni(error)

        use utils, only : str

        type(error_type), allocatable, intent(out) :: error

        integer :: i
        real(kind=wp) :: val

        do i = 1, 10000
            val = ranu(0._wp, 100._wp)
            if(val < 0. .or. val > 100.)then
                call test_failed(error, "ranu check failed!", "Expected a value between [0, 100]! Got "//str(val, 5))
            end if
            if(allocated(error))return
        end do
    end subroutine test_uni

    subroutine test_ran2(error)

        use utils, only : str

        type(error_type), allocatable, intent(out) :: error

        integer :: i
        real(kind=wp) :: val

        do i = 1, 10000
            val = ran2()
            if(val < 0. .or. val > 1.)then
                call test_failed(error, "Ran2 check failed!", "Expected a value between [0.0, 1.0]! Got "//str(val, 5))
            end if
            if(allocated(error))return
        end do

    end subroutine test_ran2

    subroutine test_int(error)

        use utils, only : str

        type(error_type), allocatable, intent(out) :: error

        integer :: i,  val

        do i = 1, 10000
            val = randint(0, 100)
            if(val < 0 .or. val > 100)then
                call test_failed(error, "randint check failed!", "Expected a value between [0, 100]! Got "//str(val, 5))
            end if
            if(allocated(error))return
        end do

    end subroutine test_int

    subroutine test_Seq(error)

        type(error_type), allocatable, intent(out) :: error
        
        type(seq) :: seqs
        real(kind=wp) :: val, vals(4)
        integer :: index, base, i

        index = 1234
        base = 10

        vals = [4./10., 3./100., 2./1000., 1./10000.]

        seqs = seq(index, base)
        val = seqs%next()
        call check(error, val, sum(vals), thr=0.001_wp)
        if(allocated(error))return

    end subroutine test_Seq
end module testsrandom