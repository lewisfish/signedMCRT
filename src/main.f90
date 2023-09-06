program mcpolar
!! Entry point for program
    use kernels, only : weight_scatter, pathlength_scatter

    integer :: num_args, i
    character(len=64), allocatable :: args(:)

    num_args = command_argument_count()
    if(num_args == 0)then
        allocate(args(1))
        args(1) = "scat_test.toml"
    else
        allocate(args(num_args))
        do i = 1, num_args
            call get_command_argument(i, args(i))
        end do
    end if
    
    ! call weight_scatter(trim(args(1)))
    call pathlength_scatter(trim(args(1)))

end program