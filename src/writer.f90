module writer_mod
! module provides output routines in raw binary and .nrrd formats
!
!
    use constants, only : wp
    use utils,     only : str
    use dict_mod

implicit none

    interface nrrd_write
        module procedure write_3d_r8_nrrd
    end interface nrrd_write


    interface raw_write
        module procedure write_3d_r8_raw
    end interface raw_write

    private
    public :: raw_write, nrrd_write, normalise_fluence, write

    contains


        function normalise_fluence(array, nphotons) result(out)
        ! normalise fluence in the Lucy 1999 way
            
            use sim_state_mod, only : state

            implicit none

            real(kind=wp),   intent(IN) :: array(:, :, :)
            integer,         intent(IN) :: nphotons
            
            real(kind=wp), allocatable :: out(:, :, :)

            real(kind=wp) :: xmax, ymax, zmax
            integer       :: nxg, nyg, nzg

            nxg = state%grid%nxg
            nyg = state%grid%nyg
            nzg = state%grid%nzg
            xmax = state%grid%xmax
            ymax = state%grid%ymax
            zmax = state%grid%zmax

            allocate(out(size(array, 1), size(array, 2), size(array, 3)))

            out  = array * ((2._wp*xmax*2._wp*ymax)/(nphotons * (2._wp * xmax / nxg) * (2._wp * ymax / nyg) * (2._wp * zmax / nzg)))

        end function normalise_fluence


        subroutine write(array, filename, dicts)
        ! routine automatically selects which way to write ouresults based upon file extension
            implicit none
        
            real(kind=wp),          intent(IN) :: array(:,:,:)
            character(*),           intent(IN) :: filename
            type(dict_t), optional, intent(IN) :: dicts

            integer :: pos
            
            pos = index(filename, ".nrrd")
            if(pos > 0)then
                if(present(dicts))then
                    call nrrd_write(array, filename, dicts)
                else
                    call nrrd_write(array, filename)
                end if
                return
            end if

            pos = index(filename, ".raw")
            if(pos > 0)then
                call raw_write(array, filename)
                return
            end if

            pos = index(filename, ".dat")
            if(pos > 0)then
                call raw_write(array, filename)
                return
            end if

            error stop "File type not supported!"

        end subroutine write

        subroutine write_3d_r8_raw(array, filename)

            implicit none

            real(kind=wp), intent(IN) :: array(:, :, :)
            character(*),  intent(IN) :: filename

            integer :: u
            character(len=:), allocatable :: file

            if(check_file(filename))then
                file = get_new_file_name(filename)
            else
                file = filename
            end if
            open(newunit=u,file=file,access='stream',status='REPLACE',form='unformatted')
            write(u) array
            close(u)

        end subroutine write_3d_r8_raw


        function get_new_file_name(file) result(res)

            implicit none

            character(len=*), intent(IN) :: file
            character(len=:), allocatable :: res
            integer :: pos, i

            i = 1
            do
                pos = scan(trim(file), ".", back=.true.)
                res = file(1:pos-1)//" ("//str(i)//")"//file(pos:)
                if(.not. check_file(res))exit
                i = i + 1
            end do

        end function get_new_file_name

        logical function check_file(file) result(res)
            
            implicit none
        
            character(len=*), intent(IN) :: file

            inquire(file=trim(file), exist=res)
        
        end function check_file

        subroutine write_hdr(u, sizes, type)

            implicit none

            character(*), intent(IN) :: type
            integer,      intent(IN) :: sizes(:), u
            
            character(len=100) :: string
            integer :: i

            string = ""
            do i = 1, size(sizes)
                if(i == 1)then
                    string = str(sizes(i))            
                else
                    string = trim(string) // " " // str(sizes(i))
                end if
            end do

            write(u,"(A)")"NRRD0004"
            write(u,"(A)")"type: "//type
            write(u,"(A)")"dimension: "//str(size(sizes))
            write(u,"(A)")"sizes: "//trim(string)
            write(u,"(A)")"space dimension: "//str(size(sizes))
            write(u,"(A)")"encoding: raw"
            write(u,"(A)")"endian: little"


        end subroutine write_hdr

        subroutine write_3d_r8_nrrd(array, filename, dicts)
            
            implicit none
        
            character(*),           intent(IN) :: filename
            real(kind=wp),          intent(IN) :: array(:, :, :)
            type(dict_t), optional, intent(IN) :: dicts

            character(len=64) :: key
            character(len=:), allocatable :: file
            integer :: u, i

            if(check_file(filename))then
                file = get_new_file_name(filename)
            else
                file = filename
            end if

            open(newunit=u,file=file,form="formatted")
            !to do fix precision
            call write_hdr(u, [size(array, 1), size(array, 2), size(array, 3)], "double")

            if(present(dicts))then
                do i = 1, dicts%count
                    key = trim(dicts%dict(i)%key)
                    write(u, "(A)")trim(key)//": "//dicts%get_value_str(trim(key))
                end do
            end if

            write(u,"(A)")new_line("C")
            close(u)
            open(newunit=u,file=file,access="stream",form="unformatted",position="append")
            write(u)array
            close(u)
        
        end subroutine write_3d_r8_nrrd
end module writer_mod
