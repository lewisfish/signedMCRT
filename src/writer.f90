module writer_mod
! module provides output routines in raw binary and .nrrd formats
!
!
    use constants, only : wp
    use utils,     only : str

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


        function normalise_fluence(grid, array, nphotons) result(out)
        ! normalise fluence in the Lucy 1999 way
            
            use gridMod

            implicit none

            type(cart_grid), intent(IN) :: grid
            real(kind=wp),   intent(IN) :: array(:, :, :)
            integer,         intent(IN) :: nphotons
            
            real(kind=wp), allocatable :: out(:, :, :)

            real(kind=wp) :: xmax, ymax, zmax
            integer       :: nxg, nyg, nzg

            nxg = grid%nxg
            nyg = grid%nyg
            nzg = grid%nzg
            xmax = grid%xmax
            ymax = grid%ymax
            zmax = grid%zmax

            allocate(out(size(array, 1), size(array, 2), size(array, 3)))

            out  = array * ((2._wp*xmax*2._wp*ymax)/(nphotons * (2._wp * xmax / nxg) * (2._wp * ymax / nyg) * (2._wp * zmax / nzg)))

        end function normalise_fluence


        subroutine write(array, filename, dict)
        ! routine automatically selects which way to write ouresults based upon file extension
            use fhash,        only : fhash_tbl_t

            implicit none
        
            real(kind=wp),          intent(IN) :: array(:,:,:)
            character(*),           intent(IN) :: filename
            type(fhash_tbl_t), optional, intent(INOUT) :: dict

            integer :: pos
            
            pos = index(filename, ".nrrd")
            if(pos > 0)then
                if(present(dict))then
                    call nrrd_write(array, filename, dict)
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

        subroutine write_3d_r8_nrrd(array, filename, dict)
            
            use fhash,           only : fhash_tbl_t, key=>fhash_tbl_t, fhash_key_t
            use iso_fortran_env, only : int32, int64, real32, real64

            implicit none
        
            character(*),                intent(IN)    :: filename
            real(kind=wp),               intent(IN)    :: array(:, :, :)
            type(fhash_tbl_t), optional, intent(INOUT) :: dict

            class(fhash_key_t), pointer   :: keyd
            character(len=:), allocatable :: key_out, file, type_str
            integer :: u, i, stat

            integer(kind=int32) :: val_i32
            integer(kind=int64) :: val_i64
            real(kind=real32)   :: val_r32
            real(kind=real64)   :: val_r64
            character(len=:), allocatable :: val_char
            logical :: val_bool

            if(check_file(filename))then
                file = get_new_file_name(filename)
            else
                file = filename
            end if

            open(newunit=u,file=file,form="formatted")
            !to do fix precision
            call write_hdr(u, [size(array, 1), size(array, 2), size(array, 3)], "double")

            if(present(dict))then
                do
                    call dict%next(keyd, key_out, type_str, stat)
                    if(stat < 0)exit
                    select case(type_str)
                        case("integer32")
                            call dict%get(keyd, val_i32)
                            write(u, "(A)")trim(key_out)//": "//str(val_i32)
                        case("integer64")
                            call dict%get(keyd, val_i64)
                            write(u, "(A)")trim(key_out)//": "//str(val_i64)
                        case("real32")
                            call dict%get(keyd, val_r32)
                            write(u, "(A)")trim(key_out)//": "//str(val_r32)
                        case("real64")
                            call dict%get(keyd, val_r64)
                            write(u, "(A)")trim(key_out)//": "//str(val_r64)
                        case("character*")
                            call dict%get(keyd, val_char)
                            write(u, "(A)")trim(key_out)//": "//val_char
                        case("logical")
                            call dict%get(keyd, val_bool)
                            write(u, "(A)")trim(key_out)//": "//str(val_bool)
                        case default
                            continue
                    end select
                end do
            end if
            write(u,"(A)")new_line("C")
            close(u)
            open(newunit=u,file=file,access="stream",form="unformatted",position="append")
            write(u)array
            close(u)
        
        end subroutine write_3d_r8_nrrd
end module writer_mod
