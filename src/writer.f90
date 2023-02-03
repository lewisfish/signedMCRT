module writer_mod
! module provides output routines in raw binary and .nrrd formats
!
!
    use constants, only : wp

    implicit none

    interface nrrd_write
        module procedure write_3d_r8_nrrd
    end interface nrrd_write

    interface raw_write
        module procedure write_3d_r8_raw
    end interface raw_write

    private
    public :: normalise_fluence, write_fluence, write_detected_photons

    contains
        function normalise_fluence(grid, array, nphotons) result(out)
        ! normalise fluence in the Lucy 1999 way
            
            use gridMod

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


        subroutine write_detected_photons(detectors)

            use detector_mod
            use constants, only: fileplace
            use utils, only : str

            type(dect_array), intent(in) :: detectors(:)

            integer :: i, j, u
            character(len=:), allocatable :: hdr

            do i = 1, size(detectors)
                open(newunit=u, file=trim(fileplace)//"detectors/detector_"//str(i)//".dat")
                associate(x => detectors(i)%p)
                    select type(x)
                    type is(circle_dect)
                        ! hdr = "# pos, layer, nbins, bin_wid, radius"//new_line("a")//str(x%pos)//","//str(x%layer)//","//str(x%nbins)//","//str(x%bin_wid)//","//str(x%radius)
                        ! write(u, "(a)")hdr
                        ! write(u, "(a)")"#data:"
                        do j = 1, x%nbins
                            write(u,*)real(j,kind=wp) * x%bin_wid, x%data(j)
                        end do
                    type is(annulus_dect)
                        ! hdr = "#pos, layer, nbins, bin_wid, radius1, radius2"//new_line("a")//str(x%pos)//","//str(x%layer)//","//str(x%nbins)//","//str(x%bin_wid)//","//str(x%r1)//","//str(x%r2)
                    type is(camera)
                        print*,"Warning not yet implmented!"
                    end select
                    end associate
                close(u)
            end do

        end subroutine write_detected_photons


        subroutine write_fluence(array, filename, state, dict, overwrite)
        ! routine automatically selects which way to write out results based upon file extension
            
            use sim_state_mod, only : settings_t
            use tomlf,         only : toml_table, get_value

            type(settings_t),           intent(IN)    :: state
            real(kind=wp),              intent(IN)    :: array(:,:,:)
            character(*),               intent(IN)    :: filename
            type(toml_table), optional, intent(INOUT) :: dict
            logical,          optional, intent(IN)    :: overwrite

            Logical :: over_write
            integer :: pos
            
            if(present(overwrite))then
                over_write = overwrite
            else
                over_write = state%overwrite
            end if

            pos = index(filename, ".nrrd")
            if(pos > 0)then
                if(present(dict))then
                    call nrrd_write(array, filename, over_write, dict)
                else
                    call nrrd_write(array, filename, over_write)
                end if
                return
            end if

            pos = index(filename, ".raw")
            if(pos > 0)then
                call raw_write(array, filename, over_write)
                return
            end if

            pos = index(filename, ".dat")
            if(pos > 0)then
                call raw_write(array, filename, over_write)
                return
            end if

            error stop "File type not supported!"

        end subroutine write_fluence

        subroutine write_3d_r8_raw(array, filename, overwrite)

            real(kind=wp), intent(IN) :: array(:, :, :)
            character(*),  intent(IN) :: filename
            logical,       intent(IN) :: overwrite

            integer :: u
            character(len=:), allocatable :: file

            if(check_file(filename) .and. .not. overwrite)then
                file = get_new_file_name(filename)
            else
                file = filename
            end if
            open(newunit=u,file=file,access='stream',status='REPLACE',form='unformatted')
            write(u) array
            close(u)

        end subroutine write_3d_r8_raw


        function get_new_file_name(file) result(res)

            use utils, only : str

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
                    
            character(len=*), intent(IN) :: file

            inquire(file=trim(file), exist=res)
        
        end function check_file

        subroutine write_hdr(u, sizes, type)

            use utils, only : str

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

        subroutine write_3d_r8_nrrd(array, filename, overwrite, dict)
            
            use tomlf,           only : toml_table, toml_serializer
            use iso_fortran_env, only : int32, int64, real32, real64
            use utils,    only : str
        
            character(*),               intent(IN)    :: filename
            real(kind=wp),              intent(IN)    :: array(:, :, :)
            type(toml_table), optional, intent(INOUT) :: dict
            logical,                    intent(IN)    :: overwrite

            character(len=:), allocatable :: file
            integer :: u
            type(toml_serializer) :: ser

            if(check_file(filename) .and. .not. overwrite)then
                file = get_new_file_name(filename)
            else
                file = filename
            end if

            open(newunit=u,file=file,form="formatted")
            !to do fix precision
            call write_hdr(u, [size(array, 1), size(array, 2), size(array, 3)], "double")

            if(present(dict))then
                ser = toml_serializer(u)
                call dict%accept(ser)
            end if
            write(u,"(A)")new_line("C")
            close(u)
            open(newunit=u,file=file,access="stream",form="unformatted",position="append")
            write(u)array
            close(u)
        
        end subroutine write_3d_r8_nrrd
end module writer_mod
