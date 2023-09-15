module writer_mod
!! This module defines all functions that write simulation data to the disk or pre-process data before writing.
!! normalise_fluence. Normalises fluence by number of photons run and size of each voxel. **!Does not normalise by power!**
!! write_fluence. Write out fluence in either raw or nrrd format. Default is nrrd.
!! write_detected_photons. Write out photons detected by detectors.

!! Changes should only be made here if there is a bug or new data types need to be written to disk (phase information) or new file format is needed.

    use constants, only : wp

    implicit none

    interface nrrd_write
        module procedure write_3d_r8_nrrd, write_3d_r4_nrrd
    end interface nrrd_write

    interface raw_write
        module procedure write_3d_r8_raw, write_3d_r4_raw
    end interface raw_write

    private
    public :: normalise_fluence, write_data, write_detected_photons

    contains
        subroutine normalise_fluence(grid, array, nphotons)
        !! normalise fluence in the Lucy 1999 way
            
            use gridMod
            use constants, only : sp

            !> grid class
            type(cart_grid), intent(in) :: grid
            !> array to normalise
            real(kind=sp),   intent(inout) :: array(:, :, :)
            !> number of photons run
            integer,         intent(in) :: nphotons
            
            real(kind=wp) :: xmax, ymax, zmax
            integer       :: nxg, nyg, nzg

            nxg = grid%nxg
            nyg = grid%nyg
            nzg = grid%nzg
            xmax = grid%xmax
            ymax = grid%ymax
            zmax = grid%zmax

array  = array * ((2._sp*xmax*2._sp*ymax)/(nphotons * (2._sp * xmax / nxg) * (2._sp * ymax / nyg) * (2._sp * zmax / nzg)))

        end subroutine normalise_fluence


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


        subroutine write_data(array, filename, state, dict, overwrite)
        !! routine automatically selects which way to write out results based upon file extension
            
            use sim_state_mod, only : settings_t
            use tomlf,         only : toml_table, get_value
            use constants,     only : sp

            !> simulation state
            type(settings_t),           intent(IN)    :: state
            !> array to write out
            real(kind=sp),              intent(IN)    :: array(:,:,:)
            !> filename to save array as
            character(*),               intent(IN)    :: filename
            !> dictionary of metadata
            type(toml_table), optional, intent(INOUT) :: dict
            !> overwrite flag
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

        end subroutine write_data

        subroutine write_3d_r8_raw(array, filename, overwrite)
        !! write 3D array of float64s to disk as raw binary data
            
            !> array to write to disk
            real(kind=wp), intent(IN) :: array(:, :, :)
            !> filename to save array as
            character(*),  intent(IN) :: filename
            !> overwrite flag
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

        subroutine write_3d_r4_raw(array, filename, overwrite)
        !! write 3D array of float32's to disk as raw binary data
            use constants, only : sp

            !> array to write to disk
            real(kind=sp), intent(IN) :: array(:, :, :)
            !> filename to save array as
            character(*),  intent(IN) :: filename
            !> overwrite flag
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

        end subroutine write_3d_r4_raw

        function get_new_file_name(file) result(res)
            !! If file exits, get numeral to append to filename 

            use utils, only : str

            !> file to be checked
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
            !! Functional wrapper around inquire to check if file exits

            !> file to be checked
            character(len=*), intent(IN) :: file

            inquire(file=trim(file), exist=res)
        
        end function check_file

        subroutine write_hdr(u, sizes, type)
            !! write out header information for .nrrd file format
            use utils, only : str

            !> data dtype
            character(*), intent(IN) :: type
            !> file handle
            integer,      intent(IN) :: u
            !> dimensions of data
            integer,      intent(IN) :: sizes(:)
            
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
            !! write 3D array of float64's to .nrrd fileformat

            use tomlf,           only : toml_table, toml_dump, toml_error
            use iso_fortran_env, only : int32, int64, real32, real64
            use utils,    only : str
            
            !> filename
            character(*),               intent(IN)    :: filename
            !> array to be written to disk
            real(kind=wp),              intent(IN)    :: array(:, :, :)
            !> dictionary of metadata
            type(toml_table), optional, intent(INOUT) :: dict
            !> overwrite flag
            logical,                    intent(IN)    :: overwrite

            type(toml_error), allocatable :: error
            character(len=:), allocatable :: file
            integer :: u

            if(check_file(filename) .and. .not. overwrite)then
                file = get_new_file_name(filename)
            else
                file = filename
            end if

            open(newunit=u,file=file,form="formatted")
            !to do fix precision
            call write_hdr(u, [size(array, 1), size(array, 2), size(array, 3)], "double")

            if(present(dict))then
                call toml_dump(dict, u, error)
            end if
            write(u,"(A)")new_line("C")
            close(u)
            open(newunit=u,file=file,access="stream",form="unformatted",position="append")
            write(u)array
            close(u)
        
        end subroutine write_3d_r8_nrrd

        subroutine write_3d_r4_nrrd(array, filename, overwrite, dict)
            !! write 3D array of float32's to .nrrd fileformat

            use tomlf,           only : toml_table, toml_dump, toml_error
            use iso_fortran_env, only : int32, int64, real32, real64
            use utils,           only : str
            use constants,       only : sp
            
            !> filename
            character(*),               intent(IN)    :: filename
            !> array to be written to disk
            real(kind=sp),              intent(IN)    :: array(:, :, :)
            !> dictionary of metadata
            type(toml_table), optional, intent(INOUT) :: dict
            !> overwrite flag
            logical,                    intent(IN)    :: overwrite

            type(toml_error), allocatable :: error
            character(len=:), allocatable :: file
            integer :: u

            if(check_file(filename) .and. .not. overwrite)then
                file = get_new_file_name(filename)
            else
                file = filename
            end if

            open(newunit=u,file=file,form="formatted")
            !to do fix precision
            call write_hdr(u, [size(array, 1), size(array, 2), size(array, 3)], "float")

            if(present(dict))then
                call toml_dump(dict, u, error)
            end if
            write(u,"(A)")new_line("C")
            close(u)
            open(newunit=u,file=file,access="stream",form="unformatted",position="append")
            write(u)array
            close(u)
        
        end subroutine write_3d_r4_nrrd
end module writer_mod
