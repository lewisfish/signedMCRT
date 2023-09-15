module setupMod
!! This file sets up some simulations variables and assigns the geometry for the simulation.
   
    use constants, only : wp
    use tomlf

    implicit none

    private
    public  :: setup_simulation, dealloc_array, directory

    contains

        subroutine setup_simulation(sdfarray, dict)
        !! Read in parameters
        !! Setup up various simulation parameters and routines

            use sdfs,          only : sdf
            use setupGeometry
            use sim_state_mod, only : settings => state
            use vector_class
            
            !> dictionary used to store metadata
            type(toml_table), optional, intent(INOUT) :: dict
            !> output array of geometry
            type(sdf), allocatable,     intent(OUT)   :: sdfarray(:)

            !allocate and set arrays to 0
            call alloc_array(settings%grid%nxg, settings%grid%nyg, settings%grid%nzg)
            call zarray()

            ! setup geometry using SDFs
            select case(settings%experiment)
                case("logo")
                    sdfarray = setup_logo()
                case("omg")
                    sdfarray = setup_omg_sdf()
                case("scat_test")
                    sdfarray = setup_scat_test(dict)
                case("scat_test2")
                    sdfarray = setup_scat_test2(dict)
                case("aptran")
                    sdfarray = setup_sphere()
                case("vessels")
                    sdfarray = get_vessels()
                case("sphere_scene")
                    sdfarray = setup_sphere_scene(dict)
                case("test_egg")
                    sdfarray = setup_egg()
                case default
                    error stop "no such routine"
            end select

        end subroutine setup_simulation

        subroutine directory()
        !!  subroutine defines vars to hold paths to various folders     
  
            use constants, only : homedir, fileplace, resdir

            character(len=256) :: cwd
            logical :: dataExists, jmeanExists, depositExists, detectorsExists, phasorExists

            !get current working directory
            call get_environment_variable('PWD', cwd)
  
            ! get 'home' dir from cwd
            homedir = trim(cwd)
            ! get data dir
            fileplace = trim(homedir)//'/data/'
            !check if data directory and subdirectories exists. if not create it
#ifdef __GFORTRAN__
            inquire(file=trim(fileplace)//"/.", exist=dataExists)
            inquire(file=trim(fileplace)//"/jmean/.", exist=jmeanExists)
            inquire(file=trim(fileplace)//"/deposit/.", exist=depositExists)
            inquire(file=trim(fileplace)//"/detectors/.", exist=detectorsExists)
            inquire(file=trim(fileplace)//"/phasor/.", exist=phasorExists)
#elif __INTEL_COMPILER
            inquire(directory=trim(fileplace), exist=dataExists)
            inquire(directory=trim(fileplace)//"/jmean", exist=jmeanExists)
            inquire(directory=trim(fileplace)//"/deposit", exist=depositExists)
            inquire(directory=trim(fileplace)//"/detectors", exist=detectorsExists)
            inquire(directory=trim(fileplace)//"/phasor", exist=phasorExists)
#else 
    error stop "Compiler not supported!"
#endif
            if(.not. dataExists)then
                call create_directory("", dataExists, "", .false.)
                call create_directory("jmean/", jmeanExists, "data/", .false.)
                call create_directory("deposit/", depositExists, "data/", .false.)
                call create_directory("detectors/", detectorsExists, "data/", .false.)
                call create_directory("phasor/", phasorExists, "data/", .false.)
            else
                call create_directory("jmean/", jmeanExists, "data/", .true.)
                call create_directory("deposit/", depositExists, "data/", .true.)
                call create_directory("detectors/", detectorsExists, "data/", .true.)
                call create_directory("phasor/", phasorExists, "data/", .true.)
            end if

            ! get res dir
            resdir = trim(homedir)//'/res/'

        end subroutine directory


        subroutine create_directory(name, flag, appendname, newline)
        !! create directories if they don't exist
            use constants, only : fileplace

            character(*),      intent(in) :: name, appendname
            logical,           intent(in) :: flag
            logical, optional, intent(in) :: newline

            character(len=:), allocatable :: mkdirCMD

            if(.not. flag)then
                mkdirCMD = "mkdir -p "//trim(fileplace)//name
                call execute_command_line(mkdirCMD)
                ! output correct message for base data dir
                if(len(name) == 0)then
                    mkdirCMD = "Created "//appendname//"data/"                    
                else
                    mkdirCMD = "Created "//appendname//name
                end if
                if(newline)mkdirCMD = mkdirCMD//new_line("a")
                print*,mkdirCMD
            end if

        end subroutine create_directory

        subroutine zarray
        !! zero data arrays
            use iarray

            !sets all arrays to zer

            phasor = 0._wp
            phasorGLOBAL = 0._wp
            jmean = 0._wp
            jmeanGLOBAL = 0._wp
            absorb = 0.0_wp
            absorbGLOBAL = 0.0_wp

        end subroutine zarray


        subroutine alloc_array(nxg, nyg, nzg)
        !!  subroutine allocates allocatable arrays  
  
            use iarray
            !> grid size
            integer, intent(IN) :: nxg, nyg, nzg

            allocate(phasor(nxg, nyg, nzg), phasorGLOBAL(nxg, nyg, nzg))
            allocate(jmean(nxg, nyg, nzg), jmeanGLOBAL(nxg, nyg, nzg))
            allocate(absorb(nxg, nyg, nzg), absorbGLOBAL(nxg, nyg, nzg))

        end subroutine alloc_array

        subroutine dealloc_array()
        !! deallocate data arrays
            use iarray

            deallocate(jmean)
            deallocate(jmeanGLOBAL)
            deallocate(absorb)
            deallocate(absorbGLOBAL)
            deallocate(phasor)
            deallocate(phasorGLOBAL)
        end subroutine dealloc_array
end module setupMod