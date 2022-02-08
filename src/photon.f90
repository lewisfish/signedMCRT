module photonMod
    
    use constants, only : wp
    use vector_class

    implicit none
    
    type :: photon
    
        type(vector)  :: pos                          ! position
        real(kind=wp) :: nxp, nyp, nzp                ! direction vectors
        real(kind=wp) :: sint, cost, sinp, cosp, phi  ! direction cosines
        real(kind=wp) :: wavelength                   ! Only used if tracking the phase
        integer       :: xcell, ycell, zcell          ! grid cell position
        logical       :: tflag                        ! Is photon dead?
        integer       :: layer                        ! id of sdf the packet is inside
        integer       :: id                           ! thread running packet
        integer       :: cnts                         !  number of sdf evals.

        procedure(generic_emit), pointer :: emit => null()

    end type photon

    interface photon
        module procedure init_source
    end interface photon

    abstract interface
        subroutine generic_emit(this, dict)
            
            use fhash, only : fhash_tbl_t

            import :: photon
            class(photon) :: this
            type(fhash_tbl_t), optional, intent(IN) :: dict

        end subroutine generic_emit
    end interface

    private
    public :: photon, init_source

    contains
        
        type(photon) function init_source(choice)

            implicit none

            character(*), intent(IN) :: choice

            if(choice == "uniform")then
                init_source%emit => uniform
            elseif(choice == "pencil")then
                init_source%emit => pencil
            elseif(choice == "annulus")then
                init_source%emit => annulus
            elseif(choice == "focus")then
                init_source%emit => focus
            elseif(choice == "point")then
                init_source%emit => point
            else
                error stop "No such source!"
            end if


        end function init_source

        subroutine point(this, dict)
        !isotropic point source

            use sim_state_mod, only : state
            use random,        only : ran2
            use constants,     only : twoPI
            use fhash,         only : fhash_tbl_t, key=>fhash_key

            implicit none
        
            class(photon) :: this
            type(fhash_tbl_t), optional, intent(IN) :: dict
            
            integer :: cell(3)

            call dict%get(key("pos%x"), value=this%pos%x)
            call dict%get(key("pos%y"), value=this%pos%y)
            call dict%get(key("pos%z"), value=this%pos%z)

            this%phi  = ran2()*twoPI
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)
            this%cost = 2._wp*ran2()-1._wp
            this%sint = sqrt(1._wp - this%cost**2)

            this%nxp = this%sint * this%cosp
            this%nyp = this%sint * this%sinp
            this%nzp = this%cost

            this%tflag  = .false.
            this%cnts   = 0
            this%layer  = 1

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)
        
        end subroutine point

        subroutine focus(this, dict)

            use random,        only : ranu
            use sim_state_mod, only : state
            use utils,         only : deg2rad
            use vector_class,  only : length
            use fhash,         only : fhash_tbl_t

            implicit none

            class(photon) :: this
            type(fhash_tbl_t), optional, intent(IN) :: dict

            type(vector)  :: targ, dir
            real(kind=wp) :: dist
            integer       :: cell(3)

            targ = vector(0._wp,0._wp,0._wp)

            this%pos%x = ranu(-state%grid%xmax, state%grid%xmax)
            this%pos%y = ranu(-state%grid%ymax, state%grid%ymax)
            this%pos%z = state%grid%zmax - 1e-8_wp

            dist = length(this%pos)

            dir = (-1._wp)*this%pos / dist
            dir = dir%magnitude()

            this%nxp = dir%x
            this%nyp = dir%y
            this%nzp = dir%z
            this%phi  = atan2(this%nyp, this%nxp)
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)
            this%cost = this%nzp
            this%sint = sqrt(1._wp - this%cost**2)

            this%nxp = this%sint * this%cosp
            this%nyp = this%sint * this%sinp
            this%nzp = this%cost

            this%tflag = .false.

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)

        end subroutine focus


        subroutine uniform(this, dict)
        !uniformly illuminate a surface of the simulation media
        !TODO change to user defined patch inplace of whole side
            use random,        only : ranu, ran2, randint
            use sim_state_mod, only : state
            use fhash,         only : fhash_tbl_t, key=>fhash_key

            implicit none

            class(photon) :: this
            type(fhash_tbl_t), optional, intent(IN) :: dict

            integer :: cell(3)
            character(len=:),allocatable :: dir

            call dict%get(key("dir"), dir)
            select case(dir)
            case("-z")
                this%pos%x = ranu(-state%grid%xmax, state%grid%xmax)
                this%pos%y = ranu(-state%grid%ymax, state%grid%ymax)
                this%pos%z = state%grid%zmax-1e-8_wp

                this%phi  = 0._wp
                this%cosp = 1._wp
                this%sinp = 0._wp
                this%cost = -1._wp
                this%sint = 0._wp
            case("+z")
                this%pos%x = ranu(-state%grid%xmax, state%grid%xmax)
                this%pos%y = ranu(-state%grid%ymax, state%grid%ymax)
                this%pos%z = -state%grid%zmax + 1e-8_wp

                this%phi  = 0._wp
                this%cosp = 1._wp
                this%sinp = 0._wp
                this%cost = 1._wp
                this%sint = 0._wp

            case("-x")
                this%pos%x = state%grid%xmax-1e-8_wp
                this%pos%y = ranu(-state%grid%ymax, state%grid%ymax)
                this%pos%z = ranu(-state%grid%zmax, state%grid%zmax)

                this%phi  = 0._wp
                this%cosp = 1._wp
                this%sinp = 0._wp
                this%cost = 0._wp
                this%sint = -1._wp
            case("+x")
                this%pos%x = -state%grid%xmax+1e-8_wp
                this%pos%y = ranu(-state%grid%ymax, state%grid%ymax)
                this%pos%z = ranu(-state%grid%zmax, state%grid%zmax)

                this%phi  = 0._wp
                this%cosp = 1._wp
                this%sinp = 0._wp
                this%cost = 0._wp
                this%sint = 1._wp
            case("-y")
                this%pos%x = ranu(-state%grid%xmax, state%grid%xmax)
                this%pos%y = state%grid%xmax-1e-8_wp
                this%pos%z = ranu(-state%grid%zmax, state%grid%zmax)

                this%phi  = 0._wp
                this%cosp = 0._wp
                this%sinp = -1._wp
                this%cost = 0._wp
                this%sint = 1._wp
            case("+y")
                this%pos%x = ranu(-state%grid%xmax, state%grid%xmax)
                this%pos%y = -state%grid%xmax+1e-8_wp
                this%pos%z = ranu(-state%grid%zmax, state%grid%zmax)

                this%phi  = 0._wp
                this%cosp = 0._wp
                this%sinp = 1._wp
                this%cost = 0._wp
                this%sint = 1._wp
            case default
                error stop "No such direction for uniform source"
            end select

            this%nxp = this%sint * this%cosp
            this%nyp = this%sint * this%sinp
            this%nzp = this%cost

            this%tflag = .false.
            this%cnts = 0

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)

        end subroutine uniform


        subroutine pencil(this, dict)

            use random,        only : ranu
            use sim_state_mod, only : state
            use fhash,         only : fhash_tbl_t, key=>fhash_key

            implicit none

            class(photon) :: this
            type(fhash_tbl_t), optional, intent(IN) :: dict

            integer :: cell(3)

            call dict%get(key("pos%x"), this%pos%x)
            call dict%get(key("pos%y"), this%pos%y)
            call dict%get(key("pos%z"), this%pos%z)

            this%phi = 0._wp
            this%cosp = 0._wp
            this%sinp = 0._wp     
            this%cost = -1._wp
            this%sint =  0._wp

            this%nxp = this%sint * this%cosp  
            this%nyp = this%sint * this%sinp
            this%nzp = this%cost

            this%tflag = .false.

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)

        end subroutine pencil
        
        subroutine annulus(this, dict)

            use utils,     only : deg2rad
            use random,    only : rang
            use surfaces,  only : reflect_refract
            use sim_state_mod, only : state
            use fhash,     only : fhash_tbl_t, key=>fhash_key

            implicit none

            class(photon) :: this
            type(fhash_tbl_t), optional, intent(IN) :: dict


            real(kind=wp) :: ra, ta, alpha, Ls, beam_width, rp, da, dist
            real(kind=wp) :: tana, x, y, z, total_length, axicon_n, height, k
            type(vector)  :: pos, dir, normal, centre, newpos
            logical       :: flag
            integer       :: cell(3)


            ra = 4._wp ! 12.7mm
            beam_width = .5e0_wp !100um
            do
                call rang(x, y, 0._wp, beam_width)
                rp = sqrt(x**2 + y**2)
                axicon_n = 1.45_wp
                ta = 0.0e-3_wp  ! 0.0 mm
                alpha = deg2rad(20.0_wp)
                tana = tan(alpha)
                height = tana * ra
                Ls = 10._wp  ! 100mm
                centre = vector(0._wp, 0._wp, 0._wp)
                k  = (rp / height)**2
                
                da = (ra - rp) * tana
                pos = vector(x, y, -da)

                total_length = Ls + height

                dir = vector(0._wp, 0._wp, -1._wp)
                !https://math.stackexchange.com/a/3349448
                normal =vector(2._wp*pos%x, 2._wp*pos%y, -2._wp*ra**2*(pos%z+height)/height**2)
                normal = normal%magnitude()
                call reflect_refract(dir, normal, axicon_n, 1._wp, flag)
                if(.not. flag)then
                    exit
                end if
            end do

            !move packet to L1 location
            pos = pos + dir * (-total_length / dir%z)

            !focus to a point
            pos%z = 4.e0_wp !f distance 40mm
            call rang(x, y, 0._wp, 5.e-4_wp)

            call dict%get(key("focus"), z)!not implmented yet!
            error stop "Not implmented focus yet!"

            newpos = vector(x, y, z)!.775
            dist = sqrt((pos%x - newpos%x)**2 + (pos%y - newpos%y)**2 +(pos%z - newpos%z)**2)

            dir = (-1._wp)*pos / dist
            dist = (state%grid%zmax + 5e-8_wp - pos%z) / dir%z
            pos = pos + dist * dir
            this%pos = pos

            dir = dir%magnitude()

            this%nxp = dir%x
            this%nyp = dir%y
            this%nzp = dir%z

            this%cost = this%nzp
            this%sint = sqrt(1.e0_wp - this%cost**2)

            this%phi = atan2(this%nyp, this%nxp)
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)

            this%tflag = .false.
            this%layer = 3
            this%cnts = 0

            !teleport to just inside medium
            this%pos%z = state%grid%zmax - 1e-8_wp

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)

        end subroutine annulus
end module photonMod