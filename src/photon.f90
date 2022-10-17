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
        integer       :: cnts, bounces                !  number of sdf evals.
        real(kind=wp) :: weight, step

        procedure(generic_emit), pointer :: emit => null()

    end type photon

    interface photon
        module procedure init_source
        module procedure init_photon
    end interface photon

    abstract interface
        subroutine generic_emit(this, dict)
            
            use tomlf, only : toml_table, get_value

            import :: photon
            class(photon) :: this
            type(toml_table), optional, intent(inout) :: dict

        end subroutine generic_emit
    end interface

    type(photon) :: photon_origin

    private
    public :: photon, init_source, set_photon

    contains
        
        subroutine set_photon(pos, dir)

            type(vector), intent(in) :: pos, dir

            associate(p => photon_origin)
                p%pos = pos
                p%nxp = dir%x
                p%nyp = dir%y
                p%nzp = dir%z
            end associate

        end subroutine set_photon

        type(photon) function init_photon(val)

            real(kind=wp), intent(in) :: val

            init_photon%pos = vector(val, val, val)
            init_photon%nxp = val
            init_photon%nyp = val
            init_photon%nzp = val
            init_photon%sint = val
            init_photon%cost = val
            init_photon%sinp = val
            init_photon%cosp = val
            init_photon%phi = val
            init_photon%wavelength = val
            init_photon%zcell = int(val)
            init_photon%ycell = int(val)
            init_photon%zcell = int(val)
            init_photon%tflag = .true.
            init_photon%layer = int(val)
            init_photon%id = int(val)
            init_photon%cnts = int(val)
            init_photon%bounces = int(val)
            init_photon%weight = val
            init_photon%step = val 

        end function init_photon

        type(photon) function init_source(choice)

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
            elseif(choice == "circular")then
                init_source%emit => circular
            else
                error stop "No such source!"
            end if


        end function init_source

        subroutine circular(this, dict)
            ! circular source
    
                use sim_state_mod, only : state
                use random,        only : ran2
                use constants,     only : twoPI
                use tomlf,         only : toml_table, get_value
                ! use sdfs,          only : rotationAlign, rotmat
                use mat_class,     only : invert
                use vector_class
            
                class(photon) :: this
                type(toml_table), optional, intent(inout) :: dict
                
                type(vector) :: pos, dir, tmp, axis
                integer :: cell(3)
                real(kind=wp) :: t(4,4), radius, r, theta, angle

                this%pos = photon_origin%pos
                this%nxp = photon_origin%nxp
                this%nyp = photon_origin%nyp
                this%nzp = photon_origin%nzp

                call get_value(dict, "radius", radius)

                ! !https://math.stackexchange.com/a/1681815
                dir = vector(this%nxp, this%nyp, this%nzp)
                dir = dir%magnitude()

                tmp = vector(0._wp, 0._wp, 1._wp)

                axis = tmp .cross. dir
                axis = axis%magnitude()
                angle = (dir .dot. axis) / (length(dir) * length(axis))
                ! t = rotmat(axis, angle)
                ! print*,t(:,1)
                ! print*,t(:,2)
                ! print*,t(:,3)
                ! print*,t(:,4)
! 
                print*," "
                print*,"tmp ",tmp
                print*,"dir ",dir
                print*,"orig",this%pos
                print*,"new ",this%pos .dot. invert(t)
                print*,"axis",axis
                print*,"angl",angle
                stop

                r = radius * sqrt(ran2())
                theta = ran2() * TWOPI
                pos = this%pos
                this%pos%x = this%pos%x + r * cos(theta)
                this%pos%y = this%pos%y + r * sin(theta)

                ! b = vector(this%nxp, this%nyp, this%nzp)
                ! b = b%magnitude()
                ! t = rotmat(vector(0._wp,0._wp,-1._wp), b)
                print*,this%pos
                this%pos = this%pos - pos 
                this%pos = this%pos .dot. t
                print*,this%pos + pos
                print*," "


                this%tflag  = .false.
                this%cnts   = 0
                this%bounces = 0
                this%layer  = 1
                this%weight = 1.0_wp
    
                ! Linear Grid 
                cell = state%grid%get_voxel(this%pos)
                this%xcell = cell(1)
                this%ycell = cell(2)
                this%zcell = cell(3)
            
            end subroutine circular


        subroutine point(this, dict)
        !isotropic point source

            use sim_state_mod, only : state
            use random,        only : ran2
            use constants,     only : twoPI
            use tomlf,         only : toml_table, get_value
        
            class(photon) :: this
            type(toml_table), optional, intent(inout) :: dict
            
            integer :: cell(3)

            this%pos = photon_origin%pos
            this%nxp = photon_origin%nxp
            this%nyp = photon_origin%nyp
            this%nzp = photon_origin%nzp

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
            this%bounces = 0
            this%layer  = 1
            this%weight = 1.0_wp

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
            use tomlf,         only : toml_table, get_value

            class(photon) :: this
            type(toml_table), optional, intent(inout) :: dict

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
            this%bounces = 0
            this%cnts = 0
            this%weight = 1.0_wp

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)

        end subroutine focus


        subroutine uniform(this, dict)
        !uniformly illuminate a surface of the simulation media
            use random,        only : ranu, ran2, randint
            use sim_state_mod, only : state
            use tomlf,         only : toml_table, get_value

            class(photon) :: this
            type(toml_table), optional, intent(inout) :: dict

            integer       :: cell(3)
            type(vector)  :: pos1, pos2, pos3
            real(kind=wp) :: rx, ry

            this%nxp = photon_origin%nxp
            this%nyp = photon_origin%nyp
            this%nzp = photon_origin%nzp

            this%cost = this%nzp
            this%sint = sqrt(1._wp - this%cost**2)

            this%phi = atan2(this%nyp, this%nxp)
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)

            call get_value(dict, "pos1%x", pos1%x)
            call get_value(dict, "pos1%y", pos1%y)
            call get_value(dict, "pos1%z", pos1%z)

            call get_value(dict, "pos2%x", pos2%x)
            call get_value(dict, "pos2%y", pos2%y)
            call get_value(dict, "pos2%z", pos2%z)

            call get_value(dict, "pos3%x", pos3%x)
            call get_value(dict, "pos3%y", pos3%y)
            call get_value(dict, "pos3%z", pos3%z)

            rx = ran2()
            ry = ran2()
            this%pos%x = pos1%x + rx * pos2%x + ry * pos3%x
            this%pos%y = pos1%y + rx * pos2%y + ry * pos3%y
            this%pos%z = pos1%z + rx * pos2%z + ry * pos3%z

            this%tflag = .false.
            this%cnts = 0
            this%bounces = 0
            this%weight = 1.0_wp

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)

        end subroutine uniform


        subroutine pencil(this, dict)

            use random,        only : ranu
            use sim_state_mod, only : state
            use tomlf,         only : toml_table, get_value

            implicit none

            class(photon) :: this
            type(toml_table), optional, intent(inout) :: dict

            integer :: cell(3)

            this%pos = photon_origin%pos
            this%nxp = photon_origin%nxp
            this%nyp = photon_origin%nyp
            this%nzp = photon_origin%nzp

            this%phi = atan2(this%nyp, this%nxp)
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)
            this%cost = this%nzp
            this%sint = sqrt(1._wp - this%cost**2)

            this%tflag = .false.
            this%bounces = 0
            this%cnts = 0
            this%weight = 1.0_wp

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)
        end subroutine pencil
        
        subroutine annulus(this, dict)

            use utils,         only : deg2rad
            use random,        only : rang
            use surfaces,      only : reflect_refract
            use sim_state_mod, only : state
            use tomlf,         only : toml_table, get_value

            class(photon) :: this
            type(toml_table), optional, intent(inout) :: dict

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
                call reflect_refract(dir, normal, axicon_n, 1._wp, flag, k)
                if(.not. flag)then
                    exit
                end if
            end do

            !move packet to L1 location
            pos = pos + dir * (-total_length / dir%z)

            !focus to a point
            pos%z = 4.e0_wp !f distance 40mm
            call rang(x, y, 0._wp, 5.e-4_wp)

            call get_value(dict, "focus", z)!not implmented yet!
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
            this%bounces = 0
            this%weight = 1.0_wp

            !teleport to just inside medium
            this%pos%z = state%grid%zmax - 1e-8_wp

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)

        end subroutine annulus
end module photonMod