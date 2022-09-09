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
                use fhash,         only : fhash_tbl_t, key=>fhash_key
                use sdfs, only : rotmat, rotationAlign
                use vector_class

                implicit none
            
                class(photon) :: this
                type(fhash_tbl_t), optional, intent(IN) :: dict
                
                type(vector) :: pos, b
                integer :: cell(3)
                real(kind=wp) :: t(4,4), radius,r,theta
    
                call dict%get(key("pos%x"), value=this%pos%x)
                call dict%get(key("pos%y"), value=this%pos%y)
                call dict%get(key("pos%z"), value=this%pos%z)
    
                call dict%get(key("dir%x"), this%nxp)
                call dict%get(key("dir%y"), this%nyp)
                call dict%get(key("dir%z"), this%nzp)
        
                call dict%get(key("radius"), radius)

                r = radius * sqrt(ran2())
                theta = ran2() * TWOPI
                this%pos%x = this%pos%x + r * cos(theta)
                this%pos%y = this%pos%y + r * sin(theta)

                ! b = vector(this%nxp, this%nyp, this%nzp)
                ! b = b%magnitude()
                ! t = rotmat(vector(0._wp,0._wp,-1._wp), b)

                ! pos = this%pos .dot. t

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
            this%bounces = 0
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
            this%bounces = 0
            this%cnts = 0

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

            integer        :: cell(3)
            type(vector)   :: pos1, pos2, pos3
            real(kind=wp) :: rx, ry

            call dict%get(key("dir%x"), this%nxp)
            call dict%get(key("dir%y"), this%nyp)
            call dict%get(key("dir%z"), this%nzp)

            this%cost = this%nzp
            this%sint = sqrt(1._wp - this%cost**2)

            this%phi = atan2(this%nyp, this%nxp)
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)

            call dict%get(key("pos1%x"), pos1%x)
            call dict%get(key("pos1%y"), pos1%y)
            call dict%get(key("pos1%z"), pos1%z)

            call dict%get(key("pos2%x"), pos2%x)
            call dict%get(key("pos2%y"), pos2%y)
            call dict%get(key("pos2%z"), pos2%z)

            call dict%get(key("pos3%x"), pos3%x)
            call dict%get(key("pos3%y"), pos3%y)
            call dict%get(key("pos3%z"), pos3%z)

            rx = ran2()
            ry = ran2()
            this%pos%x = pos1%x + rx * pos2%x + ry * pos3%x
            this%pos%y = pos1%y + rx * pos2%y + ry * pos3%y
            this%pos%z = pos1%z + rx * pos2%z + ry * pos3%z

            this%tflag = .false.
            this%cnts = 0
            this%bounces = 0

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

            call dict%get(key("dir%x"), this%nxp)
            call dict%get(key("dir%y"), this%nyp)
            call dict%get(key("dir%z"), this%nzp)

            this%phi = atan2(this%nyp, this%nxp)
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)
            this%cost = this%nzp
            this%sint = sqrt(1._wp - this%cost**2)

            this%tflag = .false.
            this%bounces = 0
            this%cnts = 0

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
            this%bounces = 0

            !teleport to just inside medium
            this%pos%z = state%grid%zmax - 1e-8_wp

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)

        end subroutine annulus
end module photonMod