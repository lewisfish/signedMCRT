module photonMod
    
    use constants, only : wp
    use vector_class
    use dict_mod

    implicit none
    
    type :: photon
    
        type(vector)  :: pos                     ! position
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
            
            use dict_mod

            import :: photon
            class(photon) :: this
            type(dict_t), optional, intent(IN) :: dict

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
            ! elseif(choice == "annulus")then
            !     init_source%emit => annulus
            elseif(choice == "focus")then
                init_source%emit => focus
            elseif(choice == "point")then
                init_source%emit => point
            ! else
            !     init_source%emit => circular_beam
            end if


        end function init_source

        subroutine point(this, dict)
        !isotropic point source

            use sim_state_mod, only : state
            use random,        only : ran2
            use constants,     only : twoPI

            implicit none
        
            class(photon) :: this
            type(dict_t), optional, intent(IN) :: dict
            integer :: cell(3)

            this%pos%x = dict%get_value_real("pos%x")
            this%pos%y = dict%get_value_real("pos%y")
            this%pos%z = dict%get_value_real("pos%z")

            this%phi  = ran2()*twoPI
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)
            this%cost = 2._wp*ran2()-1._wp
            this%sint = sqrt(1._wp - this%cost**2)

            this%nxp = this%sint * this%cosp
            this%nyp = this%sint * this%sinp
            this%nzp = this%cost

            this%tflag = .false.
            this%cnts = 0
            this%layer=1

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

            implicit none

            class(photon) :: this
            type(dict_t), optional, intent(IN) :: dict

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

            use random,        only : ranu, ran2, randint
            use sim_state_mod, only : state

            implicit none

            class(photon) :: this
            type(dict_t), optional, intent(IN) :: dict

            integer :: cell(3)
            character(len=2) :: dir

            dir = dict%get_value_str("dir")

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

            implicit none

            class(photon) :: this
            type(dict_t), optional, intent(IN) :: dict

            integer :: cell(3)

            this%pos%z = state%grid%zmax - epsilon(1._wp)
            this%pos%x = ranu(-state%grid%xmax/10._wp, state%grid%xmax/10._wp)
            this%pos%y = ranu(-state%grid%ymax/10._wp, state%grid%ymax/10._wp)

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
        
        ! subroutine annulus(this, dict)

        !     use gridMod
        !     use utils,     only : deg2rad
        !     use random,    only : ranu, ran2, rang
        !     use surfaces,  only : reflect_refract

        !     implicit none

        !     class(photon) :: this
        !     type(dict_t), optional, intent(IN) :: dict


        !     real :: ra, ta, alpha, Ls, beam_width, rp, da, dist
        !     real :: tana, x, y, total_length, axicon_n, height, k
        !     type(vector) :: pos, dir, normal, centre, newpos
        !     logical :: flag

        !     ra = 4. ! 12.7mm
        !     beam_width = .5d0 !100um
        !     do
        !         call rang(x, y, 0., beam_width)
        !         rp = sqrt(x**2 + y**2)
        !         axicon_n = 1.45
        !         ta = 0.0d-3  ! 0.0 mm
        !         alpha = deg2rad(20.0)
        !         tana = tan(alpha)
        !         height = tana * ra
        !         Ls = 10.  ! 100mm
        !         centre = vector(0., 0., 0.)
        !         k  = (rp / height)**2
                
        !         da = (ra - rp) * tana
        !         pos = vector(x, y, -da)

        !         total_length = Ls + height

        !         dir = vector(0., 0., -1.)
        !         !https://math.stackexchange.com/a/3349448
        !         normal =vector(2.*pos%x, 2.*pos%y, -2*ra**2*(pos%z+height)/height**2)
        !         normal = normal%magnitude()
        !         call reflect_refract(dir, normal, axicon_n, 1., flag)
        !         if(.not. flag)then
        !             exit
        !         end if
        !     end do

        !     !move packet to L1 location
        !     pos = pos + dir * (-total_length / dir%z)

        !     !focus to a point
        !     pos%z = 4.e0 !f distance 40mm
        !     call rang(x, y, 0., 5.d-4)

        !     newpos = vector(x, y, dict%get_value_real("focus"))!.775
        !     dist = sqrt((pos%x - newpos%x)**2 + (pos%y - newpos%y)**2 +(pos%z - newpos%z)**2)

        !     dir = (-1.)*pos / dist
        !     dist = (grid%zmax + 5e-8 - pos%z) / dir%z
        !     pos = pos + dist * dir
        !     this%pos = pos

        !     dir = dir%magnitude()

        !     this%nxp = dir%x
        !     this%nyp = dir%y
        !     this%nzp = dir%z

        !     this%cost = this%nzp
        !     this%sint = sqrt(1.e0 - this%cost**2)

        !     this%phi = atan2(this%nyp, this%nxp)
        !     this%cosp = cos(this%phi)
        !     this%sinp = sin(this%phi)

        !     this%tflag = .false.
        !     this%layer = 3
        !     this%cnts = 0

        !     !teleport to just inside medium
        !     this%pos%z = grid%zmax - 1e-8

        !     ! Linear Grid 
        !     this%xcell=int(grid%nxg*(this%pos%x+grid%xmax)/(2.*grid%xmax))+1
        !     this%ycell=int(grid%nyg*(this%pos%y+grid%ymax)/(2.*grid%ymax))+1
        !     this%zcell=int(grid%nzg*(this%pos%z+grid%zmax)/(2.*grid%zmax))+1


        ! end subroutine annulus
        
        ! subroutine circular_beam(this dict)

        !     use random,    only : ranu, rang
        !     use surfaces,  only : intersect_cone, reflect_refract
        !     use vector_class

        !     implicit none

        !     class(photon) :: this
        !     type(dict_t), optional,  intent(IN) :: dict


        !     real :: seperation, beam_width, radius, height, alpha, axicon_n, base_pos
        !     real :: posx, posy, t, k
        !     type(vector) :: centre, pos, dir, normal
        !     logical :: flag

        !     seperation = 2.5d-3
        !     !beam paramaters
        !     beam_width = 100d-6
        !     !axicon paramaters
        !     radius = 12.7d-3
        !     height = 1.1d-3
        !     alpha = atan(height / radius)
        !     k  = (radius / height)**2
        !     axicon_n = 1.45
        !     base_pos = grid%zmax + ((seperation + beam_width) / tan(alpha * (axicon_n -1.)))
        !     centre = vector(0., 0., base_pos)

        !     call rang(posx, posy, 0., beam_width)
        !     pos = centre + vector(posx, posy, 2*height+base_pos)!2*height as want the upper cone
        !     dir = vector(0., 0., -1.)
        !     ! cartesian equation defines these two cones. usually want the lower cone
        !     ! in this case we want the upper cone as the axicon points down
        !     ! therefore need to change the computation of the normals and initial postion of packet 
        !     !
        !     ! \      /
        !     !  \    /   upper cone
        !     !   \  /
        !     !    \/
        !     !    /\
        !     !   /  \
        !     !  /    \   lower cone
        !     ! /      \

        !     flag = intersect_cone(pos, dir, t, centre, radius, height)
        !     if(flag)then
        !         pos = pos + t*dir
        !         if(pos%z >= centre%z+height)then! >= for upper cone
        !             !derivative of the cartesian cone eqn
        !             normal = vector(2*(pos%x-centre%x) / k, 2*(pos%y-centre%y) / k, -2*(pos%z-centre%z)+2*height)
        !             normal = normal *(-1.)! upper cone so invert normals
        !             normal = normal%magnitude()
        !             call reflect_refract(dir, normal, axicon_n, 1., flag)
        !             ! move to mediums surface and step inside
        !             t = ((grid%zmax- epsilon(1.e0)) - pos%z) / dir%z
        !             pos = pos + t*dir
        !         end if
        !     end if

        !     this%pos = pos

        !     this%nxp = dir%x
        !     this%nyp = dir%y
        !     this%nzp = dir%z

        !     this%cost = this%nzp
        !     this%sint = sqrt(1.e0 - this%cost**2)

        !     this%phi = atan2(this%nyp, this%nxp)
        !     this%cosp = cos(this%phi)
        !     this%sinp = sin(this%phi)

        !     this%tflag = .false.

        !     ! Linear Grid 
        !     this%xcell=int(grid%nxg*(this%pos%x+grid%xmax)/(2.*grid%xmax))+1
        !     this%ycell=int(grid%nyg*(this%pos%y+grid%ymax)/(2.*grid%ymax))+1
        !     this%zcell=int(grid%nzg*(this%pos%z+grid%zmax)/(2.*grid%zmax))+1


        ! end subroutine circular_beam

        ! subroutine bessel(this, grid)

        !     use gridMod,   only : cart_grid
        !     use constants, only : PI, twoPI
        !     use random,    only : ranu, rang

        !     implicit none

        !     class(photon) :: this
        !     type(cart_grid), intent(IN)    :: grid

        !     real :: tana, r_pos, x0, y0, z0, dist, n


        !     real :: waist, d, raxi, alpha

        !     waist = .5d-1
        !     alpha = 5.
        !     raxi = 12.7
        !     d = 10.e0
        !     this%wavelength = 488d-9
        !     this%fact = twopi/this%wavelength

        !     this%phase = 0.e0
        !     tana = tan(alpha *pi/180.)

        !     call rang(this%xp, this%yp, 0.e0, sqrt(2.)*waist/4.)
        !     r_pos = sqrt(this%xp**2 + this%yp**2)

        !     this%zp = (raxi - r_pos) * tana
        !     this%phase = this%zp * n

        !     x0 = ranu(-grid%xmax, grid%xmax)
        !     y0 = ranu(-grid%xmax, grid%xmax)
        !     z0 = (r_pos * tana) + d + this%zp
            

        !     dist = sqrt((x0 - this%xp)**2 + (y0 - this%yp)**2 + (z0 - this%zp)**2)

        !     this%phase = this%phase + dist

        !     this%nxp = (x0 - this%xp) / dist
        !     this%nyp = (y0 - this%yp) / dist
        !     this%nzp = -(z0 - this%zp) / dist ! -ve due to way z pos is defined

        !     this%cost = this%nzp
        !     this%sint = sqrt(1.e0 - this%cost**2)

        !     this%phi = atan2(this%nyp, this%nxp)
        !     this%cosp = cos(this%phi)
        !     this%sinp = sin(this%phi)

        !     this%xp = x0
        !     this%yp = y0
        !     this%zp = grid%zmax - grid%delta
        !     this%tflag = .false.
        !     this%xcell = int(grid%nxg * (this%xp + grid%xmax) / (2. * grid%xmax)) + 1
        !     this%ycell = int(grid%nyg * (this%yp + grid%ymax) / (2. * grid%ymax)) + 1
        !     this%zcell = int(grid%nzg * (this%zp + grid%zmax) / (2. * grid%zmax)) + 1

        ! end subroutine bessel

end module photonMod