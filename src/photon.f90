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
        contains
            procedure :: scatter => scatter
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
                use sdfs,          only : rotationAlign, translate
                use mat_class,     only : invert
                use vector_class
            
                class(photon) :: this
                type(toml_table), optional, intent(inout) :: dict
                
                type(vector) :: a, b
                integer :: cell(3)
                real(kind=wp) :: t(4,4), radius, r, theta

                this%nxp = photon_origin%nxp
                this%nyp = photon_origin%nyp
                this%nzp = photon_origin%nzp

                call get_value(dict, "radius", radius)

                ! https://math.stackexchange.com/a/1681815
                r = radius * sqrt(ran2())
                theta = ran2() * TWOPI
                
                !set inital vector from which the source points
                a = vector(0._wp, 0._wp, 1._wp)
                a = a%magnitude()
                !set vector to rotate to. User defined.
                b = vector(this%nxp, this%nyp, this%nzp)
                b = b%magnitude()
                
                
                ! method fails if below condition. So change a vector to point down x-axis
                if(abs(a) == abs(b))then
                    a = vector(1._wp, 0._wp, 0._wp)
                    a = a%magnitude()
                    this%pos = vector(0._wp, radius * cos(theta), radius * sin(theta))
                else
                    this%pos = vector(radius * cos(theta), radius * sin(theta), 0._wp)
                end if

                ! get rotation matrix
                t = rotationAlign(a, b)
                ! get translation matrix
                t = t + translate(photon_origin%pos)
                ! transform point
                this%pos = this%pos .dot. t

                this%phi  = atan2(this%nyp, this%nxp)
                this%cosp = cos(this%phi)
                this%sinp = sin(this%phi)
                this%cost = this%nzp
                this%sint = sqrt(1._wp - this%cost**2)  
                
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

            use constants,     only : TWOPI 
            use utils,         only : deg2rad
            use tomlf,         only : toml_table, get_value
            use random,        only : ran2, rang
            use sim_state_mod, only : state

            class(photon) :: this
            type(toml_table), optional, intent(inout) :: dict

            character(len=:), allocatable :: beam_type
            real(kind=wp) :: beta, rlo, rhi, radius, tmp, mid, angle, x, y, z, phi, sinp, cosp
            type(vector)  :: pos
            integer       :: cell(3)

            call get_value(dict, "beta", beta)
            call get_value(dict, "radius", rlo)
            call get_value(dict, "radius_hi", rhi)
            call get_value(dict, "annulus_type", beam_type)

            if(beam_type == "tophat")then
                radius = rlo + (rhi - rlo) * sqrt(ran2())
            elseif(beam_type == "gaussian")then
                mid = (rhi - rlo) / 2.
                call rang(radius, tmp, mid, 0.04_wp)
            else
                error stop "No such beam type!"
            end if

            phi = TWOPI * ran2()

            angle = deg2rad(beta)

            cosp = cos(phi)
            sinp = sin(phi)
            x = radius * cosp
            y = radius * sinp
            z = state%grid%zmax - 1e-8_wp! just inside surface of medium. TODO make this user configurable?
            pos = vector(x, y, z)
            this%pos = pos

            this%nxp = sin(angle) * cosp
            this%nyp = sin(angle) * sinp
            this%nzp = -cos(angle)

            this%phi = phi
            this%cosp = cosp
            this%sinp = sinp
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
        end subroutine annulus

        subroutine scatter(this, hgg, g2)
        !This routine contains the scattering routines maths for isotropic and henyey-greenstein governed scattering
        ! Original is due to K. Woods
            use constants, only : PI, TWOPI, wp
            use random,    only : ran2   
         
            class(photon), intent(inout) :: this
            real(kind=wp), intent(in)    :: hgg, g2
            
            real(kind=wp) :: costp, sintp, phip, bmu, b, ri1, ri3, cosi3, sini3
            real(kind=wp) :: cosb2, sinbt, cosi2, sini1, cosi1, sini2, bott, cosdph
         
            !***** isotropic scattering if g = 0.0 ******************************
            if(hgg == 0.0_wp) then
               this%cost=2._wp * ran2() - 1._wp
               this%sint=(1._wp - this%cost**2)
               if(this%sint <= 0._wp)then
                  this%sint = 0._wp
               else
                  this%sint=sqrt(this%sint)
               endif
         
               this%phi=TWOPI*ran2()
               this%sinp=sin(this%phi)
               this%cosp=cos(this%phi)
         
               this%nxp=this%sint*this%cosp
               this%nyp=this%sint*this%sinp
               this%nzp=this%cost
         
            else
         
            !***** heyney greenstein scattering ********************************
         
               costp=this%cost
               sintp=this%sint
               phip=this%phi
         
               bmu=((1._wp+g2)-((1._wp-g2)/(1._wp-hgg+2._wp*hgg*ran2()))**2)/(2._wp*hgg)
               cosb2=bmu**2
               b=cosb2-1._wp
         
               if(abs(bmu) > 1._wp) then
                  if(bmu > 1._wp) then
                     bmu=1._wp
                     cosb2=1._wp
                     b=0._wp
                  else
                     bmu=-1._wp
                     cosb2=1._wp
                     b=0._wp
                  end if
               end if
               sinbt=sqrt(1._wp-cosb2)
               ri1=TWOPI*ran2()
         
               if(ri1 > PI) then
                  ri3=TWOPI-ri1
                  cosi3=cos(ri3)
                  sini3=sin(ri3)
         
                  if(bmu == 1._wp) then
                     goto 100
                  else
                     if(bmu == -1._wp) then
                        goto 100
                     end if
                  end if
         
                  this%cost=costp*bmu+sintp*sinbt*cosi3
                  if(abs(this%cost) < 1._wp) then
                     this%sint=abs(sqrt(1._wp - this%cost**2))
                     sini2=sini3*sintp/this%sint
                     bott=this%sint*sinbt
                     cosi2=costp/bott-this%cost*bmu/bott
                  else
                     this%sint=0.
                     sini2=0.
                     if(this%cost >= 1._wp)  cosi2=-1._wp
                     if(this%cost <= -1._wp) cosi2=1._wp
                  end if
         
                  cosdph=-cosi2*cosi3+sini2*sini3*bmu
                  if(abs(cosdph) > 1._wp) then
                     if(cosdph > 1._wp) then
                        cosdph=1._wp
                     else
                        cosdph=-1._wp
                     end if
                  end if
                  
                  this%phi=phip+acos(cosdph)
                  if(this%phi > TWOPI) this%phi=this%phi-TWOPI
                  if(this%phi.lt.0.)    this%phi=this%phi+TWOPI
         
                  !      elseif(ri1 <= PI) then
               else  
                  cosi1=cos(ri1)
                  sini1=sin(ri1)
                  if(bmu == 1.) then
                     goto 100
                  else
                     if(bmu == -1._wp) then
                        goto 100
                  end if
               end if
         
                  this%cost=costp*bmu+sintp*sinbt*cosi1
                  if(abs(this%cost).lt.1.) then
                     this%sint=abs(sqrt(1. - this%cost**2))
                     sini2=sini1*sintp/this%sint
                     bott=this%sint*sinbt
                     cosi2=costp/bott-this%cost*bmu/bott
                  else
                     this%sint=0._wp
                     sini2=0._wp
                     if(this%cost >= 1._wp)  cosi2=-1._wp
                     if(this%cost <= -1._wp) cosi2=1._wp
                  end if
         
                  cosdph=-cosi1*cosi2+sini1*sini2*bmu
                  if(abs(cosdph) > 1._wp) then
                     if(cosdph > 1._wp) then
                        cosdph=1._wp
                     else
                        cosdph=-1._wp
                     end if
                  end if
                  this%phi=phip-acos(cosdph)
                  if(this%phi > TWOPI) this%phi=this%phi-TWOPI
                  if(this%phi < 0._wp)    this%phi=this%phi+TWOPI
               end if
         
               this%cosp=cos(this%phi)
               this%sinp=sin(this%phi)
         
               this%nxp=this%sint*this%cosp
               this%nyp=this%sint*this%sinp
               this%nzp=this%cost
         
            end if
         
            100   continue
         
            end subroutine scatter
end module photonMod