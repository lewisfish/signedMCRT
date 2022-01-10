module stokes_mod
!
! This module contains the scattering routines for isotropic and henyey-greenstein governed scattering
! Original is due to K. Woods
!

implicit none

contains
   subroutine stokes(packet, hgg, g2)

   use constants, only : PI, TWOPI, wp
   use random,    only : ran2   
   use photonMod

   implicit none

   type(photon),  intent(INOUT) :: packet
   real(kind=wp), intent(IN)    :: hgg, g2
   
   real(kind=wp) :: costp, sintp, phip, bmu, b, ri1, ri3, cosi3, sini3
   real(kind=wp) :: cosb2, sinbt, cosi2, sini1, cosi1, sini2, bott, cosdph

   !***** isotropic scattering if g = 0.0 ******************************
   if(hgg == 0.0_wp) then
      packet%cost=2._wp * ran2() - 1._wp
      packet%sint=(1._wp - packet%cost**2)
      if(packet%sint <= 0._wp)then
         packet%sint = 0._wp
      else
         packet%sint=sqrt(packet%sint)
      endif

      packet%phi=TWOPI*ran2()
      packet%sinp=sin(packet%phi)
      packet%cosp=cos(packet%phi)

      packet%nxp=packet%sint*packet%cosp
      packet%nyp=packet%sint*packet%sinp
      packet%nzp=packet%cost

   else

   !***** heyney greenstein scattering ********************************

      costp=packet%cost
      sintp=packet%sint
      phip=packet%phi

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

         packet%cost=costp*bmu+sintp*sinbt*cosi3
         if(abs(packet%cost) < 1._wp) then
            packet%sint=abs(sqrt(1._wp - packet%cost**2))
            sini2=sini3*sintp/packet%sint
            bott=packet%sint*sinbt
            cosi2=costp/bott-packet%cost*bmu/bott
         else
            packet%sint=0.
            sini2=0.
            if(packet%cost >= 1._wp)  cosi2=-1._wp
            if(packet%cost <= -1._wp) cosi2=1._wp
         end if

         cosdph=-cosi2*cosi3+sini2*sini3*bmu
         if(abs(cosdph) > 1._wp) then
            if(cosdph > 1._wp) then
               cosdph=1._wp
            else
               cosdph=-1._wp
            end if
         end if
         
         packet%phi=phip+acos(cosdph)
         if(packet%phi > TWOPI) packet%phi=packet%phi-TWOPI
         if(packet%phi.lt.0.)    packet%phi=packet%phi+TWOPI

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

         packet%cost=costp*bmu+sintp*sinbt*cosi1
         if(abs(packet%cost).lt.1.) then
            packet%sint=abs(sqrt(1. - packet%cost**2))
            sini2=sini1*sintp/packet%sint
            bott=packet%sint*sinbt
            cosi2=costp/bott-packet%cost*bmu/bott
         else
            packet%sint=0._wp
            sini2=0._wp
            if(packet%cost >= 1._wp)  cosi2=-1._wp
            if(packet%cost <= -1._wp) cosi2=1._wp
         end if

         cosdph=-cosi1*cosi2+sini1*sini2*bmu
         if(abs(cosdph) > 1._wp) then
            if(cosdph > 1._wp) then
               cosdph=1._wp
            else
               cosdph=-1._wp
            end if
         end if
         packet%phi=phip-acos(cosdph)
         if(packet%phi > TWOPI) packet%phi=packet%phi-TWOPI
         if(packet%phi < 0._wp)    packet%phi=packet%phi+TWOPI
      end if

      packet%cosp=cos(packet%phi)
      packet%sinp=sin(packet%phi)

      packet%nxp=packet%sint*packet%cosp
      packet%nyp=packet%sint*packet%sinp
      packet%nzp=packet%cost

   end if

   100   continue

   end subroutine stokes
end MODULE stokes_mod
