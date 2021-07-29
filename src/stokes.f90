MODULE stokes_mod

implicit none

CONTAINS
   subroutine stokes(packet, grid)

   use constants, only : PI, TWOPI
   use random,    only : ran2   
   use photonMod
   use gridMod

   implicit none

   type(photon),    intent(INOUT) :: packet
   type(cart_grid), intent(IN)    :: grid
   
   real :: costp, sintp, phip, bmu, b, ri1, ri3, cosi3, sini3
   real :: cosb2, sinbt, cosi2, sini1, cosi1, sini2, bott, cosdph
   integer :: lay

   lay = packet%layer

   !***** isotropic scattering if g = 0.0 ******************************
   if(grid%hgg(lay) == 0.0) then
      packet%cost=2.*ran2()-1.
      packet%sint=(1. - packet%cost**2)
      if(packet%sint.le.0.)then
         packet%sint=0.
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

      bmu=((1.+grid%g2(lay))-((1.-grid%g2(lay))/(1.-grid%hgg(lay)+2.*grid%hgg(lay)*ran2()))**2)/(2.*grid%hgg(lay))
      cosb2=bmu**2
      b=cosb2-1.

      if(abs(bmu).gt.1.) then
         if(bmu.gt.1.) then
            bmu=1.
            cosb2=1.
            b=0.
         else
            bmu=-1.
            cosb2=1.
            b=0.
         end if
      end if
      sinbt=sqrt(1.-cosb2)
      ri1=TWOPI*ran2()

      if(ri1.gt.PI) then
         ri3=TWOPI-ri1
         cosi3=cos(ri3)
         sini3=sin(ri3)

         if(bmu.eq.1.) then
            goto 100
         else
            if(bmu.eq.-1.) then
               goto 100
            end if
         end if

         packet%cost=costp*bmu+sintp*sinbt*cosi3
         if(abs(packet%cost).lt.1.) then
            packet%sint=abs(sqrt(1. - packet%cost**2))
            sini2=sini3*sintp/packet%sint
            bott=packet%sint*sinbt
            cosi2=costp/bott-packet%cost*bmu/bott
         else
            packet%sint=0.
            sini2=0.
            if(packet%cost.ge.1.)  cosi2=-1.
            if(packet%cost.le.-1.) cosi2=1.
         end if

         cosdph=-cosi2*cosi3+sini2*sini3*bmu
         if(abs(cosdph).gt.1.) then
            if(cosdph.gt.1.) then
               cosdph=1.
            else
               cosdph=-1.
            end if
         end if
         
         packet%phi=phip+acos(cosdph)
         if(packet%phi.gt.TWOPI) packet%phi=packet%phi-TWOPI
         if(packet%phi.lt.0.)    packet%phi=packet%phi+TWOPI

         !      elseif(ri1.le.PI) then
      else  
         cosi1=cos(ri1)
         sini1=sin(ri1)
         if(bmu.eq.1.) then
            goto 100
         else
            if(bmu.eq.-1.) then
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
            packet%sint=0.
            sini2=0.
            if(packet%cost.ge.1.)  cosi2=-1.
            if(packet%cost.le.-1.) cosi2=1.
         end if

         cosdph=-cosi1*cosi2+sini1*sini2*bmu
         if(abs(cosdph).gt.1.) then
            if(cosdph.gt.1.) then
               cosdph=1.
            else
               cosdph=-1.
            end if
         end if
         packet%phi=phip-acos(cosdph)
         if(packet%phi.gt.TWOPI) packet%phi=packet%phi-TWOPI
         if(packet%phi.lt.0.)    packet%phi=packet%phi+TWOPI
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
