MODULE ch_opt

implicit none

CONTAINS
   
   pure subroutine init_opt1(kappa, albedo, hgg, g2)
!
!  subroutine to set tissue optical properties 630nm
!
   ! use opt_prop
   
   implicit none

   real, intent(OUT) :: kappa(:), albedo(:), hgg(:), g2(:)
   real :: mua(3), mus(3)


   !contents
   hgg(3) = .6
   g2(3)  = hgg(3)**2.
   mua(3) = 1d-17
   mus(3) = 5.!/(1.d0 - hgg)

   kappa(3) = mus(3) + mua(3)
   albedo(3) = 1.!mus(1) / kappa(1)

   !container
   hgg(2) = 0.
   g2(2)  = hgg(2)**2.
   mua(2) = 1d-17
   mus(2) = 5.d0!/(1.d0 - hgg)

   kappa(2) = mus(2) + mua(2)
   albedo(2) = 1.!mus(1) / kappa(1)
   !air
   hgg(1) = 0.
   g2(1)  = hgg(1)**2.
   mua(1) = 0.d0
   mus(1) = 0.d0!/(1.d0 - hgg)

   kappa(1) = mus(1) + mua(1)
   albedo(1) = 0.!mus(1) / kappa(1)

   end subroutine init_opt1
   
!    subroutine init_opt2
! !
! !  subroutine to set tissue optical properties 420nm
! !
!    use opt_prop

!    implicit none

!    hgg = 0.9
!    g2  = hgg**2.
!    mua = 1.8d0
!    mus = 82.d0/(1.d0 - hgg)

!    kappa  = mus + mua
!    albedo = mus / kappa

!    end subroutine init_opt2
   
!    subroutine init_opt3
! !
! !  subroutine to set tissue optical properties 705nm
! !
!    use opt_prop

!    implicit none

!    hgg = 0.9
!    g2  = hgg**2.
!    mua = 0.23d0
!    mus = 17.d0/(1.d0 - hgg)

!    kappa  = mus + mua 
!    albedo = mus / kappa

!    end subroutine init_opt3
   

   subroutine sample(array, size_of, cdf, wave, iseed)
!      
!  samples a random value from an array based upon its cdf     
!      
      implicit none
      
      integer, intent(IN)    :: iseed, size_of
      real,    intent(IN)    :: array(size_of, 2), cdf(size_of)
      real,    intent(OUT)   :: wave

      real :: ran2, value
      integer :: nlow
      
      value = ran2(iseed)
      
      call search_1D(size(cdf), cdf, nlow, value)
      call lin_inter_1D(array, cdf, value, size(cdf), nlow, wave)
   
   end subroutine sample
   
   subroutine lin_inter_1D(array, cdf, value, length, nlow, y)
!
!  linear interpolates between values for an array and its cdf
!   
      implicit none
   
      real,    intent(OUT)  :: y
      integer, intent(IN)   :: length
      real,    intent(IN)   :: value,array(length,2),cdf(length-1)
      integer, intent(IN)   :: nlow
   
      y = array(nlow+1,1) + (array(nlow+2,1) - array(nlow+1,1)) * (value - cdf(nlow))/(cdf(nlow+1) - cdf(nlow))
   
   end subroutine lin_inter_1D
   
   subroutine lin_inter_2D(array,value,length,nlow,y)
!
!  linear interpolation for an array
!
      implicit none

      real,    intent(OUT)  :: y
      integer, intent(IN)   :: length
      real,    intent(IN)   :: value,array(length,2)
      integer, intent(IN)   :: nlow
   
      y = array(nlow,2) + (array(nlow+1,2) - array(nlow,2)) * (value - array(nlow,1))/(array(nlow+1,1) - array(nlow,1))
   
   end subroutine lin_inter_2D
   
   subroutine search_1D(length,array,nlow,value)
!
!  search by bisection for 1D array
!
      implicit none
      
      integer              :: nup,length,middle
      integer, intent(OUT) :: nlow
      real,    intent(in)  :: array(length),value
      
      nup = length
      nlow = 1
      middle = int((nup+nlow)/2.)

      do while((nup - nlow).gt.1)
         middle = int((nup + nlow)/2.)
         if(value.gt.array(middle))then
            nlow = middle
         else
            nup = middle   
         end if
      end do
   end subroutine search_1D
   
   subroutine search_2D(length,array,nlow,value)
!
!  search by bisection for 2D array
!
      implicit none
      
      integer              :: nup,length,middle
      integer, intent(OUT) :: nlow
      real,    intent(in)  :: array(length,2),value
      
      nup = length
      nlow = 1
      middle = int((nup+nlow)/2.)

      do while((nup - nlow).gt.1)
         middle = int((nup + nlow)/2.)
         if(value.gt.array(middle,1))then
            nlow = middle
         else
            nup = middle   
         end if
      end do
   end subroutine search_2D
   
   subroutine mk_cdf(array,cdf,length)
!
!  subroutine that creates cdf for an array of values.
!
      implicit none

      integer, intent(IN)    :: length
      real,    intent(IN)    :: array(length,2)
      real,    intent(INOUT) :: cdf(length)
      real                   :: summ
      integer                :: i,j
   
      do j=1,length-1
         summ=0.
         do i=1,j   
            summ=summ+0.5*(array(i+1,2)+array(i,2))*(array(i+1,1)-array(i,1))
         end do
         cdf(j)=summ      
      end do
      cdf=cdf/cdf(length-1)
   
   end subroutine mk_cdf
end module ch_opt
