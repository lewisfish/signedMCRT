module surfaces

    use vector_class, only : vector
    use constants,    only : wp

    implicit none

    private
    public :: reflect_refract

    contains

    subroutine reflect_refract(I, N, n1, n2, rflag, ri)
    ! wrapper routine for fresnel calculation
    !
    !
        use random, only : ran2

        type(vector),  intent(INOUT) :: I !incident vector
        type(vector),  intent(INOUT) :: N ! normal vector
        real(kind=wp), intent(IN)    :: n1, n2 !refractive indcies
        real(kind=wp), intent(OUT)   :: Ri
        logical,       intent(OUT)   :: rflag !reflection flag

        rflag = .FALSE.

        !draw random number, if less than fresnel coefficents, then reflect, else refract
        Ri = fresnel(I, N, n1, n2)
        if(ran2() <= Ri)then
            call reflect(I, N)
            rflag = .true.
        else
            call refract(I, N, n1/n2)
        end if

    end subroutine reflect_refract

    subroutine reflect(I, N)
    !   get vector of reflected photon
    !
    !

        type(vector), intent(INOUT) :: I ! incident vector
        type(vector), intent(IN)    :: N ! normal vector

        type(vector) :: R

        R = I - 2._wp * (N .dot. I) * N
        I = R

    end subroutine reflect

    subroutine refract(I, N, eta)
    !   get vector of refracted photon
    !
    !

        type(vector),  intent(INOUT) :: I
        type(vector),  intent(IN)    :: N
        real(kind=wp), intent(IN)    :: eta

        type(vector)  :: T, Ntmp
        real(kind=wp) :: c1, c2

        Ntmp = N

        c1 = (Ntmp .dot. I)
        if(c1 < 0._wp)then
            c1 = -c1
        else
            Ntmp = (-1._wp) * N
        end if
        c2 = sqrt(1._wp - (eta)**2 * (1._wp-c1**2))

        T = eta*I + (eta * c1 - c2) * Ntmp 

        I = T

    end subroutine refract

    function fresnel(I, N, n1, n2) result (tir)
    !   calculates the fresnel coefficents
    !
    !
        use ieee_arithmetic, only : ieee_is_nan

        real(kind=wp), intent(IN) :: n1, n2
        type(vector),  intent(IN) :: I, N

        real(kind=wp) :: costt, sintt, sint2, cost2, tir, f1, f2

        costt = abs(I .dot. N)

        sintt = sqrt(1._wp - costt * costt)
        sint2 = n1/n2 * sintt
        if(sint2 > 1._wp)then
            tir = 1.0_wp
            return
        elseif(costt == 1._wp)then
            tir = 0._wp
            return
        else
            sint2 = (n1/n2)*sintt
            cost2 = sqrt(1._wp - sint2 * sint2)
            f1 = abs((n1*costt - n2*cost2) / (n1*costt + n2*cost2))**2
            f2 = abs((n1*cost2 - n2*costt) / (n1*cost2 + n2*costt))**2

            tir = 0.5_wp * (f1 + f2)
        if(ieee_is_nan(tir) .or. tir > 1._wp .or. tir < 0._wp)print*,'TIR: ', tir, f1, f2, costt,sintt,cost2,sint2
            return
        end if
    end function fresnel
end module surfaces