module surfaces

    use vector_class
    use constants,    only : wp

    implicit none

    private
    public :: intersect_sphere, intersect_cylinder, intersect_ellipse, intersect_cone
    public :: reflect_refract

    contains

    logical function intersect_sphere(orig, dir, t, centre, radius)
    ! calculates where a line, with origin:orig and direction:dir hits a sphere, centre:centre and radius:radius
    ! returns true if intersection exists
    ! returns t, the paramertised parameter of the line equation
    ! adapted from scratchapixel
        
        implicit none

        type(vector),  intent(IN)  :: dir, orig, centre
        real(kind=wp), intent(OUT) :: t
        real(kind=wp), intent(IN)  :: radius

        type(vector)  :: L
        real(kind=wp) :: t0, t1, a, b, c, tmp

        intersect_sphere = .false.

        L = orig - centre
        a = dir .dot. dir
        b = 2._wp * (dir .dot. L)
        c = (l .dot. l) - radius**2

        if(.not. solveQuadratic(a, b, c, t0, t1))return
        if(t0 > t1)then
            tmp = t1
            t1 = t0
            t0 = tmp
        end if
        if(t0 < 0._wp)then
            t0 = t1
            if(t0 < 0._wp)return
        end if

        t = t0
        intersect_sphere = .true.
        return

    end function intersect_sphere

    logical function intersect_cylinder(orig, dir, t, centre, radius)
    ! calculates where a line, with origin:orig and direction:dir hits a cylinder, centre:centre and radius:radius
    ! returns true if intersection exists
    ! returns t, the paramertised parameter of the line equation
    ! adapted from scratchapixel
    ! need to check z height after moving ray
    ! if not this is an infinte cylinder
    ! cylinder lies length ways along z-axis
        
        implicit none

        type(vector),  intent(IN)  :: dir, orig, centre
        real(kind=wp), intent(OUT) :: t
        real(kind=wp), intent(IN)  :: radius

        type(vector)  :: L
        real(kind=wp) :: t0, t1, a, b, c, tmp

        intersect_cylinder = .false.

        L = orig - centre
        a = dir%z**2 + dir%y**2
        b = 2._wp * (dir%z * L%z + dir%y * L%y)
        c = L%z**2 + L%y**2 - radius**2

        if(.not. solveQuadratic(a, b, c, t0, t1))return
        if(t0 > t1)then
            tmp = t1
            t1 = t0
            t0 = tmp
        end if
        if(t0 < 0._wp)then
            t0 = t1
            if(t0 < 0._wp)return
        end if

        t = t0
        intersect_cylinder = .true.
        return
    end function intersect_cylinder


    logical function intersect_ellipse(orig, dir, t, centre, semia, semib)
    ! calculates where a line, with origin:orig and direction:dir hits a ellipse, centre:centre and axii:semia, semib
    ! returns true if intersection exists
    ! returns t, the paramertised parameter of the line equation
    ! adapted from scratchapixel and pbrt
    ! need to check z height after moving ray
    ! if not this is an infinte ellipse-cylinder
    ! ellipse lies length ways along z-axis
    ! semia and semib are the semimajor axis which are the half width and height.
        
        implicit none

        type(vector),  intent(IN)  :: dir, orig, centre
        real(kind=wp), intent(OUT) :: t
        real(kind=wp), intent(IN)  :: semia, semib

        type(vector)  :: L
        real(kind=wp) :: t0, t1, a, b, c, tmp, semia2div, semib2div

        intersect_ellipse = .false.

        semia2div = 1._wp / semia**2
        semib2div = 1._wp / semib**2

        L = orig - centre
        a = semia2div * dir%z**2 + semib2div * dir%y**2
        b = 2._wp * (semia2div * dir%z * L%z + semib2div * dir%y * L%y)
        c = semia2div * L%z**2 + semib2div * L%y**2 - 1._wp

        if(.not. solveQuadratic(a, b, c, t0, t1))return
        if(t0 > t1)then
            tmp = t1
            t1 = t0
            t0 = tmp
        end if
        if(t0 < 0._wp)then
            t0 = t1
            if(t0 < 0._wp)return
        end if

        t = t0
        intersect_ellipse = .true.
        return
    end function intersect_ellipse


    logical function intersect_cone(orig, dir, t, centre, radius, height)
    ! calculates where a line, with origin:orig and direction:dir hits a cone, radius:radius and height:height with centre:centre.
    ! centre is the point under the apex at the cone's base.
    ! returns true if intersection exists
    ! returns t, the paramertised parameter of the line equation
    ! adapted from scratchapixel and pbrt
    ! need to check z height after moving ray
    ! if not this is an infinte cone
    ! cone lies height ways along z-axis

        implicit none

        type(vector),  intent(IN)  :: orig, dir, centre
        real(kind=wp), intent(IN)  :: radius, height
        real(kind=wp), intent(OUT) :: t

        type(vector)  :: L
        real(kind=wp) :: t0, t1, a, b, c, tmp, k


        intersect_cone = .false.
        k = radius / height
        k = k**2

        L = orig - centre
        a = dir%x**2 + dir%y**2 - (k*dir%z**2)
        b = 2._wp*((dir%x * L%x) + (dir%y * L%y) - (k*dir%z * (L%z - height)))
        c = L%x**2 + L%y**2 - (k*(L%z - height)**2)

        if(.not. solveQuadratic(a, b, c, t0, t1))return
        if(t0 > t1)then
            tmp = t1
            t1 = t0
            t0 = tmp
        end if
        if(t0 < 0._wp)then
            t0 = t1
            if(t0 < 0._wp)return
        end if

        t = t0

        intersect_cone = .true.
        return

    end function intersect_cone


    logical function solveQuadratic(a, b, c, x0, x1)
    ! solves quadratic equation given coeffs a, b, and c
    ! returns true if real soln
    ! returns x0 and x1
    ! adapted from scratchapixel

        implicit none

        real(kind=wp), intent(IN)  :: a, b, c
        real(kind=wp), intent(OUT) :: x0, x1

        real(kind=wp) :: discrim, q

        solveQuadratic = .false.

        discrim = b**2 - 4._wp * a * c
        if(discrim < 0._wp)then
            return
        elseif(discrim == 0._wp)then
            x0 = -0.5_wp*b/a
            x1 = x0
        else
            if(b > 0._wp)then
                q = -0.5_wp * (b + sqrt(discrim))
            else
                q = -0.5_wp * (b - sqrt(discrim))
            end if
            x0 = q / a
            x1 = c / q
        end if
        solveQuadratic = .true.
        return

    end function solveQuadratic



    subroutine reflect_refract(I, N, n1, n2, rflag, ri)
    ! wrapper routine for fresnel calculation
    !
    !
        use random, only : ran2

        implicit none

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

        implicit none

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

        implicit none

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

        implicit none

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