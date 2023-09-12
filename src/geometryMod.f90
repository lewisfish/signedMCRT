module geometry
!!    Defines a set of functions for intersecting a line and a surface.
!!
!!    - Circle
!!    - Plane
!!    - Cone
!!    - Cylinder
!!    - Ellipse
!!    - Sphere


    use vector_class, only : vector
    use constants,    only : wp

    implicit none

    private
    public :: intersectCircle, intersectPlane, intersectCone, intersectCylinder, intersectEllipse, intersectSphere

    contains

    logical function intersectSphere(orig, dir, t, centre, radius)
    ! calculates where a line, with origin:orig and direction:dir hits a sphere, centre:centre and radius:radius
    ! returns true if intersection exists
    ! returns t, the paramertised parameter of the line equation
    ! adapted from scratchapixel

        type(vector),  intent(IN)  :: dir, orig, centre
        real(kind=wp), intent(OUT) :: t
        real(kind=wp), intent(IN)  :: radius

        type(vector)  :: L
        real(kind=wp) :: t0, t1, a, b, c, tmp

        intersectSphere = .false.

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
        intersectSphere = .true.
        return

    end function intersectSphere

    logical function intersectCylinder(orig, dir, t, centre, radius)
    ! calculates where a line, with origin:orig and direction:dir hits a cylinder, centre:centre and radius:radius
    ! returns true if intersection exists
    ! returns t, the paramertised parameter of the line equation
    ! adapted from scratchapixel
    ! need to check z height after moving ray
    ! if not this is an infinte cylinder
    ! cylinder lies length ways along z-axis

        type(vector),  intent(IN)  :: dir, orig, centre
        real(kind=wp), intent(OUT) :: t
        real(kind=wp), intent(IN)  :: radius

        type(vector)  :: L
        real(kind=wp) :: t0, t1, a, b, c, tmp

        intersectCylinder = .false.

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
        intersectCylinder = .true.
        return
    end function intersectCylinder


    logical function intersectEllipse(orig, dir, t, centre, semia, semib)
    ! calculates where a line, with origin:orig and direction:dir hits a ellipse, centre:centre and axii:semia, semib
    ! returns true if intersection exists
    ! returns t, the paramertised parameter of the line equation
    ! adapted from scratchapixel and pbrt
    ! need to check z height after moving ray
    ! if not this is an infinte ellipse-cylinder
    ! ellipse lies length ways along z-axis
    ! semia and semib are the semimajor axis which are the half width and height.

        type(vector),  intent(IN)  :: dir, orig, centre
        real(kind=wp), intent(OUT) :: t
        real(kind=wp), intent(IN)  :: semia, semib

        type(vector)  :: L
        real(kind=wp) :: t0, t1, a, b, c, tmp, semia2div, semib2div

        intersectEllipse = .false.

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
        intersectEllipse = .true.
        return
    end function intersectEllipse


    logical function intersectCone(orig, dir, t, centre, radius, height)
    ! calculates where a line, with origin:orig and direction:dir hits a cone, radius:radius and height:height with centre:centre.
    ! centre is the point under the apex at the cone's base.
    ! returns true if intersection exists
    ! returns t, the paramertised parameter of the line equation
    ! adapted from scratchapixel and pbrt
    ! need to check z height after moving ray
    ! if not this is an infinte cone
    ! cone lies height ways along z-axis

        type(vector),  intent(IN)  :: orig, dir, centre
        real(kind=wp), intent(IN)  :: radius, height
        real(kind=wp), intent(OUT) :: t

        type(vector)  :: L
        real(kind=wp) :: t0, t1, a, b, c, tmp, k

        intersectCone = .false.
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

        intersectCone = .true.
        return

    end function intersectCone

    logical function intersectPlane(n, p0, l0, l, t)
    !https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
        type(vector),  intent(in) :: n, p0, l, l0
        real(kind=wp), intent(inout) :: t

        real(kind=wp) :: denom
        type(vector)  :: p0l0

        intersectPlane = .false.
        denom = l .dot. n
        if(denom > 1e-6_wp)then
            p0l0 = p0 - l0
            t = p0l0 .dot. n
            t = t / denom
            if(t >= 0._wp)intersectPlane=.true.
        end if
    end function intersectPlane


    logical function intersectCircle(n, p0, radius, l0, l, t)
    !https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
        type(vector),  intent(in) :: n, p0, l, l0
        real(kind=wp), intent(in) :: radius

        real(kind=wp) :: t, d2
        type(vector) :: v, p

        intersectCircle = .false.
        t = 0._wp
        if(intersectPlane(n, p0, l0, l, t))then
            p = l0 + l * t
            v = p - p0
            d2 = v .dot. v
            if(sqrt(d2) <= radius)intersectCircle=.true.
        end if
    end function intersectCircle

    logical function solveQuadratic(a, b, c, x0, x1)
    ! solves quadratic equation given coeffs a, b, and c
    ! returns true if real soln
    ! returns x0 and x1
    ! adapted from scratchapixel

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
end module geometry