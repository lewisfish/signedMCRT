module geometry

    use vector_class, only : vector
    use constants,    only : wp

    implicit none

    private
    public :: intersectCircle, intersectPlane

    contains

    logical function intersectPlane(n, p0, l0, l, t)
    !https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
        type(vector),  intent(in) :: n, p0, l, l0
        real(kind=wp), intent(inout) :: t

        real(kind=wp) :: denom
        type(vector)  :: p0l0

        intersectPlane = .false.
        denom = n .dot. l
        if(denom > 1e-6_wp)then
            p0l0 = p0 - l0
            t = p0l0 .dot. n
            t = t / denom
            if(t >= 0._wp)intersectPlane=.true.
        end if
    end function intersectPlane


    logical function intersectCircle(n, p0, radius, l0, l)
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
end module geometry