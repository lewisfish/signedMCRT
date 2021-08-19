module testy

implicit none

    contains

        real function dip_map(p)
            
            use vector_class

            implicit none

            type(vector), intent(IN) :: p

            dip_map = min(0.0, 0.09*cos(3.14*10.0*p%x)*cos(3.14*10.0*p%y)*cos(3.14*10.0*p%z));

        end function dip_map


        subroutine draw(array, source, name, extent, samples)
            
            use sdfs
            use vector_class

            implicit none
        
            type(container) :: array(1)
            class(sdf),   target, intent(IN) :: source
            character(*),         intent(IN) :: name
            type(vector), optional, intent(IN) :: extent
            integer,      optional, intent(IN) :: samples

            integer      :: s
            type(vector) :: ext

            if(present(samples))then
                s = samples
            else
                s = 200
            end if
            if(present(extent))then
                ext = extent
            else
                ext = vector(1., 1., 1.)
            end if

            allocate(array(1)%p, source=source)
            array(1)%p => source
            call render(array, ext, s, name)
        
        end subroutine draw

end module testy

program test_sdf
    
    use sdfs
    use testy
    use vector_class

    implicit none

    type(container), allocatable :: array(:), cnta(:)
    
    type(torus),     target :: tor
    type(box),       target :: boxy, sboxy
    type(cone),      target :: cony
    type(triprisim), target :: tri
    type(cylinder),  target :: cyl
    type(sphere),    target :: sph
    type(capsule),   target :: cap
    type(plane),     target :: plan

    type(displacement), target :: dp
    type(elongate),     target :: ep
    type(repeat),       target :: rp
    type(twist),        target :: tp
    type(bend),         target :: bp

    type(model), target :: mod

    allocate(array(1), cnta(2))

    boxy = box(1., 1., 2., 3., 4., 1)
    sboxy = box(vector(.4, .2, .2), 1., 2., 3., 4., 1)
    tor = torus(.5, .2, 1., 2., 3., 4., 1)
    cony = cone(vector(0., 0., -.5), vector(0., 0., .5), 1., 0.0, 2., 3., 4., 1., 1)
    tri = triprisim(1., .05, 1., 2., 3., 4., 1)
    sph = sphere(1., 1., 2., 3., 4., 1)
    cyl = cylinder(vector(0., 0., 0.), vector(0., 0., .5), .2, 1., 2., 3., 4., 1)
    cap = capsule(vector(0., 0., 0.), vector(0., 0., 0.5), .2, 1., 2., 3., 4., 1)
    plan = plane(vector(0., 0., 1.), 1., 2., 3., 4., 1)

    bp = bend(boxy, 1.)
    dp = displacement(tor, dip_map)
    ep = elongate(tor, vector(0., 0., 1.))
    tp = twist(boxy, 3.14/2.)
    rp = repeat(sboxy, 1., vector(-1., -2., -2.), vector(1., 2., 2.))

    call draw(array, cony, "../cone.dat")
    call draw(array, tri, "../triprisim.dat")
    call draw(array, tor, "../torus.dat")
    call draw(array, boxy, "../box.dat")
    call draw(array, sph, "../sphere.dat")
    call draw(array, cyl, "../cylinder.dat")
    call draw(array, cap, "../capsule.dat")
    call draw(array, plan, "../plane.dat")

    call draw(array, bp, "../bend.dat")
    call draw(array, dp, "../displacement.dat")
    call draw(array, ep, "../elongate.dat", vector(2.5, 2.5, 2.5))
    call draw(array, tp, "../twist.dat")
    call draw(array, rp, "../repeat.dat", vector(2., 4., 4.))


    sph = sphere(.5, 1., 2., 3., 4., 1)
    boxy = box(vector(1.1, 1.5, .5), 1., 2., 3., 4., 1)
    allocate(cnta(2)%p, source=boxy)
    allocate(cnta(1)%p, source=sph)
    cnta(2)%p => boxy
    cnta(1)%p => sph

    mod = model_init(cnta, union)
    call draw(array, mod, "../union.dat")
    mod = model_init(cnta, smoothunion, .5)
    call draw(array, mod, "../smoothunion.dat")
    mod = model_init(cnta, subtraction)
    call draw(array, mod, "../subtraction.dat")
    mod = model_init(cnta, intersection)
    call draw(array, mod, "../intersection.dat")

end program test_sdf

! 16

! 4x4

! cone, triprisim, torus, box
! sphere, cylinder, capsule, plane
! union(smooth), inter, sub, repeat
! bend, twist, disp, elongate