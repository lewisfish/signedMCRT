module testy

implicit none

    contains

        real function dip_map(p)
            
            use vector_class

            implicit none

            type(vector), intent(IN) :: p

            dip_map = min(0.0, 0.09*cos(3.14*10.0*p%x)*cos(3.14*10.0*p%y)*cos(3.14*10.0*p%z));

        end function dip_map
end module testy

program test_sdf
    
    use sdfs
    use testy
    use vector_class

    implicit none

    type(container), allocatable :: array(:)
    
    type(torus), target :: tor
    type(box),   target :: boxy

    type(displacement), target :: dp
    type(elongate),     target :: ep
    type(twist),        target :: tp
    type(bend),         target :: bp

    allocate(array(1))

    boxy = box(1., 1., 2., 3., 4., 1)
    tor = torus(.5, .2, 1., 2., 3., 4., 1)

    bp = bend(boxy, 1.)
    dp = displacement(tor, dip_map)
    ep = elongate(tor, vector(0., 0., 1.))
    tp = twist(boxy, 3.14/2.)

    allocate(array(1)%p, source=bp)
    array(1)%p => bp
    call render(array, vector(1., 1., 1.), 200, "../bend.dat")

    allocate(array(1)%p, source=dp)
    array(1)%p => dp
    call render(array, vector(1., 1., 1.), 200, "../displacement.dat")

    allocate(array(1)%p, source=ep)
    array(1)%p => ep
    call render(array, vector(2.5, 2.5, 2.5), 200, "../elongate.dat")

    allocate(array(1)%p, source=tp)
    array(1)%p => tp
    call render(array, vector(1., 1., 1.), 200, "../twist.dat")

end program test_sdf