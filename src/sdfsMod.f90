module sdfs

    use vector_class

    implicit none


    type :: sphere
        real :: radius
        contains
        procedure :: evaluate => eval_sphere
    end type sphere


    type :: cylinder
        real :: height, radius
        contains
        procedure :: evaluate => eval_cylinder
    end type cylinder

    contains

        real function eval_cylinder(this, pos)

            implicit none

            class(cylinder) :: this
            type(vector), intent(IN) :: pos

            eval_cylinder = abs(cylinder_fn(pos, this%radius, this%height))

        end function eval_cylinder


        real function eval_sphere(this, pos)

            implicit none

            class(sphere) :: this
            type(vector), intent(IN) :: pos

            eval_sphere = sphere_fn(pos, this%radius)

        end function eval_sphere


        real function cylinder_fn(p, r, h)

            implicit none

            type(vector), intent(IN) :: p
            real,         intent(IN) :: h, r

            type(vector) :: d

            !lies along x plane
            d = abs(vector(sqrt(p%z**2 + p%y**2), p%x, 0.)) - vector(r, h, 0.)
            cylinder_fn = min(max(d%x, d%y), 0.) + length(max(d, 0.))

        end function cylinder_fn


        real function sphere_fn(p, c)

            implicit none

            type(vector), intent(IN) :: p
            real :: c

            sphere_fn = sqrt(p%x**2+p%y**2+p%z**2) - c

        end function sphere_fn

        type(vector) function translate_fn(position, offset)

            implicit none

            type(vector) :: position, offset

            translate_fn = position - offset

        end function translate_fn

end module sdfs