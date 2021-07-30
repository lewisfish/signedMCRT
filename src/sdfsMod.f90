module sdfs

    use vector_class

    implicit none

    type :: sdf
        contains
        procedure :: evaluate => evaluate_fn
    end type sdf

    type :: container
        class(sdf), pointer :: p => null()
    end type container

    type, extends(sdf) :: sphere
        real :: radius
        contains
        procedure :: evaluate => eval_sphere
    end type sphere


    type, extends(sdf) :: cylinder
        real         :: height, radius
        type(vector) :: centre
        contains
        procedure :: evaluate => eval_cylinder
    end type cylinder

    type, extends(sdf) :: box
        real :: length
        contains
        procedure :: evaluate => eval_box
    end type box


    type, extends(sdf) :: model
        type(container), allocatable :: array(:)
        procedure(op),nopass, pointer :: func
        contains
        procedure :: evaluate => eval_model
    end type model

    interface cylinder
        module procedure cylinder_init
    end interface cylinder

    interface model
        module procedure model_init
    end interface model

    abstract interface
        real function op(d1, d2)
            implicit none
            real, intent(IN) :: d1, d2
        end function op
    end interface

    private
    public :: sdf, model, cylinder, sphere, container
    public :: union, intersection, subtraction

    contains

        function model_init(array, func) result(out)
            
            implicit none

            type(model) :: out
            procedure(op), pointer :: func
            type(container), intent(IN) :: array(:)
            
            out%array = array
            out%func => func

        end function model_init

        real function eval_model(this, pos)

            implicit none

            class(model) :: this
            type(vector), intent(IN) :: pos

            integer :: i

            eval_model = this%array(1)%p%evaluate(pos)

            do i = 2, size(this%array)
                eval_model = this%func(eval_model, this%array(i)%p%evaluate(pos))
            end do

        end function eval_model

        real function evaluate_fn(this, pos)
            
            implicit none
            class(sdf) :: this
            type(vector), intent(IN) :: pos
        end function evaluate_fn

        function cylinder_init(height, radius, c) result(out)
        
            implicit none
        
            type(cylinder) :: out
            
            real, intent(IN) :: height, radius
            type(vector), optional, intent(IN) :: c
            
            type(vector) :: centre

            if(present(c))then
                centre = c
            else
                centre = vector(0., 0., 0.)
            end if

            out%height = height
            out%radius = radius
            out%centre = centre

        end function cylinder_init

        real function eval_cylinder(this, pos)

            implicit none

            class(cylinder) :: this
            type(vector), intent(IN) :: pos

            eval_cylinder = cylinder_fn(translate_fn(pos, this%centre), this%radius, this%height)

        end function eval_cylinder


        real function eval_sphere(this, pos)

            implicit none

            class(sphere) :: this
            type(vector), intent(IN) :: pos

            eval_sphere = sphere_fn(pos, this%radius)

        end function eval_sphere


        real function eval_box(this, pos)

            implicit none

            class(box) :: this
            type(vector), intent(IN) :: pos

            eval_box = box_fn(pos, this%length)

        end function eval_box


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
            real,         intent(IN) :: c

            sphere_fn = sqrt(p%x**2+p%y**2+p%z**2) - c

        end function sphere_fn


        real function box_fn(p, b)

            implicit none

            type(vector), intent(IN) :: p
            real,         intent(IN) :: b

            type(vector) :: q

            q = abs(p) - b
            box_fn = length(max(q, 0.)) + min(max(q%x, max(q%y, q%z)), 0.)

        end function box_fn

        type(vector) function translate_fn(position, offset)

            implicit none

            type(vector) :: position, offset

            translate_fn = position - offset

        end function translate_fn


        real function union(d1, d2)

            implicit none

            real, intent(IN) :: d1, d2

            union = min(d1, d2)
        end function union

        real function subtraction(d1, d2)

            implicit none

            real, intent(IN) :: d1, d2

            subtraction = max(-d1, d2)

        end function subtraction

        real function intersection(d1, d2)

            implicit none

            real, intent(IN) :: d1, d2

            intersection = max(d1, d2)

        end function intersection

end module sdfs