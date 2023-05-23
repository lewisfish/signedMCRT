!   Module provides signed distance functions (SDFs) for various shapes 
!   and some operations to modify them
!   All SDF functions are adapted from Inigo Quilex exhaustive list at:
!   https://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm!
!   API based upon https://fortran-lang.discourse.group/t/attempting-type-erasure-in-fortran/4402/2
module sdf_baseMod

    use constants,         only : wp
    use opticalProperties, only : opticalProp_t
    use sdfHelpers,        only : identity
    use vector_class

    implicit none

    type, abstract :: sdf_base
        type(opticalProp_t) :: optProps
        real(kind=wp) :: transform(4, 4)
        integer :: layer
        contains
            procedure(evalInterface), deferred :: evaluate
    end type sdf_base
        
    type, extends(sdf_base) :: sdf
        class(sdf_base), allocatable :: value
        contains
            procedure :: getKappa
            procedure :: getAlbedo
            procedure :: getMua, gethgg, getG2, getN
            procedure :: evaluate => sdf_evaluate
            procedure, private :: sdf_assign
            generic :: assignment(=) => sdf_assign
    end type sdf

    type, extends(sdf_base) :: model
        type(sdf), allocatable         :: array(:)
        procedure(op), nopass, pointer :: func
        real(kind=wp) :: k
        contains
            procedure :: evaluate => eval_model
    end type model
!####################################################################
    abstract interface
        pure elemental function evalInterface(this, pos) result(res)
            use vector_class
            use constants, only : wp
            import sdf_base
            class(sdf_base), intent(in) :: this
            type(vector),    intent(in) :: pos
            real(kind=wp) :: res
        end function

        pure function primitive(pos) result(res)
            use vector_class, only : vector
            use constants,    only : wp
            implicit none
            type(vector), intent(IN) :: pos
            real(kind=wp) :: res
        end function primitive

        pure function op(d1, d2, k) result(res)
            use constants, only : wp
            implicit none
            real(kind=wp), intent(IN) :: d1, d2, k
            real(kind=wp) :: res
        end function op
    end interface

    interface sdf
        module procedure sdf_new
    end interface

    interface model
        module procedure model_init
    end interface

    interface render
        module procedure render_sub, render_vec
    end interface

    private
    public :: model, sdf, sdf_base, primitive, op, calcNormal, render

    contains

    function model_init(array, func, kopt) result(out)

        type(model) :: out

        procedure(op) :: func
        type(sdf),               intent(IN) :: array(:)
        real(kind=wp), optional, intent(IN) :: kopt
        integer :: i

        out%array = array

        out%func => func
        if(present(kopt))then
            out%k = kopt
        else
            out%k = 0._wp
        end if

        do i = 2, size(array)
            if(array(1)%value%optProps%value%mus /= array(i)%value%optProps%value%mus)then
                print*,"Error mismatch in model mus in object: ",i
            end if
            if(array(1)%value%optProps%value%mua /= array(i)%value%optProps%value%mua)then
                print*,"Error mismatch in model mua in object: ",i
            end if
            if(array(1)%value%optProps%value%hgg /= array(i)%value%optProps%value%hgg)then
                print*,"Error mismatch in model hgg in object: ",i
            end if
            if(array(1)%value%optProps%value%n /= array(i)%value%optProps%value%n)then
                print*,"Error mismatch in model n in object: ",i
            end if
            if(array(1)%value%layer /= array(i)%value%layer)then
                print*,"Error mismatch in model layer in object: ",i
            end if
        end do

            out%optProps = array(1)%value%optProps
            out%layer = array(1)%value%layer

    end function model_init

    pure elemental function eval_model(this, pos) result(res)

        class(model), intent(in) :: this
        type(vector), intent(in) :: pos
        real(kind=wp) :: res

        integer :: i

        res = this%array(1)%value%evaluate(pos)
        do i = 2, size(this%array)
            res = this%func(res, this%array(i)%value%evaluate(pos), this%k)
        end do

    end function eval_model

!#############################################################
!           Helpers
!#############################################################
    type(vector) function calcNormal(p, obj)

        type(vector), intent(IN) :: p
        class(sdf_base) :: obj

        real(kind=wp) :: h
        type(vector) :: xyy, yyx, yxy, xxx

        h = 1e-6_wp
        xyy = vector(1._wp, -1._wp, -1._wp)
        yyx = vector(-1._wp, -1._wp, 1._wp)
        yxy = vector(-1._wp, 1._wp, -1._wp)
        xxx = vector(1._wp, 1._wp, 1._wp)

        calcNormal = xyy*obj%evaluate(p + xyy*h) +  &
                    yyx*obj%evaluate(p + yyx*h) +  &
                    yxy*obj%evaluate(p + yxy*h) +  &
                    xxx*obj%evaluate(p + xxx*h)

        calcNormal = calcNormal%magnitude()

    end function calcNormal

    function getKappa(this) result(res)

        class(sdf) :: this
        real(kind=wp) :: res

        res = this%value%optProps%value%kappa

    end function getKappa

    function getMua(this) result(res)

        class(sdf) :: this
        real(kind=wp) :: res

        res = this%value%optProps%value%mua

    end function getMua

    function gethgg(this) result(res)

        class(sdf) :: this
        real(kind=wp) :: res

        res = this%value%optProps%value%hgg

    end function gethgg

    function getg2(this) result(res)

        class(sdf) :: this
        real(kind=wp) :: res

        res = this%value%optProps%value%g2

    end function getg2

    function getN(this) result(res)

        class(sdf) :: this
        real(kind=wp) :: res

        res = this%value%optProps%value%n

    end function getN

    function getAlbedo(this) result(res)

        class(sdf) :: this
        real(kind=wp) :: res

        res = this%value%optProps%value%albedo

    end function getAlbedo
!#########################################################################
!       SDF bound procedures
!#########################################################################
    pure elemental function sdf_evaluate(this, pos) result(res)

        class(sdf), intent(in) :: this
        type(vector), intent(in) :: pos
        real(kind=wp) :: res

        res = this%value%evaluate(pos)

    end function sdf_evaluate
 
    ! sdf initializer
    subroutine sdf_assign(lhs, rhs)

        class(sdf),      intent(inout) :: lhs
        class(sdf_base), intent(in)    :: rhs

        if (allocated(lhs%value))deallocate(lhs%value)
        ! Prevent nested derived type
        select type (rhsT=>rhs)
            class is (sdf)
                if(allocated(rhsT%value))allocate(lhs%value,source=rhsT%value)
            class default
                allocate(lhs%value,source=rhsT)
        end select
    end subroutine sdf_assign
 
    ! sdf initializer
    type(sdf) function sdf_new(rhs) result(lhs)

        class(sdf_base), intent(in) :: rhs
        allocate(lhs%value,source=rhs)

    end function sdf_new


    subroutine render_vec(cnt, state)

        use sim_state_mod, only : settings_t

        type(settings_t), intent(IN) :: state
        type(sdf),  intent(IN) :: cnt(:)
        
        type(vector):: extent

        extent = vector(state%grid%xmax, state%grid%ymax, state%grid%zmax)
    
        call render_sub(cnt, extent, state%render_size, state)

    end subroutine render_vec

    subroutine render_sub(cnt, extent, samples, state)

        use sim_state_mod, only : settings_t
        use utils,         only : pbar
        use constants,     only : fileplace, sp
        use writer_mod
                  
        type(settings_t), intent(IN) :: state
        type(sdf),  intent(IN) :: cnt(:)
        integer,          intent(IN) :: samples(3)
        type(vector),     intent(IN) :: extent

        type(vector)               :: pos, wid
        integer                    :: i, j, k, u, id
        real(kind=wp)              :: x, y, z, ds(size(cnt)), ns(3), minvalue
        real(kind=sp), allocatable :: image(:, :, :)
        type(pbar)                 :: bar

        ns = nint(samples / 2._wp)
        allocate(image(samples(1), samples(2), samples(3)))
        wid = vector(extent%x/ns(1), extent%y/ns(2), extent%z/ns(3))
        bar = pbar(samples(1))
!$omp parallel default(none) shared(cnt, ns, wid, image, samples, bar)&
!$omp private(i, x, y, z, pos, j, k, u, ds, id, minvalue)
!$omp do
        do i = 1, samples(1)
            x = (i-ns(1)) *wid%x
            do j = 1, samples(2)
                y = (j-ns(2)) *wid%y
                do k = 1, samples(3)
                    z = (k-ns(3)) * wid%z
                    pos = vector(x, y, z)
                    ds = 0._wp
                    do u = 1, size(ds)
                        ds(u) = cnt(u)%evaluate(pos)
                    end do
                    if(all(ds > 0._wp))then
                        id=0
                    else
                        if(all(ds < 0._wp))then
                            id = cnt(maxloc(ds,dim=1))%layer
                        else
                            ds = -1.*ds
                            minvalue = 1000000._wp
                            do u = 1, size(ds)
                                if(ds(u) > 0. .and. ds(u) < minvalue)then
                                    minvalue = ds(u)
                                    id = u
                                end if
                            end do
                            ! id = cnt(minloc(ds),dim=1))%p%layer
                        end if
                    end if
                    ! if(id == 0._wp)then
                    !     image(i,j,k)=-99._wp
                    ! else
                        image(i, j, k) = real(id)!cnt(id)%p%mua
                    ! end if
                end do
            end do
            call bar%progress()
        end do
!$OMP end  do
!$OMP end parallel
        call write_data(image, trim(fileplace)//state%renderfile, state, overwrite=.true.)
    end subroutine render_sub
end module sdf_baseMod