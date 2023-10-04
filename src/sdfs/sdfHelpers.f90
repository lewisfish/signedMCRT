module sdfHelpers
    !! Collection of helper functions for SDFs:

    !! This module defines transforms that can be applied to each SDF:

    !! - Rotate_{x,y,z}
    !! - Translate
    !! - RotationAlign (not tested)
    !! - RotMat (not tested)
    !! - Identity
    !! - SkewSymm
    
    use vector_class
    use constants, only : wp

    implicit none
    
    private
    public :: rotate_x, rotate_y, rotate_z, rotmat, rotationAlign, identity, skewSymm, translate

contains

    function rotate_x(angle) result(r)
       !! rotation in the x-axis function from [here](https://inspirnathan.com/posts/54-shadertoy-tutorial-part-8/)
        use utils, only : deg2rad
        
        !> Angle to rotate by
        real(kind=wp), intent(IN) :: angle
        
        real(kind=wp) :: r(4, 4), c, s, a

        a = deg2rad(angle)
        c = cos(a)
        s = sin(a)

        r(:, 1) = [1._wp, 0._wp, 0._wp, 0._wp]
        r(:, 2) = [0._wp, c,  -s,  0._wp]
        r(:, 3) = [0._wp, s,  c,  0._wp]
        r(:, 4) = [0._wp, 0._wp, 0._wp, 1._wp]

    end function rotate_x

    function rotate_y(angle) result(r)
        !! rotation in the y-axis function from [here](https://inspirnathan.com/posts/54-shadertoy-tutorial-part-8/)

        use utils, only : deg2rad
        
        !> Angle to rotate by
        real(kind=wp), intent(IN) :: angle

        real(kind=wp) :: r(4, 4), c, s, a

        a = deg2rad(angle)
        c = cos(a)
        s = sin(a)

        r(:, 1) = [c,  0._wp, s,  0._wp]
        r(:, 2) = [0._wp, 1._wp,  0._wp, 0._wp]
        r(:, 3) = [-s,  0._wp,  c,  0._wp]
        r(:, 4) = [0._wp, 0._wp,  0._wp, 1._wp]

    end function rotate_y

    function rotate_z(angle) result(r)
        !! rotation in the z-axis function from [here](https://inspirnathan.com/posts/54-shadertoy-tutorial-part-8/)

        use utils, only : deg2rad
        
        !> Angle to rotate by
        real(kind=wp), intent(IN) :: angle

        real(kind=wp) :: r(4, 4), c, s, a

        a = deg2rad(angle)
        c = cos(a)
        s = sin(a)

        r(:, 1) = [c, -s,  0._wp, 0._wp]
        r(:, 2) = [s,  c,  0._wp, 0._wp]
        r(:, 3) = [0._wp, 0._wp, 1._wp, 0._wp]
        r(:, 4) = [0._wp, 0._wp, 0._wp, 1._wp]

    end function rotate_z

    function rotmat(axis, angle)
    !! Rotate around around an axis by a given angle taken from [here](http://www.neilmendoza.com/glsl-rotation-about-an-arbitrary-axis/)
        use utils, only : deg2rad
        !> Axis to rotate around
        type(vector),  intent(in) :: axis
        !> Angle to rotate by in degrees
        real(kind=wp), intent(in) :: angle

        type(vector) :: axist

        real(kind=wp) :: rotmat(4, 4), s, c, oc, a

        axist = axis%magnitude()
        a = deg2rad(angle)

        s = sin(a)
        c = cos(a)
        oc = 1._wp - c

    rotmat(:, 1) = [oc * axist%x * axist%x + c, oc * axist%x * axist%y - axist%z * s,&
                    oc * axist%z * axist%x + axist%y * s, 0.0_wp]
    rotmat(:, 2) = [oc * axist%x * axist%y + axist%z * s, oc * axist%y * axist%y + c,&
                    oc * axist%y * axist%z - axist%x * s, 0.0_wp]
    rotmat(:, 3) = [oc * axist%z * axist%x - axist%y * s, oc * axist%y * axist%z + axist%x * s,&
                    oc * axist%z * axist%z + c, 0.0_wp]
    rotmat(:, 4) = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp]

    end function rotmat

    function rotationAlign(a, b) result(res)
        !! Calculate the rotation matrix to rotate vector a onto b
        !! [ref1](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula)
        !! [ref2](https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d)
        
        !> Vector to rotate. Unit vector
        type(vector), intent(in) :: a
        !> Vector to be rotated onto. Unit vector
        type(vector), intent(in) :: b
        
        type(vector)  :: v
        real(kind=wp) :: c, k, res(4, 4), v_x(4, 4), v_x2(4, 4)

        v = a .cross. b
        c = a .dot. b
        k = 1._wp / (1._wp + c)
        
        !skew-symmetric matrix
        v_x(:, 1) = [0._wp     , -1._wp*v%z, v%y       , 0._wp]
        v_x(:, 2) = [v%z       , 0._wp     , -1._wp*v%x, 0._wp]
        v_x(:, 3) = [-1._wp*v%y, v%x       , 0._wp     , 0._wp]
        v_x(:, 4) = [0._wp, 0._wp, 0._wp, 0._wp]
        
        v_x2 = matmul(v_x, v_x)
        res = identity() + v_x + v_x2*k

    end function rotationAlign

    function identity() result(r)
        !! Returns the identity transformation matrix

        real(kind=wp) :: r(4, 4)

        r(:, 1) = [1._wp, 0._wp, 0._wp, 0._wp]
        r(:, 2) = [0._wp, 1._wp, 0._wp, 0._wp]
        r(:, 3) = [0._wp, 0._wp, 1._wp, 0._wp]
        r(:, 4) = [0._wp, 0._wp, 0._wp, 1._wp]

    end function identity


    function skewSymm(a) result(out)
        !! Calculate the Skew Symmetric matrix for a given vector

        !> Vector to calculate the skew symmetric matrix for.
        type(vector), intent(in) :: a
        real(kind=wp) :: out(4,4)

        out(:, 1) = [0._wp, -a%z, a%y, 0._wp]
        out(:, 2) = [a%z, 0._wp, -a%x, 0._wp]
        out(:, 3) = [-a%y, a%x, 0._wp, 0._wp]
        out(:, 4) = [0._wp, 0._wp, 0._wp, 0._wp]

    end function skewSymm

    function translate(o) result(out)
        !! Returns the Translation matrix for a given vector translation.

        !> Vector to translate by.
        type(vector), intent(IN) :: o

        real(kind=wp) :: out(4, 4)

        out(:, 1) = [1._wp, 0._wp, 0._wp, o%x] 
        out(:, 2) = [0._wp, 1._wp, 0._wp, o%y] 
        out(:, 3) = [0._wp, 0._wp, 1._wp, o%z] 
        out(:, 4) = [0._wp, 0._wp, 0._wp, 1._wp] 

    end function translate
end module sdfHelpers