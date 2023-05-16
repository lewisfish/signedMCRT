module sdfHelpers

    use vector_class
    use constants, only : wp

    implicit none
    
contains

    function rotate_x(angle) result(r)
    ! rotation funcs from https://inspirnathan.com/posts/54-shadertoy-tutorial-part-8/
        use utils, only : deg2rad
        
        real(kind=wp), intent(IN) :: angle
        
        real(kind=wp) :: r(4, 4), c, s, a

        a = deg2rad(angle)
        c = cos(a)
        s = sin(a)

        r(:, 1) = [1._wp, 0._wp, 0._wp, 0._wp]
        r(:, 2) = [0._wp, c,  s,  0._wp]
        r(:, 3) = [0._wp,-s,  c,  0._wp]
        r(:, 4) = [0._wp, 0._wp, 0._wp, 1._wp]

    end function rotate_x

    function rotate_y(angle) result(r)
        
        use utils, only : deg2rad
        
        real(kind=wp), intent(IN) :: angle

        real(kind=wp) :: r(4, 4), c, s, a

        a = deg2rad(angle)
        c = cos(a)
        s = sin(a)

        r(:, 1) = [c,  0._wp, -s,  0._wp]
        r(:, 2) = [0._wp, 1._wp,  0._wp, 0._wp]
        r(:, 3) = [s,  0._wp,  c,  0._wp]
        r(:, 4) = [0._wp, 0._wp,  0._wp, 1._wp]

    end function rotate_y

    function rotate_z(angle) result(r)
        
        use utils, only : deg2rad
        
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
    ! http://www.neilmendoza.com/glsl-rotation-about-an-arbitrary-axis/

        type(vector),  intent(in) :: axis
        real(kind=wp), intent(in) :: angle

        type(vector) :: axist

        real(kind=wp) :: rotmat(4, 4), s, c, oc

        axist = axis%magnitude()
        s = sqrt(1. - angle**2)
        c = angle
        oc = 1._wp - c

    rotmat(:, 1) = [oc * axist%x * axist%x + c, oc * axist%x * axist%y - axist%z * s,&
                    oc * axist%z * axist%x + axist%y * s, 0.0_wp]
    rotmat(:, 2) = [oc * axist%x * axist%y + axist%z * s, oc * axist%y * axist%y + c,&
                    oc * axist%y * axist%z - axist%x * s, 0.0_wp]
    rotmat(:, 3) = [oc * axist%z * axist%x - axist%y * s, oc * axist%y * axist%z + axist%x * s,&
                    oc * axist%z * axist%z + c, 0.0_wp]
    rotmat(:, 4) = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp]

    end function rotmat

    ! https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    ! https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    function rotationAlign(a, b) result(res)
        
        type(vector), intent(in) :: a, b
        
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
                    
        real(kind=wp) :: r(4, 4)

        r(:, 1) = [1._wp, 0._wp, 0._wp, 0._wp]
        r(:, 2) = [0._wp, 1._wp, 0._wp, 0._wp]
        r(:, 3) = [0._wp, 0._wp, 1._wp, 0._wp]
        r(:, 4) = [0._wp, 0._wp, 0._wp, 1._wp]

    end function identity


    function skewSymm(a) result(out)

        type(vector), intent(in) :: a
        real(kind=wp) :: out(4,4)

        out(:, 1) = [0._wp, -a%z, a%y, 0._wp]
        out(:, 2) = [a%z, 0._wp, -a%x, 0._wp]
        out(:, 3) = [-a%y, a%x, 0._wp, 0._wp]
        out(:, 4) = [0._wp, 0._wp, 0._wp, 0._wp]

    end function skewSymm

    function translate(o) result(out)

        type(vector), intent(IN) :: o

        real(kind=wp) :: out(4, 4)

        out(:, 1) = [1._wp, 0._wp, 0._wp, o%x] 
        out(:, 2) = [0._wp, 1._wp, 0._wp, o%y] 
        out(:, 3) = [0._wp, 0._wp, 1._wp, o%z] 
        out(:, 4) = [0._wp, 0._wp, 0._wp, 1._wp] 

    end function translate
end module sdfHelpers