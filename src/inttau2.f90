module inttau2

    implicit none
    
    private
    public :: tauint2

    contains   

    subroutine tauint2(packet, grid, sdfs_array)
    !optical depth integration subroutine
    !
    !
        use random, only : ran2
        use photonMod
        use sdfs
        use gridMod

        use vector_class
        use stackMod
   
        implicit none

        type(cart_grid), intent(IN)    :: grid
        type(photon),    intent(INOUT) :: packet
        type(container), intent(IN) :: sdfs_array(:)

        real         :: tau, d_sdf, t_sdf, taurun, ds(size(sdfs_array)), dstmp(size(sdfs_array))
        real         :: eps, dtot, signs(size(sdfs_array)), n1, n2
        integer      :: i, cur_layer, oldlayer
        type(vector) :: pos, dir, oldpos, N
        logical :: rflag

        pos = packet%pos
        oldpos = pos
        dir = vector(packet%nxp, packet%nyp, packet%nzp)

        !setup sdf distance and current layer
        ds = 0.
        do i = 1, size(ds)
            ds(i) = abs(sdfs_array(i)%p%evaluate(pos))
        end do
        ! packet%cnts = packet%cnts + size(ds)
        d_sdf = minval(ds)

        eps = 1d-8
        tau = -log(ran2())
        taurun = 0.

        cur_layer = packet%layer
        dtot = 0.
        do
            do while(d_sdf > eps)

                t_sdf = d_sdf * sdfs_array(packet%layer)%p%kappa
                if(taurun + t_sdf <= tau)then
                    !move full distance to sdf
                    taurun = taurun + t_sdf
                    oldpos = pos
                    call update_jmean(oldpos, grid, dir, d_sdf, packet, sdfs_array(packet%layer)%p%mua)
                    pos = pos + d_sdf * dir
                    dtot = dtot + d_sdf
                else
                    !run out of tau so move remaining tau and exit
                    d_sdf = (tau - taurun) / sdfs_array(packet%layer)%p%kappa
                    dtot = dtot + d_sdf
                    taurun = tau
                    oldpos = pos
                    pos = pos + d_sdf * dir
                    call update_jmean(oldpos, grid, dir, d_sdf, packet, sdfs_array(packet%layer)%p%mua)
                    exit
                end if
                ! get distance to nearest sdf
                ds = 0.
                do i = 1, size(ds)
                    ds(i) = sdfs_array(i)%p%evaluate(pos)
                end do
                packet%cnts = packet%cnts + size(ds)

                d_sdf = minval(abs(ds),dim=1)
                !check if outside all sdfs
                if(minval(ds) >= 0.)then
                    packet%tflag = .true.
                    exit
                end if
            end do

            if(taurun >= tau .or. packet%tflag)then
                exit
            end if

            ds = 0.
            do i = 1, size(ds)
                ds(i) = sdfs_array(i)%p%evaluate(pos)
            end do
            packet%cnts = packet%cnts + size(ds)

            dstmp = ds
            ds = abs(ds)

            !step a bit into next sdf to get n2
            d_sdf = minval(ds) + 2.*eps
            pos = pos + d_sdf*dir
            ds = 0.
            do i = 1, size(ds)
                ds(i) = sdfs_array(i)%p%evaluate(pos)
            end do
            packet%cnts = packet%cnts + size(ds)

            !step back inside original sdf
            pos = pos - d_sdf*dir

            ! check going in to media i.e +ive to -ive
            signs = 0.
            do i = 1, size(ds)
                if(dstmp(i) >= 0. .and. ds(i) < 0.)then
                    signs(i) = -1.
                else
                    signs(i) = 1.
                end if
            end do

            !if all signs +ive then do additional check t0 get correct layer
            if(sum(signs) == real(size(ds)))then
                do i = 1, size(ds)
                    if(dstmp(i) <= 0. .and. ds(i) > 0. .and. i /= cur_layer .or. ds(i) < 0.)then
                        signs(i) = -1.
                    else
                        signs(i) = 1.
                    end if
                end do
            end if

            !check for fresnel reflection
            n1 = sdfs_array(cur_layer)%p%n
            oldlayer = cur_layer
            cur_layer = minloc(signs,dim=1)
            n2 = sdfs_array(cur_layer)%p%n

            if (n1 /= n2)then
                N = calcNormal(pos, sdfs_array(cur_layer)%p)
                N = N%magnitude()
                rflag = .false.
                call reflect_refract(dir, N, n1, n2, rflag)
                tau = -log(ran2())
                taurun = 0.
                if(rflag)then
                    cur_layer = oldlayer
                end if
            end if
            packet%layer = cur_layer

        end do

        packet%pos = pos
        packet%nxp = dir%x
        packet%nyp = dir%y
        packet%nzp = dir%z

        packet%phi = atan2(dir%y, dir%x)
        packet%sinp = sin(packet%phi)
        packet%cosp = cos(packet%phi)

        packet%cost = dir%z
        packet%sint = sqrt(1.-packet%cost**2)

        if(abs(packet%pos%x) > grid%xmax)then
            packet%tflag = .true.
        end if
        if(abs(packet%pos%y) > grid%ymax)then
            packet%tflag = .true.
        end if
        if(abs(packet%pos%z) > grid%zmax)then
            packet%tflag = .true.
        end if
    end subroutine tauint2
   

    subroutine update_jmean(pos, grid, dir, d_sdf, packet, mua)
        
        use vector_class
        use iarray, only: jmean
        use gridMod
        use photonMod

        implicit none
        
        type(cart_grid), intent(IN)    :: grid
        type(vector),    intent(IN)    :: dir
        real,            intent(IN)    :: d_sdf, mua
        type(vector),    intent(INOUT) :: pos
        type(photon),    intent(INOUT) :: packet

        type(vector) :: old_pos
        logical      :: ldir(3)
        integer      :: celli, cellj, cellk
        real         :: dcell, delta=1d-8, d

        old_pos = vector(pos%x+grid%xmax, pos%y+grid%ymax, pos%z+grid%zmax)
        call update_voxels(old_pos, grid, celli, cellj, cellk)
        packet%xcell = celli
        packet%ycell = cellj
        packet%zcell = cellk

        d = 0.
        if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
            packet%tflag = .true.
            pos = vector(old_pos%x-grid%xmax, old_pos%y-grid%ymax, old_pos%z-grid%zmax)
            return
        end if
        do
            ldir = (/.FALSE., .FALSE., .FALSE./)

            dcell = wall_dist(grid, celli, cellj, cellk, old_pos, dir, ldir)
            if(d + dcell > d_sdf)then
                dcell = d_sdf - d
                d = d_sdf
!$omp atomic
                jmean(celli, cellj, cellk) = jmean(celli, cellj, cellk) + dcell!*mua!real(packet%cnts)
                ! packet%cnts = 0
                call update_pos(old_pos, grid, celli, cellj, cellk, dcell, .false., dir, ldir, delta)
                exit
            else
                d = d + dcell
!$omp atomic
                jmean(celli, cellj, cellk) = jmean(celli, cellj, cellk) + dcell!*mua!real(packet%cnts)
                ! packet%cnts = 0
                call update_pos(old_pos, grid, celli, cellj, cellk, dcell, .true., dir, ldir, delta)
            end if
            if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
                packet%tflag = .true.
                exit
            end if
        end do
        pos = vector(old_pos%x-grid%xmax, old_pos%y-grid%ymax, old_pos%z-grid%zmax)
        packet%xcell = celli
        packet%ycell = cellj
        packet%zcell = cellk

    end subroutine update_jmean

    real function wall_dist(grid, celli, cellj, cellk, pos, dir, ldir)
    !funtion that returns distant to nearest wall and which wall that is (x,y or z)
    !
    !
        use vector_class
        use gridMod

        implicit none

        type(cart_grid), intent(IN) :: grid
        type(vector),    intent(IN) :: pos, dir
        logical, intent(INOUT) :: ldir(:)
        integer, intent(INOUT) :: celli, cellj, cellk
        real                   :: dx, dy, dz


        dx = -999.
        dy = -999.
        dz = -999.

        if(dir%x > 0.)then
            dx = (grid%xface(celli+1) - pos%x)/dir%x
        elseif(dir%x < 0.)then
            dx = (grid%xface(celli) - pos%x)/dir%x
        elseif(dir%x == 0.)then
            dx = 100000.
        end if

        if(dir%y > 0.)then
            dy = (grid%yface(cellj+1) - pos%y)/dir%y
        elseif(dir%y < 0.)then
            dy = (grid%yface(cellj) - pos%y)/dir%y
        elseif(dir%y == 0.)then
            dy = 100000.
        end if

        if(dir%z > 0.)then
            dz = (grid%zface(cellk+1) - pos%z)/dir%z
        elseif(dir%z < 0.)then
            dz = (grid%zface(cellk) - pos%z)/dir%z
        elseif(dir%z == 0.)then
            dz = 100000.
        end if

        wall_dist = min(dx, dy, dz)
        if(wall_dist < 0.)then
            print*,'dcell < 0.0 warning! ',wall_dist
            print*,dx,dy,dz
            print*,dir
            print*,celli,cellj,cellk
            error stop 1
        end if

        if(wall_dist == dx)ldir = [.TRUE., .FALSE., .FALSE.]
        if(wall_dist == dy)ldir = [.FALSE., .TRUE., .FALSE.]
        if(wall_dist == dz)ldir = [.FALSE., .FALSE., .TRUE.]
        if(.not.ldir(1) .and. .not.ldir(2) .and. .not.ldir(3))print*,'Error in dir flag'
      
   end function wall_dist


    subroutine update_pos(pos, grid, celli, cellj, cellk, dcell, wall_flag, dir, ldir, delta)
    !routine that upates postions of photon and calls fresnel routines if photon leaves current voxel
    !
    !
        ! use iarray,      only : xface, yface, zface
        use utils,       only : str, red, bold, colour
        use vector_class
        use gridMod

        implicit none
      
        type(cart_grid), intent(IN)    :: grid
        type(vector),    intent(IN)    :: dir
        logical,         intent(IN)    :: wall_flag, ldir(:)
        real,            intent(IN)    :: dcell, delta
        type(vector),    intent(INOUT) :: pos
        integer,         intent(INOUT) :: celli, cellj, cellk

        character(len=32)      :: tmp  

        if(wall_flag)then

            if(ldir(1))then
                if(dir%x > 0.)then
                    pos%x = grid%xface(celli+1) + delta
                elseif(dir%x < 0.)then
                    pos%x = grid%xface(celli) - delta
                else
                    print*,'Error in x ldir in update_pos', ldir, dir
                end if
                pos%y = pos%y + dir%y*dcell 
                pos%z = pos%z + dir%z*dcell
            elseif(ldir(2))then
                if(dir%y > 0.)then
                    pos%y = grid%yface(cellj+1) + delta
                elseif(dir%y < 0.)then
                    pos%y = grid%yface(cellj) - delta
                else
                    print*,'Error in y ldir in update_pos', ldir, dir
                end if
                pos%x = pos%x + dir%x*dcell
                pos%z = pos%z + dir%z*dcell
            elseif(ldir(3))then
                if(dir%z > 0.)then
                    pos%z = grid%zface(cellk+1) + delta
                elseif(dir%z < 0.)then
                    pos%z = grid%zface(cellk) - delta
                else
                    print*,'Error in z ldir in update_pos', ldir, dir
                end if
                pos%x = pos%x + dir%x*dcell
                pos%y = pos%y + dir%y*dcell 
            else
                tmp = colour('Error in update_pos... '//str(ldir), red, bold)
                error stop 1
            end if
        else
            pos%x = pos%x + dir%x*dcell
            pos%y = pos%y + dir%y*dcell 
            pos%z = pos%z + dir%z*dcell
        end if

        if(wall_flag)then
            call update_voxels(pos, grid, celli, cellj, cellk)
        end if

    end subroutine update_pos


    subroutine update_voxels(pos, grid, celli, cellj, cellk)
    !updates the current voxel based upon position
    !
    !
        use vector_class
        use gridMod

        implicit none

        type(cart_grid), intent(IN) :: grid
        type(vector),    intent(IN)    :: pos
        integer, intent(INOUT) :: celli, cellj, cellk

        celli = find(pos%x, grid%xface) 
        cellj = find(pos%y, grid%yface)
        cellk = find(pos%z, grid%zface) 
        ! celli = floor(grid%nxg * (pos%x) / (2. * grid%xmax)) + 1
        ! cellj = floor(grid%nyg * (pos%y) / (2. * grid%ymax)) + 1
        ! cellk = floor(grid%nzg * (pos%z) / (2. * grid%zmax)) + 1

        if(celli > grid%nxg .or. celli < 1)celli = -1
        if(cellj > grid%nyg .or. cellj < 1)cellj = -1
        if(cellk > grid%nzg .or. cellk < 1)cellk = -1

    end subroutine update_voxels

    integer function find(val, a)
    !searchs for bracketing indicies for a value val in an array a
    !
    !
        implicit none

        real, intent(IN) :: val, a(:)
        integer          :: n, lo, mid, hi

        n = size(a)
        lo = 0
        hi = n + 1

        if (val == a(1)) then
            find = 1
        else if (val == a(n)) then
            find = n-1
        else if((val > a(n)) .or. (val < a(1))) then
            find = -1
        else
            do
                if (hi-lo <= 1) exit
                mid = (hi+lo)/2
                if (val >= a(mid)) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            find = lo
        end if
    end function find


    subroutine reflect_refract(I, N, n1, n2, rflag)

        use random, only : ran2
        use vector_class

        implicit none

        type(vector), intent(INOUT) :: I
        type(vector), intent(INOUT) :: N
        real,         intent(IN)    :: n1, n2
        logical,      intent(OUT)   :: rflag

        rflag = .FALSE.

        if(ran2() <= fresnel(I, N, n1, n2))then
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
        use vector_class

        implicit none

        type(vector), intent(INOUT) :: I
        type(vector), intent(IN)    :: N

        type(vector) :: R

        R = I - 2. * (N .dot. I) * N
        I = R

    end subroutine reflect


    subroutine refract(I, N, eta)
    !   get vector of refracted photon
    !
    !
        use vector_class

        implicit none

        type(vector), intent(INOUT) :: I
        type(vector), intent(IN)    :: N
        real,         intent(IN)    :: eta

        type(vector) :: T, Ntmp

        real :: c1, c2

        Ntmp = N

        c1 = (Ntmp .dot. I)
        if(c1 < 0.)then
            c1 = -c1
        else
            Ntmp = (-1.) * N
        end if
        c2 = sqrt(1. - (eta)**2 * (1.-c1**2))

        T = eta*I + (eta * c1 - c2) * Ntmp 

        I = T

    end subroutine refract


    function fresnel(I, N, n1, n2) result (tir)
    !   calculates the fresnel coefficents
    !
    !
        use vector_class
        use ieee_arithmetic, only : ieee_is_nan

        implicit none

        real, intent(IN)         :: n1, n2
        type(vector), intent(IN) :: I, N

        real             ::  costt, sintt, sint2, cost2, tir, f1, f2

        costt = abs(I .dot. N)

        sintt = sqrt(1. - costt * costt)
        sint2 = n1/n2 * sintt
        if(sint2 > 1.)then
            tir = 1.0
            return
        elseif(costt == 1.)then
            tir = 0.
            return
        else
            sint2 = (n1/n2)*sintt
            cost2 = sqrt(1. - sint2 * sint2)
            f1 = abs((n1*costt - n2*cost2) / (n1*costt + n2*cost2))**2
            f2 = abs((n1*cost2 - n2*costt) / (n1*cost2 + n2*costt))**2

            tir = 0.5 * (f1 + f2)
        if(ieee_is_nan(tir) .or. tir > 1. .or. tir < 0.)print*,'TIR: ', tir, f1, f2, costt,sintt,cost2,sint2
            return
        end if
    end function fresnel
end module inttau2