module inttau2

    use constants, only : wp

    implicit none
    
    private
    public :: tauint2

    contains   

    subroutine tauint2(packet, sdfs_array)
    ! optical depth integration subroutine
    ! Moves photons to interaction location
    ! Calculated is any reflection or refraction happens whilst moving
    !
        use random, only : ran2
        use photonMod
        use sdfs
        use sim_state_mod, only : state

        use vector_class, only : vector
   
        implicit none

        type(photon),    intent(INOUT) :: packet
        type(container), intent(IN)    :: sdfs_array(:)

        real(kind=wp) :: tau, d_sdf, t_sdf, taurun, ds(size(sdfs_array)), dstmp(size(sdfs_array))
        real(kind=wp) :: eps, dtot, signs(size(sdfs_array)), n1, n2
        integer       :: i, cur_layer, oldlayer
        type(vector)  :: pos, dir, oldpos, N
        logical       :: rflag

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

        eps = 1e-8_wp
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
                    call update_jmean(oldpos, dir, d_sdf, packet, sdfs_array(packet%layer)%p%mua)
                    pos = pos + d_sdf * dir
                    dtot = dtot + d_sdf
                else
                    !run out of tau so move remaining tau and exit
                    d_sdf = (tau - taurun) / sdfs_array(packet%layer)%p%kappa
                    dtot = dtot + d_sdf
                    taurun = tau
                    oldpos = pos
                    pos = pos + d_sdf * dir
                    call update_jmean(oldpos, dir, d_sdf, packet, sdfs_array(packet%layer)%p%mua)
                    exit
                end if
                ! get distance to nearest sdf
                ds = 0._wp
                do i = 1, size(ds)
                    ds(i) = sdfs_array(i)%p%evaluate(pos)
                end do
                packet%cnts = packet%cnts + size(ds)

                d_sdf = minval(abs(ds),dim=1)
                !check if outside all sdfs
                if(minval(ds) >= 0._wp)then
                    packet%tflag = .true.
                    exit
                end if
            end do

            if(taurun >= tau .or. packet%tflag)then
                exit
            end if

            ds = 0._wp
            do i = 1, size(ds)
                ds(i) = sdfs_array(i)%p%evaluate(pos)
            end do
            packet%cnts = packet%cnts + size(ds)

            dstmp = ds
            ds = abs(ds)

            !step a bit into next sdf to get n2
            d_sdf = minval(ds) + 2._wp*eps
            pos = pos + d_sdf*dir
            ds = 0._wp
            do i = 1, size(ds)
                ds(i) = sdfs_array(i)%p%evaluate(pos)
            end do
            packet%cnts = packet%cnts + size(ds)

            !step back inside original sdf
            pos = pos - d_sdf*dir

            ! check going in to media i.e +ive to -ive
            signs = 0._wp
            do i = 1, size(ds)
                if(dstmp(i) >= 0._wp .and. ds(i) < 0._wp)then
                    signs(i) = -1._wp
                else
                    signs(i) = 1._wp
                end if
            end do

            !if all signs +ive then do additional check to get correct layer
            if(sum(signs) == real(size(ds), kind=wp))then
                do i = 1, size(ds)
                    if(dstmp(i) <= 0._wp .and. ds(i) > 0._wp .and. i /= cur_layer .or. ds(i) < 0._wp)then
                        signs(i) = -1._wp
                    else
                        signs(i) = 1._wp
                    end if
                end do
            end if

            !check for fresnel reflection
            n1 = sdfs_array(cur_layer)%p%n
            oldlayer = cur_layer
            cur_layer = minloc(signs,dim=1)
            n2 = sdfs_array(cur_layer)%p%n

            !carry out refelction/refraction
            if (n1 /= n2)then
                N = calcNormal(pos, sdfs_array(cur_layer)%p)
                N = N%magnitude()
                rflag = .false.
                call reflect_refract(dir, N, n1, n2, rflag)
                tau = -log(ran2())
                taurun = 0._wp
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
        packet%sint = sqrt(1._wp - packet%cost**2)

        if(abs(packet%pos%x) > state%grid%xmax)then
            packet%tflag = .true.
        end if
        if(abs(packet%pos%y) > state%grid%ymax)then
            packet%tflag = .true.
        end if
        if(abs(packet%pos%z) > state%grid%zmax)then
            packet%tflag = .true.
        end if
    end subroutine tauint2
   

    subroutine update_jmean(pos, dir, d_sdf, packet, mua)
    ! recored fluence using path length estimators.

        use vector_class
        use photonMod

        use iarray,        only: jmean
        use sim_state_mod, only : state
    
        implicit none
        
        type(vector),    intent(IN)    :: dir
        real(kind=wp),   intent(IN)    :: d_sdf, mua
        type(vector),    intent(INOUT) :: pos
        type(photon),    intent(INOUT) :: packet

        type(vector)  :: old_pos
        logical       :: ldir(3)
        integer       :: celli, cellj, cellk
        real(kind=wp) :: dcell, delta=1e-8_wp, d

        !convert to different coordinate system. Origin is at lower left corner of fluence grid
        old_pos = vector(pos%x+state%grid%xmax, pos%y+state%grid%ymax, pos%z+state%grid%zmax)
        call update_voxels(old_pos, celli, cellj, cellk)
        packet%xcell = celli
        packet%ycell = cellj
        packet%zcell = cellk

        d = 0._wp
        !if packet outside grid return
        if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
            packet%tflag = .true.
            pos = vector(old_pos%x-state%grid%xmax, old_pos%y-state%grid%ymax, old_pos%z-state%grid%zmax)
            return
        end if
        !move photon through grid updating path length estimators
        do
            ! print*,old_pos,celli,cellj,cellk,delta,1.d-8
            ! stop
            ldir = (/.FALSE., .FALSE., .FALSE./)

            dcell = wall_dist(celli, cellj, cellk, old_pos, dir, ldir)
            if(d + dcell > d_sdf)then
                dcell = d_sdf - d
                d = d_sdf
!$omp atomic
                jmean(celli, cellj, cellk) = jmean(celli, cellj, cellk) + dcell!*mua!real(packet%cnts)
                ! packet%cnts = 0
                call update_pos(old_pos, celli, cellj, cellk, dcell, .false., dir, ldir, delta)
                exit
            else
                d = d + dcell
!$omp atomic
                jmean(celli, cellj, cellk) = jmean(celli, cellj, cellk) + dcell!*mua!real(packet%cnts)
                ! packet%cnts = 0
                call update_pos(old_pos, celli, cellj, cellk, dcell, .true., dir, ldir, delta)
            end if
            if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
                packet%tflag = .true.
                exit
            end if
        end do
        pos = vector(old_pos%x-state%grid%xmax, old_pos%y-state%grid%ymax, old_pos%z-state%grid%zmax)
        packet%xcell = celli
        packet%ycell = cellj
        packet%zcell = cellk

    end subroutine update_jmean

    function wall_dist(celli, cellj, cellk, pos, dir, ldir) result(res)
    !funtion that returns distant to nearest wall and which wall that is (x, y, or z)
    !
    !
        use vector_class
        use sim_state_mod, only : state

        implicit none

        type(vector),    intent(IN)    :: pos, dir
        logical,         intent(INOUT) :: ldir(:)
        integer,         intent(INOUT) :: celli, cellj, cellk
        real(kind=wp) :: res

        real(kind=wp) :: dx, dy, dz

        dx = -999._wp
        dy = -999._wp
        dz = -999._wp

        if(dir%x > 0._wp)then
            dx = (state%grid%xface(celli+1) - pos%x)/dir%x
        elseif(dir%x < 0._wp)then
            dx = (state%grid%xface(celli) - pos%x)/dir%x
        elseif(dir%x == 0._wp)then
            dx = 100000._wp
        end if

        if(dir%y > 0._wp)then
            dy = (state%grid%yface(cellj+1) - pos%y)/dir%y
        elseif(dir%y < 0._wp)then
            dy = (state%grid%yface(cellj) - pos%y)/dir%y
        elseif(dir%y == 0._wp)then
            dy = 100000._wp
        end if

        if(dir%z > 0._wp)then
            dz = (state%grid%zface(cellk+1) - pos%z)/dir%z
        elseif(dir%z < 0._wp)then
            dz = (state%grid%zface(cellk) - pos%z)/dir%z
        elseif(dir%z == 0._wp)then
            dz = 100000._wp
        end if

        res = min(dx, dy, dz)
        if(res < 0._wp)then
            print*,'dcell < 0.0 warning! ',res
            print*,dx,dy,dz
            print*,dir
            print*,celli,cellj,cellk
            error stop 1
        end if

        if(res == dx)ldir = [.TRUE., .FALSE., .FALSE.]
        if(res == dy)ldir = [.FALSE., .TRUE., .FALSE.]
        if(res == dz)ldir = [.FALSE., .FALSE., .TRUE.]
        if(.not.ldir(1) .and. .not.ldir(2) .and. .not.ldir(3))print*,'Error in dir flag'
      
   end function wall_dist


    subroutine update_pos(pos, celli, cellj, cellk, dcell, wall_flag, dir, ldir, delta)
    !routine that upates postions of photon and calls fresnel routines if photon leaves current voxel
    !
    !
        use vector_class
        use sim_state_mod, only : state
        use utils,         only : str

        implicit none
      
        type(vector),    intent(IN)    :: dir
        logical,         intent(IN)    :: wall_flag, ldir(:)
        real(kind=wp),   intent(IN)    :: dcell, delta
        type(vector),    intent(INOUT) :: pos
        integer,         intent(INOUT) :: celli, cellj, cellk

        if(wall_flag)then

            if(ldir(1))then
                if(dir%x > 0._wp)then
                    pos%x = state%grid%xface(celli+1) + delta
                elseif(dir%x < 0._wp)then
                    pos%x = state%grid%xface(celli) - delta
                else
                    print*,'Error in x ldir in update_pos', ldir, dir
                end if
                pos%y = pos%y + dir%y*dcell 
                pos%z = pos%z + dir%z*dcell
            elseif(ldir(2))then
                if(dir%y > 0._wp)then
                    pos%y = state%grid%yface(cellj+1) + delta
                elseif(dir%y < 0._wp)then
                    pos%y = state%grid%yface(cellj) - delta
                else
                    print*,'Error in y ldir in update_pos', ldir, dir
                end if
                pos%x = pos%x + dir%x*dcell
                pos%z = pos%z + dir%z*dcell
            elseif(ldir(3))then
                if(dir%z > 0._wp)then
                    pos%z = state%grid%zface(cellk+1) + delta
                elseif(dir%z < 0._wp)then
                    pos%z = state%grid%zface(cellk) - delta
                else
                    print*,'Error in z ldir in update_pos', ldir, dir
                end if
                pos%x = pos%x + dir%x*dcell
                pos%y = pos%y + dir%y*dcell 
            else
                print*,'Error in update_pos... '//str(ldir)
                error stop 1
            end if
        else
            pos%x = pos%x + dir%x*dcell
            pos%y = pos%y + dir%y*dcell 
            pos%z = pos%z + dir%z*dcell
        end if

        if(wall_flag)then
            call update_voxels(pos, celli, cellj, cellk)
        end if

    end subroutine update_pos


    subroutine update_voxels(pos, celli, cellj, cellk)
    !updates the current voxel based upon position
    !
    !
        use vector_class
        use sim_state_mod, only : state

        implicit none
        
        type(vector),    intent(IN)    :: pos
        integer,         intent(INOUT) :: celli, cellj, cellk

        celli = find(pos%x, state%grid%xface) 
        cellj = find(pos%y, state%grid%yface)
        cellk = find(pos%z, state%grid%zface) 
        ! celli = floor(state%grid%nxg * (pos%x) / (2. * state%grid%xmax)) + 1
        ! cellj = floor(state%grid%nyg * (pos%y) / (2. * state%grid%ymax)) + 1
        ! cellk = floor(state%grid%nzg * (pos%z) / (2. * state%grid%zmax)) + 1

        if(celli > state%grid%nxg .or. celli < 1)celli = -1
        if(cellj > state%grid%nyg .or. cellj < 1)cellj = -1
        if(cellk > state%grid%nzg .or. cellk < 1)cellk = -1

    end subroutine update_voxels

    integer function find(val, a)
    !searchs for bracketing indicies for a value val in an array a
    !
    !
        implicit none

        real(kind=wp), intent(IN) :: val, a(:)
        integer :: n, lo, mid, hi

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
    ! wrapper routine for fresnel calculation
    !
    !
        use random, only : ran2
        use vector_class, only : vector

        implicit none

        type(vector),  intent(INOUT) :: I !incident vector
        type(vector),  intent(INOUT) :: N ! normal vector
        real(kind=wp), intent(IN)    :: n1, n2 !refractive indcies
        logical,       intent(OUT)   :: rflag !reflection flag

        rflag = .FALSE.

        !draw random number, if less than fresnel coefficents, then reflect, else refract
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

        type(vector), intent(INOUT) :: I ! incident vector
        type(vector), intent(IN)    :: N ! normal vector

        type(vector) :: R

        R = I - 2._wp * (N .dot. I) * N
        I = R

    end subroutine reflect


    subroutine refract(I, N, eta)
    !   get vector of refracted photon
    !
    !
        use vector_class

        implicit none

        type(vector),  intent(INOUT) :: I
        type(vector),  intent(IN)    :: N
        real(kind=wp), intent(IN)    :: eta

        type(vector) :: T, Ntmp

        real(kind=wp) :: c1, c2

        Ntmp = N

        c1 = (Ntmp .dot. I)
        if(c1 < 0._wp)then
            c1 = -c1
        else
            Ntmp = (-1._wp) * N
        end if
        c2 = sqrt(1._wp - (eta)**2 * (1._wp-c1**2))

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

        real(kind=wp), intent(IN) :: n1, n2
        type(vector),  intent(IN) :: I, N

        real(kind=wp) ::  costt, sintt, sint2, cost2, tir, f1, f2

        costt = abs(I .dot. N)

        sintt = sqrt(1._wp - costt * costt)
        sint2 = n1/n2 * sintt
        if(sint2 > 1._wp)then
            tir = 1.0_wp
            return
        elseif(costt == 1._wp)then
            tir = 0._wp
            return
        else
            sint2 = (n1/n2)*sintt
            cost2 = sqrt(1._wp - sint2 * sint2)
            f1 = abs((n1*costt - n2*cost2) / (n1*costt + n2*cost2))**2
            f2 = abs((n1*cost2 - n2*costt) / (n1*cost2 + n2*costt))**2

            tir = 0.5_wp * (f1 + f2)
        if(ieee_is_nan(tir) .or. tir > 1._wp .or. tir < 0._wp)print*,'TIR: ', tir, f1, f2, costt,sintt,cost2,sint2
            return
        end if
    end function fresnel
end module inttau2