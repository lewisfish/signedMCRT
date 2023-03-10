module inttau2

    use constants, only : wp

    implicit none
    
    private
    public :: tauint2

    contains   

    subroutine tauint2(grid, packet, sdfs_array)
    ! optical depth integration subroutine
    ! Moves photons to interaction location
    ! Calculated is any reflection or refraction happens whilst moving
    !
        use random,       only : ran2
        use photonMod,    only : photon
        use gridMod,      only : cart_grid
        use vector_class, only : vector
        use surfaces,     only : reflect_refract
        use sdfs
   
        type(cart_grid),   intent(in)    :: grid
        type(photon),      intent(inout) :: packet
        type(container),   intent(in)    :: sdfs_array(:)

        real(kind=wp) :: tau, d_sdf, t_sdf, taurun, ds(size(sdfs_array)), dstmp(size(sdfs_array))
        real(kind=wp) :: eps, dtot, old(size(sdfs_array)), new(size(sdfs_array)), n1, n2, Ri
        integer       :: i, oldlayer, new_layer
        type(vector)  :: pos, dir, oldpos, N
        logical       :: rflag

        !setup temp variables
        pos = packet%pos
        oldpos = pos
        dir = vector(packet%nxp, packet%nyp, packet%nzp)


        !round off distance
        eps = 1e-8_wp
        !get random tau
        tau = -log(ran2())
        taurun = 0.
        dtot = 0.
        do
            !setup sdf distance and current layer
            ds = 0.
            do i = 1, size(ds)
                ds(i) = abs(sdfs_array(i)%p%evaluate(pos))
            end do
            packet%cnts = packet%cnts + size(ds)
            d_sdf = minval(ds)

            ! if(d_sdf < eps)then
            !     packet%tflag=.true.
            !     exit
            ! end if

            do while(d_sdf > eps)
                t_sdf = d_sdf * sdfs_array(packet%layer)%p%kappa
                if(taurun + t_sdf <= tau)then
                    !move full distance to sdf surface
                    taurun = taurun + t_sdf
                    oldpos = pos
                    !comment out for phase screen
                    call update_phasor(grid, oldpos, dir, d_sdf, packet, sdfs_array(packet%layer)%p%mua)
                    pos = pos + d_sdf * dir
                    dtot = dtot + d_sdf
                else
                    !run out of tau so move remaining tau and exit
                    d_sdf = (tau - taurun) / sdfs_array(packet%layer)%p%kappa
                    dtot = dtot + d_sdf
                    taurun = tau
                    oldpos = pos
                    pos = pos + d_sdf * dir
                    !comment out for phase screen
                    call update_phasor(grid, oldpos, dir, d_sdf, packet, sdfs_array(packet%layer)%p%mua)
                    exit
                end if
                ! get distance to nearest sdf
                ds = 0._wp
                do i = 1, size(ds)
                    ds(i) = sdfs_array(i)%p%evaluate(pos)
                end do
                d_sdf = minval(abs(ds),dim=1)
                packet%cnts = packet%cnts + size(ds)

                !check if outside all sdfs
                if(minval(ds) >= 0._wp)then
                    packet%tflag = .true.
                    exit
                end if
            end do

            !exit early if conditions met
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
            oldpos = pos
            pos = pos + d_sdf*dir
            ds = 0._wp
            do i = 1, size(ds)
                ds(i) = sdfs_array(i)%p%evaluate(pos)
            end do
            packet%cnts = packet%cnts + size(ds)
            
            new = 0._wp
            old = 0._wp
            do i = 1, size(ds)
                if(dstmp(i) < 0.)then
                    old(i)=-1._wp
                    exit
                end if
            end do
            do i = 1, size(ds)
                if(ds(i) < 0.)then
                    new(i)=-1._wp     
                    exit
                end if
            end do

            !check for fresnel reflection
            n1 = sdfs_array(packet%layer)%p%n
            new_layer = minloc(new, dim=1)
            n2 = sdfs_array(new_layer)%p%n
            !carry out refelction/refraction
            if (n1 /= n2)then
                !get correct sdf normal
                if(ds(packet%layer) < 0._wp .and. ds(new_layer) < 0._wp)then
                    oldlayer = minloc(abs([ds(packet%layer), ds(new_layer)]), dim=1)
                elseif(dstmp(packet%layer) < 0._wp .and. dstmp(new_layer) < 0._wp)then
                    oldlayer=maxloc([dstmp(packet%layer), dstmp(new_layer)],dim=1)
                elseif(ds(packet%layer) > 0._wp .and. ds(new_layer) < 0._wp)then
                    oldlayer = packet%layer
                elseif(ds(packet%layer) > 0._wp .and. ds(new_layer) > 0._wp)then
                    packet%tflag = .true.
                    exit
                else
                    error stop "This should not be reached!"
                end if
                if(oldlayer == 1)then
                    oldlayer = packet%layer
                else
                    oldlayer = new_layer
                end if
                N = calcNormal(pos, sdfs_array(oldlayer)%p)

                rflag = .false.
                call reflect_refract(dir, N, n1, n2, rflag, Ri)
                packet%weight = packet%weight * Ri
                tau = -log(ran2())
                taurun = 0._wp
                if(.not.rflag)then
                    packet%layer = new_layer
                else
                    !step back inside original sdf
                    pos = oldpos
                    !reflect so incrment bounce counter
                    packet%bounces = packet%bounces + 1
                    if(packet%bounces > 1000)then
                        packet%tflag=.true.
                        exit
                    end if
                end if
            else
                packet%layer = new_layer
            end if
            if(packet%tflag)exit
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

        ! packet%step = dtot
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


    subroutine update_phasor(grid, pos, dir, d_sdf, packet, mua)
    ! record fluence using path length estimators. Uses voxel grid
    ! grid stores voxel grid information (voxel walls and etc)
    ! pos is current position with origin in centre of medium (0,0,0)
    ! dir is the current direction (0,0,1) is up
    ! d_sdf is the distance to travel in voxel grid
    ! packet stores the photon related variables

        use vector_class
        use photonMod
        use gridMod
        use iarray,        only: phasor
        
        type(cart_grid), intent(IN)    :: grid
        type(vector),    intent(IN)    :: dir
        real(kind=wp),   intent(IN)    :: d_sdf
        real(kind=wp), optional, intent(IN) :: mua
        complex(kind=wp) :: phasec
        type(vector),    intent(INOUT) :: pos
        type(photon),    intent(INOUT) :: packet

        type(vector)  :: old_pos
        logical       :: ldir(3)
        integer       :: celli, cellj, cellk
        real(kind=wp) :: dcell, delta=1e-8_wp, d, mua_real

        if(present(mua))then
            mua_real = mua
        else
            mua_real = 1._wp
        end if

        !convert to different coordinate system. Origin is at lower left corner of fluence grid
        old_pos = vector(pos%x+grid%xmax, pos%y+grid%ymax, pos%z+grid%zmax)
        call update_voxels(grid, old_pos, celli, cellj, cellk)
        packet%xcell = celli
        packet%ycell = cellj
        packet%zcell = cellk

        d = 0._wp
        !if packet outside grid return
        if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
            packet%tflag = .true.
            pos = vector(old_pos%x-grid%xmax, old_pos%y-grid%ymax, old_pos%z-grid%zmax)
            return
        end if
        !move photon through grid updating path length estimators
        do
            ldir = (/.FALSE., .FALSE., .FALSE./)

            dcell = wall_dist(grid, celli, cellj, cellk, old_pos, dir, ldir)
            if(d + dcell > d_sdf)then
                dcell = d_sdf - d
                d = d_sdf
! needs to be atomic so dont write to same array address with more than 1 thread at a time
                    packet%phase = packet%phase + dcell
                    phasec = packet%energy*cmplx(cos(packet%fact*packet%phase), sin(packet%fact*packet%phase))
!$omp atomic
                    phasor(celli,cellj,cellk) = phasor(celli,cellj,cellk) + phasec!*mua_real
                    !jmean(celli, cellj, cellk) = jmean(celli, cellj, cellk) + dcell*mua_real
                call update_pos(grid, old_pos, celli, cellj, cellk, dcell, .false., dir, ldir, delta)
                exit
            else
                d = d + dcell
                    packet%phase = packet%phase + dcell
                    phasec = packet%energy*cmplx(cos(packet%fact*packet%phase), sin(packet%fact*packet%phase))
!$omp atomic
                    phasor(celli,cellj,cellk) = phasor(celli,cellj,cellk) + phasec!*mua_real
                    !jmean(celli, cellj, cellk) = jmean(celli, cellj, cellk) + dcell*mua_real
                call update_pos(grid, old_pos, celli, cellj, cellk, dcell, .true., dir, ldir, delta)
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

    end subroutine update_phasor

    function wall_dist(grid, celli, cellj, cellk, pos, dir, ldir) result(res)
    !funtion that returns distant to nearest wall and which wall that is (x, y, or z)
    !
    !
        use vector_class
        use gridMod

        type(cart_grid), intent(IN)    :: grid
        type(vector),    intent(IN)    :: pos, dir
        logical,         intent(INOUT) :: ldir(:)
        integer,         intent(INOUT) :: celli, cellj, cellk
        real(kind=wp) :: res

        real(kind=wp) :: dx, dy, dz

        dx = -999._wp
        dy = -999._wp
        dz = -999._wp

        if(dir%x > 0._wp)then
            dx = (grid%xface(celli+1) - pos%x)/dir%x
        elseif(dir%x < 0._wp)then
            dx = (grid%xface(celli) - pos%x)/dir%x
        elseif(dir%x == 0._wp)then
            dx = 100000._wp
        end if

        if(dir%y > 0._wp)then
            dy = (grid%yface(cellj+1) - pos%y)/dir%y
        elseif(dir%y < 0._wp)then
            dy = (grid%yface(cellj) - pos%y)/dir%y
        elseif(dir%y == 0._wp)then
            dy = 100000._wp
        end if

        if(dir%z > 0._wp)then
            dz = (grid%zface(cellk+1) - pos%z)/dir%z
        elseif(dir%z < 0._wp)then
            dz = (grid%zface(cellk) - pos%z)/dir%z
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

        ldir = [res == dx, res==dy, res==dz]
        if(.not.ldir(1) .and. .not.ldir(2) .and. .not.ldir(3))print*,'Error in dir flag'
      
   end function wall_dist


    subroutine update_pos(grid, pos, celli, cellj, cellk, dcell, wall_flag, dir, ldir, delta)
    !routine that updates positions of photon and calls Fresnel routines if photon leaves current voxel
    !
    !
        use vector_class
        use gridMod
        use utils, only : str
      
        type(cart_grid), intent(IN)    :: grid
        type(vector),    intent(IN)    :: dir
        logical,         intent(IN)    :: wall_flag, ldir(:)
        real(kind=wp),   intent(IN)    :: dcell, delta
        type(vector),    intent(INOUT) :: pos
        integer,         intent(INOUT) :: celli, cellj, cellk

        if(wall_flag)then

            if(ldir(1))then
                if(dir%x > 0._wp)then
                    pos%x = grid%xface(celli+1) + delta
                elseif(dir%x < 0._wp)then
                    pos%x = grid%xface(celli) - delta
                else
                    print*,'Error in x ldir in update_pos', ldir, dir
                end if
                pos%y = pos%y + dir%y*dcell 
                pos%z = pos%z + dir%z*dcell
            elseif(ldir(2))then
                if(dir%y > 0._wp)then
                    pos%y = grid%yface(cellj+1) + delta
                elseif(dir%y < 0._wp)then
                    pos%y = grid%yface(cellj) - delta
                else
                    print*,'Error in y ldir in update_pos', ldir, dir
                end if
                pos%x = pos%x + dir%x*dcell
                pos%z = pos%z + dir%z*dcell
            elseif(ldir(3))then
                if(dir%z > 0._wp)then
                    pos%z = grid%zface(cellk+1) + delta
                elseif(dir%z < 0._wp)then
                    pos%z = grid%zface(cellk) - delta
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
            call update_voxels(grid, pos, celli, cellj, cellk)
        end if

    end subroutine update_pos


    subroutine update_voxels(grid, pos, celli, cellj, cellk)
    ! updates the current voxel based upon position
    !
    !
        use vector_class
        use gridmod
        
        type(cart_grid), intent(IN)    :: grid
        type(vector),    intent(IN)    :: pos
        integer,         intent(INOUT) :: celli, cellj, cellk

        !accurate but slow
        ! celli = find(pos%x, grid%xface) 
        ! cellj = find(pos%y, grid%yface)
        ! cellk = find(pos%z, grid%zface) 

        !fast but can be inaccurate in some cases...
        celli = floor(grid%nxg * (pos%x) / (2. * grid%xmax)) + 1
        cellj = floor(grid%nyg * (pos%y) / (2. * grid%ymax)) + 1
        cellk = floor(grid%nzg * (pos%z) / (2. * grid%zmax)) + 1

        if(celli > grid%nxg .or. celli < 1)celli = -1
        if(cellj > grid%nyg .or. cellj < 1)cellj = -1
        if(cellk > grid%nzg .or. cellk < 1)cellk = -1

    end subroutine update_voxels

    integer function find(val, a)
    ! searches for bracketing indices for a value value in an array a
    !
    !
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
end module inttau2