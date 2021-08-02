module inttau2

   implicit none
   
   private
   public :: tauint2

CONTAINS

!     subroutine tauint1(packet, grid)
!     !optical depth integration subroutine
!     !
!     !
!         use iarray, only : jmean
!         use random, only : ran2
!         use photonMod
!         use gridMod

!         use vector_class
   
!         implicit none

!         type(cart_grid), intent(IN)    :: grid
!         type(photon),    intent(INOUT) :: packet

!         real    :: tau, taurun, taucell, xcur, ycur, zcur, d, dcell
!         integer :: celli, cellj, cellk
!         logical :: dir(3)

!         xcur = packet%xp + grid%xmax
!         ycur = packet%yp + grid%ymax
!         zcur = packet%zp + grid%zmax

!         celli = packet%xcell
!         cellj = packet%ycell
!         cellk = packet%zcell

!         taurun = 0.
!         d = 0.
!         dir = (/.FALSE., .FALSE., .FALSE./)

!         !sample optical distance
!         tau = -log(ran2())
!         do
!             dir = (/.FALSE., .FALSE., .FALSE./)
!             !get distance to nearest wall in direction dir
!             dcell = wall_dist(packet, grid, celli, cellj, cellk, xcur, ycur, zcur, dir)
!             !calculate optical distnace to cell wall
!             taucell = dcell * grid%rhokap(celli,cellj,cellk)

!             if(taurun + taucell < tau)then!still some tau to move
!                 taurun = taurun + taucell
!                 d = d + dcell
! !$omp atomic
!                 jmean(celli, cellj, cellk) = jmean(celli, cellj, cellk) + dcell ! record fluence

!                 call update_pos(packet, grid, xcur, ycur, zcur, celli, cellj, cellk, dcell, .TRUE., dir, &
!                                 tau, taurun)
!             else!moved full distance

!                 dcell = (tau - taurun) / grid%rhokap(celli,cellj,cellk)
!                 d = d + dcell
! !$omp atomic                
!                 jmean(celli, cellj, cellk) = jmean(celli, cellj, cellk) + dcell ! record fluence

!                 call update_pos(packet, grid, xcur, ycur, zcur, celli, cellj, cellk, dcell, .FALSE., dir, &
!                                 tau, taurun)
!                 exit
!             end if
!             if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
!                 if(celli == -1 .or. cellj == -1)then
!                     call repeat_bounds(celli, cellj, xcur, ycur, grid%xmax, grid%ymax, grid%nxg, grid%nyg, grid%delta)
!                     packet%tflag = .false.
!                     if(celli == -1 .or. cellj == -1 .or. packet%tflag)then
!                        print*,'error',celli,cellj,packet%tflag
!                     end if
!                 else
!                     packet%tflag = .true.
!                     exit
!                 end if
!             end if
!         end do
   
!         packet%xp = xcur - grid%xmax
!         packet%yp = ycur - grid%ymax
!         packet%zp = zcur - grid%zmax
!         packet%xcell = celli
!         packet%ycell = cellj
!         packet%zcell = cellk

!     end subroutine tauint1
   

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

        type(istack) :: layer

        real    :: tau, d_sdf, t_sdf, taurun, ds(size(sdfs_array)), eps, dtot
        integer :: i, cur_layer, idxs(size(sdfs_array))
        type(vector) :: pos, dir, oldpos

        pos = packet%pos
        oldpos = pos
        dir = vector(packet%nxp, packet%nyp, packet%nzp)

        ds = 0.
        do i = 1, size(ds)
            ds(i) = abs(sdfs_array(i)%p%evaluate(pos))
        end do

        idxs = argsort(ds)
        do i = size(idxs), 1, -1
            call layer%push(idxs(i))
        end do

        eps = 1d-8
        tau = -log(ran2())
        taurun = 0.

        cur_layer = layer%pop()
        d_sdf = abs(sdfs_array(cur_layer)%p%evaluate(pos))
        dtot = 0.
        do while(packet%tflag .eqv. .false.)
            do while(d_sdf > eps)
                t_sdf = d_sdf * sdfs_array(packet%layer)%p%kappa
                if(taurun + t_sdf <= tau)then
                    taurun = taurun + t_sdf
                    pos = pos + d_sdf * dir
                    dtot = dtot + d_sdf
                else
                    d_sdf = (tau - taurun) / sdfs_array(packet%layer)%p%kappa
                    dtot = dtot + d_sdf
                    taurun = tau
                    pos = pos + d_sdf * dir
                    exit
                end if
                ds = 0.
                do i = 1, size(ds)
                    ds(i) = abs(sdfs_array(i)%p%evaluate(pos))
                end do
                d_sdf = minval(ds,dim=1)
                if(abs(pos%x) > grid%xmax)then
                    packet%tflag = .true.
                    exit
                elseif(abs(pos%y) > grid%ymax)then
                    packet%tflag = .true.
                    exit
                elseif(abs(pos%z) > grid%zmax)then
                    packet%tflag = .true.
                    exit
                end if
            end do
            ! print*,pos
            if(taurun >= tau)exit

            if(layer%peek() == -99)then
                ds = 0.
                do i = 1, size(ds)
                    ds(i) = abs(sdfs_array(i)%p%evaluate(pos))
                end do

                idxs = argsort(ds)
                do i = size(idxs), 1, -1
                    call layer%push(idxs(i))
                end do

                cur_layer = layer%pop()
                d_sdf = 10.*eps
            else
                cur_layer = layer%pop()
                ds = 0.
                do i = 1, size(ds)
                    ds(i) = abs(sdfs_array(i)%p%evaluate(pos))
                end do
                d_sdf = minval(ds,dim=1)
            end if
                packet%layer = cur_layer
        end do

        call update_jmean(oldpos, grid, dir, dtot, packet)
        packet%pos = pos

        if(abs(packet%pos%x) > grid%xmax)then
            packet%tflag = .true.
        elseif(abs(packet%pos%y) > grid%ymax)then
            packet%tflag = .true.
        elseif(abs(packet%pos%z) > grid%zmax)then
            packet%tflag = .true.
        end if
    end subroutine tauint2
   

    subroutine update_jmean(pos, grid, dir, d_sdf, packet)
        
        use vector_class
        use iarray, only: jmean
        use gridMod
        use photonMod

        implicit none
        
        type(cart_grid), intent(IN) :: grid
        type(vector), intent(INOUT) :: pos, dir
        real,         intent(IN) :: d_sdf
        type(photon), intent(INOUT) :: packet
        type(vector) :: old_pos
        logical :: ldir(3)
        integer :: celli, cellj, cellk
        real :: dcell, delta=1d-8,d

        old_pos = vector(pos%x+grid%xmax, pos%y+grid%ymax, pos%z+grid%zmax)
        celli = packet%xcell        
        cellj = packet%ycell
        cellk = packet%zcell
        d = 0.
        do
            if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
                packet%tflag = .true.
                exit
            end if
            ldir = (/.FALSE., .FALSE., .FALSE./)
            ! print*,celli,cellj,cellk
            dcell = wall_dist(grid, celli, cellj, cellk, old_pos, dir, ldir)
            if(d + dcell > d_sdf)then
                dcell = d_sdf - d
                d = d_sdf
                jmean(celli, cellj, cellk) = jmean(celli, cellj, cellk) + dcell
                call update_pos(old_pos, grid, celli, cellj, cellk, dcell, .false., dir, ldir, delta)
                exit
            else
                d = d + dcell
                jmean(celli, cellj, cellk) = jmean(celli, cellj, cellk) + dcell
                call update_pos(old_pos, grid, celli, cellj, cellk, dcell, .true., dir, ldir, delta)
            end if
        end do

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
            ! print*,'dcell < 0.0 warning! ',wall_dist,dx,dy,dz,dir
            ! error stop 1
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
      
        type(cart_grid), intent(IN) :: grid
        type(vector),    intent(INOUT) :: pos, dir
        real,    intent(IN)    :: dcell, delta
        integer, intent(INOUT) :: celli, cellj, cellk
        logical, intent(IN)    :: wall_flag, ldir(:)
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

        celli = floor(grid%nxg * (pos%x) / (2. * grid%xmax)) + 1
        cellj = floor(grid%nyg * (pos%y) / (2. * grid%ymax)) + 1
        cellk = floor(grid%nzg * (pos%z) / (2. * grid%zmax)) + 1

        if(celli > grid%nxg .or. celli < 1)celli = -1
        if(cellj > grid%nyg .or. cellj < 1)cellj = -1
        if(cellk > grid%nzg .or. cellk < 1)cellk = -1

    end subroutine update_voxels


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

    function argsort(a) result(b)

    real, intent(IN) :: a(:)   ! array of numbers
    
    integer :: b(size(a))         ! indices into the array 'a' that sort it

    integer :: N, i, imin, temp1                           ! number of numbers/vectors
    real    :: temp2, a2(size(a))

    a2 = a
    N = size(a)

    do i = 1, N
        b(i) = i
    end do

    do i = 1, N-1
        ! find ith smallest in 'a'
        imin = minloc(a2(i:),1) + i - 1
        ! swap to position i in 'a' and 'b', if not already there
        if (imin /= i) then
            temp2 = a2(i); a2(i) = a2(imin); a2(imin) = temp2
            temp1 = b(i); b(i) = b(imin); b(imin) = temp1
        end if
    end do
    end function argsort
end module inttau2