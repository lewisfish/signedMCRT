module writer_mod

implicit none

    contains
        subroutine writer(nphotons, grid, optprop)

            use constants, only : fileplace
            use iarray,    only : jmeanGLOBAL
            use utils,     only : str
            use gridMod

            implicit none

            type(cart_grid), intent(IN) :: grid
            integer,         intent(IN) :: nphotons
            real,            intent(IN) :: optprop(:)

            real :: jm(grid%nzg), xmax, ymax, zmax
            integer :: u, i, j, k
            integer :: nxg, nyg, nzg
            character(len=:), allocatable :: filename

            nxg = grid%nxg
            nyg = grid%nyg
            nzg = grid%nzg
            xmax = grid%xmax
            ymax = grid%ymax
            zmax = grid%zmax


    jmeanGLOBAL  = jmeanGLOBAL * ((2.*xmax*2.*ymax)/(nphotons * (2. * xmax / nxg) * (2. * ymax / nyg) * (2. * zmax / nzg)))

            filename = trim(fileplace)//"jmean/bunny-point.dat"

            open(newunit=u,file=filename,access='stream',status='REPLACE',form='unformatted')
            write(u) jmeanGLOBAL !/ nphotons
            close(u)

            ! jm = 0.d0
            ! do k = 1, nzg
            !     do j = 1, nyg
            !         do i = 1, nxg
            !             jm(k) = jm(k) + jmeanGLOBAL(i,j,k)
            !         end do
            !     end do
            ! end do

            ! jm = jm / (nxg*nyg)

            ! open(newunit=u,file=trim(fileplace)//"jmean/validation-420.dat",status="replace")
            ! do i = nzg,1,-1
            !     write(u,*)real(nzg-i)*(2.*zmax/(nzg)),jm(i)
            ! end do
            ! close(u)


        end subroutine writer
end module writer_mod
