module writer_mod

implicit none

    contains
        subroutine writer()

            use constants, only : fileplace
            use iarray,    only : jmeanGLOBAL
            use gridMod

            implicit none

            integer :: u

            ! jmeanGLOBAL =jmeanGLOBAL * ((2.*grid%xmax)**2./(nphotons*numproc&
            !     *(2.*grid%xmax/grid%nxg)*(2.*grid%ymax/grid%nyg)*(2.*grid%zmax/grid%nzg)))

            open(newunit=u,file=trim(fileplace)//'jmean/jmean.dat',access='stream',status='REPLACE',form='unformatted')
            write(u) jmeanGLOBAL
            close(u)

        end subroutine writer
end module writer_mod
