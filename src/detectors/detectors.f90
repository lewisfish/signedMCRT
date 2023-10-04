module detectors

    !! Module contains each detector type which inherits from the base detector class.
    !! detectors detect photon packets colliding with the detectors.

    use constants,    only : wp
    use detector_mod, only : detector, detector1D, detector2D, hit_t
    use vector_class, only : vector, length

    implicit none
    
    !> Circle detector
    type, extends(detector1D) :: circle_dect
        !> Radius of detector
        real(kind=wp) :: radius
    contains
        procedure :: check_hit  => check_hit_circle
    end type circle_dect

    interface circle_dect
        module procedure init_circle_dect
    end interface circle_dect

    !> Annuluar detector
    type, extends(detector1D) :: annulus_dect
        !> Inner radius
        real(kind=wp) :: r1
        !> Outer radius
        real(kind=wp) :: r2
        contains
        procedure :: check_hit => check_hit_annulus
    end type annulus_dect

    interface annulus_dect
        module procedure init_annulus_dect
    end interface annulus_dect

    !> Rectangular or "camera" detector
    type, extends(detector2D) :: camera
        !> Normal of the detector
        type(vector)  :: n
        !> Vector from pos (1st corner) to the 2nd corner of the detector
        type(vector)  :: p2
        !> Vector from pos (1st corner) to the 3rd corner of the detector
        type(vector)  :: p3
        !> Edge vector of detector
        type(vector)  :: e1
        !> Edge vector of detector
        type(vector)  :: e2
        !> Width of the detector
        real(kind=wp) :: width
        !> Height of the detector
        real(kind=wp) :: height
        contains
        procedure :: check_hit => check_hit_camera
    end type camera

    interface camera
        module procedure init_camera
    end interface camera
    
    !> Detector array
    type :: dect_array
        class(detector), pointer :: p => null()
    end type dect_array

    private
    public :: camera, annulus_dect, circle_dect, dect_array

    contains
    
    function init_circle_dect(pos, dir, layer, radius, nbins, maxval, trackHistory) result(out)
        !> Centre of detector
        type(vector),  intent(in) :: pos
        !> Normal of the detector
        type(vector),  intent(in) :: dir
        !> Layer ID
        integer,       intent(in) :: layer
        !> Number of bins in the detector
        integer,       intent(in) :: nbins
        !> Radius of the detector
        real(kind=wp), intent(in) :: radius
        !> Maximum value to store in bins
        real(kind=wp), intent(in) :: maxval
        !> Boolean on if to store photon's history prior to hitting the detector.
        logical,       intent(in) :: trackHistory

        type(circle_dect) :: out

        out%dir = dir
        out%pos = pos
        out%layer = layer
        !extra bin for data beyond end of array
        out%nbins = nbins + 1
        out%radius = radius
        allocate(out%data(out%nbins))
        out%data = 0.0_wp
        if(nbins == 0)then
            out%bin_wid = 1._wp
        else
            out%bin_wid = maxval / real(nbins-1, kind=wp)
        end if
        out%trackHistory = trackHistory

    end function init_circle_dect

    logical function check_hit_circle(this, hitpoint)
        
        use geometry, only : intersectCircle

        class(circle_dect), intent(INOUT) :: this
        type(hit_t),        intent(IN)    :: hitpoint
        
        real(kind=wp) :: t 

        check_hit_circle = .false.
        if(this%layer /= hitpoint%layer)return
        check_hit_circle = intersectCircle(this%dir, this%pos, this%radius, hitpoint%pos, hitpoint%dir, t)
        if(check_hit_circle)then
            if(t > 5e-3_wp)check_hit_circle=.false.
        end if
    end function check_hit_circle

    function init_annulus_dect(pos, dir, layer, r1, r2, nbins, maxval, trackHistory) result(out)
        !> Centre of detector
        type(vector),  intent(in) :: pos
        !> Normal of the detector
        type(vector),  intent(in) :: dir
        !> Layer ID
        integer,       intent(in) :: layer
        !> Inner radius
        real(kind=wp), intent(IN) :: r1
        !> Outer radius
        real(kind=wp), intent(IN) :: r2
        !> Number of bins in the detector
        integer,       intent(in) :: nbins
        !> Maximum value to store in bins
        real(kind=wp), intent(in) :: maxval
        !> Boolean on if to store photon's history prior to hitting the detector.
        logical,       intent(in) :: trackHistory

        type(annulus_dect) :: out

        out%pos = pos
        out%dir = dir
        out%layer = layer
        !extra bin for data beyond end of array
        out%nbins = nbins + 1
        out%r1 = r1
        out%r2 = r2
        allocate(out%data(out%nbins))
        out%data = 0.0_wp
        if(nbins == 0)then
            out%bin_wid = 1._wp
        else
            out%bin_wid = maxval / real(nbins, kind=wp)
        end if
        out%trackHistory = trackHistory

    end function init_annulus_dect

    logical function check_hit_annulus(this, hitpoint)

        class(annulus_dect), intent(INOUT) :: this
        type(hit_t),         intent(IN)    :: hitpoint

        real(kind=wp) :: newpos

        check_hit_annulus = .false.
        if(this%layer /= hitpoint%layer)return
        newpos = sqrt((hitpoint%pos%x - this%pos%x)**2 + (hitpoint%pos%y - this%pos%y)**2 + (hitpoint%pos%z - this%pos%z)**2)
        if(newpos >= this%r1 .and. newpos <= this%r2)then
            check_hit_annulus = .true.
        end if

    end function check_hit_annulus

    function init_camera(p1, p2, p3, layer, nbins, maxval, trackHistory) result(out)
        !> Position of the 1st corner of the detector
        type(vector),  intent(in) :: p1
        !> Distance from p1 to the 2nd corner
        type(vector),  intent(in) :: p2
        !> Distance from p1 to the 3rd corner
        type(vector),  intent(in) :: p3
        !> Layer ID
        integer,       intent(in) :: layer
        !> Number of bins in the detector
        integer,       intent(in) :: nbins
        !> Maximum value to store in bins
        real(kind=wp), intent(in) :: maxval
        !> Boolean on if to store photon's history prior to hitting the detector.
        logical,       intent(in) :: trackHistory
        type(camera) :: out

        out%pos = p1
        out%p2 = p2
        out%p3 = p3
        out%e1 = p2 - p1
        out%e2 = p3 - p1
        out%width = length(out%e1)
        out%height = length(out%e2)
        out%n = out%e2 .cross. out%e1
        out%n = out%n%magnitude()
        out%layer = layer
        !extra bin for data beyond end of array
        out%nbinsX = nbins + 1
        out%nbinsY = nbins + 1
        allocate(out%data(out%nbinsX, out%nbinsY))
        out%data = 0.0_wp
        if(nbins == 0)then
            out%bin_wid_x = 1._wp
            out%bin_wid_y = 1._wp
        else
            out%bin_wid_x = maxval / real(out%nbinsX, kind=wp)
            out%bin_wid_y = maxval / real(out%nbinsY, kind=wp)
        end if
        out%trackHistory = trackHistory

    end function init_camera

    logical function check_hit_camera(this, hitpoint)
    !! [ref](https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection)
        class(camera), intent(inout) :: this
        type(hit_t),   intent(in)    :: hitpoint

        real(kind=wp) :: t, proj1, proj2
        type(vector)  :: v

        check_hit_camera = .false.
        if(this%layer /= hitpoint%layer)return

        t = ((this%pos - hitpoint%pos) .dot. this%n) / (hitpoint%dir .dot. this%n)
        if(t >= 0._wp)then
            v = (hitpoint%pos + t * hitpoint%dir) - this%pos
            proj1 = (v .dot. this%e1) / this%width
            proj2 = (v .dot. this%e2) / this%height
            if((proj1 < this%width .and. proj1 > 0._wp) .and. (proj2 < this%height .and. proj2 > 0._wp))then
                check_hit_camera = .true.
            end if
        end if
    end function check_hit_camera

end module detectors