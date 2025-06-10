!> Module for handling problematic values in simulations, such as negative
!> pressures
module mod_small_values

  implicit none
  private

  !> How to handle small values
  character(len=20), public :: small_values_method = "error"

  !> Average over this many cells in each direction
  integer, public :: small_values_daverage = 1

  !> trace small values in the source file using traceback flag of compiler
  logical, public :: trace_small_values=.false.

  !> Whether to apply small value fixes to certain variables
  logical, public, allocatable :: small_values_fix_iw(:)

  public :: small_values_error
  public :: small_values_average

contains

  subroutine small_values_error(wprim, x, ixImin1,ixImax1, ixOmin1,ixOmax1,&
      w_flag, subname)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: wprim(ixImin1:ixImax1, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1, 1:ndim)
    logical, intent(in)          :: w_flag(ixImin1:ixImax1,1:nw)
    character(len=*), intent(in) :: subname
    integer                      :: iw,iiw,ix1

    if (.not.crash) then
      do iw=1,nw
       do ix1= ixOmin1,ixOmax1
          if(w_flag(ix1,iw)) then
            write(*,*) "Error: small value of ", trim(prim_wnames(iw)),&
               wprim(ix1,iw)," encountered when call ", subname
            write(*,*) "Iteration: ", it, " Time: ", global_time,&
                "Processor: ",mype
            write(*,*) "Location: ", x(ix1,:)
            write(*,*) "Cell number: ", ix1
            do iiw=1,nw
              write(*,*) trim(prim_wnames(iiw)),": ",wprim(ix1,iiw)
            end do
            ! use erroneous arithmetic operation to crash the run
            if(trace_small_values) write(*,*) sqrt(wprim(ix1,iw)-bigdouble)
            write(*,*) "Saving status at the previous time step"
            crash=.true.
          end if
       enddo
      end do
    end if

  end subroutine small_values_error

  subroutine small_values_average(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x,&
      w_flag, windex)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    logical, intent(in)             :: w_flag(ixImin1:ixImax1,1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1, 1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)
    integer, optional, intent(in)   :: windex
    integer                         :: iw, kxOmin1,kxOmax1, ix1, i, nwstart,&
        nwend

    if(present(windex)) then
      nwstart=windex
      nwend=windex
    else
      nwstart=1
      nwend=nw
    end if

    do iw = nwstart, nwend
     do ix1= ixOmin1,ixOmax1
      ! point with local failure identified by w_flag
        if (w_flag(ix1,iw)) then
          ! verify in cube with border width small_values_daverage the presence of
          ! cells where all went ok
          do i = 1, max(small_values_daverage, 1)
            kxOmin1= max(ix1-i, ixImin1);
            kxOmax1= min(ix1+i, ixImax1);
            ! in case cells are fine within smaller cube than 
            ! the userset small_values_daverage: use that smaller cube
            if(any(w_flag(kxOmin1:kxOmax1,iw) .eqv. .false.)) exit
          end do

          if(any(w_flag(kxOmin1:kxOmax1,iw) .eqv. .false.)) then
            ! within surrounding cube, cells without problem were found

            ! faulty cells are corrected by averaging here
            ! only average those which were ok and replace faulty cells
            if(small_values_fix_iw(iw)) then
              w(ix1, iw) = sum(w(kxOmin1:kxOmax1, iw), w_flag(kxOmin1:kxOmax1,&
                 iw) .eqv. .false.)/ count(w_flag(kxOmin1:kxOmax1,&
                 iw) .eqv. .false.)
            end if
          else
            write(*,*) "no cells without error were found in cube of size",&
                small_values_daverage
            write(*,*) "at location:", x(ix1, 1:ndim)
            write(*,*) "at index:", ix1
            write(*,*) "w numer:", iw
            !write(*,*) "Saving status at the previous time step"
            !crash=.true.
            write(*,*) "replace with small values"
            if(iw==iw_e) w(ix1, iw)=small_pressure
            if(iw==iw_rho) w(ix1, iw)=small_density
          end if
        end if
     enddo
    end do

  end subroutine small_values_average

end module mod_small_values
