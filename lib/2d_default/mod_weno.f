module mod_weno
  ! All kinds of (W)ENO schemes
  !
  ! 2019.9.19 WENO(-JS)5 transplant from the BHAC code;
  ! 2019.9.20 WENO3;
  ! 2019.9.21 WENO-Z5;
  ! 2019.9.22 WENO-Z+5 transplant from the BHAC code;
  ! 2019.10.30 WENO(-JS)7;
  ! 2019.10.31 MPWENO7;
  ! 2019.11.1 exENO7;
  ! 2019.11.7 clean up the code, comment out the interpolation variation;
  ! 2019.12.9 WENO-YC3;
  ! 2020.1.15 new WENO5 limiter WENO5NM for stretched grid.
  ! 2020.4.15 WENO5-CU6: hybrid sixth-order linear & WENO5
  ! 2021.6.12 generic treatment for fixing unphysical reconstructions
  ! 2022.10.21 remove exENO7
   
  implicit none
  private

  double precision, parameter     :: weno_eps_machine = 1.0d-18
  double precision, parameter     :: weno_dx_exp = 2.0d0/3.0d0

  public :: WENO3limiter
  public :: WENO5limiter
  public :: WENO5NMlimiter
  public :: WENO5limiterL
  public :: WENO5NMlimiterL
  public :: WENO5limiterR
  public :: WENO5NMlimiterR
  public :: TENO5ADlimiter
  public :: WENO5CU6limiter
  public :: WENO7limiter
  public :: fix_limiter
  public :: fix_limiter1
  public :: fix_onelimiter
  public :: fix_onelimiter1

contains

  subroutine fix_onelimiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
     iLmax1,iLmax2,wCin,wCout)
    use mod_global_parameters
    use mod_physics, only: phys_check_w

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2
    double precision, intent(in)    :: wCin(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    double precision, intent(inout) :: wCout(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)

    integer :: iw
    logical :: flagC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        flag(ixImin1:ixImax1,ixImin2:ixImax2)

    ! When limiter not TVD, negative pressures or densities could result.
    ! Fall back to flat interpolation 
    ! flagC(*,iw) indicates failed state (T when failed)
    ! assumes wCin contains primitive variables
    call phys_check_w(.true.,ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
       iLmax1,iLmax2,wCin,flagC)

    ! collect all failures
    flag(iLmin1:iLmax1,iLmin2:iLmax2)=.false.
    do iw=1,nwflux
       flag(iLmin1:iLmax1,iLmin2:iLmax2)=flag(iLmin1:iLmax1,&
          iLmin2:iLmax2).or.(flagC(iLmin1:iLmax1,iLmin2:iLmax2,iw))
    end do
    ! only use WENO reconstructions when no field failed
    ! in other places: do not modify the initial state in wCout
    do iw=1,nwflux
       where (flag(iLmin1:iLmax1,iLmin2:iLmax2) .eqv. .false.)
          wCout(iLmin1:iLmax1,iLmin2:iLmax2,iw)=wCin(iLmin1:iLmax1,&
             iLmin2:iLmax2,iw)
       end where
    enddo
  
  end subroutine fix_onelimiter

  subroutine fix_onelimiter1(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
     iLmax1,iLmax2,wCin,wCout)
    use mod_global_parameters
    use mod_physics, only: phys_check_w

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2
    double precision, intent(in)    :: wCin(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    double precision, intent(inout) :: wCout(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)

    integer :: iw
    logical :: flagC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    ! When limiter not TVD, negative pressures or densities could result.
    ! Fall back to flat interpolation 
    ! flagC(*,iw) indicates failed state (T when failed)
    ! assumes wCin contains primitive variables
    call phys_check_w(.true.,ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
       iLmax1,iLmax2,wCin,flagC)

    ! only use WENO reconstructions when no field failed
    ! in other places: do not modify the initial state in wCout
    do iw=1,nwflux
       where (flagC(iLmin1:iLmax1,iLmin2:iLmax2,iw) .eqv. .false.)
          wCout(iLmin1:iLmax1,iLmin2:iLmax2,iw)=wCin(iLmin1:iLmax1,&
             iLmin2:iLmax2,iw)
       end where
    enddo
  
  end subroutine fix_onelimiter1

  subroutine fix_limiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
     iLmax2,wLCin,wRCin,wLCout,wRCout)
    use mod_global_parameters
    use mod_physics, only: phys_check_w

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2
    double precision, intent(in)    :: wRCin(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw),wLCin(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) 
    double precision, intent(inout) :: wRCout(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw),wLCout(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) 

    integer :: iw
    logical :: flagL(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        flagR(ixImin1:ixImax1,ixImin2:ixImax2,1:nw), flag(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    ! When limiter not TVD, negative pressures or densities could result.
    ! Fall back to flat interpolation 
    ! flagL(*,iw) indicates failed L state (T when failed)
    ! flagR(*,iw) indicates failed R state (T when failed)
    ! assumes wLCin and wRCin contain primitive variables
    call phys_check_w(.true.,ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
       iLmax1,iLmax2,wLCin,flagL)
    call phys_check_w(.true.,ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
       iLmax1,iLmax2,wRCin,flagR)

    ! collect all failures
    flag(iLmin1:iLmax1,iLmin2:iLmax2)=.false.
    do iw=1,nwflux
       flag(iLmin1:iLmax1,iLmin2:iLmax2)=flag(iLmin1:iLmax1,&
          iLmin2:iLmax2).or.(flagL(iLmin1:iLmax1,iLmin2:iLmax2,&
          iw).or.flagR(iLmin1:iLmax1,iLmin2:iLmax2,iw))
    end do
    ! only use WENO reconstructions L and R when no neighbour field failed
    ! in other places: do not modify the initial state in wLCout/wRCout
    do iw=1,nwflux
       where (flag(iLmin1:iLmax1,iLmin2:iLmax2) .eqv. .false.)
          wLCout(iLmin1:iLmax1,iLmin2:iLmax2,iw)=wLCin(iLmin1:iLmax1,&
             iLmin2:iLmax2,iw)
          wRCout(iLmin1:iLmax1,iLmin2:iLmax2,iw)=wRCin(iLmin1:iLmax1,&
             iLmin2:iLmax2,iw)
       end where
    enddo
  
  end subroutine fix_limiter

  subroutine fix_limiter1(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
     iLmax2,wLCin,wRCin,wLCout,wRCout)
    use mod_global_parameters
    use mod_physics, only: phys_check_w

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2
    double precision, intent(in)    :: wRCin(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw),wLCin(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) 
    double precision, intent(inout) :: wRCout(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw),wLCout(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) 

    integer :: iw
    logical :: flagL(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        flagR(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    ! When limiter not TVD, negative pressures or densities could result.
    ! Fall back to flat interpolation 
    ! flagL(*,iw) indicates failed L state (T when failed)
    ! flagR(*,iw) indicates failed R state (T when failed)
    ! assumes wLCin and wRCin contain primitive variables
    call phys_check_w(.true.,ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
       iLmax1,iLmax2,wLCin,flagL)
    call phys_check_w(.true.,ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
       iLmax1,iLmax2,wRCin,flagR)

    ! only use WENO reconstructions L and R when no neighbour field failed
    ! in other places: do not modify the initial state in wLCout/wRCout
    do iw=1,nwflux
       where ((flagL(iLmin1:iLmax1,iLmin2:iLmax2,&
          iw) .eqv. .false.) .and. (flagR(iLmin1:iLmax1,iLmin2:iLmax2,&
          iw) .eqv. .false.))
          wLCout(iLmin1:iLmax1,iLmin2:iLmax2,iw)=wLCin(iLmin1:iLmax1,&
             iLmin2:iLmax2,iw)
          wRCout(iLmin1:iLmax1,iLmin2:iLmax2,iw)=wRCin(iLmin1:iLmax1,&
             iLmin2:iLmax2,iw)
       end where
    enddo
  
  end subroutine fix_limiter1

  subroutine WENO3limiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
     iLmax2,idims,dxdim,w,wLC,wRC,var)
    use mod_global_parameters
  
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2, idims, var
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw),wLC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) 
    !> local
    integer                         :: iLmmin1,iLmmin2,iLmmax1,iLmmax2,&
        iLpmin1,iLpmin2,iLpmax1,iLpmax2, iLppmin1,iLppmin2,iLppmax1,iLppmax2
    double precision                :: f_array(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw,2), d_array(2)
    double precision                :: beta(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw,2),tau(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision                :: u1_coeff(2), u2_coeff(2)
    double precision                :: alpha_array(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw,2), alpha_sum(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), flux(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer                         :: i

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)  :: wRCtmp, wLCtmp

    ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  
    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);
    d_array(1:2) = (/ 1.0d0/3.0d0, 2.0d0/3.0d0 /)
    u1_coeff(1:2) = (/ -1.d0/2.d0, 3.d0/2.d0 /)
    u2_coeff(1:2) = (/ 1.d0/2.d0, 1.d0/2.d0 /)
    
    !> left side
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = u1_coeff(1) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) + u1_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = u2_coeff(1) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)  + u2_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux)
  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,1) = (w(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux) - w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,2) = (w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,1:nwflux) - w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
  
    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,2
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) + weno_eps_machine)**2
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
      end do
    case(2)
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,2) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         1))
      do i = 1,2
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i) * (1.d0 + (tau(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) / (beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) + dxdim**2)))
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
      end do
    end select

    flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0
    do i = 1,2
      flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
    end do
  
    !> left value at right interface
    wLCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux)
  
    !> right side
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = u1_coeff(1) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) + u1_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = u2_coeff(1) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)  + u2_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,1) = (w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,1:nwflux) - w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,2) = (w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,1:nwflux) - w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
  
    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,2
       alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
          i) = d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
          i) + weno_eps_machine)**2
       alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
      end do
    case(2)
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,2) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         1))
      do i = 1,2
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i) * (1.d0 + (tau(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) / (beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) + dxdim**2)))
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
      end do
    end select

    flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0
    do i = 1,2
       flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
          iLmin2:iLmax2,1:nwflux) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
          1:nwflux,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
          i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
    end do
  
    !> right value at right interface
    wRCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux)

    call fix_limiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
       iLmax2,wLCtmp,wRCtmp,wLC,wRC)

  end subroutine WENO3limiter

  subroutine WENO5limiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
     iLmax2,idims,dxdim,w,wLC,wRC,var)
    use mod_global_parameters
  
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2, idims, var
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw),wLC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) 
    !> local
    integer                         :: iLmmin1,iLmmin2,iLmmax1,iLmmax2,&
        iLmmmin1,iLmmmin2,iLmmmax1,iLmmmax2, iLpmin1,iLpmin2,iLpmax1,iLpmax2,&
        iLppmin1,iLppmin2,iLppmax1,iLppmax2, iLpppmin1,iLpppmin2,iLpppmax1,&
       iLpppmax2
    double precision                :: f_array(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw,3), d_array(3)
    double precision                :: beta(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw,3), beta_coeff(2)
    double precision                :: tau(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), tmp(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw,3), alpha_sum(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), flux(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer                         :: i
    double precision                :: lambda

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)  :: wRCtmp, wLCtmp

    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
    iLmmmax1=iLmmax1-kr(idims,1);iLmmmax2=iLmmax2-kr(idims,2);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);
    iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
    iLpppmax1=iLppmax1+kr(idims,1);iLpppmax2=iLppmax2+kr(idims,2);
    lambda = dxdim**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
!   reconstruction variation
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
    u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
    u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
!   interpolation variation
!    d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
!    u1_coeff(1:3) = (/ 3.d0/8.d0, -10.d0/8.d0, 15.d0/8.d0 /)
!    u2_coeff(1:3) = (/ -1.d0/8.d0, 6.d0/8.d0, 3.d0/8.d0 /)
!    u3_coeff(1:3) = (/ 3.d0/8.d0, 6.d0/8.d0, -1.d0/8.d0 /)
    
    !> left side
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = u1_coeff(1) * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux) + u1_coeff(2) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) + u1_coeff(3) * w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = u2_coeff(1) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)  + u2_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)  + u2_coeff(3) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = u3_coeff(1) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)   + u3_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + u3_coeff(3) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux)  
  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = beta_coeff(1) * (w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux) + w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - 2.0d0*w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux))**2 + beta_coeff(2) * (w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux) - 4.0d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) + 3.0d0*w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = beta_coeff(1) * (w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux))**2 + beta_coeff(2) * (w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = beta_coeff(1) * (w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) - 2.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux))**2 + beta_coeff(2) * (3.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
        1:nwflux) - 4.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,1:nwflux))**2
 
    select case(var)
    ! case1 for wenojs, case2 for wenoz, case3 for wenoz+ 
    case(1)
      do i = 1,3
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) + weno_eps_machine)**2
      end do
    case(2)
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3))
      do i = 1,3
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i) * (1.d0 + (tau(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) / (beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) + weno_eps_machine))**2)
      end do
    case(3)
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3))
      do i = 1,3
        tmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = (tau(iLmin1:iLmax1,&
           iLmin2:iLmax2,1:nwflux) + weno_eps_machine) / (beta(iLmin1:iLmax1,&
           iLmin2:iLmax2,1:nwflux,i) + weno_eps_machine)
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i) * (1.0d0 + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux)**2 + lambda/tmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
      end do
    end select

    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0 
    do i = 1,3
       alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
    end do
    flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0
    do i = 1,3
       flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
          iLmin2:iLmax2,1:nwflux) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
          1:nwflux,i) *alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
          i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
    end do
    wLCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux)

    !> right side
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = u1_coeff(1) * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       1:nwflux) + u1_coeff(2) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) + u1_coeff(3) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = u2_coeff(1) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux)  + u2_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)  + u2_coeff(3) * w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = u3_coeff(1) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)   + u3_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)   + u3_coeff(3) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)  
  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = beta_coeff(1) * (w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       1:nwflux) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - 2.0d0*w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux))**2 + beta_coeff(2) * (w(iLpppmin1:iLpppmax1,&
       iLpppmin2:iLpppmax2,1:nwflux) - 4.0d0 * w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,1:nwflux) + 3.0d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = beta_coeff(1) * (w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) + w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - 2.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux))**2 + beta_coeff(2) * (w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) - w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = beta_coeff(1) * (w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux))**2 + beta_coeff(2) * (3.0d0 * w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2, 1:nwflux) - 4.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,1:nwflux))**2
  
    select case(var)
    case(1)
      do i = 1,3
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) + weno_eps_machine)**2
      end do
    case(2) 
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3))
      do i = 1,3
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i) * (1.d0 + (tau(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) / (beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) + weno_eps_machine))**2)
      end do
    case(3)
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3))
      do i = 1,3
        tmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = (tau(iLmin1:iLmax1,&
           iLmin2:iLmax2,1:nwflux) + weno_eps_machine) / (beta(iLmin1:iLmax1,&
           iLmin2:iLmax2,1:nwflux,i) + weno_eps_machine)
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i) * (1.0d0 + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux)**2 + lambda/tmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
      end do
    end select
  
    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0 
    ! do nothing, normal case
    do i = 1,3
       alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
    end do
    flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0
    do i = 1,3
       flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
          iLmin2:iLmax2,1:nwflux) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
          1:nwflux,i) *alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
          i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
    end do
    wRCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux)
    
    call fix_limiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
       iLmax2,wLCtmp,wRCtmp,wLC,wRC)

  end subroutine WENO5limiter

  subroutine TENO5ADlimiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
     iLmax1,iLmax2,idims,dxdim,w,wLC,wRC)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2, idims
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw),wLC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) 
    !> local
    integer                         :: iLmmin1,iLmmin2,iLmmax1,iLmmax2,&
        iLmmmin1,iLmmmin2,iLmmmax1,iLmmmax2, iLpmin1,iLpmin2,iLpmax1,iLpmax2,&
        iLppmin1,iLppmin2,iLppmax1,iLppmax2, iLpppmin1,iLpppmin2,iLpppmax1,&
       iLpppmax2
    double precision                :: f_array(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw,3), d_array(3)
    double precision                :: beta(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw,3), beta_coeff(2)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: gam_sum(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw),tau(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
       delta_sum(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision                :: gam(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw,3), kai(ixImin1:ixImax1,ixImin2:ixImax2,1:nw,3),&
        delta(ixImin1:ixImax1,ixImin2:ixImax2,1:nw,3)
    double precision                :: flux(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), kai1(ixImin1:ixImax1,ixImin2:ixImax2,1:nw,3),&
        theta(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer                         :: marray(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    integer                         :: i

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)  :: wRCtmp, wLCtmp

    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
    iLmmmax1=iLmmax1-kr(idims,1);iLmmmax2=iLmmax2-kr(idims,2);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);
    iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
    iLpppmax1=iLppmax1+kr(idims,1);iLpppmax2=iLppmax2+kr(idims,2);
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
    u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
    u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
    
    !> left side
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = u1_coeff(1) * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux) + u1_coeff(2) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) + u1_coeff(3) * w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = u2_coeff(1) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)  + u2_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)  + u2_coeff(3) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = u3_coeff(1) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)   + u3_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + u3_coeff(3) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux)  
  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = beta_coeff(1) * (w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux) + w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - 2.0d0*w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux))**2 + beta_coeff(2) * (w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux) - 4.0d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) + 3.0d0*w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = beta_coeff(1) * (w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux))**2 + beta_coeff(2) * (w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = beta_coeff(1) * (w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) - 2.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux))**2 + beta_coeff(2) * (3.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
        1:nwflux) - 4.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,1:nwflux))**2
 
    gam_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0
    tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3))
    do i = 1,3
      kai1(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i) = (tau(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux) / (beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i) + weno_eps_machine))
      gam(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i) = (1.d0 + kai1(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,i))**2
      gam_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = gam_sum(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux) + gam(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
    end do
    theta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = one / (one + &
       maxval(kai1(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,1:3)/10.d0,&
        dim=ndim+2))
    marray(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=-floor(4.d0 + &
       theta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)*6.d0)
    do i = 1,3    
      kai(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i) = gam(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,i) / gam_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux)
      where(kai(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i) .lt. 10**marray(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
        delta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)=zero
      else where
        delta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)=one
      end where
    end do
    delta_sum=zero
    do i = 1,3
      delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=delta_sum(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux)+delta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)*d_array(i)
    end do
    flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=0.0d0
    do i = 1,3
      flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=flux(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux)+f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)*delta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)*d_array(i)/(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
    end do
  
    !> left value at right interface
    wLCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux)
  
    !> right side
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = u1_coeff(1) * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       1:nwflux) + u1_coeff(2) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) + u1_coeff(3) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = u2_coeff(1) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux)  + u2_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)  + u2_coeff(3) * w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = u3_coeff(1) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)   + u3_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)   + u3_coeff(3) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)  
  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = beta_coeff(1) * (w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       1:nwflux) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - 2.0d0*w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux))**2 + beta_coeff(2) * (w(iLpppmin1:iLpppmax1,&
       iLpppmin2:iLpppmax2,1:nwflux) - 4.0d0 * w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,1:nwflux) + 3.0d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = beta_coeff(1) * (w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) + w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - 2.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux))**2 + beta_coeff(2) * (w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) - w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = beta_coeff(1) * (w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux))**2 + beta_coeff(2) * (3.0d0 * w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2, 1:nwflux) - 4.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,1:nwflux))**2
 

    gam_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=0.0d0
    tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=abs(beta(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux,1)-beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,3))
    do i=1,3
      kai1(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)=(tau(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)+weno_eps_machine))
      gam(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)=(1.d0+kai1(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,i))**6
      gam_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=gam_sum(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux)+gam(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
    end do
    theta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=one/(one+&
       maxval(kai1(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,1:3)/10.d0,&
       dim=ndim+2))
    marray(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=-floor(4.d0+&
       theta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)*6.d0)
    do i=1,3    
      kai(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i) = gam(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,i)/gam_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux)
      where(kai(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i) .lt. 10**marray(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
        delta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)=zero
      else where
        delta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)=one
      end where
    end do
    delta_sum=zero
    do i = 1,3
      delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=delta_sum(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux)+delta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)*d_array(i)
    end do
    flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=0.0d0
    do i = 1,3
      flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=flux(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux)+f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)*delta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)*d_array(i)/(delta_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
    end do

    !> right value at right interface
    wRCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux)

    call fix_limiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
       iLmax2,wLCtmp,wRCtmp,wLC,wRC)

  end subroutine TENO5ADlimiter

  subroutine WENO5NMlimiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
     iLmax1,iLmax2,idims,dxdim,w,wLC,wRC,var)
    use mod_global_parameters
  
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
       iLmax1,iLmax2,idims,var
    double precision, intent(in) :: dxdim
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw),wLC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) 
    !> local
    integer                         :: iLmmin1,iLmmin2,iLmmax1,iLmmax2,&
        iLmmmin1,iLmmmin2,iLmmmax1,iLmmmax2, iLpmin1,iLpmin2,iLpmax1,iLpmax2,&
        iLppmin1,iLppmin2,iLppmax1,iLppmax2, iLpppmin1,iLpppmin2,iLpppmax1,&
       iLpppmax2
    integer                         :: iMmin1,iMmin2,iMmax1,iMmax2, iMpmin1,&
       iMpmin2,iMpmax1,iMpmax2
    double precision                :: f_array(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw,3), d_array(3)
    double precision                :: beta(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw,3), beta_coeff(2)
    double precision                :: tau(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), tmp(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw,3), alpha_sum(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), flux(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision                :: wc(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), wd(ixImin1:ixImax1,ixImin2:ixImax2,1:nw), we(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw)
    integer                         :: i,j
    double precision                :: lambda(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)  :: wRCtmp, wLCtmp

    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
    iLmmmax1=iLmmax1-kr(idims,1);iLmmmax2=iLmmax2-kr(idims,2);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);
    iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
    iLpppmax1=iLppmax1+kr(idims,1);iLpppmax2=iLppmax2+kr(idims,2);
    iMmin1=iLmmin1;iMmin2=iLmmin2;
    iMmax1=iLpmax1;iMmax2=iLpmax2;
    iMpmin1=iMmin1+kr(idims,1);iMpmin2=iMmin2+kr(idims,2)
    iMpmax1=iMmax1+kr(idims,1);iMpmax2=iMmax2+kr(idims,2);
    lambda(iLmin1:iLmax1,iLmin2:iLmax2)=block%dx(iLmin1:iLmax1,iLmin2:iLmax2,&
       idims)**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    u1_coeff(1:3) = (/ -2.d0/3.d0, -1.d0/3.d0, 2.d0 /)
    u2_coeff(1:3) = (/ -1.d0/3.d0, 2.d0/3.d0, 2.d0/3.d0 /)
    u3_coeff(1:3) = (/ 2.d0/3.d0, 2.d0/3.d0, -1.d0/3.d0 /)
    do i = 1, nwflux
      wc(iMmin1:iMmax1,iMmin2:iMmax2,i) = (block%dx(iMpmin1:iMpmax1,&
         iMpmin2:iMpmax2,idims) * w(iMmin1:iMmax1,iMmin2:iMmax2,&
         i) + block%dx(iMmin1:iMmax1,iMmin2:iMmax2,idims) *  w(iMpmin1:iMpmax1,&
         iMpmin2:iMpmax2,i)) / (block%dx(iMpmin1:iMpmax1,iMpmin2:iMpmax2,&
         idims) + block%dx(iMmin1:iMmax1,iMmin2:iMmax2,idims))
      wd(iLmin1:iLmax1,iLmin2:iLmax2,i) = ((2.d0 * block%dx(iLmmin1:iLmmax1,&
         iLmmin2:iLmmax2,idims) + block%dx(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
         idims)) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         i) - block%dx(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         idims) * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
         i)) / (block%dx(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
         idims) + block%dx(iLmmin1:iLmmax1,iLmmin2:iLmmax2,idims))
      we(iLmin1:iLmax1,iLmin2:iLmax2,i) = ((2.d0 * block%dx(iLppmin1:iLppmax1,&
         iLppmin2:iLppmax2,idims) + block%dx(iLpppmin1:iLpppmax1,&
         iLpppmin2:iLpppmax2,idims)) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         i) - block%dx(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         idims) * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
         i)) / (block%dx(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
         idims) + block%dx(iLppmin1:iLppmax1,iLppmin2:iLppmax2,idims))
    enddo
    !> left side
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = u1_coeff(1) * wd(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)   + u1_coeff(2) * wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)+ u1_coeff(3) * w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = u2_coeff(1) * wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)  + u2_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)  + u2_coeff(3) * wc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = u3_coeff(1) * wc(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)   + u3_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + u3_coeff(3) * wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)  
  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = beta_coeff(1) * (wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - wd(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux))**2 + beta_coeff(2) * (2.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - wd(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = beta_coeff(1) * (wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) + wc(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux))**2 + beta_coeff(2) * (wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - wc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = beta_coeff(1) * (wc(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) + wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - 2.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux))**2 + beta_coeff(2) * (3.0d0 * wc(iLmin1:iLmax1,iLmin2:iLmax2,&
        1:nwflux) - 4.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux))**2
    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0 
    select case(var)
    ! case1 for wenojs, case2 for wenoz, case3 for wenoz+ 
    case(1)
      do i = 1,3
         alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
            i) = d_array(i)/(4.d0 * beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
            i) + weno_eps_machine)**2
         alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
            1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
            1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
      end do
    case(2)
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3)) * 4.d0
      do i = 1,3
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i) * (1.d0 + (tau(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) / (4.d0 * beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) + weno_eps_machine))**2)
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
      end do
    case(3)
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3)) * 4.d0
      do i=1,3
        do j=1,nwflux
          tmp(iLmin1:iLmax1,iLmin2:iLmax2,j) = (tau(iLmin1:iLmax1,&
             iLmin2:iLmax2,j) + weno_eps_machine) / (4.d0 * beta(iLmin1:iLmax1,&
             iLmin2:iLmax2,j,i) + weno_eps_machine)
          alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,j,&
             i) = d_array(i) * (1.0d0 + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
             j)**2 + lambda(iLmin1:iLmax1,iLmin2:iLmax2)/tmp(iLmin1:iLmax1,&
             iLmin2:iLmax2,j))
          alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,j) = alpha_sum(iLmin1:iLmax1,&
             iLmin2:iLmax2,j) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,j,i)
        end do
      end do
    end select
    flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0
    do i = 1,3
      flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
    end do
    !> left value at right interface
    wLCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux)

    !> right side
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = u1_coeff(1) * we(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)  + u1_coeff(2) * wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + u1_coeff(3) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = u2_coeff(1) * wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + u2_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + u2_coeff(3) * wc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = u3_coeff(1) * wc(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)  + u3_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)  + u3_coeff(3) * wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = beta_coeff(1) * (wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - we(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux))**2 + beta_coeff(2) * (2.d0 * w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,1:nwflux) - wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - we(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = beta_coeff(1) * (wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + wc(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - 2.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux))**2 + beta_coeff(2) * (wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - wc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = beta_coeff(1) * (wc(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) + wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux))**2 + beta_coeff(2) * (3.0d0 * wc(iLmin1:iLmax1,iLmin2:iLmax2,&
        1:nwflux) - 4.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) + wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,1:nwflux))**2
    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,3
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i)/(4.d0 * beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) + weno_eps_machine)**2
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
      end do
    case(2) 
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3)) * 4.d0
      do i = 1,3
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i) * (1.d0 + (tau(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) / (4.d0 * beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) + weno_eps_machine))**2)
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
      end do
    case(3)
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3)) * 4.d0
      do i = 1,3
        do j = 1,nwflux
          tmp(iLmin1:iLmax1,iLmin2:iLmax2,j) = (tau(iLmin1:iLmax1,&
             iLmin2:iLmax2,j) + weno_eps_machine) / (4.d0 * beta(iLmin1:iLmax1,&
             iLmin2:iLmax2,j,i) + weno_eps_machine)
          alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,j,&
             i) = d_array(i) * (1.0d0 + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
             j)**2 + lambda(iLmin1:iLmax1,iLmin2:iLmax2)/tmp(iLmin1:iLmax1,&
             iLmin2:iLmax2,j))
          alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,j) = alpha_sum(iLmin1:iLmax1,&
             iLmin2:iLmax2,j) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,j,i)
        end do
      end do
    end select
    flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0
    do i = 1,3
      flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
    end do
    !> right value at right interface
    wRCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux)

    call fix_limiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
       iLmax2,wLCtmp,wRCtmp,wLC,wRC)

  end subroutine WENO5NMlimiter

  subroutine WENO5limiterL(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
     iLmax1,iLmax2,idims,w,wLC,var)
    use mod_global_parameters
  
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2, idims
    integer, intent(in)             :: var
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw) 
    !> local
    integer                         :: iLmmin1,iLmmin2,iLmmax1,iLmmax2,&
        iLmmmin1,iLmmmin2,iLmmmax1,iLmmmax2, iLpmin1,iLpmin2,iLpmax1,iLpmax2,&
        iLppmin1,iLppmin2,iLppmax1,iLppmax2, iLpppmin1,iLpppmin2,iLpppmax1,&
       iLpppmax2
    double precision                :: f_array(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw,3), d_array(3)
    double precision                :: beta(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw,3), beta_coeff(2)
    double precision                :: tau(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), tmp(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw,3), alpha_sum(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), flux(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer                         :: i,j
    double precision                :: lambda(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)  :: wLCtmp

    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
    iLmmmax1=iLmmax1-kr(idims,1);iLmmmax2=iLmmax2-kr(idims,2);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);
    lambda=block%dx(iLmin1:iLmax1,iLmin2:iLmax2,idims)**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
!   reconstruction variation
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
    u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
    u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
!   interpolation variation
!    d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
!    u1_coeff(1:3) = (/ 3.d0/8.d0, -10.d0/8.d0, 15.d0/8.d0 /)
!    u2_coeff(1:3) = (/ -1.d0/8.d0, 6.d0/8.d0, 3.d0/8.d0 /)
!    u3_coeff(1:3) = (/ 3.d0/8.d0, 6.d0/8.d0, -1.d0/8.d0 /)
    
    !> left side
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = u1_coeff(1) * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux) + u1_coeff(2) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) + u1_coeff(3) * w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = u2_coeff(1) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)  + u2_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)  + u2_coeff(3) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = u3_coeff(1) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)   + u3_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + u3_coeff(3) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux)  
  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = beta_coeff(1) * (w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux) + w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - 2.0d0*w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux))**2 + beta_coeff(2) * (w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux) - 4.0d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) + 3.0d0*w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = beta_coeff(1) * (w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux))**2 + beta_coeff(2) * (w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = beta_coeff(1) * (w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) - 2.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux))**2 + beta_coeff(2) * (3.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
        1:nwflux) - 4.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,1:nwflux))**2
 
    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0 
    select case(var)
    ! case1 for wenojs, case2 for wenoz
    case(1)
      do i = 1,3
         alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
            i) = d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
            i) + weno_eps_machine)**2
         alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
            1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
            1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
      end do
    case(2)
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3))
      do i = 1,3
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i) * (1.d0 + (tau(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) / (beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) + weno_eps_machine))**2)
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
      end do
    case(3)
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3)) * 4.d0
      do i=1,3
        do j=1,nwflux
          tmp(iLmin1:iLmax1,iLmin2:iLmax2,j) = (tau(iLmin1:iLmax1,&
             iLmin2:iLmax2,j) + weno_eps_machine) / (4.d0 * beta(iLmin1:iLmax1,&
             iLmin2:iLmax2,j,i) + weno_eps_machine)
          alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,j,&
             i) = d_array(i) * (1.0d0 + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
             j)**2 + lambda(iLmin1:iLmax1,iLmin2:iLmax2)/tmp(iLmin1:iLmax1,&
             iLmin2:iLmax2,j))
          alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,j) = alpha_sum(iLmin1:iLmax1,&
             iLmin2:iLmax2,j) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,j,i)
        end do
      end do
    end select
    flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0
    do i = 1,3
      flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
    end do
  
    !> left value at right interface
    wLCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux)
  
    call fix_onelimiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
       iLmax2,wLCtmp,wLC)

  end subroutine WENO5limiterL

  subroutine WENO5limiterR(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
     iLmax1,iLmax2,idims,w,wRC,var)
    use mod_global_parameters
  
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2, idims
    integer, intent(in)             :: var
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    !> local
    integer                         :: iLmmin1,iLmmin2,iLmmax1,iLmmax2,&
        iLmmmin1,iLmmmin2,iLmmmax1,iLmmmax2, iLpmin1,iLpmin2,iLpmax1,iLpmax2,&
        iLppmin1,iLppmin2,iLppmax1,iLppmax2, iLpppmin1,iLpppmin2,iLpppmax1,&
       iLpppmax2
    double precision                :: f_array(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw,3), d_array(3)
    double precision                :: beta(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw,3), beta_coeff(2)
    double precision                :: tau(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), tmp(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw,3), alpha_sum(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), flux(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer                         :: i,j
    double precision                :: lambda(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)  :: wRCtmp

    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);
    iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
    iLpppmax1=iLppmax1+kr(idims,1);iLpppmax2=iLppmax2+kr(idims,2);
    lambda=block%dx(iLmin1:iLmax1,iLmin2:iLmax2,idims)**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
!   reconstruction variation
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
    u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
    u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
!   interpolation variation
!    d_array(1:3) = (/ 1.0d0/16.0d0, 10.0d0/16.0d0, 5.0d0/16.0d0 /)
!    u1_coeff(1:3) = (/ 3.d0/8.d0, -10.d0/8.d0, 15.d0/8.d0 /)
!    u2_coeff(1:3) = (/ -1.d0/8.d0, 6.d0/8.d0, 3.d0/8.d0 /)
!    u3_coeff(1:3) = (/ 3.d0/8.d0, 6.d0/8.d0, -1.d0/8.d0 /)
    
    !> right side
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = u1_coeff(1) * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       1:nwflux) + u1_coeff(2) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) + u1_coeff(3) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = u2_coeff(1) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux)  + u2_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)  + u2_coeff(3) * w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = u3_coeff(1) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)   + u3_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)   + u3_coeff(3) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)  
  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = beta_coeff(1) * (w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       1:nwflux) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - 2.0d0*w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux))**2 + beta_coeff(2) * (w(iLpppmin1:iLpppmax1,&
       iLpppmin2:iLpppmax2,1:nwflux) - 4.0d0 * w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,1:nwflux) + 3.0d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = beta_coeff(1) * (w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) + w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - 2.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux))**2 + beta_coeff(2) * (w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) - w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = beta_coeff(1) * (w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux))**2 + beta_coeff(2) * (3.0d0 * w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2, 1:nwflux) - 4.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,1:nwflux))**2
  
    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,3
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) + weno_eps_machine)**2
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
      end do
    case(2) 
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3))
      do i = 1,3
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i) * (1.d0 + (tau(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) / (beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) + weno_eps_machine))**2)
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
      end do
    case(3)
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3)) * 4.d0
      do i = 1,3
        do j = 1,nwflux
          tmp(iLmin1:iLmax1,iLmin2:iLmax2,j) = (tau(iLmin1:iLmax1,&
             iLmin2:iLmax2,j) + weno_eps_machine) / (4.d0 * beta(iLmin1:iLmax1,&
             iLmin2:iLmax2,j,i) + weno_eps_machine)
          alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,j,&
             i) = d_array(i) * (1.0d0 + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
             j)**2 + lambda(iLmin1:iLmax1,iLmin2:iLmax2)/tmp(iLmin1:iLmax1,&
             iLmin2:iLmax2,j))
          alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,j) = alpha_sum(iLmin1:iLmax1,&
             iLmin2:iLmax2,j) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,j,i)
        end do
      end do
    end select
    flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0
    do i = 1,3
      flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
    end do
  
    !> right value at right interface
    wRCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux)

    call fix_onelimiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
       iLmax2,wRCtmp,wRC)

  end subroutine WENO5limiterR

  subroutine WENO5NMlimiterL(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
     iLmax1,iLmax2,idims,w,wLC,var)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
       iLmax1,iLmax2,idims,var
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw) 
    !> local
    integer                         :: iLmmin1,iLmmin2,iLmmax1,iLmmax2,&
        iLmmmin1,iLmmmin2,iLmmmax1,iLmmmax2, iLpmin1,iLpmin2,iLpmax1,iLpmax2,&
        iLppmin1,iLppmin2,iLppmax1,iLppmax2
    integer                         :: iMmin1,iMmin2,iMmax1,iMmax2, iMpmin1,&
       iMpmin2,iMpmax1,iMpmax2
    double precision                :: f_array(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw,3), d_array(3)
    double precision                :: beta(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw,3), beta_coeff(2)
    double precision                :: tau(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), tmp(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw,3), alpha_sum(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), flux(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision                :: wc(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), wd(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer                         :: i,j
    double precision                :: lambda(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)  :: wLCtmp

    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
    iLmmmax1=iLmmax1-kr(idims,1);iLmmmax2=iLmmax2-kr(idims,2);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);
    iMmin1=iLmmin1;iMmin2=iLmmin2;
    iMmax1=iLpmax1;iMmax2=iLpmax2;
    iMpmin1=iMmin1+kr(idims,1);iMpmin2=iMmin2+kr(idims,2)
    iMpmax1=iMmax1+kr(idims,1);iMpmax2=iMmax2+kr(idims,2);
    lambda(iLmin1:iLmax1,iLmin2:iLmax2)=block%dx(iLmin1:iLmax1,iLmin2:iLmax2,&
       idims)**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    u1_coeff(1:3) = (/ -2.d0/3.d0, -1.d0/3.d0, 2.d0 /)
    u2_coeff(1:3) = (/ -1.d0/3.d0, 2.d0/3.d0, 2.d0/3.d0 /)
    u3_coeff(1:3) = (/ 2.d0/3.d0, 2.d0/3.d0, -1.d0/3.d0 /)
    do i = 1, nwflux
      wc(iMmin1:iMmax1,iMmin2:iMmax2,i) = (block%dx(iMpmin1:iMpmax1,&
         iMpmin2:iMpmax2,idims) * w(iMmin1:iMmax1,iMmin2:iMmax2,&
         i) + block%dx(iMmin1:iMmax1,iMmin2:iMmax2,idims) * w(iMpmin1:iMpmax1,&
         iMpmin2:iMpmax2,i)) / (block%dx(iMpmin1:iMpmax1,iMpmin2:iMpmax2,&
         idims) + block%dx(iMmin1:iMmax1,iMmin2:iMmax2,idims))
      wd(iLmin1:iLmax1,iLmin2:iLmax2,i) = ((2.d0 * block%dx(iLmmin1:iLmmax1,&
         iLmmin2:iLmmax2,idims) + block%dx(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
         idims)) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         i) - block%dx(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         idims) * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
         i)) / (block%dx(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
         idims) + block%dx(iLmmin1:iLmmax1,iLmmin2:iLmmax2,idims))
    enddo
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = u1_coeff(1) * wd(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)   + u1_coeff(2) * wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)+ u1_coeff(3) * w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = u2_coeff(1) * wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)  + u2_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)  + u2_coeff(3) * wc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = u3_coeff(1) * wc(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)   + u3_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + u3_coeff(3) * wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = beta_coeff(1) * (wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - wd(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux))**2 + beta_coeff(2) * (2.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - wd(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = beta_coeff(1) * (wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) + wc(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux))**2 + beta_coeff(2) * (wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - wc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = beta_coeff(1) * (wc(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) + wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - 2.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux))**2 + beta_coeff(2) * (3.0d0 * wc(iLmin1:iLmax1,iLmin2:iLmax2,&
        1:nwflux) - 4.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux))**2
    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0 
    select case(var)
    ! case1 for wenojs, case2 for wenoz, case3 for wenoz+ 
    case(1)
      do i = 1,3
         alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
            i) = d_array(i)/(4.d0 * beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
            i) + weno_eps_machine)**2
         alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
            1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
            1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
      end do
    case(2)
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3)) * 4.d0
      do i = 1,3
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i) * (1.d0 + (tau(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) / (4.d0 * beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) + weno_eps_machine))**2)
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
      end do
    case(3)
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3)) * 4.d0
      do i=1,3
        do j=1,nwflux
          tmp(iLmin1:iLmax1,iLmin2:iLmax2,j) = (tau(iLmin1:iLmax1,&
             iLmin2:iLmax2,j) + weno_eps_machine) / (4.d0 * beta(iLmin1:iLmax1,&
             iLmin2:iLmax2,j,i) + weno_eps_machine)
          alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,j,&
             i) = d_array(i) * (1.0d0 + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
             j)**2 + lambda(iLmin1:iLmax1,iLmin2:iLmax2)/tmp(iLmin1:iLmax1,&
             iLmin2:iLmax2,j))
          alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,j) = alpha_sum(iLmin1:iLmax1,&
             iLmin2:iLmax2,j) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,j,i)
        end do
      end do
    end select
    flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0
    do i = 1,3
      flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
    end do
    !> left value at right interface
    wLCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux)

    call fix_onelimiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
       iLmax2,wLCtmp,wLC)

  end subroutine WENO5NMlimiterL

  subroutine WENO5NMlimiterR(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
     iLmax1,iLmax2,idims,w,wRC,var)
    use mod_global_parameters
  
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
       iLmax1,iLmax2,idims,var
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    !> local
    integer                         :: iLmmin1,iLmmin2,iLmmax1,iLmmax2,iLpmin1,&
       iLpmin2,iLpmax1,iLpmax2,iLppmin1,iLppmin2,iLppmax1,iLppmax2,iLpppmin1,&
       iLpppmin2,iLpppmax1,iLpppmax2
    integer                         :: iMmin1,iMmin2,iMmax1,iMmax2, iMpmin1,&
       iMpmin2,iMpmax1,iMpmax2
    double precision                :: f_array(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw,3), d_array(3)
    double precision                :: beta(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw,3), beta_coeff(2)
    double precision                :: tau(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), tmp(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision                :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision                :: alpha_array(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw,3), alpha_sum(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), flux(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision                :: wc(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), we(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer                         :: i,j
    double precision                :: lambda(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)  :: wRCtmp

    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);
    iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
    iLpppmax1=iLppmax1+kr(idims,1);iLpppmax2=iLppmax2+kr(idims,2);
    iMmin1=iLmmin1;iMmin2=iLmmin2;
    iMmax1=iLpmax1;iMmax2=iLpmax2;
    iMpmin1=iMmin1+kr(idims,1);iMpmin2=iMmin2+kr(idims,2)
    iMpmax1=iMmax1+kr(idims,1);iMpmax2=iMmax2+kr(idims,2);
    lambda(iLmin1:iLmax1,iLmin2:iLmax2)=block%dx(iLmin1:iLmax1,iLmin2:iLmax2,&
       idims)**weno_dx_exp
    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    u1_coeff(1:3) = (/ -2.d0/3.d0, -1.d0/3.d0, 2.d0 /)
    u2_coeff(1:3) = (/ -1.d0/3.d0, 2.d0/3.d0, 2.d0/3.d0 /)
    u3_coeff(1:3) = (/ 2.d0/3.d0, 2.d0/3.d0, -1.d0/3.d0 /)
    do i = 1, nwflux
      wc(iMmin1:iMmax1,iMmin2:iMmax2,i) = (block%dx(iMpmin1:iMpmax1,&
         iMpmin2:iMpmax2,idims) * w(iMmin1:iMmax1,iMmin2:iMmax2,&
         i) + block%dx(iMmin1:iMmax1,iMmin2:iMmax2,idims) * w(iMpmin1:iMpmax1,&
         iMpmin2:iMpmax2,i)) / (block%dx(iMpmin1:iMpmax1,iMpmin2:iMpmax2,&
         idims) + block%dx(iMmin1:iMmax1,iMmin2:iMmax2,idims))
      we(iLmin1:iLmax1,iLmin2:iLmax2,i) = ((2.d0 * block%dx(iLppmin1:iLppmax1,&
         iLppmin2:iLppmax2,idims) + block%dx(iLpppmin1:iLpppmax1,&
         iLpppmin2:iLpppmax2,idims)) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         i) - block%dx(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         idims) * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
         i)) / (block%dx(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
         idims) + block%dx(iLppmin1:iLppmax1,iLppmin2:iLppmax2,idims))
    enddo
    !> right side
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = u1_coeff(1) * we(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)  + u1_coeff(2) * wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + u1_coeff(3) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = u2_coeff(1) * wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + u2_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + u2_coeff(3) * wc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = u3_coeff(1) * wc(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)  + u3_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)  + u3_coeff(3) * wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = beta_coeff(1) * (wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - we(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux))**2 + beta_coeff(2) * (2.d0 * w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,1:nwflux) - wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - we(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = beta_coeff(1) * (wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + wc(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - 2.0d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux))**2 + beta_coeff(2) * (wc(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - wc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = beta_coeff(1) * (wc(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) + wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - 2.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux))**2 + beta_coeff(2) * (3.0d0 * wc(iLmin1:iLmax1,iLmin2:iLmax2,&
        1:nwflux) - 4.0d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) + wc(iLmmin1:iLmmax1,iLmmin2:iLmmax2,1:nwflux))**2
    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0 
    select case(var)
    case(1)
      do i = 1,3
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i)/(4.d0 * beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) + weno_eps_machine)**2
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
      end do
    case(2) 
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3)) * 4.d0
      do i = 1,3
        alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) = d_array(i) * (1.d0 + (tau(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) / (4.d0 * beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
           i) + weno_eps_machine))**2)
        alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
           1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
      end do
    case(3)
      tau(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = abs(beta(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1) - beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3)) * 4.d0
      do i = 1,3
        do j = 1,nwflux
          tmp(iLmin1:iLmax1,iLmin2:iLmax2,j) = (tau(iLmin1:iLmax1,&
             iLmin2:iLmax2,j) + weno_eps_machine) / (4.d0 * beta(iLmin1:iLmax1,&
             iLmin2:iLmax2,j,i) + weno_eps_machine)
          alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,j,&
             i) = d_array(i) * (1.0d0 + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
             j)**2 + lambda(iLmin1:iLmax1,iLmin2:iLmax2)/tmp(iLmin1:iLmax1,&
             iLmin2:iLmax2,j))
          alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,j) = alpha_sum(iLmin1:iLmax1,&
             iLmin2:iLmax2,j) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,j,i)
        end do
      end do
    end select
    flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0
    do i = 1,3
      flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
    end do
    !> right value at right interface
    wRCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux)

    call fix_onelimiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
       iLmax2,wRCtmp,wRC)

  end subroutine WENO5NMlimiterR

  subroutine WENO5CU6limiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,&
     iLmax1,iLmax2,idims,w,wLC,wRC)
    use mod_global_parameters
  
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,iLmin2,&
       iLmax1,iLmax2, idims
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw),wLC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) 
    !> local
    integer :: iLmmin1,iLmmin2,iLmmax1,iLmmax2, iLmmmin1,iLmmmin2,iLmmmax1,&
       iLmmmax2, iLpmin1,iLpmin2,iLpmax1,iLpmax2, iLppmin1,iLppmin2,iLppmax1,&
       iLppmax2, iLpppmin1,iLpppmin2,iLpppmax1,iLpppmax2
    double precision :: f_array(ixImin1:ixImax1,ixImin2:ixImax2,1:nw,3),&
        d_array(3)
    double precision :: beta(ixImin1:ixImax1,ixImin2:ixImax2,1:nw,3),&
        beta_coeff(2)
    double precision :: u1_coeff(3), u2_coeff(3), u3_coeff(3)
    double precision :: alpha_array(ixImin1:ixImax1,ixImin2:ixImax2,1:nw,3),&
        alpha_sum(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision :: theta2(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer :: i
    double precision, parameter :: theta_limit=0.7d0

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)  :: wRCtmp, wLCtmp

    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
    iLmmmax1=iLmmax1-kr(idims,1);iLmmmax2=iLmmax2-kr(idims,2);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);
    iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
    iLpppmax1=iLppmax1+kr(idims,1);iLpppmax2=iLppmax2+kr(idims,2);

    beta_coeff(1:2) = (/ 1.0833333333333333d0, 0.25d0/)
    d_array(1:3) = (/ 1.0d0/10.0d0, 3.0d0/5.0d0, 3.0d0/10.0d0 /)
    u1_coeff(1:3) = (/ 1.d0/3.d0, -7.d0/6.d0, 11.d0/6.d0 /)
    u2_coeff(1:3) = (/ -1.d0/6.d0, 5.d0/6.d0, 1.d0/3.d0 /)
    u3_coeff(1:3) = (/ 1.d0/3.d0, 5.d0/6.d0, -1.d0/6.d0 /)
    
    !> left side
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1)=beta_coeff(1)*(w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux)+w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)-2.0d0*w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux))**2+beta_coeff(2)*(w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux)-4.0d0*w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)+3.0d0*w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2)=beta_coeff(1)*(w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)-2.0d0*w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux))**2+beta_coeff(2)*(w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)-w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3)=beta_coeff(1)*(w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)+w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux)-2.0d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux))**2+beta_coeff(2)*(3.0d0*w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)-4.0d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)+w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,1:nwflux))**2
    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=zero
    do i=1,3
      alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)=d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)+weno_eps_machine)**2
      alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=alpha_sum(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux)+alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux,i)
    end do
    do i=1,3
      alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)=alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)/alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
    end do
    theta2(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=((alpha_array(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux,1)/d_array(1)-1.d0)**2+&
       (alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2)/d_array(2)-1.d0)**2+(alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux,3)/d_array(3)-1.d0)**2)/83.d0
    where(theta2(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) .gt. theta_limit)
      f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         1)=u1_coeff(1)*w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
         1:nwflux)+u1_coeff(2)*w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         1:nwflux)+u1_coeff(3)*w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
      f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         2)=u2_coeff(1)*w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         1:nwflux)+u2_coeff(2)*w(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux)+u2_coeff(3)*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux)
      f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3)=u3_coeff(1)*w(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux)+u3_coeff(2)*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         1:nwflux)+u3_coeff(3)*w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         1:nwflux)  
      wLCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=f_array(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1)*alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux,1)+f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         2)*alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         2)+f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3)*alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,3)
    else where
      wLCtmp(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux)=1.d0/60.d0*(w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
         1:nwflux)-8.d0*w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         1:nwflux)+37.d0*w(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux)+37.d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         1:nwflux)-8.d0*w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         1:nwflux)+w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,1:nwflux))
    end where
  
    !> right side
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1)=beta_coeff(1)*(w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       1:nwflux)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)-2.0d0*w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux))**2+beta_coeff(2)*(w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       1:nwflux)-4.0d0*w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux)+3.0d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2)=beta_coeff(1)*(w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux)+w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)-2.0d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux))**2+beta_coeff(2)*(w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux)-w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3)=beta_coeff(1)*(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)+w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)-2.0d0*w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux))**2+beta_coeff(2)*(3.0d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)-4.0d0*w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)+w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,1:nwflux))**2
    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=zero
    do i=1,3
      alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)=d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)+weno_eps_machine)**2
      alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=alpha_sum(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux)+alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux,i)
    end do
    do i=1,3
      alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)=alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         i)/alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
    end do
    theta2(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=((alpha_array(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux,1)/d_array(1)-1.d0)**2+&
       (alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2)/d_array(2)-1.d0)**2+(alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux,3)/d_array(3)-1.d0)**2)/83.d0
    where(theta2(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) .gt. theta_limit)
      f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         1)=u1_coeff(1)*w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
         1:nwflux)+u1_coeff(2)*w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         1:nwflux)+u1_coeff(3)*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux)
      f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         2)=u2_coeff(1)*w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         1:nwflux)+u2_coeff(2)*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         1:nwflux)+u2_coeff(3)*w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
      f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3)=u3_coeff(1)*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         1:nwflux)+u3_coeff(2)*w(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux)+u3_coeff(3)*w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,1:nwflux)  
      wRCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)=f_array(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux,1)*alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux,1)+f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         2)*alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         2)+f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
         3)*alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,3)
    else where
      wRCtmp(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux)=1.d0/60.d0*(w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
         1:nwflux)-8.d0*w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         1:nwflux)+37.d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         1:nwflux)+37.d0*w(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux)-8.d0*w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         1:nwflux)+w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,1:nwflux))
    end where

    call fix_limiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
       iLmax2,wLCtmp,wRCtmp,wLC,wRC)

  end subroutine WENO5CU6limiter

  subroutine WENO7limiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
     iLmax2,idims,w,wLC,wRC,var)
    use mod_global_parameters
  
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,&
       iLmin2,iLmax1,iLmax2, idims, var
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw),wLC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) 
    !> local
    integer                         :: iLmmin1,iLmmin2,iLmmax1,iLmmax2,&
        iLmmmin1,iLmmmin2,iLmmmax1,iLmmmax2, iLmmmmin1,iLmmmmin2,iLmmmmax1,&
       iLmmmmax2
    integer                         :: iLpmin1,iLpmin2,iLpmax1,iLpmax2,&
        iLppmin1,iLppmin2,iLppmax1,iLppmax2, iLpppmin1,iLpppmin2,iLpppmax1,&
       iLpppmax2, iLppppmin1,iLppppmin2,iLppppmax1,iLppppmax2
    integer                         :: idmin1,idmin2,idmax1,idmax2, idpmin1,&
       idpmin2,idpmax1,idpmax2, idppmin1,idppmin2,idppmax1,idppmax2, idmmin1,&
       idmmin2,idmmax1,idmmax2, iemin1,iemin2,iemax1,iemax2, iemmin1,iemmin2,&
       iemmax1,iemmax2, iepmin1,iepmin2,iepmax1,iepmax2, ieppmin1,ieppmin2,&
       ieppmax1,ieppmax2

    double precision, dimension(4)  :: d_array, u1_coeff, u2_coeff, u3_coeff,&
        u4_coeff
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw,&
       4)  :: f_array, beta, alpha_array
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)         :: a,&
        b, c, tmp, tmp2, tmp3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)    :: alpha_sum, d, dm4
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)    :: flux, flux_min, flux_max, flux_ul, flux_md, flux_lc
    integer                         :: i,iw
    double precision, parameter     :: mpalpha = 2.d0, mpbeta = 4.d0

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)  :: wRCtmp, wLCtmp

    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmax1=iLmax1-kr(idims,1);iLmmax2=iLmax2-kr(idims,2);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
    iLmmmax1=iLmmax1-kr(idims,1);iLmmmax2=iLmmax2-kr(idims,2);
    iLmmmmin1=iLmmmin1-kr(idims,1);iLmmmmin2=iLmmmin2-kr(idims,2)
    iLmmmmax1=iLmmmax1-kr(idims,1);iLmmmmax2=iLmmmax2-kr(idims,2);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmax1=iLmax1+kr(idims,1);iLpmax2=iLmax2+kr(idims,2);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmax1=iLpmax1+kr(idims,1);iLppmax2=iLpmax2+kr(idims,2);
    iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
    iLpppmax1=iLppmax1+kr(idims,1);iLpppmax2=iLppmax2+kr(idims,2);
    iLppppmin1=iLpppmin1+kr(idims,1);iLppppmin2=iLpppmin2+kr(idims,2)
    iLppppmax1=iLpppmax1+kr(idims,1);iLppppmax2=iLpppmax2+kr(idims,2);

    d_array(1:4) = (/ 1.d0/35.d0, 12.d0/35.d0, 18.d0/35.d0, 4.d0/35.d0 /)
    u1_coeff(1:4) = (/ -1.d0/4.d0, 13.d0/12.d0, -23.d0/12.d0, 25.d0/12.d0 /)
    u2_coeff(1:4) = (/ 1.d0/12.d0, -5.d0/12.d0, 13.d0/12.d0, 1.d0/4.d0 /)
    u3_coeff(1:4) = (/ -1.d0/12.d0, 7.d0/12.d0, 7.d0/12.d0, -1.d0/12.d0 /)
    u4_coeff(1:4) = (/ 1.d0/4.d0, 13.d0/12.d0, -5.d0/12.d0, 1.d0/12.d0 /)
    
    !> left side
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = u1_coeff(1) * w(iLmmmmin1:iLmmmmax1,iLmmmmin2:iLmmmmax2,&
       1:nwflux) + u1_coeff(2) * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux) + u1_coeff(3) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) + u1_coeff(4) * w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = u2_coeff(1) * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux)  + u2_coeff(2) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)  + u2_coeff(3) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)  + u2_coeff(4) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = u3_coeff(1) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)   + u3_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)   + u3_coeff(3) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)   + u3_coeff(4) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       4) = u4_coeff(1) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)    + u4_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)    + u4_coeff(3) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux)  + u4_coeff(4) * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       1:nwflux)
  
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,1) = w(iLmmmmin1:iLmmmmax1,&
       iLmmmmin2:iLmmmmax2,1:nwflux) * (547.d0 * w(iLmmmmin1:iLmmmmax1,&
       iLmmmmin2:iLmmmmax2,1:nwflux) - 3882.d0 * w(iLmmmin1:iLmmmax1,&
       iLmmmin2:iLmmmax2,1:nwflux) + 4642.d0 * w(iLmmin1:iLmmax1,&
       iLmmin2:iLmmax2,1:nwflux) - 1854.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)) + w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux) * (7043.d0 * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux) - 17246.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) + 7042.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)) + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) * (11003.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - 9402.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)) + 2107.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,2) = w(iLmmmin1:iLmmmax1,&
       iLmmmin2:iLmmmax2,1:nwflux) * (267.d0 * w(iLmmmin1:iLmmmax1,&
       iLmmmin2:iLmmmax2,1:nwflux) - 1642.d0 * w(iLmmin1:iLmmax1,&
       iLmmin2:iLmmax2,1:nwflux) + 1602.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - 494.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux))  + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) * (2843.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - 5966.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) + 1922.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)) + w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) * (3443.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - 2522.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)) + 547.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux) ** 2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,3) = w(iLmmin1:iLmmax1,&
       iLmmin2:iLmmax2,1:nwflux) * (547.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - 2522.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) + 1922.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - 494.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux))  + w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) * (3443.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - 5966.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + 1602.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux)) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) * (2843.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - 1642.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux)) + 267.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) ** 2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,4) = w(iLmin1:iLmax1,&
       iLmin2:iLmax2,1:nwflux) * (2107.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - 9402.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + 7042.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) - 1854.d0 * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       1:nwflux))  + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) * (11003.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - 17246.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) + 4642.d0 * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       1:nwflux)) + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) * (7043.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) - 3882.d0 * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       1:nwflux)) + 547.d0 * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       1:nwflux) ** 2

    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0 
    do i = 1,4
       alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
          i) = d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
          i) + weno_eps_machine)**2
       alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
    end do
    flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0
    do i = 1,4
       flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
          iLmin2:iLmax2,1:nwflux) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
          1:nwflux,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
          i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
    end do
    
    select case(var)
    ! case1 for wenojs, case2 for mpweno
    case(1) 
      wLCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux)
    case(2)
      idmax1=iLmax1;idmax2=iLmax2; idmin1=iLmin1-kr(idims,1)
      idmin2=iLmin2-kr(idims,2);
      idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
      idmmax1=idmax1-kr(idims,1);idmmax2=idmax2-kr(idims,2);
      idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
      idpmax1=idmax1+kr(idims,1);idpmax2=idmax2+kr(idims,2);
  
      iemax1=idmax1+kr(idims,1);iemax2=idmax2+kr(idims,2); iemin1=idmin1
      iemin2=idmin2;
      iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
      iemmax1=iemax1-kr(idims,1);iemmax2=iemax2-kr(idims,2);
      iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
      iepmax1=iemax1+kr(idims,1);iepmax2=iemax2+kr(idims,2);
  
      d(iemin1:iemax1,iemin2:iemax2,1:nwflux) = w(iepmin1:iepmax1,&
         iepmin2:iepmax2,1:nwflux)-2.0d0*w(iemin1:iemax1,iemin2:iemax2,&
         1:nwflux)+w(iemmin1:iemmax1,iemmin2:iemmax2,1:nwflux)
  
      do iw=1,nwflux
         a(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idmin1:idmax1,idmin2:idmax2,&
            iw)-d(idpmin1:idpmax1,idpmin2:idpmax2,iw)
         b(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idpmin1:idpmax1,&
            idpmin2:idpmax2,iw)-d(idmin1:idmax1,idmin2:idmax2,iw)
         call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,&
            idmax2,a,b,tmp)
         a(idmin1:idmax1,idmin2:idmax2) = d(idmin1:idmax1,idmin2:idmax2,iw)
         b(idmin1:idmax1,idmin2:idmax2) = d(idpmin1:idpmax1,idpmin2:idpmax2,&
            iw)
         call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,&
            idmax2,a,b,tmp2)
         call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,&
            idmax2,tmp,tmp2,tmp3)
         dm4(idmin1:idmax1,idmin2:idmax2,iw) = tmp3(idmin1:idmax1,&
            idmin2:idmax2)
      end do

      flux_ul(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = w(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux) + mpalpha * (w(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux) - w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,1:nwflux))
      flux_md(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = half * (w(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         1:nwflux) - dm4(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
      flux_lc(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux) = half * (3.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux) - w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         1:nwflux)) + mpbeta / 3.d0 * dm4(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
         1:nwflux)
    
      flux_min(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = max(min(w(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux), w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux),&
          flux_md(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)), min(w(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux), flux_ul(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux),flux_lc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)))
  
      flux_max(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = min(max(w(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux), w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux),&
          flux_md(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)), max(w(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux), flux_ul(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux),flux_lc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)))
      do iw=1,nwflux
        a(iLmin1:iLmax1,iLmin2:iLmax2) = flux(iLmin1:iLmax1,iLmin2:iLmax2,iw)
        b(iLmin1:iLmax1,iLmin2:iLmax2) = flux_min(iLmin1:iLmax1,iLmin2:iLmax2,&
           iw)
        c(iLmin1:iLmax1,iLmin2:iLmax2) = flux_max(iLmin1:iLmax1,iLmin2:iLmax2,&
           iw)
        call median(ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,iLmin2,iLmax1,&
           iLmax2, a, b, c, tmp) 
        wLCtmp(iLmin1:iLmax1,iLmin2:iLmax2,iw) = tmp(iLmin1:iLmax1,&
           iLmin2:iLmax2)
      end do
    end select

    !> right side
    !>> mmm -> pppp
    !>> mm  -> ppp
    !>> m   -> pp
    !>> 0   -> p
    !>> p   -> 0
    !>> pp  -> m
    !>> ppp -> mm
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       1) = u1_coeff(1) * w(iLppppmin1:iLppppmax1,iLppppmin2:iLppppmax2,&
       1:nwflux) + u1_coeff(2) * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       1:nwflux) + u1_coeff(3) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) + u1_coeff(4) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       2) = u2_coeff(1) * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       1:nwflux)  + u2_coeff(2) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux)  + u2_coeff(3) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)  + u2_coeff(4) * w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       3) = u3_coeff(1) * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux)   + u3_coeff(2) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)   + u3_coeff(3) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)   + u3_coeff(4) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,1:nwflux)
    f_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
       4) = u4_coeff(1) * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)    + u4_coeff(2) * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)    + u4_coeff(3) * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)  + u4_coeff(4) * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux)

    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,1) = w(iLppppmin1:iLppppmax1,&
       iLppppmin2:iLppppmax2,1:nwflux) * (547.d0 * w(iLppppmin1:iLppppmax1,&
       iLppppmin2:iLppppmax2,1:nwflux) - 3882.d0 * w(iLpppmin1:iLpppmax1,&
       iLpppmin2:iLpppmax2,1:nwflux) + 4642.d0 * w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,1:nwflux) - 1854.d0 * w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,1:nwflux)) + w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       1:nwflux) * (7043.d0 * w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       1:nwflux) - 17246.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) + 7042.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)) + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) * (11003.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) - 9402.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux)) + 2107.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux)**2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,2) = w(iLpppmin1:iLpppmax1,&
       iLpppmin2:iLpppmax2,1:nwflux) * (267.d0 * w(iLpppmin1:iLpppmax1,&
       iLpppmin2:iLpppmax2,1:nwflux) - 1642.d0 * w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,1:nwflux) + 1602.d0 * w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,1:nwflux) - 494.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux))  + w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) * (2843.d0 * w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       1:nwflux) - 5966.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) + 1922.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) * (3443.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - 2522.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux)) + 547.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) ** 2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,3) = w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,1:nwflux) * (547.d0 * w(iLppmin1:iLppmax1,&
       iLppmin2:iLppmax2,1:nwflux) - 2522.d0 * w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,1:nwflux) + 1922.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - 494.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux))  + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) * (3443.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       1:nwflux) - 5966.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) + 1602.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)) + w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) * (2843.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - 1642.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux)) + 267.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,1:nwflux) ** 2
    beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,4) = w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,1:nwflux) * (2107.d0 * w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,1:nwflux) - 9402.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) + 7042.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - 1854.d0 * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux))  + w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) * (11003.d0 * w(iLmin1:iLmax1,iLmin2:iLmax2,&
       1:nwflux) - 17246.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) + 4642.d0 * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux)) + w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) * (7043.d0 * w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       1:nwflux) - 3882.d0 * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux)) + 547.d0 * w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       1:nwflux) ** 2

    alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0 
    do i = 1,4
       alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
          i) = d_array(i)/(beta(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
          i) + weno_eps_machine)**2
       alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          1:nwflux) = alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,&
          1:nwflux) + alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,i)
    end do
    flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = 0.0d0
    do i = 1,4
       flux(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
          iLmin2:iLmax2,1:nwflux) + f_array(iLmin1:iLmax1,iLmin2:iLmax2,&
          1:nwflux,i) * alpha_array(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux,&
          i)/(alpha_sum(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
    end do

    select case(var)
    case(1)
      wRCtmp(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = flux(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux)
    case(2)
      idmax1=iLmax1+kr(idims,1);idmax2=iLmax2+kr(idims,2); idmin1=iLmin1
      idmin2=iLmin2;
      idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
      idmmax1=idmax1-kr(idims,1);idmmax2=idmax2-kr(idims,2);
      idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
      idpmax1=idmax1+kr(idims,1);idpmax2=idmax2+kr(idims,2);
  
      iemax1=idmax1;iemax2=idmax2; iemin1=idmin1-kr(idims,1)
      iemin2=idmin2-kr(idims,2);
      iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
      iemmax1=iemax1-kr(idims,1);iemmax2=iemax2-kr(idims,2);
      iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
      iepmax1=iemax1+kr(idims,1);iepmax2=iemax2+kr(idims,2);
      ieppmin1=iepmin1+kr(idims,1);ieppmin2=iepmin2+kr(idims,2)
      ieppmax1=iepmax1+kr(idims,1);ieppmax2=iepmax2+kr(idims,2);
  
      d(iemin1:iemax1,iemin2:iemax2,1:nwflux) = w(iemin1:iemax1,iemin2:iemax2,&
         1:nwflux)-2.0d0*w(iepmin1:iepmax1,iepmin2:iepmax2,&
         1:nwflux)+w(ieppmin1:ieppmax1,ieppmin2:ieppmax2,1:nwflux)
  
      do iw = 1,nwflux
        a(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idmin1:idmax1,idmin2:idmax2,&
           iw)-d(idmmin1:idmmax1,idmmin2:idmmax2,iw)
        b(idmin1:idmax1,idmin2:idmax2) = 4.0d0*d(idmmin1:idmmax1,&
           idmmin2:idmmax2,iw)-d(idmin1:idmax1,idmin2:idmax2,iw)
        call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,&
           idmax2,a,b,tmp)
        a(idmin1:idmax1,idmin2:idmax2) = d(idmin1:idmax1,idmin2:idmax2,iw)
        b(idmin1:idmax1,idmin2:idmax2) = d(idmmin1:idmmax1,idmmin2:idmmax2,iw)
        call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,&
           idmax2,a,b,tmp2)
        call minmod(ixImin1,ixImin2,ixImax1,ixImax2,idmin1,idmin2,idmax1,&
           idmax2,tmp,tmp2,tmp3)
        dm4(idmin1:idmax1,idmin2:idmax2,iw) = tmp3(idmin1:idmax1,&
           idmin2:idmax2)
      end do
   
      flux_ul(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = w(iLpmin1:iLpmax1,&
         iLpmin2:iLpmax2,1:nwflux) + mpalpha * (w(iLpmin1:iLpmax1,&
         iLpmin2:iLpmax2,1:nwflux) - w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         1:nwflux))
      flux_md(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux) = half * (w(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux) + w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         1:nwflux) - dm4(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux))
      flux_lc(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux) = half * (3.d0 * w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         1:nwflux) - w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
         1:nwflux)) + mpbeta / 3.d0 * dm4(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         1:nwflux)
    
      flux_min(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux) = max(min(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux),&
          w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux), flux_md(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux)), min(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         1:nwflux), flux_ul(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux),&
         flux_lc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)))
  
      flux_max(iLmin1:iLmax1,iLmin2:iLmax2,&
         1:nwflux) = min(max(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,1:nwflux),&
          w(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux), flux_md(iLmin1:iLmax1,&
         iLmin2:iLmax2,1:nwflux)), max(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
         1:nwflux), flux_ul(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux),&
         flux_lc(iLmin1:iLmax1,iLmin2:iLmax2,1:nwflux)))
      do iw=1,nwflux
        a(iLmin1:iLmax1,iLmin2:iLmax2) = flux(iLmin1:iLmax1,iLmin2:iLmax2,iw)
        b(iLmin1:iLmax1,iLmin2:iLmax2) = flux_min(iLmin1:iLmax1,iLmin2:iLmax2,&
           iw)
        c(iLmin1:iLmax1,iLmin2:iLmax2) = flux_max(iLmin1:iLmax1,iLmin2:iLmax2,&
           iw)
        call median(ixImin1,ixImin2,ixImax1,ixImax2, iLmin1,iLmin2,iLmax1,&
           iLmax2, a, b, c, tmp) 
        wRCtmp(iLmin1:iLmax1,iLmin2:iLmax2,iw) = tmp(iLmin1:iLmax1,&
           iLmin2:iLmax2)
      end do
    end select

    call fix_limiter(ixImin1,ixImin2,ixImax1,ixImax2,iLmin1,iLmin2,iLmax1,&
       iLmax2,wLCtmp,wRCtmp,wLC,wRC)

  end subroutine WENO7limiter

  subroutine minmod(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,a,b,minm)

    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: a(ixImin1:ixImax1,ixImin2:ixImax2),&
        b(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(out):: minm(ixImin1:ixImax1,ixImin2:ixImax2)

    minm(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (sign(one,a(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))+sign(one,b(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)))/2.0d0 * min(abs(a(ixOmin1:ixOmax1,ixOmin2:ixOmax2)),&
       abs(b(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))

  end subroutine minmod

  subroutine median(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,a,b,c,med)

    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: a(ixImin1:ixImax1,ixImin2:ixImax2),&
        b(ixImin1:ixImax1,ixImin2:ixImax2), c(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(out):: med(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision             :: tmp1(ixImin1:ixImax1,ixImin2:ixImax2),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2)

    tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = b(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) - a(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = c(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) - a(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    med(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = a(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) + (sign(one,tmp1(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))+sign(one,tmp2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)))/2.0d0 * min(abs(tmp1(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)),abs(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))

  end subroutine median
end module mod_weno
