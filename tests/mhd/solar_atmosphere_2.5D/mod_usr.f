module mod_usr
  use mod_mhd
  implicit none
  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: usr_grav,vmax,La
  double precision :: heatunit,gzone,B0,theta,SRadius,kx,ly,bQ0,dya
  double precision, allocatable :: pa(:),ra(:)
  integer, parameter :: jmax=8000

contains

  subroutine usr_init()
    call set_coordinate_system("Cartesian_2.5D")

    unit_length        = 1.d9                                         ! cm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9 !cm-3,cm-3

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_source          => special_source
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_set_B0          => specialset_B0
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output 
    usr_init_vector_potential=>initvecpot_usr
    usr_set_wLR         => boundary_wLR
    !usr_refine_threshold=> specialthreshold
    !usr_set_electric_field => boundary_electric_field

    call mhd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr()
    heatunit=unit_pressure/unit_time !3.697693390805347E-003 erg*cm-3/s,erg*cm-3/s

    usr_grav=-2.74d4*unit_length/unit_velocity**2 ! solar gravity
    bQ0=1.d-4/heatunit ! background heating power density
    gzone=0.2d0 ! thickness of a ghostzone below the bottom boundary
    dya=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax) !cells size of high-resolution 1D solar atmosphere
    B0=Busr/unit_magneticfield ! magnetic field strength at the bottom
    theta=60.d0*dpi/180.d0 !the angle to the plane xy, 90-theta is the angle to the polarity inversion line of the arcade 
    kx=dpi/(xprobmax1-xprobmin1)
    ly=kx*dcos(theta)
    SRadius=69.61d0 ! Solar radius
    vmax=7.d5/unit_velocity ! maximal driven velocity
    La=1.5d9/unit_length
    ! hydrostatic vertical stratification of density, temperature, pressure
    call inithdstatic

  end subroutine initglobaldata_usr

  subroutine inithdstatic
    use mod_solar_atmosphere
    ! initialize the table in a vertical line through the global domain
    integer :: j,na,ibc
    double precision, allocatable :: Ta(:),gg(:),ya(:)
    double precision :: rpho,Ttop,Tpho,wtra,res,rhob,pb,htra,Ttr,Fc,invT,kappa
    double precision :: rhohc,hc
    logical :: simple_temperature_curve

    simple_temperature_curve=.true.

    allocate(ya(jmax),Ta(jmax),gg(jmax),pa(jmax),ra(jmax))

    if(simple_temperature_curve) then
      rpho=1.151d15/unit_numberdensity !number density at the bottom of height table
      Tpho=8.d3/unit_temperature ! temperature of chromosphere
      Ttop=1.5d6/unit_temperature ! estimated temperature in the top
      htra=0.2d0 ! height of initial transition region
      wtra=0.02d0 ! width of initial transition region 
      Ttr=1.6d5/unit_temperature ! lowest temperature of upper profile
      Fc=2.d5/heatunit/unit_length ! constant thermal conduction flux
      kappa=8.d-7*unit_temperature&
         **3.5d0/unit_length/unit_density/unit_velocity**3
      do j=1,jmax
         ya(j)=(dble(j)-0.5d0)*dya-gzone
         if(ya(j)>htra) then
           Ta(j)=(3.5d0*Fc/kappa*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
         else
           Ta(j)=Tpho+0.5d0*(Ttop-Tpho)*(tanh((ya(j)-htra-0.027d0)/wtra)+1.d0)
         endif
         gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
      enddo
      !! solution of hydrostatic equation 
      ra(1)=rpho
      pa(1)=rpho*Tpho
      invT=gg(1)/Ta(1)
      invT=0.d0
      do j=2,jmax
         invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
         pa(j)=pa(1)*dexp(invT*dya)
         ra(j)=pa(j)/Ta(j)
      end do
    else
      do j=1,jmax
         ! get height table
         ya(j)=(dble(j)-0.5d0)*dya-gzone
         ! get gravity table
         gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
      enddo
      ! a coronal height at 10 Mm
      hc=1.d9/unit_length
      ! the number density at the coronal height
      rhohc=6.2d8/unit_numberdensity
      ! get density and pressure table in hydrostatic state with a preset temperature table
      call get_atm_para(ya,ra,pa,gg,jmax,'AL-C7',hc,rhohc)
    end if
    deallocate(ya,gg,Ta)
    !! initialized rho and p in the fixed bottom boundary
    na=floor(gzone/dya+0.5d0)
    res=gzone-(dble(na)-0.5d0)*dya
    rhob=ra(na)+res/dya*(ra(na+1)-ra(na))
    pb=pa(na)+res/dya*(pa(na+1)-pa(na))
    allocate(rbc(nghostcells))
    allocate(pbc(nghostcells))
    do ibc=nghostcells,1,-1
      na=floor((gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+&
         1)-0.5d0))/dya+0.5d0)
      res=gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+&
         1)-0.5d0)-(dble(na)-0.5d0)*dya
      rbc(ibc)=ra(na)+res/dya*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+res/dya*(pa(na+1)-pa(na))
    end do

    if (mype==0) then
     print*,'minra',minval(ra)
     print*,'rhob',rhob
     print*,'pb',pb
    endif

  end subroutine inithdstatic

  subroutine initonegrid_usr(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x)
    ! initialize one grid
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: res
    integer :: ix1,ix2,na
    logical, save :: first=.true.

    if(first)then
      if(mype==0) then
        write(*,*)'Simulating 2.5D solar atmosphere'
      endif
      first=.false.
    endif
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
        na=floor((x(ix1,ix2,2)-xprobmin2+gzone)/dya+0.5d0)
        res=x(ix1,ix2,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dya
        w(ix1,ix2,rho_)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
        w(ix1,ix2,p_)  =pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
    end do
    end do
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(:))=zero
    if(B0field) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:))=zero
    else if(stagger_grid) then
      call b_from_vector_potential(ixGslo1,ixGslo2,ixGshi1,ixGshi2,ixImin1,&
         ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,block%ws,x)
      call mhd_face_to_center(ixOmin1,ixOmin2,ixOmax1,ixOmax2,block)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(3))=-B0*dcos(kx*x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         2))*dsin(theta)
    else
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(1))=-B0*dcos(kx*x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         2))*dcos(theta)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(2))= B0*dsin(kx*x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(3))=-B0*dcos(kx*x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         2))*dsin(theta)
    endif

    call mhd_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w,x)

  end subroutine initonegrid_usr

  subroutine initvecpot_usr(ixImin1,ixImin2,ixImax1,ixImax2, ixCmin1,ixCmin2,&
     ixCmax1,ixCmax2, xC, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixCmin1,ixCmin2,ixCmax1,ixCmax2,idir
    double precision, intent(in)       :: xC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(out)      :: A(ixImin1:ixImax1,ixImin2:ixImax2)

    if (idir==3) then
      A(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = B0/ly*dcos(kx*xC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,1))*dexp(-ly*xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2))*dcos(theta)
    else
      A(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = 0.d0
    end if

  end subroutine initvecpot_usr

  subroutine boundary_electric_field(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,qt,qdt,fE,s)
    ! specify tangential electric field at physical boundaries 
    ! to fix or drive normal magnetic field
    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: qt,qdt
    type(state)                        :: s
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndim:3)

    double precision :: xC(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
        tmp(ixImin1:ixImax1,ixImin2:ixImax2)
    integer :: idir,ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixAmin1,ixAmin2,ixAmax1,&
       ixAmax2

    if(s%is_physical_boundary(3)) then
      if(iprob<3) then
        ! to fix normal magnetic field at bottom boundary surface
        idir=3
        ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1;
        ixCmax1=ixOmax1;ixCmax2=ixOmax2;
        fE(ixCmin1:ixCmax1,ixCmin2,3)=0.d0
      else
        ! horizontal flow driven boundary
        ixCmin1=ixOmin1-1;ixCmin2=ixOmin2-1;
        ixCmax1=ixOmax1;ixCmax2=ixOmax2;
        ixAmin1=ixCmin1;ixAmin2=ixCmin2;
        ixAmax1=ixCmax1+1;ixAmax2=ixCmax2+1;
        !! cell center electric field
        tmp(ixAmin1:ixAmax1,ixAmin2)=-s%ws(ixAmin1:ixAmax1,ixAmin2,&
           2)*s%w(ixAmin1:ixAmax1,ixAmin2,mom(1))/s%w(ixAmin1:ixAmax1,ixAmin2,&
           rho_)
        ! neighor cell center average to get cell edge
        ixAmin1=ixCmin1+kr(1,1);ixAmin2=ixCmin2+kr(1,2);
        ixAmax1=ixCmax1+kr(1,1);ixAmax2=ixCmax2+kr(1,2);
        !print*,'difference',maxval(abs(fE(ixCmin2^%2ixC^S,3)-vx(ixCmin2^%2ixC^S)))
        fE(ixCmin1:ixCmax1,ixCmin2,3)=0.5d0*(tmp(ixCmin1:ixCmax1,&
           ixCmin2)+tmp(ixAmin1:ixAmax1,ixAmin2))
      end if
    end if

  end subroutine boundary_electric_field

  !> allow user to specify variables' left and right state at physical boundaries to control flux through the boundary surface 
  subroutine boundary_wLR(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,qt,wLC,wRC,wLp,wRp,s,idir)
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    type(state)                     :: s

    double precision :: vx(ixImin1:ixImax1,ixImin2:ixImax2)

    if(iprob<3) then
      if(s%is_physical_boundary(3).and.idir==2) then
        wLp(ixOmin1:ixOmax1,ixOmin2,mom(:))=0.d0
        wRp(ixOmin1:ixOmax1,ixOmin2,mom(:))=0.d0
        wLC(ixOmin1:ixOmax1,ixOmin2,mom(:))=0.d0
        wRC(ixOmin1:ixOmax1,ixOmin2,mom(:))=0.d0
      end if
    else
      if(s%is_physical_boundary(3).and.idir==2) then
        call driven_velocity(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,qt,s%x,vx)
        wLp(ixOmin1:ixOmax1,ixOmin2,mom(1))=vx(ixOmin1:ixOmax1,ixOmin2)
        wRp(ixOmin1:ixOmax1,ixOmin2,mom(1))=wRp(ixOmin1:ixOmax1,ixOmin2,&
           mom(1))
        wLC(ixOmin1:ixOmax1,ixOmin2,mom(1))=wLp(ixOmin1:ixOmax1,ixOmin2,&
           mom(1))*wLp(ixOmin1:ixOmax1,ixOmin2,rho_)
        wRC(ixOmin1:ixOmax1,ixOmin2,mom(1))=wRp(ixOmin1:ixOmax1,ixOmin2,&
           mom(1))*wLp(ixOmin1:ixOmax1,ixOmin2,rho_)
        wLp(ixOmin1:ixOmax1,ixOmin2,mom(2:3))=0.d0
        wRp(ixOmin1:ixOmax1,ixOmin2,mom(2:3))=0.d0
        wLC(ixOmin1:ixOmax1,ixOmin2,mom(2:3))=0.d0
        wRC(ixOmin1:ixOmax1,ixOmin2,mom(2:3))=0.d0
      end if
    end if

  end subroutine boundary_wLR

  subroutine driven_velocity(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,qt,x,vdr)
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(out) :: vdr(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision :: tramp1,ft

    tramp1=3.d0
    if(qt<tramp1) then
      ft=qt/tramp1
    else
      ft=1.d0
    end if
    where(abs(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1))<=La)
      vdr(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-ft*vmax*sin(dpi*x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1)/La)
    else where
      vdr(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=0.d0
    end where

  end subroutine driven_velocity

  subroutine specialbound_usr(qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixOmin1,ixOmin2,ixOmax1,ixOmax2, iB, ixImin1,&
       ixImin2,ixImax1,ixImax2
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2),&
       tmp(ixImin1:ixImax1,ixImin2:ixImax2),ggrid(ixImin1:ixImax1,&
       ixImin2:ixImax2),invT(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: Q(ixImin1:ixImax1,ixImin2:ixImax2),Qp(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    integer :: ix1,ix2,ixOsmin1,ixOsmin2,ixOsmax1,ixOsmax2,ixCmin1,ixCmin2,&
       ixCmax1,ixCmax2,hxCmin1,hxCmin2,hxCmax1,hxCmax2,jxOmin1,jxOmin2,jxOmax1,&
       jxOmax2,idir

    select case(iB)
    case(3)
      if(iprob<3) then
        !! fixed zero velocity
        do idir=1,ndir
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=-w(ixOmin1:ixOmax1,&
             ixOmax2+nghostcells:ixOmax2+1:-1,mom(idir))/w(ixOmin1:ixOmax1,&
             ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
        end do
      else
        call driven_velocity(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,qt,x,w(ixImin1:ixImax1,ixImin2:ixImax2,mom(1)))
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2))=-w(ixOmin1:ixOmax1,&
           ixOmax2+nghostcells:ixOmax2+1:-1,mom(2))/w(ixOmin1:ixOmax1,&
           ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(3))=-w(ixOmin1:ixOmax1,&
           ixOmax2+nghostcells:ixOmax2+1:-1,mom(3))/w(ixOmin1:ixOmax1,&
           ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
      end if
      !! fixed b1 b2 b3
      if(iprob==0 .or. B0field) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:))=0.d0
      else if(stagger_grid) then
        do idir=1,nws
          if(idir==2) cycle
          ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;
          ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir);
          ! 2nd order one-sided zero-gradient extrapolation
          !do ix2=ixOsmax2,ixOsmin2,-1
          !   block%ws(ix2^%2ixOs^S,idir)=1.d0/3.d0*&
          !         (-block%ws(ix2+2^%2ixOs^S,idir)&
          !     +4.d0*block%ws(ix2+1^%2ixOs^S,idir))
          !end do
          ! 4th order one-sided equal-gradient extrapolation
          do ix2=ixOsmax2,ixOsmin2,-1
            block%ws(ixOsmin1:ixOsmax1,ix2,&
               idir)= 0.12d0*block%ws(ixOsmin1:ixOsmax1,ix2+5,&
               idir) -0.76d0*block%ws(ixOsmin1:ixOsmax1,ix2+4,&
               idir) +2.08d0*block%ws(ixOsmin1:ixOsmax1,ix2+3,&
               idir) -3.36d0*block%ws(ixOsmin1:ixOsmax1,ix2+2,&
               idir) +2.92d0*block%ws(ixOsmin1:ixOsmax1,ix2+1,idir)
          end do
        end do
        ixOsmin1=ixOmin1-kr(2,1);ixOsmin2=ixOmin2-kr(2,2)
        ixOsmax1=ixOmax1-kr(2,1);ixOsmax2=ixOmax2-kr(2,2);
        jxOmin1=ixOmin1+nghostcells*kr(2,1)
        jxOmin2=ixOmin2+nghostcells*kr(2,2)
        jxOmax1=ixOmax1+nghostcells*kr(2,1)
        jxOmax2=ixOmax2+nghostcells*kr(2,2);
        block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,2)=zero
        do ix2=ixOsmax2,ixOsmin2,-1
          call get_divb(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
             ixOmax1,ixOmax2,Qp)
          block%ws(ixOsmin1:ixOsmax1,ix2,2)=Qp(ixOmin1:ixOmax1,&
             ix2+1)*block%dvolume(ixOmin1:ixOmax1,&
             ix2+1)/block%surfaceC(ixOsmin1:ixOsmax1,ix2,2)
        end do
        call mhd_face_to_center(ixOmin1,ixOmin2,ixOmax1,ixOmax2,block)
        if(iprob<3) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mag(3))=-B0*dcos(kx*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))*dsin(theta)
        else
          do ix2=ixOmax2,ixOmin2,-1
            w(ixOmin1:ixOmax1,ix2,mag(3))= 0.12d0*w(ixOmin1:ixOmax1,ix2+5,&
               mag(3)) -0.76d0*w(ixOmin1:ixOmax1,ix2+4,&
               mag(3)) +2.08d0*w(ixOmin1:ixOmax1,ix2+3,&
               mag(3)) -3.36d0*w(ixOmin1:ixOmax1,ix2+2,&
               mag(3)) +2.92d0*w(ixOmin1:ixOmax1,ix2+1,mag(3))
          end do
        end if
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(1))=-B0*dcos(kx*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))*dcos(theta)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(2))= B0*dsin(kx*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(3))=-B0*dcos(kx*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))*dsin(theta)
      endif
      !! fixed gravity stratification of density and pressure pre-determined in initial condition
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
        w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
      enddo
      call mhd_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,w,x)
    case(4)
      ixOsmin1=ixOmin1;ixOsmin2=ixOmin2;ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmin2-1;
      call mhd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOsmin1,&
         ixOsmin2,ixOsmax1,ixOsmax2,pth)
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmax2;
      call getggrav(ggrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOsmin1,ixOsmin2,&
         ixOsmax1,ixOsmax2,x)
      !> fill pth, rho ghost layers according to gravity stratification
      invT(ixOmin1:ixOmax1,ixOmin2-1)=w(ixOmin1:ixOmax1,ixOmin2-1,&
         rho_)/pth(ixOmin1:ixOmax1,ixOmin2-1)
      tmp=0.d0
      do ix2=ixOmin2,ixOmax2
        tmp(ixOmin1:ixOmax1,ixOmin2-1)=tmp(ixOmin1:ixOmax1,&
           ixOmin2-1)+0.5d0*(ggrid(ixOmin1:ixOmax1,ix2)+ggrid(ixOmin1:ixOmax1,&
           ix2-1))*invT(ixOmin1:ixOmax1,ixOmin2-1)
        w(ixOmin1:ixOmax1,ix2,p_)=pth(ixOmin1:ixOmax1,&
           ixOmin2-1)*dexp(tmp(ixOmin1:ixOmax1,ixOmin2-1)*dxlevel(2))
        w(ixOmin1:ixOmax1,ix2,rho_)=w(ixOmin1:ixOmax1,ix2,&
           p_)*invT(ixOmin1:ixOmax1,ixOmin2-1)
      enddo
      !> fixed zero velocity
      do idir=1,ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir)) =-w(ixOmin1:ixOmax1,&
           ixOmin2-1:ixOmin2-nghostcells:-1,mom(idir))/w(ixOmin1:ixOmax1,&
           ixOmin2-1:ixOmin2-nghostcells:-1,rho_)
      end do
      !> zero normal gradient extrapolation
      if(stagger_grid) then
        do idir=1,nws
          if(idir==2) cycle
          ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;
          ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir);
          do ix2=ixOsmin2,ixOsmax2
             block%ws(ixOsmin1:ixOsmax1,ix2,&
                idir)=1.d0/3.d0*(-block%ws(ixOsmin1:ixOsmax1,ix2-2,&
                idir)+4.d0*block%ws(ixOsmin1:ixOsmax1,ix2-1,idir))
          end do
        end do
        ixOsmin1=ixOmin1;ixOsmin2=ixOmin2;ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;
        jxOmin1=ixOmin1-nghostcells*kr(2,1)
        jxOmin2=ixOmin2-nghostcells*kr(2,2)
        jxOmax1=ixOmax1-nghostcells*kr(2,1)
        jxOmax2=ixOmax2-nghostcells*kr(2,2);
        block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,2)=zero
        call get_divb(w,ixImin1,ixImin2,ixImax1,ixImax2,jxOmin1,jxOmin2,&
           jxOmax1,jxOmax2,Q)
        do ix2=ixOsmin2,ixOsmax2
          call get_divb(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
             ixOmax1,ixOmax2,Qp)
          block%ws(ixOsmin1:ixOsmax1,ix2,2)=(Q(jxOmin1:jxOmax1,&
             jxOmax2)*block%dvolume(jxOmin1:jxOmax1,&
             jxOmax2)-Qp(ixOmin1:ixOmax1,ix2)*block%dvolume(ixOmin1:ixOmax1,&
             ix2))/block%surfaceC(ixOsmin1:ixOsmax1,ix2,2)
        end do
        call mhd_face_to_center(ixOmin1,ixOmin2,ixOmax1,ixOmax2,block)
        do ix2=ixOmin2,ixOmax2
          w(ixOmin1:ixOmax1,ix2,mag(3))=(1.0d0/3.0d0)* (-w(ixOmin1:ixOmax1,&
             ix2-2,mag(3))+4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(3)))
        enddo
      else
        do ix2=ixOmin2,ixOmax2
          w(ixOmin1:ixOmax1,ix2,mag(:))=(1.0d0/3.0d0)* (-w(ixOmin1:ixOmax1,&
             ix2-2,mag(:))+4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(:)))
        enddo
      end if
      call mhd_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select
    
  end subroutine specialbound_usr

  subroutine gravity(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,wCT,x,gravity_field)
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    double precision, intent(out)   :: gravity_field(ixImin1:ixImax1,&
       ixImin2:ixImax2,ndim)

    double precision                :: ggrid(ixImin1:ixImax1,ixImin2:ixImax2)

    gravity_field=0.d0
    call getggrav(ggrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,x)
    gravity_field(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)=ggrid(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

  end subroutine gravity

  subroutine getggrav(ggrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,x)
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(out)   :: ggrid(ixImin1:ixImax1,ixImin2:ixImax2)

    ggrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=usr_grav*(SRadius/(SRadius+&
       x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))**2
  end subroutine

  subroutine special_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, iwmin,iwmax
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
        wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: lQgrid(ixImin1:ixImax1,ixImin2:ixImax2),&
       bQgrid(ixImin1:ixImax1,ixImin2:ixImax2)

    ! add global background heating bQ
    call getbQ(bQgrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,qtC,wCT,x)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_)+qdt*bQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    if(iprob==2) then
      call getlQ(lQgrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,qt,wCT,x)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_)+qdt*lQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end if

  end subroutine special_source

  subroutine getbQ(bQgrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,qt,w,x)
  ! calculate background heating bQ
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim), w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: bQgrid(ixImin1:ixImax1,ixImin2:ixImax2)

    bQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=bQ0*dexp(-x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,2)/5.d0)

  end subroutine getbQ

  subroutine getlQ(lQgrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,qt,w,x)
  ! calculate localized heating lQ
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim), w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: lQgrid(ixImin1:ixImax1,ixImin2:ixImax2),&
        A(ixImin1:ixImax1,ixImin2:ixImax2), lQd, lQ0, hp, lQt,Atop,Alow

    lQ0=1.d-2/heatunit
    lQt=5.d2/unit_time
    lQd=0.2d0
    hp=0.3d0

    call initvecpot_usr(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, x, A, 3)

    lQgrid=0.d0
    Atop=B0/kx*cos(kx*1.35d0)
    Alow=B0/kx*cos(kx*1.75d0)
    where(A(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<Atop .and. A(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)>Alow)
      lQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=lQ0
    end where
    where(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)>hp)
      lQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=lQgrid(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*dexp(-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         2)-hp)**2/lQd)
    end where
    if(qt<lQt) then
      lQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=qt/lQt*lQgrid(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)
    endif

  end subroutine getlQ

  subroutine special_refine_grid(igrid,level,ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,qt,w,x,refine,coarsen)
  ! Enforce additional refinement or coarsening
  ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    integer, intent(in) :: igrid, level, ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: qt, w(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! fix the bottom layer to the highest level
    if (any(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)<=xprobmin2+0.05d0)) then
      refine=1
      coarsen=-1
    endif

  end subroutine special_refine_grid

  subroutine specialvar_output(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_radiative_cooling
    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2),&
       B2(ixImin1:ixImax1,ixImin2:ixImax2),tmp2(ixImin1:ixImax1,&
       ixImin2:ixImax2),dRdT(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: ens(ixImin1:ixImax1,ixImin2:ixImax2),&
       divb(ixImin1:ixImax1,ixImin2:ixImax2),wlocal(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw)
    double precision :: Btotal(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir),&
       curlvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)
    double precision :: Te(ixImin1:ixImax1,ixImin2:ixImax2),tco_local
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim) :: gradT, bunitvec
    integer :: idirmin,idir,ix1,ix2
    logical :: lrlt(ixImin1:ixImax1,ixImin2:ixImax2)

    wlocal(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw)
    ! output temperature
    call mhd_get_pthermal(wlocal,x,ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,&
       ixImin2,ixImax1,ixImax2,pth)
    Te(ixImin1:ixImax1,ixImin2:ixImax2)=pth(ixImin1:ixImax1,&
       ixImin2:ixImax2)/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1)=Te(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    do idir=1,ndir
      if(B0field) then
        Btotal(ixImin1:ixImax1,ixImin2:ixImax2,idir)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,mag(idir))+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,&
           idir,0)
      else
        Btotal(ixImin1:ixImax1,ixImin2:ixImax2,idir)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,mag(idir))
      endif
    end do
    ! B^2
    B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sum((Btotal(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,:))**2,dim=ndim+1)

    ! output Alfven wave speed B/sqrt(rho)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2)=dsqrt(B2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))

    ! output divB1
    call get_divb(wlocal,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,divb)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+3)=divb(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    ! output the plasma beta p*2/B**2
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+4)=pth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*two/B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    ! output heating rate
    call getbQ(ens,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,global_time,wlocal,x)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+5)=ens(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    ! store the cooling rate 
    if(mhd_radiative_cooling)call getvar_cooling(ixImin1,ixImin2,ixImax1,&
       ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wlocal,x,ens,rc_fl)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+6)=ens(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    ! store current
    call get_current(wlocal,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,idirmin,curlvec)
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+6+idir)=curlvec(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,idir)
    end do

    if(nw_extra>0) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       nw+10)=block%wextra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw_tcoff)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    if(nw_extra>0) then
      varnames='Te Alfv divB beta bQ rad j1 j2 j3 Tcutoff'
    else
      varnames='Te Alfv divB beta bQ rad j1 j2 j3'
    end if

  end subroutine specialvarnames_output

  subroutine specialset_B0(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,x,wB0)
  ! Here add a steady (time-independent) potential or 
  ! linear force-free background field
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: wB0(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndir)

    wB0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=-B0*dcos(kx*x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       2))*dcos(theta)
    wB0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)=+B0*dsin(kx*x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
    wB0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)=-B0*dcos(kx*x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       2))*dsin(theta)

  end subroutine specialset_B0

  subroutine specialthreshold(wlocal,xlocal,tolerance,qt,level)
    !PURPOSE: use different tolerance in special regions for AMR to
    !reduce/increase resolution there where nothing/something interesting happens.
    use mod_global_parameters

    double precision, intent(in) :: wlocal(1:nw),xlocal(1:ndim),qt
    double precision, intent(inout) :: tolerance
    integer, intent(in) :: level

    double precision :: bczone1,bczone2,addtol,tol_add

    tol_add=0.d0
    !amplitude of additional tolerance
    addtol=0.6d0
    ! thickness of near-boundary region
    bczone1=0.8d0
    bczone2=2.d0
    ! linear changing of additional tolerance
    if(xlocal(1)-xprobmin1 < bczone1 .or. xprobmax1-xlocal(1) < bczone1) then
      tol_add=(1.d0-min(xlocal(1)-xprobmin1,&
         xprobmax1-xlocal(1))/bczone1)*addtol
    endif
    if(xprobmax2-xlocal(2) < bczone2) then
      tol_add=(1.d0-(xprobmax2-xlocal(2))/bczone2)*addtol
    endif
    tolerance=tolerance+tol_add

  end subroutine specialthreshold

end module mod_usr
