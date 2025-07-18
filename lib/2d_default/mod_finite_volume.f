!> Module with finite volume methods for fluxes
module mod_finite_volume
#include "amrvac.h"
  implicit none
  private

  public :: finite_volume
  public :: hancock
  public :: reconstruct_LR

contains

  !> The non-conservative Hancock predictor for TVDLF
  !>
  !> on entry:
  !> input available on ixI^L=ixG^L asks for output on ixO^L=ixG^L^LSUBnghostcells
  !> one entry: (predictor): wCT -- w_n        wnew -- w_n   qdt=dt/2
  !> on exit :  (predictor): wCT -- w_n        wnew -- w_n+1/2
  subroutine hancock(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idimsmin,idimsmax,qtC,sCT,qt,snew,dxs,x)
    use mod_physics
    use mod_global_parameters
    use mod_source, only: addsource2

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, idimsmin,idimsmax
    double precision, intent(in) :: qdt, qtC, qt, dxs(ndim), x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndim)
    type(state) :: sCT, snew

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wprim,&
        wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLp,&
        wRp
    double precision, dimension(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)      :: inv_volume
    double precision :: fLC(ixImin1:ixImax1,ixImin2:ixImax2, nwflux),&
        fRC(ixImin1:ixImax1,ixImin2:ixImax2, nwflux)
    double precision :: dxinv(1:ndim)
    integer :: idims, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
       hxOmax2
    logical :: active

    associate(wCT=>sCT%w,wnew=>snew%w)
    ! Expand limits in each idims direction in which fluxes are added
    ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
    do idims= idimsmin,idimsmax
       ixmin1=ixmin1-kr(idims,1);ixmin2=ixmin2-kr(idims,2)
       ixmax1=ixmax1+kr(idims,1);ixmax2=ixmax2+kr(idims,2);
    end do
    if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2<ixmax2) &
       call mpistop("Error in Hancock: Nonconforming input limits")

    wprim=wCT
    call phys_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
       ixImax1,ixImax2,wprim,x)

    dxinv=-qdt/dxs
    if(.not.slab_uniform) inv_volume = 1.d0/block%dvolume(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    do idims= idimsmin,idimsmax
      b0i=idims
      ! Calculate w_j+g_j/2 and w_j-g_j/2
      ! First copy all variables, then upwind wLC and wRC.
      ! wLC is to the left of ixO, wRC is to the right of wCT.
      hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
      hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);

      wRp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,1:nwflux)=wprim(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1:nwflux)
      wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)=wprim(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1:nwflux)

      ! apply limited reconstruction for left and right status at cell interfaces
      call reconstruct_LR(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,hxOmin1,hxOmin2,hxOmax1,hxOmax2,idims,wprim,wLC,wRC,&
         wLp,wRp,x,dxs(idims))

      ! Calculate the fLC and fRC fluxes
      call phys_get_flux(wRC,wRp,x,ixImin1,ixImin2,ixImax1,ixImax2,hxOmin1,&
         hxOmin2,hxOmax1,hxOmax2,idims,fRC)
      call phys_get_flux(wLC,wLp,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,idims,fLC)

      ! Advect w(iw)
      if (slab_uniform) then
        do iw=1,nwflux
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw)+dxinv(idims)* (fLC(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2, iw)-fRC(hxOmin1:hxOmax1,hxOmin2:hxOmax2, iw))
        end do
      else
        do iw=1,nwflux
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,iw) - qdt * inv_volume &
             *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             idims)*fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              iw) -block%surfaceC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
             idims)*fRC(hxOmin1:hxOmax1,hxOmin2:hxOmax2, iw))
        end do
      end if
    end do ! next idims
    b0i=0

    if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,wnew,x)

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,wnew,x,&
       .false.,active,wprim)

    ! check and optionally correct unphysical values
    if(fix_small_values) then
       call phys_handle_small_values(.false.,wnew,x,ixImin1,ixImin2,ixImax1,&
          ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
          'exit hancock finite_volume')
    endif
    end associate
  end subroutine hancock

  !> finite volume method
  subroutine finite_volume(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idimsmin,idimsmax, qtC,sCT,qt,snew,sold,fC,fE,dxs,&
     x)
    use mod_physics
    use mod_variables
    use mod_global_parameters
    use mod_tvd, only:tvdlimit2
    use mod_source, only: addsource2
    use mod_usr_methods

    integer, intent(in)                                   :: method
    double precision, intent(in)                          :: qdt, qtC, qt,&
        dxs(ndim)
    integer, intent(in)                                   :: ixImin1,ixImin2,&
       ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idimsmin,idimsmax
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
        intent(in) :: x
    type(state)                                           :: sCT, snew, sold
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,&
       1:ndim)    :: fC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndim:3)         :: fE

    ! primitive w at cell center
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wprim
    ! left and right constructed status in conservative form
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLC,&
        wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLp,&
        wRp
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nwflux) :: fLC, fRC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:number_species)      :: cmaxC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:number_species)      :: cminC
    double precision, dimension(ixImin1:ixImax1,&
       ixImin2:ixImax2)      :: Hspeed
    double precision, dimension(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)      :: inv_volume
    double precision, dimension(1:ndim)     :: dxinv
    integer, dimension(ixImin1:ixImax1,ixImin2:ixImax2)               :: &
       patchf
    integer :: idims, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
       hxOmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2, ixCRmin1,ixCRmin2,ixCRmax1,&
       ixCRmax2, kxCmin1,kxCmin2,kxCmax1,kxCmax2, kxRmin1,kxRmin2,kxRmax1,&
       kxRmax2, ii
    logical :: active
    type(ct_velocity) :: vcts

    associate(wCT=>sCT%w, wnew=>snew%w, wold=>sold%w)

    fC=0.d0
    fLC=0.d0
    fRC=0.d0

    ! The flux calculation contracts by one in the idims direction it is applied.
    ! The limiter contracts the same directions by one more, so expand ixO by 2.
    ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
    do idims= idimsmin,idimsmax
       ixmin1=ixmin1-2*kr(idims,1);ixmin2=ixmin2-2*kr(idims,2)
       ixmax1=ixmax1+2*kr(idims,1);ixmax2=ixmax2+2*kr(idims,2);
    end do
    if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2<ixmax2) &
       call mpistop("Error in fv : Nonconforming input limits")

    wprim=wCT
    call phys_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
       ixImax1,ixImax2,wprim,x)

    do idims= idimsmin,idimsmax
       ! use interface value of w0 at idims
       b0i=idims

       hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
       hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);

       kxCmin1=ixImin1;kxCmin2=ixImin2; kxCmax1=ixImax1-kr(idims,1)
       kxCmax2=ixImax2-kr(idims,2);
       kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
       kxRmax1=kxCmax1+kr(idims,1);kxRmax2=kxCmax2+kr(idims,2);

       if(stagger_grid) then
         ! ct needs all transverse cells
         ixCmax1=ixOmax1+nghostcells-nghostcells*kr(idims,1)
         ixCmax2=ixOmax2+nghostcells-nghostcells*kr(idims,2);
         ixCmin1=hxOmin1-nghostcells+nghostcells*kr(idims,1)
         ixCmin2=hxOmin2-nghostcells+nghostcells*kr(idims,2);
       else
         ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
         ixCmax1=ixOmax1;ixCmax2=ixOmax2; ixCmin1=hxOmin1;ixCmin2=hxOmin2;
       end if

       ! wRp and wLp are defined at the same locations, and will correspond to
       ! the left and right reconstructed values at a cell face. Their indexing
       ! is similar to cell-centered values, but in direction idims they are
       ! shifted half a cell towards the 'lower' direction.
       wRp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nw)=wprim(kxRmin1:kxRmax1,&
          kxRmin2:kxRmax2,1:nw)
       wLp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nw)=wprim(kxCmin1:kxCmax1,&
          kxCmin2:kxCmax2,1:nw)

       ! Determine stencil size
       ixCRmin1 = max(ixCmin1 - phys_wider_stencil,ixGlo1)
       ixCRmin2 = max(ixCmin2 - phys_wider_stencil,ixGlo2)
       ixCRmax1 = min(ixCmax1 + phys_wider_stencil,ixGhi1)
       ixCRmax2 = min(ixCmax2 + phys_wider_stencil,ixGhi2)

       ! apply limited reconstruction for left and right status at cell interfaces
       call reconstruct_LR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,&
          ixCRmax1,ixCRmax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,idims,wprim,&
          wLC,wRC,wLp,wRp,x,dxs(idims))

       ! special modification of left and right status before flux evaluation
       call phys_modify_wLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,&
          ixCRmax1,ixCRmax2,qt,wLC,wRC,wLp,wRp,sCT,idims)

       ! evaluate physical fluxes according to reconstructed status
       call phys_get_flux(wLC,wLp,x,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
          ixCmin2,ixCmax1,ixCmax2,idims,fLC)
       call phys_get_flux(wRC,wRp,x,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
          ixCmin2,ixCmax1,ixCmax2,idims,fRC)

       if(H_correction) then
         call phys_get_H_speed(wprim,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
            ixOmin2,ixOmax1,ixOmax2,idims,Hspeed)
       end if
       ! estimating bounds for the minimum and maximum signal velocities
       if(method==fs_tvdlf.or.method==fs_tvdmu) then
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,ixImax1,&
            ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,Hspeed,cmaxC)
         ! index of var  velocity appears in the induction eq. 
         if(stagger_grid) call phys_get_ct_velocity(vcts,wLp,wRp,ixImin1,&
            ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,&
            cmaxC(ixImin1:ixImax1,ixImin2:ixImax2,index_v_mag))
       else
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,ixImax1,&
            ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,Hspeed,cmaxC,cminC)
         if(stagger_grid) call phys_get_ct_velocity(vcts,wLp,wRp,ixImin1,&
            ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,&
            cmaxC(ixImin1:ixImax1,ixImin2:ixImax2,index_v_mag),&
            cminC(ixImin1:ixImax1,ixImin2:ixImax2,index_v_mag))
       end if

       ! use approximate Riemann solver to get flux at interfaces
       select case(method)
       case(fs_hll)
         do ii=1,number_species
           call get_Riemann_flux_hll(start_indices(ii),stop_indices(ii))
         end do
       case(fs_hllc,fs_hllcd)
         do ii=1,number_species
           call get_Riemann_flux_hllc(start_indices(ii),stop_indices(ii))
         end do
       case(fs_hlld)
         do ii=1,number_species
           if(ii==index_v_mag) then
             call get_Riemann_flux_hlld(start_indices(ii),stop_indices(ii))
           else
             call get_Riemann_flux_hll(start_indices(ii),stop_indices(ii))
           endif   
         end do
       case(fs_tvdlf)
         do ii=1,number_species
           call get_Riemann_flux_tvdlf(start_indices(ii),stop_indices(ii))
         end do
       case(fs_tvdmu)
         call get_Riemann_flux_tvdmu()
       case default
         call mpistop('unkown Riemann flux in finite volume')
       end select

    end do ! Next idims
    b0i=0

    if(stagger_grid) call phys_update_faces(ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,qt,qdt,wprim,fC,fE,sCT,snew,vcts)

    if(slab_uniform) then
      dxinv=-qdt/dxs
      do idims= idimsmin,idimsmax
        hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
        hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);

        ! Multiply the fluxes by -dt/dx since Flux fixing expects this
        fC(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,&
           idims)=dxinv(idims)*fC(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,&
           idims)

        wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           iwstart:nwflux)=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           iwstart:nwflux)+(fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iwstart:nwflux,&
           idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iwstart:nwflux,idims))

        ! For the MUSCL scheme apply the characteristic based limiter
        if(method==fs_tvdmu) call tvdlimit2(method,qdt,ixImin1,ixImin2,ixImax1,&
           ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixOmin1,ixOmin2,ixOmax1,&
           ixOmax2,idims,wLC,wRC,wnew,x,fC,dxs)

      end do ! Next idims
    else
      inv_volume = 1.d0/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      do idims= idimsmin,idimsmax
         hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
         hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);

         if(.not. angmomfix) then ! default case
           do iw=iwstart,nwflux
             fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,&
                idims)=-qdt*fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,&
                idims)*block%surfaceC(ixImin1:ixImax1,ixImin2:ixImax2,idims)
             wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
                idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
                idims)) * inv_volume
           end do
         else
           ! If angular momentum conserving way to solve the equations,
           ! some fluxes additions need to be treated specifically
           call phys_angmomfix(fC,x,wnew,ixImin1,ixImin2,ixImax1,ixImax2,&
              ixOmin1,ixOmin2,ixOmax1,ixOmax2,idims)
         end if

         ! For the MUSCL scheme apply the characteristic based limiter
         if (method==fs_tvdmu) call tvdlimit2(method,qdt,ixImin1,ixImin2,&
            ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixOmin1,ixOmin2,&
            ixOmax1,ixOmax2,idims,wLC,wRC,wnew,x,fC,dxs)

      end do ! Next idims
    end if

    if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,wnew,x)

    if(stagger_grid) call phys_face_to_center(ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       snew)

    ! check and optionally correct unphysical values
    if(fix_small_values) then
       call phys_handle_small_values(.false.,wnew,x,ixImin1,ixImin2,ixImax1,&
          ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,'multi-D finite_volume')
    end if
 
    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,wnew,x,&
       .false.,active,wprim)

    if(phys_solve_eaux.and.levmin==levmax) then
      ! synchronize internal energy for uniform grid
      call phys_energy_synchro(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,wnew,x)
    end if

  end associate
  contains

    subroutine get_Riemann_flux_tvdmu()
      do iw=iwstart,nwflux
         ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=half*(fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2, iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw))
         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2, iw)
      end do
    end subroutine get_Riemann_flux_tvdmu

    subroutine get_Riemann_flux_tvdlf(iws,iwe)
      integer, intent(in) :: iws,iwe
      double precision :: fac(ixCmin1:ixCmax1,ixCmin2:ixCmax2)

      fac(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = &
         -0.5d0*tvdlfeps*cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ii)
      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=iws,iwe
         ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=0.5d0*(fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2, iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw))
         ! Add TVDLF dissipation to the flux
         if (flux_type(idims, iw) /= flux_no_dissipation) then
            fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=fLC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2, iw) + fac(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
         end if
         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2, iw)
      end do ! Next iw

    end subroutine get_Riemann_flux_tvdlf

    subroutine get_Riemann_flux_hll(iws,iwe)
      integer, intent(in) :: iws,iwe
      integer :: ix1,ix2

      do iw=iws,iwe
        if(flux_type(idims, iw) == flux_tvdlf) then
          ! CT MHD does not need normal B flux
          if(stagger_grid) cycle
          fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,&
             idims) = -tvdlfeps*half*max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ii),dabs(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ii))) * (wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
        else
         do ix2=ixCmin2,ixCmax2
         do ix1=ixCmin1,ixCmax1
           if(cminC(ix1,ix2,ii) >= zero) then
             fC(ix1,ix2,iw,idims)=fLC(ix1,ix2,iw)
           else if(cmaxC(ix1,ix2,ii) <= zero) then
             fC(ix1,ix2,iw,idims)=fRC(ix1,ix2,iw)
           else
             ! Add hll dissipation to the flux
             fC(ix1,ix2,iw,idims)=(cmaxC(ix1,ix2,ii)*fLC(ix1,ix2,&
                 iw)-cminC(ix1,ix2,ii)*fRC(ix1,ix2,iw)+cminC(ix1,ix2,&
                ii)*cmaxC(ix1,ix2,ii)*(wRC(ix1,ix2,iw)-wLC(ix1,ix2,&
                iw)))/(cmaxC(ix1,ix2,ii)-cminC(ix1,ix2,ii))
           end if
         end do
         end do
       endif 
      end do

    end subroutine get_Riemann_flux_hll

    subroutine get_Riemann_flux_hllc(iws,iwe)
      integer, intent(in) :: iws, iwe  
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:nwflux)     :: whll, Fhll, fCD
      double precision, dimension(ixImin1:ixImax1,&
         ixImin2:ixImax2)              :: lambdaCD

      integer  :: rho_, p_, e_, eaux_, mom(1:ndir)

      rho_ = iw_rho
      if (allocated(iw_mom)) mom(:) = iw_mom(:)
      e_ = iw_e 
      eaux_ = iw_eaux

      if(associated(phys_hllc_init_species)) then
       call phys_hllc_init_species(ii, rho_, mom(:), e_, eaux_)
      endif  

      p_ = e_

      patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  1
      where(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1) >= zero)
         patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = -2
      elsewhere(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1) <= zero)
         patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  2
      endwhere
      ! Use more diffusive scheme, is actually TVDLF and selected by patchf=4
      if(method==fs_hllcd) call phys_diffuse_hllcd(ixImin1,ixImin2,ixImax1,&
         ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,wLC,wRC,fLC,fRC,patchf)

      !---- calculate speed lambda at CD ----!
      if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==1)) call &
         phys_get_lCD(wLC,wRC,fLC,fRC,cminC(ixImin1:ixImax1,ixImin2:ixImax2,&
         ii),cmaxC(ixImin1:ixImax1,ixImin2:ixImax2,ii),idims,ixImin1,ixImin2,&
         ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2, whll,Fhll,lambdaCD,&
         patchf)

      ! now patchf may be -1 or 1 due to phys_get_lCD
      if(any(abs(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2))== 1))then
         !======== flux at intermediate state ========!
         call phys_get_wCD(wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,&
            cminC(ixImin1:ixImax1,ixImin2:ixImax2,ii),cmaxC(ixImin1:ixImax1,&
            ixImin2:ixImax2,ii),ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
            ixCmin2,ixCmax1,ixCmax2,idims,fCD)
      endif ! Calculate the CD flux

      ! use hll flux for the auxiliary internal e
      if(phys_energy.and.phys_solve_eaux .and. eaux_>0) then
        iw=eaux_
        fCD(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw) = (cmaxC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ii)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
            iw)-cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ii) * fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
            iw) +cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ii)*cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ii)*(wRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           iw)))/(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ii)-cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ii))
      end if

      do iw=iws,iwe
         if (flux_type(idims, iw) == flux_tvdlf) then
            fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw)=-tvdlfeps*half*max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ii),abs(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ii))) * (wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw) - wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
         else
            where(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==-2)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fLC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,iw)
            elsewhere(abs(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2))==1)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fCD(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==2)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fRC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==3)
               ! fallback option, reducing to HLL flux
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=Fhll(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==4)
               ! fallback option, reducing to TVDLF flux
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw) = half*((fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw)) -tvdlfeps * max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ii), dabs(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ii))) * (wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)))
            endwhere
         end if

         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2,iw)

      end do ! Next iw
    end subroutine get_Riemann_flux_hllc

    !> HLLD Riemann flux from Miyoshi 2005 JCP, 208, 315 and Guo 2016 JCP, 327, 543
    subroutine get_Riemann_flux_hlld(iws,iwe)
      integer, intent(in) :: iws, iwe
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:nwflux) :: w1R,w1L,f1R,f1L,f2R,f2L
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:nwflux) :: w2R,w2L
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: sm,s1R,&
         s1L,suR,suL,Bx
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: pts,ptR,&
         ptL,signBx,r1L,r1R,tmp
      ! velocity from the right and the left reconstruction
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ndir) :: vRC,&
          vLC
      ! magnetic field from the right and the left reconstruction
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ndir) :: BR,&
          BL
      integer :: ip1,ip2,ip3,idir,ix1,ix2
      integer  :: rho_, p_, e_, eaux_, mom(1:ndir), mag(1:ndir)

      associate (sR=>cmaxC,sL=>cminC)

      rho_ = iw_rho
      mom(:) = iw_mom(:)
      mag(:) = iw_mag(:) 
      e_ = iw_e 
      eaux_ = iw_eaux 

      p_ = e_

      f1R=0.d0
      f1L=0.d0
      f2R=0.d0
      f2L=0.d0
      w1L=0.d0
      w1R=0.d0
      w2L=0.d0
      w2R=0.d0
      ip1=idims
      ip3=3
      vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wRp(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mom(:))
      vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wLp(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mom(:))
      if(B0field) then
        BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))+block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,&
           ip1)
        BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))+block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,&
           ip1)
      else
        BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))
        BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))
      end if
      if(stagger_grid) then
        Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=block%ws(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip1)
      else
        ! HLL estimation of normal magnetic field at cell interfaces
        ! Li, Shenghai, 2005 JCP, 203, 344, equation (33)
        Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(sR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ii)*BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1)-sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ii)*BL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip1))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ii)-sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ii))
      end if
      ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=wRp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         p_)+0.5d0*sum(BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)**2,dim=ndim+1)
      ptL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=wLp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         p_)+0.5d0*sum(BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)**2,dim=ndim+1)
      if(iw_equi_rho>0) then
        suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = wRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,rho_)+ block%equi_vars(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw_equi_rho,ip1)
      else
        suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = wRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,rho_)
      endif
      suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ii)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ip1))*suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      if(iw_equi_rho>0) then
        suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = wLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,rho_)+ block%equi_vars(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw_equi_rho,ip1)
      else
        suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = wLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,rho_)
      endif
      suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ii)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ip1))*suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      ! Miyoshi equation (38) and Guo euqation (20)
      sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1)-suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*vLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip1)-ptR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)+ptL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))/(suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      ! Miyoshi equation (39) and Guo euqation (28)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip1))=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip1))=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip1))=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip1))=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      ! Guo equation (22)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip1))=Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip1))=Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      if(B0field) then
        ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=wRp(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,p_)+0.5d0*sum(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(:))**2,dim=ndim+1)
        ptL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=wLp(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,p_)+0.5d0*sum(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(:))**2,dim=ndim+1)
      end if

      ! Miyoshi equation (43) and Guo equation (27)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)=suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)=suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)/(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))

      ip2=mod(ip1+1,ndir)
      if(ip2==0) ip2=ndir
      r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))-Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)**2
      where(abs(r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2))>smalldouble)
        r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=1.d0/r1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)
      else where
        r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.d0
      end where
      r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))-Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)**2
      where(abs(r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2))>smalldouble)
        r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=1.d0/r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)
      else where
        r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.d0
      end where
      ! Miyoshi equation (44)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip2))=vRC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip2)-Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip2)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip1))*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip2))=vLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip2)-Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip2)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip1))*r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      ! partial solution for later usage
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=(suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ii)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ip1))-Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)**2)*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=(suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ii)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ip1))-Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)**2)*r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      if(ndir==3) then
        ip3=mod(ip1+2,ndir)
        if(ip3==0) ip3=ndir
        ! Miyoshi equation (46)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip3))=vRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip3)-Bx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip3)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip1))*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip3))=vLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip3)-Bx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip3)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip1))*r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        ! Miyoshi equation (47)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip3))=BR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip3)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip3))=BL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip3)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))
      end if
      ! Miyoshi equation (45)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=BR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip2)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=BL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip2)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))
      if(B0field) then
        ! Guo equation (26)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:))=w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))-block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,&
           ip1)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:))=w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))-block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,&
           ip1)
      end if
      ! equation (48)
      if(phys_energy) then
        ! Guo equation (25) equivalent to Miyoshi equation (41)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)=suR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sm(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1))+ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        !w1L(ixC^S,p_)=suL(ixC^S)*(sm(ixC^S)-vLC(ixC^S,ip1))+ptL(ixC^S)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)=w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,p_)
        if(B0field) then
          ! Guo equation (32)
          w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)=w1R(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,p_)+sum(block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             :,ip1)*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))-w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:))),dim=ndim+1)
          w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)=w1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,p_)+sum(block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             :,ip1)*(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:))),dim=ndim+1)
        end if
        ! Miyoshi equation (48) and main part of Guo euqation (31)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=((sR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ii)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1))*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)-ptR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1)+w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)*sm(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)+Bx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sum(vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           :)*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)-sum(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=((sL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ii)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1))*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)-ptL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1)+w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)*sm(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)+Bx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sum(vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           :)*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)-sum(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)))/(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
        if(B0field) then
          ! Guo equation (31)
          w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=w1R(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,e_)+(sum(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,ip1),&
             dim=ndim+1)*sm(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)-sum(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,ip1),&
             dim=ndim+1)*vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ip1))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ii)-sm(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2))
          w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=w1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,e_)+(sum(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,ip1),&
             dim=ndim+1)*sm(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)-sum(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,ip1),&
             dim=ndim+1)*vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ip1))/(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ii)-sm(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2))
        end if
        if(iw_equi_p>0) then
#if !defined(E_RM_W0) || E_RM_W0 == 1
          w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)= w1R(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,e_) + 1d0/(phys_gamma - 1) * &
             block%equi_vars(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw_equi_p,&
             ip1) * (sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ip1))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
          w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)= w1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,e_) + 1d0/(phys_gamma - 1) * &
             block%equi_vars(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw_equi_p,&
             ip1) * (sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vLC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ip1))/(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
#else
          w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)= w1R(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,e_) + phys_gamma /(phys_gamma - 1) * &
             block%equi_vars(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw_equi_p,&
             ip1) * (sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ip1))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
          w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)= w1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,e_) + phys_gamma /(phys_gamma - 1) * &
             block%equi_vars(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw_equi_p,&
             ip1) * (sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vLC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ip1))/(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
#endif
        endif
      end if

      ! Miyoshi equation (49) and Guo equation (35)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)=w1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,rho_)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)=w1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,rho_)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip1))=w1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mag(ip1))
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip1))=w1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mag(ip1))

      r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sqrt(w1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,rho_))
      r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sqrt(w1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,rho_))
      tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=1.d0/(r1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)+r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sign(1.d0,Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))
      ! Miyoshi equation (51) and Guo equation (33)
      s1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)+abs(Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))/r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      s1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-abs(Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))/r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      ! Miyoshi equation (59) and Guo equation (41)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip2))=(r1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mom(ip2))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*w1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mom(ip2))+(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip2))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip2)))*signBx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip2))=w2R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mom(ip2))
      ! Miyoshi equation (61) and Guo equation (43)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=(r1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip2))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*w1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mag(ip2))+r1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*r1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mom(ip2))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mom(ip2)))*signBx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=w2R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mag(ip2))
      if(ndir==3) then
        ! Miyoshi equation (60) and Guo equation (42)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip3))=(r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(ip3))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(ip3))+(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(ip3))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(ip3)))*signBx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip3))=w2R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(ip3))
        ! Miyoshi equation (62) and Guo equation (44)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip3))=(r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(ip3))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(ip3))+r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*r1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(ip3))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(ip3)))*signBx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip3))=w2R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(ip3))
      end if
      ! Miyoshi equation (63) and Guo equation (45)
      if(phys_energy) then
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,e_)+r1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sum(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)-sum(w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,e_)-r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sum(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)-sum(w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      end if

      ! convert velocity to momentum
      do idir=1,ndir
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(idir))=w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(idir))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           rho_)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(idir))=w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(idir))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           rho_)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(idir))=w2R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(idir))*w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           rho_)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(idir))=w2L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(idir))*w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           rho_)
      end do
      if(iw_equi_rho>0) then
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_) = w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,rho_) - block%equi_vars(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw_equi_rho,ip1)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_) = w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,rho_) - block%equi_vars(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw_equi_rho,ip1)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_) = w2R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,rho_) - block%equi_vars(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw_equi_rho,ip1)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_) = w2L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,rho_) - block%equi_vars(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw_equi_rho,ip1)
      endif
      ! get fluxes of intermedate states
      do iw=iws,iwe
        if(flux_type(idims, iw) == flux_special) then
          ! known flux (fLC=fRC) for normal B and psi_ in GLM method
          f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fLC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
          f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
          f2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
          f2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
        else if(flux_type(idims, iw) == flux_hll) then
          ! using hll flux for eaux and tracers
          f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=(sR(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ii)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              iw)-sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ii)*fRC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2, iw) +sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ii)*sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ii)*(wRC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iw)))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ii)-sL(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ii))
          f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
          f2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
          f2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
        else
          ! construct hlld flux
          f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fLC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)+sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ii)*(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)-wLC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw))
          f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fRC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)+sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ii)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)-wRC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw))
          f2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)+s1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)*(w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iw)-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
          f2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1R(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)+s1R(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)*(w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iw)-w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
        end if
      end do

      ! Miyoshi equation (66) and Guo equation (46)
     do ix2=ixCmin2,ixCmax2
     do ix1=ixCmin1,ixCmax1
        if(sL(ix1,ix2,ii)>0.d0) then
          fC(ix1,ix2,iws:iwe,ip1)=fLC(ix1,ix2,iws:iwe)
        else if(s1L(ix1,ix2)>=0.d0) then
          fC(ix1,ix2,iws:iwe,ip1)=f1L(ix1,ix2,iws:iwe)
        else if(sm(ix1,ix2)>=0.d0) then
          fC(ix1,ix2,iws:iwe,ip1)=f2L(ix1,ix2,iws:iwe)
        else if(s1R(ix1,ix2)>=0.d0) then
          fC(ix1,ix2,iws:iwe,ip1)=f2R(ix1,ix2,iws:iwe)
        else if(sR(ix1,ix2,ii)>=0.d0) then
          fC(ix1,ix2,iws:iwe,ip1)=f1R(ix1,ix2,iws:iwe)
        else if(sR(ix1,ix2,ii)<0.d0) then
          fC(ix1,ix2,iws:iwe,ip1)=fRC(ix1,ix2,iws:iwe)
        end if
     end do
     end do

      end associate
    end subroutine get_Riemann_flux_hlld

    !> HLLD Riemann flux from Miyoshi 2005 JCP, 208, 315 and Guo 2016 JCP, 327, 543
    !> https://arxiv.org/pdf/2108.04991.pdf
    subroutine get_Riemann_flux_hlld_mag2()
      !use mod_mhd_phys
      use mod_variables
      use mod_physics
      implicit none
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:nwflux) :: w1R,w1L,f1R,f1L,f2R,f2L
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:nwflux) :: w2R,w2L
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: sm,s1R,&
         s1L,suR,suL,Bx
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: pts,ptR,&
         ptL,signBx,r1L,r1R,tmp
      ! velocity from the right and the left reconstruction
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ndir) :: vRC,&
          vLC
      ! magnetic field from the right and the left reconstruction
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ndir) :: BR,&
          BL
      integer :: ip1,ip2,ip3,idir,ix1,ix2
      double precision :: phiPres, thetaSM, du, dv, dw
      integer :: ixVmin1,ixVmin2,ixVmax1,ixVmax2, ixVbmin1,ixVbmin2,ixVbmax1,&
         ixVbmax2, ixVcmin1,ixVcmin2,ixVcmax1,ixVcmax2, ixVdmin1,ixVdmin2,&
         ixVdmax1,ixVdmax2, ixVemin1,ixVemin2,ixVemax1,ixVemax2, ixVfmin1,&
         ixVfmin2,ixVfmax1,ixVfmax2
      integer  :: rho_, p_, e_, eaux_, mom(1:ndir), mag(1:ndir)
      double precision, parameter :: aParam = 4d0

      rho_ = iw_rho
      mom(:) = iw_mom(:)
      mag(:) = iw_mag(:) 
      p_ = iw_e
      e_ = iw_e 
      eaux_ = iw_eaux 

      associate (sR=>cmaxC,sL=>cminC)

      f1R=0.d0
      f1L=0.d0
      ip1=idims
      ip3=3

      vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wRp(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mom(:))
      vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wLp(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mom(:))

      ! reuse s1L s1R
      call get_hlld2_modif_c(wLp,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,s1L)
      call get_hlld2_modif_c(wRp,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,s1R)
      !phiPres = min(1, maxval(max(s1L(ixO^S),s1R(ixO^S))/cmaxC(ixO^S,1)))
      phiPres = min(1d0, maxval(max(s1L(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
         s1R(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))/maxval(cmaxC(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1)))
      phiPres = phiPres*(2D0 - phiPres)  

     !we use here not reconstructed velocity: wprim?       
     ixVmin1=ixOmin1;ixVmin2=ixOmin2;ixVmax1=ixOmax1;ixVmax2=ixOmax2;
     !first dim
     ixVmin1=ixOmin1+1  
     ixVmax1=ixOmax1+1
     du = minval(wprim(ixVmin1:ixVmax1,ixVmin2:ixVmax2,&
        mom(1))-wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1)))
     if(du>0d0) du=0d0  
     dv = 0d0 
     dw = 0d0 

     
     !second dim
     !i,j-1,k 
     ixVmin1=ixOmin1;ixVmin2=ixOmin2;ixVmax1=ixOmax1;ixVmax2=ixOmax2;
     ixVmin2=ixOmin2-1  
     ixVmax2=ixOmax2-1
 
     !i,j+1,k 
     ixVbmin1=ixOmin1;ixVbmin2=ixOmin2;ixVbmax1=ixOmax1;ixVbmax2=ixOmax2;
     ixVbmin2=ixOmin2+1  
     ixVbmax2=ixOmax2+1

     !i+1,j,k 
     ixVcmin1=ixOmin1;ixVcmin2=ixOmin2;ixVcmax1=ixOmax1;ixVcmax2=ixOmax2;
     ixVcmin1=ixOmin1+1  
     ixVcmax1=ixOmax1+1

     !i+1,j-1,k 
     ixVdmin1=ixOmin1;ixVdmin2=ixOmin2;ixVdmax1=ixOmax1;ixVdmax2=ixOmax2;
     ixVdmin1=ixOmin1+1  
     ixVdmax1=ixOmax1+1
     ixVdmin2=ixOmin2-1  
     ixVdmax2=ixOmax2-1

     !i+1,j+1,k 
     ixVemin1=ixOmin1;ixVemin2=ixOmin2;ixVemax1=ixOmax1;ixVemax2=ixOmax2;
     ixVemin1=ixOmin1+1  
     ixVemax1=ixOmax1+1
     ixVemin2=ixOmin2+1  
     ixVemax2=ixOmax2+1

     dv = minval(min(wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        mom(2))-wprim(ixVmin1:ixVmax1,ixVmin2:ixVmax2,mom(2)),&
        wprim(ixVbmin1:ixVbmax1,ixVbmin2:ixVbmax2,&
        mom(2))-wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2)),&
        wprim(ixVcmin1:ixVcmax1,ixVcmin2:ixVcmax2,&
        mom(2))-wprim(ixVdmin1:ixVdmax1,ixVdmin2:ixVdmax2,mom(2)),&
        wprim(ixVemin1:ixVemax1,ixVemin2:ixVemax2,&
        mom(2))-wprim(ixVcmin1:ixVcmax1,ixVcmin2:ixVcmax2,mom(2))))
     if(dv>0d0) dv=0d0 

      
     thetaSM = maxval(cmaxC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)) 
       
     thetaSM = (min(1d0, (thetaSM-du)/(thetaSM-min(dv,dw))))**aParam
     !print*, "HLLD2 ", du,dv,dw, thetaSM, phiPres 

      if(B0field) then
        BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))+block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,&
           ip1)
        BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))+block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,&
           ip1)
      else
        BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))
        BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))
      end if
      ! HLL estimation of normal magnetic field at cell interfaces
      Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         index_v_mag)*BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1)-sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         index_v_mag)*BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1)-fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip1))-fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip1)))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         index_v_mag)-sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,index_v_mag))
      ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=wRp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         p_)+0.5d0*sum(BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)**2,dim=ndim+1)
      ptL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=wLp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         p_)+0.5d0*sum(BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)**2,dim=ndim+1)
      suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         index_v_mag)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1))*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)
      suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         index_v_mag)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1))*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)
      ! Miyoshi equation (38) and Guo euqation (20)
      sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1)-suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*vLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip1)-thetaSM*(ptR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-ptL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)) )/(suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      ! Miyoshi equation (39) and Guo euqation (28)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip1))=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip1))=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip1))=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip1))=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      ! Guo equation (22)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip1))=Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip1))=Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      if(B0field) then
        ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=wRp(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,p_)+0.5d0*sum(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(:))**2,dim=ndim+1)
        ptL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=wLp(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,p_)+0.5d0*sum(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(:))**2,dim=ndim+1)
      end if

      ! Miyoshi equation (43) and Guo equation (27)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)=suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         index_v_mag)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)=suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)/(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         index_v_mag)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))

      ip2=mod(ip1+1,ndir)
      if(ip2==0) ip2=ndir
      r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         index_v_mag)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))-Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)**2
      where(r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/=0.d0)
        r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=1.d0/r1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)
      endwhere
      r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         index_v_mag)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))-Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)**2
      where(r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/=0.d0)
        r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=1.d0/r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)
      endwhere
      ! Miyoshi equation (44)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip2))=vRC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip2)-Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip2)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip1))*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip2))=vLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip2)-Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip2)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip1))*r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      ! partial solution for later usage
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=(suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         index_v_mag)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)**2)*r1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=(suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         index_v_mag)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)**2)*r1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      if(ndir==3) then
        ip3=mod(ip1+2,ndir)
        if(ip3==0) ip3=ndir
        ! Miyoshi equation (46)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip3))=vRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip3)-Bx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip3)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip1))*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip3))=vLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip3)-Bx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip3)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip1))*r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        ! Miyoshi equation (47)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip3))=BR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip3)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip3))=BL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip3)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))
      end if
      ! Miyoshi equation (45)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=BR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip2)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=BL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip2)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))
      if(B0field) then
        ! Guo equation (26)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:))=w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))-block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,&
           ip1)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:))=w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))-block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,&
           ip1)
      end if
      ! equation (48)
      if(phys_energy) then
        ! Guo equation (25) equivalent to Miyoshi equation (41)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)=(suR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*ptL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2) - suL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*ptR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2) +phiPres * suR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*suL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1)))/(suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-suL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2))
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)=w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,p_)
        !if(mhd_solve_eaux) then
        !  w1R(ixC^S,eaux_)=(w1R(ixC^S,p_)-half*sum(w1R(ixC^S,mag(:))**2,dim=ndim+1))/(mhd_gamma-one)
        !  w1L(ixC^S,eaux_)=(w1L(ixC^S,p_)-half*sum(w1L(ixC^S,mag(:))**2,dim=ndim+1))/(mhd_gamma-one)
        !end if
        if(B0field) then
          ! Guo equation (32)
          w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)=w1R(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,p_)+sum(block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             :,ip1)*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))-w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:))),dim=ndim+1)
          w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)=w1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,p_)+sum(block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             :,ip1)*(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:))),dim=ndim+1)
        end if
        ! Miyoshi equation (48) and main part of Guo euqation (31)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=((sR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,index_v_mag)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1))*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)-ptR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1)+w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)*sm(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)+Bx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sum(vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           :)*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)-sum(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           index_v_mag)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=((sL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,index_v_mag)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1))*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)-ptL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1)+w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)*sm(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)+Bx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sum(vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           :)*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)-sum(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)))/(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           index_v_mag)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
        if(B0field) then
          ! Guo equation (31)
          w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=w1R(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,e_)+(sum(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,ip1),&
             dim=ndim+1)*sm(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)-sum(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,ip1),&
             dim=ndim+1)*vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ip1))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             index_v_mag)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
          w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=w1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,e_)+(sum(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,ip1),&
             dim=ndim+1)*sm(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)-sum(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,ip1),&
             dim=ndim+1)*vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ip1))/(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             index_v_mag)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
        end if
      end if

      ! Miyoshi equation (49) and Guo equation (35)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)=w1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,rho_)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)=w1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,rho_)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip1))=w1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mag(ip1))
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip1))=w1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mag(ip1))

      r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sqrt(w1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,rho_))
      r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sqrt(w1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,rho_))
      tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=1.d0/(r1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)+r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sign(1.d0,Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))
      ! Miyoshi equation (51) and Guo equation (33)
      s1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)+abs(Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))/r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      s1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-abs(Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))/r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      ! Miyoshi equation (59) and Guo equation (41)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip2))=(r1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mom(ip2))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*w1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mom(ip2))+(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip2))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip2)))*signBx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip2))=w2R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mom(ip2))
      ! Miyoshi equation (61) and Guo equation (43)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=(r1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip2))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*w1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mag(ip2))+r1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*r1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mom(ip2))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mom(ip2)))*signBx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=w2R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mag(ip2))
      if(ndir==3) then
        ! Miyoshi equation (60) and Guo equation (42)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip3))=(r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(ip3))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(ip3))+(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(ip3))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(ip3)))*signBx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip3))=w2R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(ip3))
        ! Miyoshi equation (62) and Guo equation (44)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip3))=(r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(ip3))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(ip3))+r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*r1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(ip3))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(ip3)))*signBx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip3))=w2R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(ip3))
      end if
      ! Miyoshi equation (63) and Guo equation (45)
      if(phys_energy) then
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,e_)+r1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sum(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)-sum(w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,e_)-r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sum(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)-sum(w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      end if

      ! convert velocity to momentum
      do idir=1,ndir
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(idir))=w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(idir))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           rho_)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(idir))=w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(idir))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           rho_)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(idir))=w2R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(idir))*w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           rho_)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(idir))=w2L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(idir))*w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           rho_)
      end do

      ! get fluxes of intermedate states
      do iw=1,nwflux
        ! CT MHD does not need normal B flux
        if(stagger_grid .and. flux_type(idims, iw) == flux_tvdlf) cycle
        if(flux_type(idims, iw) == flux_special) then
          ! known flux (fLC=fRC) for normal B and psi_ in GLM method
          f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fLC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
          f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
          f2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
          f2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
        else if(flux_type(idims, iw) == flux_hll) then
          ! using hll flux for eaux and tracers
          f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=(sR(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,index_v_mag)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              iw)-sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             index_v_mag)*fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              iw) +sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             index_v_mag)*sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             index_v_mag)*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)))/(sR(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,index_v_mag)-sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             index_v_mag))
          f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
          f2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
          f2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
        else
          f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fLC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)+sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             index_v_mag)*(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
          f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fRC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)+sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             index_v_mag)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iw)-wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
          f2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)+s1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)*(w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iw)-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
          f2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1R(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)+s1R(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)*(w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iw)-w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
        end if
      end do

      ! Miyoshi equation (66) and Guo equation (46)
     do ix2=ixCmin2,ixCmax2
     do ix1=ixCmin1,ixCmax1
        if(sL(ix1,ix2,index_v_mag)>0.d0) then
          fC(ix1,ix2,1:nwflux,ip1)=fLC(ix1,ix2,1:nwflux)
        else if(s1L(ix1,ix2)>=0.d0) then
          fC(ix1,ix2,1:nwflux,ip1)=f1L(ix1,ix2,1:nwflux)
        else if(sm(ix1,ix2)>=0.d0) then
          fC(ix1,ix2,1:nwflux,ip1)=f2L(ix1,ix2,1:nwflux)
        else if(s1R(ix1,ix2)>=0.d0) then
          fC(ix1,ix2,1:nwflux,ip1)=f2R(ix1,ix2,1:nwflux)
        else if(sR(ix1,ix2,index_v_mag)>=0.d0) then
          fC(ix1,ix2,1:nwflux,ip1)=f1R(ix1,ix2,1:nwflux)
        else if(sR(ix1,ix2,index_v_mag)<0.d0) then
          fC(ix1,ix2,1:nwflux,ip1)=fRC(ix1,ix2,1:nwflux)
        end if
     end do
     end do

      end associate
    end subroutine get_Riemann_flux_hlld_mag2

    !> Calculate fast magnetosonic wave speed
    subroutine get_hlld2_modif_c(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,csound)
      use mod_global_parameters

      integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2
      double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw),&
          x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
      double precision, intent(out):: csound(ixImin1:ixImax1,ixImin2:ixImax2)
      double precision :: cfast2(ixImin1:ixImax1,ixImin2:ixImax2),&
          AvMinCs2(ixImin1:ixImax1,ixImin2:ixImax2), b2(ixImin1:ixImax1,&
         ixImin2:ixImax2), kmax
      double precision :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
          gamma_A2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      integer  :: rho_, p_, e_, eaux_, mom(1:ndir), mag(1:ndir)

        rho_ = iw_rho
        mom(:) = iw_mom(:)
        mag(:) = iw_mag(:) 
        p_ = iw_e
        e_ = iw_e 
        eaux_ = iw_eaux 

      inv_rho=1.d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

      ! store |B|^2 in v

      if (B0field) then
        b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = sum((w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2, mag(:))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:,&
           b0i))**2, dim=ndim+1)
      else
        b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = sum(w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2, mag(:))**2, dim=ndim+1)
      end if


      if (B0field) then
        AvMinCs2= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(idims))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims,b0i)
      else
        AvMinCs2= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mag(idims))
      end if


      csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = sum(w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, mom(:))**2, dim=ndim+1)

      cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)   = b2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2) * inv_rho+csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cfast2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)**2-4.0d0*csound(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2) * AvMinCs2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)**2 * inv_rho 

      where(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<zero)
         AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
      end where

      AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(AvMinCs2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))

     csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
        sqrt(half*(cfast2(ixOmin1:ixOmax1,&
        ixOmin2:ixOmax2)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))

    end subroutine get_hlld2_modif_c

  end subroutine finite_volume

  !> Determine the upwinded wLC(ixL) and wRC(ixR) from w.
  !> the wCT is only used when PPM is exploited.
  subroutine reconstruct_LR(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
     ixLmax1,ixLmax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,idims,w,wLC,wRC,wLp,wRp,x,&
     dxdim)
    use mod_physics
    use mod_global_parameters
    use mod_limiter

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixLmin1,ixLmin2,&
       ixLmax1,ixLmax2, ixRmin1,ixRmin2,ixRmax1,ixRmax2, idims
    double precision, intent(in) :: dxdim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: w
    ! left and right constructed status in conservative form
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLC,&
        wRC
    ! left and right constructed status in primitive form
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLp,&
        wRp
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim) :: x

    integer            :: jxRmin1,jxRmin2,jxRmax1,jxRmax2, ixCmin1,ixCmin2,&
       ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2, iw
    double precision   :: ldw(ixImin1:ixImax1,ixImin2:ixImax2),&
        rdw(ixImin1:ixImax1,ixImin2:ixImax2), dwC(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision   :: a2max

    select case (type_limiter(block%level))
    case (limiter_mp5)
       call MP5limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,ixLmax1,&
          ixLmax2,idims,w,wLp,wRp)
    case (limiter_weno3)
       call WENO3limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp,1)
    case (limiter_wenoyc3)
       call WENO3limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp,2)
    case (limiter_weno5)
       call WENO5limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp,1)
    case (limiter_weno5nm)
       call WENO5NMlimiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp,1)
    case (limiter_wenoz5)
       call WENO5limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp,2)
    case (limiter_wenoz5nm)
       call WENO5NMlimiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp,2)
    case (limiter_wenozp5)
       call WENO5limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp,3)
    case (limiter_wenozp5nm)
       call WENO5NMlimiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp,3)
    case (limiter_weno5cu6)
       call WENO5CU6limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,w,wLp,wRp)
    case (limiter_teno5ad)
       call TENO5ADlimiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp)
    case (limiter_weno7)
       call WENO7limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,w,wLp,wRp,1)
    case (limiter_mpweno7)
       call WENO7limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,w,wLp,wRp,2)
    case (limiter_venk)
       call venklimiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp) 
       if(fix_small_values) then
          call phys_handle_small_values(.true.,wLp,x,ixImin1,ixImin2,ixImax1,&
             ixImax2,ixLmin1,ixLmin2,ixLmax1,ixLmax2,'reconstruct left')
          call phys_handle_small_values(.true.,wRp,x,ixImin1,ixImin2,ixImax1,&
             ixImax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,'reconstruct right')
       end if
    case (limiter_ppm)
       call PPMlimiter(ixImin1,ixImin2,ixImax1,ixImax2,ixMlo1,ixMlo2,ixMhi1,&
          ixMhi2,idims,w,w,wLp,wRp)
       if(fix_small_values) then
          call phys_handle_small_values(.true.,wLp,x,ixImin1,ixImin2,ixImax1,&
             ixImax2,ixLmin1,ixLmin2,ixLmax1,ixLmax2,'reconstruct left')
          call phys_handle_small_values(.true.,wRp,x,ixImin1,ixImin2,ixImax1,&
             ixImax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,'reconstruct right')
       end if
    case default
       jxRmin1=ixRmin1+kr(idims,1);jxRmin2=ixRmin2+kr(idims,2)
       jxRmax1=ixRmax1+kr(idims,1);jxRmax2=ixRmax2+kr(idims,2);
       ixCmax1=jxRmax1;ixCmax2=jxRmax2; ixCmin1=ixLmin1-kr(idims,1)
       ixCmin2=ixLmin2-kr(idims,2);
       jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
       jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);
       do iw=1,nwflux
          if (loglimit(iw)) then
             w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,iw)=dlog10(w(ixCmin1:jxCmax1,&
                ixCmin2:jxCmax2,iw))
             wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
                iw)=dlog10(wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw))
             wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                iw)=dlog10(wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw))
          end if

          dwC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=w(jxCmin1:jxCmax1,&
             jxCmin2:jxCmax2,iw)-w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)
          if(need_global_a2max) then 
            a2max=a2max_global(idims)
          else
            select case(idims)
            case(1)
              a2max=schmid_rad1
            
            case(2)
              a2max=schmid_rad2
            
            case default
              call mpistop("idims is wrong in mod_limiter")
            end select
          end if
            
          ! limit flux from left and/or right
          call dwlimiter2(dwC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
             ixCmax1,ixCmax2,idims,type_limiter(block%level),ldw,rdw,&
             a2max=a2max)
          wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw)=wLp(ixLmin1:ixLmax1,&
             ixLmin2:ixLmax2,iw)+half*ldw(ixLmin1:ixLmax1,ixLmin2:ixLmax2)
          wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw)=wRp(ixRmin1:ixRmax1,&
             ixRmin2:ixRmax2,iw)-half*rdw(jxRmin1:jxRmax1,jxRmin2:jxRmax2)

          if (loglimit(iw)) then
             w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,iw)=10.0d0**w(ixCmin1:jxCmax1,&
                ixCmin2:jxCmax2,iw)
             wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
                iw)=10.0d0**wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw)
             wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                iw)=10.0d0**wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw)
          end if
       end do
       if(fix_small_values) then
          call phys_handle_small_values(.true.,wLp,x,ixImin1,ixImin2,ixImax1,&
             ixImax2,ixLmin1,ixLmin2,ixLmax1,ixLmax2,'reconstruct left')
          call phys_handle_small_values(.true.,wRp,x,ixImin1,ixImin2,ixImax1,&
             ixImax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,'reconstruct right')
       end if
    end select

    wLC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nw)=wLp(ixLmin1:ixLmax1,&
       ixLmin2:ixLmax2,1:nw)
    wRC(ixRmin1:ixRmax1,ixRmin2:ixRmax2,1:nw)=wRp(ixRmin1:ixRmax1,&
       ixRmin2:ixRmax2,1:nw)
    call phys_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
       ixLmax1,ixLmax2,wLC,x)
    call phys_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixRmin1,ixRmin2,&
       ixRmax1,ixRmax2,wRC,x)

  end subroutine reconstruct_LR

end module mod_finite_volume
