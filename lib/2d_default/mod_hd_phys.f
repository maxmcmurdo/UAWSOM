!> Hydrodynamics physics module
module mod_hd_phys
  use mod_thermal_conduction, only: tc_fluid
  use mod_radiative_cooling, only: rc_fluid
  use mod_thermal_emission, only: te_fluid
  implicit none
  private

  !> Whether an energy equation is used
  logical, public, protected              :: hd_energy = .true.

  !> Whether thermal conduction is added
  logical, public, protected              :: hd_thermal_conduction = .false.
  type(tc_fluid), allocatable, public :: tc_fl
  type(te_fluid), allocatable, public :: te_fl_hd

  !> Whether radiative cooling is added
  logical, public, protected              :: hd_radiative_cooling = .false.
  type(rc_fluid), allocatable, public :: rc_fl

  !> Whether dust is added
  logical, public, protected              :: hd_dust = .false.

  !> Whether viscosity is added
  logical, public, protected              :: hd_viscosity = .false.

  !> Whether gravity is added
  logical, public, protected              :: hd_gravity = .false.

  !> Whether particles module is added
  logical, public, protected              :: hd_particles = .false.

  !> Whether rotating frame is activated
  logical, public, protected              :: hd_rotating_frame = .false.

  !> Whether CAK radiation line force is activated
  logical, public, protected              :: hd_cak_force = .false.

  !> Number of tracer species
  integer, public, protected              :: hd_n_tracer = 0

  !> Index of the density (in the w array)
  integer, public, protected              :: rho_

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public, protected              :: p_

  !> Index of the cutoff temperature for the TRAC method
  integer, public, protected              :: Tcoff_

  !> The adiabatic index
  double precision, public                :: hd_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: hd_adiab = 1.0d0

  !> The small_est allowed energy
  double precision, protected             :: small_e

  !> Whether TRAC method is used
  logical, public, protected              :: hd_trac = .false.
  integer, public, protected              :: hd_trac_type = 1

  !> Allows overruling default corner filling (for debug mode, since otherwise corner primitives fail)
  logical, public, protected              :: hd_force_diagonal = .false.

  !> Helium abundance over Hydrogen
  double precision, public, protected  :: He_abundance=0.1d0

  ! Public methods
  public :: hd_phys_init
  public :: hd_kin_en
  public :: hd_get_pthermal
  public :: hd_get_csound2
  public :: hd_to_conserved
  public :: hd_to_primitive
  public :: hd_check_params
  public :: hd_check_w

contains

  !> Read this module's parameters from a file
  subroutine hd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /hd_list/ hd_energy, hd_n_tracer, hd_gamma, hd_adiab, hd_dust,&
        hd_thermal_conduction, hd_radiative_cooling, hd_viscosity, hd_gravity,&
        He_abundance, SI_unit, hd_particles, hd_rotating_frame, hd_trac,&
        hd_force_diagonal, hd_trac_type, hd_cak_force

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, hd_list, end=111)
111    close(unitpar)
    end do

  end subroutine hd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine hd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = hd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine hd_write_info

  !> Add fluxes in an angular momentum conserving way
  subroutine hd_angmomfix(fC,x,wnew,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim)
    use mod_global_parameters
    use mod_dust, only: dust_n_species, dust_mom
    use mod_geometry
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout)    :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nwflux,1:ndim),  wnew(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    integer, intent(in)                :: idim
    integer                            :: hxOmin1,hxOmin2,hxOmax1,hxOmax2,&
        kxCmin1,kxCmin2,kxCmax1,kxCmax2, iw
    double precision                   :: inv_volume(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    logical isangmom

    ! shifted indexes
    hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
    hxOmax1=ixOmax1-kr(idim,1);hxOmax2=ixOmax2-kr(idim,2);
    ! all the indexes
    kxCmin1=hxOmin1;kxCmin2=hxOmin2;
    kxCmax1=ixOmax1;kxCmax2=ixOmax2;

    inv_volume(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
       1.0d0/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    select case(coordinate)
    case (cylindrical)
       do iw=1,nwflux
        isangmom = (iw==iw_mom(phi_))
        if (hd_dust) isangmom = (isangmom .or. any(dust_mom(phi_,&
           1:dust_n_species) == iw))
        if (idim==r_ .and. isangmom) then
          fC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,iw,idim)= fC(kxCmin1:kxCmax1,&
             kxCmin2:kxCmax2,iw,idim)*(x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,&
             r_)+half*block%dx(kxCmin1:kxCmax1,kxCmin2:kxCmax2,idim))
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
             idim)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
             idim)) * (inv_volume(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim))
        else
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
             idim)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
             idim)) * inv_volume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
        endif
      enddo
     case (spherical)
      if (hd_dust) call mpistop("Error: hd_angmomfix is not implemented &
      
      
        &with dust and coordinate=='spherical'")
      do iw=1,nwflux
        if     (idim==r_ .and. (iw==iw_mom(2) .or. iw==iw_mom(phi_))) then
          fC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,iw,idim)= fC(kxCmin1:kxCmax1,&
             kxCmin2:kxCmax2,iw,idim)*(x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,&
             idim)+half*block%dx(kxCmin1:kxCmax1,kxCmin2:kxCmax2,idim))
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
             idim)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
             idim)) * (inv_volume(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim))
        elseif (idim==2  .and. iw==iw_mom(phi_)) then
          fC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,iw,idim)=fC(kxCmin1:kxCmax1,&
             kxCmin2:kxCmax2,iw,idim)*sin(x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,&
             idim)+half*block%dx(kxCmin1:kxCmax1,kxCmin2:kxCmax2,idim)) !(x(4,3,1)-x(3,3,1)))
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
             idim)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
             idim)) * (inv_volume(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)/sin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim)))
        else
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
             idim)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
             idim)) * inv_volume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
        endif
      enddo

    end select

  end subroutine hd_angmomfix

  !> Initialize the module
  subroutine hd_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_dust, only: dust_init
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_particles, only: particles_init
    use mod_rotating_frame, only:rotating_frame_init
    use mod_cak_force, only: cak_init
    use mod_physics
    use mod_supertimestepping, only: sts_init, add_sts_method,&
       set_conversion_methods_to_head, set_error_handling_to_head

    integer :: itr, idir

    call hd_read_params(par_files)

    physics_type = "hd"
    phys_energy  = hd_energy
    phys_total_energy  = hd_energy
    phys_gamma = hd_gamma

    phys_trac=hd_trac
    if(phys_trac) then
      if(ndim .eq. 1) then
        if(hd_trac_type .gt. 2) then
          hd_trac_type=1
          if(mype==0) write(*,*) 'WARNING: set hd_trac_type=1'
        end if
        phys_trac_type=hd_trac_type
      else
        phys_trac=.false.
        if(mype==0) write(*,*) 'WARNING: set hd_trac=F when ndim>=2'
      end if
    end if

    ! set default gamma for polytropic/isothermal process
    if(.not.hd_energy) then
      if(hd_thermal_conduction) then
        hd_thermal_conduction=.false.
        if(mype==0) write(*,*)&
            'WARNING: set hd_thermal_conduction=F when hd_energy=F'
      end if
      if(hd_radiative_cooling) then
        hd_radiative_cooling=.false.
        if(mype==0) write(*,*)&
            'WARNING: set hd_radiative_cooling=F when hd_energy=F'
      end if
    end if
    use_particles = hd_particles

    allocate(start_indices(number_species),stop_indices(number_species))

    ! set the index of the first flux variable for species 1
    start_indices(1)=1
    ! Determine flux variables
    rho_ = var_set_rho()

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    if (hd_energy) then
       e_ = var_set_energy()
       p_ = e_
    else
       e_ = -1
       p_ = -1
    end if

    phys_get_dt              => hd_get_dt
    phys_get_cmax            => hd_get_cmax
    phys_get_a2max           => hd_get_a2max
    phys_get_tcutoff         => hd_get_tcutoff
    phys_get_cbounds         => hd_get_cbounds
    phys_get_flux            => hd_get_flux
    phys_add_source_geom     => hd_add_source_geom
    phys_add_source          => hd_add_source
    phys_to_conserved        => hd_to_conserved
    phys_to_primitive        => hd_to_primitive
    !phys_ei_to_e             => hd_ei_to_e
    !phys_e_to_ei             => hd_e_to_ei
    phys_check_params        => hd_check_params
    phys_check_w             => hd_check_w
    phys_get_pthermal        => hd_get_pthermal
    phys_get_v               => hd_get_v
    phys_write_info          => hd_write_info
    phys_handle_small_values => hd_handle_small_values
    phys_angmomfix           => hd_angmomfix

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.

    ! derive units from basic units
    call hd_physical_units()

    if (hd_dust) then
        call dust_init(rho_, mom(:), e_)
    endif

    if (hd_force_diagonal) then
       ! ensure corners are filled, otherwise divide by zero when getting primitives
       !  --> only for debug purposes
       phys_req_diagonal = .true.
    endif

    allocate(tracer(hd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, hd_n_tracer
       tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do

    ! set number of variables which need update ghostcells
    nwgc=nwflux

    ! set the index of the last flux variable for species 1
    stop_indices(1)=nwflux

    if(hd_trac) then
      Tcoff_ = var_set_wextra()
      iw_tcoff=Tcoff_
    else
      Tcoff_ = -1
    end if

    ! initialize thermal conduction module
    if (hd_thermal_conduction) then
      if (.not. hd_energy) call mpistop("thermal conduction needs hd_energy=T")

      phys_req_diagonal = .true.

      call sts_init()
      call tc_init_params(hd_gamma)

      allocate(tc_fl)
      call tc_get_hd_params(tc_fl,tc_params_read_hd)
      tc_fl%get_temperature_from_conserved => hd_get_temperature_from_etot
      call add_sts_method(hd_get_tc_dt_hd,hd_sts_set_source_tc_hd,e_,1,e_,1,&
         .false.)
      call set_conversion_methods_to_head(hd_e_to_ei, hd_ei_to_e)
      call set_error_handling_to_head(hd_tc_handle_small_e)
      tc_fl%get_temperature_from_eint => hd_get_temperature_from_eint
      tc_fl%get_rho => hd_get_rho
      tc_fl%e_ = e_
    end if

    ! Initialize radiative cooling module
    if (hd_radiative_cooling) then
      if (.not. hd_energy) call mpistop("radiative cooling needs hd_energy=T")
      call radiative_cooling_init_params(hd_gamma,He_abundance)
      allocate(rc_fl)
      call radiative_cooling_init(rc_fl,rc_params_read)
      rc_fl%get_rho => hd_get_rho
      rc_fl%get_pthermal => hd_get_pthermal
      rc_fl%e_ = e_
      rc_fl%Tcoff_ = Tcoff_
    end if
    allocate(te_fl_hd)
    te_fl_hd%get_rho=> hd_get_rho
    te_fl_hd%get_pthermal=> hd_get_pthermal
    te_fl_hd%Rfactor = 1d0

    ! Initialize viscosity module
    if (hd_viscosity) call viscosity_init(phys_wider_stencil,&
       phys_req_diagonal)

    ! Initialize gravity module
    if (hd_gravity) call gravity_init()

    ! Initialize rotating_frame module
    if (hd_rotating_frame) call rotating_frame_init()

    ! Initialize CAK radiation force module
    if (hd_cak_force) call cak_init(hd_gamma)

    ! Initialize particles module
    if (hd_particles) then
       call particles_init()
       phys_req_diagonal = .true.
    end if

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1

  end subroutine hd_phys_init


!!start th cond
  ! wrappers for STS functions in thermal_conductivity module
  ! which take as argument the tc_fluid (defined in the physics module)
  subroutine  hd_sts_set_source_tc_hd(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_thermal_conduction, only: sts_set_source_tc_hd
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, igrid, nflux
    double precision, intent(in) ::  x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) ::  wres(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    call sts_set_source_tc_hd(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux,tc_fl)
  end subroutine hd_sts_set_source_tc_hd


  function hd_get_tc_dt_hd(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,dx1,dx2,x) result(dtnew)
    !Check diffusion time limit dt < dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    use mod_thermal_conduction, only: get_tc_dt_hd
 
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: dx1,dx2, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision :: dtnew

    dtnew=get_tc_dt_hd(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,dx1,dx2,x,tc_fl)
  end function hd_get_tc_dt_hd

  
  subroutine hd_tc_handle_small_e(w, x, ixImin1,ixImin2,ixImax1,ixImax2,&
      ixOmin1,ixOmin2,ixOmax1,ixOmax2, step)
    ! move this in a different  routine as in mhd if needed in more places
    use mod_global_parameters
    use mod_small_values

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    integer, intent(in)    :: step

    integer :: idir
    logical :: flag(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    character(len=140) :: error_msg

    flag=.false.
    where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)<small_e) flag(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,e_)=.true.
    if(any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_))) then
      select case (small_values_method)
      case ("replace")
        where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)) w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,e_)=small_e
      case ("average")
        call small_values_average(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2, w, x, flag, e_)
      case default
        ! small values error shows primitive variables
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,e_)*(hd_gamma - 1.0d0)
        do idir = 1, ndir
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               iw_mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               iw_mom(idir))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
        end do
        write(error_msg,"(a,i3)") "Thermal conduction step ", step
        call small_values_error(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2, flag, error_msg)
      end select
    end if
  end subroutine hd_tc_handle_small_e

    ! fill in tc_fluid fields from namelist
    subroutine tc_params_read_hd(fl)
      use mod_global_parameters, only: unitpar,par_files
      use mod_global_parameters, only: unitpar
      type(tc_fluid), intent(inout) :: fl
      integer                      :: n
      logical :: tc_saturate=.false.
      double precision :: tc_k_para=0d0

      namelist /tc_list/ tc_saturate, tc_k_para

      do n = 1, size(par_files)
         open(unitpar, file=trim(par_files(n)), status="old")
         read(unitpar, tc_list, end=111)
111      close(unitpar)
      end do
      fl%tc_saturate = tc_saturate
      fl%tc_k_para = tc_k_para

    end subroutine tc_params_read_hd

  subroutine hd_get_rho(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,rho)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
       x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out) :: rho(ixImin1:ixImax1,ixImin2:ixImax2)

    rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_)

  end subroutine hd_get_rho

!!end th cond
!!rad cool
    subroutine rc_params_read(fl)
      use mod_global_parameters, only: unitpar,par_files
      use mod_constants, only: bigdouble
      use mod_basic_types, only: std_len
      type(rc_fluid), intent(inout) :: fl
      integer                      :: n
      ! list parameters
      integer :: ncool = 4000
      double precision :: cfrac=0.1d0
    
      !> Name of cooling curve
      character(len=std_len)  :: coolcurve='JCcorona'
    
      !> Name of cooling method
      character(len=std_len)  :: coolmethod='exact'
    
      !> Fixed temperature not lower than tlow
      logical    :: Tfix=.false.
    
      !> Lower limit of temperature
      double precision   :: tlow=bigdouble
    
      !> Add cooling source in a split way (.true.) or un-split way (.false.)
      logical    :: rc_split=.false.


      namelist /rc_list/ coolcurve, coolmethod, ncool, cfrac, tlow, Tfix,&
          rc_split
  
      do n = 1, size(par_files)
        open(unitpar, file=trim(par_files(n)), status="old")
        read(unitpar, rc_list, end=111)
111     close(unitpar)
      end do

      fl%ncool=ncool
      fl%coolcurve=coolcurve
      fl%coolmethod=coolmethod
      fl%tlow=tlow
      fl%Tfix=Tfix
      fl%rc_split=rc_split
      fl%cfrac=cfrac
    end subroutine rc_params_read
!! end rad cool

  subroutine hd_check_params
    use mod_global_parameters
    use mod_dust, only: dust_check_params, dust_implicit_update,&
        dust_evaluate_implicit
    use mod_physics, only: phys_implicit_update, phys_evaluate_implicit

    if (.not. hd_energy) then
       if (hd_gamma <= 0.0d0) call mpistop ("Error: hd_gamma <= 0")
       if (hd_adiab < 0.0d0) call mpistop  ("Error: hd_adiab < 0")
       small_pressure= hd_adiab*small_density**hd_gamma
    else
       if (hd_gamma <= 0.0d0 .or. hd_gamma == 1.0d0) call mpistop &
          ("Error: hd_gamma <= 0 or hd_gamma == 1.0")
       small_e = small_pressure/(hd_gamma - 1.0d0)
    end if

    if (hd_dust) call dust_check_params()
    if(use_imex_scheme) then
        ! implicit dust update
        phys_implicit_update => dust_implicit_update
        phys_evaluate_implicit => dust_evaluate_implicit
    endif  

  end subroutine hd_check_params

  subroutine hd_physical_units
    use mod_global_parameters
    double precision :: mp,kB
    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
    else
      mp=mp_cgs
      kB=kB_cgs
    end if
    if(unit_density/=1.d0) then
      unit_numberdensity=unit_density/((1.d0+4.d0*He_abundance)*mp)
    else
      ! unit of numberdensity is independent by default
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
    end if
    if(unit_velocity/=1.d0) then
      unit_pressure=unit_density*unit_velocity**2
      unit_temperature=unit_pressure/((2.d0+&
         3.d0*He_abundance)*unit_numberdensity*kB)
    else if(unit_pressure/=1.d0) then
      unit_temperature=unit_pressure/((2.d0+&
         3.d0*He_abundance)*unit_numberdensity*kB)
      unit_velocity=sqrt(unit_pressure/unit_density)
    else
      ! unit of temperature is independent by default
      unit_pressure=(2.d0+3.d0*He_abundance)&
         *unit_numberdensity*kB*unit_temperature
      unit_velocity=sqrt(unit_pressure/unit_density)
    end if
    if(unit_time/=1.d0) then
      unit_length=unit_time*unit_velocity
    else
      ! unit of length is independent by default
      unit_time=unit_length/unit_velocity
    end if
    unit_mass = unit_density * unit_length**3

  end subroutine hd_physical_units

  !> Returns logical argument flag where values are ok
  subroutine hd_check_w(primitive, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, flag)
    use mod_global_parameters
    use mod_dust, only: dust_check_w

    logical, intent(in)          :: primitive
    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    logical, intent(inout)       :: flag(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision             :: tmp(ixImin1:ixImax1,ixImin2:ixImax2)

    flag=.false.

    if (hd_energy) then
       if (primitive) then
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              e_) < small_pressure) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             e_) = .true.
       else
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (hd_gamma - &
             1.0d0)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) - hd_kin_en(w,&
              ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,&
             ixOmax2))
          where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) < small_pressure) &
             flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = .true.
       endif
    end if

    where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        rho_) < small_density) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_) = .true.

    if(hd_dust) call dust_check_w(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,w,flag)

  end subroutine hd_check_w

  !> Transform primitive variables into conservative ones
  subroutine hd_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, w, x)
    use mod_global_parameters
    use mod_dust, only: dust_to_conserved
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    double precision                :: invgam
    integer                         :: idir, itr

    !!if (fix_small_values) then
    !!  call hd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'hd_to_conserved')
    !!end if

    if (hd_energy) then
       invgam = 1.d0/(hd_gamma - 1.0d0)
       ! Calculate total energy from pressure and kinetic energy
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, e_) * invgam + 0.5d0 * sum(w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, mom(:))**2, dim=ndim+1) * w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, rho_)
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, rho_) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir))
    end do

    if (hd_dust) then
      call dust_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2, w, x)
    end if

  end subroutine hd_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine hd_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, w, x)
    use mod_global_parameters
    use mod_dust, only: dust_to_primitive
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    integer                         :: itr, idir
    double precision                :: inv_rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    if (fix_small_values) then
      call hd_handle_small_values(.false., w, x, ixImin1,ixImin2,ixImax1,&
         ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, 'hd_to_primitive')
    end if

    inv_rho = 1.0d0 / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)

    if (hd_energy) then
       ! Compute pressure
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           e_) = (hd_gamma - 1.0d0) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           e_) - hd_kin_en(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2, inv_rho))
    end if

    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, mom(idir)) * inv_rho
    end do

    ! Convert dust momentum to dust velocity
    if (hd_dust) then
      call dust_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2, w, x)
    end if

  end subroutine hd_to_primitive

  !> Transform internal energy to total energy
  subroutine hd_ei_to_e(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)

    ! Calculate total energy from internal and kinetic energy
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_)+hd_kin_en(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2)

  end subroutine hd_ei_to_e

  !> Transform total energy to internal energy
  subroutine hd_e_to_ei(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)

    ! Calculate ei = e - ek
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_)-hd_kin_en(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2)

  end subroutine hd_e_to_ei

  subroutine e_to_rhos(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, w, x)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision             :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)

    if (hd_energy) then
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           e_) = (hd_gamma - 1.0d0) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           rho_)**(1.0d0 - hd_gamma) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           e_) - hd_kin_en(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2))
    else
       call mpistop("energy from entropy can not be used with -eos = iso !")
    end if
  end subroutine e_to_rhos

  subroutine rhos_to_e(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, w, x)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision             :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)

    if (hd_energy) then
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, rho_)**(hd_gamma - 1.0d0) * w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, e_) / (hd_gamma - 1.0d0) + hd_kin_en(w, ixImin1,&
          ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2)
    else
       call mpistop("entropy from energy can not be used with -eos = iso !")
    end if
  end subroutine rhos_to_e

  !> Calculate v_i = m_i / rho within ixO^L
  subroutine hd_get_v_idim(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, idim, v)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2)

    v(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        mom(idim)) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)
  end subroutine hd_get_v_idim

  !> Calculate velocity vector v_i = m_i / rho within ixO^L
  subroutine hd_get_v(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:2)
    double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)

    integer :: idir

    do idir=1,ndir
      v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, mom(idir)) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          rho_)
    end do

  end subroutine hd_get_v

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine hd_get_cmax(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, idim, cmax)
    use mod_global_parameters
    use mod_dust, only: dust_get_cmax

    integer, intent(in)                       :: ixImin1,ixImin2,ixImax1,&
       ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)              :: w(ixImin1:ixImax1,&
       ixImin2:ixImax2, nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(inout)           :: cmax(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision                          :: csound(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision                          :: v(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    call hd_get_v_idim(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, idim, v)
    call hd_get_csound2(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,csound)
    csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = dsqrt(csound(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))

    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = dabs(v(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))+csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    if (hd_dust) then
      call dust_get_cmax(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2, idim, cmax)
    end if
  end subroutine hd_get_cmax

  subroutine hd_get_a2max(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,a2max)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: a2max(ndim)
    double precision :: a2(ixImin1:ixImax1,ixImin2:ixImax2,ndim,nw)
    integer :: gxOmin1,gxOmin2,gxOmax1,gxOmax2,hxOmin1,hxOmin2,hxOmax1,hxOmax2,&
       jxOmin1,jxOmin2,jxOmax1,jxOmax2,kxOmin1,kxOmin2,kxOmax1,kxOmax2,i,j

    a2=zero
    do i = 1,ndim
      !> 4th order
      hxOmin1=ixOmin1-kr(i,1);hxOmin2=ixOmin2-kr(i,2);hxOmax1=ixOmax1-kr(i,1)
      hxOmax2=ixOmax2-kr(i,2);
      gxOmin1=hxOmin1-kr(i,1);gxOmin2=hxOmin2-kr(i,2);gxOmax1=hxOmax1-kr(i,1)
      gxOmax2=hxOmax2-kr(i,2);
      jxOmin1=ixOmin1+kr(i,1);jxOmin2=ixOmin2+kr(i,2);jxOmax1=ixOmax1+kr(i,1)
      jxOmax2=ixOmax2+kr(i,2);
      kxOmin1=jxOmin1+kr(i,1);kxOmin2=jxOmin2+kr(i,2);kxOmax1=jxOmax1+kr(i,1)
      kxOmax2=jxOmax2+kr(i,2);
      a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i,1:nw)=dabs(-w(kxOmin1:kxOmax1,&
         kxOmin2:kxOmax2,1:nw)+16.d0*w(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
         1:nw)-30.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1:nw)+16.d0*w(hxOmin1:hxOmax1,hxOmin2:hxOmax2,1:nw)-w(gxOmin1:gxOmax1,&
         gxOmin2:gxOmax2,1:nw))
      a2max(i)=maxval(a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i,&
         1:nw))/12.d0/dxlevel(i)**2
    end do
  end subroutine hd_get_a2max

  !> get adaptive cutoff temperature for TRAC (Johnston 2019 ApJL, 873, L22)
  subroutine hd_get_tcutoff(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,tco_local,Tmax_local)
    use mod_global_parameters
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(out) :: tco_local, Tmax_local

    double precision, parameter :: trac_delta=0.25d0
    double precision :: tmp1(ixImin1:ixImax1,ixImin2:ixImax2),&
       Te(ixImin1:ixImax1,ixImin2:ixImax2),lts(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision :: ltr(ixImin1:ixImax1,ixImin2:ixImax2),ltrc,ltrp,&
       tcoff(ixImin1:ixImax1,ixImin2:ixImax2)
    integer :: jxOmin1,jxOmin2,jxOmax1,jxOmax2,hxOmin1,hxOmin2,hxOmax1,hxOmax2
    integer :: jxPmin1,jxPmin2,jxPmax1,jxPmax2,hxPmin1,hxPmin2,hxPmax1,hxPmax2,&
       ixPmin1,ixPmin2,ixPmax1,ixPmax2
    logical :: lrlt(ixImin1:ixImax1,ixImin2:ixImax2)

    
  end subroutine hd_get_tcutoff

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine hd_get_cbounds(wLC, wRC, wLp, wRp, x, ixImin1,ixImin2,ixImax1,&
     ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim,Hspeed,cmax, cmin)
    use mod_global_parameters
    use mod_dust, only: dust_get_cmax
    use mod_variables

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    ! conservative left and right status
    double precision, intent(in)    :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
        nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    ! primitive left and right status
    double precision, intent(in)    :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
        nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:number_species)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:number_species)
    double precision, intent(in)    :: Hspeed(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:number_species)

    double precision :: wmean(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: umean,&
        dmean, csoundL, csoundR, tmp1,tmp2,tmp3
    integer :: ix1,ix2

    select case(boundspeed)
    case (1)
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.

      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsqrt(wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,rho_))
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsqrt(wRp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,rho_))
      tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=1.d0/(dsqrt(wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,rho_))+dsqrt(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)))
      umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idim))*tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idim))*tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

      if(hd_energy) then
        csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=hd_gamma*wLp(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,p_)/wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
        csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=hd_gamma*wRp(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,p_)/wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
      else
        call hd_get_csound2(wLC,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,csoundL)
        call hd_get_csound2(wRC,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,csoundR)
      end if

      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*csoundL(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)+tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*csoundR(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)) * tmp3(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2) + 0.5d0*tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*tmp3(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)**2 * (wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idim))-wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idim)))**2

      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsqrt(dmean(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=umean(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)-dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=umean(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
        if(H_correction) then
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,1)=sign(one,cmin(ix1,ix2,1))*max(abs(cmin(ix1,ix2,1)),&
               Hspeed(ix1,ix2,1))
            cmax(ix1,ix2,1)=sign(one,cmax(ix1,ix2,1))*max(abs(cmax(ix1,ix2,1)),&
               Hspeed(ix1,ix2,1))
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=dabs(umean(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2))+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      end if

      if (hd_dust) then
        wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1:nwflux)=0.5d0*(wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1:nwflux)+wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux))
        call dust_get_cmax(wmean, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2, idim, cmax, cmin)
      end if

    case (2)
      wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1:nwflux)=0.5d0*(wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1:nwflux)+wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux))
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=wmean(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idim))/wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)
      call hd_get_csound2(wmean,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,csoundR)
      csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = dsqrt(csoundR(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))

      if(present(cmin)) then
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=max(tmp1(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2),zero)
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=min(tmp1(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)-csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2),zero)
        if(H_correction) then
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,1)=sign(one,cmin(ix1,ix2,1))*max(abs(cmin(ix1,ix2,1)),&
               Hspeed(ix1,ix2,1))
            cmax(ix1,ix2,1)=sign(one,cmax(ix1,ix2,1))*max(abs(cmax(ix1,ix2,1)),&
               Hspeed(ix1,ix2,1))
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=dabs(tmp1(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2))+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      end if

      if (hd_dust) then
        call dust_get_cmax(wmean, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2, idim, cmax, cmin)
      end if
    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      if(hd_energy) then
        csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=hd_gamma*wLp(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,p_)/wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
        csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=hd_gamma*wRp(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,p_)/wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
      else
        call hd_get_csound2(wLC,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,csoundL)
        call hd_get_csound2(wRC,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,csoundR)
      end if
      csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(dsqrt(csoundL(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2)),dsqrt(csoundR(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=min(wLp(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mom(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idim)))-csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=max(wLp(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mom(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idim)))+csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
        if(H_correction) then
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,1)=sign(one,cmin(ix1,ix2,1))*max(abs(cmin(ix1,ix2,1)),&
               Hspeed(ix1,ix2,1))
            cmax(ix1,ix2,1)=sign(one,cmax(ix1,ix2,1))*max(abs(cmax(ix1,ix2,1)),&
               Hspeed(ix1,ix2,1))
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=max(wLp(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mom(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idim)))+csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      end if
      if (hd_dust) then
        wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1:nwflux)=0.5d0*(wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1:nwflux)+wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux))
        call dust_get_cmax(wmean, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2, idim, cmax, cmin)
      end if
    end select

  end subroutine hd_get_cbounds

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamma*p/rho
  subroutine hd_get_csound2(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(out)   :: csound2(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    call hd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,csound2)
    csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=hd_gamma*csound2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

  end subroutine hd_get_csound2

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho) within ixO^L
  subroutine hd_get_pthermal(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, pth)
    use mod_global_parameters
    use mod_usr_methods, only: usr_set_pthermal
    use mod_small_values, only: trace_small_values

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(out):: pth(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                      :: iw, ix1,ix2

    if (hd_energy) then
       pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (hd_gamma - 1.0d0) * &
          (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) - hd_kin_en(w, ixImin1,&
          ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2))
    else
       if (.not. associated(usr_set_pthermal)) then
          pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = hd_adiab * w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2, rho_)**hd_gamma
       else
          call usr_set_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
             ixOmin2,ixOmax1,ixOmax2,pth)
       end if
    end if

    if (check_small_values) then
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(pth(ix1,ix2)<small_pressure) then
           write(*,*) "Error: small value of gas pressure",pth(ix1,ix2),&
              " encountered when call hd_get_pthermal"
           write(*,*) "Iteration: ", it, " Time: ", global_time
           write(*,*) "Location: ", x(ix1,ix2,:)
           write(*,*) "Cell number: ", ix1,ix2
           do iw=1,nw
             write(*,*) trim(cons_wnames(iw)),": ",w(ix1,ix2,iw)
           end do
           ! use erroneous arithmetic operation to crash the run
           if(trace_small_values) write(*,*) dsqrt(pth(ix1,ix2)-bigdouble)
           write(*,*) "Saving status at the previous time step"
           crash=.true.
         end if
      enddo
      enddo
    end if

    if (fix_small_values) then
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(pth(ix1,ix2)<small_pressure) then
            pth(ix1,ix2)=small_pressure
         endif
      enddo
      enddo
    endif

  end subroutine hd_get_pthermal

  !> Calculate temperature=p/rho when in e_ the  total energy is stored
  subroutine hd_get_temperature_from_etot(w, x, ixImin1,ixImin2,ixImax1,&
     ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2)

    call hd_get_pthermal(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, res)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=res(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
  end subroutine hd_get_temperature_from_etot

  
  !> Calculate temperature=p/rho when in e_ the  internal energy is stored
  subroutine hd_get_temperature_from_eint(w, x, ixImin1,ixImin2,ixImax1,&
     ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (hd_gamma - 1.0d0) * &
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) /w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_)
  end subroutine hd_get_temperature_from_eint

  !these are very similar to the subroutines without 1, used in mod_thermal_conductivity
  !but no check on whether energy variable is present
  subroutine hd_ei_to_e1(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)

    ! Calculate total energy from internal and kinetic energy
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_)+hd_kin_en(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2)

  end subroutine hd_ei_to_e1

  !> Transform total energy to internal energy
  !but no check on whether energy variable is present
  subroutine hd_e_to_ei1(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)

    ! Calculate ei = e - ek
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_)-hd_kin_en(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2)

  end subroutine hd_e_to_ei1

  ! Calculate flux f_idim[iw]
  subroutine hd_get_flux_cons(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, idim, f)
    use mod_global_parameters
    use mod_dust, only: dust_get_flux
    use mod_rotating_frame, only: rotating_frame_velocity

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(out)   :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
        nwflux)
    double precision                :: pth(ixImin1:ixImax1,ixImin2:ixImax2),&
        v(ixImin1:ixImax1,ixImin2:ixImax2),frame_vel(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    integer                         :: idir, itr

    call hd_get_pthermal(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, pth)
    call hd_get_v_idim(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, idim, v)

    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_) = v(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)

    ! Momentum flux is v_i*m_i, +p in direction idim
    do idir = 1, ndir
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = v(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir))
       if (hd_rotating_frame .and. angmomfix .and. idir==phi_) then
          call rotating_frame_velocity(x,ixImin1,ixImin2,ixImax1,ixImax2,&
             ixOmin1,ixOmin2,ixOmax1,ixOmax2,frame_vel)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = f(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2, mom(idir)) + v(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) * frame_vel(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
       end if
    end do

    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idim)) = f(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, mom(idim)) + pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    if(hd_energy) then
      ! Energy flux is v_i*(e + p)
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) = v(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          e_) + pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    end if

    do itr = 1, hd_n_tracer
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, tracer(itr)) = v(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, tracer(itr))
    end do

    ! Dust fluxes
    if (hd_dust) then
      call dust_get_flux(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2, idim, f)
    end if

  end subroutine hd_get_flux_cons

  ! Calculate flux f_idim[iw]
  subroutine hd_get_flux(wC, w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, idim, f)
    use mod_global_parameters
    use mod_dust, only: dust_get_flux_prim
    use mod_viscosity, only: visc_get_flux_prim ! viscInDiv
    use mod_rotating_frame, only: rotating_frame_velocity

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    ! conservative w
    double precision, intent(in)    :: wC(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    ! primitive w
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    double precision, intent(out)   :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
        nwflux)
    double precision                :: pth(ixImin1:ixImax1,ixImin2:ixImax2),&
       frame_vel(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                         :: idir, itr

    if (hd_energy) then
       pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,p_)
    else
       call hd_get_pthermal(wC, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
          ixOmin2,ixOmax1,ixOmax2, pth)
    end if

    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(idim)) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)

    ! Momentum flux is v_i*m_i, +p in direction idim
    do idir = 1, ndir
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(idim)) * wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir))
       if (hd_rotating_frame .and. angmomfix .and. idir==phi_) then
          call mpistop("hd_rotating_frame not implemented yet with angmomfix")
          !One have to compute the frame velocity on cell edge (but we dont know if right of left edge here!!!)
          call rotating_frame_velocity(x,ixImin1,ixImin2,ixImax1,ixImax2,&
             ixOmin1,ixOmin2,ixOmax1,ixOmax2,frame_vel)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = f(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2, mom(idir)) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mom(idim))* wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              rho_) * frame_vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       end if
    end do

    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idim)) = f(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, mom(idim)) + pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    if(hd_energy) then
      ! Energy flux is v_i*(e + p)
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idim)) * (wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          e_) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_))
    end if

    do itr = 1, hd_n_tracer
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, tracer(itr)) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(idim)) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           tracer(itr))
    end do

    ! Dust fluxes
    if (hd_dust) then
      call dust_get_flux_prim(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2, idim, f)
    end if

    ! Viscosity fluxes - viscInDiv
    if (hd_viscosity) then
      call visc_get_flux_prim(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2, idim, f, hd_energy)
    endif

  end subroutine hd_get_flux

  !> Add geometrical source terms to w
  !>
  !> Notice that the expressions of the geometrical terms depend only on ndir,
  !> not ndim. Eg, they are the same in 2.5D and in 3D, for any geometry.
  !>
  !> Ileyk : to do :
  !>     - address the source term for the dust in case (coordinate == spherical)
  subroutine hd_add_source_geom(qdt, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, wCT, w, x)
    use mod_global_parameters
    use mod_usr_methods, only: usr_set_surface
    use mod_viscosity, only: visc_add_source_geom ! viscInDiv
    use mod_rotating_frame, only: rotating_frame_add_source
    use mod_dust, only: dust_n_species, dust_mom, dust_rho
    use mod_geometry
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2, 1:nw)
    ! to change and to set as a parameter in the parfile once the possibility to
    ! solve the equations in an angular momentum conserving form has been
    ! implemented (change tvdlf.t eg)
    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2),&
        source(ixImin1:ixImax1,ixImin2:ixImax2), minrho
    integer                         :: iw,idir, h1xmin1,h1xmin2,h1xmax1,&
       h1xmax2, h2xmin1,h2xmin2,h2xmax1,h2xmax2
    integer :: mr_,mphi_ ! Polar var. names
    integer :: irho, ifluid, n_fluids
    double precision :: exp_factor(ixImin1:ixImax1,ixImin2:ixImax2),&
        del_exp_factor(ixImin1:ixImax1,ixImin2:ixImax2),&
        exp_factor_primitive(ixImin1:ixImax1,ixImin2:ixImax2)

    if (hd_dust) then
       n_fluids = 1 + dust_n_species
    else
       n_fluids = 1
    end if

    select case (coordinate)

    case(Cartesian_expansion)
      !the user provides the functions of exp_factor and del_exp_factor
      if(associated(usr_set_surface)) call usr_set_surface(ixImin1,ixImin2,&
         ixImax1,ixImax2,x,block%dx,exp_factor,del_exp_factor,&
         exp_factor_primitive)
      call hd_get_pthermal(wCT, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2, source)
      source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = source(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*del_exp_factor(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)/exp_factor(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1)) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(1)) + qdt*source(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    case (cylindrical)
       do ifluid = 0, n_fluids-1
          ! s[mr]=(pthermal+mphi**2/rho)/radius
          if (ifluid == 0) then
             ! gas
             irho  = rho_
             mr_   = mom(r_)
             if(phi_>0) mphi_ = mom(phi_)
             call hd_get_pthermal(wCT, x, ixImin1,ixImin2,ixImax1,ixImax2,&
                 ixOmin1,ixOmin2,ixOmax1,ixOmax2, source)
             minrho = 0.0d0
          else
             ! dust : no pressure
             irho  = dust_rho(ifluid)
             mr_   = dust_mom(r_, ifluid)
             if(phi_>0) mphi_ = dust_mom(phi_, ifluid)
             source(ixImin1:ixImax1,ixImin2:ixImax2) = zero
             minrho = 0.0d0
          end if
          if (phi_ > 0) then
             where (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, irho) > minrho)
                source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
                   source(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2) + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                    mphi_)**2 / wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, irho)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mr_) = w(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2, mr_) + qdt * source(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2, r_)
             end where
             ! s[mphi]=(-mphi*mr/rho)/radius
             if(.not. angmomfix) then
                where (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, irho) > minrho)
                   source(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2) = -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       mphi_) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       mr_) / wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, irho)
                   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       mphi_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       mphi_) + qdt * source(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       r_)
                end where
             end if
          else
             ! s[mr]=2pthermal/radius
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mr_) = w(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2, mr_) + qdt * source(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2, r_)
          end if
       end do
    case (spherical)
       if (hd_dust) then
          call mpistop&
("Dust geom source terms not implemented yet with spherical geometries")
       end if
       mr_   = mom(r_)
       if(phi_>0) mphi_ = mom(phi_)
       h1xmin1=ixOmin1-kr(1,1);h1xmin2=ixOmin2-kr(1,2)
       h1xmax1=ixOmax1-kr(1,1);h1xmax2=ixOmax2-kr(1,2)
       h2xmin1=ixOmin1-kr(2,1);h2xmin2=ixOmin2-kr(2,2)
       h2xmax1=ixOmax1-kr(2,1);h2xmax2=ixOmax2-kr(2,2);
       ! s[mr]=((mtheta**2+mphi**2)/rho+2*p)/r
       call hd_get_pthermal(wCT, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
          ixOmin2,ixOmax1,ixOmax2, pth)
       source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = pth(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1) - block%surfaceC(h1xmin1:h1xmax1,h1xmin2:h1xmax2,&
           1)) /block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       if (ndir > 1) then
         do idir = 2, ndir
           source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = source(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2) + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               mom(idir))**2 / wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)
         end do
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mr_) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, mr_) + qdt * source(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2, 1)

       
       ! s[mtheta]=-(mr*mtheta/rho)/r+cot(theta)*(mphi**2/rho+p)/r
       source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = pth(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1) * (block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           2) - block%surfaceC(h2xmin1:h2xmax1,h2xmin2:h2xmax2,&
           2)) / block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       if (ndir == 3) then
          source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = source(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) + (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(3))**2 / wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              rho_)) / tan(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2, 2))
       end if
       if (.not. angmomfix) then
          source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = source(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) - (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(2)) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mr_)) / wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(2)) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, mom(2)) + qdt * source(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2, 1)

       if (ndir == 3) then
         ! s[mphi]=-(mphi*mr/rho)/r-cot(theta)*(mtheta*mphi/rho)/r
         if (.not. angmomfix) then
           source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = -(wCT(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2, mom(3)) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               mr_)) / wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               rho_)- (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               mom(2)) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               mom(3))) / wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               rho_) / tan(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2, 2))
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(3)) = w(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2, mom(3)) + qdt * source(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2, 1)
         end if
       end if
      
    end select

    if (hd_viscosity) call visc_add_source_geom(qdt,ixImin1,ixImin2,ixImax1,&
       ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x)

    if (hd_rotating_frame) then
       if (hd_dust) then
          call mpistop("Rotating frame not implemented yet with dust")
       else
          call rotating_frame_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,&
             ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
       end if
    end if

  end subroutine hd_add_source_geom

  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine hd_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,wCT,w,x,qsourcesplit,active,wCTprim)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_dust, only: dust_add_source, dust_mom, dust_rho, dust_n_species
    use mod_viscosity, only: viscosity_add_source
    use mod_usr_methods, only: usr_gravity
    use mod_gravity, only: gravity_add_source, grav_split
    use mod_cak_force, only: cak_add_source

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active
    double precision, intent(in),optional :: wCTprim(ixImin1:ixImax1,&
       ixImin2:ixImax2, 1:nw)

    double precision :: gravity_field(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    integer :: idust, idim

    if(hd_dust .and. .not. use_imex_scheme) then
      call dust_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,wCT,w,x,qsourcesplit,active)
    end if

    if(hd_radiative_cooling) then
      call radiative_cooling_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,&
         ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x,qsourcesplit,active, rc_fl)
    end if

    if(hd_viscosity) then
      call viscosity_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,wCT,w,x,hd_energy,qsourcesplit,active)
    end if

    if (hd_gravity) then
      call gravity_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,wCT,w,x,hd_energy,qsourcesplit,active)

      if (hd_dust .and. qsourcesplit .eqv. grav_split) then
         active = .true.

         call usr_gravity(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
            ixOmax1,ixOmax2, wCT, x, gravity_field)
         do idust = 1, dust_n_species
            do idim = 1, ndim
               w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idim,&
                   idust)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idim,&
                   idust)) + qdt * gravity_field(ixOmin1:ixOmax1,&
                  ixOmin2:ixOmax2, idim) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   dust_rho(idust))
            end do
         end do
      end if
    end if

    if (hd_cak_force) then
      call cak_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,wCT,w,x,hd_energy,qsourcesplit,active)
    end if

  end subroutine hd_add_source

  subroutine hd_get_dt(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, dtnew, dx1,dx2, x)
    use mod_global_parameters
    use mod_dust, only: dust_get_dt
    use mod_radiative_cooling, only: cooling_get_dt
    use mod_viscosity, only: viscosity_get_dt
    use mod_gravity, only: gravity_get_dt
    use mod_cak_force, only: cak_get_dt

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: dx1,dx2, x(ixImin1:ixImax1,&
       ixImin2:ixImax2, 1:2)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble

    if(hd_dust) then
      call dust_get_dt(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2, dtnew, dx1,dx2, x)
    end if

    if(hd_radiative_cooling) then
      call cooling_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,dtnew,dx1,dx2,x,rc_fl)
    end if

    if(hd_viscosity) then
      call viscosity_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
    end if

    if(hd_gravity) then
      call gravity_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
   end if

   if (hd_cak_force) then
     call cak_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
        ixOmax2,dtnew,dx1,dx2,x)
   end if

  end subroutine hd_get_dt

  function hd_kin_en(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, inv_rho) result(ke)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)                    :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)           :: w(ixImin1:ixImax1,&
       ixImin2:ixImax2, nw)
    double precision                       :: ke(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    double precision, intent(in), optional :: inv_rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    if (present(inv_rho)) then
       ke = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(:))**2,&
           dim=ndim+1) * inv_rho
    else
       ke = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(:))**2,&
           dim=ndim+1) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)
    end if
  end function hd_kin_en

  function hd_inv_rho(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2) result(inv_rho)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision              :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    ! Can make this more robust
    inv_rho = 1.0d0 / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)
  end function hd_inv_rho

  subroutine hd_handle_small_values(primitive, w, x, ixImin1,ixImin2,ixImax1,&
     ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, subname)
    ! handles hydro (density,pressure,velocity) bootstrapping
    ! any negative dust density is flagged as well (and throws an error)
    ! small_values_method=replace also for dust
    use mod_global_parameters
    use mod_small_values
    use mod_dust, only: dust_n_species, dust_mom, dust_rho
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: n,idir
    logical :: flag(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    if (small_values_method == "ignore") return

    call hd_check_w(primitive, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, w, flag)

    if (any(flag)) then
      select case (small_values_method)
      case ("replace")
        where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)) w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,rho_) = small_density
        do idir = 1, ndir
          if(small_values_fix_iw(mom(idir))) then
            where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               rho_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = 0.0d0
          end if
        end do
        if(hd_energy)then
          if(small_values_fix_iw(e_)) then
            if(primitive) then
              where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 rho_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                  p_) = small_pressure
            else
              where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 rho_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                  e_) = small_e + hd_kin_en(w,ixImin1,ixImin2,ixImax1,ixImax2,&
                 ixOmin1,ixOmin2,ixOmax1,ixOmax2)
            endif
          end if
        endif

        if(hd_energy) then
          if(primitive) then
            where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)) w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,p_) = small_pressure
          else
            where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_))
              ! Add kinetic energy
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = small_e + hd_kin_en(w,&
                 ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
                 ixOmax2)
            end where
          end if
        end if

        if(hd_dust)then
           do n=1,dust_n_species
              where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 dust_rho(n))) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 dust_rho(n)) = 0.0d0
              do idir = 1, ndir
                  where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     dust_rho(n))) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     dust_mom(idir,n)) = 0.0d0
              enddo
           enddo
        endif
      case ("average")
        if(primitive)then
           ! averaging for all primitive fields, including dust
           call small_values_average(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
              ixOmin2,ixOmax1,ixOmax2, w, x, flag)
        else
           ! do averaging of density
           call small_values_average(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
              ixOmin2,ixOmax1,ixOmax2, w, x, flag, rho_)
           if(hd_energy) then
             ! do averaging of pressure
             w(ixImin1:ixImax1,ixImin2:ixImax2,&
                p_)=(hd_gamma-1.d0)*(w(ixImin1:ixImax1,ixImin2:ixImax2,&
                e_) -0.5d0*sum(w(ixImin1:ixImax1,ixImin2:ixImax2, mom(:))**2,&
                 dim=ndim+1)/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_))
             call small_values_average(ixImin1,ixImin2,ixImax1,ixImax2,&
                 ixOmin1,ixOmin2,ixOmax1,ixOmax2, w, x, flag, p_)
             w(ixImin1:ixImax1,ixImin2:ixImax2,e_)=w(ixImin1:ixImax1,&
                ixImin2:ixImax2,p_)/(hd_gamma-1.d0) &
                +0.5d0*sum(w(ixImin1:ixImax1,ixImin2:ixImax2, mom(:))**2,&
                 dim=ndim+1)/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
           end if
           if(hd_dust)then
              do n=1,dust_n_species
                 where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                    dust_rho(n))) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                    dust_rho(n)) = 0.0d0
                 do idir = 1, ndir
                    where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       dust_rho(n))) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       dust_mom(idir,n)) = 0.0d0
                 enddo
              enddo
          endif
        endif
      case default
        if(.not.primitive) then
          !convert w to primitive
          ! Calculate pressure = (gamma-1) * (e-ek)
          if(hd_energy) then
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               p_)=(hd_gamma-1.d0)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               e_)-hd_kin_en(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
               ixOmax1,ixOmax2))
          end if
          ! Convert gas momentum to velocity
          do idir = 1, ndir
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2, mom(idir))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               rho_)
          end do
        end if
        ! NOTE: dust entries may still have conserved values here
        call small_values_error(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2, flag, subname)
      end select
    end if
  end subroutine hd_handle_small_values

end module mod_hd_phys
