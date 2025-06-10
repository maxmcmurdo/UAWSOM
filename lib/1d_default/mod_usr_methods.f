!> Module with all the methods that users can customize in AMRVAC
!>
!> Each procedure pointer can be initialized in a user's mod_usr.t
module mod_usr_methods

  implicit none
  public

  !> Initialize the user's settings (after initializing amrvac)
  procedure(p_no_args), pointer       :: usr_set_parameters   => null()
  !> Initialize earch grid block data
  procedure(init_one_grid), pointer   :: usr_init_one_grid    => null()

  ! Boundary condition related
  procedure(special_bc), pointer      :: usr_special_bc       => null()
  procedure(special_mg_bc), pointer   :: usr_special_mg_bc    => null()

  procedure(internal_bc), pointer     :: usr_internal_bc      => null()

  ! Output related
  procedure(p_no_args), pointer       :: usr_print_log        => null()
  procedure(p_no_args), pointer       :: usr_write_analysis   => null()
  procedure(transform_w), pointer     :: usr_transform_w      => null()
  procedure(aux_output), pointer      :: usr_aux_output       => null()
  procedure(add_aux_names), pointer   :: usr_add_aux_names    => null()
  procedure(sub_modify_io), pointer   :: usr_modify_output    => null()
  procedure(special_convert), pointer :: usr_special_convert  => null()

  ! Called at the beginning of every time step (after determining dt)
  procedure(process_grid), pointer    :: usr_process_grid     => null()
  procedure(process_global), pointer  :: usr_process_global   => null()

  ! Called every time step just after advance (with w^(n+1), it^n, global_time^n)
  procedure(process_adv_grid), pointer   :: usr_process_adv_grid   => null()
  procedure(process_adv_global), pointer :: usr_process_adv_global => null()

  ! Called after initial condition before the start of the simulation
  procedure(p_no_args), pointer       :: usr_improve_initial_condition => &
     null()

  ! Called before the start of the simulation
  procedure(p_no_args), pointer       :: usr_before_main_loop => null()

  ! Source terms
  procedure(source), pointer          :: usr_source           => null()
  procedure(get_dt), pointer          :: usr_get_dt           => null()
  procedure(phys_gravity), pointer    :: usr_gravity          => null()

  ! Usr defined dust drag force
  procedure(phys_dust_get_dt), pointer :: usr_dust_get_dt     => null()
  procedure(phys_dust_get_3d_dragforce), pointer :: usr_get_3d_dragforce => &
     null()

  ! Usr defined space varying viscosity
  procedure(phys_visco), pointer      :: usr_setvisco         => null()

  ! Usr defined thermal pressure for hydro & energy=.False.
  procedure(hd_pthermal), pointer     :: usr_set_pthermal     => null()

  ! Refinement related procedures
  procedure(refine_grid), pointer     :: usr_refine_grid      => null()
  procedure(var_for_errest), pointer  :: usr_var_for_errest   => null()
  procedure(a_refine_threshold), pointer :: usr_refine_threshold => null()
  procedure(flag_grid), pointer       :: usr_flag_grid        => null()

  ! Set time-independent magnetic field for B0 splitting
  procedure(set_B0), pointer          :: usr_set_B0           => null()
  ! Set time-independent variables for equilibrium splitting, except for B0
  procedure(set_equi_vars), pointer   :: usr_set_equi_vars           => null()
  ! Set time-independent current density for B0 splitting
  procedure(set_J0), pointer          :: usr_set_J0           => null()
  procedure(special_resistivity), pointer :: usr_special_resistivity => null()

  ! Particle module related
  procedure(update_payload), pointer    :: usr_update_payload    => null()
  procedure(create_particles), pointer  :: usr_create_particles  => null()
  procedure(check_particle), pointer    :: usr_check_particle    => null()
  procedure(particle_fields), pointer   :: usr_particle_fields   => null()
  procedure(particle_analytic), pointer :: usr_particle_analytic => null()
  procedure(particle_position), pointer :: usr_particle_position => null()

  ! Radiation quantity related
  procedure(special_opacity), pointer   :: usr_special_opacity => null()
  procedure(special_aniso_opacity), pointer   :: usr_special_aniso_opacity => &
     null()
  procedure(special_opacity_qdot), pointer   :: usr_special_opacity_qdot => &
     null()
  procedure(special_fluxlimiter), pointer   :: usr_special_fluxlimiter => &
     null()
  procedure(special_diffcoef), pointer   :: usr_special_diffcoef => null()

  ! Called after the mesh has been adjuste
  procedure(after_refine), pointer      :: usr_after_refine => null()

  ! initialize vector potential on cell edges for magnetic field
  procedure(init_vector_potential), pointer :: usr_init_vector_potential => &
     null()

  ! allow user to change inductive electric field, especially for boundary driven applications
  procedure(set_electric_field), pointer :: usr_set_electric_field => null()

  ! allow user to specify variables at physical boundaries
  procedure(set_wLR), pointer :: usr_set_wLR => null()

  ! allow user to specify the expansion function for the surface of a cross sectional
  ! area of a 1D prominence, along with the analytical derivative of that function and its
  ! primitive shape evaluated in the boundaries \int_(x_i-dx_i/2)^(x_i+dx_i/2) A(s) ds
  procedure(set_surface), pointer  :: usr_set_surface          => null()

  ! for tracing field. allow user to specify variables and field
  procedure(set_field_w), pointer :: usr_set_field_w => null()
  procedure(set_field), pointer :: usr_set_field => null()

  abstract interface

    subroutine p_no_args()
    end subroutine p_no_args

    !> Initialize one grid
    subroutine init_one_grid(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
      double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    end subroutine init_one_grid

     !> special boundary types, users must assign conservative
     !> variables in boundaries
     subroutine special_bc(qt,ixImin1,ixImax1,ixOmin1,ixOmax1,iB,w,x)
       use mod_global_parameters
       !> Shape of input arrays
       integer, intent(in)             :: ixImin1,ixImax1
       !> Region where boundary values have to be set
       integer, intent(in)             :: ixOmin1,ixOmax1
       !> Integer indicating direction of boundary
       integer, intent(in)             :: iB
       double precision, intent(in)    :: qt, x(ixImin1:ixImax1,1:ndim)
       double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
     end subroutine special_bc

     !> Special boundary type for radiation hydrodynamics module, only used to
     !> set the boundary conditions for the radiation energy.
     subroutine special_mg_bc(iB)
       use mod_global_parameters
       integer, intent(in)             :: iB
     end subroutine special_mg_bc

    !> internal boundary, user defined
    !> This subroutine can be used to artificially overwrite ALL conservative
    !> variables in a user-selected region of the mesh, and thereby act as
    !> an internal boundary region. It is called just before external (ghost cell)
    !> boundary regions will be set by the BC selection. Here, you could e.g.
    !> want to introduce an extra variable (nwextra, to be distinguished from nwaux)
    !> which can be used to identify the internal boundary region location.
    !> Its effect should always be local as it acts on the mesh.
    subroutine internal_bc(level,qt,ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1,level
      double precision, intent(in)    :: qt
      double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
      double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    end subroutine internal_bc

    !> this subroutine is ONLY to be used for computing auxiliary variables
    !> which happen to be non-local (like div v), and are in no way used for
    !> flux computations. As auxiliaries, they are also not advanced
    subroutine process_grid(igrid,level,ixImin1,ixImax1,ixOmin1,ixOmax1,qt,w,&
       x)
      use mod_global_parameters
      integer, intent(in)             :: igrid,level,ixImin1,ixImax1,ixOmin1,&
         ixOmax1
      double precision, intent(in)    :: qt,x(ixImin1:ixImax1,1:ndim)
      double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    end subroutine process_grid

    !> If defined, this routine is called before writing output, and it can
    !> set/modify the variables in the w array.
    subroutine sub_modify_io(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
      double precision, intent(in)    :: qt,x(ixImin1:ixImax1,1:ndim)
      double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    end subroutine sub_modify_io

    !> This subroutine is called at the beginning of each time step
    !> by each processor. No communication is specified, so the user
    !> has to implement MPI routines if information has to be shared
    subroutine process_global(iit,qt)
      use mod_global_parameters
      integer, intent(in)          :: iit
      double precision, intent(in) :: qt
    end subroutine process_global

    !> for processing after the advance (PIC-MHD, e.g.)
    subroutine process_adv_grid(igrid,level,ixImin1,ixImax1,ixOmin1,ixOmax1,qt,&
       w,x)
      use mod_global_parameters
      integer, intent(in)             :: igrid,level,ixImin1,ixImax1,ixOmin1,&
         ixOmax1
      double precision, intent(in)    :: qt,x(ixImin1:ixImax1,1:ndim)
      double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    end subroutine process_adv_grid

    !> for processing after the advance (PIC-MHD, e.g.)
    subroutine process_adv_global(iit,qt)
      use mod_global_parameters
      integer, intent(in)          :: iit
      double precision, intent(in) :: qt
    end subroutine process_adv_global

    !> this subroutine can be used in convert, to add auxiliary variables to the
    !> converted output file, for further analysis using tecplot, paraview, ....
    !> these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    !
    !> the array normconv can be filled in the (nw+1:nw+nwauxio) range with
    !> corresponding normalization values (default value 1)
    subroutine aux_output(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,normconv)
      use mod_global_parameters
      integer, intent(in)          :: ixImin1,ixImax1,ixOmin1,ixOmax1
      double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
      double precision             :: w(ixImin1:ixImax1,nw+nwauxio)
      double precision             :: normconv(0:nw+nwauxio)
    end subroutine aux_output

    !> Add names for the auxiliary variables
    subroutine add_aux_names(varnames)
      use mod_global_parameters
      character(len=*) :: varnames
    end subroutine add_aux_names

    !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
    !> iw=iwmin...iwmax.  wCT is at time qCT
    subroutine source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,wCT,&
       qt,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1,&
          iwmin,iwmax
      double precision, intent(in)    :: qdt, qtC, qt
      double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw),&
          x(ixImin1:ixImax1,1:ndim)
      double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    end subroutine source

    !> Limit "dt" further if necessary, e.g. due to the special source terms.
    !> The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
    !> module have already been called.
    subroutine get_dt(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dtnew,dx1,x)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in)    :: dx1, x(ixImin1:ixImax1,1:ndim)
      double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)
      double precision, intent(inout) :: dtnew
    end subroutine get_dt

    !> Calculate gravitational acceleration in each dimension
    subroutine phys_gravity(ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,x,&
       gravity_field)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
      double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw)
      double precision, intent(out)   :: gravity_field(ixImin1:ixImax1,ndim)
    end subroutine phys_gravity

    !> Calculate the 3d drag force of gas onto dust
    subroutine phys_dust_get_3d_dragforce(ixImin1,ixImax1, ixOmin1,ixOmax1, w,&
        x, fdrag, ptherm, vgas,dust_n_species)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1,&
          dust_n_species
      double precision, intent(in)    :: x(ixImin1:ixImax1, 1:ndim)
      double precision, intent(in)    :: w(ixImin1:ixImax1, 1:nw)
      double precision, intent(out)   :: fdrag(ixImin1:ixImax1, 1:ndir,&
          1:dust_n_species)
      double precision, intent(in)    :: ptherm(ixImin1:ixImax1),&
          vgas(ixImin1:ixImax1, ndir)
    end subroutine phys_dust_get_3d_dragforce

    !> Calculate the time step associated with the usr drag force
    subroutine phys_dust_get_dt(w, ixImin1,ixImax1, ixOmin1,ixOmax1, dtdust,&
        dx1, x, dust_n_species)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1,&
          dust_n_species
      double precision, intent(in)    :: dx1, x(ixImin1:ixImax1,1:ndim)
      double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)
      double precision, intent(inout) :: dtdust(1:dust_n_species)
    end subroutine phys_dust_get_dt

    !>Calculation anormal viscosity depending on space
    subroutine phys_visco(ixImin1,ixImax1,ixOmin1,ixOmax1,x,w,mu)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
      double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)
      double precision, intent(out)   :: mu(ixImin1:ixImax1)
    end subroutine phys_visco

    !>Calculation anormal pressure for hd & energy=.False.
    subroutine hd_pthermal(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,pth)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
      double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)
      double precision, intent(out)   :: pth(ixImin1:ixImax1)
    end subroutine hd_pthermal

    !> Set the "eta" array for resistive MHD based on w or the
    !> "current" variable which has components between idirmin and 3.
    subroutine special_resistivity(w,ixImin1,ixImax1,ixOmin1,ixOmax1,idirmin,x,&
       current,eta)
      use mod_global_parameters
      integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1,&
          idirmin
      double precision, intent(in) :: w(ixImin1:ixImax1,nw), x(ixImin1:ixImax1,&
         1:ndim)
      double precision             :: current(ixImin1:ixImax1,7-2*ndir:3),&
          eta(ixImin1:ixImax1)
    end subroutine special_resistivity


    !> Set user defined opacity for use in diffusion coeff, heating and cooling, and radiation force
    subroutine special_opacity(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,kappa)
      use mod_global_parameters
      integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in) :: w(ixImin1:ixImax1,1:nw),&
          x(ixImin1:ixImax1,1:ndim)
      double precision, intent(out):: kappa(ixOmin1:ixOmax1)
    end subroutine special_opacity

    !> Set user defined, anisotropic opacity for use in diffusion coeff, heating and cooling, and radiation force
    subroutine special_aniso_opacity(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,kappa,&
       idir)
      use mod_global_parameters
      integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1, idir
      double precision, intent(in) :: w(ixImin1:ixImax1,1:nw),&
          x(ixImin1:ixImax1,1:ndim)
      double precision, intent(out):: kappa(ixOmin1:ixOmax1)
    end subroutine special_aniso_opacity

    !> Set user defined opacity for use in diffusion coeff, heating and cooling, and radiation force. Overwrites special_opacity
    subroutine special_opacity_qdot(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,kappa)
      use mod_global_parameters
      integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in) :: w(ixImin1:ixImax1,1:nw),&
          x(ixImin1:ixImax1,1:ndim)
      double precision, intent(out):: kappa(ixOmin1:ixOmax1)
    end subroutine special_opacity_qdot

    !> Set user defined FLD flux limiter, lambda
    subroutine special_fluxlimiter(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,&
       fld_lambda,fld_R)
      use mod_global_parameters
      integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in) :: w(ixImin1:ixImax1,1:nw),&
          x(ixImin1:ixImax1,1:ndim)
      double precision, intent(out):: fld_lambda(ixImin1:ixImax1),&
         fld_R(ixImin1:ixImax1)
    end subroutine special_fluxlimiter

    !> Set user defined FLD diffusion coefficient
    subroutine special_diffcoef(w, wCT, x, ixImin1,ixImax1, ixOmin1,ixOmax1)
      use mod_global_parameters
      integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(inout) :: w(ixImin1:ixImax1, 1:nw)
      double precision, intent(in) :: wCT(ixImin1:ixImax1, 1:nw)
      double precision, intent(in) :: x(ixImin1:ixImax1, 1:ndim)
    end subroutine special_diffcoef

    !> Enforce additional refinement or coarsening
    !> One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    !> you must set consistent values for integers refine/coarsen:
    !> refine = -1 enforce to not refine
    !> refine =  0 doesn't enforce anything
    !> refine =  1 enforce refinement
    !> coarsen = -1 enforce to not coarsen
    !> coarsen =  0 doesn't enforce anything
    !> coarsen =  1 enforce coarsen
    !> e.g. refine for negative first coordinate x < 0 as
    !> if (any(x(ix^S,1) < zero)) refine=1
    subroutine refine_grid(igrid,level,ixImin1,ixImax1,ixOmin1,ixOmax1,qt,w,x,&
       refine,coarsen)
      use mod_global_parameters
      integer, intent(in)          :: igrid, level, ixImin1,ixImax1, ixOmin1,&
         ixOmax1
      double precision, intent(in) :: qt, w(ixImin1:ixImax1,1:nw),&
          x(ixImin1:ixImax1,1:ndim)
      integer, intent(inout)       :: refine, coarsen
    end subroutine refine_grid

    !> this is the place to compute a local auxiliary variable to be used
    !> as refinement criterion for the Lohner error estimator only
    !>  -->it is then requiring and iflag>nw
    !> note that ixO=ixI=ixG, hence the term local (gradients need special attention!)
    subroutine var_for_errest(ixImin1,ixImax1,ixOmin1,ixOmax1,iflag,w,x,var)
      use mod_global_parameters
      integer, intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1,iflag
      double precision, intent(in)  :: w(ixImin1:ixImax1,1:nw),&
          x(ixImin1:ixImax1,1:ndim)
      double precision, intent(out) :: var(ixImin1:ixImax1)
    end subroutine var_for_errest

    !> Here one can add a steady (time-independent) potential background field
    subroutine set_B0(ixImin1,ixImax1,ixOmin1,ixOmax1,x,wB0)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
      double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
      double precision, intent(inout) :: wB0(ixImin1:ixImax1,1:ndir)
    end subroutine set_B0

    !> Here one can add a time-independent background current density
    subroutine set_J0(ixImin1,ixImax1,ixOmin1,ixOmax1,x,wJ0)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
      double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
      double precision, intent(inout) :: wJ0(ixImin1:ixImax1,7-2*ndir:ndir)
    end subroutine set_J0

    !> Here one can add a steady (time-independent) equi vars
    subroutine set_equi_vars(ixImin1,ixImax1,ixOmin1,ixOmax1,x,w0)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
      double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
      double precision, intent(inout) :: w0(ixImin1:ixImax1,&
         1:number_equi_vars)
    end subroutine set_equi_vars

    !> adjust w when restart from dat file with different w variables
    subroutine transform_w(ixImin1,ixImax1,ixOmin1,ixOmax1,nw_in,w_in,x,w_out)
      use mod_global_parameters
      integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1, nw_in
      double precision, intent(in)  :: w_in(ixImin1:ixImax1,1:nw_in)
      double precision, intent(in)  :: x(ixImin1:ixImax1, 1:ndim)
      double precision, intent(out) :: w_out(ixImin1:ixImax1,1:nw)
    end subroutine transform_w

    !> use different threshold in special regions for AMR to
    !> reduce/increase resolution there where nothing/something interesting happens.
    subroutine a_refine_threshold(wlocal,xlocal,threshold,qt,level)
      use mod_global_parameters
      double precision, intent(in)    :: wlocal(1:nw),xlocal(1:ndim),qt
      double precision, intent(inout) :: threshold
      integer, intent(in) :: level
    end subroutine a_refine_threshold

    !> Allow user to use their own data-postprocessing procedures
    subroutine special_convert(qunitconvert)
      use mod_global_parameters
      integer, intent(in) :: qunitconvert
      character(len=20)   :: userconvert_type
    end subroutine special_convert

    !> flag=-1 : Treat all cells active, omit deactivation (onentry, default)
    !> flag=0  : Treat as normal domain
    !> flag=1  : Treat as passive, but reduce by safety belt
    !> flag=2  : Always treat as passive
    subroutine flag_grid(qt,ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,flag)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
      integer, intent(inout)          :: flag
      double precision, intent(in)    :: qt
      double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
      double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    end subroutine flag_grid

    !> Update payload of particles
    subroutine update_payload(igrid,w,wold,xgrid,x,u,q,m,mypayload,mynpayload,&
       particle_time)
      use mod_global_parameters
      integer, intent(in)           :: igrid,mynpayload
      double precision, intent(in)  :: w(ixGlo1:ixGhi1,1:nw),&
         wold(ixGlo1:ixGhi1,1:nw)
      double precision, intent(in)  :: xgrid(ixGlo1:ixGhi1,1:ndim),x(1:ndir),&
         u(1:ndir),q,m,particle_time
      double precision, intent(out) :: mypayload(mynpayload)
    end subroutine update_payload

    !> Create particles
    subroutine create_particles(n_particles, x, v, q, m, follow)
      integer, intent(in)           :: n_particles
      double precision, intent(out) :: x(3, n_particles)
      double precision, intent(out) :: v(3, n_particles)
      double precision, intent(out) :: q(n_particles)
      double precision, intent(out) :: m(n_particles)
      logical, intent(out)          :: follow(n_particles)
    end subroutine create_particles

    !> Check arbitrary particle conditions or modifications
    subroutine check_particle(igrid,x,v,q,m,follow,check)
      use mod_global_parameters
      integer, intent(in)           :: igrid
      double precision, intent(inout) :: x(1:ndir)
      double precision, intent(inout) :: v(1:ndir),q,m
      logical, intent(inout) :: follow
      logical, intent(out)   :: check
    end subroutine check_particle

    !> Associate fields to particle
    subroutine particle_fields(w, x, E, B)
      use mod_global_parameters
      double precision, intent(in)  :: w(ixGlo1:ixGhi1,1:nw)
      double precision, intent(in)  :: x(ixGlo1:ixGhi1,1:ndim)
      double precision, intent(out) :: E(ixGlo1:ixGhi1, ndir)
      double precision, intent(out) :: B(ixGlo1:ixGhi1, ndir)
    end subroutine particle_fields

    subroutine particle_analytic(ix, x, tloc, vec)
      use mod_global_parameters
      integer, intent(in)           :: ix(ndir) !< Indices in gridvars
      double precision, intent(in)  :: x(ndir)
      double precision, intent(in)  :: tloc
      double precision, intent(out) :: vec(ndir)
    end subroutine particle_analytic

    !> User-defined particle movement
    subroutine particle_position(x, n, tloc, tlocnew)
      use mod_global_parameters
      integer, intent(in)             :: n
      double precision, intent(inout) :: x(3)
      double precision, intent(in)    :: tloc, tlocnew
    end subroutine particle_position

    subroutine after_refine(n_coarsen, n_refine)
      integer, intent(in) :: n_coarsen
      integer, intent(in) :: n_refine
    end subroutine after_refine

    !> initialize vector potential on cell edges for magnetic field
    subroutine init_vector_potential(ixImin1,ixImax1, ixCmin1,ixCmax1, xC, A,&
        idir)
      use mod_global_parameters

      integer, intent(in)                :: ixImin1,ixImax1, ixCmin1,ixCmax1,&
          idir
      double precision, intent(in)       :: xC(ixImin1:ixImax1,1:ndim)
      double precision, intent(out)      :: A(ixImin1:ixImax1)

    end subroutine init_vector_potential

    ! allow user to change inductive electric field, especially for boundary driven applications
    subroutine set_electric_field(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,qdt,fE,s)
      use mod_global_parameters
      integer, intent(in)                :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in)       :: qt, qdt
      type(state)                        :: s
      double precision, intent(inout)    :: fE(ixImin1:ixImax1,7-2*ndim:3)

      !integer :: ixC^L,ixA^L
      ! For example, to set inductive electric field at bottom boundary in a 3D box for induction equation
      ! v and b are  from observational data for data-driven application

      !associate(w=>s%w,ws=>s%ws)

      !if(s%is_physical_boundary(5)) then
      !  ixCmin^D=ixOmin^D-1;
      !  ixCmax^D=ixOmax^D;
      !  ixAmin^D=ixCmin^D;
      !  ixAmax^D=ixCmax^D+1;
      !  fE(nghostcells^%3ixA^S,1)=-ws(nghostcells^%3ixA^S,3)*w(nghostcells^%3ixA^S,mom(2))
      !  fE(nghostcells^%3ixA^S,2)= ws(nghostcells^%3ixA^S,3)*w(nghostcells^%3ixA^S,mom(1))
      !  ixAmin^D=ixCmin^D+kr(2,^D);
      !  ixAmax^D=ixCmax^D+kr(2,^D);
      !  fE(nghostcells^%3ixC^S,1)=0.5d0*(fE(nghostcells^%3ixC^S,1)+fE(nghostcells^%3ixA^S,1))*&
      !                            qdt*s%dsC(nghostcells^%3ixC^S,1)
      !  ixAmin^D=ixCmin^D+kr(1,^D);
      !  ixAmax^D=ixCmax^D+kr(1,^D);
      !  fE(nghostcells^%3ixC^S,2)=0.5d0*(fE(nghostcells^%3ixC^S,2)+fE(nghostcells^%3ixA^S,2))*&
      !                            qdt*s%dsC(nghostcells^%3ixC^S,2)
      !end if

      !end associate

    end subroutine set_electric_field

    !> allow user to specify 'variables' left and right state at physical boundaries to control flux through the boundary surface
    subroutine set_wLR(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,wLC,wRC,wLp,wRp,s,&
       idir)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1,&
          idir
      double precision, intent(in)    :: qt
      double precision, intent(inout) :: wLC(ixImin1:ixImax1,1:nw),&
          wRC(ixImin1:ixImax1,1:nw)
      double precision, intent(inout) :: wLp(ixImin1:ixImax1,1:nw),&
          wRp(ixImin1:ixImax1,1:nw)
      type(state)                     :: s

      !if(s%is_physical_boundary(3).and.idir==2) then
      !  wLp(ixOmin2^%2ixO^S,mom(1))=1.d0
      !  wRp(ixOmin2^%2ixO^S,mom(1))=wRp(ixOmin2^%2ixO^S,mom(1))
      !  wLC(ixOmin2^%2ixO^S,mom(1))=wLp(ixOmin2^%2ixO^S,mom(1))*wLp(ixOmin2^%2ixO^S,rho_)
      !  wRC(ixOmin2^%2ixO^S,mom(1))=wRp(ixOmin2^%2ixO^S,mom(1))*wRp(ixOmin2^%2ixO^S,rho_)
      !end if
    end subroutine set_wLR

    subroutine set_surface(ixImin1,ixImax1,x,delx,exp_factor,del_exp_factor,&
       exp_factor_primitive)
      use mod_global_parameters
      integer, intent(in)              :: ixImin1,ixImax1
      double precision, intent(in)     :: delx(ixImin1:ixImax1,1:ndim),&
          x(ixImin1:ixImax1,1:ndim)
      double precision, intent(out)    :: exp_factor(ixImin1:ixImax1),&
          del_exp_factor(ixImin1:ixImax1)
      double precision, intent(out)    :: &
         exp_factor_primitive(ixImin1:ixImax1)

    end subroutine set_surface

    subroutine set_field_w(igrid,ip,xf,wP,wL,numP,nwP,nwL,dL,forward,ftype,&
       tcondi)
      use mod_global_parameters
      !use mod_point_searching

      integer, intent(in)                 :: igrid,ip,numP,nwP,nwL
      double precision, intent(in)        :: xf(numP,ndim)
      double precision, intent(inout)     :: wP(numP,nwP),wL(1+nwL)
      double precision, intent(in)        :: dL
      logical, intent(in)                 :: forward
      character(len=std_len), intent(in)  :: ftype,tcondi

      !double precision :: xpp(1:ndim),wpp(1:nw)

      !! nwP=2,nwL=0. get rho/T at line
      !if (tcondi=='user') then
      !  xpp(1:ndim)=xf(ip,1:ndim)
      !  call get_point_w_ingrid(igrid,xpp,wpp,'primitive')
      !  wP(ip,1)=wpp(rho_)
      !  wP(ip,2)=wpp(p_)/wpp(rho_)
      !endif

    end subroutine set_field_w

    subroutine set_field(xfn,igrid,field,ftype)
      use mod_global_parameters

      integer,intent(in)                  :: igrid
      double precision, intent(in)        :: xfn(ndim)
      double precision, intent(inout)     :: field(ndim)
      character(len=std_len), intent(in)  :: ftype

      !if (ftype='xdir') then
      !  field(:)=zero
      !  field(1)=1.d0
      !endif

    end subroutine set_field

  end interface

end module mod_usr_methods
