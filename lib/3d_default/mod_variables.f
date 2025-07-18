module mod_variables
  use mod_basic_types

  implicit none
  public

  !> Number of flux variables
  integer           :: nwflux = 0

  !> Number of flux variables which need user to specify boundary type
  integer           :: nwfluxbc = 0

  !> Number of auxiliary variables in w
  integer           :: nwaux = 0

  !> Number of extra variables in w
  integer           :: nwextra = 0

  !> Number of extra variables in wextra seperated from w
  integer           :: nw_extra = 0

  !> Total number of variables
  integer           :: nw = 0

  !> Total number of stagger variables
  integer           :: nws = 0

  !> Number of variables which need to be updated in ghost cells
  integer           :: nwgc = 0

  !> Number of vector variables (used for writing output)
  integer           :: nvector = 0

  !> Indices of vector variables
  integer, dimension(:), allocatable :: iw_vector

  ! the number of the first w variable to exchange ghost cells
  integer            :: iwstart=1

  !> Maximum number of variables
  integer, parameter :: max_nw = 50

  !> Primitive variable names
  character(len=name_len) :: prim_wnames(max_nw)

  !> Conservative variable names
  character(len=name_len) :: cons_wnames(max_nw)

  ! Global indices of variables that are often used

  !> Index of the (gas) density
  integer :: iw_rho = -1

  !> Indices of the momentum density
  integer, allocatable :: iw_mom(:)

  !> Index of the energy density
  integer :: iw_e = -1

  !> Index of the radiation energy density
  integer :: iw_r_e = -1

  !> Index of the internal energy density
  integer  :: iw_eaux = -1

  !> Indices of the magnetic field components
  integer, allocatable, protected :: iw_mag(:)

  !> Index of wplus. This is used when no AWs in the system and we revert back to an older version of mod_uawsom_phys.t
  integer :: iw_wplus = -1

  !> Index of wminus. This is used when no AWs in the system and we revert back to an older version of mod_uawsom_phys.t
  integer :: iw_wminus = -1

  !> Index of wkplus
  integer :: iw_wkplus = -1
  
  !> Index of wkminus
  integer :: iw_wkminus = -1

  !> Index of wAplus
  integer :: iw_wAplus = -1
  
  !> Index of wAminus
  integer :: iw_wAminus = -1

  !> Index of the cutoff temperature for the TRAC method
  integer :: iw_tcoff = -1

  !> number of species: each species has different characterictic speeds and should
  !> be used accordingly in mod_finite_volume and mod_finite_difference
  integer :: number_species = 1

  !> index of the var
  !> whose velocity appears in the induction eq.
  integer :: index_v_mag = 1

  !> the indices in 1:nwflux array are assumed consecutive for each species
  !> this array should be of size number_species and contain the first index in the array of
  !> the number_species
  integer, allocatable :: start_indices(:)
  !> the indices in 1:nwflux array are assumed consecutive for each species
  !> this array should be of size number_species and contain the last index in the array of
  !> the first number_species, the last index for the last one is nwflux
  integer, allocatable :: stop_indices(:)

  ! indices of equi for the species index_v_mag
  ! these are needed for hlld solver, TODO: consider moving in a separate file
  integer :: iw_equi_rho = -1
  integer :: iw_equi_p = -1

contains

  !> Set generic flux variable
  function var_set_fluxvar(name_cons, name_prim, ix, need_bc) result(iw)
    character(len=*), intent(in)  :: name_cons !< Conservative name
    character(len=*), intent(in)  :: name_prim !< Primitive name
    integer, intent(in), optional :: ix !< Optional index (to make var1, var2, ...)
    logical, intent(in), optional :: need_bc !< Require boundary condition (default: true)
    integer                       :: iw
    logical                       :: add_bc

    nwflux = nwflux + 1
    nw     = nw + 1
    iw     = nwflux

    add_bc = .true.
    if (present(need_bc)) add_bc = need_bc
    if (add_bc) nwfluxbc = nwfluxbc + 1

    if (.not. present(ix)) then
      prim_wnames(nwflux) = name_cons
      cons_wnames(nwflux) = name_prim
    else
      write(cons_wnames(nwflux),"(A,I0)") name_cons, ix
      write(prim_wnames(nwflux),"(A,I0)") name_prim, ix
    end if
  end function var_set_fluxvar

  !> Set extra variable in w, which is not advected and has no boundary conditions.
  !> This has to be done after defining flux variables and auxiliary variables.
  function var_set_extravar(name_cons, name_prim, ix) result(iw)
    character(len=*), intent(in)  :: name_cons, name_prim
    integer, intent(in), optional :: ix
    integer                       :: iw

    nwextra = nwextra + 1
    nw      = nw + 1
    iw      = nw

    if (.not. present(ix)) then
      prim_wnames(iw) = name_cons
      cons_wnames(iw) = name_prim
    else
      write(cons_wnames(iw),"(A,I0)") name_cons, ix
      write(prim_wnames(iw),"(A,I0)") name_prim, ix
    end if
  end function var_set_extravar

  !> Set extra variable in wextra, which is not advected and has no boundary conditions and not output in dat.
  !> This has to be done after defining flux variables and auxiliary variables.
  function var_set_wextra() result(iw)
    integer :: iw

    nw_extra = nw_extra + 1
    iw      = nw_extra

  end function var_set_wextra

  !> Set auxiliary variable, which is not advected but has boundary conditions.
  !> This has to be done after defining flux variables.
  function var_set_auxvar(name_cons, name_prim, ix) result(iw)
    character(len=*), intent(in)  :: name_cons, name_prim
    integer, intent(in), optional :: ix
    integer                       :: iw

    nwaux   = nwaux + 1
    nw      = nw + 1
    iw      = nw

    if (.not. present(ix)) then
      prim_wnames(iw) = name_cons
      cons_wnames(iw) = name_prim
    else
      write(cons_wnames(iw),"(A,I0)") name_cons, ix
      write(prim_wnames(iw),"(A,I0)") name_prim, ix
    end if
  end function var_set_auxvar

  !> Set density variable
  function var_set_rho() result(iw)
    integer :: iw

    nwflux              = nwflux + 1
    nwfluxbc            = nwfluxbc + 1
    nw                  = nw + 1
    iw_rho              = nwflux
    iw                  = nwflux
    prim_wnames(nwflux) = 'rho'
    cons_wnames(nwflux) = 'rho'
  end function var_set_rho

  !> Set momentum variables
  function var_set_momentum(ndir) result(iw)
    integer, intent(in) :: ndir
    integer             :: iw(ndir), idir

    if (allocated(iw_mom)) call mpistop("Error: set_mom was already called")
    allocate(iw_mom(ndir))

    do idir = 1, ndir
      nwflux       = nwflux + 1
      nwfluxbc     = nwfluxbc + 1
      nw           = nw + 1
      iw_mom(idir) = nwflux
      iw(idir)     = nwflux
      write(cons_wnames(nwflux),"(A1,I1)") "m", idir
      write(prim_wnames(nwflux),"(A1,I1)") "v", idir
    end do
  end function var_set_momentum

  !> Set energy variable
  function var_set_energy() result(iw)
    integer :: iw

    nwflux              = nwflux + 1
    nwfluxbc            = nwfluxbc + 1
    nw                  = nw + 1
    iw_e                = nwflux
    iw                  = nwflux
    cons_wnames(nwflux) = 'e'
    prim_wnames(nwflux) = 'p'
  end function var_set_energy

  function var_set_radiation_energy() result(iw)
    integer :: iw

    nwflux              = nwflux + 1
    nwfluxbc            = nwfluxbc + 1
    nw                  = nw + 1
    iw_r_e              = nwflux
    iw                  = nwflux
    cons_wnames(nwflux) = 'r_e'
    prim_wnames(nwflux) = 'r_e'
  end function var_set_radiation_energy

  !> Set magnetic field variables
  function var_set_bfield(ndir) result(iw)
    integer, intent(in) :: ndir
    integer             :: iw(ndir), idir

    if (allocated(iw_mag)) call mpistop("Error: set_mag was already called")
    allocate(iw_mag(ndir))

    do idir = 1, ndir
      nwflux       = nwflux + 1
      nwfluxbc     = nwfluxbc + 1
      nw           = nw + 1
      iw_mag(idir) = nwflux
      iw(idir)     = nwflux
      write(cons_wnames(nwflux),"(A1,I1)") "b", idir
      write(prim_wnames(nwflux),"(A1,I1)") "b", idir
    end do
  end function var_set_bfield

  !> Set internal energy variable
  function var_set_internal_energy() result(iw)
    integer :: iw

    nwflux              = nwflux + 1
    nw                  = nw + 1
    iw_eaux             = nwflux
    iw                  = nwflux
    cons_wnames(nwflux) = 'eaux'
    prim_wnames(nwflux) = 'paux'
  end function var_set_internal_energy

  !> For the UAWSOM physics module, set W+ and W-

  function var_set_wplus() result(iw)
    integer :: iw

    nwflux              = nwflux + 1
    nwfluxbc            = nwfluxbc + 1
    nw                  = nw + 1
    iw_wplus            = nwflux
    iw                  = nwflux
    prim_wnames(nwflux) = 'wplus'
    cons_wnames(nwflux) = 'wplus'
  end function var_set_wplus
  
  function var_set_wminus() result(iw)
    integer :: iw

    nwflux              = nwflux + 1
    nwfluxbc            = nwfluxbc + 1
    nw                  = nw + 1
    iw_wminus           = nwflux
    iw                  = nwflux
    prim_wnames(nwflux) = 'wminus'
    cons_wnames(nwflux) = 'wminus'
  end function var_set_wminus

  !> For the UAWSOM physics module, set WA+ and WA- for AWs and Wk+ and Wk- for the kink waves. With the AWs included need to differentiate between AWs and Kink waves
  
  function var_set_wkplus() result(iw)
    integer :: iw

    nwflux              = nwflux + 1
    nwfluxbc            = nwfluxbc + 1
    nw                  = nw + 1
    iw_wkplus            = nwflux
    iw                  = nwflux
    prim_wnames(nwflux) = 'wkplus'
    cons_wnames(nwflux) = 'wkplus'
  end function var_set_wkplus
  
  function var_set_wkminus() result(iw)
    integer :: iw

    nwflux              = nwflux + 1
    nwfluxbc            = nwfluxbc + 1
    nw                  = nw + 1
    iw_wkminus           = nwflux
    iw                  = nwflux
    prim_wnames(nwflux) = 'wkminus'
    cons_wnames(nwflux) = 'wkminus'
  end function var_set_wkminus

  function var_set_wAplus() result(iw)
    integer :: iw

    nwflux              = nwflux + 1
    nwfluxbc            = nwfluxbc + 1
    nw                  = nw + 1
    iw_wAplus            = nwflux
    iw                  = nwflux
    prim_wnames(nwflux) = 'wAplus'
    cons_wnames(nwflux) = 'wAplus'
  end function var_set_wAplus
  
  function var_set_wAminus() result(iw)
    integer :: iw

    nwflux              = nwflux + 1
    nwfluxbc            = nwfluxbc + 1
    nw                  = nw + 1
    iw_wAminus           = nwflux
    iw                  = nwflux
    prim_wnames(nwflux) = 'wAminus'
    cons_wnames(nwflux) = 'wAminus'
  end function var_set_wAminus

end module mod_variables
