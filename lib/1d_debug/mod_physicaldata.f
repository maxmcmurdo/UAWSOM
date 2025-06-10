module mod_physicaldata
  implicit none
  save

  type state
     !> ID of a grid block
     integer :: igrid=-1
     !> index range of block array in cell centers
     integer :: ixGmin1,ixGmax1
     !> index range of block array in cell faces
     integer :: ixGsmin1,ixGsmax1
     !> level of AMR
     integer :: level
     !> If it face a physical boundary
     logical, dimension(:), pointer :: is_physical_boundary(:) =>Null()
     !> Variables, normally cell center conservative values
     double precision, dimension(:,:), allocatable :: w
     !> Variables, cell face values
     double precision, dimension(:,:), allocatable :: ws
     !> Variables, cell edge values
     double precision, dimension(:,:), allocatable :: we
     !> Variables, cell corner values
     double precision, dimension(:,:), allocatable :: wc
     !> extra variables do not need ghost cell and equation flux
     double precision, dimension(:,:), pointer :: wextra=>Null()
     !> Time-independent magnetic field at cell center and cell interface
     double precision, dimension(:,:,:), pointer :: B0=>Null()
     !> Time-independent electric current density at cell center
     double precision, dimension(:,:), pointer :: J0=>Null()
     !> Time-independent equi vars (B0 is not into this array)
     double precision, dimension(:,:,:), pointer :: equi_vars=>Null()
     !> Cell-center positions
     double precision, dimension(:,:), pointer :: x=>Null()
     !> Cell sizes in coordinate units
     double precision, dimension(:,:), pointer :: dx=>Null()
     !> Cell sizes at cell center in length unit
     double precision, dimension(:,:), pointer :: ds=>Null()
     !> Cell sizes at cell face in length unit
     double precision, dimension(:,:), pointer :: dsC=>Null()
     !> Volumes of a cell
     double precision, dimension(:), pointer :: dvolume=>Null()
     !> Areas of cell-center surfaces
     double precision, dimension(:,:), pointer :: surface=>Null()
     !> Areas of cell-face surfaces
     double precision, dimension(:,:), pointer :: surfaceC=>Null()
     !> special values for a block
     double precision, dimension(:), pointer :: special_values=>Null()
  end type state



  type state_sub
     !> ID of a grid block
     integer :: igrid=-1
     !> Variables, normally center
     double precision, dimension(:), allocatable :: w
     !> Variables for the cornerpositions on the slice 
     double precision, dimension(:), allocatable :: wC
     !> Variables, normally center, one level coarser representative
     double precision, dimension(:), allocatable :: wcoarse
     !> Cell-center positions
     double precision, dimension(:), allocatable :: x
     !> Corner positions on the slice
     double precision, dimension(:), allocatable :: xC
     !> Cell-center positions, one level coarser representative
     double precision, dimension(:), allocatable :: xcoarse
  end type state_sub

  type grid_field
     !> Variables new state
     double precision, dimension(:,:), allocatable :: w
     !> Variables old state
     double precision, dimension(:,:), allocatable :: wold
  end type grid_field
  !> buffer for pole boundary
  type(state) :: pole_buf

  !> array of physical states for all blocks on my processor
  type(state), dimension(:), allocatable, target :: ps
  !> array of physical states, temp 1 for multi-step time integrator
  type(state), dimension(:), allocatable, target :: ps1
  !> array of physical states, temp 2 for multi-step time integrator
  type(state), dimension(:), allocatable, target :: ps2
  !> array of physical states, temp 3 for multi-step time integrator
  type(state), dimension(:), allocatable, target :: ps3
  !> array of physical states, temp 4 for multi-step time integrator
  type(state), dimension(:), allocatable, target :: ps4
  !> array of physical states, at the beginning of each iteration
  type(state), dimension(:), allocatable, target :: pso
  !> array of physical blocks, one level coarser representative
  type(state), dimension(:), allocatable, target :: psc

  !> array of physical blocks in reduced dimension
  type(state_sub), dimension(:), allocatable, target :: ps_sub


  double precision, dimension(:), allocatable :: collapsedData


  !> array of physical blocks of meshed fields for particles
  type(grid_field), dimension(:), allocatable, target :: gridvars

  !> velocities store for constrained transport
  type ct_velocity
    double precision, dimension(:,:), allocatable :: vnorm,cbarmin,cbarmax
    double precision, dimension(:,:,:), allocatable :: vbarC,vbarLC,vbarRC
  end type ct_velocity

end module mod_physicaldata
