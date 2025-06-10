!> This module handles the initialization of various components of amrvac
module mod_initialize

  implicit none
  private

  logical :: initialized_already = .false.

  ! Public methods
  public :: initialize_amrvac

contains

  !> Initialize amrvac: read par files and initialize variables
  subroutine initialize_amrvac()
    use mod_global_parameters
    use mod_input_output
    use mod_physics, only: phys_check, phys_check_params
    use mod_usr_methods, only: usr_set_parameters
    use mod_bc_data, only: bc_data_init
    use mod_trac, only:init_trac_line, init_trac_block
    use mod_init_datafromfile, only: read_data_init

    if (initialized_already) return

    ! add auxiliary variable(s) to update boundary ghost cells
    nwgc=nwgc+nwaux

    ! Check whether the user has loaded a physics module
    call phys_check()

    ! Read input files
    call read_par_files()
    call initialize_vars()
    call init_comm_types()

    ! Possibly load boundary condition data or initial data
    call bc_data_init()
    call read_data_init()

    if(associated(usr_set_parameters)) call usr_set_parameters()

    call phys_check_params()

    initialized_already = .true.
  end subroutine initialize_amrvac

  !> Initialize (and allocate) simulation and grid variables
  subroutine initialize_vars
    use mod_forest
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve, only: pflux
    use mod_amr_fct, only: pface, fine_neighbors, old_neighbor
    use mod_geometry

    integer :: igrid, level, ipe, ig1
    logical :: ok

    allocate(ps(max_blocks))
    allocate(ps1(max_blocks))
    allocate(ps2(max_blocks))
    allocate(ps3(max_blocks))
    allocate(ps4(max_blocks))
    allocate(pso(max_blocks))
    allocate(psc(max_blocks))
    allocate(ps_sub(max_blocks))
    allocate(neighbor(2,-1:1,max_blocks),neighbor_child(2,0:3,max_blocks))
    allocate(neighbor_type(-1:1,max_blocks),neighbor_active(-1:1,max_blocks))
    allocate(neighbor_pole(-1:1,max_blocks))
    allocate(igrids(max_blocks),igrids_active(max_blocks),&
       igrids_passive(max_blocks))
    allocate(rnode(rnodehi,max_blocks),rnode_sub(rnodehi,max_blocks))
    allocate(node(nodehi,max_blocks),node_sub(nodehi,max_blocks),&
       phyboundblock(max_blocks))
    allocate(pflux(2,1,max_blocks))
    ! allocate mesh for particles
    if(use_particles) allocate(gridvars(max_blocks))
    if(stagger_grid) then
      allocate(pface(2,1,max_blocks),fine_neighbors(2,1,max_blocks))
      allocate(old_neighbor(2,-1:11,max_blocks))
    end if

    it=it_init
    global_time=time_init

    dt=zero

    ! no poles initially
    neighbor_pole=0

    ! check resolution
    if (mod(ixGhi1,2)/=0) then
       call mpistop("mesh widths must give even number grid points")
    end if
    ixMlo1=ixGlo1+nghostcells;ixMhi1=ixGhi1-nghostcells;

    if (nbufferx1>(ixMhi1-ixMlo1+1)) then
       write(unitterm,*) "nbufferx^D bigger than mesh size makes no sense."
       write(unitterm,*) "Decrease nbufferx or increase mesh size"
       call mpistop("")
    end if

    ! initialize dx arrays on finer (>1) levels
    do level=2,refine_max_level
       dx(1,level) = dx(1,level-1) * half  ! refine ratio 2
    end do

    ! domain decomposition
    ! physical extent of a grid block at level 1, per dimension
    dg1(1)=dx(1,1)*dble(block_nx1)
    ! number of grid blocks at level 1 in simulation domain, per dimension
    ng1(1)=nint((xprobmax1-xprobmin1)/dg1(1))
    ! total number of grid blocks at level 1
    nglev1=ng1(1)

    do level=2,refine_max_level
       dg1(level)=half*dg1(level-1);
       ng1(level)=ng1(level-1)*2;
    end do

    ! check that specified stepsize correctly divides domain
    ok=((abs(dble(ng1(1))*dg1(1)-(xprobmax1-xprobmin1))<=smalldouble))
    if (.not.ok) then
       write(unitterm,*)"domain cannot be divided by meshes of given gridsize"
       call mpistop("domain cannot be divided by meshes of given gridsize")
    end if

    poleB=.false.
    if (.not.slab) call set_pole

    ! number of grid blocks at level 1 along a dimension, which does not have a pole or periodic boundary, 
    ! must be larger than 1 for a rectangular AMR mesh
    if((ng1(1)/=1).and.refine_max_level>1) then
      
      if(ng1(1)==1.and..not.poleB(1,1).and..not.poleB(2,&
         1).and..not.periodB(1).and..not.aperiodB(1)) then
        write(unitterm,"(a,i2,a)")&
            "number of grid blocks at level 1 in dimension",1,&
           " be larger than 1 for a rectangular AMR mesh!"
        write(unitterm,"(a,i1)") "increase domain_nx",1
        call mpistop("")
      end if
      
    end if

    ! initialize connectivity data
    igridstail=0

    ! allocate memory for forest data structures
    allocate(level_head(refine_max_level),level_tail(refine_max_level))
    do level=1,refine_max_level
       nullify(level_head(level)%node,level_tail(level)%node)
    end do

    allocate(igrid_to_node(max_blocks,0:npe-1))
    do ipe=0,npe-1
       do igrid=1,max_blocks
          nullify(igrid_to_node(igrid,ipe)%node)
       end do
    end do

    allocate(sfc(1:3,max_blocks*npe))

    allocate(igrid_to_sfc(max_blocks))

    sfc=0
    allocate(Morton_start(0:npe-1),Morton_stop(0:npe-1))
    allocate(Morton_sub_start(0:npe-1),Morton_sub_stop(0:npe-1))

    allocate(nleafs_level(1:nlevelshi))

    allocate(coarsen(max_blocks,0:npe-1),refine(max_blocks,0:npe-1))
    coarsen=.false.
    refine=.false.
    if (nbufferx1/=0) then
       allocate(buffer(max_blocks,0:npe-1))
       buffer=.false.
    end if
    allocate(igrid_inuse(max_blocks,0:npe-1))
    igrid_inuse=.false.

    allocate(tree_root(1:ng1(1)))
    do ig1=1,ng1(1)
    nullify(tree_root(ig1)%node)
    end do

    ! define index ranges and MPI send/receive derived datatype for ghost-cell swap
    call init_bc()
    type_send_srl=>type_send_srl_f
    type_recv_srl=>type_recv_srl_f
    type_send_r=>type_send_r_f
    type_recv_r=>type_recv_r_f
    type_send_p=>type_send_p_f
    type_recv_p=>type_recv_p_f
    call create_bc_mpi_datatype(iwstart,nwgc)

  end subroutine initialize_vars


end module mod_initialize
