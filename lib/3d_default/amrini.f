!> Generate and initialize all grids at the coarsest level (level one)
subroutine initlevelone
  use mod_global_parameters
  use mod_ghostcells_update

  integer :: iigrid, igrid

  levmin=1
  levmax=1

  call init_forest_root

  call getigrids
  call build_connectivity

  ! fill solution space of all root grids
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     call alloc_node(igrid)
     ! in case gradient routine used in initial condition, ensure geometry known
     call initial_condition(igrid)
  end do
  

  ! update ghost cells
  call getbc(global_time,0.d0,ps,iwstart,nwgc)

end subroutine initlevelone

!> fill in initial condition
subroutine initial_condition(igrid)
  ! Need only to set the mesh values (can leave ghost cells untouched)
  use mod_usr_methods, only: usr_init_one_grid
  use mod_global_parameters

  integer, intent(in) :: igrid

  ! in case gradient routine used in initial condition, ensure geometry known
  block=>ps(igrid)
  dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
  dxlevel(3)=rnode(rpdx3_,igrid);

  if (.not. associated(usr_init_one_grid)) then
     call mpistop("usr_init_one_grid not defined")
  else
     call usr_init_one_grid(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,&
        ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,ps(igrid)%w,ps(igrid)%x)
  end if

end subroutine initial_condition

!> modify initial condition
subroutine modify_IC
  use mod_usr_methods, only: usr_init_one_grid
  use mod_global_parameters

  integer :: iigrid, igrid

  do iigrid=1,igridstail; igrid=igrids(iigrid);
     block=>ps(igrid)
     dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
     dxlevel(3)=rnode(rpdx3_,igrid);

     if (.not. associated(usr_init_one_grid)) then
        call mpistop("usr_init_one_grid not defined")
     else
        call usr_init_one_grid(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
           ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,ps(igrid)%w,ps(igrid)%x)
     end if
  end do

end subroutine modify_IC


!> improve initial condition after initialization
subroutine improve_initial_condition()
  use mod_global_parameters
  use mod_usr_methods
  use mod_constrained_transport
  use mod_multigrid_coupling
  use mod_physics
  use mod_ghostcells_update

  logical :: active

  if(associated(usr_improve_initial_condition)) then
    call usr_improve_initial_condition
  else if(stagger_grid) then
    if(associated(usr_init_vector_potential)) then
      ! re-calculate magnetic field from the vector potential in a 
      ! completely divergency free way for AMR mesh in 3D
      if(levmax>levmin.and.ndim==3) call recalculateB
    end if
    if(slab_uniform.and.associated(phys_clean_divb)) then
      ! Project out the divB using multigrid poisson solver 
      ! if not initialised from vector potential
      if(.not.use_multigrid) call mg_setup_multigrid()
      call phys_clean_divb(global_time,0.d0,active)
      call getbc(global_time,0.d0,ps,iwstart,nwgc)
    end if
  end if

end subroutine improve_initial_condition

