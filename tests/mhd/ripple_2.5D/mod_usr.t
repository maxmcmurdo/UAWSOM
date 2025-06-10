module mod_usr
  use mod_mhd
  implicit none

contains

  subroutine usr_init()
    usr_init_one_grid => initonegrid_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output 
    usr_set_B0          => specialset_B0

    call set_coordinate_system("Cartesian_2.5D")

    call mhd_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
  ! initialize one grid
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: ixC^L

    w(ixO^S,mom(1))=0.d0
    w(ixO^S,mom(2))=0.d0
    w(ixO^S,mom(3))=0.01d0*dexp(-(^D&x(ixO^S,^D)**2+)/1d-4)
    w(ixO^S,p_)    =1.d0+0.1*dexp(-(^D&x(ixO^S,^D)**2+)/1d-4)
    w(ixO^S,rho_)  =1.d0+0.1*dexp(-(^D&x(ixO^S,^D)**2+)/1d-4)
    if(stagger_grid) then
      ixCmax^D=ixOmax^D;
      ixCmin^D=ixOmin^D-kr(1,^D);
      block%ws(ixC^S,1)=1.d0
      call mhd_face_to_center(ixO^L,block)
    else if(B0field)then
      w(ixO^S,mag(:))=0.d0
    else
      w(ixO^S,mag(1))=1.d0
      w(ixO^S,mag(2))=0.d0
      w(ixO^S,mag(3))=0.d0
    end if
    
    call mhd_to_conserved(ixI^L,ixO^L,w,x)
  
  end subroutine initonegrid_usr

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: tmp(ixI^S)

    call mhd_get_pthermal(w,x,ixI^L,ixO^L,tmp)
    ! output the temperature
    w(ixO^S,nw+1)=tmp(ixO^S)/w(ixO^S,rho_)
    ! output the plasma beta p*2/B**2
    if(B0field)then
      w(ixO^S,nw+2)=tmp(ixO^S)*two/sum((w(ixO^S,mag(:))+&
                    block%B0(ixO^S,:,0))**2,dim=ndim+1)
    else
      w(ixO^S,nw+2)=tmp(ixO^S)*two/sum(w(ixO^S,mag(:))**2,dim=ndim+1)
    endif
    ! output divB1
    call get_divb(w,ixI^L,ixO^L,tmp)
    w(ixO^S,nw+3)=tmp(ixO^S)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames
    varnames='Te beta divb'

  end subroutine specialvarnames_output

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a steady (time-independent) potential or 
  ! linear force-free background field
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)

    wB0(ixO^S,1)=1.d0
    wB0(ixO^S,2)=0.d0
    wB0(ixO^S,3)=0.d0

  end subroutine specialset_B0

end module mod_usr
