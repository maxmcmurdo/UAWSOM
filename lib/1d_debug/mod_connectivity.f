!> This module contains variables that describe the connectivity of the mesh and
!> also data structures for connectivity-related communication.
module mod_connectivity
   implicit none
   save

   integer, parameter :: neighbor_boundary = 1
   integer, parameter :: neighbor_coarse = 2
   integer, parameter :: neighbor_sibling = 3
   integer, parameter :: neighbor_fine = 4

   integer, dimension(:,:,:), allocatable :: neighbor
   integer, dimension(:,:,:), allocatable :: neighbor_child
   integer, dimension(:,:), allocatable :: neighbor_type
   logical, dimension(:,:), allocatable :: neighbor_active
   integer, dimension(:,:), allocatable :: neighbor_pole

   ! grid number array per processor
   integer, dimension(:), allocatable :: igrids
   integer, dimension(:), allocatable :: igrids_active
   integer, dimension(:), allocatable :: igrids_passive
   ! number of grids on current processor
   integer :: igridstail
   integer :: igridstail_active
   integer :: igridstail_passive

   integer, dimension(1) :: nrecv_fc, nsend_fc
   ! cc for corner coarse
   integer, dimension(1) :: nrecv_cc, nsend_cc


end module mod_connectivity
