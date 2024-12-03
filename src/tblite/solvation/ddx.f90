! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/solvation/ddx.f90
!> Provides a polarizable continuum model

!> Implicit solvation model based on a polarizable dielectric continuum
module tblite_solvation_ddx
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   ! use tblite_blas, only : dot, gemv
   ! use tblite_container_cache, only : container_cache
   use tblite_mesh_lebedev, only : grid_size, get_angular_grid, list_bisection
   ! use tblite_scf_info, only : scf_info, atom_resolved
   ! use tblite_scf_potential, only : potential_type
   ! use tblite_wavefunction_type, only : wavefunction_type
   use tblite_solvation_data, only : get_vdw_rad_cosmo
   use tblite_solvation_type, only : solvation_type

   use ddx, only: ddinit, ddx_type, ddx_error_type, check_error, allocate_state, ddx_state_type


   implicit none
   private

   public :: new_ddx, ddx_input

   !> Input for ddx solvation
   type :: ddx_input
      !> Dielectric constant
      real(wp) :: dielectric_const
      !> Number of grid points for each atom
      integer :: nang = grid_size(6)
      !> Scaling of van-der-Waals radii
      real(wp) :: rscale = 1.0_wp
      !> Accuracy for iterative solver
      real(wp) :: conv = 1.0e-8_wp
      !> Regularization parameter
      real(wp) :: eta = 2.0_wp
      !> Maximum angular momentum of basis functions
      integer :: lmax = 6
      !> Van-der-Waals radii for all species
      real(wp), allocatable :: rvdw(:)
   end type ddx_input

   ! !> Definition of polarizable continuum model
   ! type, extends(solvation_type) :: ddx_solvation
   !    !> Dielectric function
   !    real(wp) :: keps
   !    !> Van-der-Waal radii for all atoms
   !    real(wp), allocatable :: rvdw(:)
   ! contains
   !    ! !> Update cache from container
   !    ! procedure :: update
   !    ! !> Return dependency on density
   !    ! procedure :: variable_info
   !    ! !> Get electric field energy
   !    ! procedure :: get_energy
   !    ! !> Get electric field potential
   !    ! procedure :: get_potential
   !    ! !> Get electric field gradient
   !    ! procedure :: get_gradient
   ! end type ddx_solvation

contains

!> Create new electric field container
subroutine new_ddx(mol)
   !! @param[in] model: 1 for COSMO, 2 for PCM and 3 for LPB
   !! @param[in] nsph: number of spheres n > 0
   !! @param[in] coords: coordinates of the spheres, size (3, nsph)
   !! @param[in] radii: radii of the spheres, size (nsph)
   !! @param[in] eps: relative dielectric permittivity eps > 1
   !! @param[out] ddx_data: Container for ddX data structures

   type(structure_type), intent(in) :: mol
   integer :: model
   real(wp) :: eps
   real(wp), allocatable :: coords(:,:)
   real(wp), allocatable :: radii(:)
   integer :: nsph
   type(ddx_type) :: ddx_solvation
   type(ddx_error_type) :: ddx_error

   type(ddx_state_type) :: ddx_state

   integer :: i
   

   write(*,*) 'new_ddx'

   model = 1
   nsph = mol%nat ! number of spheres = number of atoms, 3 for h2o
   allocate(coords(3, nsph), radii(nsph))
   do i = 1, nsph
      coords(:,i) = mol%xyz(:,i)
      radii(i) = get_vdw_rad_cosmo(mol%num(i))
   end do
   eps = 78.5_wp


   call ddinit(model, nsph, coords, radii, eps, ddx_solvation, ddx_error)
   call check_error(ddx_error)



   call allocate_state(ddx_solvation%params, ddx_solvation%constants,  &
      ddx_state, ddx_error)
   call check_error(ddx_error)






end subroutine new_ddx

end module tblite_solvation_ddx
