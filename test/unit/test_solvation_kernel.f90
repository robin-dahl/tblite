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

module test_solvation_kernel
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_solvation_born, only : born_integrator, new_born_integrator
   use tblite_solvation_data, only : get_vdw_rad_cosmo, get_vdw_rad_d3
   use tblite_solvation_kernel, only : kernel_type, new_kernel, kernel_enum
   use tblite_solvation_alpb, only : alpb_solvation, alpb_input, alpb_cache
   use tblite_solvation_data_alpb, only : get_alpb_param
   use tblite_solvation_data, only : solvent_data, get_vdw_rad_d3, get_solvent_data

   use mctc_io_structure, only: new_structure

   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_multipole, only : damped_multipole, new_damped_multipole
   use tblite_container_cache, only : container_cache
   implicit none
   private

   public :: collect_solvation_kernel

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   abstract interface
      subroutine multipole_maker(multipole, mol, error)
         import :: damped_multipole, structure_type, error_type
         type(damped_multipole), intent(out) :: multipole
         type(structure_type), intent(in) :: mol
         type(error_type), allocatable, intent(out) :: error
      end subroutine multipole_maker
   end interface

contains

!> Collect all exported unit tests
subroutine collect_solvation_kernel(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("amat-still", test_amat_still), &
      new_unittest("amat-p16", test_amat_p16), &
      new_unittest("amat-coulomb", test_amat_coulomb), &
      new_unittest("kernel-gradient-still", test_kernel_gradient_still), &
      new_unittest("kernel-gradient-born-still", test_kernel_dborn_still), &
      new_unittest("kernel-gradient-bornrad-still", test_kernel_gradient_dborn_still), &
      new_unittest("kernel-hessian-still", test_kernel_hessian_still), &
      new_unittest("kernel-hessian-bornrad-still", test_kernel_hessian_dborn_still), &
      new_unittest("kernel-third-still", test_kernel_third_still), &
      new_unittest("kernel-third-bornrad-still", test_kernel_third_dborn_still), &
      new_unittest("kernel-fourth-still", test_kernel_fourth_still), &
      new_unittest("kernel-fourth-bornrad-still", test_kernel_fourth_dborn_still), &
      new_unittest("kernel-fifth-still", test_kernel_fifth_still), &
      new_unittest("kernel-gradient-p16", test_kernel_gradient_p16), &
      new_unittest("kernel-gradient-bornrad-p16", test_kernel_gradient_dborn_p16), &
      new_unittest("kernel-hessian-p16", test_kernel_hessian_p16), &
      new_unittest("kernel-hessian-bornrad-p16", test_kernel_hessian_dborn_p16), &
      new_unittest("kernel-third-p16", test_kernel_third_p16), &
      new_unittest("kernel-third-bornrad-p16", test_kernel_third_dborn_p16), &
      new_unittest("kernel-fourth-p16", test_kernel_fourth_p16), &
      new_unittest("kernel-fourth-bornrad-p16", test_kernel_fourth_dborn_p16), &
      new_unittest("kernel-fifth-p16", test_kernel_fifth_p16), &
      new_unittest("kernel-gradient-coulomb", test_kernel_gradient_coulomb), &
      new_unittest("kernel-hessian-coulomb", test_kernel_hessian_coulomb), &
      new_unittest("kernel-third-coulomb", test_kernel_third_coulomb), &
      new_unittest("kernel-fourth-coulomb", test_kernel_fourth_coulomb), &
      new_unittest("kernel-fifth-coulomb", test_kernel_fifth_coulomb) &
      ]

end subroutine collect_solvation_kernel


!> Test Still kernel gradient against numerical derivative
subroutine test_kernel_gradient_still(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 0.5_wp
   
   call get_structure(mol, "MB16-43", "01")
   call test_numg(error, mol, kernel_enum%still, keps)

end subroutine test_kernel_gradient_still

!> Test gradient of Still kernel gradient wrt Born radii against numerical derivative
subroutine test_kernel_dborn_still(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 0.5_wp

   call get_structure(mol, "MB16-43", "01")
   call test_num_dborn_value(error, mol, kernel_enum%still, keps)

end subroutine test_kernel_dborn_still

!> Test gradient of Still kernel gradient wrt Born radii against numerical derivative
subroutine test_kernel_gradient_dborn_still(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 0.5_wp

   call get_structure(mol, "MB16-43", "01")
   call test_numg_dborn(error, mol, kernel_enum%still, keps)

end subroutine test_kernel_gradient_dborn_still

!> Test Still kernel Hessian against numerical derivative
subroutine test_kernel_hessian_still(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 0.5_wp
  
   call get_structure(mol, "MB16-43", "01")
   call test_numh(error, mol, kernel_enum%still, keps)

end subroutine test_kernel_hessian_still

!> Test gradient of Still kernel Hessian wrt Born radii against numerical derivative
subroutine test_kernel_hessian_dborn_still(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 0.5_wp
   
   call get_structure(mol, "MB16-43", "01")
   call test_numh_dborn(error, mol, kernel_enum%still, keps)

end subroutine test_kernel_hessian_dborn_still


!> Test Still kernel third derivative against numerical derivative
subroutine test_kernel_third_still(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 0.5_wp

   call get_structure(mol, "MB16-43", "01")
   call test_numt(error, mol, kernel_enum%still, keps)

end subroutine test_kernel_third_still


!> Test gradient of Still kernel third derivative wrt Born radii against numerical derivative
subroutine test_kernel_third_dborn_still(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 0.5_wp

   call get_structure(mol, "MB16-43", "01")
   call test_numt_dborn(error, mol, kernel_enum%still, keps)

end subroutine test_kernel_third_dborn_still


!> Test Still kernel fourth derivative against numerical derivative
subroutine test_kernel_fourth_still(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 0.5_wp

   call get_structure(mol, "MB16-43", "01")
   call test_numq(error, mol, kernel_enum%still, keps)

end subroutine test_kernel_fourth_still

!> Test gradient of Still kernel fourth derivative wrt Born radii against numerical derivative
subroutine test_kernel_fourth_dborn_still(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 0.5_wp

   call get_structure(mol, "MB16-43", "01")
   call test_numq_dborn(error, mol, kernel_enum%still, keps)

end subroutine test_kernel_fourth_dborn_still

!> Test Still kernel fifth derivative against numerical derivative
subroutine test_kernel_fifth_still(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 0.5_wp

   call get_structure(mol, "MB16-43", "01")
   call test_num5(error, mol, kernel_enum%still, keps)

end subroutine test_kernel_fifth_still


!> Test P16 kernel gradient against numerical derivative
subroutine test_kernel_gradient_p16(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: keps = 0.5_wp

   call get_structure(mol, "MB16-43", "01")
   call test_numg(error, mol, kernel_enum%p16, keps)

end subroutine test_kernel_gradient_p16

!> Test gradient of P16 kernel wrt Born radii against numerical derivative
subroutine test_kernel_gradient_dborn_p16(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: keps = 0.5_wp

   call get_structure(mol, "MB16-43", "01")
   call test_numg_dborn(error, mol, kernel_enum%p16, keps)

end subroutine test_kernel_gradient_dborn_p16


!> Test P16 kernel Hessian against numerical derivative
subroutine test_kernel_hessian_p16(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: keps = 0.5_wp

   call get_structure(mol, "MB16-43", "01")
   call test_numh(error, mol, kernel_enum%p16, keps)

end subroutine test_kernel_hessian_p16

!> Test gradient of P16 kernel Hessian wrt Born radii against numerical derivative
subroutine test_kernel_hessian_dborn_p16(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: keps = 0.5_wp

   call get_structure(mol, "MB16-43", "01")
   call test_numh_dborn(error, mol, kernel_enum%p16, keps)

end subroutine test_kernel_hessian_dborn_p16


!> Test P16 kernel third derivative against numerical derivative
subroutine test_kernel_third_p16(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: keps = 0.5_wp

   call get_structure(mol, "MB16-43", "01")
   call test_numt(error, mol, kernel_enum%p16, keps)

end subroutine test_kernel_third_p16

!> Test gradient of P16 kernel third derivative wrt Born radii against numerical derivative
subroutine test_kernel_third_dborn_p16(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: keps = 0.5_wp

   call get_structure(mol, "MB16-43", "01")
   call test_numt_dborn(error, mol, kernel_enum%p16, keps)

end subroutine test_kernel_third_dborn_p16

!> Test P16 kernel fourth derivative against numerical derivative
subroutine test_kernel_fourth_p16(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: keps = 0.5_wp

   call get_structure(mol, "MB16-43", "01")
   call test_numq(error, mol, kernel_enum%p16, keps)

end subroutine test_kernel_fourth_p16

!> Test gradient of P16 kernel fourth derivative wrt Born radii against numerical derivative
subroutine test_kernel_fourth_dborn_p16(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: keps = 0.5_wp

   call get_structure(mol, "MB16-43", "01")
   call test_numq_dborn(error, mol, kernel_enum%p16, keps)

end subroutine test_kernel_fourth_dborn_p16

!> Test P16 kernel fifth derivative against numerical derivative
subroutine test_kernel_fifth_p16(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: keps = 0.5_wp

   call get_structure(mol, "MB16-43", "01")
   call test_num5(error, mol, kernel_enum%p16, keps)

end subroutine test_kernel_fifth_p16

!> Test Coulomb kernel gradient against numerical derivative
subroutine test_kernel_gradient_coulomb(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 1.0_wp

   call get_structure(mol, "MB16-43", "01")
   call test_numg(error, mol, kernel_enum%coulomb, keps)

end subroutine test_kernel_gradient_coulomb

!> Test Coulomb kernel Hessian against numerical derivative
subroutine test_kernel_hessian_coulomb(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 1.0_wp

   call get_structure(mol, "MB16-43", "01")
   call test_numh(error, mol, kernel_enum%coulomb, keps)

end subroutine test_kernel_hessian_coulomb

!> Test Coulomb kernel third derivative against numerical derivative
subroutine test_kernel_third_coulomb(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 1.0_wp

   call get_structure(mol, "MB16-43", "01")
   call test_numt(error, mol, kernel_enum%coulomb, keps)

end subroutine test_kernel_third_coulomb

!> Test Coulomb kernel fourth derivative against numerical derivative
subroutine test_kernel_fourth_coulomb(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 1.0_wp

   call get_structure(mol, "MB16-43", "01")
   call test_numq(error, mol, kernel_enum%coulomb, keps)

end subroutine test_kernel_fourth_coulomb

!> Test Coulomb kernel fifth derivative against numerical derivative
subroutine test_kernel_fifth_coulomb(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 1.0_wp

   call get_structure(mol, "MB16-43", "01")
   call test_num5(error, mol, kernel_enum%coulomb, keps)

end subroutine test_kernel_fifth_coulomb

!> Test construction of multipole interaction matrices
subroutine test_amat_still(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(alpb_input) :: input
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 0.5_wp

   real(wp), allocatable :: rad(:), ds(:)

   real(wp) :: amat_sd_ref(3,3,3), amat_dd_ref(3,3,3,3), &
      & amat_sq_ref(6,3,3), amat_dq_ref(3,3,6,3), amat_qq_ref(6,3,6,3)

   call get_amat_ref_still(amat_sd_ref, amat_dd_ref, &
      & amat_sq_ref, amat_dq_ref, amat_qq_ref)

   solvent = get_solvent_data("water")
   input = alpb_input(solvent%eps, solvent=solvent%solvent, &
         & kernel=1, alpb=.true., do_multipoles=.true.)
   call get_structure(mol, "MB16-43", "BeH2")
   call test_amat(error, mol, keps, input, &
      & amat_sd_ref, amat_dd_ref, amat_sq_ref, amat_dq_ref, amat_qq_ref) 

end subroutine test_amat_still

!> Test construction of multipole interaction matrices
subroutine test_amat_p16(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(alpb_input) :: input
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 0.5_wp

   real(wp), allocatable :: rad(:), ds(:)

   real(wp) :: amat_sd_ref(3,3,3), amat_dd_ref(3,3,3,3), &
      & amat_sq_ref(6,3,3), amat_dq_ref(3,3,6,3), amat_qq_ref(6,3,6,3)

   call get_amat_ref_p16(amat_sd_ref, amat_dd_ref, &
      & amat_sq_ref, amat_dq_ref, amat_qq_ref)

   solvent = get_solvent_data("water")
   input = alpb_input(solvent%eps, solvent=solvent%solvent, &
         & kernel=2, alpb=.true., do_multipoles=.true.)
   call get_structure(mol, "MB16-43", "BeH2")
   call test_amat(error, mol, keps, input, &
      & amat_sd_ref, amat_dd_ref, amat_sq_ref, amat_dq_ref, amat_qq_ref) 

end subroutine test_amat_p16

!> Test construction of multipole interaction matrices
subroutine test_amat_coulomb(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(alpb_input) :: input
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 0.5_wp

   real(wp), allocatable :: rad(:), ds(:)

   real(wp) :: amat_sd_ref(3,3,3), amat_dd_ref(3,3,3,3), &
      & amat_sq_ref(6,3,3), amat_dq_ref(3,3,6,3), amat_qq_ref(6,3,6,3)

   call get_amat_ref_coulomb(amat_sd_ref, amat_dd_ref, &
      & amat_sq_ref, amat_dq_ref, amat_qq_ref)

   solvent = get_solvent_data("water")
   input = alpb_input(solvent%eps, solvent=solvent%solvent, &
         & kernel=3, alpb=.true., do_multipoles=.true.)
   call get_structure(mol, "MB16-43", "BeH2")
   call test_amat(error, mol, keps, input, &
      & amat_sd_ref, amat_dd_ref, amat_sq_ref, amat_dq_ref, amat_qq_ref) 

end subroutine test_amat_coulomb


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> Test the kernel gradient against numerical derivative
subroutine test_numg(error, mol, kernel_id, keps)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data (modified during finite difference)
   type(structure_type), intent(inout) :: mol
   !> Kernel identifier (still, p16, or coulomb)
   integer, intent(in) :: kernel_id
   !> Dielectric screening factor
   real(wp), intent(in) :: keps

   !> Born radii integrator
   type(born_integrator) :: gbobc

   !> Kernel instance
   class(kernel_type), allocatable :: kernel
   !> Van der Waals radii for all atoms
   real(wp), allocatable :: rvdw(:)
   !> Born radii for all atoms
   real(wp), allocatable :: brad(:)
   !> Kernel matrix with positive displacement
   real(wp), allocatable :: kernel_r(:, :)
   !> Kernel matrix with negative displacement
   real(wp), allocatable :: kernel_l(:, :)
   !> Numerical gradient: ∂K_mn/∂R_k,α (3, nat, nat, nat)
   real(wp), allocatable :: numg_kernel(:, :, :, :)
   !> Analytical gradient for a single pair (3, nat)
   real(wp), allocatable :: anag_pair_kernel(:, :)

   !> Finite difference step size
   real(wp), parameter :: step = 1.0e-6_wp
   !> Loop indices for atoms and Cartesian directions
   integer :: iat, jat, ic, jc, k, alpha
   !> Maximum difference between analytical and numerical gradients
   real(wp) :: maxdiff
   !> Individual gradient element
   real(wp) :: dKdr_elem
   !> Location of maximum difference
   integer :: loc(2)

   !> Gradient for atom pair interaction
   real(wp) :: grad_m(3)

   kernel = new_kernel(kernel_id, keps)
   
   allocate(rvdw(mol%nat), brad(mol%nat))
   allocate(kernel_r(mol%nat, mol%nat), kernel_l(mol%nat, mol%nat))
   allocate(numg_kernel(3, mol%nat, mol%nat, mol%nat))
   allocate(anag_pair_kernel(3, mol%nat))

   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   call gbobc%get_rad(mol, brad)

   ! Numerical derivative
   numg_kernel(:, :, :, :) = 0.0_wp

   do iat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         kernel_r(:, :) = 0.0_wp
         call kernel%add_kernel_mat(mol%nat, mol%xyz, brad, kernel_r)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         kernel_l(:, :) = 0.0_wp
         call kernel%add_kernel_mat(mol%nat, mol%xyz, brad, kernel_l)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step

         do jat = 1, mol%nat
            do jc = 1, mol%nat
               numg_kernel(ic, iat, jat, jc) = 0.5_wp * (kernel_r(jat, jc) - kernel_l(jat, jc)) / step
            end do
         end do
      end do
   end do

   ! Compare for each (m,n): build analytical dK_mn/dR_k,alpha element-by-element
   do jat = 1, mol%nat
      do jc = 1, mol%nat

         anag_pair_kernel(:, :) = 0.0_wp

         call kernel%kernel_d1_pair(mol%xyz(:, jat), mol%xyz(:, jc), brad(jat), brad(jc), grad_m)

         do k = 1, mol%nat
            do alpha = 1, 3
               dKdr_elem = 0.0_wp

               if (jat /= jc) then
                  if (k == jat) then
                     dKdr_elem = grad_m(alpha)
                  else if (k == jc) then
                     dKdr_elem = -grad_m(alpha)
                  end if
               end if

               anag_pair_kernel(alpha, k) = dKdr_elem
            end do
         end do

         maxdiff = maxval(abs(anag_pair_kernel(:, :) - numg_kernel(:, :, jat, jc)))
         if (maxdiff > thr2) then
            loc = maxloc(abs(anag_pair_kernel(:, :) - numg_kernel(:, :, jat, jc)))
            call test_failed(error, "Analytical gradient does not match finite difference solution")
            print '(a,2i6, a, es20.13)', "Mismatch at (m,n)=(", jat, jc, "), max|diff|=", maxdiff
            print '(a,i2,a,i6)', "Worst entry at alpha=", loc(1), " k=", loc(2)
            print '(a,3es20.13)', "anag_pair_kernel(:,k)=", anag_pair_kernel(:, loc(2))
            print '(a,3es20.13)', "numg   (:,k)=", numg_kernel(:, loc(2), jat, jc)
            print '(a,3es20.13)', "diff   (:,k)=", anag_pair_kernel(:, loc(2)) - numg_kernel(:, loc(2), jat, jc)
            return
         end if

      end do
   end do

end subroutine test_numg


!> Test kernel Born-radius gradient (∂K/∂born) against numerical derivative
!> Test kernel Born-radius gradient (∂K/∂born) against numerical derivative
subroutine test_num_dborn_value(error, mol, kernel_id, keps)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol
   !> Kernel identifier (still, p16, or coulomb)
   integer, intent(in) :: kernel_id
   !> Dielectric screening factor
   real(wp), intent(in) :: keps

   !> Born radii integrator
   type(born_integrator) :: gbobc
   !> Kernel instance
   class(kernel_type), allocatable :: kernel

   !> Van der Waals radii for all atoms
   real(wp), allocatable :: rvdw(:)
   !> Born radii for all atoms
   real(wp), allocatable :: brad(:)
   !> Temporary Born radii (for FD)
   real(wp), allocatable :: brtmp(:)

   !> Kernel matrix with positive/negative Born displacement
   real(wp), allocatable :: kernel_r(:, :)
   real(wp), allocatable :: kernel_l(:, :)

   !> Finite difference step size
   real(wp), parameter :: step = 1.0e-6_wp

   integer :: jat, jc
   real(wp) :: bornA0, bornB0
   real(wp) :: num_dK_bA, num_dK_bB
   real(wp) :: ana_dK_bA, ana_dK_bB
   real(wp) :: diffA, diffB, maxdiff

   kernel = new_kernel(kernel_id, keps)

   allocate(rvdw(mol%nat), brad(mol%nat), brtmp(mol%nat))
   allocate(kernel_r(mol%nat, mol%nat), kernel_l(mol%nat, mol%nat))

   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   call gbobc%get_rad(mol, brad)

   do jat = 1, mol%nat
      do jc = 1, mol%nat

         ! Same convention as your dborn test: skip self-pair
         if (jat == jc) cycle

         bornA0 = brad(jat)
         bornB0 = brad(jc)

         ! ---------------------------
         ! Numerical: ∂K(jat,jc)/∂bornA  (perturb brad(jat))
         ! ---------------------------
         brtmp(:) = brad(:)

         brtmp(jat) = bornA0 + step
         kernel_r(:, :) = 0.0_wp
         call kernel%add_kernel_mat(mol%nat, mol%xyz, brtmp, kernel_r)

         brtmp(jat) = bornA0 - step
         kernel_l(:, :) = 0.0_wp
         call kernel%add_kernel_mat(mol%nat, mol%xyz, brtmp, kernel_l)

         num_dK_bA = 0.5_wp * (kernel_r(jat, jc) - kernel_l(jat, jc)) / step

         ! ---------------------------
         ! Numerical: ∂K(jat,jc)/∂bornB  (perturb brad(jc))
         ! ---------------------------
         brtmp(:) = brad(:)

         brtmp(jc) = bornB0 + step
         kernel_r(:, :) = 0.0_wp
         call kernel%add_kernel_mat(mol%nat, mol%xyz, brtmp, kernel_r)

         brtmp(jc) = bornB0 - step
         kernel_l(:, :) = 0.0_wp
         call kernel%add_kernel_mat(mol%nat, mol%xyz, brtmp, kernel_l)

         num_dK_bB = 0.5_wp * (kernel_r(jat, jc) - kernel_l(jat, jc)) / step

         ! ---------------------------
         ! Analytical: call kernel_pair_dborn
         ! ---------------------------
         call kernel%kernel_pair_dborn( mol%xyz(:, jat), mol%xyz(:, jc), bornA0, bornB0, &
                                        ana_dK_bA, ana_dK_bB )

         diffA = ana_dK_bA - num_dK_bA
         diffB = ana_dK_bB - num_dK_bB
         maxdiff = max(abs(diffA), abs(diffB))

         if (maxdiff > thr2) then
            call test_failed(error, "Analytical kernel Born-gradient does not match finite difference solution")
            print '(a,2i6,a,es20.13)', "Mismatch at (A,B)=(", jat, jc, "), max|diff|=", maxdiff
            print '(a,es20.13,a,es20.13)', "ana_dK_bA=", ana_dK_bA, "  num_dK_bA=", num_dK_bA
            print '(a,es20.13,a,es20.13)', "ana_dK_bB=", ana_dK_bB, "  num_dK_bB=", num_dK_bB
            print '(a,es20.13,a,es20.13)', "diffA=", diffA, "  diffB=", diffB
            return
         end if

      end do
   end do

end subroutine test_num_dborn_value




!> Test kernel gradient Born radius derivative against numerical derivative
subroutine test_numg_dborn(error, mol, kernel_id, keps)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol
   !> Kernel identifier (still, p16, or coulomb)
   integer, intent(in) :: kernel_id
   !> Dielectric screening factor
   real(wp), intent(in) :: keps

   !> Born radii integrator
   type(born_integrator) :: gbobc
   !> Kernel instance
   class(kernel_type), allocatable :: kernel

   !> Van der Waals radii for all atoms
   real(wp), allocatable :: rvdw(:)
   !> Born radii for all atoms
   real(wp), allocatable :: brad(:)

   !> Finite difference step size
   real(wp), parameter :: step = 1.0e-6_wp
   !> Loop indices for atom pairs
   integer :: jat, jc
   !> Maximum difference for Born radius derivatives
   real(wp) :: maxdiffA, maxdiffB
   !> Location of maximum difference
   integer :: locA(1), locB(1)

   !> Original Born radii before perturbation
   real(wp) :: bornA0, bornB0
   !> Gradient with positive/negative Born radius displacement
   real(wp) :: grad_r(3), grad_l(3)
   !> Numerical Born radius derivatives
   real(wp) :: num_bA(3), num_bB(3)

   !> Analytical Born radius derivatives from kernel_d1_pair_dborn
   real(wp) :: ana_bA(3), ana_bB(3)

   kernel = new_kernel(kernel_id, keps)

   allocate(rvdw(mol%nat), brad(mol%nat))
   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   call gbobc%get_rad(mol, brad)

   ! Loop over all ordered pairs (jat, jc) like your original test.
   do jat = 1, mol%nat
      do jc = 1, mol%nat

         ! Self-pair has r=0 => d1_pair = 0 by construction; dborn should also be 0.
         if (jat == jc) cycle

         bornA0 = brad(jat)
         bornB0 = brad(jc)

         ! ---------------------------
         ! Numerical: d(d1_pair)/d(bornA)
         ! ---------------------------
         brad(jat) = bornA0 + step
         call kernel%kernel_d1_pair(mol%xyz(:, jat), mol%xyz(:, jc), brad(jat), brad(jc), grad_r)

         brad(jat) = bornA0 - step
         call kernel%kernel_d1_pair(mol%xyz(:, jat), mol%xyz(:, jc), brad(jat), brad(jc), grad_l)

         brad(jat) = bornA0
         num_bA(:) = 0.5_wp * (grad_r(:) - grad_l(:)) / step

         ! ---------------------------
         ! Numerical: d(d1_pair)/d(bornB)
         ! ---------------------------
         brad(jc) = bornB0 + step
         call kernel%kernel_d1_pair(mol%xyz(:, jat), mol%xyz(:, jc), brad(jat), brad(jc), grad_r)

         brad(jc) = bornB0 - step
         call kernel%kernel_d1_pair(mol%xyz(:, jat), mol%xyz(:, jc), brad(jat), brad(jc), grad_l)

         brad(jc) = bornB0
         num_bB(:) = 0.5_wp * (grad_r(:) - grad_l(:)) / step

         ! ---------------------------
         ! Analytical: kernel_d1_pair_dborn
         ! ---------------------------
         call kernel%kernel_d1_pair_dborn( mol%xyz(:, jat), mol%xyz(:, jc), bornA0, bornB0, &
                                           ana_bA, ana_bB )

         ! ---------------------------
         ! Compare + diagnostics
         ! ---------------------------
         maxdiffA = maxval(abs(ana_bA(:) - num_bA(:)))
         maxdiffB = maxval(abs(ana_bB(:) - num_bB(:)))

         if (maxdiffA > thr2 .or. maxdiffB > thr2) then
            call test_failed(error, "Analytical d1_pair Born-derivative does not match finite difference solution")

            if (maxdiffA >= maxdiffB) then
               locA = maxloc(abs(ana_bA(:) - num_bA(:)))
               print '(a,2i6,a,es20.13)', "Mismatch d(d1)/d(bornA) at (A,B)=(", jat, jc, "), max|diff|=", maxdiffA
               print '(a,i2)', "Worst component alpha=", locA(1)
               print '(a,3es20.13)', "ana_bA =", ana_bA(:)
               print '(a,3es20.13)', "num_bA =", num_bA(:)
               print '(a,3es20.13)', "diff  =", ana_bA(:) - num_bA(:)
            else
               locB = maxloc(abs(ana_bB(:) - num_bB(:)))
               print '(a,2i6,a,es20.13)', "Mismatch d(d1)/d(bornB) at (A,B)=(", jat, jc, "), max|diff|=", maxdiffB
               print '(a,i2)', "Worst component alpha=", locB(1)
               print '(a,3es20.13)', "ana_bB =", ana_bB(:)
               print '(a,3es20.13)', "num_bB =", num_bB(:)
               print '(a,3es20.13)', "diff  =", ana_bB(:) - num_bB(:)
            end if

            return
         end if

      end do
   end do

end subroutine test_numg_dborn


!> Test the kernel Hessian against numerical derivative
subroutine test_numh(error, mol, kernel_id, keps)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol
   !> Kernel identifier (still, p16, or coulomb)
   integer, intent(in) :: kernel_id
   !> Dielectric screening factor
   real(wp), intent(in) :: keps

   !> Kernel instance
   class(kernel_type), allocatable :: kernel

   !> Born radii integrator
   type(born_integrator) :: gbobc

   !> Finite difference step size
   real(wp), parameter :: step = 1.0e-6_wp

   !> Loop indices for atom pairs and Cartesian directions
   integer :: jat, jc, i, j
   !> Van der Waals radii for all atoms
   real(wp), allocatable  :: rvdw(:)
   !> Born radii for all atoms
   real(wp), allocatable :: brad(:)
   !> Atomic coordinates and perturbed coordinates
   real(wp) :: rA(3), rB(3), rA_r(3), rA_l(3)
   !> Interatomic distance
   real(wp) :: r
   !> Gradient with positive/negative displacement
   real(wp) :: dkernel_l(3), dkernel_r(3)
   !> Analytical and numerical Hessian matrices
   real(wp) :: anah_kernel(3,3), numh_kernel(3,3)
   !> Maximum difference between analytical and numerical Hessian
   real(wp) :: maxdiff
   !> Location of maximum difference
   integer :: loc(2)

   kernel = new_kernel(kernel_id, keps)

   allocate(rvdw(mol%nat), brad(mol%nat))
   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   call gbobc%get_rad(mol, brad)

   do jat = 1, mol%nat
      do jc = 1, mol%nat
         if (jat == jc) cycle

         rA = mol%xyz(:, jat)
         rB = mol%xyz(:, jc)
         r  = norm2(rA - rB)

         ! analytic Hessian
         call kernel%kernel_d2_pair(rA, rB, brad(jat), brad(jc), anah_kernel)

         ! numerical Hessian from analytic gradient:
         ! n2(i,j) = d/dR_A,j [ d1(i) ]
         numh_kernel(:,:) = 0.0_wp
         do j = 1, 3
            rA_r = rA; rA_r(j) = rA_r(j) + step
            rA_l = rA; rA_l(j) = rA_l(j) - step

            call kernel%kernel_d1_pair(rA_r, rB, brad(jat), brad(jc), dkernel_r)
            call kernel%kernel_d1_pair(rA_l, rB, brad(jat), brad(jc), dkernel_l)

            do i = 1, 3
               numh_kernel(i,j) = 0.5_wp * (dkernel_r(i) - dkernel_l(i)) / step
            end do
         end do

         maxdiff = maxval(abs(anah_kernel - numh_kernel))
         if (maxdiff > thr2) then
            loc = maxloc(abs(anah_kernel - numh_kernel))
            call test_failed(error, "Analytical Hessian does not match finite difference solution")
            print '(a,2i6,a,es20.13)', "Mismatch at pair (A,B)=(", jat, jc, "), max|diff|=", maxdiff
            print '(a,2i2)', "Worst entry (i,j)=", loc(1), loc(2)
            print '(a,es20.13)', "anah_kernel(i,j)=", anah_kernel(loc(1),loc(2))
            print '(a,es20.13)', "numh_kernel(i,j)=", numh_kernel(loc(1),loc(2))
            print '(a,es20.13)', "diff        =", anah_kernel(loc(1),loc(2)) - numh_kernel(loc(1),loc(2))
            return
         end if

      end do
   end do
end subroutine test_numh

!> Test kernel_d2_pair_dborn against numerical derivative of kernel_d2_pair wrt Born radii
subroutine test_numh_dborn(error, mol, kernel_id, keps)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol
   !> Kernel identifier (still, p16, or coulomb)
   integer, intent(in) :: kernel_id
   !> Dielectric screening factor
   real(wp), intent(in) :: keps

   !> Kernel instance
   class(kernel_type), allocatable :: kernel
   !> Born radii integrator
   type(born_integrator) :: gbobc

   !> Finite difference step size
   real(wp), parameter :: step = 1.0e-6_wp

   !> Loop indices for atom pairs
   integer :: jat, jc
   !> Van der Waals radii for all atoms
   real(wp), allocatable :: rvdw(:)
   !> Born radii for all atoms
   real(wp), allocatable :: brad(:)

   !> Atomic coordinates
   real(wp) :: rA(3), rB(3)
   !> Original Born radii before perturbation
   real(wp) :: bornA0, bornB0

   !> Hessian with positive/negative Born radius displacement
   real(wp) :: d2_r(3,3), d2_l(3,3)
   !> Numerical Born radius derivatives of Hessian
   real(wp) :: num_bA(3,3), num_bB(3,3)
   !> Analytical Born radius derivatives from kernel_d2_pair_dborn
   real(wp) :: ana_bA(3,3), ana_bB(3,3)

   !> Maximum differences for Born radius derivatives
   real(wp) :: maxdiffA, maxdiffB
   !> Locations of maximum differences
   integer :: locA(2), locB(2)

   kernel = new_kernel(kernel_id, keps)

   allocate(rvdw(mol%nat), brad(mol%nat))
   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   call gbobc%get_rad(mol, brad)

   do jat = 1, mol%nat
      do jc = 1, mol%nat
         if (jat == jc) cycle

         rA = mol%xyz(:, jat)
         rB = mol%xyz(:, jc)

         bornA0 = brad(jat)
         bornB0 = brad(jc)

         ! ---- analytic d/d(bornA), d/d(bornB)
         call kernel%kernel_d2_pair_dborn(rA, rB, bornA0, bornB0, ana_bA, ana_bB)

         ! ---- numerical d/d(bornA) of d2_pair
         call kernel%kernel_d2_pair(rA, rB, bornA0 + step, bornB0, d2_r)
         call kernel%kernel_d2_pair(rA, rB, bornA0 - step, bornB0, d2_l)
         num_bA(:,:) = 0.5_wp * (d2_r(:,:) - d2_l(:,:)) / step

         ! ---- numerical d/d(bornB) of d2_pair
         call kernel%kernel_d2_pair(rA, rB, bornA0, bornB0 + step, d2_r)
         call kernel%kernel_d2_pair(rA, rB, bornA0, bornB0 - step, d2_l)
         num_bB(:,:) = 0.5_wp * (d2_r(:,:) - d2_l(:,:)) / step

         maxdiffA = maxval(abs(ana_bA - num_bA))
         maxdiffB = maxval(abs(ana_bB - num_bB))

         if (maxdiffA > thr2 .or. maxdiffB > thr2) then
            call test_failed(error, "kernel_d2_pair_dborn does not match FD(d2_pair)")

            if (maxdiffA >= maxdiffB) then
               locA = maxloc(abs(ana_bA - num_bA))
               print '(a,2i6,a,es20.13)', "Mismatch d(d2)/d(bornA) at pair (A,B)=(", jat, jc, "), max|diff|=", maxdiffA
               print '(a,2i2)', "Worst entry (i,j)=", locA(1), locA(2)
               print '(a,es20.13)', "ana_bA(i,j)=", ana_bA(locA(1),locA(2))
               print '(a,es20.13)', "num_bA(i,j)=", num_bA(locA(1),locA(2))
               print '(a,es20.13)', "diff       =", ana_bA(locA(1),locA(2)) - num_bA(locA(1),locA(2))
            else
               locB = maxloc(abs(ana_bB - num_bB))
               print '(a,2i6,a,es20.13)', "Mismatch d(d2)/d(bornB) at pair (A,B)=(", jat, jc, "), max|diff|=", maxdiffB
               print '(a,2i2)', "Worst entry (i,j)=", locB(1), locB(2)
               print '(a,es20.13)', "ana_bB(i,j)=", ana_bB(locB(1),locB(2))
               print '(a,es20.13)', "num_bB(i,j)=", num_bB(locB(1),locB(2))
               print '(a,es20.13)', "diff       =", ana_bB(locB(1),locB(2)) - num_bB(locB(1),locB(2))
            end if

            return
         end if

      end do
   end do
end subroutine test_numh_dborn


!> Test the kernel third derivative against numerical derivative
subroutine test_numt(error, mol, kernel_id, keps)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol
   !> Kernel identifier (still, p16, or coulomb)
   integer, intent(in) :: kernel_id
   !> Dielectric screening factor
   real(wp), intent(in) :: keps
   !> Kernel instance
   class(kernel_type), allocatable :: kernel

   !> Born radii integrator
   type(born_integrator) :: gbobc

   !> Finite difference step size
   real(wp), parameter :: step = 1.0e-6_wp

   !> Loop indices for atom pairs and Cartesian directions
   integer :: jat, jc, i, j, k
   !> Van der Waals radii for all atoms
   real(wp), allocatable  :: rvdw(:)
   !> Born radii for all atoms
   real(wp), allocatable :: brad(:)
   !> Atomic coordinates and perturbed coordinates
   real(wp) :: rA(3), rB(3), rA_r(3), rA_l(3)
   !> Interatomic distance
   real(wp) :: r
   !> Hessian with positive/negative displacement
   real(wp) :: d2kernel_r(3,3), d2kernel_l(3,3)
   !> Analytical and numerical third derivative tensors
   real(wp) :: anat_kernel(3,3,3), numt_kernel(3,3,3)
   !> Maximum difference between analytical and numerical third derivatives
   real(wp) :: maxdiff
   !> Location of maximum difference
   integer :: loc(3)

   kernel = new_kernel(kernel_id, keps)

   allocate(rvdw(mol%nat), brad(mol%nat))
   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   call gbobc%get_rad(mol, brad)

   do jat = 1, mol%nat
      do jc = 1, mol%nat
         if (jat == jc) cycle

         rA = mol%xyz(:, jat)
         rB = mol%xyz(:, jc)
         r  = norm2(rA - rB)

         ! analytic 3rd derivative
         call kernel%kernel_d3_pair(rA, rB, brad(jat), brad(jc), anat_kernel)

         ! numerical 3rd derivative from analytic Hessian:
         ! n3(i,j,k) = d/dR_A,k [ d2(i,j) ]
         numt_kernel(:,:,:) = 0.0_wp
         do k = 1, 3
            rA_r = rA; rA_r(k) = rA_r(k) + step
            rA_l = rA; rA_l(k) = rA_l(k) - step

            call kernel%kernel_d2_pair(rA_r, rB, brad(jat), brad(jc), d2kernel_r)
            call kernel%kernel_d2_pair(rA_l, rB, brad(jat), brad(jc), d2kernel_l)

            do i = 1, 3
               do j = 1, 3
                  numt_kernel(i,j,k) = 0.5_wp * (d2kernel_r(i,j) - d2kernel_l(i,j)) / step
               end do
            end do
         end do

         maxdiff = maxval(abs(anat_kernel - numt_kernel))
         if (maxdiff > thr2) then
            loc = maxloc(abs(anat_kernel - numt_kernel))
            call test_failed(error, "Coulomb d3 (pair) does not match FD(d2)")
            print '(a,2i6,a,es20.13)', "Mismatch at pair (A,B)=(", jat, jc, "), max|diff|=", maxdiff
            print '(a,3i2)', "Worst entry (i,j,k)=", loc(1), loc(2), loc(3)
            print '(a,es20.13)', "ana d3(i,j,k)=", anat_kernel(loc(1),loc(2),loc(3))
            print '(a,es20.13)', "num d3(i,j,k)=", numt_kernel(loc(1),loc(2),loc(3))
            print '(a,es20.13)', "diff          =", anat_kernel(loc(1),loc(2),loc(3)) - numt_kernel(loc(1),loc(2),loc(3))
            return
         end if

      end do
   end do
end subroutine test_numt

!> Test kernel_d3_pair_dborn against numerical derivative of kernel_d3_pair wrt Born radii
subroutine test_numt_dborn(error, mol, kernel_id, keps)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol
   !> Kernel identifier (still, p16, or coulomb)
   integer, intent(in) :: kernel_id
   !> Dielectric screening factor
   real(wp), intent(in) :: keps

   !> Kernel instance
   class(kernel_type), allocatable :: kernel
   !> Born radii integrator
   type(born_integrator) :: gbobc

   !> Finite difference step size
   real(wp), parameter :: step = 1.0e-6_wp

   !> Loop indices for atom pairs
   integer :: jat, jc
   !> Van der Waals radii for all atoms
   real(wp), allocatable :: rvdw(:)
   !> Born radii for all atoms
   real(wp), allocatable :: brad(:)

   !> Atomic coordinates
   real(wp) :: rA(3), rB(3)
   !> Original Born radii before perturbation
   real(wp) :: bornA0, bornB0

   !> Third derivative with positive/negative Born radius displacement
   real(wp) :: d3_r(3,3,3), d3_l(3,3,3)
   !> Numerical Born radius derivatives of third derivative
   real(wp) :: num_bA(3,3,3), num_bB(3,3,3)
   !> Analytical Born radius derivatives from kernel_d3_pair_dborn
   real(wp) :: ana_bA(3,3,3), ana_bB(3,3,3)

   !> Maximum differences for Born radius derivatives
   real(wp) :: maxdiffA, maxdiffB
   !> Locations of maximum differences
   integer :: locA(3), locB(3)

   kernel = new_kernel(kernel_id, keps)

   allocate(rvdw(mol%nat), brad(mol%nat))
   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   call gbobc%get_rad(mol, brad)

   do jat = 1, mol%nat
      do jc = 1, mol%nat
         if (jat == jc) cycle

         rA = mol%xyz(:, jat)
         rB = mol%xyz(:, jc)

         bornA0 = brad(jat)
         bornB0 = brad(jc)

         ! ---- analytic d/d(bornA), d/d(bornB)
         call kernel%kernel_d3_pair_dborn(rA, rB, bornA0, bornB0, ana_bA, ana_bB)

         ! ---- numerical d/d(bornA) of d3_pair
         call kernel%kernel_d3_pair(rA, rB, bornA0 + step, bornB0, d3_r)
         call kernel%kernel_d3_pair(rA, rB, bornA0 - step, bornB0, d3_l)
         num_bA(:,:,:) = 0.5_wp * (d3_r(:,:,:) - d3_l(:,:,:)) / step

         ! ---- numerical d/d(bornB) of d3_pair
         call kernel%kernel_d3_pair(rA, rB, bornA0, bornB0 + step, d3_r)
         call kernel%kernel_d3_pair(rA, rB, bornA0, bornB0 - step, d3_l)
         num_bB(:,:,:) = 0.5_wp * (d3_r(:,:,:) - d3_l(:,:,:)) / step

         maxdiffA = maxval(abs(ana_bA - num_bA))
         maxdiffB = maxval(abs(ana_bB - num_bB))

         if (maxdiffA > thr2 .or. maxdiffB > thr2) then
            call test_failed(error, "kernel_d3_pair_dborn does not match FD(d3_pair)")

            if (maxdiffA >= maxdiffB) then
               locA = maxloc(abs(ana_bA - num_bA))
               print '(a,2i6,a,es20.13)', "Mismatch d(d3)/d(bornA) at pair (A,B)=(", jat, jc, "), max|diff|=", maxdiffA
               print '(a,3i2)', "Worst entry (i,j,k)=", locA(1), locA(2), locA(3)
               print '(a,es20.13)', "ana_bA(i,j,k)=", ana_bA(locA(1),locA(2),locA(3))
               print '(a,es20.13)', "num_bA(i,j,k)=", num_bA(locA(1),locA(2),locA(3))
               print '(a,es20.13)', "diff         =", ana_bA(locA(1),locA(2),locA(3)) - num_bA(locA(1),locA(2),locA(3))
            else
               locB = maxloc(abs(ana_bB - num_bB))
               print '(a,2i6,a,es20.13)', "Mismatch d(d3)/d(bornB) at pair (A,B)=(", jat, jc, "), max|diff|=", maxdiffB
               print '(a,3i2)', "Worst entry (i,j,k)=", locB(1), locB(2), locB(3)
               print '(a,es20.13)', "ana_bB(i,j,k)=", ana_bB(locB(1),locB(2),locB(3))
               print '(a,es20.13)', "num_bB(i,j,k)=", num_bB(locB(1),locB(2),locB(3))
               print '(a,es20.13)', "diff         =", ana_bB(locB(1),locB(2),locB(3)) - num_bB(locB(1),locB(2),locB(3))
            end if

            return
         end if

      end do
   end do
end subroutine test_numt_dborn


!> Test the kernel fourth derivative against numerical derivative
subroutine test_numq(error, mol, kernel_id, keps)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol
   !> Kernel identifier (still, p16, or coulomb)
   integer, intent(in) :: kernel_id
   !> Dielectric screening factor
   real(wp), intent(in) :: keps

   !> Kernel instance
   class(kernel_type), allocatable :: kernel

   !> Born radii integrator
   type(born_integrator) :: gbobc

   !> Finite difference step size
   real(wp), parameter :: step = 1.0e-6_wp

   !> Loop indices for atom pairs and Cartesian directions
   integer :: jat, jc, i, j, k, l
   !> Van der Waals radii for all atoms
   real(wp), allocatable :: rvdw(:)
   !> Born radii for all atoms
   real(wp), allocatable :: brad(:)
   !> Atomic coordinates and perturbed coordinates
   real(wp) :: rA(3), rB(3), rA_r(3), rA_l(3)
   !> Interatomic distance
   real(wp) :: r
   !> Third derivative with positive/negative displacement
   real(wp) :: d3kernel_r(3,3,3), d3kernel_l(3,3,3)
   !> Analytical and numerical fourth derivative tensors
   real(wp) :: anaq_kernel(3,3,3,3), numq_kernel(3,3,3,3)
   !> Maximum difference between analytical and numerical fourth derivatives
   real(wp) :: maxdiff
   !> Location of maximum difference
   integer :: loc(4)

   kernel = new_kernel(kernel_id, keps)

   allocate(rvdw(mol%nat), brad(mol%nat))
   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   call gbobc%get_rad(mol, brad)

   do jat = 1, mol%nat
      do jc = 1, mol%nat
         if (jat == jc) cycle

         rA = mol%xyz(:, jat)
         rB = mol%xyz(:, jc)
         r  = norm2(rA - rB)

         ! analytic 4th derivative
         call kernel%kernel_d4_pair(rA, rB, brad(jat), brad(jc), anaq_kernel)

         ! numerical 4th derivative from analytic 3rd derivative:
         ! n4(i,j,k,l) = d/dR_A,l [ d3(i,j,k) ]
         numq_kernel(:,:,:,:) = 0.0_wp
         do l = 1, 3
            rA_r = rA; rA_r(l) = rA_r(l) + step
            rA_l = rA; rA_l(l) = rA_l(l) - step

            call kernel%kernel_d3_pair(rA_r, rB, brad(jat), brad(jc), d3kernel_r)
            call kernel%kernel_d3_pair(rA_l, rB, brad(jat), brad(jc), d3kernel_l)

            do i = 1, 3
               do j = 1, 3
                  do k = 1, 3
                     numq_kernel(i,j,k,l) = 0.5_wp * (d3kernel_r(i,j,k) - d3kernel_l(i,j,k)) / step
                  end do
               end do
            end do
         end do

         maxdiff = maxval(abs(anaq_kernel - numq_kernel))
         if (maxdiff > thr2) then
            loc = maxloc(abs(anaq_kernel - numq_kernel))
            call test_failed(error, "Coulomb d4 (pair) does not match FD(d3)")
            print '(a,2i6,a,es20.13)', "Mismatch at pair (A,B)=(", jat, jc, "), max|diff|=", maxdiff
            print '(a,4i2)', "Worst entry (i,j,k,l)=", loc(1), loc(2), loc(3), loc(4)
            print '(a,es20.13)', "ana d4(i,j,k,l)=", anaq_kernel(loc(1),loc(2),loc(3),loc(4))
            print '(a,es20.13)', "num d4(i,j,k,l)=", numq_kernel(loc(1),loc(2),loc(3),loc(4))
            print '(a,es20.13)', "diff            =", anaq_kernel(loc(1),loc(2),loc(3),loc(4)) - numq_kernel(loc(1),loc(2),loc(3),loc(4))
            return
         end if

      end do
   end do
end subroutine test_numq

!> Test kernel_d4_pair_dborn against numerical derivative of kernel_d4_pair wrt Born radii
subroutine test_numq_dborn(error, mol, kernel_id, keps)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol
   !> Kernel identifier (still, p16, or coulomb)
   integer, intent(in) :: kernel_id
   !> Dielectric screening factor
   real(wp), intent(in) :: keps

   !> Kernel instance
   class(kernel_type), allocatable :: kernel
   !> Born radii integrator
   type(born_integrator) :: gbobc

   !> Finite difference step size
   real(wp), parameter :: step = 1.0e-6_wp

   !> Loop indices for atom pairs
   integer :: jat, jc
   !> Van der Waals radii for all atoms
   real(wp), allocatable :: rvdw(:)
   !> Born radii for all atoms
   real(wp), allocatable :: brad(:)

   !> Atomic coordinates
   real(wp) :: rA(3), rB(3)
   !> Original Born radii before perturbation
   real(wp) :: bornA0, bornB0

   !> Fourth derivative with positive/negative Born radius displacement
   real(wp) :: d4_r(3,3,3,3), d4_l(3,3,3,3)
   !> Numerical Born radius derivatives of fourth derivative
   real(wp) :: num_bA(3,3,3,3), num_bB(3,3,3,3)
   !> Analytical Born radius derivatives from kernel_d4_pair_dborn
   real(wp) :: ana_bA(3,3,3,3), ana_bB(3,3,3,3)

   !> Maximum differences for Born radius derivatives
   real(wp) :: maxdiffA, maxdiffB
   !> Locations of maximum differences
   integer :: locA(4), locB(4)

   kernel = new_kernel(kernel_id, keps)

   allocate(rvdw(mol%nat), brad(mol%nat))
   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   call gbobc%get_rad(mol, brad)

   do jat = 1, mol%nat
      do jc = 1, mol%nat
         if (jat == jc) cycle

         rA = mol%xyz(:, jat)
         rB = mol%xyz(:, jc)

         bornA0 = brad(jat)
         bornB0 = brad(jc)

         ! ---- analytic d/d(bornA), d/d(bornB)
         call kernel%kernel_d4_pair_dborn(rA, rB, bornA0, bornB0, ana_bA, ana_bB)

         ! ---- numerical d/d(bornA) of d4_pair
         call kernel%kernel_d4_pair(rA, rB, bornA0 + step, bornB0, d4_r)
         call kernel%kernel_d4_pair(rA, rB, bornA0 - step, bornB0, d4_l)
         num_bA(:,:,:,:) = 0.5_wp * (d4_r(:,:,:,:) - d4_l(:,:,:,:)) / step

         ! ---- numerical d/d(bornB) of d4_pair
         call kernel%kernel_d4_pair(rA, rB, bornA0, bornB0 + step, d4_r)
         call kernel%kernel_d4_pair(rA, rB, bornA0, bornB0 - step, d4_l)
         num_bB(:,:,:,:) = 0.5_wp * (d4_r(:,:,:,:) - d4_l(:,:,:,:)) / step

         maxdiffA = maxval(abs(ana_bA - num_bA))
         maxdiffB = maxval(abs(ana_bB - num_bB))

         if (maxdiffA > thr2 .or. maxdiffB > thr2) then
            call test_failed(error, "kernel_d4_pair_dborn does not match FD(d4_pair)")

            if (maxdiffA >= maxdiffB) then
               locA = maxloc(abs(ana_bA - num_bA))
               print '(a,2i6,a,es20.13)', "Mismatch d(d4)/d(bornA) at pair (A,B)=(", jat, jc, "), max|diff|=", maxdiffA
               print '(a,4i2)', "Worst entry (i,j,k,l)=", locA(1), locA(2), locA(3), locA(4)
               print '(a,es20.13)', "ana_bA(i,j,k,l)=", ana_bA(locA(1),locA(2),locA(3),locA(4))
               print '(a,es20.13)', "num_bA(i,j,k,l)=", num_bA(locA(1),locA(2),locA(3),locA(4))
               print '(a,es20.13)', "diff          =", ana_bA(locA(1),locA(2),locA(3),locA(4)) - num_bA(locA(1),locA(2),locA(3),locA(4))
            else
               locB = maxloc(abs(ana_bB - num_bB))
               print '(a,2i6,a,es20.13)', "Mismatch d(d4)/d(bornB) at pair (A,B)=(", jat, jc, "), max|diff|=", maxdiffB
               print '(a,4i2)', "Worst entry (i,j,k,l)=", locB(1), locB(2), locB(3), locB(4)
               print '(a,es20.13)', "ana_bB(i,j,k,l)=", ana_bB(locB(1),locB(2),locB(3),locB(4))
               print '(a,es20.13)', "num_bB(i,j,k,l)=", num_bB(locB(1),locB(2),locB(3),locB(4))
               print '(a,es20.13)', "diff          =", ana_bB(locB(1),locB(2),locB(3),locB(4)) - num_bB(locB(1),locB(2),locB(3),locB(4))
            end if

            return
         end if

      end do
   end do
end subroutine test_numq_dborn



!> Test the kernel fifth derivative against numerical derivative
subroutine test_num5(error, mol, kernel_id, keps)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol
   !> Kernel identifier (still, p16, or coulomb)
   integer, intent(in) :: kernel_id
   !> Dielectric screening factor
   real(wp), intent(in) :: keps

   !> Kernel instance
   class(kernel_type), allocatable :: kernel

   !> Born radii integrator
   type(born_integrator) :: gbobc

   !> Finite difference step size
   real(wp), parameter :: step = 1.0e-6_wp

   !> Loop indices for atom pairs and Cartesian directions
   integer :: jat, jc, i, j, k, l, m
   !> Van der Waals radii for all atoms
   real(wp), allocatable :: rvdw(:)
   !> Born radii for all atoms
   real(wp), allocatable :: brad(:)
   !> Atomic coordinates and perturbed coordinates
   real(wp) :: rA(3), rB(3), rA_r(3), rA_l(3)
   !> Interatomic distance
   real(wp) :: r
   !> Fourth derivative with positive/negative displacement
   real(wp) :: d4kernel_r(3,3,3,3), d4kernel_l(3,3,3,3)
   !> Analytical and numerical fifth derivative tensors
   real(wp) :: ana5_kernel(3,3,3,3,3), num5_kernel(3,3,3,3,3)
   !> Maximum difference between analytical and numerical fifth derivatives
   real(wp) :: maxdiff
   !> Location of maximum difference
   integer :: loc(5)

   kernel = new_kernel(kernel_id, keps)

   allocate(rvdw(mol%nat), brad(mol%nat))
   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   call gbobc%get_rad(mol, brad)

   do jat = 1, mol%nat
      do jc = 1, mol%nat
         if (jat == jc) cycle

         rA = mol%xyz(:, jat)
         rB = mol%xyz(:, jc)
         r  = norm2(rA - rB)

         ! analytic 5th derivative
         call kernel%kernel_d5_pair(rA, rB, brad(jat), brad(jc), ana5_kernel)

         ! numerical 5th derivative from analytic 4th derivative:
         ! n4(i,j,k,l) = d/dR_A,l [ d3(i,j,k) ]
         num5_kernel(:,:,:,:,:) = 0.0_wp
         do l = 1, 3
            rA_r = rA; rA_r(l) = rA_r(l) + step
            rA_l = rA; rA_l(l) = rA_l(l) - step

            call kernel%kernel_d4_pair(rA_r, rB, brad(jat), brad(jc), d4kernel_r)
            call kernel%kernel_d4_pair(rA_l, rB, brad(jat), brad(jc), d4kernel_l)

            do i = 1, 3
               do j = 1, 3
                  do k = 1, 3
                     do m = 1,3
                        num5_kernel(i,j,k,m,l) = 0.5_wp * (d4kernel_r(i,j,k,m) - d4kernel_l(i,j,k,m)) / step
                     end do 
                  end do
               end do
            end do
         end do

         maxdiff = maxval(abs(ana5_kernel - num5_kernel))
         if (maxdiff > thr2) then
            loc = maxloc(abs(ana5_kernel - num5_kernel))
            call test_failed(error, "Coulomb d5 (pair) does not match FD(d3)")
            print '(a,2i6,a,es20.13)', "Mismatch at pair (A,B)=(", jat, jc, "), max|diff|=", maxdiff
            print '(a,4i2)', "Worst entry (i,j,k,l)=", loc(1), loc(2), loc(3), loc(4), loc(5)
            print '(a,es20.13)', "ana d5(i,j,k,l)=", ana5_kernel(loc(1),loc(2),loc(3),loc(4),loc(5))
            print '(a,es20.13)', "num d5(i,j,k,l)=", num5_kernel(loc(1),loc(2),loc(3),loc(4),loc(5))
            print '(a,es20.13)', "diff            =", ana5_kernel(loc(1),loc(2),loc(3),loc(4),loc(5)) - num5_kernel(loc(1),loc(2),loc(3),loc(4),loc(5))
            return
         end if

      end do
   end do
end subroutine test_num5



!> Test setting up the kernel interacion matrices based on kernel derivatives
subroutine test_amat(error, mol, keps, input, &
   & amat_sd_ref, amat_dd_ref, amat_sq_ref, amat_dq_ref, amat_qq_ref)
   !> Error handler for test failures
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol
   !> Dielectric screening factor (1 - 1/epsilon)
   real(wp), intent(in) :: keps
   !> Reference charge-dipole interaction matrix (3,3,3)
   real(wp), intent(in) :: amat_sd_ref(:,:,:)
   !> Reference dipole-dipole interaction matrix (3,3,3,3)
   real(wp), intent(in) :: amat_dd_ref(:,:,:,:)
   !> Reference charge-quadrupole interaction matrix (6,3,3)
   real(wp), intent(in) :: amat_sq_ref(:,:,:)
   !> Reference dipole-quadrupole interaction matrix (3,3,6,3)
   real(wp), intent(in) :: amat_dq_ref(:,:,:,:)
   !> Reference quadrupole-quadrupole interaction matrix (6,3,6,3)
   real(wp), intent(in) :: amat_qq_ref(:,:,:,:)

   !> Factory to create new electrostatic objects
   type(alpb_input), intent(in) :: input

   !> Container cache for solvation data
   type(container_cache) :: cache

   !> Kernel type for solvation interactions
   class(kernel_type), allocatable :: kernel
   !> ALPB solvation model object
   type(alpb_solvation) :: solv
   !> Temporary copy of ALPB input parameters
   type(alpb_input), allocatable :: scratch_input
   !> Pointer to cached ALPB data structure containing multipole matrices
   type(alpb_cache), pointer :: ptr

   !> Born radii array (unused but allocated)
   !> Born radii derivative array (unused but allocated)
   real(wp), allocatable :: rad(:), draddr(:,:,:)

   integer :: i, j, k, l, ii, jj

   scratch_input = input
   call get_alpb_param(scratch_input, mol, 'gfn2', error)
   solv = alpb_solvation(mol, scratch_input, 'gfn2')
   
   allocate(rad(mol%nat), source=0.0_wp)
   allocate(draddr(3, mol%nat, mol%nat), source=0.0_wp)

   call taint(cache, ptr)
   call solv%update(mol, cache)
   call view(cache, ptr)


   if (abs(maxval(ptr%amat_sd - amat_sd_ref)) > thr2) then
      call test_failed(error, "Monopole-dipole interaction matrices do no match!")
      print '(a,es20.13)', "Max difference amat_sd: ", maxval(ptr%amat_sd - amat_sd_ref)
   end if
   if (abs(maxval(ptr%amat_dd - amat_dd_ref)) > thr2) then
      call test_failed(error, "Dipole-dipole interaction matrices do no match!")
      print '(a,es20.13)', "Max difference amat_dd: ", maxval(ptr%amat_dd - amat_dd_ref)
   end if
   if (abs(maxval(ptr%amat_sq - amat_sq_ref)) > thr2) then
      call test_failed(error, "Monopole-quadrupole interaction matrices do no match!")
      print '(a,es20.13)', "Max difference amat_sq: ", maxval(ptr%amat_sq - amat_sq_ref)
   end if
   if (abs(maxval(ptr%amat_dq - amat_dq_ref)) > thr2) then
      call test_failed(error, "Dipole-quadrupole interaction matrices do no match!")
      print '(a,es20.13)', "Max difference amat_dq: ", maxval(ptr%amat_dq - amat_dq_ref)
   end if
   if (abs(maxval(ptr%amat_qq - amat_qq_ref)) > thr2) then
      call test_failed(error, "Quadrupole-quadrupole interaction matrices do no match!")
      print *, "Max difference amat_qq: ", maxval(ptr%amat_qq - amat_qq_ref)
   end if


end subroutine test_amat


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine get_amat_ref_still(amat_sd_ref, amat_dd_ref, amat_sq_ref, &
   & amat_dq_ref, amat_qq_ref)
   real(wp), intent(out) :: amat_sd_ref(3,3,3)
   real(wp), intent(out) :: amat_dd_ref(3,3,3,3)
   real(wp), intent(out) :: amat_sq_ref(6,3,3)
   real(wp), intent(out) :: amat_dq_ref(3,3,6,3)
   real(wp), intent(out) :: amat_qq_ref(6,3,6,3)

   amat_sd_ref = reshape([ &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.49759202499840E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.49759202499840E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.49759202499840E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.58870044522578E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.49759202499840E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.58870044522578E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp &
      ], [3,3,3])



   amat_dd_ref = reshape([ &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-9.88004455944290E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-9.88004455944290E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-9.88004455944290E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-9.88004455944290E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.42990642057011E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.42990642057011E-3_wp,&
      &-9.88004455944290E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-5.12022690132844E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-9.88004455944290E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-5.12022690132844E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.42990642057011E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 5.16579034330314E-3_wp,&
      &-9.88004455944290E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-5.12022690132844E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-9.88004455944290E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-5.12022690132844E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.42990642057011E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 5.16579034330314E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp &
      ], [3,3,3,3])

    amat_sq_ref = reshape([ &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.48337937962426E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.48337937962426E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.48337937962426E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.42867241487719E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.48337937962426E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.42867241487719E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp &
      ], [6,3,3])

   amat_dq_ref = reshape([ &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-9.82382177838076E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.82382177838076E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-9.82382177838076E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.82382177838076E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.96476435567615E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.96476435567615E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.96476435567615E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.96476435567615E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.71978808702705E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.71978808702705E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.82382177838076E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 6.78161923557974E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.82382177838076E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 6.78161923557974E-4_wp,&
      & 1.96476435567615E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.35632384711595E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.96476435567615E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.35632384711595E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.71978808702705E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.73000361720903E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-9.82382177838076E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-6.78161923557974E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-9.82382177838076E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-6.78161923557974E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.96476435567615E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.35632384711595E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.96476435567615E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.35632384711595E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.71978808702705E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.73000361720903E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp &
      ], [3,3,6,3])

   amat_qq_ref = reshape([ &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.74876073331153E-3_wp, 0.00000000000000E+0_wp,-5.82920244437178E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-9.73983489342622E-5_wp,&
      &-1.74876073331153E-3_wp, 0.00000000000000E+0_wp,-5.82920244437178E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-9.73983489342622E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.33168097774871E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.33168097774871E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-5.82920244437178E-4_wp, 0.00000000000000E+0_wp,-1.74876073331153E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-9.73983489342622E-5_wp,&
      &-5.82920244437178E-4_wp, 0.00000000000000E+0_wp,-1.74876073331153E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-9.73983489342622E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-3.89593395737049E-4_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-3.89593395737049E-4_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-3.89593395737049E-4_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-3.89593395737049E-4_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-9.73983489342622E-5_wp, 0.00000000000000E+0_wp,-9.73983489342622E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 3.11932087590825E-4_wp,&
      &-9.73983489342622E-5_wp, 0.00000000000000E+0_wp,-9.73983489342622E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 3.11932087590825E-4_wp,&
      &-1.74876073331153E-3_wp, 0.00000000000000E+0_wp,-5.82920244437178E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-9.73983489342622E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-6.03605689058357E-4_wp, 0.00000000000000E+0_wp,-2.01201896352786E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.35419880424931E-4_wp,&
      & 0.00000000000000E+0_wp,-2.33168097774871E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-8.04807585411143E-4_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-5.82920244437178E-4_wp, 0.00000000000000E+0_wp,-1.74876073331153E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-9.73983489342622E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.01201896352786E-4_wp, 0.00000000000000E+0_wp,-6.03605689058357E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.35419880424931E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-3.89593395737049E-4_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.41679521699725E-4_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-3.89593395737049E-4_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.41679521699725E-4_wp, 0.00000000000000E+0_wp,&
      &-9.73983489342622E-5_wp, 0.00000000000000E+0_wp,-9.73983489342622E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 3.11932087590825E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 2.35419880424931E-4_wp, 0.00000000000000E+0_wp, 2.35419880424931E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 3.34121574004972E-5_wp,&
      &-1.74876073331153E-3_wp, 0.00000000000000E+0_wp,-5.82920244437178E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-9.73983489342622E-5_wp,&
      &-6.03605689058357E-4_wp, 0.00000000000000E+0_wp,-2.01201896352786E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.35419880424931E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.33168097774871E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-8.04807585411143E-4_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-5.82920244437178E-4_wp, 0.00000000000000E+0_wp,-1.74876073331153E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-9.73983489342622E-5_wp,&
      &-2.01201896352786E-4_wp, 0.00000000000000E+0_wp,-6.03605689058357E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.35419880424931E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-3.89593395737049E-4_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.41679521699725E-4_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-3.89593395737049E-4_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.41679521699725E-4_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-9.73983489342622E-5_wp, 0.00000000000000E+0_wp,-9.73983489342622E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 3.11932087590825E-4_wp,&
      & 2.35419880424931E-4_wp, 0.00000000000000E+0_wp, 2.35419880424931E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 3.34121574004972E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp &
      ], [6,3,6,3])


end subroutine get_amat_ref_still



subroutine get_amat_ref_p16(amat_sd_ref, amat_dd_ref, amat_sq_ref, &
   & amat_dq_ref, amat_qq_ref)
   real(wp), intent(out) :: amat_sd_ref(3,3,3)
   real(wp), intent(out) :: amat_dd_ref(3,3,3,3)
   real(wp), intent(out) :: amat_sq_ref(6,3,3)
   real(wp), intent(out) :: amat_dq_ref(3,3,6,3)
   real(wp), intent(out) :: amat_qq_ref(6,3,6,3)


   amat_sd_ref = reshape([ &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.65921845083053E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.65921845083053E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.65921845083053E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.48236785603179E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.65921845083053E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.48236785603179E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp &
      ], [3,3,3])


   amat_dd_ref = reshape([ &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.05194109063970E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.05194109063970E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.05194109063970E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.05194109063970E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.06516102610275E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.06516102610275E-3_wp,&
      &-1.05194109063970E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-4.90991018249639E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.05194109063970E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-4.90991018249639E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.06516102610275E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 5.30193154514398E-3_wp,&
      &-1.05194109063970E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-4.90991018249639E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.05194109063970E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-4.90991018249639E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.06516102610275E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 5.30193154514398E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp &
      ], [3,3,3,3])

    amat_sq_ref = reshape([ &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.15141662676476E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.15141662676476E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.15141662676476E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.40394724254679E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.15141662676476E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.40394724254679E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp &
    ], [6,3,3])

   amat_dq_ref = reshape([ &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.24664622509055E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.24664622509055E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.24664622509055E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.24664622509055E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.49329245018109E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 2.49329245018109E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.49329245018109E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 2.49329245018109E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.72933625417785E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.72933625417785E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.24664622509055E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 6.73271497060787E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.24664622509055E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 6.73271497060787E-4_wp,&
      & 2.49329245018109E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.34654299412157E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 2.49329245018109E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.34654299412157E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.72933625417785E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.91929119230186E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.24664622509055E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-6.73271497060787E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.24664622509055E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-6.73271497060787E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.49329245018109E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.34654299412157E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.49329245018109E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.34654299412157E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.72933625417785E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.91929119230186E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp &
      ], [3,3,6,3])

   amat_qq_ref = reshape([ &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.21918314068676E-3_wp, 0.00000000000000E+0_wp,-7.39727713562253E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 5.56320294583537E-5_wp,&
      &-2.21918314068676E-3_wp, 0.00000000000000E+0_wp,-7.39727713562253E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 5.56320294583537E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.95891085424901E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.95891085424901E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-7.39727713562253E-4_wp, 0.00000000000000E+0_wp,-2.21918314068676E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 5.56320294583537E-5_wp,&
      &-7.39727713562253E-4_wp, 0.00000000000000E+0_wp,-2.21918314068676E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 5.56320294583537E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 2.22528117833415E-4_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 2.22528117833415E-4_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 2.22528117833415E-4_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 2.22528117833415E-4_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 5.56320294583537E-5_wp, 0.00000000000000E+0_wp, 5.56320294583537E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 5.92849366434424E-4_wp,&
      & 5.56320294583537E-5_wp, 0.00000000000000E+0_wp, 5.56320294583537E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 5.92849366434424E-4_wp,&
      &-2.21918314068676E-3_wp, 0.00000000000000E+0_wp,-7.39727713562253E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 5.56320294583537E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-5.99252909651137E-4_wp, 0.00000000000000E+0_wp,-1.99750969883712E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.57492039628580E-4_wp,&
      & 0.00000000000000E+0_wp,-2.95891085424901E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-7.99003879534850E-4_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-7.39727713562253E-4_wp, 0.00000000000000E+0_wp,-2.21918314068676E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 5.56320294583537E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.99750969883712E-4_wp, 0.00000000000000E+0_wp,-5.99252909651137E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.57492039628580E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 2.22528117833415E-4_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.02996815851432E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 2.22528117833415E-4_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.02996815851432E-3_wp, 0.00000000000000E+0_wp,&
      & 5.56320294583537E-5_wp, 0.00000000000000E+0_wp, 5.56320294583537E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 5.92849366434424E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 2.57492039628580E-4_wp, 0.00000000000000E+0_wp, 2.57492039628580E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-6.20121728894948E-5_wp,&
      &-2.21918314068676E-3_wp, 0.00000000000000E+0_wp,-7.39727713562253E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 5.56320294583537E-5_wp,&
      &-5.99252909651137E-4_wp, 0.00000000000000E+0_wp,-1.99750969883712E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.57492039628580E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.95891085424901E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-7.99003879534850E-4_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-7.39727713562253E-4_wp, 0.00000000000000E+0_wp,-2.21918314068676E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 5.56320294583537E-5_wp,&
      &-1.99750969883712E-4_wp, 0.00000000000000E+0_wp,-5.99252909651137E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.57492039628580E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 2.22528117833415E-4_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.02996815851432E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 2.22528117833415E-4_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.02996815851432E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 5.56320294583537E-5_wp, 0.00000000000000E+0_wp, 5.56320294583537E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 5.92849366434424E-4_wp,&
      & 2.57492039628580E-4_wp, 0.00000000000000E+0_wp, 2.57492039628580E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-6.20121728894948E-5_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp &
      ], [6,3,6,3])

end subroutine get_amat_ref_p16




subroutine get_amat_ref_coulomb(amat_sd_ref, amat_dd_ref, amat_sq_ref, &
   & amat_dq_ref, amat_qq_ref)
   real(wp), intent(out) :: amat_sd_ref(3,3,3)
   real(wp), intent(out) :: amat_dd_ref(3,3,3,3)
   real(wp), intent(out) :: amat_sq_ref(6,3,3)
   real(wp), intent(out) :: amat_dq_ref(3,3,6,3)
   real(wp), intent(out) :: amat_qq_ref(6,3,6,3)

   amat_sd_ref = reshape([ &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.56485754645450E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.56485754645450E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.56485754645450E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.91214386613625E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.56485754645450E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 3.91214386613625E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp &
      ], [3,3,3])



   amat_dd_ref = reshape([ &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 6.19030735740793E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 6.19030735740793E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 6.19030735740793E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 6.19030735740793E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.23806147148159E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.23806147148159E-1_wp,&
      & 6.19030735740793E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 7.73788419675991E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 6.19030735740793E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 7.73788419675991E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.23806147148159E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.54757683935198E-2_wp,&
      & 6.19030735740793E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 7.73788419675991E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 6.19030735740793E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 7.73788419675991E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.23806147148159E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.54757683935198E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp &
      ], [3,3,3,3])

   amat_sq_ref = reshape([ &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 6.19030735740793E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 6.19030735740793E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 6.19030735740793E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 7.73788419675991E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 6.19030735740793E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 7.73788419675991E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp &
    ], [6,3,3])

   amat_dq_ref = reshape([ &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.44877914069560E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.44877914069560E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.44877914069560E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.44877914069560E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.89755828139120E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-4.89755828139120E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 4.89755828139120E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-4.89755828139120E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-4.89755828139120E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 4.89755828139120E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.44877914069560E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.53048696293475E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.44877914069560E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.53048696293475E-3_wp,&
      &-4.89755828139120E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-3.06097392586950E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-4.89755828139120E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-3.06097392586950E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 4.89755828139120E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 3.06097392586950E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.44877914069560E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.53048696293475E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.44877914069560E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.53048696293475E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.89755828139120E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 3.06097392586950E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 4.89755828139120E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 3.06097392586950E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-4.89755828139120E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.06097392586950E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp &
   ], [3,3,6,3])

   amat_qq_ref = reshape([ &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.35912713240068E-2_wp, 0.00000000000000E+0_wp, 1.45304237746689E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.39043221408942E-2_wp,&
      & 4.35912713240068E-2_wp, 0.00000000000000E+0_wp, 1.45304237746689E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.39043221408942E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 5.81216950986757E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 5.81216950986757E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.45304237746689E-2_wp, 0.00000000000000E+0_wp, 4.35912713240068E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.39043221408942E-2_wp,&
      & 1.45304237746689E-2_wp, 0.00000000000000E+0_wp, 4.35912713240068E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.39043221408942E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.35617288563577E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.35617288563577E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.35617288563577E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.35617288563577E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-3.39043221408942E-2_wp, 0.00000000000000E+0_wp,-3.39043221408942E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.20260172395699E-2_wp,&
      &-3.39043221408942E-2_wp, 0.00000000000000E+0_wp,-3.39043221408942E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.20260172395699E-2_wp,&
      & 4.35912713240068E-2_wp, 0.00000000000000E+0_wp, 1.45304237746689E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.39043221408942E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.36222722887521E-3_wp, 0.00000000000000E+0_wp, 4.54075742958404E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.05951006690294E-3_wp,&
      & 0.00000000000000E+0_wp, 5.81216950986757E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.81630297183361E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.45304237746689E-2_wp, 0.00000000000000E+0_wp, 4.35912713240068E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.39043221408942E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.54075742958404E-4_wp, 0.00000000000000E+0_wp, 1.36222722887521E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.05951006690294E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.35617288563577E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-4.23804026761177E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.35617288563577E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-4.23804026761177E-3_wp, 0.00000000000000E+0_wp,&
      &-3.39043221408942E-2_wp, 0.00000000000000E+0_wp,-3.39043221408942E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.20260172395699E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.05951006690294E-3_wp, 0.00000000000000E+0_wp,-1.05951006690294E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.87581303873656E-3_wp,&
      & 4.35912713240068E-2_wp, 0.00000000000000E+0_wp, 1.45304237746689E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.39043221408942E-2_wp,&
      & 1.36222722887521E-3_wp, 0.00000000000000E+0_wp, 4.54075742958404E-4_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.05951006690294E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 5.81216950986757E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.81630297183361E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.45304237746689E-2_wp, 0.00000000000000E+0_wp, 4.35912713240068E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-3.39043221408942E-2_wp,&
      & 4.54075742958404E-4_wp, 0.00000000000000E+0_wp, 1.36222722887521E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.05951006690294E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.35617288563577E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-4.23804026761177E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.35617288563577E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-4.23804026761177E-3_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-3.39043221408942E-2_wp, 0.00000000000000E+0_wp,-3.39043221408942E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.20260172395699E-2_wp,&
      &-1.05951006690294E-3_wp, 0.00000000000000E+0_wp,-1.05951006690294E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.87581303873656E-3_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp &
      ], [6,3,6,3])

end subroutine get_amat_ref_coulomb





subroutine taint(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(alpb_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(alpb_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

subroutine view(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(alpb_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(alpb_cache)
      ptr => target
   end select
end subroutine view




end module test_solvation_kernel
