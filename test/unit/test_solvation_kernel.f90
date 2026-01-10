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
   use tblite_solvation_kernel, only : kernel_type, new_kernel, kernel_enum, compute_kernel_dkdr
   implicit none
   private

   public :: collect_solvation_kernel

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains

!> Collect all exported unit tests
subroutine collect_solvation_kernel(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("kernel-gradient-still", test_kernel_gradient_still), &
      new_unittest("kernel-gradient-p16", test_kernel_gradient_p16) &
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
   real(wp), parameter :: qat(*) = [&
      & -2.11018727757438E-1_wp, -6.04389222813257E-2_wp, -1.90601159250311E-1_wp, &
      &  1.49237694872530E-1_wp,  1.35835820853652E-1_wp,  1.27732431639016E-1_wp, &
      &  1.78559147201780E-1_wp,  1.42324484825195E-1_wp,  1.92106458233743E-1_wp, &
      &  1.45841758574287E-1_wp,  1.56456166394024E-1_wp,  1.59746890863949E-1_wp, &
      & -2.70765876809499E-1_wp, -3.27435355522312E-1_wp, -4.70046325670683E-2_wp, &
      &  1.10838969762146E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_kernel_numg(error, mol, kernel_enum%still, keps, qat)

end subroutine test_kernel_gradient_still


!> Test P16 kernel gradient against numerical derivative
subroutine test_kernel_gradient_p16(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: keps = 0.5_wp
   real(wp), parameter :: qat(*) = [&
      & -2.11018727757438E-1_wp, -6.04389222813257E-2_wp, -1.90601159250311E-1_wp, &
      &  1.49237694872530E-1_wp,  1.35835820853652E-1_wp,  1.27732431639016E-1_wp, &
      &  1.78559147201780E-1_wp,  1.42324484825195E-1_wp,  1.92106458233743E-1_wp, &
      &  1.45841758574287E-1_wp,  1.56456166394024E-1_wp,  1.59746890863949E-1_wp, &
      & -2.70765876809499E-1_wp, -3.27435355522312E-1_wp, -4.70046325670683E-2_wp, &
      &  1.10838969762146E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_kernel_numg(error, mol, kernel_enum%p16, keps, qat)

end subroutine test_kernel_gradient_p16


!> Test kernel gradient against numerical derivative by finite difference
subroutine test_kernel_numg(error, mol, kernel_id, keps, qat)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol
   !> Kernel type identifier
   integer, intent(in) :: kernel_id
   !> Dielectric screening parameter
   real(wp), intent(in) :: keps
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), allocatable :: rvdw(:), rad(:), draddr(:, :, :)
   real(wp), allocatable :: amat_r(:, :), amat_l(:, :), amat_0(:, :)
   real(wp), allocatable :: kernel_grad_spatial(:, :, :), kernel_grad_born(:, :)
   real(wp), allocatable :: numg_kernel(:, :, :, :), ana_kernel(:, :, :, :)
   real(wp), allocatable :: gradient_numerical(:, :), gradient_analytical(:, :)
   real(wp), parameter :: step = 1.0e-6_wp
   integer :: iat, jat, ic, jc

   ! Initialize Born integrator and kernel
   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   kernel = new_kernel(kernel_id, keps)

   ! Allocate arrays
   allocate(rad(mol%nat), draddr(3, mol%nat, mol%nat))
   allocate(amat_r(mol%nat, mol%nat), amat_l(mol%nat, mol%nat), amat_0(mol%nat, mol%nat))
   allocate(kernel_grad_spatial(3, mol%nat, mol%nat))
   allocate(kernel_grad_born(mol%nat, mol%nat))
   allocate(numg_kernel(3, mol%nat, mol%nat, mol%nat))
   allocate(ana_kernel(3, mol%nat, mol%nat, mol%nat))

   ! Compute numerical gradient of kernel by finite difference
   ! When we displace atom k and recompute Born radii and kernel,
   ! we get dK_mn/dr_k - how each kernel element K_mn changes when atom k moves
   numg_kernel(:, :, :, :) = 0.0_wp
   
   do iat = 1, mol%nat  ! iat = k, the atom we're displacing
      do ic = 1, 3  ! ic = alpha, the direction
         ! Right displacement
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call gbobc%get_rad(mol, rad)
         amat_r(:, :) = 0.0_wp
         call kernel%add_kernel_mat(mol%nat, mol%xyz, rad, amat_r)

         ! Left displacement
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call gbobc%get_rad(mol, rad)
         amat_l(:, :) = 0.0_wp
         call kernel%add_kernel_mat(mol%nat, mol%xyz, rad, amat_l)

         ! Restore coordinate
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step

         ! Numerical derivative: dK_mn/dr_k,alpha
         ! amat_r(m, n) and amat_l(m, n) are kernel matrices
         ! numg_kernel(alpha, k, m, n) = dK_mn/dr_k,alpha
         do jat = 1, mol%nat  ! jat = m (row index)
            do jc = 1, mol%nat  ! jc = n (column index)
               numg_kernel(ic, iat, jat, jc) = 0.5_wp * (amat_r(jat, jc) - amat_l(jat, jc)) / step
            end do
         end do
      end do
   end do

   ! Get Born radii at original geometry and compute analytical kernel gradient
   call gbobc%get_rad(mol, rad, draddr)
   call kernel%add_kernel_mat(mol%nat, mol%xyz, rad, amat_r)
   
   ! Compute analytical kernel gradient using compute_kernel_deriv
   call compute_kernel_dkdr(kernel_id, keps, mol%nat, mol%xyz, rad, draddr, ana_kernel)
   

   ! Compare analytical and numerical energy gradients
   if (any(abs(ana_kernel - numg_kernel) > thr2)) then
      call test_failed(error, "Kernel gradient does not match finite difference solution")
      print '(a)', "Analytical energy gradient:"
      print '(3es20.13)', ana_kernel
      print '(a)', "Numerical energy gradient:"
      print '(3es20.13)', numg_kernel
      print '(a)', "Difference:"
      print '(3es20.13)', ana_kernel - numg_kernel
   end if

end subroutine test_kernel_numg


end module test_solvation_kernel
