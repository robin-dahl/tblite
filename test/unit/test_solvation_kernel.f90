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
   use tblite_solvation_kernel, only : kernel_type, new_kernel, kernel_enum, compute_kernel_dKdr_ij, compute_kernel_d2Kdr2_ij, compute_kernel_d3Kdr3_ij, &
      & compute_kernel_d4Kdr4_ij
   use tblite_solvation_alpb, only : alpb_solvation, alpb_input, get_multipole_matrix, alpb_cache
   use tblite_solvation_data_alpb, only : get_alpb_param
   use tblite_solvation_data, only : solvent_data, get_vdw_rad_d3, get_solvent_data


   use mctc_io_structure, only: new_structure

   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_multipole, only : damped_multipole, new_damped_multipole, get_multipole_matrix_0d
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
      new_unittest("amat-coulomb", test_amat_coulomb) &
      !   new_unittest("amat-higher-order-coulomb", test_amat_higher_order_coulomb) &
      ! new_unittest("kernel-gradient-still", test_kernel_gradient_still), &
      ! new_unittest("kernel-hessian-still", test_kernel_hessian_still), &
      ! new_unittest("kernel-third-still", test_kernel_third_still), &
      !   new_unittest("kernel-fourth-still", test_kernel_fourth_still), &
      ! new_unittest("kernel-gradient-p16", test_kernel_gradient_p16), &
      ! new_unittest("kernel-hessian-p16", test_kernel_hessian_p16), &
      ! new_unittest("kernel-third-p16", test_kernel_third_p16), &
      !   new_unittest("kernel-fourth-p16", test_kernel_fourth_p16) &
      ! new_unittest("kernel-gradient-coulomb", test_kernel_gradient_coulomb), &
      ! new_unittest("kernel-hessian-coulomb", test_kernel_hessian_coulomb), &
      ! new_unittest("kernel-third-coulomb", test_kernel_third_coulomb) &
      !   new_unittest("kernel-fourth-coulomb", test_kernel_fourth_coulomb) &
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

!> Test Still kernel Hessian against numerical derivative
subroutine test_kernel_hessian_still(error)
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
   call test_kernel_numh(error, mol, kernel_enum%still, keps, qat)

end subroutine test_kernel_hessian_still


!> Test Still kernel third derivative against numerical derivative
subroutine test_kernel_third_still(error)
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
   call test_kernel_numt(error, mol, kernel_enum%still, keps, qat)

end subroutine test_kernel_third_still



!> Test Still kernel fourth derivative against numerical derivative
!  Turned off due to long runtimes
subroutine test_kernel_fourth_still(error)
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
   call test_kernel_numq(error, mol, kernel_enum%still, keps, qat)

end subroutine test_kernel_fourth_still


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


!> Test P16 kernel Hessian against numerical derivative
subroutine test_kernel_hessian_p16(error)
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
   call test_kernel_numh(error, mol, kernel_enum%p16, keps, qat)

end subroutine test_kernel_hessian_p16


!> Test P16 kernel third derivative against numerical derivative
subroutine test_kernel_third_p16(error)
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
   call test_kernel_numt(error, mol, kernel_enum%p16, keps, qat)

end subroutine test_kernel_third_p16

!> Test P16 kernel fourth derivative against numerical derivative
!  Turned off due to long runtimes
subroutine test_kernel_fourth_p16(error)
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
   call test_kernel_numq(error, mol, kernel_enum%p16, keps, qat)

end subroutine test_kernel_fourth_p16


!> Test Coulomb kernel gradient against numerical derivative
subroutine test_kernel_gradient_coulomb(error)
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

end subroutine test_kernel_gradient_coulomb

!> Test Coulomb kernel Hessian against numerical derivative
subroutine test_kernel_hessian_coulomb(error)
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
   call test_kernel_numh(error, mol, kernel_enum%coulomb, keps, qat)

end subroutine test_kernel_hessian_coulomb


!> Test Coulomb kernel third derivative against numerical derivative
subroutine test_kernel_third_coulomb(error)
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
   call test_kernel_numt(error, mol, kernel_enum%coulomb, keps, qat)

end subroutine test_kernel_third_coulomb

!> Test Coulomb kernel fourth derivative against numerical derivative
!  Turned off due to long runtimes
subroutine test_kernel_fourth_coulomb(error)
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
   call test_kernel_numq(error, mol, kernel_enum%coulomb, keps, qat)

end subroutine test_kernel_fourth_coulomb



!> Test if amat construction based on derivative routines works
subroutine test_amat_coulomb(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(alpb_input) :: input
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

   real(wp), allocatable :: rad(:), ds(:)

   solvent = get_solvent_data("water")
   input = alpb_input(solvent%eps, solvent=solvent%solvent, &
         & kernel=3, alpb=.true.)
   call get_structure(mol, "MB16-43", "01")
   call test_amat(error, mol, kernel_enum%coulomb, keps, qat, make_multipole2, input) 

end subroutine test_amat_coulomb

subroutine test_amat_higher_order_coulomb(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   character(len=*), parameter :: sym(2) =  [ 'H', 'H' ]
   integer, parameter :: num(2) =  [ 1, 1 ]
   real(wp) :: xyz(3, 2)

   xyz(:,1) = [0.0_wp, 0.0_wp, 0.0_wp]
   xyz(:,2) = [0.0_wp, 0.0_wp, 0.74_wp]

   call new_structure(mol, num, sym, xyz)

   call test_amat_ho(error, mol)

end subroutine test_amat_higher_order_coulomb



!> Test kernel gradient against numerical one.
subroutine test_kernel_numg(error, mol, kernel_id, keps, qat)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type), intent(inout) :: mol
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   real(wp), intent(in) :: qat(:)

   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), allocatable :: rvdw(:), rad(:), draddr(:, :, :)
   real(wp), allocatable :: amat_r(:, :), amat_l(:, :)
   real(wp), allocatable :: numg_kernel(:, :, :, :)
   real(wp), allocatable :: ana_ij(:, :)
   real(wp), parameter :: step = 1.0e-6_wp
   integer :: iat, jat, ic, jc
   real(wp) :: maxdiff

   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   kernel = new_kernel(kernel_id, keps)

   allocate(rad(mol%nat), draddr(3, mol%nat, mol%nat))
   allocate(amat_r(mol%nat, mol%nat), amat_l(mol%nat, mol%nat))
   allocate(numg_kernel(3, mol%nat, mol%nat, mol%nat))
   allocate(ana_ij(3, mol%nat))

   ! --- Numerical derivative: numg_kernel(alpha, k, m, n) = dK_mn / dr_k,alpha
   numg_kernel(:, :, :, :) = 0.0_wp

   do iat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call gbobc%get_rad(mol, rad)
         amat_r(:, :) = 0.0_wp
         call kernel%add_kernel_mat(mol%nat, mol%xyz, rad, amat_r)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call gbobc%get_rad(mol, rad)
         amat_l(:, :) = 0.0_wp
         call kernel%add_kernel_mat(mol%nat, mol%xyz, rad, amat_l)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step

         do jat = 1, mol%nat
            do jc = 1, mol%nat
               numg_kernel(ic, iat, jat, jc) = 0.5_wp * (amat_r(jat, jc) - amat_l(jat, jc)) / step
            end do
         end do
      end do
   end do

   ! --- Analytical: get brad + brdr once at reference geometry
   call gbobc%get_rad(mol, rad, draddr)

   ! Compare: for each (m,n), compute ana_ij(:,k)=dK_mn/dr_k and compare to numg_kernel(:,:,m,n)
   do jat = 1, mol%nat
      do jc = 1, mol%nat
         call compute_kernel_dkdr_ij(kernel_id, keps, mol%nat, mol%xyz, rad, draddr, jat, jc, ana_ij)

         maxdiff = maxval(abs(ana_ij(:, :) - numg_kernel(:, :, jat, jc)))
         if (maxdiff > thr2) then
            call test_failed(error, "Kernel gradient does not match finite difference solution")
            print '(a,2i6, a, es20.13)', "Mismatch at (m,n)=(", jat, jc, "), max|diff|=", maxdiff
            print '(a)', "Analytical dK_mn/dr_k (3,nat):"
            print '(3es20.13)', ana_ij
            print '(a)', "Numerical dK_mn/dr_k (3,nat):"
            print '(3es20.13)', numg_kernel(:, :, jat, jc)
            print '(a)', "Difference (ana - num):"
            print '(3es20.13)', ana_ij - numg_kernel(:, :, jat, jc)
            return
         end if
      end do
   end do
end subroutine test_kernel_numg

!> Test kernel 2nd derivative against numerical finite difference of the 
!> analytical 1st derivative.
subroutine test_kernel_numh(error, mol, kernel_id, keps, qat)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type), intent(inout) :: mol
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   real(wp), intent(in) :: qat(:)

   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel

   real(wp), allocatable :: rvdw(:), rad(:)
   real(wp), allocatable :: draddr(:, :, :)
   real(wp), allocatable :: draddr2(:, :, :, :, :)

   real(wp), allocatable :: dkdr_p_ij(:, :)     ! (3,nat)
   real(wp), allocatable :: dkdr_m_ij(:, :)     ! (3,nat)
   real(wp), allocatable :: num2(:, :, :, :)    ! (3,nat,3,nat) for one (m,n)
   real(wp), allocatable :: ana2(:, :, :, :)    ! (3,nat,3,nat) for one (m,n)

   real(wp), parameter :: step = 1.0e-6_wp
   integer :: nat, m, n, l, beta, k, alpha
   real(wp) :: diff, maxdiff
   integer :: imax_alpha, imax_k, imax_beta, imax_l, imax_m, imax_n

   nat = mol%nat

   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   kernel = new_kernel(kernel_id, keps)

   allocate(rad(nat))
   allocate(draddr(3, nat, nat))
   allocate(draddr2(3, nat, 3, nat, nat))

   allocate(dkdr_p_ij(3, nat), dkdr_m_ij(3, nat))
   allocate(num2(3, nat, 3, nat))
   allocate(ana2(3, nat, 3, nat))

   ! Analytical Born radii derivatives at reference geometry (need brdr2!)
   call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2)

   maxdiff = 0.0_wp
   imax_alpha=1; imax_k=1; imax_beta=1; imax_l=1; imax_m=1; imax_n=1

   do m = 1, nat
      do n = 1, nat

         ! ---- Analytical slice for this (m,n)
         call compute_kernel_d2kdr2_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, m, n, ana2)

         ! ---- Numerical slice for this (m,n) by FD of dkdr_ij
         num2(:, :, :, :) = 0.0_wp

         do l = 1, nat
            do beta = 1, 3

               ! +step
               mol%xyz(beta, l) = mol%xyz(beta, l) + step
               call gbobc%get_rad(mol, rad, draddr)
               call compute_kernel_dkdr_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, m, n, dkdr_p_ij)

               ! -step
               mol%xyz(beta, l) = mol%xyz(beta, l) - 2.0_wp*step
               call gbobc%get_rad(mol, rad, draddr)
               call compute_kernel_dkdr_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, m, n, dkdr_m_ij)

               ! restore
               mol%xyz(beta, l) = mol%xyz(beta, l) + step

               do k = 1, nat
                  do alpha = 1, 3
                     num2(alpha, k, beta, l) = 0.5_wp * (dkdr_p_ij(alpha, k) - dkdr_m_ij(alpha, k)) / step
                  end do
               end do

            end do
         end do

         ! ---- Compare this slice
         do k = 1, nat
            do alpha = 1, 3
               do l = 1, nat
                  do beta = 1, 3
                     diff = abs(ana2(alpha, k, beta, l) - num2(alpha, k, beta, l))
                     if (diff > maxdiff) then
                        maxdiff = diff
                        imax_alpha = alpha
                        imax_k     = k
                        imax_beta  = beta
                        imax_l     = l
                        imax_m     = m
                        imax_n     = n
                     end if
                  end do
               end do
            end do
         end do

      end do
   end do

   if (maxdiff > thr2) then
      call test_failed(error, "Kernel second derivative does not match finite difference solution")
      print '(a,es20.13)', "Max |d2K/dr2| difference: ", maxdiff
      print '(a,6(i0,1x))', "At indices (alpha,k,beta,l,m,n) = ", imax_alpha, imax_k, imax_beta, imax_l, imax_m, imax_n
      ! Recompute the offending slice and report the single entry
      call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2)
      call compute_kernel_d2kdr2_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, imax_m, imax_n, ana2)

      mol%xyz(imax_beta, imax_l) = mol%xyz(imax_beta, imax_l) + step
      call gbobc%get_rad(mol, rad, draddr)
      call compute_kernel_dkdr_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, imax_m, imax_n, dkdr_p_ij)

      mol%xyz(imax_beta, imax_l) = mol%xyz(imax_beta, imax_l) - 2.0_wp*step
      call gbobc%get_rad(mol, rad, draddr)
      call compute_kernel_dkdr_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, imax_m, imax_n, dkdr_m_ij)

      mol%xyz(imax_beta, imax_l) = mol%xyz(imax_beta, imax_l) + step

      num2(imax_alpha, imax_k, imax_beta, imax_l) = 0.5_wp * (dkdr_p_ij(imax_alpha, imax_k) - dkdr_m_ij(imax_alpha, imax_k)) / step

      print '(a,es20.13)', "Analytical value: ", ana2(imax_alpha, imax_k, imax_beta, imax_l)
      print '(a,es20.13)', "Numerical  value: ", num2(imax_alpha, imax_k, imax_beta, imax_l)
      print '(a,es20.13)', "Difference       : ", ana2(imax_alpha, imax_k, imax_beta, imax_l) - num2(imax_alpha, imax_k, imax_beta, imax_l)
   end if

end subroutine test_kernel_numh



!> Test kernel 3rd derivative against numerical finite difference of the 
!> analytical 2nd derivative.
subroutine test_kernel_numt(error, mol, kernel_id, keps, qat)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type), intent(inout) :: mol
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   real(wp), intent(in) :: qat(:)

   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel

   real(wp), allocatable :: rvdw(:), rad(:)
   real(wp), allocatable :: draddr(:, :, :)
   real(wp), allocatable :: draddr2(:, :, :, :, :)
   real(wp), allocatable :: draddr3(:, :, :, :, :, :, :)

   real(wp), allocatable :: d2p(:, :, :, :)      ! (3,nat,3,nat) for one (i,j)
   real(wp), allocatable :: d2m(:, :, :, :)      ! (3,nat,3,nat) for one (i,j)
   real(wp), allocatable :: num3(:, :, :, :, :, :) ! (3,nat,3,nat,3,nat) for one (i,j)
   real(wp), allocatable :: ana3(:, :, :, :, :, :) ! (3,nat,3,nat,3,nat) for one (i,j)

   real(wp), parameter :: step = 1.0e-6_wp
   integer :: nat, i, j, k, l, m, alpha, beta, gamma
   real(wp) :: diff, maxdiff
   integer :: ia, ik, ib, il, ig, im, ii, ij  ! index record for reporting

   nat = mol%nat
   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   kernel = new_kernel(kernel_id, keps)

   allocate(rad(nat))
   allocate(draddr(3,nat,nat))
   allocate(draddr2(3,nat,3,nat,nat))
   allocate(draddr3(3,nat,3,nat,3,nat,nat))

   allocate(d2p(3,nat,3,nat), d2m(3,nat,3,nat))
   allocate(num3(3,nat,3,nat,3,nat))
   allocate(ana3(3,nat,3,nat,3,nat))

   ! reference geometry: need up to brdr3 for analytical
   call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2, dradd3r=draddr3)

   maxdiff = 0.0_wp
   ia=1;ik=1;ib=1;il=1;ig=1;im=1;ii=1;ij=1

   do i = 1, nat
      do j = 1, nat

         ! ---- analytical 3rd-derivative slab for this (i,j)
         call compute_kernel_d3kdr3_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, draddr3, i, j, ana3)

         ! ---- numerical slab via FD of the ij Hessian
         num3(:, :, :, :, :, :) = 0.0_wp

         do m = 1, nat
            do gamma = 1, 3

               mol%xyz(gamma, m) = mol%xyz(gamma, m) + step
               call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2)
               call compute_kernel_d2kdr2_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, i, j, d2p)

               mol%xyz(gamma, m) = mol%xyz(gamma, m) - 2.0_wp*step
               call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2)
               call compute_kernel_d2kdr2_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, i, j, d2m)

               mol%xyz(gamma, m) = mol%xyz(gamma, m) + step

               do k = 1, nat
                  do alpha = 1, 3
                     do l = 1, nat
                        do beta = 1, 3
                           num3(alpha,k,beta,l,gamma,m) = 0.5_wp * (d2p(alpha,k,beta,l) - d2m(alpha,k,beta,l)) / step
                        end do
                     end do
                  end do
               end do

            end do
         end do

         ! ---- compare this slab; track worst entry globally
         do k = 1, nat
            do alpha = 1, 3
               do l = 1, nat
                  do beta = 1, 3
                     do m = 1, nat
                        do gamma = 1, 3
                           diff = abs(ana3(alpha,k,beta,l,gamma,m) - num3(alpha,k,beta,l,gamma,m))
                           if (diff > maxdiff) then
                              maxdiff = diff
                              ia=alpha; ik=k; ib=beta; il=l; ig=gamma; im=m; ii=i; ij=j
                           end if
                        end do
                     end do
                  end do
               end do
            end do
         end do

      end do
   end do

   if (maxdiff > thr2) then
      call test_failed(error, "Kernel third derivative does not match finite difference solution")
      print '(a,es20.13)', "Max |d3K/dr3| difference: ", maxdiff
      print '(a,8(i0,1x))', "At indices (a,k,b,l,g,m,i,j) = ", ia,ik,ib,il,ig,im,ii,ij
   end if

end subroutine test_kernel_numt


! subroutine test_kernel_numq(error, mol, kernel_id, keps, qat)
!    type(error_type), allocatable, intent(out) :: error
!    type(structure_type), intent(inout) :: mol
!    integer, intent(in) :: kernel_id
!    real(wp), intent(in) :: keps
!    real(wp), intent(in) :: qat(:)

!    type(born_integrator) :: gbobc
!    class(kernel_type), allocatable :: kernel

!    real(wp), allocatable :: rvdw(:), rad(:)
!    real(wp), allocatable :: draddr(:, :, :)
!    real(wp), allocatable :: draddr2(:, :, :, :, :)
!    real(wp), allocatable :: draddr3(:, :, :, :, :, :, :)
!    real(wp), allocatable :: draddr4(:, :, :, :, :, :, :, :, :)

!    real(wp), allocatable :: d3p(:, :, :, :, :, :)      ! (3,nat,3,nat,3,nat) for one (i,j)
!    real(wp), allocatable :: d3m(:, :, :, :, :, :)      ! (3,nat,3,nat,3,nat) for one (i,j)
!    real(wp), allocatable :: num4(:, :, :, :, :, :, :, :) ! (3,nat,3,nat,3,nat,3,nat) for one (i,j)
!    real(wp), allocatable :: ana4(:, :, :, :, :, :, :, :) ! (3,nat,3,nat,3,nat,3,nat) for one (i,j)

!    real(wp), parameter :: step = 1.0e-6_wp
!    integer :: nat, i, j, k, l, m, n, a, b, c, d
!    real(wp) :: diff, maxdiff
!    integer :: ia, ik, ib, il, ic, im, id, in, ii, ij  ! index record for reporting

!    nat = mol%nat
!    rvdw = get_vdw_rad_d3(mol%num)
!    call new_born_integrator(gbobc, mol, rvdw)
!    kernel = new_kernel(kernel_id, keps)

!    allocate(rad(nat))
!    allocate(draddr(3,nat,nat))
!    allocate(draddr2(3,nat,3,nat,nat))
!    allocate(draddr3(3,nat,3,nat,3,nat,nat))
!    allocate(draddr4(3,nat,3,nat,3,nat,3,nat,nat))

!    allocate(d3p(3,nat,3,nat,3,nat), d3m(3,nat,3,nat,3,nat))
!    allocate(num4(3,nat,3,nat,3,nat,3,nat))
!    allocate(ana4(3,nat,3,nat,3,nat,3,nat))

!    ! reference geometry: need up to brdr4 for analytical
!    call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2, dradd3r=draddr3, dradd4r=draddr4)

!    maxdiff = 0.0_wp
!    ia=1;ik=1;ib=1;il=1;ic=1;im=1;id=1;in=1;ii=1;ij=1

!    do i = 1, nat
!       do j = 1, nat

!          ! ---- analytical 4th-derivative slab for this (i,j)
!          call compute_kernel_d4kdr4_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, draddr3, draddr4, i, j, ana4)

!          ! ---- numerical slab via FD of the ij third derivative
!          num4(:, :, :, :, :, :, :, :) = 0.0_wp

!          do n = 1, nat
!             do d = 1, 3

!                mol%xyz(d, n) = mol%xyz(d, n) + step
!                call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2, dradd3r=draddr3)
!                call compute_kernel_d3kdr3_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, draddr3, i, j, d3p)

!                mol%xyz(d, n) = mol%xyz(d, n) - 2.0_wp*step
!                call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2, dradd3r=draddr3)
!                call compute_kernel_d3kdr3_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, draddr3, i, j, d3m)

!                mol%xyz(d, n) = mol%xyz(d, n) + step

!                do k = 1, nat
!                   do a = 1, 3
!                      do l = 1, nat
!                         do b = 1, 3
!                            do m = 1, nat
!                               do c = 1, 3
!                                  num4(a,k,b,l,c,m,d,n) = 0.5_wp * (d3p(a,k,b,l,c,m) - d3m(a,k,b,l,c,m)) / step
!                               end do
!                            end do
!                         end do
!                      end do
!                   end do
!                end do

!             end do
!          end do

!          ! ---- compare this slab; track worst entry globally
!          do k = 1, nat
!             do a = 1, 3
!                do l = 1, nat
!                   do b = 1, 3
!                      do m = 1, nat
!                         do c = 1, 3
!                            do n = 1, nat
!                               do d = 1, 3
!                                  diff = abs(ana4(a,k,b,l,c,m,d,n) - num4(a,k,b,l,c,m,d,n))
!                                  if (diff > maxdiff) then
!                                     maxdiff = diff
!                                     ia=a; ik=k; ib=b; il=l; ic=c; im=m; id=d; in=n; ii=i; ij=j
!                                     print *, ana4(a,k,b,l,c,m,d,n), num4(a,k,b,l,c,m,d,n)
!                                  end if
!                               end do
!                            end do
!                         end do
!                      end do
!                   end do
!                end do
!             end do
!          end do

!       end do
!    end do

!    if (maxdiff > thr2) then
!       call test_failed(error, "Kernel fourth derivative does not match finite difference solution")
!       print '(a,es20.13)', "Max |d4K/dr4| difference: ", maxdiff
!       print '(a,10(i0,1x))', "At indices (a,k,b,l,c,m,d,n,i,j) = ", ia,ik,ib,il,ic,im,id,in,ii,ij
!    end if

! end subroutine test_kernel_numq


subroutine test_kernel_numq(error, mol, kernel_id, keps, qat)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type), intent(inout) :: mol
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   real(wp), intent(in) :: qat(:)

   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel

   real(wp), allocatable :: rvdw(:), rad(:)
   real(wp), allocatable :: draddr(:, :, :)
   real(wp), allocatable :: draddr2(:, :, :, :, :)
   real(wp), allocatable :: draddr3(:, :, :, :, :, :, :)
   real(wp), allocatable :: draddr4(:, :, :, :, :, :, :, :, :)

   ! 3rd-derivative slabs at shifted coordinates (for one (i,j))
   real(wp), allocatable :: d3pp(:, :, :, :, :, :)   ! x + 2h
   real(wp), allocatable :: d3p (:, :, :, :, :, :)   ! x + h
   real(wp), allocatable :: d3m (:, :, :, :, :, :)   ! x - h
   real(wp), allocatable :: d3mm(:, :, :, :, :, :)   ! x - 2h

   real(wp), allocatable :: num4(:, :, :, :, :, :, :, :) ! (3,nat,3,nat,3,nat,3,nat) for one (i,j)
   real(wp), allocatable :: ana4(:, :, :, :, :, :, :, :) ! (3,nat,3,nat,3,nat,3,nat) for one (i,j)

   real(wp), parameter :: step = 1.0e-6_wp
   integer :: nat, i, j, k, l, m, n, a, b, c, d
   real(wp) :: diff, maxdiff
   real(wp) :: x0
   integer :: ia, ik, ib, il, ic, im, id, in, ii, ij  ! index record for reporting

   nat = mol%nat
   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   kernel = new_kernel(kernel_id, keps)

   allocate(rad(nat))
   allocate(draddr(3,nat,nat))
   allocate(draddr2(3,nat,3,nat,nat))
   allocate(draddr3(3,nat,3,nat,3,nat,nat))
   allocate(draddr4(3,nat,3,nat,3,nat,3,nat,nat))

   allocate(d3pp(3,nat,3,nat,3,nat))
   allocate(d3p (3,nat,3,nat,3,nat))
   allocate(d3m (3,nat,3,nat,3,nat))
   allocate(d3mm(3,nat,3,nat,3,nat))

   allocate(num4(3,nat,3,nat,3,nat,3,nat))
   allocate(ana4(3,nat,3,nat,3,nat,3,nat))

   ! reference geometry: need up to brdr4 for analytical
   call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2, dradd3r=draddr3, dradd4r=draddr4)

   maxdiff = 0.0_wp
   ia=1;ik=1;ib=1;il=1;ic=1;im=1;id=1;in=1;ii=1;ij=1

   do i = 1, nat
      do j = 1, nat

         ! ---- analytical 4th-derivative slab for this (i,j)
         call compute_kernel_d4kdr4_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, draddr3, draddr4, i, j, ana4)

         ! ---- numerical slab via 4-point central FD of the ij third derivative
         num4(:, :, :, :, :, :, :, :) = 0.0_wp

         do n = 1, nat
            do d = 1, 3

               x0 = mol%xyz(d, n)

               ! x + 2h
               mol%xyz(d, n) = x0 + 2.0_wp*step
               call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2, dradd3r=draddr3)
               call compute_kernel_d3kdr3_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, draddr3, i, j, d3pp)

               ! x + h
               mol%xyz(d, n) = x0 + 1.0_wp*step
               call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2, dradd3r=draddr3)
               call compute_kernel_d3kdr3_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, draddr3, i, j, d3p)

               ! x - h
               mol%xyz(d, n) = x0 - 1.0_wp*step
               call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2, dradd3r=draddr3)
               call compute_kernel_d3kdr3_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, draddr3, i, j, d3m)

               ! x - 2h
               mol%xyz(d, n) = x0 - 2.0_wp*step
               call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2, dradd3r=draddr3)
               call compute_kernel_d3kdr3_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, draddr3, i, j, d3mm)

               ! restore
               mol%xyz(d, n) = x0

               ! 4-point (4th-order) central difference for first derivative:
               ! f'(x) ≈ (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / (12h)
               do k = 1, nat
                  do a = 1, 3
                     do l = 1, nat
                        do b = 1, 3
                           do m = 1, nat
                              do c = 1, 3
                                 num4(a,k,b,l,c,m,d,n) = (-d3pp(a,k,b,l,c,m) + 8.0_wp*d3p(a,k,b,l,c,m) &
                                                          -8.0_wp*d3m(a,k,b,l,c,m) + d3mm(a,k,b,l,c,m)) &
                                                          / (12.0_wp*step)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do

            end do
         end do

         ! ---- compare this slab; track worst entry globally
         do k = 1, nat
            do a = 1, 3
               do l = 1, nat
                  do b = 1, 3
                     do m = 1, nat
                        do c = 1, 3
                           do n = 1, nat
                              do d = 1, 3
                                 diff = abs(ana4(a,k,b,l,c,m,d,n) - num4(a,k,b,l,c,m,d,n))
                                 if (diff > maxdiff) then
                                    maxdiff = diff
                                    ia=a; ik=k; ib=b; il=l; ic=c; im=m; id=d; in=n; ii=i; ij=j
                                    print *, ana4(a,k,b,l,c,m,d,n), num4(a,k,b,l,c,m,d,n)
                                 end if
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do

      end do
   end do

   if (maxdiff > thr2) then
      call test_failed(error, "Kernel fourth derivative does not match finite difference solution")
      print '(a,es20.13)', "Max |d4K/dr4| difference: ", maxdiff
      print '(a,10(i0,1x))', "At indices (a,k,b,l,c,m,d,n,i,j) = ", ia,ik,ib,il,ic,im,id,in,ii,ij
   end if

end subroutine test_kernel_numq











!> Test setting up the kernel interacion matrices based on kernel derivatives
!> For this test to work, damping needs to be disabled in the coulomb/multipole.f90
!  Therefore turned off for now.
subroutine test_amat(error, mol, kernel_id, keps, qat, make_multipole, input)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type), intent(inout) :: mol
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   real(wp), intent(in) :: qat(:)

   !> Factory to create new electrostatic objects
   procedure(multipole_maker) :: make_multipole

   type(alpb_input), intent(in) :: input

   type(container_cache) :: cache
   type(coulomb_cache), pointer :: c_cache
   type(damped_multipole) :: d_multipole

   class(kernel_type), allocatable :: kernel
   type(alpb_solvation) :: solv
   type(alpb_input), allocatable :: scratch_input
   type(alpb_cache) :: a_cache



   real(wp), allocatable :: amat_sd_mp(:,:,:), amat_dd_mp(:,:,:,:), amat_sq_mp(:,:,:), amat_dq_mp(:,:,:,:), amat_qq_mp(:,:,:,:)
   real(wp), allocatable :: amat_sd_alpb(:,:,:), amat_dd_alpb(:,:,:,:), amat_sq_alpb(:,:,:), amat_dq_alpb(:,:,:,:), amat_qq_alpb(:,:,:,:)
   real(wp), allocatable :: rad(:), draddr(:,:,:)

   call taint(cache, c_cache)
   call c_cache%update(mol)
   call make_multipole(d_multipole, mol, error)

   scratch_input = input
   call get_alpb_param(scratch_input, mol, 'gfn2', error)
   solv = alpb_solvation(mol, scratch_input, 'gnf2')
   

   allocate(amat_sd_mp(3, mol%nat, mol%nat), source=0.0_wp)
   allocate(amat_dd_mp(3, mol%nat, 3, mol%nat), source=0.0_wp)
   allocate(amat_sq_mp(6, mol%nat, mol%nat), source=0.0_wp)
   allocate(amat_dq_mp(3, mol%nat, 6, mol%nat), source=0.0_wp)
   allocate(amat_qq_mp(6, mol%nat, 6, mol%nat), source=0.0_wp)

   allocate(amat_sd_alpb(3, mol%nat, mol%nat), source=0.0_wp)
   allocate(amat_dd_alpb(3, mol%nat, 3, mol%nat), source=0.0_wp)
   allocate(amat_sq_alpb(6, mol%nat, mol%nat), source=0.0_wp)
   allocate(amat_dq_alpb(3, mol%nat, 6, mol%nat), source=0.0_wp)
   allocate(amat_qq_alpb(6, mol%nat, 6, mol%nat), source=0.0_wp)

   allocate(rad(mol%nat), source=0.0_wp)
   allocate(draddr(3, mol%nat, mol%nat), source=0.0_wp)

   call d_multipole%update(mol, cache)

   amat_sd_mp = c_cache%amat_sd
   amat_dd_mp = c_cache%amat_dd
   amat_sq_mp = c_cache%amat_sq
   amat_dq_mp = c_cache%amat_dq
   amat_qq_mp = c_cache%amat_qq

   call solv%update(mol, cache)


   call get_multipole_matrix(solv, mol, mol%xyz, solv%keps, rad, draddr, &
      & amat_sd_alpb, amat_dd_alpb, amat_sq_alpb, amat_dq_alpb, amat_qq_alpb)

   ! If all fail, this just reports amat_sq
   ! While all maxdiffs are printed, this should be redone some time (if kept in for merge)
    if (maxval(amat_sd_alpb - amat_sd_mp) > thr2) then
      call test_failed(error, "Monopole-dipole interaction matrices do no match!")
      print '(a,es20.13)', "Max difference amat_sd: ", maxval(amat_sd_alpb - amat_sd_mp)
   end if
   if (maxval(amat_dd_alpb - amat_dd_mp) > thr2) then
      call test_failed(error, "Dipole-dipole interaction matrices do no match!")
      print '(a,es20.13)', "Max difference amat_dd: ", maxval(amat_dd_alpb - amat_dd_mp)
   end if
   if (maxval(amat_sq_alpb - amat_sq_mp) > thr2) then
      call test_failed(error, "Monopole-quadrupole interaction matrices do no match!")
      print '(a,es20.13)', "Max difference amat_sq: ", maxval(amat_sq_alpb - amat_sq_mp)
   end if
   if (maxval(amat_dq_alpb - amat_dq_mp) > thr2) then
      call test_failed(error, "Dipole-quadrupole interaction matrices do no match!")
      print '(a,es20.13)', "Max difference amat_dq: ", maxval(amat_dq_alpb - amat_dq_mp)
   end if
   if (maxval(amat_qq_alpb - amat_qq_mp) > thr2) then
      call test_failed(error, "Quadrupole-quadrupole interaction matrices do no match!")
      print *, "Max difference amat_qq: ", maxval(amat_qq_alpb - amat_qq_mp)
   end if


end subroutine test_amat




subroutine test_amat_ho(error, mol)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type), intent(inout) :: mol

   integer :: nat
   real(wp), allocatable :: amat_sd_mp(:,:,:), amat_dd_mp(:,:,:,:), amat_sq_mp(:,:,:)
   real(wp), allocatable :: amat_dq_mp(:,:,:,:), amat_qq_mp(:,:,:,:)
   real(wp), allocatable :: rad(:)

   ! reference blocks for the one pair (jat=2, iat=1)
   real(wp) :: dq_ref(3,6), qq_ref(6,6)

   ! multipole moments for energy checks
   real(wp) :: mu(3,2)     ! dipoles
   real(wp) :: Q6(6,2)     ! packed quadrupoles (xx,2xy,yy,2xz,2yz,zz)

   real(wp) :: E_dq_mp, E_dq_ref
   real(wp) :: E_qq_mp, E_qq_ref

   real(wp) :: R(3)
   real(wp) :: tol
   real(wp) :: maxerr_dq, maxerr_qq

   integer :: iat, jat

   tol = 1.0e-10_wp
   nat = mol%nat

   allocate(rad(nat), source=0.0_wp)

   allocate(amat_sd_mp(3, nat, nat), source=0.0_wp)
   allocate(amat_dd_mp(3, nat, 3, nat), source=0.0_wp)
   allocate(amat_sq_mp(6, nat, nat), source=0.0_wp)
   allocate(amat_dq_mp(3, nat, 6, nat), source=0.0_wp)
   allocate(amat_qq_mp(6, nat, 6, nat), source=0.0_wp)

   ! Build matrices (Coulomb; damping effectively disabled in your current version)
   call get_multipole_matrix_0d(mol, rad, 1.0_wp, 1.0_wp, &
      & amat_sd_mp, amat_dd_mp, amat_sq_mp, amat_dq_mp, amat_qq_mp)

   ! -----------------------
   ! Define arbitrary moments
   ! -----------------------
   mu(:,:) = 0.0_wp
   Q6(:,:) = 0.0_wp

   ! Dipole on atom 2: along +z
   mu(:,2) = [0.0_wp, 0.0_wp, 1.0_wp]

   ! Traceless axial quadrupole on atom 1 and 2 (aligned with z):
   ! Q = diag(-1/2, -1/2, 1)
   Q6(:,1) = [-0.5_wp, 0.0_wp, -0.5_wp, 0.0_wp, 0.0_wp, 1.0_wp]
   Q6(:,2) = [-0.5_wp, 0.0_wp, -0.5_wp, 0.0_wp, 0.0_wp, 1.0_wp]

   ! Choose the pair orientation consistent with your matrix storage:
   ! In your get_multipole_matrix_0d: vec = xyz(:,iat) - xyz(:,jat),
   ! and then stored at ( :, jat, ..., iat ).
   iat = 1
   jat = 2
   R(:) = mol%xyz(:, iat) - mol%xyz(:, jat)

   ! -----------------------
   ! Build analytic references
   ! -----------------------
   call build_coulomb_dq_block(R, dq_ref)   ! 3x6 for (dipole on jat) vs (quad on iat)
   call build_coulomb_qq_block(R, qq_ref)   ! 6x6 for (quad on jat) vs (quad on iat)

   ! -----------------------
   ! Matrix-entry comparisons
   ! -----------------------
   maxerr_dq = maxval(abs(amat_dq_mp(:, jat, :, iat) - dq_ref(:,:)))
   maxerr_qq = maxval(abs(amat_qq_mp(:, jat, :, iat) - qq_ref(:,:)))

   if (maxerr_dq > tol) then
      print *, 'dq', amat_dq_mp(:, jat, :, iat), dq_ref(:,:), maxerr_dq
      return
   end if
   if (maxerr_qq > tol) then
      print *, 'qq', maxerr_qq
      return
   end if

   ! -----------------------
   ! Energy comparisons (same contraction on ref vs mp)
   ! -----------------------
   ! Dipole–quadrupole energy for this pair
   E_dq_mp  = dot_product(mu(:,jat), matmul(amat_dq_mp(:,jat,:,iat), Q6(:,iat)))
   E_dq_ref = dot_product(mu(:,jat), matmul(dq_ref(:,:),                Q6(:,iat)))

   if (abs(E_dq_mp - E_dq_ref) < tol) then
      print *, E_dq_mp, E_dq_ref
      return
   end if

   ! Quadrupole–quadrupole energy for this pair
   E_qq_mp  = dot_product(Q6(:,jat), matmul(amat_qq_mp(:,jat,:,iat), Q6(:,iat)))
   E_qq_ref = dot_product(Q6(:,jat), matmul(qq_ref(:,:),             Q6(:,iat)))

   if (abs(E_qq_mp - E_qq_ref) > tol) then
      call test_failed(error, "Monopole-quadrupole interaction matrices do no match!")
      return
   end if


contains

   subroutine build_coulomb_dq_block(R, dq)
      real(wp), intent(in)  :: R(3)
      real(wp), intent(out) :: dq(3,6)

      real(wp) :: r1, r2, g5, g7
      real(wp) :: I3(3,3)
      real(wp) :: U(3,3,3)
      integer  :: a,b,c

      I3 = 0.0_wp
      I3(1,1)=1.0_wp; I3(2,2)=1.0_wp; I3(3,3)=1.0_wp

      r2 = dot_product(R,R)
      r1 = sqrt(r2)
      g5 = 1.0_wp/(r1**5)
      g7 = 1.0_wp/(r1**7)

      ! U_{a,bc} = -5 R_a R_b R_c / R^7 + (δ_ab R_c + δ_ac R_b + δ_bc R_a) / R^5
      U(:,:,:) = 0.0_wp
      do a=1,3
         do b=1,3
            do c=1,3
               U(a,b,c) = -5.0_wp * R(a)*R(b)*R(c) * g7 &
                        + ( I3(a,b)*R(c) + I3(a,c)*R(b) + I3(b,c)*R(a) ) * g5
            end do
         end do
      end do

      ! Pack bc -> p in your convention: (xx, 2xy, yy, 2xz, 2yz, zz)
      do a=1,3
         dq(a,1) = U(a,1,1)
         dq(a,2) = 2.0_wp*U(a,1,2)
         dq(a,3) = U(a,2,2)
         dq(a,4) = 2.0_wp*U(a,1,3)
         dq(a,5) = 2.0_wp*U(a,2,3)
         dq(a,6) = U(a,3,3)
      end do
   end subroutine build_coulomb_dq_block


   subroutine build_coulomb_qq_block(R, qq)
      real(wp), intent(in)  :: R(3)
      real(wp), intent(out) :: qq(6,6)

      real(wp) :: r1, r2, r4, g9
      real(wp) :: I3(3,3)
      real(wp) :: sym1, sym2, W
      integer  :: p,q,a,b,c,d
      integer, parameter :: pa(6) = [1, 1, 2, 1, 2, 3]
      integer, parameter :: pb(6) = [1, 2, 2, 3, 3, 3]
      integer, parameter :: pf(6) = [1, 2, 1, 2, 2, 1]  ! 1 for diag, 2 for offdiag

      I3 = 0.0_wp
      I3(1,1)=1.0_wp; I3(2,2)=1.0_wp; I3(3,3)=1.0_wp

      r2 = dot_product(R,R)
      r1 = sqrt(r2)
      r4 = r2*r2
      g9 = 1.0_wp/(r1**9)

      ! Traceless-Theta Coulomb QQ tensor:
      ! W_abcd = [ 35 R_a R_b R_c R_d
      !          - 5 R^2 * sym(δ R R)
      !          + (3/2) R^4 * sym(δδ) ] / R^9
      do p=1,6
         a = pa(p); b = pb(p)
         do q=1,6
            c = pa(q); d = pb(q)

            sym1 = I3(a,b)*R(c)*R(d) + I3(a,c)*R(b)*R(d) + I3(a,d)*R(b)*R(c) &
                 + I3(b,c)*R(a)*R(d) + I3(b,d)*R(a)*R(c) + I3(c,d)*R(a)*R(b)

            sym2 = I3(a,b)*I3(c,d) + I3(a,c)*I3(b,d) + I3(a,d)*I3(b,c)

            W = ( 35.0_wp * R(a)*R(b)*R(c)*R(d) &
                -  5.0_wp * r2 * sym1 &
                +  1.5_wp * r4 * sym2 ) * g9

            qq(p,q) = real(pf(p)*pf(q),wp) * W
         end do
      end do
   end subroutine build_coulomb_qq_block


end subroutine test_amat_ho





!> Factory to create electrostatic objects based on GFN2-xTB values
subroutine make_multipole2(multipole, mol, error)

   !> New electrostatic object
   type(damped_multipole), intent(out) :: multipole

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: kdmp3 = 3.0_wp, kdmp5 = 3.0_wp
   real(wp), parameter :: shift = 1.2_wp, kexp = 4.0_wp, rmax = 5.0_wp
   !> Dipole exchange-correlation kernel
   real(wp), parameter :: p_dkernel(20) = 0.01_wp * [&
      & 5.563889_wp,-1.000000_wp,-0.500000_wp,-0.613341_wp,-0.481186_wp, &
      &-0.411674_wp, 3.521273_wp,-4.935670_wp,-8.339183_wp,10.000000_wp, &
      & 0.000000_wp,-0.082005_wp, 2.633341_wp,-0.025750_wp, 2.110225_wp, &
      &-0.151117_wp,-2.536958_wp,-2.077329_wp,-0.103383_wp,-0.236675_wp]
   !> Quadrupole exchange-correlation kernel
   real(wp), parameter :: p_qkernel(20) = 0.01_wp * [&
      & 0.027431_wp,-0.337528_wp, 0.020000_wp,-0.058586_wp,-0.058228_wp, &
      & 0.213583_wp, 2.026786_wp,-0.310828_wp,-0.245955_wp,-0.500000_wp, &
      & 0.020000_wp,-0.005516_wp,-0.021887_wp,-0.080000_wp, 0.028679_wp, &
      & 0.442859_wp, 0.122783_wp,-1.083404_wp, 0.025000_wp, 0.010000_wp]
   real(wp), parameter :: p_rad(20) = [&
      & 1.4_wp, 3.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 3.0_wp, 1.9_wp, 1.8_wp, 2.4_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 3.9_wp, 2.1_wp, 3.1_wp, 2.5_wp, 5.0_wp, 5.0_wp, 5.0_wp]
   real(wp), parameter :: p_vcn(20) = [&
      & 1.0_wp, 1.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 2.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 2.0_wp]
   real(wp), allocatable :: dkernel(:), qkernel(:), rad(:), vcn(:)

   dkernel = p_dkernel(mol%num)
   qkernel = p_qkernel(mol%num)
   rad = p_rad(mol%num)
   vcn = p_vcn(mol%num)

   call new_damped_multipole(multipole, mol, kdmp3, kdmp5, dkernel, qkernel, &
      & shift, kexp, rmax, rad, vcn, error)

end subroutine make_multipole2


!> Inspect container cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(coulomb_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(coulomb_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

!> Return reference to container cache after resolving its type
subroutine view(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(coulomb_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(coulomb_cache)
      ptr => target
   end select
end subroutine view





end module test_solvation_kernel
