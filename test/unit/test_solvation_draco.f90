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

module test_solvation_draco
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use tblite_solvation_data, only : get_vdw_rad_cosmo
   use tblite_solvation_radii_scaling, only : draco
   implicit none
   private

   public :: collect_solvation_draco

   real(wp), parameter :: thr = 5.0e-5_wp

contains


!> Collect all exported unit tests
subroutine collect_solvation_draco(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("gradient-water-cosmo", test_g_water_cosmo) &
      ]

end subroutine collect_solvation_draco


subroutine test_g_water_cosmo(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   integer, parameter :: num(3) = [8, 1, 1]
   real(wp), parameter :: xyz(3, 3) = reshape([&
      &  0.000000000000_wp,  0.000000000000_wp,  0.000000000000_wp, &
      &  1.515263215189_wp,  0.000000000000_wp, -1.058898509481_wp, &
      & -1.515263215189_wp,  0.000000000000_wp, -1.058898509481_wp], [3, 3])
   real(wp), allocatable :: radii_in(:), radii_out(:), drdr(:, :, :), numg(:, :, :)
   real(wp), allocatable :: right(:), left(:)
   real(wp) :: step, maxdiff
   integer :: iat, jat, ic

   call new(mol, num, xyz)

   allocate(radii_in(mol%nat), radii_out(mol%nat), right(mol%nat), left(mol%nat))
   allocate(drdr(3, mol%nat, mol%nat), numg(3, mol%nat, mol%nat))

   do iat = 1, mol%nat
      radii_in(iat) = get_vdw_rad_cosmo(mol%num(mol%id(iat)))
   end do

   call draco(mol, radii_in, radii_out, "water", "cosmo", drdr=drdr)

   step = 1.0e-5_wp
   do jat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, jat) = mol%xyz(ic, jat) + step
         call draco(mol, radii_in, right, "water", "cosmo")

         mol%xyz(ic, jat) = mol%xyz(ic, jat) - 2*step
         call draco(mol, radii_in, left, "water", "cosmo")

         mol%xyz(ic, jat) = mol%xyz(ic, jat) + step

         numg(ic, jat, :) = 0.5_wp * (right - left) / step
      end do
   end do

   maxdiff = maxval(abs(numg - drdr))
   if (maxdiff > thr) then
      call test_failed(error, "DRACO radii derivative does not match finite difference")
      print *, "max difference:", maxdiff
      print *, "analytical:"
      do iat = 1, mol%nat
         print '(3es20.13)', drdr(:, :, iat)
      end do
      print *, "numerical:"
      do iat = 1, mol%nat
         print '(3es20.13)', numg(:, :, iat)
      end do
      print *, "difference:"
      do iat = 1, mol%nat
         print '(3es20.13)', drdr(:, :, iat) - numg(:, :, iat)
      end do
      return
   end if

   call check(error, 0.0_wp, maxdiff, thr=thr)
end subroutine test_g_water_cosmo


end module test_solvation_draco
