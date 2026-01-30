#!/bin/bash
# Script to split kernel.f90 into separate files

cd "$(dirname "$0")"

echo "Splitting kernel implementations..."

# Extract Still kernel (lines 346-1838)
{
cat << 'HEADER'
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

!> @file tblite/solvation/kernel/still.f90
!> Provides the Still kernel implementation for GBSA solvation.

module tblite_solvation_kernel_still
   use mctc_env, only : wp
   use tblite_blas, only : gemv
   use tblite_solvation_kernel_kernel, only : kernel_type
   implicit none
   private

   public :: still_kernel

   type, extends(kernel_type) :: still_kernel
   contains
      procedure :: add_kernel_mat   => add_still_mat
      procedure :: add_kernel_deriv => add_still_deriv
      procedure :: compute_kernel_dKdr => compute_still_dKdr
      procedure :: compute_kernel_d2Kdr2 => compute_still_d2Kdr2
      procedure :: compute_kernel_d3Kdr3_ij => compute_still_d3Kdr3_ij
      procedure :: compute_kernel_d3Kdr3 => compute_still_d3Kdr3
      procedure :: compute_kernel_d4Kdr4 => compute_still_d4Kdr4
   end type still_kernel

contains

HEADER
sed -n '346,1838p' kernel.f90.bak
echo ""
echo "end module tblite_solvation_kernel_still"
} > kernel/still.f90

echo "Created kernel/still.f90"

# Extract P16 kernel (lines 1839-3617)
{
cat << 'HEADER'
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

!> @file tblite/solvation/kernel/p16.f90
!> Provides the P16 kernel implementation for ALPB solvation.

module tblite_solvation_kernel_p16
   use mctc_env, only : wp
   use tblite_blas, only : gemv
   use tblite_solvation_kernel_kernel, only : kernel_type
   implicit none
   private

   public :: p16_kernel

   type, extends(kernel_type) :: p16_kernel
   contains
      procedure :: add_kernel_mat   => add_p16_mat
      procedure :: add_kernel_deriv => add_p16_deriv
      procedure :: compute_kernel_dKdr => compute_p16_dKdr
      procedure :: compute_kernel_d2Kdr2 => compute_p16_d2Kdr2
      procedure :: compute_kernel_d3Kdr3 => compute_p16_d3Kdr3
      procedure :: compute_kernel_d3Kdr3_ij => compute_p16_d3Kdr3_ij
      procedure :: compute_kernel_d4Kdr4 => compute_p16_d4Kdr4
   end type p16_kernel

   real(wp), parameter :: zetaP16    = 1.028_wp
   real(wp), parameter :: zetaP16o16 = zetaP16 / 16.0_wp

contains

HEADER
sed -n '1839,3617p' kernel.f90.bak
echo ""
echo "end module tblite_solvation_kernel_p16"
} > kernel/p16.f90

echo "Created kernel/p16.f90"

# Extract Coulomb kernel (lines 3618-4107)
{
cat << 'HEADER'
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

!> @file tblite/solvation/kernel/coulomb.f90
!> Provides the Coulomb kernel implementation.

module tblite_solvation_kernel_coulomb
   use mctc_env, only : wp
   use tblite_blas, only : gemv
   use tblite_solvation_kernel_kernel, only : kernel_type
   implicit none
   private

   public :: coulomb_kernel

   type, extends(kernel_type) :: coulomb_kernel
   contains
      procedure :: add_kernel_mat => add_coulomb_mat
      procedure :: add_kernel_deriv => add_coulomb_deriv
      procedure :: compute_kernel_dKdr => compute_coulomb_dKdr
      procedure :: compute_kernel_d2Kdr2 => compute_coulomb_d2Kdr2
      procedure :: compute_kernel_d3Kdr3_ij => compute_coulomb_d3Kdr3_ij
      procedure :: compute_kernel_d3Kdr3 => compute_coulomb_d3Kdr3
      procedure :: compute_kernel_d4Kdr4 => compute_coulomb_d4Kdr4
   end type

contains

HEADER
sed -n '3618,4107p' kernel.f90.bak
echo ""
echo "end module tblite_solvation_kernel_coulomb"
} > kernel/coulomb.f90

echo "Created kernel/coulomb.f90"

echo "Done! Created 3 kernel implementation files in kernel/ directory."
echo "Files created:"
echo "  - kernel/kernel.f90 (base types and interfaces)"
echo "  - kernel/still.f90"
echo "  - kernel/p16.f90"
echo "  - kernel/coulomb.f90"
echo ""
echo "Next steps:"
echo "1. Create new main kernel.f90 wrapper module"
echo "2. Update meson.build or CMakeLists.txt"
