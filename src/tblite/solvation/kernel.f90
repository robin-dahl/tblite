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

!> @file tblite/solvation/kernel.f90
!> Provides the generalized Born interaction kernels used in the ALPB and GBSA implicit solvation models.
!> This module serves as the main interface, re-exporting types and functions from submodules.

module tblite_solvation_kernel
   use mctc_env, only : wp
   use tblite_solvation_kernel_type, only : kernel_type
   use tblite_solvation_kernel_still, only : still_kernel
   use tblite_solvation_kernel_p16, only : p16_kernel
   use tblite_solvation_kernel_coulomb, only : coulomb_kernel
   implicit none
   private

   public :: kernel_type, new_kernel
   public :: kernel_enum, kernel_enum_type

   !> Possible interaction kernels
   type :: kernel_enum_type
      !> Still interaction kernel
      integer :: still = 1
      !> P16 interaction kernel
      integer :: p16   = 2
      !> Classic Coulomb interaction kernel
      integer :: coulomb = 3
   end type kernel_enum_type

   !> Actual enumerator for possible interaction kernels
   type(kernel_enum_type), parameter :: kernel_enum = kernel_enum_type()

contains

function new_kernel(kernel_id, keps) result(kernel)
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   class(kernel_type), allocatable :: kernel

   select case(kernel_id)
   case(kernel_enum%still)
      allocate(still_kernel :: kernel)
   case(kernel_enum%p16)
      allocate(p16_kernel :: kernel)
   case(kernel_enum%coulomb)
      allocate(coulomb_kernel :: kernel)
   case default
      allocate(p16_kernel :: kernel)
   end select

   kernel%keps = keps
end function new_kernel


end module tblite_solvation_kernel
