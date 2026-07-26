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

!> @file tblite/integral/shell.f90
!> Provides dimensions and component ordering for Gaussian shells.

!> Shell dimensions and Cartesian component ordering used by integral routines.
module tblite_integral_shell
   implicit none
   private

   public :: maxl, maxl2
   public :: msao, mlao, smap, lmap, sdim, lx

   !> Maximum supported angular momentum
   integer, parameter :: maxl = 6
   !> Maximum combined angular momentum
   integer, parameter :: maxl2 = maxl*2
   !> Number of spherical functions for each angular momentum
   integer, parameter :: msao(0:maxl) = [1, 3, 5, 7, 9, 11, 13]
   !> Number of Cartesian functions for each angular momentum
   integer, parameter :: mlao(0:maxl) = [1, 3, 6, 10, 15, 21, 28]
   !> Offset of each angular momentum in a spherical shell block
   integer, parameter :: smap(0:maxl) = [0, 1, 4, 9, 16, 25, 36]
   !> Offset of each angular momentum in the Cartesian exponent table
   integer, parameter :: lmap(0:maxl) = [0, 1, 4, 10, 20, 35, 56]
   !> Size of a spherical shell block through each angular momentum
   integer, parameter :: sdim(0:maxl) = [1, 4, 9, 16, 25, 36, 49]

   !> Cartesian exponents in CCA ordering. For angular momentum l, components
   !> are ordered by decreasing lx and, for equal lx, by decreasing ly, with
   !> lz = l - lx - ly.
   integer, parameter :: lx(3, 84) = reshape([&
      & 0, &
      & 1,0,0, &
      & 2,1,1,0,0,0, &
      & 3,2,2,1,1,1,0,0,0,0, &
      & 4,3,3,2,2,2,1,1,1,1,0,0,0,0,0, &
      & 5,4,4,3,3,3,2,2,2,2,1,1,1,1,1,0,0,0,0,0,0, &
      & 6,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1,0,0,0,0,0,0,0, &
      & 0, &
      & 0,1,0, &
      & 0,1,0,2,1,0, &
      & 0,1,0,2,1,0,3,2,1,0, &
      & 0,1,0,2,1,0,3,2,1,0,4,3,2,1,0, &
      & 0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0, &
      & 0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0,6,5,4,3,2,1,0, &
      & 0, &
      & 0,0,1, &
      & 0,0,1,0,1,2, &
      & 0,0,1,0,1,2,0,1,2,3, &
      & 0,0,1,0,1,2,0,1,2,3,0,1,2,3,4, &
      & 0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,5, &
      & 0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,5,0,1,2,3,4,5,6], &
      & shape(lx), order=[2, 1])

end module tblite_integral_shell
