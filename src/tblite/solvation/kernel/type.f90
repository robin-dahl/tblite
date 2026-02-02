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

!> @file tblite/solvation/kernel/type.f90
!> Provides the base type definitions and interfaces for generalized Born interaction kernels.

module tblite_solvation_kernel_type
   use mctc_env, only : wp
   implicit none
   private

   public :: kernel_type

   ! Abstract base class for kernel types
   type, abstract :: kernel_type
      real(wp) :: keps
   contains
      procedure(kernel_K_interface),                 deferred :: kernel_K
      procedure(kernel_dKdborn_interface),           deferred :: kernel_dKdborn
      procedure(kernel_dKdr_interface),              deferred :: kernel_dKdr
      procedure(kernel_d_dKdr_dborn_interface),      deferred :: kernel_d_dKdr_dborn
      procedure(kernel_d2Kdr2_interface),            deferred :: kernel_d2Kdr2
      procedure(kernel_d_d2Kdr2_dborn_interface),    deferred :: kernel_d_d2Kdr2_dborn
      procedure(kernel_d3Kdr3_interface),            deferred :: kernel_d3Kdr3
      procedure(kernel_d_d3Kdr3_dborn_interface),    deferred :: kernel_d_d3Kdr3_dborn
      procedure(kernel_d4Kdr4_interface),            deferred :: kernel_d4Kdr4
      procedure(kernel_d_d4Kdr4_dborn_interface),    deferred :: kernel_d_d4Kdr4_dborn
      procedure(kernel_d5Kdr5_interface),            deferred :: kernel_d5Kdr5
   end type kernel_type

   abstract interface
      !> Add basic kernel contributions to interaction matrix
      subroutine kernel_K_interface(self, nat, xyz, brad, amat)
         import :: kernel_type, wp
         !> Instance of interaction kernel
         class(kernel_type), intent(in) :: self
         !> Number of atoms
         integer, intent(in) :: nat
         !> Cartesian coordinates
         real(wp), intent(in) :: xyz(:, :)
         !> Born radii
         real(wp), intent(in) :: brad(:)
         !> Charge-charge interaction amtrix 
         real(wp), intent(inout) :: amat(:, :)
      end subroutine kernel_K_interface

      !> Pairwise kernel gradient (interface) 
      subroutine kernel_dKdr_interface(self, rA, rB, bornA, bornB, d1)
         import :: kernel_type, wp
         !> Instance of interaction kernel
         class(kernel_type), intent(in) :: self
         !> Cartesian coordinates of atom A
         real(wp), intent(in) :: rA(3)
         !> Cartesian coordinates of atom B
         real(wp), intent(in) :: rB(3)
         !> Born radius of atom A
         real(wp), intent(in) :: bornA
         !> Born radius of atom B
         real(wp), intent(in) :: bornB
         !> First derivative of kernel
         real(wp), intent(out) :: d1(3)
      end subroutine kernel_dKdr_interface

      !> Gradient of pairwise kernel wrt Born radii (dK/dborn) 
      subroutine kernel_dKdborn_interface(self, rA, rB, bornA, bornB, dk_bA, dk_bB)
         import :: kernel_type, wp
         !> Instance of interaction kernel
         class(kernel_type), intent(in) :: self
         !> Cartesian coordinates of atom A
         real(wp), intent(in) :: rA(3)
         !> Cartesian coordinates of atom B
         real(wp), intent(in) :: rB(3)
         !> Born radius of atom A
         real(wp), intent(in) :: bornA
         !> Born radius of atom B
         real(wp), intent(in) :: bornB
         !> Derivative wrt Born radius of atom A
         real(wp), intent(out) :: dk_bA
         !> Derivative wrt Born radius of atom B
         real(wp), intent(out) :: dk_bB
      end subroutine kernel_dKdborn_interface

      !> Gradient of pairwise kernel gradient wrt Born radii (d(dK/dr)/dborn) 
      subroutine kernel_d_dKdr_dborn_interface(self, rA, rB, bornA, bornB, d1_bA, d1_bB)
         import :: kernel_type, wp
         !> Instance of interaction kernel
         class(kernel_type), intent(in) :: self
         !> Cartesian coordinates of atom A
         real(wp), intent(in) :: rA(3)
         !> Cartesian coordinates of atom B
         real(wp), intent(in) :: rB(3)
         !> Born radius of atom A
         real(wp), intent(in) :: bornA
         !> Born radius of atom B
         real(wp), intent(in) :: bornB
         !> Derivative wrt Born radius of atom A
         real(wp), intent(out) :: d1_bA(3)
         !> Derivative wrt Born radius of atom B
         real(wp), intent(out) :: d1_bB(3)
      end subroutine kernel_d_dKdr_dborn_interface

      !> Pairwise kernel Hessian  
      subroutine kernel_d2Kdr2_interface(self, rA, rB, bornA, bornB, d2)
         import :: kernel_type, wp
         !> Instance of interaction kernel
         class(kernel_type), intent(in) :: self
         !> Cartesian coordinates of atom A
         real(wp), intent(in) :: rA(3)
         !> Cartesian coordinates of atom B
         real(wp), intent(in) :: rB(3)
         !> Born radius of atom A
         real(wp), intent(in) :: bornA
         !> Born radius of atom B
         real(wp), intent(in) :: bornB
         !> Second derivative of kernel
         real(wp), intent(out) :: d2(3, 3)
      end subroutine kernel_d2Kdr2_interface

      !> Gradient of pairwise kernel Hessian wrt Born radii (d(d2K/dr2)/dborn)
      subroutine kernel_d_d2Kdr2_dborn_interface(self, rA, rB, bornA, bornB, d2_bA, d2_bB)
         import :: kernel_type, wp
         !> Instance of interaction kernel
         class(kernel_type), intent(in) :: self
         !> Cartesian coordinates of atom A
         real(wp), intent(in) :: rA(3)
         !> Cartesian coordinates of atom B
         real(wp), intent(in) :: rB(3)
         !> Born radius of atom A
         real(wp), intent(in) :: bornA
         !> Born radius of atom B
         real(wp), intent(in) :: bornB
         !> Derivative wrt Born radius of atom A
         real(wp), intent(out) :: d2_bA(3, 3)
         !> Derivative wrt Born radius of atom B
         real(wp), intent(out) :: d2_bB(3, 3)
      end subroutine kernel_d_d2Kdr2_dborn_interface

      !> Pairwise kernel third derivative  
      subroutine kernel_d3Kdr3_interface(self, rA, rB, bornA, bornB, d3)
         import :: kernel_type, wp
         !> Instance of interaction kernel
         class(kernel_type), intent(in) :: self
         !> Cartesian coordinates of atom A
         real(wp), intent(in) :: rA(3)
         !> Cartesian coordinates of atom B
         real(wp), intent(in) :: rB(3)
         !> Born radius of atom A
         real(wp), intent(in) :: bornA
         !> Born radius of atom B
         real(wp), intent(in) :: bornB
         !> Third derivative of kernel
         real(wp), intent(out) :: d3(3, 3, 3)
      end subroutine kernel_d3Kdr3_interface

      !> Gradient of pairwise kernel third derivative wrt Born radii (d(d3K/dr3)/dborn) 
      subroutine kernel_d_d3Kdr3_dborn_interface(self, rA, rB, bornA, bornB, d3_bA, d3_bB)
         import :: kernel_type, wp
         !> Instance of interaction kernel
         class(kernel_type), intent(in) :: self
         !> Cartesian coordinates of atom A
         real(wp), intent(in) :: rA(3)
         !> Cartesian coordinates of atom B
         real(wp), intent(in) :: rB(3)
         !> Born radius of atom A
         real(wp), intent(in) :: bornA
         !> Born radius of atom B
         real(wp), intent(in) :: bornB
         !> Derivative wrt Born radius of atom A
         real(wp), intent(out) :: d3_bA(3, 3, 3)
         !> Derivative wrt Born radius of atom B
         real(wp), intent(out) :: d3_bB(3, 3, 3)
      end subroutine kernel_d_d3Kdr3_dborn_interface

      !> Pairwise kernel fourth derivative  
      subroutine kernel_d4Kdr4_interface(self, rA, rB, bornA, bornB, d4)
         import :: kernel_type, wp
         !> Instance of interaction kernel
         class(kernel_type), intent(in) :: self
         !> Cartesian coordinates of atom A
         real(wp), intent(in) :: rA(3)
         !> Cartesian coordinates of atom B
         real(wp), intent(in) :: rB(3)
         !> Born radius of atom A
         real(wp), intent(in) :: bornA
         !> Born radius of atom B
         real(wp), intent(in) :: bornB
         !> Fourth derivative of kernel
         real(wp), intent(out) :: d4(3, 3, 3, 3)
      end subroutine kernel_d4Kdr4_interface

      !> Gradient of pairwise kernel fourth derivative wrt Born radii (d(d4K/dr4)/dborn) 
      subroutine kernel_d_d4Kdr4_dborn_interface(self, rA, rB, bornA, bornB, d4_bA, d4_bB)
         import :: kernel_type, wp
         !> Instance of interaction kernel
         class(kernel_type), intent(in) :: self
         !> Cartesian coordinates of atom A
         real(wp), intent(in) :: rA(3)
         !> Cartesian coordinates of atom B
         real(wp), intent(in) :: rB(3)
         !> Born radius of atom A
         real(wp), intent(in) :: bornA
         !> Born radius of atom B
         real(wp), intent(in) :: bornB
         !> Derivative wrt Born radius of atom A
         real(wp), intent(out) :: d4_bA(3, 3, 3, 3)
         !> Derivative wrt Born radius of atom B
         real(wp), intent(out) :: d4_bB(3, 3, 3, 3)
      end subroutine kernel_d_d4Kdr4_dborn_interface

      !> Pairwise kernel fifth derivative 
      subroutine kernel_d5Kdr5_interface(self, rA, rB, bornA, bornB, d5)
         import :: kernel_type, wp
         !> Instance of interaction kernel
         class(kernel_type), intent(in) :: self
         !> Cartesian coordinates of atom A
         real(wp), intent(in) :: rA(3)
         !> Cartesian coordinates of atom B
         real(wp), intent(in) :: rB(3)
         !> Born radius of atom A
         real(wp), intent(in) :: bornA
         !> Born radius of atom B
         real(wp), intent(in) :: bornB
         !> Fifth derivative of kernel
         real(wp), intent(out) :: d5(3, 3, 3, 3, 3)
      end subroutine kernel_d5Kdr5_interface
      
   end interface

end module tblite_solvation_kernel_type
