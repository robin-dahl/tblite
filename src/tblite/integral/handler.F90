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

!> @dir tblite/integral
!> Contains the integral evaluation implementations

!> @file tblite/integral/handler.f90
!> Provides an abstract interface for Gaussian integral evaluators,
!> either native or via libcint, and their nuclear derivatives.

!> Abstract interface for Gaussian integral evaluators.
module tblite_integral_handler
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type, cgto_type
   implicit none
   private

   public :: enum_integral_handler, integral_handler, msao

   integer, parameter :: maxl = 6
   integer, parameter :: msao(0:maxl) = [1, 3, 5, 7, 9, 11, 13]

   !> Possible Gaussian integral evaluation implementations
   type :: enum_integral_handler_type
      !> Native integral evaluator
      integer :: native = 1
      !> Libcint integral evaluator
      integer :: libcint = 2
   end type enum_integral_handler_type

   type(enum_integral_handler_type), parameter :: enum_integral_handler = &
      & enum_integral_handler_type()

   !> Common interface for evaluating shell-pair integral blocks.
   type, abstract :: integral_handler
   contains
      procedure(initialize_integral), deferred :: initialize_integral
      procedure(multipole_cgto), deferred :: multipole_cgto
      procedure(multipole_grad_cgto), deferred :: multipole_grad_cgto
      procedure(dipole_cgto), deferred :: dipole_cgto
   end type integral_handler

   abstract interface
      !> Initialize implementation-specific data for a structure and basis.
      subroutine initialize_integral(self, mol, basis)
         import :: integral_handler, structure_type, basis_type
         !> Integral evaluator
         class(integral_handler), intent(inout) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Basis set information
         type(basis_type), intent(in) :: basis
      end subroutine initialize_integral

      !> Evaluate overlap, dipole and quadrupole integrals for a shell pair.
      !>
      !> The CGTOs provide the shell data used by native evaluators, while the
      !> matching global shell indices identify cached integral representations.
      subroutine multipole_cgto(self, cgtoj, cgtoi, jsh, ish, r2, vec, intcut, overlap, &
            & dpint, qpint)
         import :: integral_handler, cgto_type, msao, wp
         !> Integral evaluator
         class(integral_handler), intent(in) :: self
         !> Description of contracted Gaussian function on center j
         type(cgto_type), intent(in) :: cgtoj
         !> Description of contracted Gaussian function on center i
         type(cgto_type), intent(in) :: cgtoi
         !> Global shell index of the contracted Gaussian function on center j
         integer, intent(in) :: jsh
         !> Global shell index of the contracted Gaussian function on center i
         integer, intent(in) :: ish
         !> Square distance between center i and j
         real(wp), intent(in) :: r2
         !> Distance vector between center i and j, ri - rj
         real(wp), intent(in) :: vec(3)
         !> Maximum value of integral prefactor to consider
         real(wp), intent(in) :: intcut
         !> Overlap integrals for the given pair i and j
         real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))
         !> Dipole moment integrals for the given pair i and j
         real(wp), intent(out) :: dpint(3, msao(cgtoj%ang), msao(cgtoi%ang))
         !> Quadrupole moment integrals for the given pair i and j
         real(wp), intent(out) :: qpint(6, msao(cgtoj%ang), msao(cgtoi%ang))
      end subroutine multipole_cgto

      !> Evaluate multipole integrals and their nuclear derivatives.
      subroutine multipole_grad_cgto(self, cgtoj, cgtoi, jsh, ish, r2, vec, intcut, &
            & overlap, dpint, qpint, doverlap, ddpintj, dqpintj, &
            & ddpinti, dqpinti)
         import :: integral_handler, cgto_type, msao, wp
         !> Integral evaluator
         class(integral_handler), intent(in) :: self
         !> Description of contracted Gaussian function on center j
         type(cgto_type), intent(in) :: cgtoj
         !> Description of contracted Gaussian function on center i
         type(cgto_type), intent(in) :: cgtoi
         !> Global shell index of the contracted Gaussian function on center j
         integer, intent(in) :: jsh
         !> Global shell index of the contracted Gaussian function on center i
         integer, intent(in) :: ish
         !> Square distance between center i and j
         real(wp), intent(in) :: r2
         !> Distance vector between center i and j, ri - rj
         real(wp), intent(in) :: vec(3)
         !> Maximum value of integral prefactor to consider
         real(wp), intent(in) :: intcut
         !> Overlap integrals for the given pair i and j
         real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))
         !> Dipole moment integrals for the given pair i and j
         real(wp), intent(out) :: dpint(3, msao(cgtoj%ang), msao(cgtoi%ang))
         !> Quadrupole moment integrals for the given pair i and j
         real(wp), intent(out) :: qpint(6, msao(cgtoj%ang), msao(cgtoi%ang))
         !> Overlap integral gradient for the given pair i and j
         real(wp), intent(out) :: doverlap(3, msao(cgtoj%ang), msao(cgtoi%ang))
         !> Dipole moment integral gradient with respect to center j
         real(wp), intent(out) :: ddpintj(3, 3, msao(cgtoj%ang), msao(cgtoi%ang))
         !> Quadrupole moment integral gradient with respect to center j
         real(wp), intent(out) :: dqpintj(3, 6, msao(cgtoj%ang), msao(cgtoi%ang))
         !> Dipole moment integral gradient with respect to center i
         real(wp), intent(out) :: ddpinti(3, 3, msao(cgtoj%ang), msao(cgtoi%ang))
         !> Quadrupole moment integral gradient with respect to center i
         real(wp), intent(out) :: dqpinti(3, 6, msao(cgtoj%ang), msao(cgtoi%ang))
      end subroutine multipole_grad_cgto

      !> Evaluate overlap and dipole integrals for a shell pair.
      subroutine dipole_cgto(self, cgtoj, cgtoi, jsh, ish, r2, vec, intcut, &
         & overlap, dpint)
         import :: integral_handler, cgto_type, msao, wp
         !> Integral evaluator
         class(integral_handler), intent(in) :: self
         !> Description of contracted Gaussian function on center j
         type(cgto_type), intent(in) :: cgtoj
         !> Description of contracted Gaussian function on center i
         type(cgto_type), intent(in) :: cgtoi
         !> Global shell index of the contracted Gaussian function on center j
         integer, intent(in) :: jsh
         !> Global shell index of the contracted Gaussian function on center i
         integer, intent(in) :: ish
         !> Square distance between center i and j
         real(wp), intent(in) :: r2
         !> Distance vector between center i and j, ri - rj
         real(wp), intent(in) :: vec(3)
         !> Maximum value of integral prefactor to consider
         real(wp), intent(in) :: intcut
         !> Overlap integrals for the given pair i and j
         real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))
         !> Dipole moment integrals for the given pair i and j
         real(wp), intent(out) :: dpint(3, msao(cgtoj%ang), msao(cgtoi%ang))
      end subroutine dipole_cgto
   end interface

end module tblite_integral_handler
