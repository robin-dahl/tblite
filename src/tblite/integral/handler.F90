! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later

!> Abstract interface for Gaussian integral evaluators.
module tblite_integral_handler
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type, cgto_type
   implicit none
   private

   public :: integral_handler

   !> Common interface for evaluating shell-pair integral blocks.
   type, abstract :: integral_handler
   contains
      procedure(initialize_interface), deferred :: initialize_integral
      procedure(multipole_interface), deferred :: multipole_cgto
      procedure(multipole_gradient_interface), deferred :: multipole_grad_cgto
      procedure(dipole_interface), deferred :: dipole_cgto
   end type integral_handler

   abstract interface
      !> Initialize implementation-specific data for a structure and basis.
      subroutine initialize_interface(self, mol, basis)
         import :: integral_handler, structure_type, basis_type
         !> Integral evaluator
         class(integral_handler), intent(inout) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Basis set information
         type(basis_type), intent(in) :: basis
      end subroutine initialize_interface

      !> Evaluate overlap, dipole and quadrupole integrals for a shell pair.
      !>
      !> The CGTOs provide the shell data used by native evaluators, while the
      !> matching global shell indices identify cached integral representations.
      subroutine multipole_interface(self, cgtoj, cgtoi, jsh, ish, r2, vec, intcut, overlap, &
            & dipole, quadrupole)
         import :: integral_handler, cgto_type, wp
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
         real(wp), intent(out) :: overlap(:)
         !> Dipole moment integrals for the given pair i and j
         real(wp), intent(out) :: dipole(:, :)
         !> Quadrupole moment integrals for the given pair i and j
         real(wp), intent(out) :: quadrupole(:, :)
      end subroutine multipole_interface

      !> Evaluate multipole integrals and their nuclear derivatives.
      subroutine multipole_gradient_interface(self, cgtoj, cgtoi, jsh, ish, r2, vec, intcut, &
            & overlap, dipole, quadrupole, doverlap, ddipole_j, dquadrupole_j, &
            & ddipole_i, dquadrupole_i)
         import :: integral_handler, cgto_type, wp
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
         real(wp), intent(out) :: overlap(:)
         !> Dipole moment integrals for the given pair i and j
         real(wp), intent(out) :: dipole(:, :)
         !> Quadrupole moment integrals for the given pair i and j
         real(wp), intent(out) :: quadrupole(:, :)
         !> Overlap integral gradient for the given pair i and j
         real(wp), intent(out) :: doverlap(:, :)
         !> Dipole moment integral gradient with respect to center j
         real(wp), intent(out) :: ddipole_j(:, :, :)
         !> Quadrupole moment integral gradient with respect to center j
         real(wp), intent(out) :: dquadrupole_j(:, :, :)
         !> Dipole moment integral gradient with respect to center i
         real(wp), intent(out) :: ddipole_i(:, :, :)
         !> Quadrupole moment integral gradient with respect to center i
         real(wp), intent(out) :: dquadrupole_i(:, :, :)
      end subroutine multipole_gradient_interface

      !> Evaluate overlap and dipole integrals for a shell pair.
      subroutine dipole_interface(self, cgtoj, cgtoi, jsh, ish, r2, vec, intcut, overlap, dipole)
         import :: integral_handler, cgto_type, wp
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
         real(wp), intent(out) :: overlap(:)
         !> Dipole moment integrals for the given pair i and j
         real(wp), intent(out) :: dipole(:, :)
      end subroutine dipole_interface
   end interface

end module tblite_integral_handler
