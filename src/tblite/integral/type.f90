! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later

!> Integral storage and backend interface.
module tblite_integral_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   implicit none
   private

   public :: integral_type, new_integral

   !> Integral container and abstract evaluator interface.
   type, abstract :: integral_type
      real(wp), allocatable :: hamiltonian(:, :)
      real(wp), allocatable :: overlap(:, :)
      real(wp), allocatable :: dipole(:, :, :)
      real(wp), allocatable :: quadrupole(:, :, :)
      real(wp), allocatable :: overlap_diat(:, :)
   contains
      procedure(initialize_interface), deferred :: initialize_integral
      procedure(multipole_interface), deferred :: multipole_integral
      procedure(multipole_gradient_interface), deferred :: multipole_gradient_integral
      procedure(dipole_interface), deferred :: dipole_integral
   end type integral_type

   abstract interface
      subroutine initialize_interface(self, mol, basis)
         import :: integral_type, structure_type, basis_type
         class(integral_type), intent(inout) :: self
         type(structure_type), intent(in) :: mol
         type(basis_type), intent(in) :: basis
      end subroutine initialize_interface

      subroutine multipole_interface(self, mol, basis, jsh, ish, r2, vec, overlap, &
            & dipole, quadrupole)
         import :: integral_type, structure_type, basis_type, wp
         class(integral_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(basis_type), intent(in) :: basis
         integer, intent(in) :: jsh, ish
         real(wp), intent(in) :: r2, vec(3)
         real(wp), intent(out) :: overlap(:), dipole(:, :), quadrupole(:, :)
      end subroutine multipole_interface

      subroutine multipole_gradient_interface(self, mol, basis, jsh, ish, r2, vec, &
            & overlap, dipole, quadrupole, doverlap, ddipole_j, dquadrupole_j, &
            & ddipole_i, dquadrupole_i)
         import :: integral_type, structure_type, basis_type, wp
         class(integral_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(basis_type), intent(in) :: basis
         integer, intent(in) :: jsh, ish
         real(wp), intent(in) :: r2, vec(3)
         real(wp), intent(out) :: overlap(:), dipole(:, :), quadrupole(:, :)
         real(wp), intent(out) :: doverlap(:, :)
         real(wp), intent(out) :: ddipole_j(:, :, :), dquadrupole_j(:, :, :)
         real(wp), intent(out) :: ddipole_i(:, :, :), dquadrupole_i(:, :, :)
      end subroutine multipole_gradient_interface

      subroutine dipole_interface(self, mol, basis, jsh, ish, r2, vec, overlap, dipole)
         import :: integral_type, structure_type, basis_type, wp
         class(integral_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(basis_type), intent(in) :: basis
         integer, intent(in) :: jsh, ish
         real(wp), intent(in) :: r2, vec(3)
         real(wp), intent(out) :: overlap(:), dipole(:, :)
      end subroutine dipole_interface
   end interface

contains

!> Clone a calculator's selected backend and allocate its matrix storage.
subroutine new_integral(self, backend, nao)
   class(integral_type), allocatable, intent(out) :: self
   class(integral_type), intent(in) :: backend
   integer, intent(in) :: nao

   allocate(self, source=backend)
   allocate(self%hamiltonian(nao, nao), source=0.0_wp)
   allocate(self%overlap(nao, nao), source=0.0_wp)
   allocate(self%dipole(3, nao, nao), source=0.0_wp)
   allocate(self%quadrupole(6, nao, nao), source=0.0_wp)
   allocate(self%overlap_diat(nao, nao), source=0.0_wp)
end subroutine new_integral

end module tblite_integral_type
