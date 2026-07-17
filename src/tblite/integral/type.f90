! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later

!> Storage for integral matrices used during a calculation.
module tblite_integral_type
   use mctc_env, only : wp
   implicit none
   private

   public :: integral_type, new_integral

   !> Cache of integral matrices and the core Hamiltonian.
   type :: integral_type
      real(wp), allocatable :: hamiltonian(:, :)
      real(wp), allocatable :: overlap(:, :)
      real(wp), allocatable :: dipole(:, :, :)
      real(wp), allocatable :: quadrupole(:, :, :)
      real(wp), allocatable :: overlap_diat(:, :)
   end type integral_type

contains

!> Allocate integral storage for a given number of atomic orbitals.
subroutine new_integral(self, nao)
   type(integral_type), intent(out) :: self
   integer, intent(in) :: nao

   allocate(self%hamiltonian(nao, nao), source=0.0_wp)
   allocate(self%overlap(nao, nao), source=0.0_wp)
   allocate(self%dipole(3, nao, nao), source=0.0_wp)
   allocate(self%quadrupole(6, nao, nao), source=0.0_wp)
   allocate(self%overlap_diat(nao, nao), source=0.0_wp)
end subroutine new_integral

end module tblite_integral_type
