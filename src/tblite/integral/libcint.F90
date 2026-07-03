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

!> @file tblite/integral/libcint.f90
!> Provides an optional libcint-backed Gaussian integral interface.
module tblite_integral_libcint
   use, intrinsic :: iso_c_binding, only : c_double, c_int, c_null_ptr, c_ptr
   use mctc_env, only : wp
   implicit none
   private

   public :: LIBCINT_1E_OVERLAP, LIBCINT_1E_KINETIC, LIBCINT_1E_NUCLEAR
   public :: LIBCINT_CARTESIAN, LIBCINT_SPHERICAL
   public :: CHARGE_OF, PTR_COORD, NUC_MOD_OF, PTR_ZETA, ATM_SLOTS
   public :: ATOM_OF, ANG_OF, NPRIM_OF, NCTR_OF, KAPPA_OF, PTR_EXP, PTR_COEFF, BAS_SLOTS
   public :: PTR_ENV_START
   public :: libcint_cgto_cart, libcint_cgto_spheric
   public :: libcint_tot_cgto_cart, libcint_tot_cgto_spheric
   public :: libcint_gto_norm
   public :: libcint_shell_size
   public :: libcint_eval_1e
   public :: libcint_eval_eri

   ! libcint C arrays are flattened as slot + slots * item.  These Fortran
   ! constants are shifted by one so normal arrays can be declared as
   ! atm(ATM_SLOTS,natm), bas(BAS_SLOTS,nbas).
   integer, parameter :: CHARGE_OF = 1
   integer, parameter :: PTR_COORD = 2
   integer, parameter :: NUC_MOD_OF = 3
   integer, parameter :: PTR_ZETA = 4
   integer, parameter :: ATM_SLOTS = 6

   integer, parameter :: ATOM_OF = 1
   integer, parameter :: ANG_OF = 2
   integer, parameter :: NPRIM_OF = 3
   integer, parameter :: NCTR_OF = 4
   integer, parameter :: KAPPA_OF = 5
   integer, parameter :: PTR_EXP = 6
   integer, parameter :: PTR_COEFF = 7
   integer, parameter :: BAS_SLOTS = 8

   ! env offsets stored in atm/bas are zero-based libcint offsets.  Fortran
   ! code should write env(offset+1:offset+n) for C locations offset:offset+n-1.
   integer, parameter :: PTR_ENV_START = 20

   integer, parameter :: LIBCINT_1E_OVERLAP = 1
   integer, parameter :: LIBCINT_1E_KINETIC = 2
   integer, parameter :: LIBCINT_1E_NUCLEAR = 3
   integer, parameter :: LIBCINT_CARTESIAN = 1
   integer, parameter :: LIBCINT_SPHERICAL = 2

   interface
      function cint_cgto_cart(bas_id, bas) bind(C, name="CINTcgto_cart") result(nao)
         import :: c_int
         integer(c_int), value :: bas_id
         integer(c_int), intent(in) :: bas(*)
         integer(c_int) :: nao
      end function cint_cgto_cart

      function cint_cgto_spheric(bas_id, bas) bind(C, name="CINTcgto_spheric") result(nao)
         import :: c_int
         integer(c_int), value :: bas_id
         integer(c_int), intent(in) :: bas(*)
         integer(c_int) :: nao
      end function cint_cgto_spheric

      function cint_tot_cgto_cart(bas, nbas) bind(C, name="CINTtot_cgto_cart") result(nao)
         import :: c_int
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value :: nbas
         integer(c_int) :: nao
      end function cint_tot_cgto_cart

      function cint_tot_cgto_spheric(bas, nbas) bind(C, name="CINTtot_cgto_spheric") result(nao)
         import :: c_int
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value :: nbas
         integer(c_int) :: nao
      end function cint_tot_cgto_spheric

      function cint_gto_norm(n, alpha) bind(C, name="CINTgto_norm") result(norm)
         import :: c_double, c_int
         integer(c_int), value :: n
         real(c_double), value :: alpha
         real(c_double) :: norm
      end function cint_gto_norm

      function int1e_ovlp_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_ovlp_cart") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         type(c_ptr), value :: dims
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_ovlp_cart

      function int1e_ovlp_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_ovlp_sph") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         type(c_ptr), value :: dims
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_ovlp_sph

      function int1e_kin_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_kin_cart") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         type(c_ptr), value :: dims
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_kin_cart

      function int1e_kin_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_kin_sph") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         type(c_ptr), value :: dims
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_kin_sph

      function int1e_nuc_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_nuc_cart") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         type(c_ptr), value :: dims
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_nuc_cart

      function int1e_nuc_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_nuc_sph") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         type(c_ptr), value :: dims
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_nuc_sph

      function cint2e_cart(out, shls, atm, natm, bas, nbas, env, opt) &
            & bind(C, name="cint2e_cart") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt
         integer(c_int) :: stat
      end function cint2e_cart

      function cint2e_sph(out, shls, atm, natm, bas, nbas, env, opt) &
            & bind(C, name="cint2e_sph") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt
         integer(c_int) :: stat
      end function cint2e_sph
   end interface

contains

function libcint_cgto_cart(shell, bas) result(nao)
   integer, intent(in) :: shell
   integer(c_int), contiguous, intent(in) :: bas(:, :)
   integer :: nao

   nao = int(cint_cgto_cart(int(shell, c_int), bas))
end function libcint_cgto_cart

function libcint_cgto_spheric(shell, bas) result(nao)
   integer, intent(in) :: shell
   integer(c_int), contiguous, intent(in) :: bas(:, :)
   integer :: nao

   nao = int(cint_cgto_spheric(int(shell, c_int), bas))
end function libcint_cgto_spheric

function libcint_tot_cgto_cart(bas) result(nao)
   integer(c_int), contiguous, intent(in) :: bas(:, :)
   integer :: nao

   nao = int(cint_tot_cgto_cart(bas, int(size(bas, 2), c_int)))
end function libcint_tot_cgto_cart

function libcint_tot_cgto_spheric(bas) result(nao)
   integer(c_int), contiguous, intent(in) :: bas(:, :)
   integer :: nao

   nao = int(cint_tot_cgto_spheric(bas, int(size(bas, 2), c_int)))
end function libcint_tot_cgto_spheric

function libcint_gto_norm(ang_mom, exponent) result(norm)
   integer, intent(in) :: ang_mom
   real(wp), intent(in) :: exponent
   real(wp) :: norm

   norm = real(cint_gto_norm(int(ang_mom, c_int), real(exponent, c_double)), wp)
end function libcint_gto_norm

function libcint_shell_size(shell, bas, representation) result(nao)
   integer, intent(in) :: shell
   integer(c_int), contiguous, intent(in) :: bas(:, :)
   integer, intent(in) :: representation
   integer :: nao

   select case(representation)
   case(LIBCINT_CARTESIAN)
      nao = libcint_cgto_cart(shell, bas)
   case(LIBCINT_SPHERICAL)
      nao = libcint_cgto_spheric(shell, bas)
   case default
      nao = -1
   end select
end function libcint_shell_size

function libcint_eval_1e(kind, representation, out, shls, atm, bas, env) result(stat)
   integer, intent(in) :: kind
   integer, intent(in) :: representation
   real(c_double), contiguous, intent(out) :: out(:, :)
   integer, intent(in) :: shls(2)
   integer(c_int), contiguous, intent(in) :: atm(:, :)
   integer(c_int), contiguous, intent(in) :: bas(:, :)
   real(c_double), contiguous, intent(in) :: env(:)
   integer :: stat

   integer(c_int) :: cshls(2)
   integer :: di, dj

   cshls = int(shls, c_int)
   di = libcint_shell_size(shls(1), bas, representation)
   dj = libcint_shell_size(shls(2), bas, representation)
   if (di < 0 .or. dj < 0 .or. size(out, 1) < di .or. size(out, 2) < dj) then
      stat = -1
      return
   end if

   out(:, :) = 0.0_c_double
   select case(representation)
   case(LIBCINT_CARTESIAN)
      select case(kind)
      case(LIBCINT_1E_OVERLAP)
         stat = int(int1e_ovlp_cart(out, c_null_ptr, cshls, atm, int(size(atm, 2), c_int), &
            & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
      case(LIBCINT_1E_KINETIC)
         stat = int(int1e_kin_cart(out, c_null_ptr, cshls, atm, int(size(atm, 2), c_int), &
            & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
      case(LIBCINT_1E_NUCLEAR)
         stat = int(int1e_nuc_cart(out, c_null_ptr, cshls, atm, int(size(atm, 2), c_int), &
            & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
      case default
         stat = -2
      end select
   case(LIBCINT_SPHERICAL)
      select case(kind)
      case(LIBCINT_1E_OVERLAP)
         stat = int(int1e_ovlp_sph(out, c_null_ptr, cshls, atm, int(size(atm, 2), c_int), &
            & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
      case(LIBCINT_1E_KINETIC)
         stat = int(int1e_kin_sph(out, c_null_ptr, cshls, atm, int(size(atm, 2), c_int), &
            & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
      case(LIBCINT_1E_NUCLEAR)
         stat = int(int1e_nuc_sph(out, c_null_ptr, cshls, atm, int(size(atm, 2), c_int), &
            & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
      case default
         stat = -2
      end select
   case default
      stat = -2
   end select
end function libcint_eval_1e

function libcint_eval_eri(representation, out, shls, atm, bas, env) result(stat)
   integer, intent(in) :: representation
   real(c_double), contiguous, intent(out) :: out(:, :, :, :)
   integer, intent(in) :: shls(4)
   integer(c_int), contiguous, intent(in) :: atm(:, :)
   integer(c_int), contiguous, intent(in) :: bas(:, :)
   real(c_double), contiguous, intent(in) :: env(:)
   integer :: stat

   integer(c_int) :: cshls(4)
   integer :: di, dj, dk, dl

   cshls = int(shls, c_int)
   di = libcint_shell_size(shls(1), bas, representation)
   dj = libcint_shell_size(shls(2), bas, representation)
   dk = libcint_shell_size(shls(3), bas, representation)
   dl = libcint_shell_size(shls(4), bas, representation)
   if (di < 0 .or. dj < 0 .or. dk < 0 .or. dl < 0 .or. &
      & size(out, 1) < di .or. size(out, 2) < dj .or. &
      & size(out, 3) < dk .or. size(out, 4) < dl) then
      stat = -1
      return
   end if

   out(:, :, :, :) = 0.0_c_double
   select case(representation)
   case(LIBCINT_CARTESIAN)
      stat = int(cint2e_cart(out, cshls, atm, int(size(atm, 2), c_int), &
         & bas, int(size(bas, 2), c_int), env, c_null_ptr))
   case(LIBCINT_SPHERICAL)
      stat = int(cint2e_sph(out, cshls, atm, int(size(atm, 2), c_int), &
         & bas, int(size(bas, 2), c_int), env, c_null_ptr))
   case default
      stat = -2
   end select
end function libcint_eval_eri

end module tblite_integral_libcint
