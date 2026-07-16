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
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use tblite_basis_type, only : basis_type
   use tblite_integral_type, only : integral_type
   implicit none
   private

   public :: LIBCINT_1E_OVERLAP, LIBCINT_1E_KINETIC, LIBCINT_1E_NUCLEAR
   public :: LIBCINT_CARTESIAN, LIBCINT_SPHERICAL
   public :: CHARGE_OF, PTR_COORD, NUC_MOD_OF, PTR_ZETA, PTR_FRAC_CHARGE, ATM_SLOTS
   public :: ATOM_OF, ANG_OF, NPRIM_OF, NCTR_OF, KAPPA_OF, PTR_EXP, PTR_COEFF, BAS_SLOTS
   public :: PTR_GRIDS, PTR_ENV_START
   public :: POINT_NUC, GAUSSIAN_NUC, FRAC_CHARGE_NUC
   public :: libcint_cgto_cart, libcint_cgto_spheric
   public :: libcint_tot_cgto_cart, libcint_tot_cgto_spheric
   public :: libcint_gto_norm
   public :: libcint_shell_size
   public :: libcint_eval_1e
   public :: libcint_eval_dipole, libcint_eval_quadrupole
   public :: libcint_eval_dipole_gradient, libcint_eval_quadrupole_gradient
   public :: libcint_eval_overlap_gradient
   public :: libcint_eval_1e_grids
   public :: libcint_eval_eri
   public :: libcint_eval_3c2e, libcint_eval_3c1e_rinv
   public :: libcint_basis_type, libcint_integral_type, new_libcint_basis

   !> Libcint representation of a molecular Gaussian basis.  The integer
   !> tables contain zero-based offsets into env, as required by libcint.
   type :: libcint_basis_type
      integer(c_int), allocatable :: atm(:, :)
      integer(c_int), allocatable :: bas(:, :)
      real(c_double), allocatable :: env(:)
   end type libcint_basis_type

   !> Libcint-backed integral evaluator and cached basis conversion.
   type, extends(integral_type) :: libcint_integral_type
      type(libcint_basis_type) :: basis
   contains
      procedure :: initialize_integral => initialize_libcint
      procedure :: multipole_integral => multipole_libcint
      procedure :: multipole_gradient_integral => multipole_gradient_libcint
      procedure :: dipole_integral => dipole_libcint
   end type libcint_integral_type

   ! libcint C arrays are flattened as slot + slots * item.  These Fortran
   ! constants are shifted by one so normal arrays can be declared as
   ! atm(ATM_SLOTS,natm), bas(BAS_SLOTS,nbas).
   integer, parameter :: CHARGE_OF = 1
   integer, parameter :: PTR_COORD = 2
   integer, parameter :: NUC_MOD_OF = 3
   integer, parameter :: PTR_ZETA = 4
   integer, parameter :: PTR_FRAC_CHARGE = 5
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
   integer, parameter :: PTR_GRIDS = 12
   integer, parameter :: PTR_ENV_START = 20

   integer, parameter :: POINT_NUC = 1
   integer, parameter :: GAUSSIAN_NUC = 2
   integer, parameter :: FRAC_CHARGE_NUC = 3

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
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_ovlp_cart

      function int1e_ovlp_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_ovlp_sph") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_ovlp_sph

      function int1e_kin_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_kin_cart") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_kin_cart

      function int1e_kin_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_kin_sph") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_kin_sph

      function int1e_nuc_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_nuc_cart") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_nuc_cart

      function int1e_nuc_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_nuc_sph") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_nuc_sph

      function int1e_grids_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_grids_sph") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_grids_sph

      function int1e_r_origj_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_r_origj_sph") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_r_origj_sph

      function int1e_rr_origj_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_rr_origj_sph") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_rr_origj_sph

      function int1e_ovlpip_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_ovlpip_sph") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_ovlpip_sph

      function int1e_ipr_origj_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_ipr_origj_sph") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_ipr_origj_sph

      function int1e_r_origj_ip_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_r_origj_ip_sph") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_r_origj_ip_sph

      function int1e_iprr_origj_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_iprr_origj_sph") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_iprr_origj_sph

      function int1e_rr_origj_ip_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_rr_origj_ip_sph") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int1e_rr_origj_ip_sph

      function int3c2e_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int3c2e_cart") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int3c2e_cart

      function int3c2e_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int3c2e_sph") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int3c2e_sph

      function int3c1e_rinv_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int3c1e_rinv_cart") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int3c1e_rinv_cart

      function int3c1e_rinv_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int3c1e_rinv_sph") result(stat)
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*), shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: stat
      end function int3c1e_rinv_sph

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

subroutine initialize_libcint(self, mol, basis)
   class(libcint_integral_type), intent(inout) :: self
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: basis

   call new_libcint_basis(self%basis, mol, basis)
end subroutine initialize_libcint

subroutine dipole_libcint(self, mol, basis, jsh, ish, r2, vec, overlap, dipole)
   class(libcint_integral_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: basis
   integer, intent(in) :: jsh, ish
   real(wp), intent(in) :: r2, vec(3)
   real(wp), intent(out) :: overlap(:), dipole(:, :)
   real(c_double), allocatable :: cdp(:, :, :), cov(:, :)
   integer :: dj, di, jao, iao, ij, stat

   dj = basis%nao_sh(jsh); di = basis%nao_sh(ish)
   allocate(cdp(dj, di, 3), cov(dj, di))
   stat = libcint_eval_1e(LIBCINT_1E_OVERLAP, LIBCINT_SPHERICAL, cov, &
      & [jsh-1, ish-1], self%basis%atm, self%basis%bas, self%basis%env)
   stat = libcint_eval_dipole(cdp, [jsh-1, ish-1], self%basis%atm, &
      & self%basis%bas, self%basis%env)
   do iao = 1, di
      do jao = 1, dj
         ij = jao + dj*(iao-1)
         overlap(ij) = real(cov(jao, iao), wp)
         dipole(:, ij) = real(cdp(jao, iao, :), wp)
      end do
   end do
end subroutine dipole_libcint

subroutine multipole_libcint(self, mol, basis, jsh, ish, r2, vec, overlap, &
      & dipole, quadrupole)
   class(libcint_integral_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: basis
   integer, intent(in) :: jsh, ish
   real(wp), intent(in) :: r2, vec(3)
   real(wp), intent(out) :: overlap(:), dipole(:, :), quadrupole(:, :)
   real(c_double), allocatable :: cqp(:, :, :)
   real(wp) :: raw(6), trace
   integer, parameter :: qmap(6) = [1, 2, 5, 3, 6, 9]
   integer :: dj, di, jao, iao, ij, stat

   call self%dipole_integral(mol, basis, jsh, ish, r2, vec, overlap, dipole)
   dj = basis%nao_sh(jsh); di = basis%nao_sh(ish)
   allocate(cqp(dj, di, 9))
   stat = libcint_eval_quadrupole(cqp, [jsh-1, ish-1], self%basis%atm, &
      & self%basis%bas, self%basis%env)
   do iao = 1, di
      do jao = 1, dj
         ij = jao + dj*(iao-1)
         raw = real(cqp(jao, iao, qmap), wp)
         trace = 0.5_wp*(raw(1) + raw(3) + raw(6))
         quadrupole(:, ij) = 1.5_wp*raw
         quadrupole([1, 3, 6], ij) = quadrupole([1, 3, 6], ij) - trace
      end do
   end do
end subroutine multipole_libcint

subroutine multipole_gradient_libcint(self, mol, basis, jsh, ish, r2, vec, &
      & overlap, dipole, quadrupole, doverlap, ddipole_j, dquadrupole_j, &
      & ddipole_i, dquadrupole_i)
   class(libcint_integral_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: basis
   integer, intent(in) :: jsh, ish
   real(wp), intent(in) :: r2, vec(3)
   real(wp), intent(out) :: overlap(:), dipole(:, :), quadrupole(:, :)
   real(wp), intent(out) :: doverlap(:, :)
   real(wp), intent(out) :: ddipole_j(:, :, :), dquadrupole_j(:, :, :)
   real(wp), intent(out) :: ddipole_i(:, :, :), dquadrupole_i(:, :, :)
   real(c_double), allocatable :: covg(:, :, :)
   real(c_double), allocatable :: cdj(:, :, :, :), cdi(:, :, :, :)
   real(c_double), allocatable :: cqj(:, :, :, :), cqi(:, :, :, :)
   real(c_double), allocatable :: sdj(:, :, :, :), sdi(:, :, :, :)
   real(c_double), allocatable :: sqj(:, :, :, :), sqi(:, :, :, :)
   real(wp) :: raw(6), trace
   integer, parameter :: qmap(6) = [1, 2, 5, 3, 6, 9]
   integer :: dj, di, jao, iao, ij, ic, ider, stat

   call self%multipole_integral(mol, basis, jsh, ish, r2, vec, overlap, dipole, quadrupole)
   dj = basis%nao_sh(jsh); di = basis%nao_sh(ish)
   allocate(covg(dj, di, 3), cdj(dj, di, 3, 3), cdi(dj, di, 3, 3), &
      & cqj(dj, di, 9, 3), cqi(dj, di, 9, 3), &
      & sdj(di, dj, 3, 3), sdi(di, dj, 3, 3), &
      & sqj(di, dj, 9, 3), sqi(di, dj, 9, 3))
   stat = libcint_eval_overlap_gradient(covg, [jsh-1, ish-1], self%basis%atm, &
      & self%basis%bas, self%basis%env)
   stat = libcint_eval_dipole_gradient(cdj, cdi, [jsh-1, ish-1], self%basis%atm, &
      & self%basis%bas, self%basis%env)
   stat = libcint_eval_quadrupole_gradient(cqj, cqi, [jsh-1, ish-1], &
      & self%basis%atm, self%basis%bas, self%basis%env)
   stat = libcint_eval_dipole_gradient(sdi, sdj, [ish-1, jsh-1], self%basis%atm, &
      & self%basis%bas, self%basis%env)
   stat = libcint_eval_quadrupole_gradient(sqi, sqj, [ish-1, jsh-1], &
      & self%basis%atm, self%basis%bas, self%basis%env)
   do iao = 1, di
      do jao = 1, dj
         ij = jao + dj*(iao-1)
         do ider = 1, 3
            doverlap(ider, ij) = real(covg(jao, iao, ider), wp)
            ddipole_i(ider, :, ij) = real(cdi(jao, iao, :, ider), wp)
            ddipole_j(ider, :, ij) = real(sdi(iao, jao, :, ider), wp)
            do ic = 1, 6
               raw(ic) = real(cqi(jao, iao, qmap(ic), ider), wp)
            end do
            trace = 0.5_wp*(raw(1) + raw(3) + raw(6))
            dquadrupole_i(ider, :, ij) = 1.5_wp*raw
            dquadrupole_i(ider, [1, 3, 6], ij) = &
               & dquadrupole_i(ider, [1, 3, 6], ij) - trace
            do ic = 1, 6
               raw(ic) = real(sqi(iao, jao, qmap(ic), ider), wp)
            end do
            trace = 0.5_wp*(raw(1) + raw(3) + raw(6))
            dquadrupole_j(ider, :, ij) = 1.5_wp*raw
            dquadrupole_j(ider, [1, 3, 6], ij) = &
               & dquadrupole_j(ider, [1, 3, 6], ij) - trace
         end do
      end do
   end do
end subroutine multipole_gradient_libcint

!> Pack a tblite structure and basis into libcint's atm/bas/env format.
subroutine new_libcint_basis(self, mol, basis)
   type(libcint_basis_type), intent(out) :: self
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: basis

   integer :: iat, isp, ish, lsh, ip, off, nenv, lang
   real(wp) :: spherical_norm

   nenv = PTR_ENV_START + 3*mol%nat
   do ish = 1, basis%nsh
      iat = basis%sh2at(ish)
      isp = mol%id(iat)
      lsh = ish - basis%ish_at(iat)
      nenv = nenv + 2*basis%cgto(lsh, isp)%nprim
   end do

   allocate(self%atm(ATM_SLOTS, mol%nat), source=0_c_int)
   allocate(self%bas(BAS_SLOTS, basis%nsh), source=0_c_int)
   allocate(self%env(nenv), source=0.0_c_double)

   off = PTR_ENV_START
   do iat = 1, mol%nat
      isp = mol%id(iat)
      self%atm(CHARGE_OF, iat) = int(mol%num(isp), c_int)
      self%atm(PTR_COORD, iat) = int(off, c_int)
      self%atm(NUC_MOD_OF, iat) = POINT_NUC
      self%env(off+1:off+3) = real(mol%xyz(:, iat), c_double)
      off = off + 3
   end do

   do ish = 1, basis%nsh
      iat = basis%sh2at(ish)
      isp = mol%id(iat)
      lsh = ish - basis%ish_at(iat)

      self%bas(ATOM_OF, ish) = int(iat-1, c_int)
      self%bas(ANG_OF, ish) = int(basis%cgto(lsh, isp)%ang, c_int)
      self%bas(NPRIM_OF, ish) = int(basis%cgto(lsh, isp)%nprim, c_int)
      self%bas(NCTR_OF, ish) = 1_c_int
      self%bas(KAPPA_OF, ish) = 0_c_int

      self%bas(PTR_EXP, ish) = int(off, c_int)
      do ip = 1, basis%cgto(lsh, isp)%nprim
         self%env(off+ip) = real(basis%cgto(lsh, isp)%alpha(ip), c_double)
      end do
      off = off + basis%cgto(lsh, isp)%nprim

      self%bas(PTR_COEFF, ish) = int(off, c_int)
      lang = basis%cgto(lsh, isp)%ang
      ! Tblite uses unnormalised real solid harmonics, whereas libcint's
      ! spherical functions contain a normalized angular factor.
      spherical_norm = sqrt(4.0_wp*pi/real(2*lang+1, wp))
      do ip = 1, basis%cgto(lsh, isp)%nprim
         self%env(off+ip) = real(spherical_norm*basis%cgto(lsh, isp)%coeff(ip), c_double)
      end do
      off = off + basis%cgto(lsh, isp)%nprim

   end do
end subroutine new_libcint_basis

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

   integer(c_int) :: cshls(2), dims(2)
   integer :: di, dj

   cshls = int(shls, c_int) ! convert to C-int for libcint
   dims = int([size(out, 1), size(out, 2)], c_int)
   di = libcint_shell_size(shls(1), bas, representation) ! determine output dimensions
   dj = libcint_shell_size(shls(2), bas, representation) ! (how many spherical basis functions in each shell)
   if (di < 0 .or. dj < 0 .or. size(out, 1) < di .or. size(out, 2) < dj) then
      stat = -1
      return
   end if

   out(:, :) = 0.0_c_double
   select case(representation)
   case(LIBCINT_CARTESIAN)
      select case(kind)
      case(LIBCINT_1E_OVERLAP)
         stat = int(int1e_ovlp_cart(out, dims, cshls, atm, int(size(atm, 2), c_int), &
            & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
      case(LIBCINT_1E_KINETIC)
         stat = int(int1e_kin_cart(out, dims, cshls, atm, int(size(atm, 2), c_int), &
            & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
      case(LIBCINT_1E_NUCLEAR)
         stat = int(int1e_nuc_cart(out, dims, cshls, atm, int(size(atm, 2), c_int), &
            & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
      case default
         stat = -2
      end select
   case(LIBCINT_SPHERICAL)
      select case(kind)
      case(LIBCINT_1E_OVERLAP)
         stat = int(int1e_ovlp_sph(out, dims, cshls, atm, int(size(atm, 2), c_int), &
            & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
      case(LIBCINT_1E_KINETIC)
         stat = int(int1e_kin_sph(out, dims, cshls, atm, int(size(atm, 2), c_int), &
            & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
      case(LIBCINT_1E_NUCLEAR)
         stat = int(int1e_nuc_sph(out, dims, cshls, atm, int(size(atm, 2), c_int), &
            & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
      case default
         stat = -2
      end select
   case default
      stat = -2
   end select
end function libcint_eval_1e

function libcint_eval_dipole(out, shls, atm, bas, env) result(stat)
   real(c_double), contiguous, intent(out) :: out(:, :, :)
   integer, intent(in) :: shls(2)
   integer(c_int), contiguous, intent(in) :: atm(:, :), bas(:, :)
   real(c_double), contiguous, intent(in) :: env(:)
   integer :: stat
   integer(c_int) :: cshls(2), dims(2)
   integer :: di, dj

   di = libcint_shell_size(shls(1), bas, LIBCINT_SPHERICAL)
   dj = libcint_shell_size(shls(2), bas, LIBCINT_SPHERICAL)
   if (size(out, 1) < di .or. size(out, 2) < dj .or. size(out, 3) < 3) then
      stat = -1
      return
   end if
   cshls = int(shls, c_int)
   dims = int([size(out, 1), size(out, 2)], c_int)
   out = 0.0_c_double
   stat = int(int1e_r_origj_sph(out, dims, cshls, atm, int(size(atm, 2), c_int), &
      & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
end function libcint_eval_dipole

function libcint_eval_quadrupole(out, shls, atm, bas, env) result(stat)
   real(c_double), contiguous, intent(out) :: out(:, :, :)
   integer, intent(in) :: shls(2)
   integer(c_int), contiguous, intent(in) :: atm(:, :), bas(:, :)
   real(c_double), contiguous, intent(in) :: env(:)
   integer :: stat
   integer(c_int) :: cshls(2), dims(2)
   integer :: di, dj

   di = libcint_shell_size(shls(1), bas, LIBCINT_SPHERICAL)
   dj = libcint_shell_size(shls(2), bas, LIBCINT_SPHERICAL)
   if (size(out, 1) < di .or. size(out, 2) < dj .or. size(out, 3) < 9) then
      stat = -1
      return
   end if
   cshls = int(shls, c_int)
   dims = int([size(out, 1), size(out, 2)], c_int)
   out = 0.0_c_double
   stat = int(int1e_rr_origj_sph(out, dims, cshls, atm, int(size(atm, 2), c_int), &
      & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
end function libcint_eval_quadrupole

function libcint_eval_dipole_gradient(out_bra, out_ket, shls, atm, bas, env) result(stat)
   !> Complete nuclear-center derivatives of a ket-centered dipole integral.
   real(c_double), contiguous, intent(out) :: out_bra(:, :, :, :), out_ket(:, :, :, :)
   integer, intent(in) :: shls(2)
   integer(c_int), contiguous, intent(in) :: atm(:, :), bas(:, :)
   real(c_double), contiguous, intent(in) :: env(:)
   integer :: stat
   integer(c_int) :: cshls(2), dims(2)
   integer :: di, dj, stat_bra

   di = libcint_shell_size(shls(1), bas, LIBCINT_SPHERICAL)
   dj = libcint_shell_size(shls(2), bas, LIBCINT_SPHERICAL)
   if (size(out_bra, 1) < di .or. size(out_bra, 2) < dj .or. &
      & size(out_bra, 3) < 3 .or. size(out_bra, 4) < 3 .or. &
      & any(shape(out_ket) < shape(out_bra))) then
      stat = -1
      return
   end if
   cshls = int(shls, c_int)
   dims = int([size(out_bra, 1), size(out_bra, 2)], c_int)
   out_bra = 0.0_c_double
   out_ket = 0.0_c_double
   stat_bra = int(int1e_ipr_origj_sph(out_bra, dims, cshls, atm, &
      & int(size(atm, 2), c_int), bas, int(size(bas, 2), c_int), env, &
      & c_null_ptr, c_null_ptr))
   stat = int(int1e_r_origj_ip_sph(out_ket, dims, cshls, atm, &
      & int(size(atm, 2), c_int), bas, int(size(bas, 2), c_int), env, &
      & c_null_ptr, c_null_ptr))
   if (stat_bra < 0) stat = stat_bra
end function libcint_eval_dipole_gradient

function libcint_eval_quadrupole_gradient(out_bra, out_ket, shls, atm, bas, env) result(stat)
   !> Complete nuclear-center derivatives of a ket-centered Cartesian second moment.
   real(c_double), contiguous, intent(out) :: out_bra(:, :, :, :), out_ket(:, :, :, :)
   integer, intent(in) :: shls(2)
   integer(c_int), contiguous, intent(in) :: atm(:, :), bas(:, :)
   real(c_double), contiguous, intent(in) :: env(:)
   integer :: stat
   integer(c_int) :: cshls(2), dims(2)
   integer :: di, dj, stat_bra

   di = libcint_shell_size(shls(1), bas, LIBCINT_SPHERICAL)
   dj = libcint_shell_size(shls(2), bas, LIBCINT_SPHERICAL)
   if (size(out_bra, 1) < di .or. size(out_bra, 2) < dj .or. &
      & size(out_bra, 3) < 9 .or. size(out_bra, 4) < 3 .or. &
      & any(shape(out_ket) < shape(out_bra))) then
      stat = -1
      return
   end if
   cshls = int(shls, c_int)
   dims = int([size(out_bra, 1), size(out_bra, 2)], c_int)
   out_bra = 0.0_c_double
   out_ket = 0.0_c_double
   stat_bra = int(int1e_iprr_origj_sph(out_bra, dims, cshls, atm, &
      & int(size(atm, 2), c_int), bas, int(size(bas, 2), c_int), env, &
      & c_null_ptr, c_null_ptr))
   stat = int(int1e_rr_origj_ip_sph(out_ket, dims, cshls, atm, &
      & int(size(atm, 2), c_int), bas, int(size(bas, 2), c_int), env, &
      & c_null_ptr, c_null_ptr))
   if (stat_bra < 0) stat = stat_bra
end function libcint_eval_quadrupole_gradient

function libcint_eval_overlap_gradient(out, shls, atm, bas, env) result(stat)
   !> Overlap derivative with respect to the nuclear centre of shls(2).
   real(c_double), contiguous, intent(out) :: out(:, :, :)
   integer, intent(in) :: shls(2)
   integer(c_int), contiguous, intent(in) :: atm(:, :), bas(:, :)
   real(c_double), contiguous, intent(in) :: env(:)
   integer :: stat
   integer(c_int) :: cshls(2), dims(2)
   integer :: di, dj

   di = libcint_shell_size(shls(1), bas, LIBCINT_SPHERICAL)
   dj = libcint_shell_size(shls(2), bas, LIBCINT_SPHERICAL)
   if (size(out, 1) < di .or. size(out, 2) < dj .or. size(out, 3) < 3) then
      stat = -1
      return
   end if
   cshls = int(shls, c_int)
   dims = int([size(out, 1), size(out, 2)], c_int)
   out = 0.0_c_double
   stat = int(int1e_ovlpip_sph(out, dims, cshls, atm, int(size(atm, 2), c_int), &
      & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
   ! Libcint differentiates the electronic coordinate of the basis function;
   ! tblite differentiates its nuclear centre, which has the opposite sign.
   if (stat >= 0) out = -out
end function libcint_eval_overlap_gradient

function libcint_eval_1e_grids(out, shls, grid_range, atm, bas, env) result(stat)
   real(c_double), contiguous, intent(out) :: out(:, :, :)
   integer, intent(in) :: shls(2)
   integer, intent(in) :: grid_range(2)
   integer(c_int), contiguous, intent(in) :: atm(:, :)
   integer(c_int), contiguous, intent(in) :: bas(:, :)
   real(c_double), contiguous, intent(in) :: env(:)
   integer :: stat

   integer(c_int) :: cshls(4), dims(3)
   integer :: di, dj, ngrids

   di = libcint_shell_size(shls(1), bas, LIBCINT_SPHERICAL)
   dj = libcint_shell_size(shls(2), bas, LIBCINT_SPHERICAL)
   ngrids = grid_range(2) - grid_range(1)
   if (di < 0 .or. dj < 0 .or. ngrids < 0 .or. &
      & size(out, 1) < ngrids .or. size(out, 2) < di .or. size(out, 3) < dj) then
      stat = -1
      return
   end if

   cshls = int([shls(1), shls(2), grid_range(1), grid_range(2)], c_int)
   dims = int([di, dj, max(ngrids, 1)], c_int)
   out(:, :, :) = 0.0_c_double
   stat = int(int1e_grids_sph(out, dims, cshls, atm, int(size(atm, 2), c_int), &
      & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
end function libcint_eval_1e_grids

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

function libcint_eval_3c2e(representation, out, shls, atm, bas, env) result(stat)
   integer, intent(in) :: representation
   real(c_double), contiguous, intent(out) :: out(:, :, :)
   integer, intent(in) :: shls(3)
   integer(c_int), contiguous, intent(in) :: atm(:, :), bas(:, :)
   real(c_double), contiguous, intent(in) :: env(:)
   integer :: stat
   integer(c_int) :: cshls(3), dims(3)
   integer :: i, shell_dims(3)

   do i = 1, 3
      shell_dims(i) = libcint_shell_size(shls(i), bas, representation)
   end do
   if (any(shell_dims < 0) .or. any(shape(out) < shell_dims)) then
      stat = -1
      return
   end if

   cshls = int(shls, c_int)
   dims = int(shape(out), c_int)
   out = 0.0_c_double
   select case(representation)
   case(LIBCINT_CARTESIAN)
      stat = int(int3c2e_cart(out, dims, cshls, atm, int(size(atm, 2), c_int), &
         & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
   case(LIBCINT_SPHERICAL)
      stat = int(int3c2e_sph(out, dims, cshls, atm, int(size(atm, 2), c_int), &
         & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
   case default
      stat = -2
   end select
end function libcint_eval_3c2e

function libcint_eval_3c1e_rinv(representation, out, shls, atm, bas, env) result(stat)
   integer, intent(in) :: representation
   real(c_double), contiguous, intent(out) :: out(:, :, :)
   integer, intent(in) :: shls(3)
   integer(c_int), contiguous, intent(in) :: atm(:, :), bas(:, :)
   real(c_double), contiguous, intent(in) :: env(:)
   integer :: stat
   integer(c_int) :: cshls(3), dims(3)
   integer :: i, shell_dims(3)

   do i = 1, 3
      shell_dims(i) = libcint_shell_size(shls(i), bas, representation)
   end do
   if (any(shell_dims < 0) .or. any(shape(out) < shell_dims)) then
      stat = -1
      return
   end if

   cshls = int(shls, c_int)
   dims = int(shape(out), c_int)
   out = 0.0_c_double
   select case(representation)
   case(LIBCINT_CARTESIAN)
      stat = int(int3c1e_rinv_cart(out, dims, cshls, atm, int(size(atm, 2), c_int), &
         & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
   case(LIBCINT_SPHERICAL)
      stat = int(int3c1e_rinv_sph(out, dims, cshls, atm, int(size(atm, 2), c_int), &
         & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
   case default
      stat = -2
   end select
end function libcint_eval_3c1e_rinv

end module tblite_integral_libcint
