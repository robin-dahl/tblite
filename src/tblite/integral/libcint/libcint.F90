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

!> @file tblite/integral/libcint/libcint.f90
!> Provides an optional libcint-backed Gaussian integral interface.

!> Interface to Libcint Gaussian integral evaluation.
module tblite_integral_libcint
   use, intrinsic :: iso_c_binding, only : c_double, c_int, c_null_ptr, c_ptr, c_size_t
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use tblite_basis_type, only : basis_type, cgto_type
   use tblite_integral_handler, only : integral_handler
   use tblite_integral_shell, only : msao
   implicit none
   private

   public :: CHARGE_OF, PTR_COORD, NUC_MOD_OF, PTR_ZETA, PTR_FRAC_CHARGE, ATM_SLOTS
   public :: ATOM_OF, ANG_OF, NPRIM_OF, NCTR_OF, KAPPA_OF, PTR_EXP, PTR_COEFF, BAS_SLOTS
   public :: PTR_GRIDS, PTR_ENV_START
   public :: POINT_NUC, GAUSSIAN_NUC, FRAC_CHARGE_NUC
   public :: libcint_eval_1e
   public :: libcint_eval_dipole, libcint_eval_quadrupole
   public :: libcint_eval_dipole_gradient, libcint_eval_quadrupole_gradient
   public :: libcint_eval_overlap_gradient
   public :: libcint_eval_3c2e
   public :: libcint_basis_type, libcint_integral_type

   !> Libcint representation of a molecular Gaussian basis.  The integer
   !> tables contain zero-based offsets into env, as required by libcint
   type :: libcint_basis_type
      !> Libcint atom table
      integer(c_int), allocatable :: atm(:, :)
      !> Libcint basis-shell table
      integer(c_int), allocatable :: bas(:, :)
      !> Libcint floating-point environment array
      real(c_double), allocatable :: env(:)
   end type libcint_basis_type

   !> Libcint-backed integral evaluator and cached basis conversion
   type, extends(integral_handler) :: libcint_integral_type
      !> Molecular basis represented in libcint data structures
      type(libcint_basis_type) :: basis
   contains
      !> Convert and cache the molecular basis for libcint evaluations
      procedure :: initialize_integral => initialize_libcint
      !> Evaluate overlap, dipole, and quadrupole integrals
      procedure :: multipole_cgto => multipole_libcint
      !> Evaluate multipole integrals and their nuclear derivatives
      procedure :: multipole_grad_cgto => multipole_grad_libcint
      !> Evaluate overlap and dipole integrals
      procedure :: dipole_cgto => dipole_libcint
      !> Evaluate AO-pair Coulomb integrals with normalized Gaussian charges
      procedure :: surface_3c2e => surface_3c2e_libcint
   end type libcint_integral_type

   ! libcint C arrays are flattened as slot + slots * item.  These Fortran
   ! constants are shifted by one so normal arrays can be declared as
   ! atm(ATM_SLOTS,natm), bas(BAS_SLOTS,nbas)
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
   ! code should write env(offset+1:offset+n) for C locations offset:offset+n-1
   integer, parameter :: PTR_GRIDS = 12
   integer, parameter :: PTR_ENV_START = 20

   integer, parameter :: POINT_NUC = 1
   integer, parameter :: GAUSSIAN_NUC = 2
   integer, parameter :: FRAC_CHARGE_NUC = 3

   interface
      !> Return the number of spherical atomic orbitals in a libcint shell
      function cint_cgto_spheric(bas_id, bas) bind(C, name="CINTcgto_spheric") result(nao)
         import :: c_int
         !> Zero-based index of the shell
         integer(c_int), value :: bas_id
         !> Flattened libcint basis-shell table
         integer(c_int), intent(in) :: bas(*)
         !> Number of spherical atomic orbitals in the shell
         integer(c_int) :: nao
      end function cint_cgto_spheric

      !> Evaluate <i| OVLP |j>
      function int1e_ovlp_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_ovlp_sph") result(stat)
         import :: c_double, c_int, c_ptr
         !> Flattened output buffer
         real(c_double), intent(out) :: out(*)
         !> Leading dimensions of the output buffer
         integer(c_int), intent(in) :: dims(*)
         !> Zero-based indices of the bra and ket shells
         integer(c_int), intent(in) :: shls(*)
         !> Flattened libcint atom table
         integer(c_int), intent(in) :: atm(*)
         !> Number of atoms in the atom table
         integer(c_int), value :: natm
         !> Flattened libcint basis-shell table
         integer(c_int), intent(in) :: bas(*)
         !> Number of shells in the basis-shell table
         integer(c_int), value :: nbas
         !> Libcint floating-point environment array
         real(c_double), intent(in) :: env(*)
         !> Libcint optimizer, or a null pointer
         type(c_ptr), value :: opt
         !> Libcint workspace cache, or a null pointer
         type(c_ptr), value :: cache
         !> Libcint return status
         integer(c_int) :: stat
      end function int1e_ovlp_sph

      !> Evaluate (i j|k) three-center two-electron repulsion integrals
      function int3c2e_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int3c2e_sph") result(stat)
         import :: c_double, c_int, c_ptr, c_size_t
         real(c_double), intent(out) :: out(*)
         integer(c_int), intent(in) :: dims(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value :: nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value :: opt
         type(c_ptr), value :: cache
         integer(c_size_t) :: stat
      end function int3c2e_sph

      !> Evaluate <i| R |j>
      function int1e_r_origj_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_r_origj_sph") result(stat)
         import :: c_double, c_int, c_ptr
         !> Flattened output buffer
         real(c_double), intent(out) :: out(*)
         !> Leading dimensions of the output buffer
         integer(c_int), intent(in) :: dims(*)
         !> Zero-based indices of the bra and ket shells
         integer(c_int), intent(in) :: shls(*)
         !> Flattened libcint atom table
         integer(c_int), intent(in) :: atm(*)
         !> Number of atoms in the atom table
         integer(c_int), value :: natm
         !> Flattened libcint basis-shell table
         integer(c_int), intent(in) :: bas(*)
         !> Number of shells in the basis-shell table
         integer(c_int), value :: nbas
         !> Libcint floating-point environment array
         real(c_double), intent(in) :: env(*)
         !> Libcint optimizer, or a null pointer
         type(c_ptr), value :: opt
         !> Libcint workspace cache, or a null pointer
         type(c_ptr), value :: cache
         !> Libcint return status
         integer(c_int) :: stat
      end function int1e_r_origj_sph

      !> Evaluate <i| R R |j>
      function int1e_rr_origj_sph(out, dims, shls, atm, natm, bas, nbas, &
            & env, opt, cache) &
            & bind(C, name="int1e_rr_origj_sph") result(stat)
         import :: c_double, c_int, c_ptr
         !> Flattened output buffer
         real(c_double), intent(out) :: out(*)
         !> Leading dimensions of the output buffer
         integer(c_int), intent(in) :: dims(*)
         !> Zero-based indices of the bra and ket shells
         integer(c_int), intent(in) :: shls(*)
         !> Flattened libcint atom table
         integer(c_int), intent(in) :: atm(*)
         !> Number of atoms in the atom table
         integer(c_int), value :: natm
         !> Flattened libcint basis-shell table
         integer(c_int), intent(in) :: bas(*)
         !> Number of shells in the basis-shell table
         integer(c_int), value :: nbas
         !> Libcint floating-point environment array
         real(c_double), intent(in) :: env(*)
         !> Libcint optimizer, or a null pointer
         type(c_ptr), value :: opt
         !> Libcint workspace cache, or a null pointer
         type(c_ptr), value :: cache
         !> Libcint return status
         integer(c_int) :: stat
      end function int1e_rr_origj_sph

      !> Evaluate <i| OVLP |NABLA j>
      function int1e_ovlpip_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            & bind(C, name="int1e_ovlpip_sph") result(stat)
         import :: c_double, c_int, c_ptr
         !> Flattened output buffer
         real(c_double), intent(out) :: out(*)
         !> Leading dimensions of the output buffer
         integer(c_int), intent(in) :: dims(*)
         !> Zero-based indices of the bra and ket shells
         integer(c_int), intent(in) :: shls(*)
         !> Flattened libcint atom table
         integer(c_int), intent(in) :: atm(*)
         !> Number of atoms in the atom table
         integer(c_int), value :: natm
         !> Flattened libcint basis-shell table
         integer(c_int), intent(in) :: bas(*)
         !> Number of shells in the basis-shell table
         integer(c_int), value :: nbas
         !> Libcint floating-point environment array
         real(c_double), intent(in) :: env(*)
         !> Libcint optimizer, or a null pointer
         type(c_ptr), value :: opt
         !> Libcint workspace cache, or a null pointer
         type(c_ptr), value :: cache
         !> Libcint return status
         integer(c_int) :: stat
      end function int1e_ovlpip_sph

      !> Evaluate <NABLA i| OVLP |R j>
      function int1e_ipr_origj_sph(out, dims, shls, atm, natm, bas, nbas, &
            & env, opt, cache) &
            & bind(C, name="int1e_ipr_origj_sph") result(stat)
         import :: c_double, c_int, c_ptr
         !> Flattened output buffer
         real(c_double), intent(out) :: out(*)
         !> Leading dimensions of the output buffer
         integer(c_int), intent(in) :: dims(*)
         !> Zero-based indices of the bra and ket shells
         integer(c_int), intent(in) :: shls(*)
         !> Flattened libcint atom table
         integer(c_int), intent(in) :: atm(*)
         !> Number of atoms in the atom table
         integer(c_int), value :: natm
         !> Flattened libcint basis-shell table
         integer(c_int), intent(in) :: bas(*)
         !> Number of shells in the basis-shell table
         integer(c_int), value :: nbas
         !> Libcint floating-point environment array
         real(c_double), intent(in) :: env(*)
         !> Libcint optimizer, or a null pointer
         type(c_ptr), value :: opt
         !> Libcint workspace cache, or a null pointer
         type(c_ptr), value :: cache
         !> Libcint return status
         integer(c_int) :: stat
      end function int1e_ipr_origj_sph

      !> Evaluate <i| OVLP |NABLA R j>
      function int1e_r_origj_ip_sph(out, dims, shls, atm, natm, bas, nbas, &
            & env, opt, cache) &
            & bind(C, name="int1e_r_origj_ip_sph") result(stat)
         import :: c_double, c_int, c_ptr
         !> Flattened output buffer
         real(c_double), intent(out) :: out(*)
         !> Leading dimensions of the output buffer
         integer(c_int), intent(in) :: dims(*)
         !> Zero-based indices of the bra and ket shells
         integer(c_int), intent(in) :: shls(*)
         !> Flattened libcint atom table
         integer(c_int), intent(in) :: atm(*)
         !> Number of atoms in the atom table
         integer(c_int), value :: natm
         !> Flattened libcint basis-shell table
         integer(c_int), intent(in) :: bas(*)
         !> Number of shells in the basis-shell table
         integer(c_int), value :: nbas
         !> Libcint floating-point environment array
         real(c_double), intent(in) :: env(*)
         !> Libcint optimizer, or a null pointer
         type(c_ptr), value :: opt
         !> Libcint workspace cache, or a null pointer
         type(c_ptr), value :: cache
         !> Libcint return status
         integer(c_int) :: stat
      end function int1e_r_origj_ip_sph

      !> Evaluate <NABLA i| OVLP |R R j>
      function int1e_iprr_origj_sph(out, dims, shls, atm, natm, bas, nbas, &
            & env, opt, cache) &
            & bind(C, name="int1e_iprr_origj_sph") result(stat)
         import :: c_double, c_int, c_ptr
         !> Flattened output buffer
         real(c_double), intent(out) :: out(*)
         !> Leading dimensions of the output buffer
         integer(c_int), intent(in) :: dims(*)
         !> Zero-based indices of the bra and ket shells
         integer(c_int), intent(in) :: shls(*)
         !> Flattened libcint atom table
         integer(c_int), intent(in) :: atm(*)
         !> Number of atoms in the atom table
         integer(c_int), value :: natm
         !> Flattened libcint basis-shell table
         integer(c_int), intent(in) :: bas(*)
         !> Number of shells in the basis-shell table
         integer(c_int), value :: nbas
         !> Libcint floating-point environment array
         real(c_double), intent(in) :: env(*)
         !> Libcint optimizer, or a null pointer
         type(c_ptr), value :: opt
         !> Libcint workspace cache, or a null pointer
         type(c_ptr), value :: cache
         !> Libcint return status
         integer(c_int) :: stat
      end function int1e_iprr_origj_sph

      !> Evaluate <i| OVLP |NABLA R R j>
      function int1e_rr_origj_ip_sph(out, dims, shls, atm, natm, bas, nbas, &
            & env, opt, cache) &
            & bind(C, name="int1e_rr_origj_ip_sph") result(stat)
         import :: c_double, c_int, c_ptr
         !> Flattened output buffer
         real(c_double), intent(out) :: out(*)
         !> Leading dimensions of the output buffer
         integer(c_int), intent(in) :: dims(*)
         !> Zero-based indices of the bra and ket shells
         integer(c_int), intent(in) :: shls(*)
         !> Flattened libcint atom table
         integer(c_int), intent(in) :: atm(*)
         !> Number of atoms in the atom table
         integer(c_int), value :: natm
         !> Flattened libcint basis-shell table
         integer(c_int), intent(in) :: bas(*)
         !> Number of shells in the basis-shell table
         integer(c_int), value :: nbas
         !> Libcint floating-point environment array
         real(c_double), intent(in) :: env(*)
         !> Libcint optimizer, or a null pointer
         type(c_ptr), value :: opt
         !> Libcint workspace cache, or a null pointer
         type(c_ptr), value :: cache
         !> Libcint return status
         integer(c_int) :: stat
      end function int1e_rr_origj_ip_sph

   end interface

contains

!> Convert and cache the complete basis for libcint evaluations
subroutine initialize_libcint(self, mol, basis)
   !> Libcint integral evaluator
   class(libcint_integral_type), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
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

   allocate(self%basis%atm(ATM_SLOTS, mol%nat), source=0_c_int)
   allocate(self%basis%bas(BAS_SLOTS, basis%nsh), source=0_c_int)
   allocate(self%basis%env(nenv), source=0.0_c_double)

   off = PTR_ENV_START
   do iat = 1, mol%nat
      isp = mol%id(iat)
      self%basis%atm(CHARGE_OF, iat) = int(mol%num(isp), c_int)
      self%basis%atm(PTR_COORD, iat) = int(off, c_int)
      self%basis%atm(NUC_MOD_OF, iat) = POINT_NUC
      self%basis%env(off+1:off+3) = real(mol%xyz(:, iat), c_double)
      off = off + 3
   end do

   do ish = 1, basis%nsh
      iat = basis%sh2at(ish)
      isp = mol%id(iat)
      lsh = ish - basis%ish_at(iat)

      self%basis%bas(ATOM_OF, ish) = int(iat-1, c_int)
      self%basis%bas(ANG_OF, ish) = int(basis%cgto(lsh, isp)%ang, c_int)
      self%basis%bas(NPRIM_OF, ish) = int(basis%cgto(lsh, isp)%nprim, c_int)
      self%basis%bas(NCTR_OF, ish) = 1_c_int
      self%basis%bas(KAPPA_OF, ish) = 0_c_int

      self%basis%bas(PTR_EXP, ish) = int(off, c_int)
      do ip = 1, basis%cgto(lsh, isp)%nprim
         self%basis%env(off+ip) = real(basis%cgto(lsh, isp)%alpha(ip), c_double)
      end do
      off = off + basis%cgto(lsh, isp)%nprim

      self%basis%bas(PTR_COEFF, ish) = int(off, c_int)
      lang = basis%cgto(lsh, isp)%ang
      ! Tblite uses unnormalised real solid harmonics, whereas libcint's
      ! spherical functions contain a normalized angular factor.
      spherical_norm = sqrt(4.0_wp*pi/real(2*lang+1, wp))
      do ip = 1, basis%cgto(lsh, isp)%nprim
         self%basis%env(off+ip) = real(spherical_norm*basis%cgto(lsh, isp)%coeff(ip), c_double)
      end do
      off = off + basis%cgto(lsh, isp)%nprim
   end do
end subroutine initialize_libcint

!> Evaluate B(i,mu,nu)=(mu nu|g_i) for unit-charge Gaussian surface functions.
subroutine surface_3c2e_libcint(self, mol, basis, xyz, xi, bmat)
   class(libcint_integral_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: basis
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: xi(:)
   real(wp), allocatable, intent(out) :: bmat(:, :, :)

   integer(c_int), allocatable :: atm(:, :), bas(:, :)
   real(c_double), allocatable :: env(:), block(:, :, :)
   integer :: natm, nbas, nenv, old_env, off, igrid, ish, jsh
   integer :: iao, jao, ni, nj, ii, jj, stat
   real(wp) :: alpha, coeff

   natm = size(self%basis%atm, 2) + size(xi)
   nbas = size(self%basis%bas, 2) + size(xi)
   old_env = size(self%basis%env)
   nenv = old_env + 5*size(xi)
   allocate(atm(ATM_SLOTS, natm), source=0_c_int)
   allocate(bas(BAS_SLOTS, nbas), source=0_c_int)
   allocate(env(nenv), source=0.0_c_double)
   atm(:, :mol%nat) = self%basis%atm
   bas(:, :basis%nsh) = self%basis%bas
   env(:old_env) = self%basis%env

   off = old_env
   do igrid = 1, size(xi)
      atm(CHARGE_OF, mol%nat+igrid) = 0_c_int
      atm(PTR_COORD, mol%nat+igrid) = int(off, c_int)
      atm(NUC_MOD_OF, mol%nat+igrid) = POINT_NUC
      env(off+1:off+3) = real(xyz(:, igrid), c_double)
      off = off + 3

      ! Moist defines xi as the square root of the Gaussian exponent:
      ! g(r)=(xi^2/pi)^(3/2)*exp(-xi^2*r^2).
      alpha = xi(igrid)**2
      ! Tblite multiplies every libcint spherical shell coefficient by its
      ! angular normalization conversion.  Apply the same l=0 conversion.
      coeff = sqrt(4.0_wp*pi)*(alpha/pi)**1.5_wp
      bas(ATOM_OF, basis%nsh+igrid) = int(mol%nat+igrid-1, c_int)
      bas(ANG_OF, basis%nsh+igrid) = 0_c_int
      bas(NPRIM_OF, basis%nsh+igrid) = 1_c_int
      bas(NCTR_OF, basis%nsh+igrid) = 1_c_int
      bas(KAPPA_OF, basis%nsh+igrid) = 0_c_int
      bas(PTR_EXP, basis%nsh+igrid) = int(off, c_int)
      env(off+1) = real(alpha, c_double)
      off = off + 1
      bas(PTR_COEFF, basis%nsh+igrid) = int(off, c_int)
      env(off+1) = real(coeff, c_double)
      off = off + 1
   end do

   allocate(bmat(size(xi), basis%nao, basis%nao), source=0.0_wp)
   do igrid = 1, size(xi)
      do ish = 1, basis%nsh
         ni = basis%nao_sh(ish)
         ii = basis%iao_sh(ish)
         do jsh = 1, basis%nsh
            nj = basis%nao_sh(jsh)
            jj = basis%iao_sh(jsh)
            allocate(block(nj, ni, 1), source=0.0_c_double)
            stat = libcint_eval_3c2e(block, &
               & [jsh-1, ish-1, basis%nsh+igrid-1], atm, bas, env)
            if (stat > 0) then
               do iao = 1, ni
                  do jao = 1, nj
                     bmat(igrid, ii+iao, jj+jao) = real(block(jao, iao, 1), wp)
                  end do
               end do
            end if
            deallocate(block)
         end do
      end do
   end do
end subroutine surface_3c2e_libcint

!> Evaluate overlap and dipole integrals for a shell pair
subroutine dipole_libcint(self, cgtoj, cgtoi, jsh, ish, r2, vec, intcut, overlap, dpint)
   !> Libcint integral evaluator
   class(libcint_integral_type), intent(in) :: self
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
   real(c_double), allocatable :: cdp(:, :, :), cov(:, :)
   integer :: dj, di, jao, iao, stat

   ! Query libcint so scratch dimensions agree with its cached shell representation
   dj = int(cint_cgto_spheric(int(jsh - 1, c_int), self%basis%bas))
   di = int(cint_cgto_spheric(int(ish - 1, c_int), self%basis%bas))
   allocate(cdp(dj, di, 3), cov(dj, di))
   stat = libcint_eval_1e(cov, &
      & [jsh-1, ish-1], self%basis%atm, self%basis%bas, self%basis%env)
   stat = libcint_eval_dipole(cdp, [jsh-1, ish-1], self%basis%atm, &
      & self%basis%bas, self%basis%env)
   do iao = 1, di
      do jao = 1, dj
         overlap(jao, iao) = real(cov(jao, iao), wp)
         dpint(:, jao, iao) = real(cdp(jao, iao, :), wp)
      end do
   end do
end subroutine dipole_libcint

!> Evaluate overlap, dipole, and quadrupole integrals for a shell pair
subroutine multipole_libcint(self, cgtoj, cgtoi, jsh, ish, r2, vec, intcut, overlap, &
      & dpint, qpint)
   !> Libcint integral evaluator
   class(libcint_integral_type), intent(in) :: self
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
   real(c_double), allocatable :: cqp(:, :, :)
   real(wp) :: raw(6), trace
   integer, parameter :: qmap(6) = [1, 2, 5, 3, 6, 9]
   integer :: dj, di, jao, iao, stat

   call self%dipole_cgto(cgtoj, cgtoi, jsh, ish, r2, vec, intcut, overlap, dpint)
   dj = int(cint_cgto_spheric(int(jsh - 1, c_int), self%basis%bas))
   di = int(cint_cgto_spheric(int(ish - 1, c_int), self%basis%bas))

   allocate(cqp(dj, di, 9))
   stat = libcint_eval_quadrupole(cqp, [jsh-1, ish-1], self%basis%atm, &
      & self%basis%bas, self%basis%env)
   do iao = 1, di
      do jao = 1, dj
         raw = real(cqp(jao, iao, qmap), wp)
         trace = 0.5_wp*(raw(1) + raw(3) + raw(6))
         qpint(:, jao, iao) = 1.5_wp*raw
         qpint([1, 3, 6], jao, iao) = qpint([1, 3, 6], jao, iao) - trace
      end do
   end do
end subroutine multipole_libcint

!> Evaluate multipole integrals and their nuclear derivatives for a shell pair
subroutine multipole_grad_libcint(self, cgtoj, cgtoi, jsh, ish, r2, vec, intcut, &
      & overlap, dpint, qpint, doverlap, ddpintj, dqpintj, ddpinti, dqpinti)
   !> Libcint integral evaluator
   class(libcint_integral_type), intent(in) :: self
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
   real(c_double), allocatable :: covg(:, :, :)
   real(c_double), allocatable :: cdj(:, :, :, :), cdi(:, :, :, :)
   real(c_double), allocatable :: cqj(:, :, :, :), cqi(:, :, :, :)
   real(c_double), allocatable :: sdj(:, :, :, :), sdi(:, :, :, :)
   real(c_double), allocatable :: sqj(:, :, :, :), sqi(:, :, :, :)
   real(wp) :: raw(6), trace
   integer, parameter :: qmap(6) = [1, 2, 5, 3, 6, 9]
   integer :: dj, di, jao, iao, ic, ider, stat

   call self%multipole_cgto(cgtoj, cgtoi, jsh, ish, r2, vec, intcut, &
      & overlap, dpint, qpint)
   dj = int(cint_cgto_spheric(int(jsh - 1, c_int), self%basis%bas))
   di = int(cint_cgto_spheric(int(ish - 1, c_int), self%basis%bas))
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
         do ider = 1, 3
            doverlap(ider, jao, iao) = real(covg(jao, iao, ider), wp)
            ddpinti(ider, :, jao, iao) = real(cdi(jao, iao, :, ider), wp)
            ddpintj(ider, :, jao, iao) = real(sdi(iao, jao, :, ider), wp)
            do ic = 1, 6
               raw(ic) = real(cqi(jao, iao, qmap(ic), ider), wp)
            end do
            trace = 0.5_wp*(raw(1) + raw(3) + raw(6))
            dqpinti(ider, :, jao, iao) = 1.5_wp*raw
            dqpinti(ider, [1, 3, 6], jao, iao) = &
               & dqpinti(ider, [1, 3, 6], jao, iao) - trace
            do ic = 1, 6
               raw(ic) = real(sqi(iao, jao, qmap(ic), ider), wp)
            end do
            trace = 0.5_wp*(raw(1) + raw(3) + raw(6))
            dqpintj(ider, :, jao, iao) = 1.5_wp*raw
            dqpintj(ider, [1, 3, 6], jao, iao) = &
               & dqpintj(ider, [1, 3, 6], jao, iao) - trace
         end do
      end do
   end do
end subroutine multipole_grad_libcint

!> Evaluate <i| OVLP |j>
function libcint_eval_1e(out, shls, atm, bas, env) result(stat)
   !> Overlap-integral block
   real(c_double), contiguous, intent(out) :: out(:, :)
   !> Zero-based indices of the bra and ket shells
   integer, intent(in) :: shls(2)
   !> Libcint atom table
   integer(c_int), contiguous, intent(in) :: atm(:, :)
   !> Libcint basis-shell table
   integer(c_int), contiguous, intent(in) :: bas(:, :)
   !> Libcint floating-point environment array
   real(c_double), contiguous, intent(in) :: env(:)
   !> Libcint return status
   integer :: stat

   integer(c_int) :: cshls(2), dims(2)
   integer :: di, dj

   cshls = int(shls, c_int) ! convert to C-int for libcint
   dims = int([size(out, 1), size(out, 2)], c_int)
   ! Determine how many spherical basis functions are in each shell
   di = int(cint_cgto_spheric(int(shls(1), c_int), bas))
   dj = int(cint_cgto_spheric(int(shls(2), c_int), bas))
   if (di < 0 .or. dj < 0 .or. size(out, 1) < di .or. size(out, 2) < dj) then
      stat = -1
      return
   end if

   out(:, :) = 0.0_c_double
   stat = int(int1e_ovlp_sph(out, dims, cshls, atm, int(size(atm, 2), c_int), &
      & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
end function libcint_eval_1e

!> Evaluate (i j|k) for two AO shells and one auxiliary Gaussian shell.
function libcint_eval_3c2e(out, shls, atm, bas, env) result(stat)
   real(c_double), contiguous, intent(out) :: out(:, :, :)
   integer, intent(in) :: shls(3)
   integer(c_int), contiguous, intent(in) :: atm(:, :)
   integer(c_int), contiguous, intent(in) :: bas(:, :)
   real(c_double), contiguous, intent(in) :: env(:)
   integer :: stat

   integer(c_int) :: cshls(3), dims(3)
   integer :: di, dj, dk

   di = int(cint_cgto_spheric(int(shls(1), c_int), bas))
   dj = int(cint_cgto_spheric(int(shls(2), c_int), bas))
   dk = int(cint_cgto_spheric(int(shls(3), c_int), bas))
   if (di < 0 .or. dj < 0 .or. dk < 0 .or. size(out, 1) < di .or. &
      & size(out, 2) < dj .or. size(out, 3) < dk) then
      stat = -1
      return
   end if

   cshls = int(shls, c_int)
   dims = int([size(out, 1), size(out, 2), size(out, 3)], c_int)
   out = 0.0_c_double
   stat = int(int3c2e_sph(out, dims, cshls, atm, int(size(atm, 2), c_int), &
      & bas, int(size(bas, 2), c_int), env, c_null_ptr, c_null_ptr))
end function libcint_eval_3c2e

!> Evaluate <i| R |j>
function libcint_eval_dipole(out, shls, atm, bas, env) result(stat)
   !> Dipole-integral block, with Cartesian operator component last
   real(c_double), contiguous, intent(out) :: out(:, :, :)
   !> Zero-based indices of the bra and ket shells
   integer, intent(in) :: shls(2)
   !> Libcint atom table
   integer(c_int), contiguous, intent(in) :: atm(:, :)
   !> Libcint basis-shell table
   integer(c_int), contiguous, intent(in) :: bas(:, :)
   !> Libcint floating-point environment array
   real(c_double), contiguous, intent(in) :: env(:)
   !> Libcint return status
   integer :: stat
   integer(c_int) :: cshls(2), dims(2)
   integer :: di, dj

   di = int(cint_cgto_spheric(int(shls(1), c_int), bas))
   dj = int(cint_cgto_spheric(int(shls(2), c_int), bas))
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

!> Evaluate <i| R R |j>
function libcint_eval_quadrupole(out, shls, atm, bas, env) result(stat)
   !> Quadrupole-integral block, with Cartesian operator component last
   real(c_double), contiguous, intent(out) :: out(:, :, :)
   !> Zero-based indices of the bra and ket shells
   integer, intent(in) :: shls(2)
   !> Libcint atom table
   integer(c_int), contiguous, intent(in) :: atm(:, :)
   !> Libcint basis-shell table
   integer(c_int), contiguous, intent(in) :: bas(:, :)
   !> Libcint floating-point environment array
   real(c_double), contiguous, intent(in) :: env(:)
   !> Libcint return status
   integer :: stat
   integer(c_int) :: cshls(2), dims(2)
   integer :: di, dj

   di = int(cint_cgto_spheric(int(shls(1), c_int), bas))
   dj = int(cint_cgto_spheric(int(shls(2), c_int), bas))
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

!> Evaluate <NABLA i| OVLP |R j> and <i| OVLP |NABLA R j>
function libcint_eval_dipole_gradient(out_bra, out_ket, shls, atm, bas, env) result(stat)
   !> Dipole derivative with respect to the bra center
   real(c_double), contiguous, intent(out) :: out_bra(:, :, :, :)
   !> Dipole derivative with respect to the ket center
   real(c_double), contiguous, intent(out) :: out_ket(:, :, :, :)
   !> Zero-based indices of the bra and ket shells
   integer, intent(in) :: shls(2)
   !> Libcint atom table
   integer(c_int), contiguous, intent(in) :: atm(:, :)
   !> Libcint basis-shell table
   integer(c_int), contiguous, intent(in) :: bas(:, :)
   !> Libcint floating-point environment array
   real(c_double), contiguous, intent(in) :: env(:)
   !> Libcint return status
   integer :: stat
   integer(c_int) :: cshls(2), dims(2)
   integer :: di, dj, stat_bra

   di = int(cint_cgto_spheric(int(shls(1), c_int), bas))
   dj = int(cint_cgto_spheric(int(shls(2), c_int), bas))
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

!> Evaluate <NABLA i| OVLP |R R j> and <i| OVLP |NABLA R R j>
function libcint_eval_quadrupole_gradient(out_bra, out_ket, shls, atm, bas, &
      & env) result(stat)
   !> Quadrupole derivative with respect to the bra center
   real(c_double), contiguous, intent(out) :: out_bra(:, :, :, :)
   !> Quadrupole derivative with respect to the ket center
   real(c_double), contiguous, intent(out) :: out_ket(:, :, :, :)
   !> Zero-based indices of the bra and ket shells
   integer, intent(in) :: shls(2)
   !> Libcint atom table
   integer(c_int), contiguous, intent(in) :: atm(:, :)
   !> Libcint basis-shell table
   integer(c_int), contiguous, intent(in) :: bas(:, :)
   !> Libcint floating-point environment array
   real(c_double), contiguous, intent(in) :: env(:)
   !> Libcint return status
   integer :: stat
   integer(c_int) :: cshls(2), dims(2)
   integer :: di, dj, stat_bra

   di = int(cint_cgto_spheric(int(shls(1), c_int), bas))
   dj = int(cint_cgto_spheric(int(shls(2), c_int), bas))
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

!> Evaluate <i| OVLP |NABLA j>
function libcint_eval_overlap_gradient(out, shls, atm, bas, env) result(stat)
   !> Overlap derivative with respect to the ket center
   real(c_double), contiguous, intent(out) :: out(:, :, :)
   !> Zero-based indices of the bra and ket shells
   integer, intent(in) :: shls(2)
   !> Libcint atom table
   integer(c_int), contiguous, intent(in) :: atm(:, :)
   !> Libcint basis-shell table
   integer(c_int), contiguous, intent(in) :: bas(:, :)
   !> Libcint floating-point environment array
   real(c_double), contiguous, intent(in) :: env(:)
   !> Libcint return status
   integer :: stat
   integer(c_int) :: cshls(2), dims(2)
   integer :: di, dj

   di = int(cint_cgto_spheric(int(shls(1), c_int), bas))
   dj = int(cint_cgto_spheric(int(shls(2), c_int), bas))
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

end module tblite_integral_libcint
