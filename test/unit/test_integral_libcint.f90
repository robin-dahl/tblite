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

module test_integral_libcint
   use, intrinsic :: iso_c_binding, only : c_double
   use mctc_env, only : wp
   use mctc_env_testing, only : check, error_type, new_unittest, test_failed, &
      & unittest_type
   use mctc_io, only : new, structure_type
   use mstore, only : get_structure
   use tblite_basis_type, only : basis_type, cgto_type, get_cutoff, new_basis, &
      & new_cgto
   use tblite_context_type, only : context_type
   use tblite_features, only : tblite_use_libcint
   use tblite_integral_handler, only : enum_integral_handler
   use tblite_wavefunction, only : new_wavefunction, wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   use tblite_integral_libcint
   use tblite_integral_native, only : native_integral_type
   use tblite_integral_native_integrals, only : get_overlap, overlap_grad_cgto
   implicit none
   private

   public :: collect_integral_libcint

   real(wp), parameter :: acc = 0.01_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp
   integer, parameter :: max_shell = 7
   integer, parameter :: max_block = max_shell**2

contains


!> Collect all exported libcint unit tests
subroutine collect_integral_libcint(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   if (tblite_use_libcint) then
      testsuite = [ &
         new_unittest("integral-handler-consistency", test_integral_handler_consistency), &
         new_unittest("energy-consistency", test_energy_consistency), &
         new_unittest("overlap-consistency", test_overlap_consistency), &
         new_unittest("overlap-gradient-consistency", &
            & test_overlap_gradient_consistency), &
         new_unittest("dipole-consistency", test_dipole_consistency), &
         new_unittest("dipole-gradient-consistency", &
            & test_dipole_gradient_consistency), &
         new_unittest("quadrupole-consistency", test_quadrupole_consistency), &
         new_unittest("quadrupole-gradient-consistency", &
            & test_quadrupole_gradient_consistency) &
         ]
   else
      testsuite = [new_unittest("disabled", test_disabled)]
   end if
end subroutine collect_integral_libcint

!> Construct a molecular basis containing s, p, d, and f shells for comparisons
subroutine make_comparison_basis(mol, basis, libcint)
   !> Molecular structure data
   type(structure_type), intent(out) :: mol
   !> Basis set information
   type(basis_type), intent(out) :: basis
   !> Libcint integral evaluator
   type(libcint_integral_type), intent(out) :: libcint
   type(cgto_type), allocatable :: cgto(:, :)
   integer, allocatable :: nshell(:)

   call new(mol, [6, 6], reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.8_wp, -0.5_wp, 1.1_wp], [3, 2]))
   allocate(nshell(mol%nid), cgto(4, mol%nid))
   nshell = 4
   call new_cgto(cgto(1, 1), 2, 0, [1.4_wp, 0.35_wp], &
      & [0.65_wp, 0.45_wp], .true.)
   call new_cgto(cgto(2, 1), 2, 1, [1.1_wp, 0.28_wp], &
      & [0.60_wp, 0.50_wp], .true.)
   call new_cgto(cgto(3, 1), 2, 2, [0.9_wp, 0.22_wp], &
      & [0.55_wp, 0.52_wp], .true.)
   call new_cgto(cgto(4, 1), 2, 3, [0.8_wp, 0.20_wp], &
      & [0.50_wp, 0.48_wp], .true.)
   call new_basis(basis, mol, nshell, cgto, 1.0_wp)
   call libcint%initialize_integral(mol, basis)
end subroutine make_comparison_basis

!> Compare native and libcint evaluators through the common handler interface
subroutine test_integral_handler_consistency(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   type(basis_type) :: basis
   type(cgto_type) :: cgtoj, cgtoi
   type(native_integral_type) :: nativeint
   type(libcint_integral_type) :: libcint
   real(wp) :: sn(max_block), sc(max_block)
   real(wp) :: dn(3, max_block), dc(3, max_block)
   real(wp) :: qn(6, max_block), qc(6, max_block), vec(3), r2
   real(wp) :: dsn(3, max_block), dsc(3, max_block)
   real(wp) :: ddjn(3, 3, max_block), ddjc(3, 3, max_block)
   real(wp) :: ddin(3, 3, max_block), ddic(3, 3, max_block)
   real(wp) :: dqjn(3, 6, max_block), dqjc(3, 6, max_block)
   real(wp) :: dqin(3, 6, max_block), dqic(3, 6, max_block)
   integer :: iat, jat, ish, jsh, n

   call make_comparison_basis(mol, basis, libcint)
   call nativeint%initialize_integral(mol, basis)

   do ish = 1, basis%nsh
      iat = basis%sh2at(ish)
      cgtoi = basis%cgto(ish-basis%ish_at(iat), mol%id(iat))
      do jsh = 1, basis%nsh
         jat = basis%sh2at(jsh)
         cgtoj = basis%cgto(jsh-basis%ish_at(jat), mol%id(jat))
         n = basis%nao_sh(ish)*basis%nao_sh(jsh)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = sum(vec**2)

         call nativeint%dipole_cgto(cgtoj, cgtoi, jsh, ish, r2, vec, &
            & basis%intcut, sn, dn)
         call libcint%dipole_cgto(cgtoj, cgtoi, jsh, ish, r2, vec, &
            & basis%intcut, sc, dc)
         call check(error, all(abs(sn(:n) - sc(:n)) < thr), &
            & message="Overlap integrals do not match")
         if (allocated(error)) return
         call check(error, all(abs(dn(:, :n) - dc(:, :n)) < thr), &
            & message="Dipole integrals do not match")
         if (allocated(error)) return

         call nativeint%multipole_cgto(cgtoj, cgtoi, jsh, ish, r2, vec, &
            & basis%intcut, sn, dn, qn)
         call libcint%multipole_cgto(cgtoj, cgtoi, jsh, ish, r2, vec, &
            & basis%intcut, sc, dc, qc)
         call check(error, all(abs(sn(:n) - sc(:n)) < thr), &
            & message="Overlap integrals do not match")
         if (allocated(error)) return
         call check(error, all(abs(dn(:, :n) - dc(:, :n)) < thr), &
            & message="Dipole integrals do not match")
         if (allocated(error)) return
         call check(error, all(abs(qn(:, :n) - qc(:, :n)) < thr), &
            & message="Quadrupole integrals do not match")
         if (allocated(error)) return

         call nativeint%multipole_grad_cgto(cgtoj, cgtoi, jsh, ish, r2, vec, &
            & basis%intcut, sn, dn, qn, dsn, ddjn, dqjn, ddin, dqin)
         call libcint%multipole_grad_cgto(cgtoj, cgtoi, jsh, ish, r2, vec, &
            & basis%intcut, sc, dc, qc, dsc, ddjc, dqjc, ddic, dqic)
         call check(error, all(abs(dsn(:, :n) - dsc(:, :n)) < thr), &
            & message="Overlap gradients do not match")
         if (allocated(error)) return
         call check(error, all(abs(ddjn(:, :, :n) - ddjc(:, :, :n)) < thr), &
            & message="Dipole gradients on center j do not match")
         if (allocated(error)) return
         call check(error, all(abs(ddin(:, :, :n) - ddic(:, :, :n)) < thr), &
            & message="Dipole gradients on center i do not match")
         if (allocated(error)) return
         call check(error, all(abs(dqjn(:, :, :n) - dqjc(:, :, :n)) < thr), &
            & message="Quadrupole gradients on center j do not match")
         if (allocated(error)) return
         call check(error, all(abs(dqin(:, :, :n) - dqic(:, :, :n)) < thr), &
            & message="Quadrupole gradients on center i do not match")
         if (allocated(error)) return
      end do
   end do
end subroutine test_integral_handler_consistency

!> Compare complete GFN2-xTB energies using native and libcint integral handlers
subroutine test_energy_consistency(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx_native
   type(context_type) :: ctx_libcint
   type(structure_type) :: mol
   type(xtb_calculator) :: calc_native
   type(xtb_calculator) :: calc_libcint
   type(wavefunction_type) :: wfn_native
   type(wavefunction_type) :: wfn_libcint
   real(wp) :: energy_native
   real(wp) :: energy_libcint

   call get_structure(mol, "MB16-43", "01")

   ! Construct independent calculators so both SCF calculations start identically.
   call new_gfn2_calculator(calc_native, mol, error)
   if (allocated(error)) return
   call calc_native%set_integral_handler(mol, error, enum_integral_handler%native)
   if (allocated(error)) return
   call new_wavefunction(wfn_native, mol%nat, calc_native%bas%nsh, &
      & calc_native%bas%nao, 1, kt)

   energy_native = 0.0_wp
   call xtb_singlepoint(ctx_native, mol, calc_native, wfn_native, acc, &
      & energy_native, verbosity=0)
   call check(error, .not. ctx_native%failed(), &
      & message="Native integral handler calculation failed")
   if (allocated(error)) return

   call new_gfn2_calculator(calc_libcint, mol, error)
   if (allocated(error)) return
   call calc_libcint%set_integral_handler(mol, error, enum_integral_handler%libcint)
   if (allocated(error)) return
   call new_wavefunction(wfn_libcint, mol%nat, calc_libcint%bas%nsh, &
      & calc_libcint%bas%nao, 1, kt)

   energy_libcint = 0.0_wp
   call xtb_singlepoint(ctx_libcint, mol, calc_libcint, wfn_libcint, acc, &
      & energy_libcint, verbosity=0)
   call check(error, .not. ctx_libcint%failed(), &
      & message="Libcint integral handler calculation failed")
   if (allocated(error)) return

   if (abs(energy_libcint - energy_native) > thr2) then
      call test_failed(error, &
         & "GFN2-xTB energies from native and libcint integrals do not match")
      print '("Native energy:  ", es21.14)', energy_native
      print '("Libcint energy: ", es21.14)', energy_libcint
      print '("Difference:     ", es21.14)', energy_libcint-energy_native
   end if
end subroutine test_energy_consistency

!> Compare the complete native and libcint overlap matrices
subroutine test_overlap_consistency(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(basis_type) :: basis
   type(libcint_integral_type) :: libcint
   real(wp), allocatable :: ref(:, :), cint_overlap(:, :)
   real(wp) :: trans(3, 1), cutoff
   real(c_double) :: block(max_shell, max_shell)
   integer :: ish, jsh, ii, jj, di, dj, iao, jao, stat

   call make_comparison_basis(mol, basis, libcint)

   call check(error, size(libcint%basis%atm, 2), mol%nat, &
      & message="Number of libcint atoms does not match tblite basis")
   if (allocated(error)) return
   call check(error, size(libcint%basis%bas, 2), basis%nsh, &
      & message="Number of libcint shells does not match tblite basis")
   if (allocated(error)) return

   cutoff = get_cutoff(basis)
   trans = 0.0_wp
   allocate(ref(basis%nao, basis%nao), cint_overlap(basis%nao, basis%nao), &
      & source=0.0_wp)
   call get_overlap(mol, trans, cutoff, basis, ref)

   do ish = 1, basis%nsh
      ii = basis%iao_sh(ish)
      di = basis%nao_sh(ish)
      do jsh = 1, basis%nsh
         jj = basis%iao_sh(jsh)
         dj = basis%nao_sh(jsh)
         stat = libcint_eval_1e(&
            & block, [ish-1, jsh-1], libcint%basis%atm, libcint%basis%bas, &
            & libcint%basis%env)
         call check(error, stat >= 0, message="Libcint overlap evaluation failed")
         if (allocated(error)) return
         do iao = 1, di
            do jao = 1, dj
               cint_overlap(ii+iao, jj+jao) = real(block(iao, jao), wp)
            end do
         end do
      end do
   end do

   do ii = 1, basis%nao
      do jj = 1, basis%nao
         call check(error, cint_overlap(jj, ii), ref(jj, ii), thr=thr, &
            & message="Libcint overlap does not match native overlap")
         if (allocated(error)) return
      end do
   end do
end subroutine test_overlap_consistency

!> Compare native and libcint overlap gradients for every shell pair
subroutine test_overlap_gradient_consistency(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   type(basis_type) :: basis
   type(libcint_integral_type) :: libcint
   real(wp) :: ref_s(max_shell, max_shell), ref_g(3, max_shell, max_shell)
   real(wp) :: vec(3), r2
   real(c_double) :: cint_g(max_shell, max_shell, 3)
   integer :: ish, jsh, iat, jat, isp, jsp, ilsh, jlsh
   integer :: ii, jj, di, dj, iao, jao, ic, stat

   call make_comparison_basis(mol, basis, libcint)
   do ish = 1, basis%nsh
      iat = basis%sh2at(ish)
      isp = mol%id(iat)
      ilsh = ish - basis%ish_at(iat)
      ii = basis%iao_sh(ish)
      di = basis%nao_sh(ish)
      do jsh = 1, basis%nsh
         jat = basis%sh2at(jsh)
         jsp = mol%id(jat)
         jlsh = jsh - basis%ish_at(jat)
         jj = basis%iao_sh(jsh)
         dj = basis%nao_sh(jsh)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = sum(vec**2)
         call overlap_grad_cgto(basis%cgto(jlsh, jsp), basis%cgto(ilsh, isp), &
            & r2, vec, basis%intcut, ref_s(1:dj, 1:di), ref_g(:, 1:dj, 1:di))
         ! The reversed order gives libcint output (j,i), with shell i in the
         ! second position differentiated by int1e_ovlpip_sph.
         stat = libcint_eval_overlap_gradient(cint_g, [jsh-1, ish-1], &
            & libcint%basis%atm, libcint%basis%bas, libcint%basis%env)
         call check(error, stat >= 0, &
            & message="Libcint overlap gradient evaluation failed")
         if (allocated(error)) return
         do ic = 1, 3
            do iao = 1, di
               do jao = 1, dj
                  call check(error, real(cint_g(jao, iao, ic), wp), &
                     & ref_g(ic, jao, iao), thr=thr, &
                     & message="Libcint overlap gradient does not match native result")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end do
end subroutine test_overlap_gradient_consistency

!> Compare native and libcint dipole integrals for every shell pair
subroutine test_dipole_consistency(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   type(basis_type) :: basis
   type(libcint_integral_type) :: libcint
   type(native_integral_type) :: nativeint
   real(wp) :: ref_s(max_shell, max_shell), ref_d(3, max_shell, max_shell)
   real(wp) :: vec(3), r2
   real(c_double) :: cint_d(max_shell, max_shell, 3)
   integer :: ish, jsh, iat, jat, isp, jsp, ilsh, jlsh
   integer :: ii, jj, di, dj, iao, jao, ic, stat

   call make_comparison_basis(mol, basis, libcint)
   do ish = 1, basis%nsh
      iat = basis%sh2at(ish)
      isp = mol%id(iat)
      ilsh = ish - basis%ish_at(iat)
      ii = basis%iao_sh(ish)
      di = basis%nao_sh(ish)
      do jsh = 1, basis%nsh
         jat = basis%sh2at(jsh)
         jsp = mol%id(jat)
         jlsh = jsh - basis%ish_at(jat)
         jj = basis%iao_sh(jsh)
         dj = basis%nao_sh(jsh)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = sum(vec**2)
         call nativeint%dipole_cgto(basis%cgto(jlsh, jsp), &
            & basis%cgto(ilsh, isp), jsh, ish, r2, vec, basis%intcut, &
            & ref_s(1:dj, 1:di), ref_d(:, 1:dj, 1:di))
         ! Reversing the libcint shells makes origj the tblite ket centre i.
         stat = libcint_eval_dipole(cint_d, [jsh-1, ish-1], &
            & libcint%basis%atm, libcint%basis%bas, libcint%basis%env)
         call check(error, stat >= 0, message="Libcint dipole evaluation failed")
         if (allocated(error)) return
         do ic = 1, 3
            do iao = 1, di
               do jao = 1, dj
                  call check(error, real(cint_d(jao, iao, ic), wp), &
                     & ref_d(ic, jao, iao), thr=thr, &
                     & message="Libcint dipole integral does not match native result")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end do
end subroutine test_dipole_consistency

!> Compare native and libcint dipole gradients for both Gaussian centers
subroutine test_dipole_gradient_consistency(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   type(basis_type) :: basis
   type(libcint_integral_type) :: libcint
   type(native_integral_type) :: nativeint
   real(wp) :: ref_s(max_shell, max_shell), ref_d(3, max_shell, max_shell)
   real(wp) :: ref_q(6, max_shell, max_shell), ref_g(3, max_shell, max_shell)
   real(wp) :: ref_dj(3, 3, max_shell, max_shell), ref_di(3, 3, max_shell, max_shell)
   real(wp) :: ref_qj(3, 6, max_shell, max_shell), ref_qi(3, 6, max_shell, max_shell)
   real(wp) :: vec(3), r2
   real(c_double) :: cint_j(max_shell, max_shell, 3, 3)
   real(c_double) :: cint_i(max_shell, max_shell, 3, 3)
   real(c_double) :: swap_i(max_shell, max_shell, 3, 3)
   real(c_double) :: swap_j(max_shell, max_shell, 3, 3)
   integer :: ish, jsh, iat, jat, isp, jsp, ilsh, jlsh
   integer :: di, dj, iao, jao, im, ider, stat

   call make_comparison_basis(mol, basis, libcint)
   do ish = 1, basis%nsh
      iat = basis%sh2at(ish)
      isp = mol%id(iat)
      ilsh = ish - basis%ish_at(iat)
      di = basis%nao_sh(ish)
      do jsh = 1, basis%nsh
         jat = basis%sh2at(jsh)
         jsp = mol%id(jat)
         jlsh = jsh - basis%ish_at(jat)
         dj = basis%nao_sh(jsh)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = sum(vec**2)
         call nativeint%multipole_grad_cgto(basis%cgto(jlsh, jsp), &
            & basis%cgto(ilsh, isp), jsh, ish, r2, vec, basis%intcut, &
            & ref_s(1:dj, 1:di), ref_d(:, 1:dj, 1:di), ref_q(:, 1:dj, 1:di), &
            & ref_g(:, 1:dj, 1:di), ref_dj(:, :, 1:dj, 1:di), &
            & ref_qj(:, :, 1:dj, 1:di), ref_di(:, :, 1:dj, 1:di), &
            & ref_qi(:, :, 1:dj, 1:di))
         ! Evaluate both shell orders because the moment origin follows the ket shell.
         stat = libcint_eval_dipole_gradient(cint_j, cint_i, [jsh-1, ish-1], &
            & libcint%basis%atm, libcint%basis%bas, libcint%basis%env)
         call check(error, stat >= 0, &
            & message="Libcint dipole gradient evaluation failed")
         if (allocated(error)) return
         stat = libcint_eval_dipole_gradient(swap_i, swap_j, [ish-1, jsh-1], &
            & libcint%basis%atm, libcint%basis%bas, libcint%basis%env)
         call check(error, stat >= 0, &
            & message="Libcint reversed dipole gradient evaluation failed")
         if (allocated(error)) return
         do ider = 1, 3
            do im = 1, 3
               do iao = 1, di
                  do jao = 1, dj
            ! For origj=i, d/dRi=d/dvec and d/dRj=-d/dvec.
            call check(error, real(cint_j(jao, iao, im, ider), wp), &
               & -ref_di(ider, im, jao, iao), thr=thr, &
               & message="Libcint bra-center dipole gradient does not match")
            if (allocated(error)) return
            call check(error, real(cint_i(jao, iao, im, ider), wp), &
               & ref_di(ider, im, jao, iao), thr=thr, &
               & message="Libcint ket-center dipole gradient does not match")
            if (allocated(error)) return
            ! For origj=j the reversed block provides tblite's j-centred operator.
            call check(error, real(swap_i(iao, jao, im, ider), wp), &
               & ref_dj(ider, im, jao, iao), thr=thr, &
               & message="Libcint reversed ket-center dipole gradient does not match")
            if (allocated(error)) return
            call check(error, real(swap_j(iao, jao, im, ider), wp), &
               & -ref_dj(ider, im, jao, iao), thr=thr, &
               & message="Libcint reversed bra-center dipole gradient does not match")
            if (allocated(error)) return
                  end do
               end do
            end do
         end do
      end do
   end do
end subroutine test_dipole_gradient_consistency

!> Compare native and libcint traceless quadrupole integrals
subroutine test_quadrupole_consistency(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   type(basis_type) :: basis
   type(libcint_integral_type) :: libcint
   type(native_integral_type) :: nativeint
   real(wp) :: ref_s(max_shell, max_shell), ref_d(3, max_shell, max_shell)
   real(wp) :: ref_q(6, max_shell, max_shell)
   real(wp) :: vec(3), r2, raw(6), quad(6), trace
   real(c_double) :: cint_q(max_shell, max_shell, 9)
   integer :: ish, jsh, iat, jat, isp, jsp, ilsh, jlsh
   integer :: ii, jj, di, dj, iao, jao, ic, stat
   ! Map libcint's full Cartesian tensor to tblite's xx, xy, yy, xz, yz, zz order.
   integer, parameter :: qmap(6) = [1, 2, 5, 3, 6, 9]

   call make_comparison_basis(mol, basis, libcint)
   do ish = 1, basis%nsh
      iat = basis%sh2at(ish)
      isp = mol%id(iat)
      ilsh = ish - basis%ish_at(iat)
      ii = basis%iao_sh(ish)
      di = basis%nao_sh(ish)
      do jsh = 1, basis%nsh
         jat = basis%sh2at(jsh)
         jsp = mol%id(jat)
         jlsh = jsh - basis%ish_at(jat)
         jj = basis%iao_sh(jsh)
         dj = basis%nao_sh(jsh)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = sum(vec**2)
         call nativeint%multipole_cgto(basis%cgto(jlsh, jsp), &
            & basis%cgto(ilsh, isp), jsh, ish, r2, vec, basis%intcut, &
            & ref_s(1:dj, 1:di), ref_d(:, 1:dj, 1:di), ref_q(:, 1:dj, 1:di))
         stat = libcint_eval_quadrupole(cint_q, [jsh-1, ish-1], &
            & libcint%basis%atm, libcint%basis%bas, libcint%basis%env)
         call check(error, stat >= 0, message="Libcint quadrupole evaluation failed")
         if (allocated(error)) return
         do iao = 1, di
            do jao = 1, dj
               do ic = 1, 6
                  raw(ic) = real(cint_q(jao, iao, qmap(ic)), wp)
               end do
               ! Apply the same traceless Cartesian transformation as tblite.
               trace = 0.5_wp*(raw(1) + raw(3) + raw(6))
               quad = 1.5_wp*raw
               quad([1, 3, 6]) = quad([1, 3, 6]) - trace
               do ic = 1, 6
                  call check(error, quad(ic), ref_q(ic, jao, iao), thr=thr, &
                     & message="Libcint quadrupole integral does not match native result")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end do
end subroutine test_quadrupole_consistency

!> Compare native and libcint quadrupole gradients for both Gaussian centers
subroutine test_quadrupole_gradient_consistency(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   type(basis_type) :: basis
   type(libcint_integral_type) :: libcint
   type(native_integral_type) :: nativeint
   real(wp) :: ref_s(max_shell, max_shell), ref_d(3, max_shell, max_shell)
   real(wp) :: ref_q(6, max_shell, max_shell), ref_g(3, max_shell, max_shell)
   real(wp) :: ref_dj(3, 3, max_shell, max_shell), ref_di(3, 3, max_shell, max_shell)
   real(wp) :: ref_qj(3, 6, max_shell, max_shell), ref_qi(3, 6, max_shell, max_shell)
   real(wp) :: raw(6), quad(6), trace
   real(wp) :: vec(3), r2
   real(c_double) :: cint_j(max_shell, max_shell, 9, 3)
   real(c_double) :: cint_i(max_shell, max_shell, 9, 3)
   real(c_double) :: swap_i(max_shell, max_shell, 9, 3)
   real(c_double) :: swap_j(max_shell, max_shell, 9, 3)
   integer :: ish, jsh, iat, jat, isp, jsp, ilsh, jlsh
   integer :: di, dj, iao, jao, im, ider, stat, center
   ! Map libcint's full Cartesian tensor to tblite's xx, xy, yy, xz, yz, zz order.
   integer, parameter :: qmap(6) = [1, 2, 5, 3, 6, 9]

   call make_comparison_basis(mol, basis, libcint)
   do ish = 1, basis%nsh
      iat = basis%sh2at(ish)
      isp = mol%id(iat)
      ilsh = ish - basis%ish_at(iat)
      di = basis%nao_sh(ish)
      do jsh = 1, basis%nsh
         jat = basis%sh2at(jsh)
         jsp = mol%id(jat)
         jlsh = jsh - basis%ish_at(jat)
         dj = basis%nao_sh(jsh)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = sum(vec**2)
         call nativeint%multipole_grad_cgto(basis%cgto(jlsh, jsp), &
            & basis%cgto(ilsh, isp), jsh, ish, r2, vec, basis%intcut, &
            & ref_s(1:dj, 1:di), ref_d(:, 1:dj, 1:di), ref_q(:, 1:dj, 1:di), &
            & ref_g(:, 1:dj, 1:di), ref_dj(:, :, 1:dj, 1:di), &
            & ref_qj(:, :, 1:dj, 1:di), ref_di(:, :, 1:dj, 1:di), &
            & ref_qi(:, :, 1:dj, 1:di))
         ! Evaluate both shell orders because the moment origin follows the ket shell.
         stat = libcint_eval_quadrupole_gradient(cint_j, cint_i, [jsh-1, ish-1], &
            & libcint%basis%atm, libcint%basis%bas, libcint%basis%env)
         call check(error, stat >= 0, &
            & message="Libcint quadrupole gradient evaluation failed")
         if (allocated(error)) return
         stat = libcint_eval_quadrupole_gradient(swap_i, swap_j, [ish-1, jsh-1], &
            & libcint%basis%atm, libcint%basis%bas, libcint%basis%env)
         call check(error, stat >= 0, &
            & message="Libcint reversed quadrupole gradient evaluation failed")
         if (allocated(error)) return
         do center = 1, 4
            do ider = 1, 3
               do iao = 1, di
                  do jao = 1, dj
            do im = 1, 6
               if (center == 1) then
                  raw(im) = real(cint_j(jao, iao, qmap(im), ider), wp)
               else if (center == 2) then
                  raw(im) = real(cint_i(jao, iao, qmap(im), ider), wp)
               else if (center == 3) then
                  raw(im) = real(swap_i(iao, jao, qmap(im), ider), wp)
               else
                  raw(im) = real(swap_j(iao, jao, qmap(im), ider), wp)
               end if
            end do
            ! Apply the same traceless Cartesian transformation as tblite.
            trace = 0.5_wp*(raw(1) + raw(3) + raw(6))
            quad = 1.5_wp*raw
            quad([1, 3, 6]) = quad([1, 3, 6]) - trace
            do im = 1, 6
               if (center == 1) then
                  call check(error, quad(im), -ref_qi(ider, im, jao, iao), &
                     & thr=thr, message= &
                     & "Libcint bra-center quadrupole gradient does not match")
               else if (center == 2) then
                  call check(error, quad(im), ref_qi(ider, im, jao, iao), &
                     & thr=thr, message= &
                     & "Libcint ket-center quadrupole gradient does not match")
               else if (center == 3) then
                  call check(error, quad(im), ref_qj(ider, im, jao, iao), &
                     & thr=thr, message= &
                     & "Libcint reversed ket-center quadrupole gradient does not match")
               else
                  call check(error, quad(im), -ref_qj(ider, im, jao, iao), &
                     & thr=thr, message= &
                     & "Libcint reversed bra-center quadrupole gradient does not match")
               end if
               if (allocated(error)) return
            end do
                  end do
               end do
            end do
         end do
      end do
   end do
end subroutine test_quadrupole_gradient_consistency

!> Check that libcint support is reported as disabled when unavailable
subroutine test_disabled(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check(error, .not. tblite_use_libcint, &
      & message="Build unexpectedly reports libcint support")
end subroutine test_disabled

end module test_integral_libcint
