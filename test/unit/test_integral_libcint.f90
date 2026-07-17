! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later

module test_integral_libcint
   use, intrinsic :: iso_c_binding, only : c_double
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_basis_type, only : basis_type, cgto_type, new_basis, new_cgto, get_cutoff
   use tblite_context_type, only : context_type
   use tblite_features, only : tblite_use_libcint
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator, integral_handler_native, &
      & integral_handler_libcint
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
#if TBLITE_HAS_LIBCINT
   use tblite_integral_libcint
   use tblite_integral_native, only : native_integral_type, dipole_cgto, &
      & multipole_cgto, multipole_grad_cgto
   use tblite_integral_overlap, only : get_overlap, overlap_grad_cgto
#endif
   implicit none
   private

   public :: collect_integral_libcint

   real(wp), parameter :: acc = 0.01_wp
   real(wp), parameter :: thr = 1.0e-10_wp
   real(wp), parameter :: thr2 = 1.0e-11_wp
   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp

contains

!> Collect all exported libcint unit tests
subroutine collect_integral_libcint(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

#if TBLITE_HAS_LIBCINT
   testsuite = [ &
      new_unittest("tblite-overlap-consistency", test_tblite_overlap_consistency), &
      new_unittest("tblite-overlap-gradient-consistency", test_tblite_overlap_gradient_consistency), &
      new_unittest("integral-handler-consistency", test_integral_handler_consistency), &
      new_unittest("integral-energy-consistency", test_integral_energy_consistency), &
      new_unittest("tblite-dipole-consistency", test_tblite_dipole_consistency), &
      new_unittest("tblite-quadrupole-consistency", test_tblite_quadrupole_consistency), &
      new_unittest("tblite-dipole-gradient-consistency", test_tblite_dipole_gradient_consistency), &
      new_unittest("tblite-quadrupole-gradient-consistency", test_tblite_quadrupole_gradient_consistency) &
      ]
#else
   testsuite = [ &
      new_unittest("disabled", test_disabled) &
      ]
#endif
end subroutine collect_integral_libcint

#if TBLITE_HAS_LIBCINT
!> Construct a molecular basis containing s, p, and d shells for comparisons
subroutine make_comparison_basis(mol, basis, cbasis)
   type(structure_type), intent(out) :: mol
   type(basis_type), intent(out) :: basis
   type(libcint_basis_type), intent(out) :: cbasis
   type(cgto_type), allocatable :: cgto(:, :)
   integer, allocatable :: nshell(:)

   call new(mol, [6, 6], reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.8_wp, -0.5_wp, 1.1_wp], [3, 2]))
   allocate(nshell(mol%nid), cgto(3, mol%nid))
   nshell = 3
   call new_cgto(cgto(1, 1), 2, 0, [1.4_wp, 0.35_wp], &
      & [0.65_wp, 0.45_wp], .true.)
   call new_cgto(cgto(2, 1), 2, 1, [1.1_wp, 0.28_wp], &
      & [0.60_wp, 0.50_wp], .true.)
   call new_cgto(cgto(3, 1), 2, 2, [0.9_wp, 0.22_wp], &
      & [0.55_wp, 0.52_wp], .true.)
   call new_basis(basis, mol, nshell, cgto, 1.0_wp)
   call new_libcint_basis(cbasis, mol, basis)
end subroutine make_comparison_basis

!> Compare the complete native and libcint overlap matrices
subroutine test_tblite_overlap_consistency(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(basis_type) :: basis
   type(libcint_basis_type) :: cbasis
   real(wp), allocatable :: ref(:, :), cint_overlap(:, :)
   real(wp) :: trans(3, 1), cutoff
   real(c_double) :: block(5, 5)
   integer :: ish, jsh, ii, jj, di, dj, iao, jao, stat

   call make_comparison_basis(mol, basis, cbasis)

   call check(error, size(cbasis%atm, 2), mol%nat, &
      & message="Number of libcint atoms does not match tblite basis")
   if (allocated(error)) return
   call check(error, size(cbasis%bas, 2), basis%nsh, &
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
         stat = libcint_eval_1e(LIBCINT_1E_OVERLAP, LIBCINT_SPHERICAL, &
            & block, [ish-1, jsh-1], cbasis%atm, cbasis%bas, cbasis%env)
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
end subroutine test_tblite_overlap_consistency

!> Compare native and libcint evaluators through the common handler interface
subroutine test_integral_handler_consistency(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   type(basis_type) :: basis
   type(libcint_basis_type) :: cbasis
   type(cgto_type) :: cgtoj, cgtoi
   type(native_integral_type) :: native
   type(libcint_integral_type) :: cint
   real(wp) :: sn(25), sc(25), dn(3, 25), dc(3, 25)
   real(wp) :: qn(6, 25), qc(6, 25), vec(3), r2
   real(wp) :: dsn(3, 25), dsc(3, 25)
   real(wp) :: ddjn(3, 3, 25), ddjc(3, 3, 25)
   real(wp) :: ddin(3, 3, 25), ddic(3, 3, 25)
   real(wp) :: dqjn(3, 6, 25), dqjc(3, 6, 25)
   real(wp) :: dqin(3, 6, 25), dqic(3, 6, 25)
   integer :: iat, jat, ish, jsh, n

   call make_comparison_basis(mol, basis, cbasis)
   call native%initialize_integral(mol, basis)
   call cint%initialize_integral(mol, basis)

   ! Select a p-d shell pair and recover the corresponding CGTO descriptions.
   ish = 2
   jsh = 6
   iat = basis%sh2at(ish)
   jat = basis%sh2at(jsh)
   cgtoi = basis%cgto(ish-basis%ish_at(iat), mol%id(iat))
   cgtoj = basis%cgto(jsh-basis%ish_at(jat), mol%id(jat))
   n = basis%nao_sh(ish)*basis%nao_sh(jsh)
   vec = mol%xyz(:, basis%sh2at(ish)) - mol%xyz(:, basis%sh2at(jsh))
   r2 = sum(vec**2)
   call native%multipole_cgto(cgtoj, cgtoi, jsh, ish, r2, vec, basis%intcut, &
      & sn, dn, qn)
   call cint%multipole_cgto(cgtoj, cgtoi, jsh, ish, r2, vec, basis%intcut, &
      & sc, dc, qc)
   call check(error, all(abs(sn(:n) - sc(:n)) < thr), &
      & message="Handler overlap integrals do not match")
   if (allocated(error)) return
   call check(error, all(abs(dn(:, :n) - dc(:, :n)) < thr), &
      & message="Handler dipole integrals do not match")
   if (allocated(error)) return
   call check(error, all(abs(qn(:, :n) - qc(:, :n)) < thr), &
      & message="Handler quadrupole integrals do not match")
   if (allocated(error)) return

   call native%multipole_grad_cgto(cgtoj, cgtoi, jsh, ish, r2, vec, basis%intcut, sn, dn, qn, &
      & dsn, ddjn, dqjn, ddin, dqin)
   call cint%multipole_grad_cgto(cgtoj, cgtoi, jsh, ish, r2, vec, basis%intcut, sc, dc, qc, &
      & dsc, ddjc, dqjc, ddic, dqic)
   call check(error, all(abs(dsn(:, :n) - dsc(:, :n)) < thr), &
      & message="Handler overlap gradients do not match")
   if (allocated(error)) return
   call check(error, all(abs(ddjn(:, :, :n) - ddjc(:, :, :n)) < thr), &
      & message="Handler dipole gradients on center j do not match")
   if (allocated(error)) return
   call check(error, all(abs(ddin(:, :, :n) - ddic(:, :, :n)) < thr), &
      & message="Handler dipole gradients on center i do not match")
   if (allocated(error)) return
   call check(error, all(abs(dqjn(:, :, :n) - dqjc(:, :, :n)) < thr), &
      & message="Handler quadrupole gradients on center j do not match")
   if (allocated(error)) return
   call check(error, all(abs(dqin(:, :, :n) - dqic(:, :, :n)) < thr), &
      & message="Handler quadrupole gradients on center i do not match")
end subroutine test_integral_handler_consistency

!> Compare complete GFN2-xTB energies using native and libcint integral handlers
subroutine test_integral_energy_consistency(error)
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
   call calc_native%set_integral_handler(mol, error, integral_handler_native)
   if (allocated(error)) return
   call new_wavefunction(wfn_native, mol%nat, calc_native%bas%nsh, &
      & calc_native%bas%nao, 1, kt)

   call new_gfn2_calculator(calc_libcint, mol, error)
   if (allocated(error)) return
   call calc_libcint%set_integral_handler(mol, error, integral_handler_libcint)
   if (allocated(error)) return
   call new_wavefunction(wfn_libcint, mol%nat, calc_libcint%bas%nsh, &
      & calc_libcint%bas%nao, 1, kt)

   energy_native = 0.0_wp
   call xtb_singlepoint(ctx_native, mol, calc_native, wfn_native, acc, &
      & energy_native, verbosity=0)
   call check(error, .not. ctx_native%failed(), &
      & message="Native integral handler calculation failed")
   if (allocated(error)) return

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
end subroutine test_integral_energy_consistency

!> Compare native and libcint overlap gradients for every shell pair
subroutine test_tblite_overlap_gradient_consistency(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   type(basis_type) :: basis
   type(libcint_basis_type) :: cbasis
   real(wp) :: ref_s(5, 5), ref_g(3, 5, 5), vec(3), r2
   real(c_double) :: cint_g(5, 5, 3)
   integer :: ish, jsh, iat, jat, isp, jsp, ilsh, jlsh
   integer :: ii, jj, di, dj, iao, jao, ic, stat

   call make_comparison_basis(mol, basis, cbasis)
   do ish = 1, basis%nsh
      iat = basis%sh2at(ish); isp = mol%id(iat)
      ilsh = ish - basis%ish_at(iat)
      ii = basis%iao_sh(ish); di = basis%nao_sh(ish)
      do jsh = 1, basis%nsh
         jat = basis%sh2at(jsh); jsp = mol%id(jat)
         jlsh = jsh - basis%ish_at(jat)
         jj = basis%iao_sh(jsh); dj = basis%nao_sh(jsh)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = sum(vec**2)
         call overlap_grad_cgto(basis%cgto(jlsh, jsp), basis%cgto(ilsh, isp), &
            & r2, vec, basis%intcut, ref_s(1:dj, 1:di), ref_g(:, 1:dj, 1:di))
         ! The reversed order gives libcint output (j,i), with shell i in the
         ! second position differentiated by int1e_ovlpip_sph.
         stat = libcint_eval_overlap_gradient(cint_g, [jsh-1, ish-1], &
            & cbasis%atm, cbasis%bas, cbasis%env)
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
end subroutine test_tblite_overlap_gradient_consistency

!> Compare native and libcint dipole integrals for every shell pair
subroutine test_tblite_dipole_consistency(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   type(basis_type) :: basis
   type(libcint_basis_type) :: cbasis
   real(wp) :: ref_s(5, 5), ref_d(3, 5, 5), vec(3), r2
   real(c_double) :: cint_d(5, 5, 3)
   integer :: ish, jsh, iat, jat, isp, jsp, ilsh, jlsh
   integer :: ii, jj, di, dj, iao, jao, ic, stat

   call make_comparison_basis(mol, basis, cbasis)
   do ish = 1, basis%nsh
      iat = basis%sh2at(ish); isp = mol%id(iat)
      ilsh = ish - basis%ish_at(iat)
      ii = basis%iao_sh(ish); di = basis%nao_sh(ish)
      do jsh = 1, basis%nsh
         jat = basis%sh2at(jsh); jsp = mol%id(jat)
         jlsh = jsh - basis%ish_at(jat)
         jj = basis%iao_sh(jsh); dj = basis%nao_sh(jsh)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = sum(vec**2)
         call dipole_cgto(basis%cgto(jlsh, jsp), basis%cgto(ilsh, isp), &
            & r2, vec, basis%intcut, ref_s(1:dj, 1:di), ref_d(:, 1:dj, 1:di))
         ! Reversing the libcint shells makes origj the tblite ket centre i.
         stat = libcint_eval_dipole(cint_d, [jsh-1, ish-1], &
            & cbasis%atm, cbasis%bas, cbasis%env)
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
end subroutine test_tblite_dipole_consistency

!> Compare native and libcint traceless quadrupole integrals
subroutine test_tblite_quadrupole_consistency(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   type(basis_type) :: basis
   type(libcint_basis_type) :: cbasis
   real(wp) :: ref_s(5, 5), ref_d(3, 5, 5), ref_q(6, 5, 5)
   real(wp) :: vec(3), r2, raw(6), quad(6), trace
   real(c_double) :: cint_q(5, 5, 9)
   integer :: ish, jsh, iat, jat, isp, jsp, ilsh, jlsh
   integer :: ii, jj, di, dj, iao, jao, ic, stat
   ! Map libcint's full Cartesian tensor to tblite's xx, xy, yy, xz, yz, zz order.
   integer, parameter :: qmap(6) = [1, 2, 5, 3, 6, 9]

   call make_comparison_basis(mol, basis, cbasis)
   do ish = 1, basis%nsh
      iat = basis%sh2at(ish); isp = mol%id(iat)
      ilsh = ish - basis%ish_at(iat)
      ii = basis%iao_sh(ish); di = basis%nao_sh(ish)
      do jsh = 1, basis%nsh
         jat = basis%sh2at(jsh); jsp = mol%id(jat)
         jlsh = jsh - basis%ish_at(jat)
         jj = basis%iao_sh(jsh); dj = basis%nao_sh(jsh)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = sum(vec**2)
         call multipole_cgto(basis%cgto(jlsh, jsp), basis%cgto(ilsh, isp), &
            & r2, vec, basis%intcut, ref_s(1:dj, 1:di), ref_d(:, 1:dj, 1:di), &
            & ref_q(:, 1:dj, 1:di))
         stat = libcint_eval_quadrupole(cint_q, [jsh-1, ish-1], &
            & cbasis%atm, cbasis%bas, cbasis%env)
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
end subroutine test_tblite_quadrupole_consistency

!> Compare native and libcint dipole gradients for both Gaussian centers
subroutine test_tblite_dipole_gradient_consistency(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   type(basis_type) :: basis
   type(libcint_basis_type) :: cbasis
   real(wp) :: ref_s(5, 5), ref_d(3, 5, 5), ref_q(6, 5, 5), ref_g(3, 5, 5)
   real(wp) :: ref_dj(3, 3, 5, 5), ref_di(3, 3, 5, 5)
   real(wp) :: ref_qj(3, 6, 5, 5), ref_qi(3, 6, 5, 5), vec(3), r2
   real(c_double) :: cint_j(5, 5, 3, 3), cint_i(5, 5, 3, 3)
   real(c_double) :: swap_i(5, 5, 3, 3), swap_j(5, 5, 3, 3)
   integer :: ish, jsh, iat, jat, isp, jsp, ilsh, jlsh
   integer :: di, dj, iao, jao, im, ider, stat

   call make_comparison_basis(mol, basis, cbasis)
   do ish = 1, basis%nsh
      iat = basis%sh2at(ish); isp = mol%id(iat)
      ilsh = ish - basis%ish_at(iat); di = basis%nao_sh(ish)
      do jsh = 1, basis%nsh
         jat = basis%sh2at(jsh); jsp = mol%id(jat)
         jlsh = jsh - basis%ish_at(jat); dj = basis%nao_sh(jsh)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat); r2 = sum(vec**2)
         call multipole_grad_cgto(basis%cgto(jlsh, jsp), basis%cgto(ilsh, isp), &
            & r2, vec, basis%intcut, ref_s(1:dj, 1:di), ref_d(:, 1:dj, 1:di), &
            & ref_q(:, 1:dj, 1:di), ref_g(:, 1:dj, 1:di), &
            & ref_dj(:, :, 1:dj, 1:di), ref_qj(:, :, 1:dj, 1:di), &
            & ref_di(:, :, 1:dj, 1:di), ref_qi(:, :, 1:dj, 1:di))
         ! Evaluate both shell orders because the moment origin follows the ket shell.
         stat = libcint_eval_dipole_gradient(cint_j, cint_i, [jsh-1, ish-1], &
            & cbasis%atm, cbasis%bas, cbasis%env)
         call check(error, stat >= 0, &
            & message="Libcint dipole gradient evaluation failed")
         if (allocated(error)) return
         stat = libcint_eval_dipole_gradient(swap_i, swap_j, [ish-1, jsh-1], &
            & cbasis%atm, cbasis%bas, cbasis%env)
         call check(error, stat >= 0, &
            & message="Libcint reversed dipole gradient evaluation failed")
         if (allocated(error)) return
         do ider = 1, 3; do im = 1, 3; do iao = 1, di; do jao = 1, dj
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
         end do; end do; end do; end do
      end do
   end do
end subroutine test_tblite_dipole_gradient_consistency

!> Compare native and libcint quadrupole gradients for both Gaussian centers
subroutine test_tblite_quadrupole_gradient_consistency(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   type(basis_type) :: basis
   type(libcint_basis_type) :: cbasis
   real(wp) :: ref_s(5, 5), ref_d(3, 5, 5), ref_q(6, 5, 5), ref_g(3, 5, 5)
   real(wp) :: ref_dj(3, 3, 5, 5), ref_di(3, 3, 5, 5)
   real(wp) :: ref_qj(3, 6, 5, 5), ref_qi(3, 6, 5, 5), raw(6), quad(6), trace
   real(wp) :: vec(3), r2
   real(c_double) :: cint_j(5, 5, 9, 3), cint_i(5, 5, 9, 3)
   real(c_double) :: swap_i(5, 5, 9, 3), swap_j(5, 5, 9, 3)
   integer :: ish, jsh, iat, jat, isp, jsp, ilsh, jlsh
   integer :: di, dj, iao, jao, im, ider, stat, center
   ! Map libcint's full Cartesian tensor to tblite's xx, xy, yy, xz, yz, zz order.
   integer, parameter :: qmap(6) = [1, 2, 5, 3, 6, 9]

   call make_comparison_basis(mol, basis, cbasis)
   do ish = 1, basis%nsh
      iat = basis%sh2at(ish); isp = mol%id(iat)
      ilsh = ish - basis%ish_at(iat); di = basis%nao_sh(ish)
      do jsh = 1, basis%nsh
         jat = basis%sh2at(jsh); jsp = mol%id(jat)
         jlsh = jsh - basis%ish_at(jat); dj = basis%nao_sh(jsh)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat); r2 = sum(vec**2)
         call multipole_grad_cgto(basis%cgto(jlsh, jsp), basis%cgto(ilsh, isp), &
            & r2, vec, basis%intcut, ref_s(1:dj, 1:di), ref_d(:, 1:dj, 1:di), &
            & ref_q(:, 1:dj, 1:di), ref_g(:, 1:dj, 1:di), &
            & ref_dj(:, :, 1:dj, 1:di), ref_qj(:, :, 1:dj, 1:di), &
            & ref_di(:, :, 1:dj, 1:di), ref_qi(:, :, 1:dj, 1:di))
         ! Evaluate both shell orders because the moment origin follows the ket shell.
         stat = libcint_eval_quadrupole_gradient(cint_j, cint_i, [jsh-1, ish-1], &
            & cbasis%atm, cbasis%bas, cbasis%env)
         call check(error, stat >= 0, &
            & message="Libcint quadrupole gradient evaluation failed")
         if (allocated(error)) return
         stat = libcint_eval_quadrupole_gradient(swap_i, swap_j, [ish-1, jsh-1], &
            & cbasis%atm, cbasis%bas, cbasis%env)
         call check(error, stat >= 0, &
            & message="Libcint reversed quadrupole gradient evaluation failed")
         if (allocated(error)) return
         do center = 1, 4; do ider = 1, 3; do iao = 1, di; do jao = 1, dj
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
                     & thr=thr, message="Libcint bra-center quadrupole gradient does not match")
               else if (center == 2) then
                  call check(error, quad(im), ref_qi(ider, im, jao, iao), &
                     & thr=thr, message="Libcint ket-center quadrupole gradient does not match")
               else if (center == 3) then
                  call check(error, quad(im), ref_qj(ider, im, jao, iao), &
                     & thr=thr, message="Libcint reversed ket-center quadrupole gradient does not match")
               else
                  call check(error, quad(im), -ref_qj(ider, im, jao, iao), &
                     & thr=thr, message="Libcint reversed bra-center quadrupole gradient does not match")
               end if
               if (allocated(error)) return
            end do
         end do; end do; end do; end do
      end do
   end do
end subroutine test_tblite_quadrupole_gradient_consistency

#else
!> Check that libcint support is reported as disabled when unavailable
subroutine test_disabled(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check(error, .not. tblite_use_libcint, &
      & message="Build unexpectedly reports libcint support")
end subroutine test_disabled
#endif

end module test_integral_libcint
