! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later

module test_integral_libcint
   use, intrinsic :: iso_c_binding, only : c_double
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check
   use mctc_io, only : structure_type, new
   use tblite_basis_type, only : basis_type, cgto_type, new_basis, new_cgto, get_cutoff
   use tblite_features, only : tblite_use_libcint
#if TBLITE_HAS_LIBCINT
   use tblite_integral_libcint
   use tblite_integral_dipole, only : dipole_cgto
   use tblite_integral_multipole, only : multipole_cgto, multipole_grad_cgto
   use tblite_integral_overlap, only : get_overlap, overlap_grad_cgto, &
      & get_cartesian_exponents
#endif
   implicit none
   private

   public :: collect_integral_libcint

   real(wp), parameter :: thr = 1.0e-10_wp

contains

subroutine collect_integral_libcint(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

#if TBLITE_HAS_LIBCINT
   testsuite = [ &
      new_unittest("tblite-overlap-consistency", test_tblite_overlap_consistency), &
      new_unittest("tblite-overlap-gradient-consistency", &
         & test_tblite_overlap_gradient_consistency), &
      new_unittest("cca-cartesian-ordering", test_cca_cartesian_ordering), &
      new_unittest("tblite-dipole-consistency", test_tblite_dipole_consistency), &
      new_unittest("tblite-quadrupole-consistency", test_tblite_quadrupole_consistency), &
      new_unittest("tblite-dipole-gradient-consistency", &
         & test_tblite_dipole_gradient_consistency), &
      new_unittest("tblite-quadrupole-gradient-consistency", &
         & test_tblite_quadrupole_gradient_consistency) &
      ]
#else
   testsuite = [ &
      new_unittest("disabled", test_disabled) &
      ]
#endif
end subroutine collect_integral_libcint

#if TBLITE_HAS_LIBCINT
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

subroutine test_tblite_overlap_consistency(error)
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(basis_type) :: basis
   type(libcint_basis_type) :: cbasis
   real(wp), allocatable :: ref(:, :), cint_overlap(:, :)
   real(wp) :: trans(3, 1), cutoff
   real(c_double) :: block(5, 5)
   integer :: ish, jsh, ii, jj, di, dj, iao, jao, stat

   call make_comparison_basis(mol, basis, cbasis)

   call check(error, size(cbasis%atm, 2), mol%nat)
   if (allocated(error)) return
   call check(error, size(cbasis%bas, 2), basis%nsh)
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
         call check(error, stat >= 0)
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
         call check(error, cint_overlap(jj, ii), ref(jj, ii), thr=thr)
         if (allocated(error)) return
      end do
   end do
end subroutine test_tblite_overlap_consistency

subroutine test_tblite_overlap_gradient_consistency(error)
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
         call check(error, stat >= 0)
         if (allocated(error)) return
         do ic = 1, 3
            do iao = 1, di
               do jao = 1, dj
                  call check(error, real(cint_g(jao, iao, ic), wp), &
                     & ref_g(ic, jao, iao), thr=thr)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end do
end subroutine test_tblite_overlap_gradient_consistency

subroutine test_cca_cartesian_ordering(error)
   type(error_type), allocatable, intent(out) :: error
   integer :: exponents(3, 15)

   call get_cartesian_exponents(1, exponents(:, 1:3))
   call check(error, all(exponents(:, 1:3) == reshape([ &
      & 1,0,0, 0,1,0, 0,0,1], [3, 3])))
   if (allocated(error)) return
   call get_cartesian_exponents(2, exponents(:, 1:6))
   call check(error, all(exponents(:, 1:6) == reshape([ &
      & 2,0,0, 1,1,0, 1,0,1, 0,2,0, 0,1,1, 0,0,2], [3, 6])))
   if (allocated(error)) return
   call get_cartesian_exponents(3, exponents(:, 1:10))
   call check(error, all(exponents(:, 1:10) == reshape([ &
      & 3,0,0, 2,1,0, 2,0,1, 1,2,0, 1,1,1, 1,0,2, &
      & 0,3,0, 0,2,1, 0,1,2, 0,0,3], [3, 10])))
   if (allocated(error)) return
   call get_cartesian_exponents(4, exponents)
   call check(error, all(exponents == reshape([ &
      & 4,0,0, 3,1,0, 3,0,1, 2,2,0, 2,1,1, 2,0,2, &
      & 1,3,0, 1,2,1, 1,1,2, 1,0,3, 0,4,0, 0,3,1, &
      & 0,2,2, 0,1,3, 0,0,4], [3, 15])))
end subroutine test_cca_cartesian_ordering

subroutine test_tblite_dipole_consistency(error)
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
         call check(error, stat >= 0)
         if (allocated(error)) return
         do ic = 1, 3
            do iao = 1, di
               do jao = 1, dj
                  call check(error, real(cint_d(jao, iao, ic), wp), &
                     & ref_d(ic, jao, iao), thr=thr)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end do
end subroutine test_tblite_dipole_consistency

subroutine test_tblite_quadrupole_consistency(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   type(basis_type) :: basis
   type(libcint_basis_type) :: cbasis
   real(wp) :: ref_s(5, 5), ref_d(3, 5, 5), ref_q(6, 5, 5)
   real(wp) :: vec(3), r2, raw(6), quad(6), trace
   real(c_double) :: cint_q(5, 5, 9)
   integer :: ish, jsh, iat, jat, isp, jsp, ilsh, jlsh
   integer :: ii, jj, di, dj, iao, jao, ic, stat
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
         call check(error, stat >= 0)
         if (allocated(error)) return
         do iao = 1, di
            do jao = 1, dj
               do ic = 1, 6
                  raw(ic) = real(cint_q(jao, iao, qmap(ic)), wp)
               end do
               trace = 0.5_wp*(raw(1) + raw(3) + raw(6))
               quad = 1.5_wp*raw
               quad([1, 3, 6]) = quad([1, 3, 6]) - trace
               do ic = 1, 6
                  call check(error, quad(ic), ref_q(ic, jao, iao), thr=thr)
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end do
end subroutine test_tblite_quadrupole_consistency

subroutine test_tblite_dipole_gradient_consistency(error)
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
         stat = libcint_eval_dipole_gradient(cint_j, cint_i, [jsh-1, ish-1], &
            & cbasis%atm, cbasis%bas, cbasis%env)
         call check(error, stat >= 0)
         if (allocated(error)) return
         stat = libcint_eval_dipole_gradient(swap_i, swap_j, [ish-1, jsh-1], &
            & cbasis%atm, cbasis%bas, cbasis%env)
         call check(error, stat >= 0)
         if (allocated(error)) return
         do ider = 1, 3; do im = 1, 3; do iao = 1, di; do jao = 1, dj
            ! For origj=i, d/dRi=d/dvec and d/dRj=-d/dvec.
            call check(error, real(cint_j(jao, iao, im, ider), wp), &
               & -ref_di(ider, im, jao, iao), thr=thr)
            if (allocated(error)) return
            call check(error, real(cint_i(jao, iao, im, ider), wp), &
               & ref_di(ider, im, jao, iao), thr=thr)
            if (allocated(error)) return
            ! For origj=j the reversed block provides tblite's j-centred operator.
            call check(error, real(swap_i(iao, jao, im, ider), wp), &
               & ref_dj(ider, im, jao, iao), thr=thr)
            if (allocated(error)) return
            call check(error, real(swap_j(iao, jao, im, ider), wp), &
               & -ref_dj(ider, im, jao, iao), thr=thr)
            if (allocated(error)) return
         end do; end do; end do; end do
      end do
   end do
end subroutine test_tblite_dipole_gradient_consistency

subroutine test_tblite_quadrupole_gradient_consistency(error)
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
         stat = libcint_eval_quadrupole_gradient(cint_j, cint_i, [jsh-1, ish-1], &
            & cbasis%atm, cbasis%bas, cbasis%env)
         call check(error, stat >= 0)
         if (allocated(error)) return
         stat = libcint_eval_quadrupole_gradient(swap_i, swap_j, [ish-1, jsh-1], &
            & cbasis%atm, cbasis%bas, cbasis%env)
         call check(error, stat >= 0)
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
            trace = 0.5_wp*(raw(1) + raw(3) + raw(6))
            quad = 1.5_wp*raw
            quad([1, 3, 6]) = quad([1, 3, 6]) - trace
            do im = 1, 6
               if (center == 1) then
                  call check(error, quad(im), -ref_qi(ider, im, jao, iao), thr=thr)
               else if (center == 2) then
                  call check(error, quad(im), ref_qi(ider, im, jao, iao), thr=thr)
               else if (center == 3) then
                  call check(error, quad(im), ref_qj(ider, im, jao, iao), thr=thr)
               else
                  call check(error, quad(im), -ref_qj(ider, im, jao, iao), thr=thr)
               end if
               if (allocated(error)) return
            end do
         end do; end do; end do; end do
      end do
   end do
end subroutine test_tblite_quadrupole_gradient_consistency

#else
subroutine test_disabled(error)
   type(error_type), allocatable, intent(out) :: error

   call check(error, .not. tblite_use_libcint)
end subroutine test_disabled
#endif

end module test_integral_libcint
