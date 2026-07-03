! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later

module test_integral_libcint
   use, intrinsic :: iso_c_binding, only : c_double, c_int
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check
   use tblite_features, only : tblite_use_libcint
#if TBLITE_HAS_LIBCINT
   use tblite_integral_libcint
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
      new_unittest("overlap-kinetic-nuclear", test_one_electron), &
      new_unittest("eri-shell", test_eri_shell), &
      new_unittest("cartesian-spheric-dimensions", test_cart_sph_dimensions) &
      ]
#else
   testsuite = [ &
      new_unittest("disabled", test_disabled) &
      ]
#endif
end subroutine collect_integral_libcint

#if TBLITE_HAS_LIBCINT
subroutine make_h2_basis(atm, bas, env)
   integer(c_int), allocatable, intent(out) :: atm(:, :)
   integer(c_int), allocatable, intent(out) :: bas(:, :)
   real(c_double), allocatable, intent(out) :: env(:)

   integer :: off

   allocate(atm(ATM_SLOTS, 2), source=0_c_int)
   allocate(bas(BAS_SLOTS, 2), source=0_c_int)
   allocate(env(64), source=0.0_c_double)

   off = PTR_ENV_START
   atm(CHARGE_OF, 1) = 1_c_int
   atm(PTR_COORD, 1) = int(off, c_int)
   env(off+1) = 0.0_c_double
   env(off+2) = 0.0_c_double
   env(off+3) = -0.7_c_double
   off = off + 3

   atm(CHARGE_OF, 2) = 1_c_int
   atm(PTR_COORD, 2) = int(off, c_int)
   env(off+1) = 0.0_c_double
   env(off+2) = 0.0_c_double
   env(off+3) = 0.7_c_double
   off = off + 3

   bas(ATOM_OF, 1) = 0_c_int
   bas(ANG_OF, 1) = 0_c_int
   bas(NPRIM_OF, 1) = 1_c_int
   bas(NCTR_OF, 1) = 1_c_int
   bas(PTR_EXP, 1) = int(off, c_int)
   env(off+1) = 1.0_c_double
   off = off + 1
   bas(PTR_COEFF, 1) = int(off, c_int)
   env(off+1) = real(libcint_gto_norm(0, 1.0_wp), c_double)
   off = off + 1

   bas(ATOM_OF, 2) = 1_c_int
   bas(ANG_OF, 2) = bas(ANG_OF, 1)
   bas(NPRIM_OF, 2) = bas(NPRIM_OF, 1)
   bas(NCTR_OF, 2) = bas(NCTR_OF, 1)
   bas(PTR_EXP, 2) = bas(PTR_EXP, 1)
   bas(PTR_COEFF, 2) = bas(PTR_COEFF, 1)
end subroutine make_h2_basis

subroutine test_one_electron(error)
   type(error_type), allocatable, intent(out) :: error

   integer(c_int), allocatable :: atm(:, :), bas(:, :)
   real(c_double), allocatable :: env(:)
   real(c_double) :: buf(1, 1)
   integer :: stat

   call make_h2_basis(atm, bas, env)

   stat = libcint_eval_1e(LIBCINT_1E_OVERLAP, LIBCINT_SPHERICAL, buf, [0, 0], atm, bas, env)
   call check(error, stat >= 0)
   if (allocated(error)) return
   call check(error, real(buf(1, 1), wp), 1.0_wp, thr=thr)
   if (allocated(error)) return

   stat = libcint_eval_1e(LIBCINT_1E_OVERLAP, LIBCINT_SPHERICAL, buf, [0, 1], atm, bas, env)
   call check(error, stat >= 0)
   if (allocated(error)) return
   call check(error, all(abs(buf) < huge(1.0_c_double)))
   if (allocated(error)) return
   call check(error, buf(1, 1) > 0.0_c_double)
   if (allocated(error)) return

   stat = libcint_eval_1e(LIBCINT_1E_KINETIC, LIBCINT_SPHERICAL, buf, [0, 1], atm, bas, env)
   call check(error, stat >= 0)
   if (allocated(error)) return
   call check(error, all(abs(buf) < huge(1.0_c_double)))
   if (allocated(error)) return

   stat = libcint_eval_1e(LIBCINT_1E_NUCLEAR, LIBCINT_SPHERICAL, buf, [0, 1], atm, bas, env)
   call check(error, stat >= 0)
   if (allocated(error)) return
   call check(error, all(abs(buf) < huge(1.0_c_double)))
end subroutine test_one_electron

subroutine test_eri_shell(error)
   type(error_type), allocatable, intent(out) :: error

   integer(c_int), allocatable :: atm(:, :), bas(:, :)
   real(c_double), allocatable :: env(:)
   real(c_double) :: eri_0101(1, 1, 1, 1), eri_1010(1, 1, 1, 1)
   integer :: stat

   call make_h2_basis(atm, bas, env)

   stat = libcint_eval_eri(LIBCINT_SPHERICAL, eri_0101, [0, 1, 0, 1], atm, bas, env)
   call check(error, stat >= 0)
   if (allocated(error)) return
   call check(error, all(abs(eri_0101) < huge(1.0_c_double)))
   if (allocated(error)) return
   call check(error, eri_0101(1, 1, 1, 1) > 0.0_c_double)
   if (allocated(error)) return

   stat = libcint_eval_eri(LIBCINT_SPHERICAL, eri_1010, [1, 0, 1, 0], atm, bas, env)
   call check(error, stat >= 0)
   if (allocated(error)) return
   call check(error, real(eri_0101(1, 1, 1, 1) - eri_1010(1, 1, 1, 1), wp), 0.0_wp, thr=thr)
end subroutine test_eri_shell

subroutine test_cart_sph_dimensions(error)
   type(error_type), allocatable, intent(out) :: error

   integer(c_int), allocatable :: atm(:, :), bas(:, :)
   real(c_double), allocatable :: env(:)

   call make_h2_basis(atm, bas, env)

   call check(error, libcint_cgto_cart(0, bas), 1)
   if (allocated(error)) return
   call check(error, libcint_cgto_spheric(0, bas), 1)
   if (allocated(error)) return
   call check(error, libcint_tot_cgto_cart(bas), 2)
   if (allocated(error)) return
   call check(error, libcint_tot_cgto_spheric(bas), 2)
end subroutine test_cart_sph_dimensions

#else
subroutine test_disabled(error)
   type(error_type), allocatable, intent(out) :: error

   call check(error, .not. tblite_use_libcint)
end subroutine test_disabled
#endif

end module test_integral_libcint
