! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later

module test_solvation_cosmo
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, test_failed
   use tblite_solvation_cosmo, only : get_gaussian_multipole_kernels
   implicit none
   private

   public :: collect_solvation_cosmo

contains

subroutine collect_solvation_cosmo(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("gaussian-multipole-finite-difference", &
         test_gaussian_multipole_finite_difference), &
      new_unittest("gaussian-multipole-point-limit", &
         test_gaussian_multipole_point_limit) &
      ]
end subroutine collect_solvation_cosmo

subroutine test_gaussian_multipole_finite_difference(error)
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: xis(3) = [0.45_wp, 1.20_wp, 2.75_wp]
   real(wp), parameter :: vecs(3, 3) = reshape([ &
      1.30_wp, -0.70_wp, 2.10_wp, &
     -1.80_wp,  1.10_wp, 0.60_wp, &
      0.55_wp,  1.65_wp, -1.25_wp], [3, 3])
   real(wp), parameter :: qmat(3, 3) = reshape([ &
       0.70_wp,  0.20_wp, -0.10_wp, &
       0.20_wp, -0.40_wp,  0.30_wp, &
      -0.10_wp,  0.30_wp, -0.30_wp], [3, 3])
   real(wp), parameter :: qvec(6) = [ &
      qmat(1, 1), qmat(1, 2), qmat(2, 2), &
      qmat(1, 3), qmat(2, 3), qmat(3, 3)]
   real(wp), parameter :: h = 5.0e-4_wp
   real(wp), parameter :: dipole_tol = 2.0e-9_wp
   real(wp), parameter :: quadrupole_tol = 2.0e-7_wp

   real(wp) :: center(3), shifted(3), c0, c1(3), c2(6)
   real(wp) :: dnum(3), hnum(3, 3), f0, qnum, qimpl
   integer :: icase, ia, ib

   do icase = 1, size(xis)
      ! Put the surface Gaussian at the origin, so s_i-R_A=vecs(:,icase).
      center = -vecs(:, icase)
      call get_gaussian_multipole_kernels(vecs(:, icase), xis(icase), &
         c0, c1, c2)
      f0 = gaussian_monopole(center, xis(icase))

      do ia = 1, 3
         shifted = center
         shifted(ia) = center(ia) - 2.0_wp*h
         dnum(ia) = gaussian_monopole(shifted, xis(icase))
         shifted(ia) = center(ia) - h
         dnum(ia) = dnum(ia) - 8.0_wp*gaussian_monopole(shifted, xis(icase))
         shifted(ia) = center(ia) + h
         dnum(ia) = dnum(ia) + 8.0_wp*gaussian_monopole(shifted, xis(icase))
         shifted(ia) = center(ia) + 2.0_wp*h
         dnum(ia) = (dnum(ia) - gaussian_monopole(shifted, xis(icase))) &
            & /(12.0_wp*h)

         if (abs(c1(ia) - dnum(ia)) > dipole_tol) then
            call test_failed(error, &
               "Gaussian dipole kernel disagrees with finite differences")
            return
         end if

         shifted = center
         shifted(ia) = center(ia) - 2.0_wp*h
         hnum(ia, ia) = -gaussian_monopole(shifted, xis(icase))
         shifted(ia) = center(ia) - h
         hnum(ia, ia) = hnum(ia, ia) &
            & + 16.0_wp*gaussian_monopole(shifted, xis(icase))
         shifted(ia) = center(ia) + h
         hnum(ia, ia) = hnum(ia, ia) &
            & + 16.0_wp*gaussian_monopole(shifted, xis(icase))
         shifted(ia) = center(ia) + 2.0_wp*h
         hnum(ia, ia) = (hnum(ia, ia) &
            & - gaussian_monopole(shifted, xis(icase)) - 30.0_wp*f0) &
            & /(12.0_wp*h*h)
      end do

      do ia = 1, 3
         do ib = ia + 1, 3
            shifted = center
            shifted(ia) = center(ia) + h
            shifted(ib) = center(ib) + h
            hnum(ia, ib) = gaussian_monopole(shifted, xis(icase))
            shifted(ib) = center(ib) - h
            hnum(ia, ib) = hnum(ia, ib) &
               & - gaussian_monopole(shifted, xis(icase))
            shifted(ia) = center(ia) - h
            hnum(ia, ib) = hnum(ia, ib) &
               & + gaussian_monopole(shifted, xis(icase))
            shifted(ib) = center(ib) + h
            hnum(ia, ib) = hnum(ia, ib) &
               & - gaussian_monopole(shifted, xis(icase))
            hnum(ia, ib) = hnum(ia, ib)/(4.0_wp*h*h)
            hnum(ib, ia) = hnum(ia, ib)
         end do
      end do

      qnum = sum(qmat*hnum)/3.0_wp
      qimpl = dot_product(qvec, c2)
      if (abs(qimpl - qnum) > quadrupole_tol) then
         call test_failed(error, &
            "Gaussian quadrupole kernel disagrees with finite differences")
         return
      end if
   end do
end subroutine test_gaussian_multipole_finite_difference

subroutine test_gaussian_multipole_point_limit(error)
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: vec(3) = [1.25_wp, -0.85_wp, 1.70_wp]
   real(wp), parameter :: xi = 30.0_wp
   real(wp), parameter :: tol = 5.0e-13_wp
   real(wp) :: c0, c1(3), c2(6), r2, r1, ref2(6)

   call get_gaussian_multipole_kernels(vec, xi, c0, c1, c2)
   r2 = sum(vec**2)
   r1 = sqrt(r2)
   ref2 = [vec(1)*vec(1), 2.0_wp*vec(1)*vec(2), &
      & vec(2)*vec(2), 2.0_wp*vec(1)*vec(3), &
      & 2.0_wp*vec(2)*vec(3), vec(3)*vec(3)]/(r1*r2*r2)

   if (abs(c0 - 1.0_wp/r1) > tol &
      & .or. maxval(abs(c1 - vec/(r1*r2))) > tol &
      & .or. maxval(abs(c2 - ref2)) > tol) then
      call test_failed(error, &
         "Gaussian multipole kernels do not recover the point limit")
   end if
end subroutine test_gaussian_multipole_point_limit

pure function gaussian_monopole(center, xi) result(value)
   real(wp), intent(in) :: center(3)
   real(wp), intent(in) :: xi
   real(wp) :: value, distance

   distance = sqrt(sum(center**2))
   value = erf(xi*distance)/distance
end function gaussian_monopole

end module test_solvation_cosmo
