! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later

!> Test for kernel gradient computation
module test_kernel_gradient
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check
   use mctc_io, only : structure_type
   use tblite_solvation_alpb, only : alpb_solvation
   use tblite_test_utils, only : make_mol
   implicit none
   private

   public :: collect_kernel_gradient

contains

!> Collect all exported unit tests
subroutine collect_kernel_gradient(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("kernel-still", test_kernel_still_spatial) &
      ]
end subroutine collect_kernel_gradient


!> Test spatial kernel gradient for Still kernel
subroutine test_kernel_still_spatial(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: xyz(:, :), brad(:), brdr(:, :, :)
   real(wp), allocatable :: kernel_grad_spatial(:, :, :), kernel_grad_born(:, :)
   real(wp), allocatable :: num_grad(:, :, :)
   real(wp) :: keps, delta, kernel_ref, kernel_pert, max_error
   integer :: nat, i, j, k

   ! Simple test system
   nat = 3
   allocate(xyz(3, nat))
   allocate(brad(nat), brdr(3, nat, nat))
   allocate(kernel_grad_spatial(3, nat, nat))
   allocate(kernel_grad_born(nat, nat))
   allocate(num_grad(3, nat, nat))

   ! Geometry
   xyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
   xyz(:, 2) = [3.0_wp, 0.0_wp, 0.0_wp]
   xyz(:, 3) = [0.0_wp, 3.0_wp, 0.0_wp]

   ! Born radii
   brad(1) = 2.5_wp
   brad(2) = 2.0_wp
   brad(3) = 2.0_wp

   ! No Born radii position dependence
   brdr(:, :, :) = 0.0_wp

   ! Dielectric screening
   keps = 0.5_wp

   ! Compute analytical kernel gradients
   call compute_kernel_deriv_still(nat, xyz, keps, brad, brdr, &
      & kernel_grad_spatial, kernel_grad_born)

   ! Numerical gradient
   delta = 1.0e-6_wp
   num_grad(:, :, :) = 0.0_wp
   
   do i = 1, nat
      do k = 1, 3
         xyz(k, i) = xyz(k, i) + delta
         kernel_pert = evaluate_kernel_sum(nat, xyz, brad, keps)
         
         xyz(k, i) = xyz(k, i) - 2.0_wp * delta
         kernel_ref = evaluate_kernel_sum(nat, xyz, brad, keps)
         
         xyz(k, i) = xyz(k, i) + delta
         
         do j = 1, nat
            if (i /= j) then
               if (i < j) then
                  num_grad(k, i, j) = (kernel_pert - kernel_ref) / (2.0_wp * delta)
               end if
            end if
         end do
      end do
   end do

   ! Symmetrize
   do i = 1, nat
      do j = i + 1, nat
         num_grad(:, j, i) = -num_grad(:, i, j)
      end do
   end do

   ! Check agreement
   max_error = 0.0_wp
   do i = 1, nat
      do j = i + 1, nat
         do k = 1, 3
            max_error = max(max_error, abs(kernel_grad_spatial(k, i, j) - num_grad(k, i, j)))
         end do
      end do
   end do

   call check(error, max_error, 0.0_wp, thr=1.0e-5_wp)

end subroutine test_kernel_still_spatial


!> Evaluate kernel sum
function evaluate_kernel_sum(n, coords, radii, kappa) result(kernel_sum)
   integer, intent(in) :: n
   real(wp), intent(in) :: coords(:, :)
   real(wp), intent(in) :: radii(:)
   real(wp), intent(in) :: kappa
   real(wp) :: kernel_sum
   
   integer :: ii, jj
   real(wp) :: vec(3), r1, r2, aa, dd, expd, fgb2, dfgb
   real(wp), parameter :: a4 = 0.25_wp
   
   kernel_sum = 0.0_wp
   
   do ii = 1, n
      do jj = 1, ii - 1
         vec(:) = coords(:, ii) - coords(:, jj)
         r1 = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
         r2 = r1 * r1
         
         aa = radii(ii) * radii(jj)
         dd = a4 * r2 / aa
         expd = exp(-dd)
         fgb2 = r2 + aa * expd
         dfgb = 1.0_wp / sqrt(fgb2)
         
         kernel_sum = kernel_sum + kappa * dfgb
      end do
   end do
   
   do ii = 1, n
      kernel_sum = kernel_sum + 0.5_wp * kappa / radii(ii)
   end do
   
end function evaluate_kernel_sum


!> Compute analytical kernel derivatives
subroutine compute_kernel_deriv_still(n, coords, kappa, radii, radii_deriv, &
      & grad_spatial, grad_born)
   integer, intent(in) :: n
   real(wp), intent(in) :: coords(:, :)
   real(wp), intent(in) :: kappa
   real(wp), intent(in) :: radii(:)
   real(wp), intent(in) :: radii_deriv(:, :, :)
   real(wp), intent(out) :: grad_spatial(:, :, :)
   real(wp), intent(out) :: grad_born(:, :)

   integer :: ii, jj
   real(wp), parameter :: a4 = 0.25_wp
   real(wp) :: aa, r2, fgb2, dd, expd, dfgb, dfgb2, dfgb3
   real(wp) :: ap, bp, r1, vec(3)

   grad_spatial(:, :, :) = 0.0_wp
   grad_born(:, :) = 0.0_wp

   do ii = 1, n
      do jj = 1, ii - 1
         vec(:) = coords(:, ii) - coords(:, jj)
         r1 = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
         r2 = r1 * r1

         aa = radii(ii) * radii(jj)
         dd = a4 * r2 / aa
         expd = exp(-dd)
         fgb2 = r2 + aa * expd
         dfgb2 = 1._wp / fgb2
         dfgb = sqrt(dfgb2)
         dfgb3 = dfgb2 * dfgb * kappa

         ap = (1._wp - a4 * expd) * dfgb3
         
         grad_spatial(:, ii, jj) = ap * vec
         grad_spatial(:, jj, ii) = -ap * vec

         bp = -0.5_wp * expd * (1._wp + dd) * dfgb3
         
         grad_born(ii, jj) = bp * radii(jj)
         grad_born(jj, ii) = bp * radii(ii)
      end do
   end do

end subroutine compute_kernel_deriv_still


end module test_kernel_gradient
