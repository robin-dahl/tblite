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

module test_solvation_kernel
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_solvation_born, only : born_integrator, new_born_integrator
   use tblite_solvation_data, only : get_vdw_rad_cosmo, get_vdw_rad_d3
   use tblite_solvation_kernel, only : kernel_type, new_kernel, kernel_enum, compute_kernel_dkdr_ij, compute_kernel_d2kdr2_ij, compute_kernel_d3kdr3_ij
   implicit none
   private

   public :: collect_solvation_kernel

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains

!> Collect all exported unit tests
subroutine collect_solvation_kernel(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("kernel-gradient-still", test_kernel_gradient_still), &
      new_unittest("kernel-hessian-still", test_kernel_hessian_still), &
      new_unittest("kernel-third-still", test_kernel_third_still), &
      new_unittest("kernel-gradient-p16", test_kernel_gradient_p16), &
      new_unittest("kernel-hessian-p16", test_kernel_hessian_p16), &
      new_unittest("kernel-third-p16", test_kernel_third_p16) &
      ]

end subroutine collect_solvation_kernel


!> Test Still kernel gradient against numerical derivative
subroutine test_kernel_gradient_still(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 0.5_wp
   real(wp), parameter :: qat(*) = [&
      & -2.11018727757438E-1_wp, -6.04389222813257E-2_wp, -1.90601159250311E-1_wp, &
      &  1.49237694872530E-1_wp,  1.35835820853652E-1_wp,  1.27732431639016E-1_wp, &
      &  1.78559147201780E-1_wp,  1.42324484825195E-1_wp,  1.92106458233743E-1_wp, &
      &  1.45841758574287E-1_wp,  1.56456166394024E-1_wp,  1.59746890863949E-1_wp, &
      & -2.70765876809499E-1_wp, -3.27435355522312E-1_wp, -4.70046325670683E-2_wp, &
      &  1.10838969762146E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_kernel_numg(error, mol, kernel_enum%still, keps, qat)

end subroutine test_kernel_gradient_still

!> Test Still kernel Hessian against numerical derivative
subroutine test_kernel_hessian_still(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 0.5_wp
   real(wp), parameter :: qat(*) = [&
      & -2.11018727757438E-1_wp, -6.04389222813257E-2_wp, -1.90601159250311E-1_wp, &
      &  1.49237694872530E-1_wp,  1.35835820853652E-1_wp,  1.27732431639016E-1_wp, &
      &  1.78559147201780E-1_wp,  1.42324484825195E-1_wp,  1.92106458233743E-1_wp, &
      &  1.45841758574287E-1_wp,  1.56456166394024E-1_wp,  1.59746890863949E-1_wp, &
      & -2.70765876809499E-1_wp, -3.27435355522312E-1_wp, -4.70046325670683E-2_wp, &
      &  1.10838969762146E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_kernel_numh(error, mol, kernel_enum%still, keps, qat)

end subroutine test_kernel_hessian_still


!> Test Still kernel third derivative against numerical derivative
subroutine test_kernel_third_still(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), parameter :: keps = 0.5_wp
   real(wp), parameter :: qat(*) = [&
      & -2.11018727757438E-1_wp, -6.04389222813257E-2_wp, -1.90601159250311E-1_wp, &
      &  1.49237694872530E-1_wp,  1.35835820853652E-1_wp,  1.27732431639016E-1_wp, &
      &  1.78559147201780E-1_wp,  1.42324484825195E-1_wp,  1.92106458233743E-1_wp, &
      &  1.45841758574287E-1_wp,  1.56456166394024E-1_wp,  1.59746890863949E-1_wp, &
      & -2.70765876809499E-1_wp, -3.27435355522312E-1_wp, -4.70046325670683E-2_wp, &
      &  1.10838969762146E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_kernel_numt(error, mol, kernel_enum%still, keps, qat)

end subroutine test_kernel_third_still


!> Test P16 kernel gradient against numerical derivative
subroutine test_kernel_gradient_p16(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: keps = 0.5_wp
   real(wp), parameter :: qat(*) = [&
      & -2.11018727757438E-1_wp, -6.04389222813257E-2_wp, -1.90601159250311E-1_wp, &
      &  1.49237694872530E-1_wp,  1.35835820853652E-1_wp,  1.27732431639016E-1_wp, &
      &  1.78559147201780E-1_wp,  1.42324484825195E-1_wp,  1.92106458233743E-1_wp, &
      &  1.45841758574287E-1_wp,  1.56456166394024E-1_wp,  1.59746890863949E-1_wp, &
      & -2.70765876809499E-1_wp, -3.27435355522312E-1_wp, -4.70046325670683E-2_wp, &
      &  1.10838969762146E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_kernel_numg(error, mol, kernel_enum%p16, keps, qat)

end subroutine test_kernel_gradient_p16


!> Test P16 kernel Hessian against numerical derivative
subroutine test_kernel_hessian_p16(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: keps = 0.5_wp
   real(wp), parameter :: qat(*) = [&
      & -2.11018727757438E-1_wp, -6.04389222813257E-2_wp, -1.90601159250311E-1_wp, &
      &  1.49237694872530E-1_wp,  1.35835820853652E-1_wp,  1.27732431639016E-1_wp, &
      &  1.78559147201780E-1_wp,  1.42324484825195E-1_wp,  1.92106458233743E-1_wp, &
      &  1.45841758574287E-1_wp,  1.56456166394024E-1_wp,  1.59746890863949E-1_wp, &
      & -2.70765876809499E-1_wp, -3.27435355522312E-1_wp, -4.70046325670683E-2_wp, &
      &  1.10838969762146E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_kernel_numh(error, mol, kernel_enum%p16, keps, qat)

end subroutine test_kernel_hessian_p16


!> Test P16 kernel Hessian against numerical derivative
subroutine test_kernel_third_p16(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: keps = 0.5_wp
   real(wp), parameter :: qat(*) = [&
      & -2.11018727757438E-1_wp, -6.04389222813257E-2_wp, -1.90601159250311E-1_wp, &
      &  1.49237694872530E-1_wp,  1.35835820853652E-1_wp,  1.27732431639016E-1_wp, &
      &  1.78559147201780E-1_wp,  1.42324484825195E-1_wp,  1.92106458233743E-1_wp, &
      &  1.45841758574287E-1_wp,  1.56456166394024E-1_wp,  1.59746890863949E-1_wp, &
      & -2.70765876809499E-1_wp, -3.27435355522312E-1_wp, -4.70046325670683E-2_wp, &
      &  1.10838969762146E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_kernel_numt(error, mol, kernel_enum%p16, keps, qat)

end subroutine test_kernel_third_p16

!> Test kernel gradient against numerical one.
subroutine test_kernel_numg(error, mol, kernel_id, keps, qat)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type), intent(inout) :: mol
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   real(wp), intent(in) :: qat(:)

   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel
   real(wp), allocatable :: rvdw(:), rad(:), draddr(:, :, :)
   real(wp), allocatable :: amat_r(:, :), amat_l(:, :)
   real(wp), allocatable :: numg_kernel(:, :, :, :)
   real(wp), allocatable :: ana_ij(:, :)
   real(wp), parameter :: step = 1.0e-6_wp
   integer :: iat, jat, ic, jc
   real(wp) :: maxdiff

   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   kernel = new_kernel(kernel_id, keps)

   allocate(rad(mol%nat), draddr(3, mol%nat, mol%nat))
   allocate(amat_r(mol%nat, mol%nat), amat_l(mol%nat, mol%nat))
   allocate(numg_kernel(3, mol%nat, mol%nat, mol%nat))
   allocate(ana_ij(3, mol%nat))

   ! --- Numerical derivative: numg_kernel(alpha, k, m, n) = dK_mn / dr_k,alpha
   numg_kernel(:, :, :, :) = 0.0_wp

   do iat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call gbobc%get_rad(mol, rad)
         amat_r(:, :) = 0.0_wp
         call kernel%add_kernel_mat(mol%nat, mol%xyz, rad, amat_r)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call gbobc%get_rad(mol, rad)
         amat_l(:, :) = 0.0_wp
         call kernel%add_kernel_mat(mol%nat, mol%xyz, rad, amat_l)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step

         do jat = 1, mol%nat
            do jc = 1, mol%nat
               numg_kernel(ic, iat, jat, jc) = 0.5_wp * (amat_r(jat, jc) - amat_l(jat, jc)) / step
            end do
         end do
      end do
   end do

   ! --- Analytical: get brad + brdr once at reference geometry
   call gbobc%get_rad(mol, rad, draddr)

   ! Compare: for each (m,n), compute ana_ij(:,k)=dK_mn/dr_k and compare to numg_kernel(:,:,m,n)
   do jat = 1, mol%nat
      do jc = 1, mol%nat
         call compute_kernel_dkdr_ij(kernel_id, keps, mol%nat, mol%xyz, rad, draddr, jat, jc, ana_ij)

         maxdiff = maxval(abs(ana_ij(:, :) - numg_kernel(:, :, jat, jc)))
         if (maxdiff > thr2) then
            call test_failed(error, "Kernel gradient does not match finite difference solution")
            print '(a,2i6, a, es20.13)', "Mismatch at (m,n)=(", jat, jc, "), max|diff|=", maxdiff
            print '(a)', "Analytical dK_mn/dr_k (3,nat):"
            print '(3es20.13)', ana_ij
            print '(a)', "Numerical dK_mn/dr_k (3,nat):"
            print '(3es20.13)', numg_kernel(:, :, jat, jc)
            print '(a)', "Difference (ana - num):"
            print '(3es20.13)', ana_ij - numg_kernel(:, :, jat, jc)
            return
         end if
      end do
   end do
end subroutine test_kernel_numg

!> Test kernel 2nd derivative against numerical finite difference of the 
!> analytical 1st derivative.
subroutine test_kernel_numh(error, mol, kernel_id, keps, qat)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type), intent(inout) :: mol
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   real(wp), intent(in) :: qat(:)

   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel

   real(wp), allocatable :: rvdw(:), rad(:)
   real(wp), allocatable :: draddr(:, :, :)
   real(wp), allocatable :: draddr2(:, :, :, :, :)

   real(wp), allocatable :: dkdr_p_ij(:, :)     ! (3,nat)
   real(wp), allocatable :: dkdr_m_ij(:, :)     ! (3,nat)
   real(wp), allocatable :: num2(:, :, :, :)    ! (3,nat,3,nat) for one (m,n)
   real(wp), allocatable :: ana2(:, :, :, :)    ! (3,nat,3,nat) for one (m,n)

   real(wp), parameter :: step = 1.0e-6_wp
   integer :: nat, m, n, l, beta, k, alpha
   real(wp) :: diff, maxdiff
   integer :: imax_alpha, imax_k, imax_beta, imax_l, imax_m, imax_n

   nat = mol%nat

   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   kernel = new_kernel(kernel_id, keps)

   allocate(rad(nat))
   allocate(draddr(3, nat, nat))
   allocate(draddr2(3, nat, 3, nat, nat))

   allocate(dkdr_p_ij(3, nat), dkdr_m_ij(3, nat))
   allocate(num2(3, nat, 3, nat))
   allocate(ana2(3, nat, 3, nat))

   ! Analytical Born radii derivatives at reference geometry (need brdr2!)
   call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2)

   maxdiff = 0.0_wp
   imax_alpha=1; imax_k=1; imax_beta=1; imax_l=1; imax_m=1; imax_n=1

   do m = 1, nat
      do n = 1, nat

         ! ---- Analytical slice for this (m,n)
         call compute_kernel_d2kdr2_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, m, n, ana2)

         ! ---- Numerical slice for this (m,n) by FD of dkdr_ij
         num2(:, :, :, :) = 0.0_wp

         do l = 1, nat
            do beta = 1, 3

               ! +step
               mol%xyz(beta, l) = mol%xyz(beta, l) + step
               call gbobc%get_rad(mol, rad, draddr)
               call compute_kernel_dkdr_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, m, n, dkdr_p_ij)

               ! -step
               mol%xyz(beta, l) = mol%xyz(beta, l) - 2.0_wp*step
               call gbobc%get_rad(mol, rad, draddr)
               call compute_kernel_dkdr_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, m, n, dkdr_m_ij)

               ! restore
               mol%xyz(beta, l) = mol%xyz(beta, l) + step

               do k = 1, nat
                  do alpha = 1, 3
                     num2(alpha, k, beta, l) = 0.5_wp * (dkdr_p_ij(alpha, k) - dkdr_m_ij(alpha, k)) / step
                  end do
               end do

            end do
         end do

         ! ---- Compare this slice
         do k = 1, nat
            do alpha = 1, 3
               do l = 1, nat
                  do beta = 1, 3
                     diff = abs(ana2(alpha, k, beta, l) - num2(alpha, k, beta, l))
                     if (diff > maxdiff) then
                        maxdiff = diff
                        imax_alpha = alpha
                        imax_k     = k
                        imax_beta  = beta
                        imax_l     = l
                        imax_m     = m
                        imax_n     = n
                     end if
                  end do
               end do
            end do
         end do

      end do
   end do

   if (maxdiff > thr2) then
      call test_failed(error, "Kernel second derivative does not match finite difference solution")
      print '(a,es20.13)', "Max |d2K/dr2| difference: ", maxdiff
      print '(a,6(i0,1x))', "At indices (alpha,k,beta,l,m,n) = ", imax_alpha, imax_k, imax_beta, imax_l, imax_m, imax_n
      ! Recompute the offending slice and report the single entry
      call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2)
      call compute_kernel_d2kdr2_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, imax_m, imax_n, ana2)

      mol%xyz(imax_beta, imax_l) = mol%xyz(imax_beta, imax_l) + step
      call gbobc%get_rad(mol, rad, draddr)
      call compute_kernel_dkdr_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, imax_m, imax_n, dkdr_p_ij)

      mol%xyz(imax_beta, imax_l) = mol%xyz(imax_beta, imax_l) - 2.0_wp*step
      call gbobc%get_rad(mol, rad, draddr)
      call compute_kernel_dkdr_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, imax_m, imax_n, dkdr_m_ij)

      mol%xyz(imax_beta, imax_l) = mol%xyz(imax_beta, imax_l) + step

      num2(imax_alpha, imax_k, imax_beta, imax_l) = 0.5_wp * (dkdr_p_ij(imax_alpha, imax_k) - dkdr_m_ij(imax_alpha, imax_k)) / step

      print '(a,es20.13)', "Analytical value: ", ana2(imax_alpha, imax_k, imax_beta, imax_l)
      print '(a,es20.13)', "Numerical  value: ", num2(imax_alpha, imax_k, imax_beta, imax_l)
      print '(a,es20.13)', "Difference       : ", ana2(imax_alpha, imax_k, imax_beta, imax_l) - num2(imax_alpha, imax_k, imax_beta, imax_l)
   end if

end subroutine test_kernel_numh



!> Test kernel 3rd derivative against numerical finite difference of the 
!> analytical 2nd derivative.
subroutine test_kernel_numt(error, mol, kernel_id, keps, qat)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type), intent(inout) :: mol
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   real(wp), intent(in) :: qat(:)

   type(born_integrator) :: gbobc
   class(kernel_type), allocatable :: kernel

   real(wp), allocatable :: rvdw(:), rad(:)
   real(wp), allocatable :: draddr(:, :, :)
   real(wp), allocatable :: draddr2(:, :, :, :, :)
   real(wp), allocatable :: draddr3(:, :, :, :, :, :, :)

   real(wp), allocatable :: d2p(:, :, :, :)      ! (3,nat,3,nat) for one (i,j)
   real(wp), allocatable :: d2m(:, :, :, :)      ! (3,nat,3,nat) for one (i,j)
   real(wp), allocatable :: num3(:, :, :, :, :, :) ! (3,nat,3,nat,3,nat) for one (i,j)
   real(wp), allocatable :: ana3(:, :, :, :, :, :) ! (3,nat,3,nat,3,nat) for one (i,j)

   real(wp), parameter :: step = 1.0e-6_wp
   integer :: nat, i, j, k, l, m, alpha, beta, gamma
   real(wp) :: diff, maxdiff
   integer :: ia, ik, ib, il, ig, im, ii, ij  ! index record for reporting

   nat = mol%nat
   rvdw = get_vdw_rad_d3(mol%num)
   call new_born_integrator(gbobc, mol, rvdw)
   kernel = new_kernel(kernel_id, keps)

   allocate(rad(nat))
   allocate(draddr(3,nat,nat))
   allocate(draddr2(3,nat,3,nat,nat))
   allocate(draddr3(3,nat,3,nat,3,nat,nat))

   allocate(d2p(3,nat,3,nat), d2m(3,nat,3,nat))
   allocate(num3(3,nat,3,nat,3,nat))
   allocate(ana3(3,nat,3,nat,3,nat))

   ! reference geometry: need up to brdr3 for analytical
   call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2, dradd3r=draddr3)

   maxdiff = 0.0_wp
   ia=1;ik=1;ib=1;il=1;ig=1;im=1;ii=1;ij=1

   do i = 1, nat
      do j = 1, nat

         ! ---- analytical 3rd-derivative slab for this (i,j)
         call compute_kernel_d3kdr3_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, draddr3, i, j, ana3)

         ! ---- numerical slab via FD of the ij Hessian
         num3(:, :, :, :, :, :) = 0.0_wp

         do m = 1, nat
            do gamma = 1, 3

               mol%xyz(gamma, m) = mol%xyz(gamma, m) + step
               call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2)
               call compute_kernel_d2kdr2_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, i, j, d2p)

               mol%xyz(gamma, m) = mol%xyz(gamma, m) - 2.0_wp*step
               call gbobc%get_rad(mol, rad, draddr, dradd2r=draddr2)
               call compute_kernel_d2kdr2_ij(kernel_id, keps, nat, mol%xyz, rad, draddr, draddr2, i, j, d2m)

               mol%xyz(gamma, m) = mol%xyz(gamma, m) + step

               do k = 1, nat
                  do alpha = 1, 3
                     do l = 1, nat
                        do beta = 1, 3
                           num3(alpha,k,beta,l,gamma,m) = 0.5_wp * (d2p(alpha,k,beta,l) - d2m(alpha,k,beta,l)) / step
                        end do
                     end do
                  end do
               end do

            end do
         end do

         ! ---- compare this slab; track worst entry globally
         do k = 1, nat
            do alpha = 1, 3
               do l = 1, nat
                  do beta = 1, 3
                     do m = 1, nat
                        do gamma = 1, 3
                           diff = abs(ana3(alpha,k,beta,l,gamma,m) - num3(alpha,k,beta,l,gamma,m))
                           if (diff > maxdiff) then
                              maxdiff = diff
                              ia=alpha; ik=k; ib=beta; il=l; ig=gamma; im=m; ii=i; ij=j
                           end if
                        end do
                     end do
                  end do
               end do
            end do
         end do

      end do
   end do

   if (maxdiff > thr2) then
      call test_failed(error, "Kernel third derivative does not match finite difference solution")
      print '(a,es20.13)', "Max |d3K/dr3| difference: ", maxdiff
      print '(a,8(i0,1x))', "At indices (a,k,b,l,g,m,i,j) = ", ia,ik,ib,il,ig,im,ii,ij
   end if

end subroutine test_kernel_numt





end module test_solvation_kernel
