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

!> @file tblite/solvation/kernel/still.f90
!> Provides the Still kernel implementation for GBSA solvation.

module tblite_solvation_kernel_still
   use mctc_env, only: wp
   use tblite_blas, only: gemv
   use tblite_solvation_kernel_type, only: kernel_type
   implicit none
   private

   public :: still_kernel

   type, extends(kernel_type) :: still_kernel
   contains
      procedure :: kernel_K => still_K
      procedure :: kernel_dKdr => still_dKdr
      procedure :: kernel_dKdborn => still_dKdborn
      procedure :: kernel_d_dKdr_dborn => still_d_dKdr_dborn
      procedure :: kernel_d2Kdr2 => still_d2Kdr2
      procedure :: kernel_d_d2Kdr2_dborn => still_d_d2Kdr2_dborn
      procedure :: kernel_d3Kdr3 => still_d3Kdr3
      procedure :: kernel_d_d3Kdr3_dborn => still_d_d3Kdr3_dborn
      procedure :: kernel_d4Kdr4 => still_d4Kdr4
      procedure :: kernel_d_d4Kdr4_dborn => still_d_d4Kdr4_dborn
      procedure :: kernel_d5Kdr5 => still_d5Kdr5
   end type still_kernel

contains

   pure subroutine still_K(self, nat, xyz, brad, Amat)
      !> Instance of Still kernel
      class(still_kernel), intent(in) :: self
      !> Number of atoms
      integer, intent(in) :: nat
      !> Cartesian coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Born radii
      real(wp), intent(in) :: brad(:)
      !> Charge-charge interaction matrix
      real(wp), intent(inout) :: Amat(:, :)

      integer  :: i, j
      real(wp), parameter :: a4 = 0.25_wp
      real(wp) :: aa, vec(3), r1, r2, bp
      real(wp) :: dd, expd, fgb2, dfgb

      Amat = 0.0_wp

      do i = 1, nat
         do j = 1, i-1
            vec(:) = xyz(:, i)-xyz(:, j)
            r1 = norm2(vec)
            r2 = r1*r1

            aa = brad(i)*brad(j)
            dd = a4*r2/aa
            expd = exp(-dd)
            fgb2 = r2+aa*expd
            dfgb = 1.0_wp/sqrt(fgb2)

            Amat(i, j) = self%keps*dfgb+Amat(i, j)
            Amat(j, i) = self%keps*dfgb+Amat(j, i)
         end do

         bp = 1.0_wp/brad(i)
         Amat(i, i) = Amat(i, i)+self%keps*bp
      end do
   end subroutine still_K

 
   subroutine still_dKdr(self, rA, rB, bornA, bornB, d1)
   !! d1(i) = ∂K/∂R_A,i where K = 1/sqrt(r^2 + a*exp(-r^2/(4a))), a=bornA*bornB
      !> Instance of Still kernel
      class(still_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in) :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in) :: rB(3)
      !> Born radius of atom A
      real(wp), intent(in) :: bornA
      !> Born radius of atom B
      real(wp), intent(in) :: bornB
      !> First derivative of Still kernel
      real(wp), intent(out) :: d1(3)

      real(wp) :: rvec(3), s, a, e, u
      real(wp) :: invu, invsqrtu, um3
      real(wp) :: u1, u2, u3, u4
      real(wp) :: Acoef, Bcoef, Ccoef, Dcoef

      call still_scalar_coeffs(rA, rB, bornA, bornB, rvec, s, a, e, u, invu, invsqrtu, um3, &
                               u1, u2, u3, u4, Acoef, Bcoef, Ccoef, Dcoef)

      ! ∂K/∂r_i = 2 r_i * dK/ds
      d1 = self%keps*2.0_wp*rvec*Acoef
   end subroutine still_dKdr

    
   pure subroutine still_dKdborn(self, rA, rB, bornA, bornB, dK_bA, dK_bB)
      !! dK_bA = ∂K/∂bornA, dK_bB = ∂K/∂bornB (coordinates held fixed)
      !! K = keps / sqrt( s + a*exp(-s/(4a)) ),  a = bornA*bornB, s = |rA-rB|^2
      class(still_kernel), intent(in) :: self
      real(wp), intent(in)  :: rA(3), rB(3)
      real(wp), intent(in)  :: bornA, bornB
      real(wp), intent(out) :: dK_bA, dK_bB

      real(wp) :: rvec(3), s, a, e, u
      real(wp) :: invu, invsqrtu, um3
      real(wp) :: du_da, dK_da

      rvec = rA - rB
      s = dot_product(rvec, rvec)

      a = bornA*bornB
      if (a <= 0.0_wp) then
         dK_bA = 0.0_wp
         dK_bB = 0.0_wp
         return
      end if

      ! Handle s==0 robustly (coincident coordinates):
      ! u = a, K = keps / sqrt(a), dK/da = -0.5 * keps * a^(-3/2)
      if (s == 0.0_wp) then
         dK_da = self%keps * (-0.5_wp) / (a*sqrt(a))
         dK_bA = bornB * dK_da
         dK_bB = bornA * dK_da
         return
      end if

      e = exp(-s/(4.0_wp*a))
      u = s + a*e

      invu     = 1.0_wp/u
      invsqrtu = 1.0_wp/sqrt(u)
      um3      = invu*invsqrtu   ! u^(-3/2)

      ! du/da = e*(1 + s/(4a))
      du_da = e * (1.0_wp + s/(4.0_wp*a))

      ! dK/da = keps * (-1/2) * u^(-3/2) * du/da
      dK_da = self%keps * (-0.5_wp) * um3 * du_da

      dK_bA = bornB * dK_da
      dK_bB = bornA * dK_da
   end subroutine still_dKdborn

   subroutine still_d_dKdr_dborn(self, rA, rB, bornA, bornB, d1_bA, d1_bB)
      !> Instance of Still kernel
      class(still_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in) :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in) :: rB(3)
      !> Born radius of atom A
      real(wp), intent(in) :: bornA
      !> Born radius of atom B
      real(wp), intent(in) :: bornB
      !> Derivative wrt Born radius of atom A
      real(wp), intent(out) :: d1_bA(3)
      !> Derivative wrt Born radius of atom B
      real(wp), intent(out) :: d1_bB(3)
      

      real(wp) :: rvec(3), s, a, e, u
      real(wp) :: invu, invsqrtu, um3
      real(wp) :: u1, u2, u3, u4
      real(wp) :: Acoef, Bcoef, Ccoef, Dcoef
      real(wp) :: dA_da, dB_da, dC_da, dD_da
      real(wp) :: dA_bA, dA_bB

      call still_scalar_coeffs(rA, rB, bornA, bornB, rvec, s, a, e, u, invu, invsqrtu, um3, &
                               u1, u2, u3, u4, Acoef, Bcoef, Ccoef, Dcoef)

      if (a <= 0.0_wp .or. s == 0.0_wp) then
         d1_bA = 0.0_wp
         d1_bB = 0.0_wp
         return
      end if

      call still_scalar_coeffs_da(s, a, e, u, u1, u2, u3, u4, dA_da, dB_da, dC_da, dD_da)

      dA_bA = bornB*dA_da
      dA_bB = bornA*dA_da

      d1_bA = self%keps*2.0_wp*rvec*dA_bA
      d1_bB = self%keps*2.0_wp*rvec*dA_bB
   end subroutine still_d_dKdr_dborn

   subroutine still_d2Kdr2(self, rA, rB, bornA, bornB, d2)
   !! d2(i,j) = ∂²K/∂R_A,i ∂R_A,j
      !> Instance of Still kernel
      class(still_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in) :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in) :: rB(3)
      !> Born radius of atom A
      real(wp), intent(in) :: bornA
      !> Born radius of atom B
      real(wp), intent(in) :: bornB
      !> Second derivative of Still kernel
      real(wp), intent(out) :: d2(3, 3)

      real(wp) :: rvec(3), s, a, e, u
      real(wp) :: invu, invsqrtu, um3
      real(wp) :: u1, u2, u3, u4
      real(wp) :: Acoef, Bcoef, Ccoef, Dcoef
      integer :: i, j

      call still_scalar_coeffs(rA, rB, bornA, bornB, rvec, s, a, e, u, invu, invsqrtu, um3, &
                               u1, u2, u3, u4, Acoef, Bcoef, Ccoef, Dcoef)

      ! For radial K(s): ∂i∂j K = 2 δij K'(s) + 4 r_i r_j K''(s)
      do i = 1, 3
         do j = 1, 3
            d2(i, j) = self%keps*4.0_wp*rvec(i)*rvec(j)*Bcoef
            if (i == j) d2(i, j) = d2(i, j)+self%keps*2.0_wp*Acoef
         end do
      end do
   end subroutine still_d2Kdr2

   subroutine still_d_d2Kdr2_dborn(self, rA, rB, bornA, bornB, d2_bA, d2_bB)
      !> Instance of Still kernel
      class(still_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in) :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in) :: rB(3)
      !> Born radius of atom A
      real(wp), intent(in) :: bornA
      !> Born radius of atom B
      real(wp), intent(in) :: bornB
      !> Derivative wrt Born radius of atom A
      real(wp), intent(out) :: d2_bA(3, 3)
      !> Derivative wrt Born radius of atom B
      real(wp), intent(out) :: d2_bB(3, 3)

      real(wp) :: rvec(3), s, a, e, u
      real(wp) :: invu, invsqrtu, um3
      real(wp) :: u1, u2, u3, u4
      real(wp) :: Acoef, Bcoef, Ccoef, Dcoef
      real(wp) :: dA_da, dB_da, dC_da, dD_da
      real(wp) :: dA_bA, dB_bA, dA_bB, dB_bB
      integer :: i, j

      call still_scalar_coeffs(rA, rB, bornA, bornB, rvec, s, a, e, u, invu, invsqrtu, um3, &
                               u1, u2, u3, u4, Acoef, Bcoef, Ccoef, Dcoef)

      if (a <= 0.0_wp .or. s == 0.0_wp) then
         d2_bA = 0.0_wp
         d2_bB = 0.0_wp
         return
      end if

      call still_scalar_coeffs_da(s, a, e, u, u1, u2, u3, u4, dA_da, dB_da, dC_da, dD_da)

      dA_bA = bornB*dA_da
      dB_bA = bornB*dB_da
      dA_bB = bornA*dA_da
      dB_bB = bornA*dB_da

      do i = 1, 3
         do j = 1, 3
            d2_bA(i, j) = self%keps*4.0_wp*rvec(i)*rvec(j)*dB_bA
            d2_bB(i, j) = self%keps*4.0_wp*rvec(i)*rvec(j)*dB_bB
            if (i == j) then
               d2_bA(i, j) = d2_bA(i, j)+self%keps*2.0_wp*dA_bA
               d2_bB(i, j) = d2_bB(i, j)+self%keps*2.0_wp*dA_bB
            end if
         end do
      end do
   end subroutine still_d_d2Kdr2_dborn

   subroutine still_d3Kdr3(self, rA, rB, bornA, bornB, d3)
    !! d3(i,j,k) = ∂³K/∂R_A,i ∂R_A,j ∂R_A,k
      !> Instance of Still kernel
      class(still_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in)  :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in)  :: rB(3)
      !> Born radius of atom A
      real(wp), intent(in)  :: bornA
      !> Born radius of atom B
      real(wp), intent(in)  :: bornB
      !> Third derivative of Still kernel
      real(wp), intent(out) :: d3(3, 3, 3)

      real(wp) :: rvec(3), s, a, e, u
      real(wp) :: invu, invsqrtu, um3
      real(wp) :: u1, u2, u3, u4
      real(wp) :: Acoef, Bcoef, Ccoef, Dcoef
      integer :: i, j, k

    call still_scalar_coeffs(rA, rB, bornA, bornB, rvec, s, a, e, u, invu, invsqrtu, um3, u1, u2, u3, u4, Acoef, Bcoef, Ccoef, Dcoef)

      ! ∂i∂j∂k K =
      !   4(δij r_k + δik r_j + δjk r_i) K''(s) + 8 r_i r_j r_k K'''(s)
      do i = 1, 3
         do j = 1, 3
            do k = 1, 3
               d3(i, j, k) = self%keps*8.0_wp*rvec(i)*rvec(j)*rvec(k)*Ccoef
               if (i == j) d3(i, j, k) = d3(i, j, k)+self%keps*4.0_wp*rvec(k)*Bcoef
               if (i == k) d3(i, j, k) = d3(i, j, k)+self%keps*4.0_wp*rvec(j)*Bcoef
               if (j == k) d3(i, j, k) = d3(i, j, k)+self%keps*4.0_wp*rvec(i)*Bcoef
            end do
         end do
      end do
   end subroutine still_d3Kdr3

   subroutine still_d_d3Kdr3_dborn(self, rA, rB, bornA, bornB, d3_bA, d3_bB)
      !> Instance of Still kernel
      class(still_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in)  :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in)  :: rB(3)
      !> Born radius of atom A
      real(wp), intent(in)  :: bornA
      !> Born radius of atom B
      real(wp), intent(in)  :: bornB
      !> Derivative wrt Born radius of atom A
      real(wp), intent(out) :: d3_bA(3, 3, 3)
      !> Derivative wrt Born radius of atom B
      real(wp), intent(out) :: d3_bB(3, 3, 3)

      real(wp) :: rvec(3), s, a, e, u
      real(wp) :: invu, invsqrtu, um3
      real(wp) :: u1, u2, u3, u4
      real(wp) :: Acoef, Bcoef, Ccoef, Dcoef
      real(wp) :: dA_da, dB_da, dC_da, dD_da
      real(wp) :: dB_bA, dC_bA, dB_bB, dC_bB
      integer :: i, j, k

      call still_scalar_coeffs(rA, rB, bornA, bornB, rvec, s, a, e, u, invu, invsqrtu, um3, &
                               u1, u2, u3, u4, Acoef, Bcoef, Ccoef, Dcoef)

      if (a <= 0.0_wp .or. s == 0.0_wp) then
         d3_bA = 0.0_wp; d3_bB = 0.0_wp
         return
      end if

      call still_scalar_coeffs_da(s, a, e, u, u1, u2, u3, u4, dA_da, dB_da, dC_da, dD_da)

      dB_bA = bornB*dB_da
      dC_bA = bornB*dC_da
      dB_bB = bornA*dB_da
      dC_bB = bornA*dC_da

      do i = 1, 3
         do j = 1, 3
            do k = 1, 3
               d3_bA(i, j, k) = self%keps*8.0_wp*rvec(i)*rvec(j)*rvec(k)*dC_bA
               d3_bB(i, j, k) = self%keps*8.0_wp*rvec(i)*rvec(j)*rvec(k)*dC_bB
               if (i == j) then
                  d3_bA(i, j, k) = d3_bA(i, j, k)+self%keps*4.0_wp*rvec(k)*dB_bA
                  d3_bB(i, j, k) = d3_bB(i, j, k)+self%keps*4.0_wp*rvec(k)*dB_bB
               end if
               if (i == k) then
                  d3_bA(i, j, k) = d3_bA(i, j, k)+self%keps*4.0_wp*rvec(j)*dB_bA
                  d3_bB(i, j, k) = d3_bB(i, j, k)+self%keps*4.0_wp*rvec(j)*dB_bB
               end if
               if (j == k) then
                  d3_bA(i, j, k) = d3_bA(i, j, k)+self%keps*4.0_wp*rvec(i)*dB_bA
                  d3_bB(i, j, k) = d3_bB(i, j, k)+self%keps*4.0_wp*rvec(i)*dB_bB
               end if
            end do
         end do
      end do
   end subroutine still_d_d3Kdr3_dborn

   subroutine still_d4Kdr4(self, rA, rB, bornA, bornB, d4)
    !! d4(i,j,k,l) = ∂⁴K/∂R_A,i ∂R_A,j ∂R_A,k ∂R_A,l
      !> Instance of Still kernel
      class(still_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in)  :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in)  :: rB(3)
      !> Born radius of atom A
      real(wp), intent(in)  :: bornA
      !> Born radius of atom B
      real(wp), intent(in)  :: bornB
      !> Fourth derivative of Still kernel
      real(wp), intent(out) :: d4(3, 3, 3, 3)

      real(wp) :: rvec(3), s, a, e, u
      real(wp) :: invu, invsqrtu, um3
      real(wp) :: u1, u2, u3, u4
      real(wp) :: Acoef, Bcoef, Ccoef, Dcoef
      integer :: i, j, k, l

    call still_scalar_coeffs(rA, rB, bornA, bornB, rvec, s, a, e, u, invu, invsqrtu, um3, u1, u2, u3, u4, Acoef, Bcoef, Ccoef, Dcoef)

      ! ∂i∂j∂k∂l K =
      !   4(δijδkl + δikδjl + δilδjk) K''(s)
      ! + 8(δij r_k r_l + δik r_j r_l + δil r_j r_k + δjk r_i r_l + δjl r_i r_k + δkl r_i r_j) K'''(s)
      ! + 16 r_i r_j r_k r_l K''''(s)
      do i = 1, 3
         do j = 1, 3
            do k = 1, 3
               do l = 1, 3
                  d4(i, j, k, l) = self%keps*16.0_wp*rvec(i)*rvec(j)*rvec(k)*rvec(l)*Dcoef

                  ! B terms
                  if (i == j .and. k == l) d4(i, j, k, l) = d4(i, j, k, l)+self%keps*4.0_wp*Bcoef
                  if (i == k .and. j == l) d4(i, j, k, l) = d4(i, j, k, l)+self%keps*4.0_wp*Bcoef
                  if (i == l .and. j == k) d4(i, j, k, l) = d4(i, j, k, l)+self%keps*4.0_wp*Bcoef

                  ! C terms (6 permutations)
                  if (i == j) d4(i, j, k, l) = d4(i, j, k, l)+self%keps*8.0_wp*rvec(k)*rvec(l)*Ccoef
                  if (i == k) d4(i, j, k, l) = d4(i, j, k, l)+self%keps*8.0_wp*rvec(j)*rvec(l)*Ccoef
                  if (i == l) d4(i, j, k, l) = d4(i, j, k, l)+self%keps*8.0_wp*rvec(j)*rvec(k)*Ccoef
                  if (j == k) d4(i, j, k, l) = d4(i, j, k, l)+self%keps*8.0_wp*rvec(i)*rvec(l)*Ccoef
                  if (j == l) d4(i, j, k, l) = d4(i, j, k, l)+self%keps*8.0_wp*rvec(i)*rvec(k)*Ccoef
                  if (k == l) d4(i, j, k, l) = d4(i, j, k, l)+self%keps*8.0_wp*rvec(i)*rvec(j)*Ccoef
               end do
            end do
         end do
      end do
   end subroutine still_d4Kdr4

   subroutine still_d_d4Kdr4_dborn(self, rA, rB, bornA, bornB, d4_bA, d4_bB)
      !> Instance of Still kernel
      class(still_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in)  :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in)  :: rB(3)
      !> Born radius of atom A
      real(wp), intent(in)  :: bornA
      !> Born radius of atom B
      real(wp), intent(in)  :: bornB
      !> Derivative wrt Born radius of atom A
      real(wp), intent(out) :: d4_bA(3, 3, 3, 3)
      !> Derivative wrt Born radius of atom B
      real(wp), intent(out) :: d4_bB(3, 3, 3, 3)

      real(wp) :: rvec(3), s, a, e, u
      real(wp) :: invu, invsqrtu, um3
      real(wp) :: u1, u2, u3, u4
      real(wp) :: Acoef, Bcoef, Ccoef, Dcoef
      real(wp) :: dA_da, dB_da, dC_da, dD_da
      real(wp) :: dB_bA, dC_bA, dD_bA, dB_bB, dC_bB, dD_bB
      integer :: i, j, k, l

      call still_scalar_coeffs(rA, rB, bornA, bornB, rvec, s, a, e, u, invu, invsqrtu, um3, &
                               u1, u2, u3, u4, Acoef, Bcoef, Ccoef, Dcoef)

      if (a <= 0.0_wp .or. s == 0.0_wp) then
         d4_bA = 0.0_wp; d4_bB = 0.0_wp
         return
      end if

      call still_scalar_coeffs_da(s, a, e, u, u1, u2, u3, u4, dA_da, dB_da, dC_da, dD_da)

      dB_bA = bornB*dB_da
      dC_bA = bornB*dC_da
      dD_bA = bornB*dD_da

      dB_bB = bornA*dB_da
      dC_bB = bornA*dC_da
      dD_bB = bornA*dD_da

      do i = 1, 3
         do j = 1, 3
            do k = 1, 3
               do l = 1, 3
                  d4_bA(i, j, k, l) = self%keps*16.0_wp*rvec(i)*rvec(j)*rvec(k)*rvec(l)*dD_bA
                  d4_bB(i, j, k, l) = self%keps*16.0_wp*rvec(i)*rvec(j)*rvec(k)*rvec(l)*dD_bB

                  ! B terms
                  if (i == j .and. k == l) then
                     d4_bA(i, j, k, l) = d4_bA(i, j, k, l)+self%keps*4.0_wp*dB_bA
                     d4_bB(i, j, k, l) = d4_bB(i, j, k, l)+self%keps*4.0_wp*dB_bB
                  end if
                  if (i == k .and. j == l) then
                     d4_bA(i, j, k, l) = d4_bA(i, j, k, l)+self%keps*4.0_wp*dB_bA
                     d4_bB(i, j, k, l) = d4_bB(i, j, k, l)+self%keps*4.0_wp*dB_bB
                  end if
                  if (i == l .and. j == k) then
                     d4_bA(i, j, k, l) = d4_bA(i, j, k, l)+self%keps*4.0_wp*dB_bA
                     d4_bB(i, j, k, l) = d4_bB(i, j, k, l)+self%keps*4.0_wp*dB_bB
                  end if

                  ! C terms (6 permutations)
                  if (i == j) then
                     d4_bA(i, j, k, l) = d4_bA(i, j, k, l)+self%keps*8.0_wp*rvec(k)*rvec(l)*dC_bA
                     d4_bB(i, j, k, l) = d4_bB(i, j, k, l)+self%keps*8.0_wp*rvec(k)*rvec(l)*dC_bB
                  end if
                  if (i == k) then
                     d4_bA(i, j, k, l) = d4_bA(i, j, k, l)+self%keps*8.0_wp*rvec(j)*rvec(l)*dC_bA
                     d4_bB(i, j, k, l) = d4_bB(i, j, k, l)+self%keps*8.0_wp*rvec(j)*rvec(l)*dC_bB
                  end if
                  if (i == l) then
                     d4_bA(i, j, k, l) = d4_bA(i, j, k, l)+self%keps*8.0_wp*rvec(j)*rvec(k)*dC_bA
                     d4_bB(i, j, k, l) = d4_bB(i, j, k, l)+self%keps*8.0_wp*rvec(j)*rvec(k)*dC_bB
                  end if
                  if (j == k) then
                     d4_bA(i, j, k, l) = d4_bA(i, j, k, l)+self%keps*8.0_wp*rvec(i)*rvec(l)*dC_bA
                     d4_bB(i, j, k, l) = d4_bB(i, j, k, l)+self%keps*8.0_wp*rvec(i)*rvec(l)*dC_bB
                  end if
                  if (j == l) then
                     d4_bA(i, j, k, l) = d4_bA(i, j, k, l)+self%keps*8.0_wp*rvec(i)*rvec(k)*dC_bA
                     d4_bB(i, j, k, l) = d4_bB(i, j, k, l)+self%keps*8.0_wp*rvec(i)*rvec(k)*dC_bB
                  end if
                  if (k == l) then
                     d4_bA(i, j, k, l) = d4_bA(i, j, k, l)+self%keps*8.0_wp*rvec(i)*rvec(j)*dC_bA
                     d4_bB(i, j, k, l) = d4_bB(i, j, k, l)+self%keps*8.0_wp*rvec(i)*rvec(j)*dC_bB
                  end if
               end do
            end do
         end do
      end do
   end subroutine still_d_d4Kdr4_dborn

   subroutine still_d5Kdr5(self, rA, rB, bornA, bornB, d5)
    !! d5(i,j,k,l,m) = ∂⁵K/∂R_A,i ∂R_A,j ∂R_A,k ∂R_A,l ∂R_A,m
    !!
    !! For radial K(s), s = r·r:
    !!   ∂⁵K =
    !!     8 * Sym(δδ r)   * K'''(s)
    !!   +16 * Sym(δ rrr)  * K''''(s)
    !!   +32 * rrrrr       * K'''''(s)
      !> Instance of Still kernel
      class(still_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in)  :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in)  :: rB(3)
      !> Born radius of atom A
      real(wp), intent(in)  :: bornA
      !> Born radius of atom B
      real(wp), intent(in)  :: bornB
      !> Fifth derivative of Still kernel
      real(wp), intent(out) :: d5(3, 3, 3, 3, 3)

      real(wp) :: rvec(3), s, a, e, u
      real(wp) :: invu, invsqrtu, um3
      real(wp) :: u1, u2, u3, u4
      real(wp) :: Acoef, Bcoef, Ccoef, Dcoef

      ! for K'''''(s)
      real(wp) :: um5, um7, um9, um11
      real(wp) :: u5, Ecoef
      integer :: i, j, k, l, m

      call still_scalar_coeffs(rA, rB, bornA, bornB, rvec, s, a, e, u, invu, invsqrtu, um3, &
                               u1, u2, u3, u4, Acoef, Bcoef, Ccoef, Dcoef)

      ! If helper short-circuited (a<=0 or s==0), everything is zero.
      if (a <= 0.0_wp .or. s == 0.0_wp) then
         d5 = 0.0_wp
         return
      end if

      ! Build u^(-p/2) ladder from invu and um3 = u^(-3/2)
      um5 = invu*um3
      um7 = invu*um5
      um9 = invu*um7
      um11 = invu*um9

      ! u(s) = s + a*exp(-s/(4a)), already have u4. Next derivative:
      u5 = -(1.0_wp/(1024.0_wp*a*a*a*a))*e

      ! K(s) = u(s)^(-1/2). 5th derivative w.r.t. s:
      ! Ecoef = K'''''(s)
      Ecoef = -0.5_wp*um3*u5 &
              +(15.0_wp/4.0_wp)*um5*(u1*u4+2.0_wp*u2*u3) &
              -(75.0_wp/8.0_wp)*um7*(2.0_wp*u1*u1*u3+3.0_wp*u1*u2*u2) &
              +(525.0_wp/8.0_wp)*um9*(u1*u1*u1*u2) &
              -(945.0_wp/32.0_wp)*um11*(u1*u1*u1*u1*u1)

      d5 = 0.0_wp

      do i = 1, 3
         do j = 1, 3
            do k = 1, 3
               do l = 1, 3
                  do m = 1, 3

                     ! 32 r_i r_j r_k r_l r_m * K'''''(s)
                     d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*32.0_wp* &
                                         rvec(i)*rvec(j)*rvec(k)*rvec(l)*rvec(m)*Ecoef

                     ! 16 * Sym(δ r r r) * K''''(s)  (10 permutations)
                     if (i == j) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*16.0_wp*rvec(k)*rvec(l)*rvec(m)*Dcoef
                     if (i == k) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*16.0_wp*rvec(j)*rvec(l)*rvec(m)*Dcoef
                     if (i == l) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*16.0_wp*rvec(j)*rvec(k)*rvec(m)*Dcoef
                     if (i == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*16.0_wp*rvec(j)*rvec(k)*rvec(l)*Dcoef
                     if (j == k) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*16.0_wp*rvec(i)*rvec(l)*rvec(m)*Dcoef
                     if (j == l) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*16.0_wp*rvec(i)*rvec(k)*rvec(m)*Dcoef
                     if (j == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*16.0_wp*rvec(i)*rvec(k)*rvec(l)*Dcoef
                     if (k == l) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*16.0_wp*rvec(i)*rvec(j)*rvec(m)*Dcoef
                     if (k == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*16.0_wp*rvec(i)*rvec(j)*rvec(l)*Dcoef
                     if (l == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*16.0_wp*rvec(i)*rvec(j)*rvec(k)*Dcoef

                     ! 8 * Sym(δδ r) * K'''(s)  (15 permutations)
                     if (i == j .and. k == l) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*8.0_wp*rvec(m)*Ccoef
                     if (i == j .and. k == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*8.0_wp*rvec(l)*Ccoef
                     if (i == j .and. l == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*8.0_wp*rvec(k)*Ccoef

                     if (i == k .and. j == l) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*8.0_wp*rvec(m)*Ccoef
                     if (i == k .and. j == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*8.0_wp*rvec(l)*Ccoef
                     if (i == k .and. l == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*8.0_wp*rvec(j)*Ccoef

                     if (i == l .and. j == k) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*8.0_wp*rvec(m)*Ccoef
                     if (i == l .and. j == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*8.0_wp*rvec(k)*Ccoef
                     if (i == l .and. k == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*8.0_wp*rvec(j)*Ccoef

                     if (i == m .and. j == k) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*8.0_wp*rvec(l)*Ccoef
                     if (i == m .and. j == l) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*8.0_wp*rvec(k)*Ccoef
                     if (i == m .and. k == l) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*8.0_wp*rvec(j)*Ccoef

                     if (j == k .and. l == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*8.0_wp*rvec(i)*Ccoef
                     if (j == l .and. k == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*8.0_wp*rvec(i)*Ccoef
                     if (j == m .and. k == l) d5(i, j, k, l, m) = d5(i, j, k, l, m)+self%keps*8.0_wp*rvec(i)*Ccoef

                  end do
               end do
            end do
         end do
      end do
   end subroutine still_d5Kdr5

   ! ---------------- internal helper: compute scalar coefficients K'(s),K''(s),K'''(s),K''''(s) ----------------

  pure subroutine still_scalar_coeffs(rA, rB, bornA, bornB, rvec, s, a, e, u, invu, invsqrtu, um3, u1, u2, u3, u4, Acoef, Bcoef, Ccoef, Dcoef)
      real(wp), intent(in)  :: rA(3), rB(3), bornA, bornB
      real(wp), intent(out) :: rvec(3), s, a, e, u
      real(wp), intent(out) :: invu, invsqrtu, um3
      real(wp), intent(out) :: u1, u2, u3, u4
      real(wp), intent(out) :: Acoef, Bcoef, Ccoef, Dcoef

      real(wp) :: um5, um7, um9
      real(wp) :: u1sq, u1cu, u1qu

      rvec = rA-rB
      s = dot_product(rvec, rvec)

      a = bornA*bornB
      if (a <= 0.0_wp .or. s == 0.0_wp) then
         e = 0.0_wp; u = 0.0_wp
         invu = 0.0_wp; invsqrtu = 0.0_wp
         um3 = 0.0_wp
         u1 = 0.0_wp; u2 = 0.0_wp; u3 = 0.0_wp; u4 = 0.0_wp
         Acoef = 0.0_wp; Bcoef = 0.0_wp; Ccoef = 0.0_wp; Dcoef = 0.0_wp
         return
      end if

      e = exp(-s/(4.0_wp*a))
      u = s+a*e

      invu = 1.0_wp/u
      invsqrtu = 1.0_wp/sqrt(u)

      ! u^(-3/2), u^(-5/2), u^(-7/2), u^(-9/2)
      um3 = invu*invsqrtu
      um5 = invu*invu*invsqrtu
      um7 = invu*invu*invu*invsqrtu
      um9 = invu*invu*invu*invu*invsqrtu

      ! derivatives of u(s) = s + a exp(-s/(4a)) w.r.t. s
      u1 = 1.0_wp-0.25_wp*e
      u2 = (1.0_wp/(16.0_wp*a))*e
      u3 = -(1.0_wp/(64.0_wp*a*a))*e
      u4 = (1.0_wp/(256.0_wp*a*a*a))*e

      u1sq = u1*u1
      u1cu = u1sq*u1
      u1qu = u1sq*u1sq

      ! K(s) = u(s)^(-1/2).  Acoef..Dcoef are K'(s)..K''''(s).
      Acoef = -0.5_wp*um3*u1

      Bcoef = 0.75_wp*um5*u1sq-0.5_wp*um3*u2

      Ccoef = (-15.0_wp/8.0_wp)*um7*u1cu+(9.0_wp/4.0_wp)*um5*(u1*u2)-0.5_wp*um3*u3

      Dcoef = (105.0_wp/16.0_wp)*um9*u1qu &
              -(45.0_wp/4.0_wp)*um7*(u1sq*u2) &
              +(9.0_wp/4.0_wp)*um5*(u2*u2) &
              +3.0_wp*um5*(u1*u3) &
              -0.5_wp*um3*u4

   end subroutine still_scalar_coeffs

   pure subroutine still_scalar_coeffs_da(s, a, e, u, u1, u2, u3, u4, dA_da, dB_da, dC_da, dD_da)
      real(wp), intent(in)  :: s, a, e, u, u1, u2, u3, u4
      real(wp), intent(out) :: dA_da, dB_da, dC_da, dD_da

      real(wp) :: invu, invsqrtu
      real(wp) :: um3, um5, um7, um9, um11
      real(wp) :: du_da, de_da
      real(wp) :: du1_da, du2_da, du3_da, du4_da
      real(wp) :: dum3_da, dum5_da, dum7_da, dum9_da
      real(wp) :: u1sq, u1cu, u1qu

      if (a <= 0.0_wp) then
         dA_da = 0.0_wp; dB_da = 0.0_wp; dC_da = 0.0_wp; dD_da = 0.0_wp
         return
      end if

      invu = 1.0_wp/u
      invsqrtu = 1.0_wp/sqrt(u)

      um3 = invu*invsqrtu
      um5 = invu*um3
      um7 = invu*um5
      um9 = invu*um7
      um11 = invu*um9

      ! de/da for e = exp(-s/(4a))
      de_da = e*(s/(4.0_wp*a*a))

      ! u = s + a e
      du_da = e+a*de_da   ! = e*(1 + s/(4a))

      ! d(u^{-p/2})/da
      dum3_da = -1.5_wp*um5*du_da
      dum5_da = -2.5_wp*um7*du_da
      dum7_da = -3.5_wp*um9*du_da
      dum9_da = -4.5_wp*um11*du_da

      ! u1..u4 dependence on a
      du1_da = -0.25_wp*de_da

      du2_da = -(1.0_wp/(16.0_wp*a*a))*e+(1.0_wp/(16.0_wp*a))*de_da
      du3_da = -(1.0_wp/(64.0_wp*a*a))*de_da+(2.0_wp/(64.0_wp*a*a*a))*e
      du4_da = -(3.0_wp/(256.0_wp*a**4))*e+(1.0_wp/(256.0_wp*a**3))*de_da

      u1sq = u1*u1
      u1cu = u1sq*u1
      u1qu = u1sq*u1sq

      ! A = -1/2 um3 u1
      dA_da = -0.5_wp*(dum3_da*u1+um3*du1_da)

      ! B = 3/4 um5 u1^2 - 1/2 um3 u2
      dB_da = 0.75_wp*(dum5_da*u1sq+um5*2.0_wp*u1*du1_da) &
              -0.5_wp*(dum3_da*u2+um3*du2_da)

      ! C = -15/8 um7 u1^3 + 9/4 um5 (u1 u2) - 1/2 um3 u3
      dC_da = (-15.0_wp/8.0_wp)*(dum7_da*u1cu+um7*3.0_wp*u1sq*du1_da) &
              +(9.0_wp/4.0_wp)*(dum5_da*(u1*u2)+um5*(du1_da*u2+u1*du2_da)) &
              -0.5_wp*(dum3_da*u3+um3*du3_da)

      ! D = 105/16 um9 u1^4
      !   - 45/4  um7 (u1^2 u2)
      !   + 9/4   um5 u2^2
      !   + 3     um5 (u1 u3)
      !   - 1/2   um3 u4
      dD_da = (105.0_wp/16.0_wp)*(dum9_da*u1qu+um9*4.0_wp*u1cu*du1_da) &
              -(45.0_wp/4.0_wp)*(dum7_da*(u1sq*u2)+um7*(2.0_wp*u1*du1_da*u2+u1sq*du2_da)) &
              +(9.0_wp/4.0_wp)*(dum5_da*(u2*u2)+um5*2.0_wp*u2*du2_da) &
              +3.0_wp*(dum5_da*(u1*u3)+um5*(du1_da*u3+u1*du3_da)) &
              -0.5_wp*(dum3_da*u4+um3*du4_da)
   end subroutine still_scalar_coeffs_da

end module tblite_solvation_kernel_still
