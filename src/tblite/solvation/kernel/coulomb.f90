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

!> @file tblite/solvation/kernel/coulomb.f90
!> Provides the Coulomb kernel implementation.

module tblite_solvation_kernel_coulomb
   use mctc_env, only: wp
   use tblite_blas, only: gemv
   use tblite_solvation_kernel_type, only: kernel_type
   implicit none
   private

   public :: coulomb_kernel

   type, extends(kernel_type) :: coulomb_kernel
   contains
      procedure :: add_kernel_mat => add_coulomb_mat
      procedure :: add_kernel_deriv => add_coulomb_deriv
      procedure :: add_kernel_deriv_multipole_contributions => add_coulomb_deriv_multipole_contributions
      procedure :: kernel_d1_pair => coulomb_d1_pair
      procedure :: kernel_d1_pair_dborn => coulomb_d1_pair_dborn
      procedure :: kernel_d2_pair => coulomb_d2_pair
      procedure :: kernel_d2_pair_dborn => coulomb_d2_pair_dborn
      procedure :: kernel_d3_pair => coulomb_d3_pair
      procedure :: kernel_d3_pair_dborn => coulomb_d3_pair_dborn
      procedure :: kernel_d4_pair => coulomb_d4_pair
      procedure :: kernel_d4_pair_dborn => coulomb_d4_pair_dborn
      procedure :: kernel_d5_pair => coulomb_d5_pair
   end type

contains

   pure subroutine add_coulomb_mat(self, nat, xyz, brad, Amat)
      !> Instance of Coulomb kernel
      class(coulomb_kernel), intent(in) :: self
      !> Number of atoms
      integer, intent(in) :: nat
      !> Cartesian coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Born radii (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in) :: brad(:)
      !> Charge-charge interaction matrix
      real(wp), intent(inout) :: Amat(:, :)

      integer  :: i, j
      real(wp) :: vec(3), r1, invr

      Amat = 0.0_wp

      ! Classic Coulomb interaction kernel: 1 / R
      do i = 1, nat
         do j = 1, i-1
            vec(:) = xyz(:, i)-xyz(:, j)
            r1 = norm2(vec)

            invr = 1.0_wp/r1

            Amat(i, j) = Amat(i, j)+self%keps*invr
            Amat(j, i) = Amat(j, i)+self%keps*invr
         end do

      end do
   end subroutine add_coulomb_mat

   subroutine add_coulomb_deriv(self, nat, xyz, qat, brad, brdr, energy, gradient)
      !> Instance of Coulomb kernel
      class(coulomb_kernel), intent(in) :: self
      !> Number of atoms
      integer, intent(in) :: nat
      !> Cartesian coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Atomic partial charges
      real(wp), intent(in) :: qat(:)
      !> Born radii (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in) :: brad(:)
      !> Born radii derivatives (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), contiguous, intent(in) :: brdr(:, :, :)
      !> Solvation energy
      real(wp), intent(out) :: energy
      !> Molecular gradient
      real(wp), contiguous, intent(inout) :: gradient(:, :)

      integer :: i, j
      real(wp) :: vec(3), r1, r2, invr, invr3
      real(wp) :: qq
      real(wp) :: dr(3)
      real(wp) :: e_coul
      real(wp), allocatable :: grddb(:)

      ! Keep for interface compatibility (unused for pure Coulomb kernel)
      allocate (grddb(nat), source=0.0_wp)

      e_coul = 0.0_wp
      grddb(:) = 0.0_wp

      do i = 1, nat
         do j = 1, i-1
            vec(:) = xyz(:, i)-xyz(:, j)
            r1 = norm2(vec)
            r2 = r1*r1

            invr = 1.0_wp/r1
            invr3 = invr/r2    ! = 1 / r^3

            qq = qat(i)*qat(j)

            ! Energy contribution: keps * q_i q_j / r_ij
            e_coul = e_coul+self%keps*qq*invr

            ! d/dr (1/r) = - r_vec / r^3
            dr = self%keps*invr3*vec

            ! Gradient on coordinates
            gradient(:, i) = gradient(:, i)-dr*qq
            gradient(:, j) = gradient(:, j)+dr*qq
         end do

      end do

      ! Keep call for interface compatibility; grddb is zero so this is a no-op.
      call gemv(brdr, grddb, gradient, beta=1.0_wp)

      energy = e_coul
   end subroutine add_coulomb_deriv

   subroutine add_coulomb_deriv_multipole_contributions(self, nat, xyz, q_at, mu_at, q_at2, brad, brdr, gradient)
      use mctc_env, only: wp
      use tblite_blas, only: gemv
      implicit none

      class(coulomb_kernel), intent(in) :: self
      integer, intent(in) :: nat
      real(wp), intent(in) :: xyz(:, :)          ! (3,nat)
      real(wp), intent(in) :: q_at(:)               ! (nat)
      real(wp), intent(in) :: mu_at(:, :)           ! (3,nat)
      real(wp), intent(in) :: q_at2(:, :)         ! (6,nat)  (lower-tri: xx,xy,yy,xz,yz,zz), moments NOT doubled
      real(wp), intent(in) :: brad(:)            ! (nat)
      real(wp), contiguous, intent(in) :: brdr(:, :, :)   ! (3,nat,nat) -> gemv-compatible like in add_still_deriv
      real(wp), contiguous, intent(inout) :: gradient(:, :) ! (3,nat)

      ! Empty subroutine - Coulomb kernel has no multipole contributions
   end subroutine add_coulomb_deriv_multipole_contributions

   subroutine coulomb_d1_pair(self, rA, rB, bornA, bornB, d1)
   !! d1(i) = ∂(1/r)/∂r_i (gradient w.r.t. r = rA-rB)
      !> Instance of Coulomb kernel
      class(coulomb_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in) :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in) :: rB(3)
      !> Born radius of atom A (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in) :: bornA
      !> Born radius of atom B (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in) :: bornB
      !> First derivative of Coulomb kernel
      real(wp), intent(out) :: d1(3)

      real(wp) :: rvec(3), r, invr, invr3

      rvec = rA-rB
      r = norm2(rvec)
      if (r == 0.0_wp) then
         d1 = 0.0_wp
         return
      end if

      invr = 1.0_wp/r
      invr3 = invr*invr*invr

      d1 = -rvec*invr3
   end subroutine coulomb_d1_pair

   subroutine coulomb_d1_pair_dborn(self, rA, rB, bornA, bornB, d1_bA, d1_bB)
      !> Instance of Coulomb kernel
      class(coulomb_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in) :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in) :: rB(3)
      !> Born radius of atom A (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in) :: bornA
      !> Born radius of atom B (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in) :: bornB
      !> Derivative wrt Born radius of atom A (zero for Coulomb kernel)
      real(wp), intent(out) :: d1_bA(3)
      !> Derivative wrt Born radius of atom B (zero for Coulomb kernel)
      real(wp), intent(out) :: d1_bB(3)

      d1_bA = 0.0_wp
      d1_bB = 0.0_wp
   end subroutine coulomb_d1_pair_dborn

   subroutine coulomb_d2_pair(self, rA, rB, bornA, bornB, d2)
    !! d2(i,j) = ∂²(1/r)/∂r_i∂r_j  (Hessian)
      !> Instance of Coulomb kernel
      class(coulomb_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in)  :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in)  :: rB(3)
      !> Born radius of atom A (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in)  :: bornA
      !> Born radius of atom B (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in)  :: bornB
      !> Second derivative of Coulomb kernel
      real(wp), intent(out) :: d2(3, 3)
      real(wp) :: rvec(3), r, invr, invr3, invr5
      integer :: i, j

      rvec = rA-rB
      r = norm2(rvec)
      if (r == 0.0_wp) then
         d2 = 0.0_wp
         return
      end if

      invr = 1.0_wp/r
      invr3 = invr**3
      invr5 = invr3*invr**2

      do i = 1, 3
         do j = 1, 3
            d2(i, j) = 3.0_wp*rvec(i)*rvec(j)*invr5
            if (i == j) d2(i, j) = d2(i, j)-invr3
         end do
      end do
   end subroutine coulomb_d2_pair

   subroutine coulomb_d2_pair_dborn(self, rA, rB, bornA, bornB, d2_bA, d2_bB)
      !> Instance of Coulomb kernel
      class(coulomb_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in)  :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in)  :: rB(3)
      !> Born radius of atom A (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in)  :: bornA
      !> Born radius of atom B (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in)  :: bornB
      !> Derivative wrt Born radius of atom A (zero for Coulomb kernel)
      real(wp), intent(out) :: d2_bA(3, 3)
      !> Derivative wrt Born radius of atom B (zero for Coulomb kernel)
      real(wp), intent(out) :: d2_bB(3, 3)
      d2_bA = 0.0_wp
      d2_bB = 0.0_wp
   end subroutine coulomb_d2_pair_dborn

   subroutine coulomb_d3_pair(self, rA, rB, bornA, bornB, d3)
    !! d3(i,j,k) = ∂³(1/r)/∂r_i∂r_j∂r_k
    !!
    !! d3 = 3(δ_ij r_k + δ_ik r_j + δ_jk r_i)/r^5 - 15 r_i r_j r_k / r^7
      !> Instance of Coulomb kernel
      class(coulomb_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in)  :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in)  :: rB(3)
      !> Born radius of atom A (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in)  :: bornA
      !> Born radius of atom B (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in)  :: bornB
      !> Third derivative of Coulomb kernel
      real(wp), intent(out) :: d3(3, 3, 3)
      real(wp) :: rvec(3), r, invr, invr5, invr7
      integer :: i, j, k

      rvec = rA-rB
      r = norm2(rvec)
      if (r == 0.0_wp) then
         d3 = 0.0_wp
         return
      end if

      invr = 1.0_wp/r
      invr5 = invr**5
      invr7 = invr**7

      do i = 1, 3
         do j = 1, 3
            do k = 1, 3
               d3(i, j, k) = -15.0_wp*rvec(i)*rvec(j)*rvec(k)*invr7
               if (i == j) d3(i, j, k) = d3(i, j, k)+3.0_wp*rvec(k)*invr5
               if (i == k) d3(i, j, k) = d3(i, j, k)+3.0_wp*rvec(j)*invr5
               if (j == k) d3(i, j, k) = d3(i, j, k)+3.0_wp*rvec(i)*invr5
            end do
         end do
      end do
   end subroutine coulomb_d3_pair

   subroutine coulomb_d3_pair_dborn(self, rA, rB, bornA, bornB, d3_bA, d3_bB)
      !> Instance of Coulomb kernel
      class(coulomb_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in)  :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in)  :: rB(3)
      !> Born radius of atom A (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in)  :: bornA
      !> Born radius of atom B (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in)  :: bornB
      !> Derivative wrt Born radius of atom A (zero for Coulomb kernel)
      real(wp), intent(out) :: d3_bA(3, 3, 3)
      !> Derivative wrt Born radius of atom B (zero for Coulomb kernel)
      real(wp), intent(out) :: d3_bB(3, 3, 3)
      d3_bA = 0.0_wp
      d3_bB = 0.0_wp
   end subroutine coulomb_d3_pair_dborn

   subroutine coulomb_d4_pair(self, rA, rB, bornA, bornB, d4)
    !! d4(i,j,k,l) = ∂⁴(1/r)/∂r_i∂r_j∂r_k∂r_l
    !!
    !! d4 = 105 r_i r_j r_k r_l / r^9
    !!    - 15/r^7 * sum_pairs δ_(ab) r_c r_d
    !!    +  3/r^5 * (δ_ij δ_kl + δ_ik δ_jl + δ_il δ_jk)
      !> Instance of Coulomb kernel
      class(coulomb_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in)  :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in)  :: rB(3)
      !> Born radius of atom A (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in)  :: bornA
      !> Born radius of atom B (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in)  :: bornB
      !> Fourth derivative of Coulomb kernel
      real(wp), intent(out) :: d4(3, 3, 3, 3)
      real(wp) :: rvec(3), r, invr, invr5, invr7, invr9
      integer :: i, j, k, l

      rvec = rA-rB
      r = norm2(rvec)
      if (r == 0.0_wp) then
         d4 = 0.0_wp
         return
      end if

      invr = 1.0_wp/r
      invr5 = invr**5
      invr7 = invr**7
      invr9 = invr**9

      do i = 1, 3
         do j = 1, 3
            do k = 1, 3
               do l = 1, 3
                  d4(i, j, k, l) = 105.0_wp*rvec(i)*rvec(j)*rvec(k)*rvec(l)*invr9

                  if (i == j) d4(i, j, k, l) = d4(i, j, k, l)-15.0_wp*rvec(k)*rvec(l)*invr7
                  if (i == k) d4(i, j, k, l) = d4(i, j, k, l)-15.0_wp*rvec(j)*rvec(l)*invr7
                  if (i == l) d4(i, j, k, l) = d4(i, j, k, l)-15.0_wp*rvec(j)*rvec(k)*invr7
                  if (j == k) d4(i, j, k, l) = d4(i, j, k, l)-15.0_wp*rvec(i)*rvec(l)*invr7
                  if (j == l) d4(i, j, k, l) = d4(i, j, k, l)-15.0_wp*rvec(i)*rvec(k)*invr7
                  if (k == l) d4(i, j, k, l) = d4(i, j, k, l)-15.0_wp*rvec(i)*rvec(j)*invr7

                  if (i == j .and. k == l) d4(i, j, k, l) = d4(i, j, k, l)+3.0_wp*invr5
                  if (i == k .and. j == l) d4(i, j, k, l) = d4(i, j, k, l)+3.0_wp*invr5
                  if (i == l .and. j == k) d4(i, j, k, l) = d4(i, j, k, l)+3.0_wp*invr5
               end do
            end do
         end do
      end do
   end subroutine coulomb_d4_pair

   subroutine coulomb_d4_pair_dborn(self, rA, rB, bornA, bornB, d4_bA, d4_bB)
      !> Instance of Coulomb kernel
      class(coulomb_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in)  :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in)  :: rB(3)
      !> Born radius of atom A (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in)  :: bornA
      !> Born radius of atom B (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in)  :: bornB
      !> Derivative wrt Born radius of atom A (zero for Coulomb kernel)
      real(wp), intent(out) :: d4_bA(3, 3, 3, 3)
      !> Derivative wrt Born radius of atom B (zero for Coulomb kernel)
      real(wp), intent(out) :: d4_bB(3, 3, 3, 3)
      d4_bA = 0.0_wp
      d4_bB = 0.0_wp
   end subroutine coulomb_d4_pair_dborn

   subroutine coulomb_d5_pair(self, rA, rB, bornA, bornB, d5)
    !! d5(i,j,k,l,m) = ∂⁵(1/r)/∂r_i∂r_j∂r_k∂r_l∂r_m
    !!
    !! d5 = -945 r_i r_j r_k r_l r_m / r^11
    !!    + 105/r^9 * sum_{one delta} (δ_(ab) r_c r_d r_e)
    !!    -  15/r^7 * sum_{two deltas} (δ_(ab) δ_(cd) r_e)
      !> Instance of Coulomb kernel
      class(coulomb_kernel), intent(in) :: self
      !> Cartesian coordinates of atom A
      real(wp), intent(in)  :: rA(3)
      !> Cartesian coordinates of atom B
      real(wp), intent(in)  :: rB(3)
      !> Born radius of atom A (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in)  :: bornA
      !> Born radius of atom B (unused for Coulomb kernel)
      !  Keep for interface compatibility
      real(wp), intent(in)  :: bornB
      !> Fifth derivative of Coulomb kernel
      real(wp), intent(out) :: d5(3, 3, 3, 3, 3)

      real(wp) :: rvec(3), r, invr, invr7, invr9, invr11
      integer :: i, j, k, l, m

      rvec = rA-rB
      r = norm2(rvec)
      if (r == 0.0_wp) then
         d5 = 0.0_wp
         return
      end if

      invr = 1.0_wp/r
      invr7 = invr**7
      invr9 = invr**9
      invr11 = invr**11

      do i = 1, 3
         do j = 1, 3
            do k = 1, 3
               do l = 1, 3
                  do m = 1, 3

                     ! Leading fully anisotropic term
                     d5(i, j, k, l, m) = -945.0_wp*rvec(i)*rvec(j)*rvec(k)*rvec(l)*rvec(m)*invr11

                     ! +105/r^9 * (single-delta) terms: 10 permutations
                     if (i == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)+105.0_wp*rvec(j)*rvec(k)*rvec(l)*invr9
                     if (j == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)+105.0_wp*rvec(i)*rvec(k)*rvec(l)*invr9
                     if (k == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)+105.0_wp*rvec(i)*rvec(j)*rvec(l)*invr9
                     if (l == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)+105.0_wp*rvec(i)*rvec(j)*rvec(k)*invr9

                     if (i == j) d5(i, j, k, l, m) = d5(i, j, k, l, m)+105.0_wp*rvec(k)*rvec(l)*rvec(m)*invr9
                     if (i == k) d5(i, j, k, l, m) = d5(i, j, k, l, m)+105.0_wp*rvec(j)*rvec(l)*rvec(m)*invr9
                     if (i == l) d5(i, j, k, l, m) = d5(i, j, k, l, m)+105.0_wp*rvec(j)*rvec(k)*rvec(m)*invr9
                     if (j == k) d5(i, j, k, l, m) = d5(i, j, k, l, m)+105.0_wp*rvec(i)*rvec(l)*rvec(m)*invr9
                     if (j == l) d5(i, j, k, l, m) = d5(i, j, k, l, m)+105.0_wp*rvec(i)*rvec(k)*rvec(m)*invr9
                     if (k == l) d5(i, j, k, l, m) = d5(i, j, k, l, m)+105.0_wp*rvec(i)*rvec(j)*rvec(m)*invr9

                     ! -15/r^7 * (double-delta) terms: 15 permutations
                     if (i == j .and. k == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)-15.0_wp*rvec(l)*invr7
                     if (i == j .and. l == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)-15.0_wp*rvec(k)*invr7

                     if (i == k .and. j == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)-15.0_wp*rvec(l)*invr7
                     if (i == k .and. l == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)-15.0_wp*rvec(j)*invr7

                     if (i == l .and. j == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)-15.0_wp*rvec(k)*invr7
                     if (i == l .and. k == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)-15.0_wp*rvec(j)*invr7

                     if (j == k .and. i == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)-15.0_wp*rvec(l)*invr7
                     if (j == k .and. l == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)-15.0_wp*rvec(i)*invr7

                     if (j == l .and. i == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)-15.0_wp*rvec(k)*invr7
                     if (j == l .and. k == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)-15.0_wp*rvec(i)*invr7

                     if (k == l .and. i == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)-15.0_wp*rvec(j)*invr7
                     if (k == l .and. j == m) d5(i, j, k, l, m) = d5(i, j, k, l, m)-15.0_wp*rvec(i)*invr7

                     if (i == j .and. k == l) d5(i, j, k, l, m) = d5(i, j, k, l, m)-15.0_wp*rvec(m)*invr7
                     if (i == k .and. j == l) d5(i, j, k, l, m) = d5(i, j, k, l, m)-15.0_wp*rvec(m)*invr7
                     if (i == l .and. j == k) d5(i, j, k, l, m) = d5(i, j, k, l, m)-15.0_wp*rvec(m)*invr7

                  end do
               end do
            end do
         end do
      end do
   end subroutine coulomb_d5_pair

end module tblite_solvation_kernel_coulomb
