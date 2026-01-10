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

!> @file tblite/solvation/alpb.f90
!> Provides the generalized Born interaction kernels used in the ALPB and GBSA implicit solvation models.

module tblite_solvation_kernel
   use mctc_env, only : wp
   use tblite_blas, only : gemv
   implicit none
   private

   public :: kernel_type, new_kernel
   public :: still_kernel, p16_kernel
   public :: kernel_enum, kernel_enum_type
   public :: compute_kernel_dkdr  ! convenience dispatcher by enum

   type :: kernel_enum_type
      integer :: still = 1
      integer :: p16   = 2
   end type kernel_enum_type

   type(kernel_enum_type), parameter :: kernel_enum = kernel_enum_type()

   ! Abstract base class for kernel types
   type, abstract :: kernel_type
      real(wp) :: keps
   contains
      procedure(add_kernel_mat_interface),   deferred :: add_kernel_mat
      procedure(add_kernel_deriv_interface), deferred :: add_kernel_deriv
      procedure(compute_kernel_dkdr_interface), deferred :: compute_kernel_dkdr
   end type kernel_type

   abstract interface
      !> Add kernel contributions to interaction matrix
      subroutine add_kernel_mat_interface(self, nat, xyz, brad, amat)
         import :: kernel_type, wp
         class(kernel_type), intent(in) :: self
         integer, intent(in) :: nat
         real(wp), intent(in) :: xyz(:, :)
         real(wp), intent(in) :: brad(:)
         real(wp), intent(inout) :: amat(:, :)
      end subroutine add_kernel_mat_interface

      !> Add kernel derivative contributions to energy and gradient
      subroutine add_kernel_deriv_interface(self, nat, xyz, qat, brad, brdr, energy, gradient)
         import :: kernel_type, wp
         class(kernel_type), intent(in) :: self
         integer, intent(in) :: nat
         real(wp), intent(in) :: xyz(:, :)
         real(wp), intent(in) :: qat(:)
         real(wp), intent(in) :: brad(:)
         real(wp), contiguous, intent(in) :: brdr(:, :, :)
         real(wp), intent(out) :: energy
         real(wp), contiguous, intent(inout) :: gradient(:, :)
      end subroutine add_kernel_deriv_interface

      !> Full element-wise derivative tensor of the actual kernel matrix:
      !> dKdr(:, k, i, j) = ∂K_ij / ∂r_k(:)
      subroutine compute_kernel_dkdr_interface(self, nat, xyz, brad, brdr, dKdr)
         import :: kernel_type, wp
         class(kernel_type), intent(in) :: self
         integer, intent(in) :: nat
         real(wp), intent(in) :: xyz(:, :)
         real(wp), intent(in) :: brad(:)
         real(wp), contiguous, intent(in) :: brdr(:, :, :)
         real(wp), contiguous, intent(out) :: dKdr(:, :, :, :)
      end subroutine compute_kernel_dkdr_interface
   end interface

   type, extends(kernel_type) :: still_kernel
   contains
      procedure :: add_kernel_mat   => add_still_mat
      procedure :: add_kernel_deriv => add_still_deriv
      procedure :: compute_kernel_dkdr => compute_still_dkdr_full
   end type still_kernel

   type, extends(kernel_type) :: p16_kernel
   contains
      procedure :: add_kernel_mat   => add_p16_mat
      procedure :: add_kernel_deriv => add_p16_deriv
      procedure :: compute_kernel_dkdr => compute_p16_dkdr_full
   end type p16_kernel

   real(wp), parameter :: zetaP16    = 1.028_wp
   real(wp), parameter :: zetaP16o16 = zetaP16 / 16.0_wp

contains

function new_kernel(kernel_id, keps) result(kernel)
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   class(kernel_type), allocatable :: kernel

   select case(kernel_id)
   case(kernel_enum%still)
      allocate(still_kernel :: kernel)
   case(kernel_enum%p16)
      allocate(p16_kernel :: kernel)
   case default
      allocate(p16_kernel :: kernel)
   end select

   kernel%keps = keps
end function new_kernel

!> Convenience dispatcher (switches by kernel enum, returns full 4D derivative tensor)
subroutine compute_kernel_dkdr(kernel_id, keps, nat, xyz, brad, brdr, dKdr)
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: brad(:)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)
   real(wp), contiguous, intent(out) :: dKdr(:, :, :, :)

   class(kernel_type), allocatable :: kernel

   kernel = new_kernel(kernel_id, keps)
   call kernel%compute_kernel_dkdr(nat, xyz, brad, brdr, dKdr)
end subroutine compute_kernel_dkdr

!==============================================================================
! P16 kernel
!==============================================================================

subroutine add_p16_mat(self, nat, xyz, brad, Amat)
   class(p16_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: brad(:)
   real(wp), intent(inout) :: Amat(:, :)

   integer :: iat, jat
   real(wp) :: r1, ab, arg, fgb, dfgb, bp, vec(3)

   do iat = 1, nat
      do jat = 1, iat - 1
         vec(:) = xyz(:, iat) - xyz(:, jat)
         r1 = norm2(vec)

         ab  = sqrt(brad(iat) * brad(jat))
         arg = ab / (ab + zetaP16o16*r1)
         arg = arg * arg
         arg = arg * arg
         arg = arg * arg
         arg = arg * arg
         fgb  = r1 + ab*arg
         dfgb = 1.0_wp / fgb

         Amat(iat, jat) = self%keps*dfgb + Amat(iat, jat)
         Amat(jat, iat) = self%keps*dfgb + Amat(jat, iat)
      enddo
      bp = 1.0_wp/brad(iat)
      Amat(iat, iat) = Amat(iat, iat) + self%keps*bp
   enddo
end subroutine add_p16_mat

subroutine add_p16_deriv(self, nat, xyz, qat, brad, brdr, energy, gradient)
   class(p16_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: qat(:)
   real(wp), intent(in) :: brad(:)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)
   real(wp), intent(out) :: energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)

   integer :: iat, jat
   real(wp) :: vec(3), r2, r1, ab, arg1, arg16, qq, fgb, dfgb, dfgb2, egb
   real(wp) :: dEdbri, dEdbrj, dG(3), ap, bp
   real(wp), allocatable :: dEdbr(:)

   allocate(dEdbr(nat), source = 0.0_wp)

   egb = 0.0_wp
   dEdbr(:) = 0.0_wp

   do iat = 1, nat
      do jat = 1, iat - 1
         vec(:) = xyz(:, iat) - xyz(:, jat)
         r1 = norm2(vec)
         r2 = r1*r1

         qq = qat(iat)*qat(jat)

         ab = sqrt(brad(iat) * brad(jat))
         arg1  = ab / (ab + zetaP16o16*r1)
         arg16 = arg1 * arg1
         arg16 = arg16 * arg16
         arg16 = arg16 * arg16
         arg16 = arg16 * arg16

         fgb   = r1 + ab*arg16
         dfgb  = 1.0_wp / fgb
         dfgb2 = dfgb * dfgb

         egb = egb + qq*self%keps*dfgb

         ap = (1.0_wp - zetaP16 * arg1 * arg16) * dfgb2
         dG(:) = ap * vec * self%keps / r1 * qq
         gradient(:, iat) = gradient(:, iat) - dG
         gradient(:, jat) = gradient(:, jat) + dG

         bp = -0.5_wp*(r1 * zetaP16 / ab * arg1 + 1.0_wp) / ab * arg16 * dfgb2
         dEdbri = brad(jat) * bp * self%keps * qq
         dEdbrj = brad(iat) * bp * self%keps * qq
         dEdbr(iat) = dEdbr(iat) + dEdbri
         dEdbr(jat) = dEdbr(jat) + dEdbrj
      end do

      bp = 1.0_wp/brad(iat)
      qq = qat(iat)*bp
      egb = egb + 0.5_wp*qat(iat)*qq*self%keps
      dEdbri = -0.5_wp*self%keps*qq*bp
      dEdbr(iat) = dEdbr(iat) + dEdbri*qat(iat)
   enddo

   call gemv(brdr, dEdbr, gradient, beta=1.0_wp)
   energy = egb
end subroutine add_p16_deriv

!> Full dKdr for the actual P16 kernel matrix (including diagonal self term)
subroutine compute_p16_dkdr_full(self, nat, xyz, brad, brdr, dKdr)
   class(p16_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                 ! (3,nat)
   real(wp), intent(in) :: brad(:)                   ! (nat)
   real(wp), contiguous, intent(in) :: brdr(:, :, :) ! (3,nat,nat)
   real(wp), contiguous, intent(out) :: dKdr(:, :, :, :) ! (3,nat,nat,nat)

   integer :: i, j, k
   real(wp) :: rvec(3), r, r2
   real(wp) :: ai, aj, ab
   real(wp) :: a1, a16
   real(wp) :: fgb, invfgb2
   real(wp) :: coef_pos
   real(wp) :: bp, dK_dai, dK_daj
   real(wp), parameter :: tiny_r = 1.0e-14_wp

   dKdr(:, :, :, :) = 0.0_wp

   ! Off-diagonal (compute i>j, then mirror)
   do i = 1, nat
      ai = brad(i)
      do j = 1, i-1
         aj = brad(j)

         rvec(:) = xyz(:, i) - xyz(:, j)
         r2      = dot_product(rvec, rvec)
         r       = sqrt(r2)

         ab = sqrt(ai * aj)

         a1  = ab / (ab + zetaP16o16 * r)
         a16 = a1 * a1
         a16 = a16 * a16
         a16 = a16 * a16
         a16 = a16 * a16

         fgb     = r + ab * a16
         invfgb2 = 1.0_wp / (fgb * fgb)

         ! Explicit coordinate part (radii held fixed):
         ! ∂K/∂r_i = -keps*(1 - zeta*a1*a16)/fgb^2 * rvec/r
         if (r > tiny_r) then
            coef_pos = -self%keps * (1.0_wp - zetaP16 * a1 * a16) * invfgb2 / r
            dKdr(:, i, i, j) = dKdr(:, i, i, j) + coef_pos * rvec(:)  ! k=i
            dKdr(:, j, i, j) = dKdr(:, j, i, j) - coef_pos * rvec(:)  ! k=j
         end if

         ! Born radii partials for chain rule:
         ! bp as in tblite add_p16_deriv (but without qq):
         bp = -0.5_wp * ( (r * zetaP16 / ab) * a1 + 1.0_wp ) / ab * a16 * invfgb2
         dK_dai = self%keps * aj * bp
         dK_daj = self%keps * ai * bp

         do k = 1, nat
            dKdr(:, k, i, j) = dKdr(:, k, i, j) &
               + dK_dai * brdr(:, k, i) &
               + dK_daj * brdr(:, k, j)
         end do

         ! Mirror symmetry
         dKdr(:, :, j, i) = dKdr(:, :, i, j)
      end do
   end do

   ! Diagonal self terms: K_ii = keps/a_i
   do i = 1, nat
      do k = 1, nat
         dKdr(:, k, i, i) = dKdr(:, k, i, i) &
            + (-self%keps / (brad(i)*brad(i))) * brdr(:, k, i)
      end do
   end do
end subroutine compute_p16_dkdr_full

!==============================================================================
! Still kernel
!==============================================================================

pure subroutine add_still_mat(self, nat, xyz, brad, Amat)
   class(still_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: brad(:)
   real(wp), intent(inout) :: Amat(:, :)

   integer  :: i, j
   real(wp), parameter :: a4 = 0.25_wp
   real(wp) :: aa, vec(3), r1, r2, bp
   real(wp) :: dd, expd, fgb2, dfgb

   do i = 1, nat
      do j = 1, i - 1
         vec(:) = xyz(:, i) - xyz(:, j)
         r1 = norm2(vec)
         r2 = r1*r1

         aa   = brad(i)*brad(j)
         dd   = a4*r2/aa
         expd = exp(-dd)
         fgb2 = r2 + aa*expd
         dfgb = 1.0_wp/sqrt(fgb2)

         Amat(i, j) = self%keps*dfgb + Amat(i, j)
         Amat(j, i) = self%keps*dfgb + Amat(j, i)
      enddo

      bp = 1.0_wp/brad(i)
      Amat(i, i) = Amat(i, i) + self%keps*bp
   enddo
end subroutine add_still_mat

subroutine add_still_deriv(self, nat, xyz, qat, brad, brdr, energy, gradient)
   class(still_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: qat(:)
   real(wp), intent(in) :: brad(:)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)
   real(wp), intent(out) :: energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)

   integer :: i, j
   real(wp), parameter :: a4 = 0.25_wp
   real(wp) :: aa, r2, fgb2
   real(wp) :: qq, dd, expd, dfgb, dfgb2, dfgb3, egb, ap, bp
   real(wp) :: grddbi, grddbj
   real(wp) :: dr(3), r1, vec(3)
   real(wp), allocatable :: grddb(:)

   allocate(grddb(nat), source = 0.0_wp)

   egb = 0.0_wp
   grddb(:) = 0.0_wp

   do i = 1, nat
      do j = 1, i - 1
         vec(:) = xyz(:, i) - xyz(:, j)
         r1 = norm2(vec)
         r2 = r1*r1

         qq = qat(i)*qat(j)
         aa   = brad(i)*brad(j)
         dd   = a4*r2/aa
         expd = exp(-dd)
         fgb2 = r2 + aa*expd
         dfgb2 = 1.0_wp/fgb2
         dfgb  = sqrt(dfgb2)
         dfgb3 = dfgb2*dfgb*self%keps

         egb = egb + qq*self%keps*dfgb

         ap = (1.0_wp - a4*expd)*dfgb3
         dr = ap*vec
         gradient(:, i) = gradient(:, i) - dr*qq
         gradient(:, j) = gradient(:, j) + dr*qq

         bp = -0.5_wp*expd*(1.0_wp+dd)*dfgb3
         grddbi = brad(j)*bp
         grddbj = brad(i)*bp
         grddb(i) = grddb(i) + grddbi*qq
         grddb(j) = grddb(j) + grddbj*qq
      enddo

      bp = 1.0_wp/brad(i)
      qq = qat(i)*bp
      egb = egb + 0.5_wp*qat(i)*qq*self%keps
      grddbi = -0.5_wp*self%keps*qq*bp
      grddb(i) = grddb(i) + grddbi*qat(i)
   enddo

   call gemv(brdr, grddb, gradient, beta=1.0_wp)
   energy = egb
end subroutine add_still_deriv

!> Full dKdr for the actual Still kernel matrix (including diagonal self term)
subroutine compute_still_dkdr_full(self, nat, xyz, brad, brdr, dKdr)
   class(still_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                 ! (3,nat)
   real(wp), intent(in) :: brad(:)                   ! (nat)
   real(wp), contiguous, intent(in) :: brdr(:, :, :) ! (3,nat,nat)
   real(wp), contiguous, intent(out) :: dKdr(:, :, :, :) ! (3,nat,nat,nat)

   integer :: i, j, k
   real(wp), parameter :: a4 = 0.25_wp
   real(wp) :: rvec(3), r2
   real(wp) :: A, d, E, f2, invf, invf3
   real(wp) :: pref_pos
   real(wp) :: dK_dai, dK_daj

   dKdr(:, :, :, :) = 0.0_wp

   ! Off-diagonal (compute i>j, then mirror)
   do i = 1, nat
      do j = 1, i-1
         rvec(:) = xyz(:, i) - xyz(:, j)
         r2      = dot_product(rvec, rvec)

         A = brad(i) * brad(j)
         d = a4 * r2 / A
         E = exp(-d)

         f2    = r2 + A * E
         invf  = 1.0_wp / sqrt(f2)
         invf3 = invf * invf * invf

         ! Explicit coordinate part (radii held fixed):
         ! ∂K/∂r_i = -keps*(1 - 1/4*E) * rvec / f^3
         pref_pos = -(1.0_wp - a4*E) * (self%keps * invf3)

         dKdr(:, i, i, j) = dKdr(:, i, i, j) + pref_pos * rvec(:)  ! k=i
         dKdr(:, j, i, j) = dKdr(:, j, i, j) - pref_pos * rvec(:)  ! k=j

         ! Born radii partials for chain rule:
         dK_dai = -0.5_wp * self%keps * E * (1.0_wp + d) * brad(j) * invf3
         dK_daj = -0.5_wp * self%keps * E * (1.0_wp + d) * brad(i) * invf3

         do k = 1, nat
            dKdr(:, k, i, j) = dKdr(:, k, i, j) &
               + dK_dai * brdr(:, k, i) &
               + dK_daj * brdr(:, k, j)
         end do

         ! Mirror symmetry
         dKdr(:, :, j, i) = dKdr(:, :, i, j)
      end do
   end do

   ! Diagonal self terms: K_ii = keps/a_i
   do i = 1, nat
      do k = 1, nat
         dKdr(:, k, i, i) = dKdr(:, k, i, i) &
            + (-self%keps / (brad(i)*brad(i))) * brdr(:, k, i)
      end do
   end do
end subroutine compute_still_dkdr_full

end module tblite_solvation_kernel
