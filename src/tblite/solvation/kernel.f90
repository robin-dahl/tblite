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
   public :: compute_kernel_dkdr_ij, compute_kernel_d2kdr2_ij, compute_kernel_d3kdr3_ij ! convenience dispatcher by enum

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
      procedure(compute_kernel_dkdr_ij_interface), deferred :: compute_kernel_dkdr_ij
      procedure(compute_kernel_d2kdr2_ij_interface), deferred :: compute_kernel_d2kdr2_ij
      procedure(compute_kernel_d3kdr3_ij_interface), deferred :: compute_kernel_d3kdr3_ij
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

      subroutine compute_kernel_dkdr_ij_interface(self, nat, xyz, brad, brdr, i, j, dKdr_ij)
         import :: kernel_type, wp
         class(kernel_type), intent(in) :: self
         integer, intent(in) :: nat
         real(wp), intent(in) :: xyz(:, :)
         real(wp), intent(in) :: brad(:)
         real(wp), contiguous, intent(in) :: brdr(:, :, :)
         integer, intent(in) :: i, j
         real(wp), contiguous, intent(out) :: dKdr_ij(:, :)   ! (3,nat)
      end subroutine compute_kernel_dkdr_ij_interface

      subroutine compute_kernel_d2kdr2_ij_interface(self, nat, xyz, brad, brdr, brdr2, i, j, d2Kdr2_ij)
         import :: kernel_type, wp
         class(kernel_type), intent(in) :: self
         integer, intent(in) :: nat
         real(wp), intent(in) :: xyz(:, :)
         real(wp), intent(in) :: brad(:)
         real(wp), contiguous, intent(in) :: brdr(:, :, :)
         real(wp), contiguous, intent(in) :: brdr2(:, :, :, :, :)
         integer, intent(in) :: i, j
         real(wp), contiguous, intent(out) :: d2Kdr2_ij(:, :, :, :)
      end subroutine compute_kernel_d2kdr2_ij_interface

      subroutine compute_kernel_d3kdr3_ij_interface(self, nat, xyz, brad, brdr, brdr2, brdr3, i, j, d3Kdr3_ij)
         import :: kernel_type, wp
         class(kernel_type), intent(in) :: self
         integer, intent(in) :: nat
         real(wp), intent(in) :: xyz(:, :)
         real(wp), intent(in) :: brad(:)
         real(wp), contiguous, intent(in) :: brdr(:, :, :)
         real(wp), contiguous, intent(in) :: brdr2(:, :, :, :, :)
         real(wp), contiguous, intent(in) :: brdr3(:, :, :, :, :, :, :)
         integer, intent(in) :: i, j
         real(wp), contiguous, intent(out) :: d3Kdr3_ij(:, :, :, :, :, :)
      end subroutine compute_kernel_d3kdr3_ij_interface
   end interface

   type, extends(kernel_type) :: still_kernel
   contains
      procedure :: add_kernel_mat   => add_still_mat
      procedure :: add_kernel_deriv => add_still_deriv
      procedure :: compute_kernel_dkdr_ij => compute_still_dkdr_ij
      procedure :: compute_kernel_d2kdr2_ij => compute_still_d2kdr2_ij
      procedure :: compute_kernel_d3kdr3_ij => compute_still_d3kdr3_ij
   end type still_kernel

   type, extends(kernel_type) :: p16_kernel
   contains
      procedure :: add_kernel_mat   => add_p16_mat
      procedure :: add_kernel_deriv => add_p16_deriv
      procedure :: compute_kernel_dkdr_ij => compute_p16_dkdr_ij
      procedure :: compute_kernel_d2kdr2_ij => compute_p16_d2kdr2_ij
      procedure :: compute_kernel_d3kdr3_ij => compute_p16_d3kdr3_ij
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

!> Convenience dispatcher (switches by kernel enum, returns full derivative tensor)
subroutine compute_kernel_dkdr_ij(kernel_id, keps, nat, xyz, brad, brdr, i, j, dKdr_ij)
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: brad(:)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)
   integer, intent(in) :: i, j
   real(wp), contiguous, intent(out) :: dKdr_ij(:, :)

   class(kernel_type), allocatable :: kernel

   kernel = new_kernel(kernel_id, keps)
   call kernel%compute_kernel_dkdr_ij(nat, xyz, brad, brdr, i, j, dKdr_ij)
end subroutine compute_kernel_dkdr_ij

!> Convenience dispatcher (switches by kernel enum, returns full derivative tensor)
subroutine compute_kernel_d2kdr2_ij(kernel_id, keps, nat, xyz, brad, brdr, brdr2, i, j, d2Kdr2_ij)
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: brad(:)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)
   real(wp), contiguous, intent(in) :: brdr2(:, :, :, :, :)
   integer, intent(in) :: i, j
   real(wp), contiguous, intent(out) :: d2Kdr2_ij(:, :, :, :)

   class(kernel_type), allocatable :: kernel

   kernel = new_kernel(kernel_id, keps)
   call kernel%compute_kernel_d2kdr2_ij(nat, xyz, brad, brdr, brdr2, i, j, d2Kdr2_ij)
end subroutine compute_kernel_d2kdr2_ij

subroutine compute_kernel_d3kdr3_ij(kernel_id, keps, nat, xyz, brad, brdr, brdr2, brdr3, i, j, d3Kdr3_ij)
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: brad(:)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)
   real(wp), contiguous, intent(in) :: brdr2(:, :, :, :, :)
   real(wp), contiguous, intent(in) :: brdr3(:, :, :, :, :, :, :)
   integer, intent(in) :: i, j   
   real(wp), contiguous, intent(out) :: d3Kdr3_ij(:, :, :, :, :, :)

   class(kernel_type), allocatable :: kernel

   kernel = new_kernel(kernel_id, keps)
   call kernel%compute_kernel_d3kdr3_ij(nat, xyz, brad, brdr, brdr2, brdr3, i, j, d3Kdr3_ij)
end subroutine compute_kernel_d3kdr3_ij



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

subroutine compute_p16_dkdr_ij(self, nat, xyz, brad, brdr, i, j, dKdr_ij)
   class(p16_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                 ! (3,nat)
   real(wp), intent(in) :: brad(:)                   ! (nat)
   real(wp), contiguous, intent(in) :: brdr(:, :, :) ! (3,nat,nat)
   integer, intent(in) :: i, j
   real(wp), contiguous, intent(out) :: dKdr_ij(:, :) ! (3,nat)

   integer :: k
   real(wp) :: rvec(3), r, r2
   real(wp) :: ai, aj, ab
   real(wp) :: a1, a16
   real(wp) :: fgb, invfgb2
   real(wp) :: coef_pos
   real(wp) :: bp, dK_dai, dK_daj
   real(wp), parameter :: tiny_r = 1.0e-14_wp

   dKdr_ij(:, :) = 0.0_wp

   ! Diagonal: K_ii = keps/a_i
   if (i == j) then
      do k = 1, nat
         dKdr_ij(:, k) = dKdr_ij(:, k) + (-self%keps / (brad(i)*brad(i))) * brdr(:, k, i)
      end do
      return
   end if

   ! Off-diagonal element (i,j)
   ai = brad(i)
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

   ! Explicit coordinate part (radii held fixed)
   if (r > tiny_r) then
      coef_pos = -self%keps * (1.0_wp - zetaP16 * a1 * a16) * invfgb2 / r
      dKdr_ij(:, i) = dKdr_ij(:, i) + coef_pos * rvec(:)  ! k=i
      dKdr_ij(:, j) = dKdr_ij(:, j) - coef_pos * rvec(:)  ! k=j
   end if

   ! Born radii partials for chain rule
   bp = -0.5_wp * ( (r * zetaP16 / ab) * a1 + 1.0_wp ) / ab * a16 * invfgb2
   dK_dai = self%keps * aj * bp
   dK_daj = self%keps * ai * bp

   do k = 1, nat
      dKdr_ij(:, k) = dKdr_ij(:, k) &
         + dK_dai * brdr(:, k, i) &
         + dK_daj * brdr(:, k, j)
   end do
end subroutine compute_p16_dkdr_ij


subroutine compute_p16_d2kdr2_ij(self, nat, xyz, brad, brdr, brdr2, i, j, d2Kdr2_ij)
   class(p16_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                         ! (3,nat)
   real(wp), intent(in) :: brad(:)                           ! (nat)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)         ! (3,nat,nat)
   real(wp), contiguous, intent(in) :: brdr2(:, :, :, :, :)  ! (3,nat,3,nat,nat)
   integer, intent(in) :: i, j
   real(wp), contiguous, intent(out) :: d2Kdr2_ij(:, :, :, :) ! (3,nat,3,nat)

   integer :: k, l
   integer :: delk, dell
   real(wp) :: v(3), r2, r, invr, invr2
   real(wp) :: ai, aj, u, t, cp16
   real(wp) :: a1, a16, a17
   real(wp) :: g, invg, invg2, invg3
   real(wp) :: gr, grr
   real(wp) :: Kr, Krr, C
   real(wp) :: gu, guu
   real(wp) :: u_ai, u_aj, u_aiai, u_ajaj, u_aiaj
   real(wp) :: g_ai, g_aj, g_aiai, g_ajaj, g_aiaj
   real(wp) :: gr_u, gr_ai, gr_aj
   real(wp) :: dK_dai, dK_daj
   real(wp) :: d2K_dai2, d2K_daj2, d2K_daida_j
   real(wp) :: Kr_ai, Kr_aj
   real(wp) :: dC_dai, dC_daj
   real(wp) :: I3(3,3), vv(3,3), Hvv(3,3), M(3,3)
   real(wp) :: dk_i(3), dk_j(3), dl_i(3), dl_j(3)
   real(wp) :: coef1, coef2
   real(wp), parameter :: tiny_r = 1.0e-14_wp

   ! identity
   I3 = 0.0_wp
   I3(1,1) = 1.0_wp; I3(2,2) = 1.0_wp; I3(3,3) = 1.0_wp

   d2Kdr2_ij(:, :, :, :) = 0.0_wp
   cp16 = zetaP16o16


   ! -------------------------
   ! Diagonal self terms: K_ii = keps / a_i
   ! -------------------------
   if (i == j) then
      ai = brad(i)

      coef1 = -self%keps / (ai*ai)
      coef2 =  2.0_wp * self%keps / (ai*ai*ai)

      do k = 1, nat
         dk_i = brdr(:, k, i)
         do l = 1, nat
            dl_i = brdr(:, l, i)

            M = coef2 * (spread(dk_i,2,3)*spread(dl_i,1,3)) + coef1 * brdr2(:, k, :, l, i)
            d2Kdr2_ij(:, k, :, l) = d2Kdr2_ij(:, k, :, l) + M
         end do
      end do
      return
   end if

   ! -------------------------
   ! Off-diagonal element (i,j)
   ! -------------------------
   ai = brad(i)
   aj = brad(j)

   v(:) = xyz(:, i) - xyz(:, j)
   r2   = dot_product(v, v)
   r    = sqrt(r2)
   if (r <= tiny_r) return

   invr  = 1.0_wp / r
   invr2 = invr * invr

   u = sqrt(ai * aj)
   t = u + cp16 * r

   ! a1 = u/t; a16 = a1^16; a17 = a1^17
   a1  = u / t
   a16 = a1 * a1
   a16 = a16 * a16
   a16 = a16 * a16
   a16 = a16 * a16
   a17 = a1 * a16

   g     = r + u * a16
   invg  = 1.0_wp / g
   invg2 = invg * invg
   invg3 = invg2 * invg

   ! r-derivatives (radii held fixed)
   gr  = 1.0_wp - zetaP16 * a17
   grr = 17.0_wp * zetaP16 * cp16 * a17 / t

   ! K_r and K_rr (radii held fixed)
   Kr  = -self%keps * gr * invg2
   Krr = -self%keps * ( grr * invg2 - 2.0_wp * (gr*gr) * invg3 )

   ! Gradient wrt v is (Kr/r)*v
   C = Kr * invr

   ! Hessian wrt v (radii held fixed)
   vv  = spread(v,2,3) * spread(v,1,3)
   Hvv = C * I3 + (Krr - C) * (vv * invr2)

   ! u-derivatives for radii chain terms (r fixed)
   gu  = a16 * (u + 17.0_wp*cp16*r) / t
   guu = 272.0_wp * cp16*cp16 * r2 * a16 / (u * t*t)

   ! u derivatives wrt ai/aj
   u_ai   = u / (2.0_wp * ai)
   u_aj   = u / (2.0_wp * aj)
   u_aiai = -u / (4.0_wp * ai*ai)
   u_ajaj = -u / (4.0_wp * aj*aj)
   u_aiaj =  u / (4.0_wp * ai*aj)

   ! g radii partials (r fixed)
   g_ai   = gu  * u_ai
   g_aj   = gu  * u_aj
   g_aiai = guu * u_ai*u_ai + gu * u_aiai
   g_ajaj = guu * u_aj*u_aj + gu * u_ajaj
   g_aiaj = guu * u_ai*u_aj + gu * u_aiaj

   ! First radii partials of K
   dK_dai = -self%keps * g_ai * invg2
   dK_daj = -self%keps * g_aj * invg2

   ! Second radii partials of K
   d2K_dai2     = -self%keps * ( g_aiai * invg2 - 2.0_wp * (g_ai*g_ai) * invg3 )
   d2K_daj2     = -self%keps * ( g_ajaj * invg2 - 2.0_wp * (g_aj*g_aj) * invg3 )
   d2K_daida_j  = -self%keps * ( g_aiaj * invg2 - 2.0_wp * (g_ai*g_aj) * invg3 )

   ! Need ∂C/∂a_i, ∂C/∂a_j for mixed v–a terms:
   gr_u  = -17.0_wp * zetaP16 * cp16 * r * a16 / (t*t)
   gr_ai = gr_u * u_ai
   gr_aj = gr_u * u_aj

   Kr_ai = -self%keps * ( gr_ai * invg2 - 2.0_wp * gr * g_ai * invg3 )
   Kr_aj = -self%keps * ( gr_aj * invg2 - 2.0_wp * gr * g_aj * invg3 )

   dC_dai = Kr_ai * invr
   dC_daj = Kr_aj * invr

   ! Assemble coordinate Hessian for this (i,j) only: d2K_ij / dr_k dr_l
   do k = 1, nat
      dk_i = brdr(:, k, i)
      dk_j = brdr(:, k, j)

      delk = 0
      if (k == i) delk = delk + 1
      if (k == j) delk = delk - 1

      do l = 1, nat
         dl_i = brdr(:, l, i)
         dl_j = brdr(:, l, j)

         dell = 0
         if (l == i) dell = dell + 1
         if (l == j) dell = dell - 1

         M = 0.0_wp

         ! (1) Explicit vv part
         if (delk /= 0 .and. dell /= 0) then
            M = M + real(delk*dell, wp) * Hvv
         end if

         ! (2) Mixed v–a parts
         if (delk /= 0) then
            M = M + real(delk, wp) * ( dC_dai * (spread(v,2,3)*spread(dl_i,1,3)) &
                                     + dC_daj * (spread(v,2,3)*spread(dl_j,1,3)) )
         end if

         if (dell /= 0) then
            M = M + real(dell, wp) * ( dC_dai * (spread(dk_i,2,3)*spread(v,1,3)) &
                                     + dC_daj * (spread(dk_j,2,3)*spread(v,1,3)) )
         end if

         ! (3) a–a parts
         M = M + d2K_dai2    * (spread(dk_i,2,3)*spread(dl_i,1,3))
         M = M + d2K_daj2    * (spread(dk_j,2,3)*spread(dl_j,1,3))
         M = M + d2K_daida_j * ( (spread(dk_i,2,3)*spread(dl_j,1,3)) &
                               + (spread(dk_j,2,3)*spread(dl_i,1,3)) )

         ! (4) brdr2 terms
         M = M + dK_dai * brdr2(:, k, :, l, i) + dK_daj * brdr2(:, k, :, l, j)

         d2Kdr2_ij(:, k, :, l) = d2Kdr2_ij(:, k, :, l) + M
      end do
   end do

end subroutine compute_p16_d2kdr2_ij



!> Full d3Kdr3 for the actual P16 kernel matrix (including diagonal self term)
subroutine compute_p16_d3kdr3_full(self, nat, xyz, brad, brdr, brdr2, brdr3, d3Kdr3)
   class(p16_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                         ! (3,nat)
   real(wp), intent(in) :: brad(:)                           ! (nat)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)         ! (3,nat,nat)
   real(wp), contiguous, intent(in) :: brdr2(:, :, :, :, :)  ! (3,nat,3,nat,nat)
   real(wp), contiguous, intent(in) :: brdr3(:, :, :, :, :, :, :)  ! (3,nat,3,nat,3,nat,nat)
   real(wp), contiguous, intent(out) :: d3Kdr3(:, :, :, :, :, :, :, :) ! (3,nat,3,nat,3,nat,nat,nat)

   integer :: i, j, k, l, m
   integer :: alpha, beta, gamma
   integer :: delk, dell, delm

   real(wp), parameter :: tiny_r = 1.0e-14_wp
   real(wp) :: c

   real(wp) :: v(3), r2, r, invr, invr2
   real(wp) :: ai, aj, u, t
   real(wp) :: a1, a16, a17
   real(wp) :: g, invg, invg2, invg3, invg4

   ! g-derivatives wrt r (radii fixed)
   real(wp) :: gr, grr, grrr

   ! g-derivatives wrt u (r fixed)
   real(wp) :: gu, guu, guuu

   ! mixed g_r,u etc (only need gr_u, gr_uu, grr_u)
   real(wp) :: gr_u, gr_uu, grr_u

   ! u-derivatives wrt ai,aj
   real(wp) :: u_ai, u_aj
   real(wp) :: u_aiai, u_ajaj, u_aiaj
   real(wp) :: u_aiaiai, u_ajajaj
   real(wp) :: u_aiaiaj, u_aiajaj

   ! convert g-derivatives to ai/aj derivatives
   real(wp) :: g_ai, g_aj
   real(wp) :: g_aiai, g_ajaj, g_aiaj
   real(wp) :: g_aiaiai, g_ajajaj, g_aiaiaj, g_aiajaj

   real(wp) :: gr_ai, gr_aj
   real(wp) :: gr_aiai, gr_ajaj, gr_aiaj
   real(wp) :: grr_ai, grr_aj

   ! f(g)=keps/g derivatives
   real(wp) :: f1, f2, f3

   ! scalar K partials (r, ai, aj)
   real(wp) :: Kr, Krr, Krrr
   real(wp) :: Kai, Kaj
   real(wp) :: Kaiai, Kajaj, Kaiaj
   real(wp) :: Kr_ai, Kr_aj
   real(wp) :: Krr_ai, Krr_aj
   real(wp) :: Kr_aiai, Kr_ajaj, Kr_aiaj
   real(wp) :: Kaiaiai, Kajajaj, Kaiaiaj, Kaiajaj

   ! radial tensors (v-derivatives holding radii fixed)
   real(wp) :: e(3)
   real(wp) :: I3(3,3)
   real(wp) :: A, B, Bp, Ap
   real(wp) :: Kvv_ai(3,3), Kvv_aj(3,3)
   real(wp) :: Kvvv(3,3,3)
   real(wp) :: Kv_ai(3), Kv_aj(3)
   real(wp) :: Kv_aiai(3), Kv_ajaj(3), Kv_aiaj(3)

   ! coordinate-to-radius scalars
   real(wp) :: dai_k, dai_l, dai_m
   real(wp) :: daj_k, daj_l, daj_m
   real(wp) :: d2ai_kl, d2ai_km, d2ai_lm
   real(wp) :: d2aj_kl, d2aj_km, d2aj_lm
   real(wp) :: d3ai_klm, d3aj_klm

   real(wp) :: term
   real(wp) :: self_f1, self_f2, self_f3

   I3 = 0.0_wp
   I3(1,1)=1.0_wp; I3(2,2)=1.0_wp; I3(3,3)=1.0_wp

   d3Kdr3(:, :, :, :, :, :, :, :) = 0.0_wp

   c = zetaP16o16   ! = zeta/16

   ! =========================
   ! Off-diagonal: i>j, mirror
   ! =========================
   do i = 1, nat
      ai = brad(i)
      do j = 1, i-1
         aj = brad(j)

         v(:) = xyz(:, i) - xyz(:, j)
         r2   = dot_product(v, v)
         r    = sqrt(r2)
         if (r <= tiny_r) cycle

         invr  = 1.0_wp / r
         invr2 = invr * invr
         e(:)  = v(:) * invr

         u = sqrt(ai*aj)
         t = u + c*r

         a1  = u / t
         a16 = a1*a1
         a16 = a16*a16
         a16 = a16*a16
         a16 = a16*a16
         a17 = a1 * a16

         g    = r + u * a16         ! = r + u^17 / t^16
         invg = 1.0_wp / g
         invg2 = invg*invg
         invg3 = invg2*invg
         invg4 = invg2*invg2

         ! --- g derivatives wrt r (radii fixed) ---
         gr   = 1.0_wp - zetaP16 * a17
         grr  = 17.0_wp * zetaP16 * c * a17 / t
         grrr = -4896.0_wp * c*c*c * a17 / (t*t)   ! -4896 c^3 a17 / t^2

         ! --- g derivatives wrt u (r fixed) ---
         ! gu  = a16*(u + 17 c r)/t
         gu   = a16 * (u + 17.0_wp*c*r) / t
         ! guu = 272 c^2 r^2 u^15 / t^18 = 272 c^2 r^2 * a16 / (u t^2)
         guu  = 272.0_wp * c*c * r2 * a16 / (u * t*t)
         ! guuu = 816 c^2 r^2 u^14 (5 c r - u)/t^19
         guuu = 816.0_wp * c*c * r2 * (u**14) * (5.0_wp*c*r - u) / (t**19)

         ! --- mixed derivatives needed for Kr_ai etc (r fixed in ai-derivs) ---
         ! gr_u  = -17 zeta c r * a16 / t^2
         gr_u  = -17.0_wp * zetaP16 * c * r * a16 / (t*t)
         ! gr_uu = -zeta * q_uu,  q_uu = 34 c r u^15 (8 c r - u)/t^19
         gr_uu = -zetaP16 * (34.0_wp * c * r * (u**15) * (8.0_wp*c*r - u) / (t**19))
         ! grr_u = 17 zeta c u^16 (17 c r - u)/t^19
         grr_u = 17.0_wp * zetaP16 * c * (u**16) * (17.0_wp*c*r - u) / (t**19)

         ! --- u derivatives wrt ai,aj ---
         u_ai   = u / (2.0_wp*ai)
         u_aj   = u / (2.0_wp*aj)
         u_aiai = -u / (4.0_wp*ai*ai)
         u_ajaj = -u / (4.0_wp*aj*aj)
         u_aiaj =  u / (4.0_wp*ai*aj)

         u_aiaiai =  3.0_wp*u / (8.0_wp*ai**3)
         u_ajajaj =  3.0_wp*u / (8.0_wp*aj**3)
         u_aiaiaj = -u / (8.0_wp*ai*ai*aj)
         u_aiajaj = -u / (8.0_wp*ai*aj*aj)

         ! --- g radii derivatives (r fixed) ---
         g_ai   = gu * u_ai
         g_aj   = gu * u_aj
         g_aiai = guu*u_ai*u_ai + gu*u_aiai
         g_ajaj = guu*u_aj*u_aj + gu*u_ajaj
         g_aiaj = guu*u_ai*u_aj + gu*u_aiaj

         g_aiaiai = guuu*u_ai**3 + 3.0_wp*guu*u_ai*u_aiai + gu*u_aiaiai
         g_ajajaj = guuu*u_aj**3 + 3.0_wp*guu*u_aj*u_ajaj + gu*u_ajajaj
         g_aiaiaj = guuu*u_ai*u_ai*u_aj + guu*u_aiai*u_aj + 2.0_wp*guu*u_ai*u_aiaj + gu*u_aiaiaj
         g_aiajaj = guuu*u_aj*u_aj*u_ai + guu*u_ajaj*u_ai + 2.0_wp*guu*u_aj*u_aiaj + gu*u_aiajaj

         ! --- gr radii derivatives (r fixed) ---
         gr_ai   = gr_u * u_ai
         gr_aj   = gr_u * u_aj
         gr_aiai = gr_uu*u_ai*u_ai + gr_u*u_aiai
         gr_ajaj = gr_uu*u_aj*u_aj + gr_u*u_ajaj
         gr_aiaj = gr_uu*u_ai*u_aj + gr_u*u_aiaj

         grr_ai = grr_u * u_ai
         grr_aj = grr_u * u_aj

         ! --- f(g)=keps/g derivatives ---
         f1 = -self%keps * invg2
         f2 =  2.0_wp * self%keps * invg3
         f3 = -6.0_wp * self%keps * invg4

         ! --- scalar K derivatives ---
         Kr   = f1 * gr
         Krr  = f2 * gr*gr + f1 * grr
         Krrr = f3 * gr*gr*gr + 3.0_wp*f2*gr*grr + f1*grrr

         Kai = f1 * g_ai
         Kaj = f1 * g_aj

         Kaiai = f2*g_ai*g_ai + f1*g_aiai
         Kajaj = f2*g_aj*g_aj + f1*g_ajaj
         Kaiaj = f2*g_ai*g_aj + f1*g_aiaj

         ! mixed with r
         Kr_ai = f2*g_ai*gr + f1*gr_ai
         Kr_aj = f2*g_aj*gr + f1*gr_aj

         Krr_ai = f3*g_ai*gr*gr + 2.0_wp*f2*gr*gr_ai + f2*g_ai*grr + f1*grr_ai
         Krr_aj = f3*g_aj*gr*gr + 2.0_wp*f2*gr*gr_aj + f2*g_aj*grr + f1*grr_aj

         Kr_aiai = f3*gr*(g_ai*g_ai) + f2*gr*g_aiai + 2.0_wp*f2*g_ai*gr_ai + f1*gr_aiai
         Kr_ajaj = f3*gr*(g_aj*g_aj) + f2*gr*g_ajaj + 2.0_wp*f2*g_aj*gr_aj + f1*gr_ajaj
         Kr_aiaj = f3*gr*(g_ai*g_aj) + f2*gr*g_aiaj + f2*(g_ai*gr_aj + g_aj*gr_ai) + f1*gr_aiaj

         ! radii-only third derivatives
         Kaiaiai = f3*g_ai**3 + 3.0_wp*f2*g_ai*g_aiai + f1*g_aiaiai
         Kajajaj = f3*g_aj**3 + 3.0_wp*f2*g_aj*g_ajaj + f1*g_ajajaj
         Kaiaiaj = f3*(g_ai*g_ai*g_aj) + f2*(g_aiai*g_aj + 2.0_wp*g_ai*g_aiaj) + f1*g_aiaiaj
         Kaiajaj = f3*(g_aj*g_aj*g_ai) + f2*(g_ajaj*g_ai + 2.0_wp*g_aj*g_aiaj) + f1*g_aiajaj

         ! --- build radial tensors in v-space (radii fixed) ---
         ! Hessian form uses A,B:
         A  = Krr - Kr*invr
         B  = Kr*invr
         Bp = (Krr*r - Kr) * invr2            ! d/dr (Kr/r)
         Ap = Krrr - Bp                        ! dA/dr = Krrr - d(Kr/r)/dr

         ! Kvv_ai / Kvv_aj (two v-derivatives + one radius)
         Kvv_ai(:,:) = 0.0_wp
         Kvv_aj(:,:) = 0.0_wp
         do alpha = 1,3
            do beta = 1,3
               Kvv_ai(alpha,beta) = (Krr_ai - Kr_ai*invr) * e(alpha)*e(beta) + (Kr_ai*invr) * I3(alpha,beta)
               Kvv_aj(alpha,beta) = (Krr_aj - Kr_aj*invr) * e(alpha)*e(beta) + (Kr_aj*invr) * I3(alpha,beta)
            end do
         end do

         ! Kv_ai etc (one v-derivative + radii derivatives)
         Kv_ai(:)   = Kr_ai   * e(:)
         Kv_aj(:)   = Kr_aj   * e(:)
         Kv_aiai(:) = Kr_aiai * e(:)
         Kv_ajaj(:) = Kr_ajaj * e(:)
         Kv_aiaj(:) = Kr_aiaj * e(:)

         ! third v-derivative tensor Kvvv
         Kvvv(:,:,:) = 0.0_wp
         do alpha = 1,3
            do beta = 1,3
               do gamma = 1,3
                  Kvvv(alpha,beta,gamma) = Ap * e(alpha)*e(beta)*e(gamma) &
                     + (A*invr) * ( I3(alpha,gamma)*e(beta) + I3(beta,gamma)*e(alpha) - 2.0_wp*e(alpha)*e(beta)*e(gamma) ) &
                     + Bp * I3(alpha,beta) * e(gamma)
               end do
            end do
         end do

         ! ==========================================================
         ! Assemble full coordinate third derivative for all k,l,m
         ! ==========================================================
         do k = 1, nat
            delk = 0; if (k==i) delk=delk+1; if (k==j) delk=delk-1
            do l = 1, nat
               dell = 0; if (l==i) dell=dell+1; if (l==j) dell=dell-1
               do m = 1, nat
                  delm = 0; if (m==i) delm=delm+1; if (m==j) delm=delm-1

                  do alpha = 1,3
                     dai_k = brdr(alpha,k,i); daj_k = brdr(alpha,k,j)
                     do beta = 1,3
                        dai_l = brdr(beta,l,i); daj_l = brdr(beta,l,j)
                        do gamma = 1,3
                           dai_m = brdr(gamma,m,i); daj_m = brdr(gamma,m,j)

                           d2ai_kl = brdr2(alpha,k,beta,l,i)
                           d2ai_km = brdr2(alpha,k,gamma,m,i)
                           d2ai_lm = brdr2(beta,l,gamma,m,i)

                           d2aj_kl = brdr2(alpha,k,beta,l,j)
                           d2aj_km = brdr2(alpha,k,gamma,m,j)
                           d2aj_lm = brdr2(beta,l,gamma,m,j)

                           d3ai_klm = brdr3(alpha,k,beta,l,gamma,m,i)
                           d3aj_klm = brdr3(alpha,k,beta,l,gamma,m,j)

                           term = 0.0_wp

                           ! (1) vvv term
                           if (delk/=0 .and. dell/=0 .and. delm/=0) then
                              term = term + real(delk*dell*delm,wp) * Kvvv(alpha,beta,gamma)
                           end if

                           ! (2) vv–a terms (three placements)
                           if (delk/=0 .and. dell/=0) then
                              term = term + real(delk*dell,wp) * ( Kvv_ai(alpha,beta)*dai_m + Kvv_aj(alpha,beta)*daj_m )
                           end if
                           if (delk/=0 .and. delm/=0) then
                              term = term + real(delk*delm,wp) * ( Kvv_ai(alpha,gamma)*dai_l + Kvv_aj(alpha,gamma)*daj_l )
                           end if
                           if (dell/=0 .and. delm/=0) then
                              term = term + real(dell*delm,wp) * ( Kvv_ai(beta,gamma)*dai_k + Kvv_aj(beta,gamma)*daj_k )
                           end if

                           ! (3) v–aa / v–bb / v–ab terms (three placements)
                           if (delk/=0) then
                              term = term + real(delk,wp) * ( Kv_aiai(alpha)*(dai_l*dai_m) + Kv_ajaj(alpha)*(daj_l*daj_m) &
                                 + Kv_aiaj(alpha)*(dai_l*daj_m + daj_l*dai_m) )
                           end if
                           if (dell/=0) then
                              term = term + real(dell,wp) * ( Kv_aiai(beta)*(dai_k*dai_m) + Kv_ajaj(beta)*(daj_k*daj_m) &
                                 + Kv_aiaj(beta)*(dai_k*daj_m + daj_k*dai_m) )
                           end if
                           if (delm/=0) then
                              term = term + real(delm,wp) * ( Kv_aiai(gamma)*(dai_k*dai_l) + Kv_ajaj(gamma)*(daj_k*daj_l) &
                                 + Kv_aiaj(gamma)*(dai_k*daj_l + daj_k*dai_l) )
                           end if

                           ! (4) radii-only cubic terms
                           term = term + Kaiaiai*(dai_k*dai_l*dai_m) + Kajajaj*(daj_k*daj_l*daj_m)
                           term = term + Kaiaiaj*(dai_k*dai_l*daj_m + dai_k*daj_l*dai_m + daj_k*dai_l*dai_m)
                           term = term + Kaiajaj*(daj_k*daj_l*dai_m + daj_k*dai_l*daj_m + dai_k*daj_l*daj_m)

                           ! (5) K_pq terms with brdr2 (three pairings)
                           ! (k,l) paired, m as remaining
                           if (delm/=0) then
                              term = term + real(delm,wp) * ( d2ai_kl * Kv_ai(gamma) + d2aj_kl * Kv_aj(gamma) )
                           end if
                           term = term + d2ai_kl * (Kaiai*dai_m + Kaiaj*daj_m) + d2aj_kl * (Kaiaj*dai_m + Kajaj*daj_m)

                           ! (k,m) paired, l remaining
                           if (dell/=0) then
                              term = term + real(dell,wp) * ( d2ai_km * Kv_ai(beta) + d2aj_km * Kv_aj(beta) )
                           end if
                           term = term + d2ai_km * (Kaiai*dai_l + Kaiaj*daj_l) + d2aj_km * (Kaiaj*dai_l + Kajaj*daj_l)

                           ! (l,m) paired, k remaining
                           if (delk/=0) then
                              term = term + real(delk,wp) * ( d2ai_lm * Kv_ai(alpha) + d2aj_lm * Kv_aj(alpha) )
                           end if
                           term = term + d2ai_lm * (Kaiai*dai_k + Kaiaj*daj_k) + d2aj_lm * (Kaiaj*dai_k + Kajaj*daj_k)

                           ! (6) K_p * brdr3 term
                           term = term + Kai*d3ai_klm + Kaj*d3aj_klm

                           d3Kdr3(alpha,k,beta,l,gamma,m,i,j) = d3Kdr3(alpha,k,beta,l,gamma,m,i,j) + term
                           d3Kdr3(alpha,k,beta,l,gamma,m,j,i) = d3Kdr3(alpha,k,beta,l,gamma,m,j,i) + term
                        end do
                     end do
                  end do

               end do
            end do
         end do

      end do
   end do

   ! =========================
   ! Diagonal self term: K_ii = keps / a_i
   ! =========================
   do i = 1, nat
      ai = brad(i)

      self_f1 = -self%keps / (ai*ai)
      self_f2 =  2.0_wp * self%keps / (ai*ai*ai)
      self_f3 = -6.0_wp * self%keps / (ai**4)

      do k = 1, nat
         do l = 1, nat
            do m = 1, nat
               do alpha = 1,3
                  dai_k = brdr(alpha,k,i)
                  do beta = 1,3
                     dai_l   = brdr(beta,l,i)
                     d2ai_kl = brdr2(alpha,k,beta,l,i)
                     do gamma = 1,3
                        dai_m   = brdr(gamma,m,i)
                        d2ai_km = brdr2(alpha,k,gamma,m,i)
                        d2ai_lm = brdr2(beta,l,gamma,m,i)
                        d3ai_klm = brdr3(alpha,k,beta,l,gamma,m,i)

                        term = 0.0_wp
                        term = term + self_f3 * (dai_k*dai_l*dai_m)
                        term = term + self_f2 * ( d2ai_kl*dai_m + d2ai_km*dai_l + d2ai_lm*dai_k )
                        term = term + self_f1 * d3ai_klm

                        d3Kdr3(alpha,k,beta,l,gamma,m,i,i) = d3Kdr3(alpha,k,beta,l,gamma,m,i,i) + term
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do

end subroutine compute_p16_d3kdr3_full


subroutine compute_p16_d3kdr3_ij(self, nat, xyz, brad, brdr, brdr2, brdr3, i, j, d3Kdr3_ij)
   class(p16_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                         ! (3,nat)
   real(wp), intent(in) :: brad(:)                           ! (nat)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)         ! (3,nat,nat)
   real(wp), contiguous, intent(in) :: brdr2(:, :, :, :, :)  ! (3,nat,3,nat,nat)
   real(wp), contiguous, intent(in) :: brdr3(:, :, :, :, :, :, :)  ! (3,nat,3,nat,3,nat,nat)
   integer, intent(in) :: i, j
   real(wp), contiguous, intent(out) :: d3Kdr3_ij(:, :, :, :, :, :) ! (3,nat,3,nat,3,nat)

   integer :: ip, jp
   integer :: k, l, m
   integer :: alpha, beta, gamma
   integer :: delk, dell, delm

   real(wp), parameter :: tiny_r = 1.0e-14_wp
   real(wp) :: c

   real(wp) :: v(3), r2, r, invr, invr2
   real(wp) :: ai, aj, u, t
   real(wp) :: a1, a16, a17
   real(wp) :: g, invg, invg2, invg3, invg4

   ! g-derivatives wrt r (radii fixed)
   real(wp) :: gr, grr, grrr

   ! g-derivatives wrt u (r fixed)
   real(wp) :: gu, guu, guuu

   ! mixed g_r,u etc (only need gr_u, gr_uu, grr_u)
   real(wp) :: gr_u, gr_uu, grr_u

   ! u-derivatives wrt ai,aj
   real(wp) :: u_ai, u_aj
   real(wp) :: u_aiai, u_ajaj, u_aiaj
   real(wp) :: u_aiaiai, u_ajajaj
   real(wp) :: u_aiaiaj, u_aiajaj

   ! convert g-derivatives to ai/aj derivatives
   real(wp) :: g_ai, g_aj
   real(wp) :: g_aiai, g_ajaj, g_aiaj
   real(wp) :: g_aiaiai, g_ajajaj, g_aiaiaj, g_aiajaj

   real(wp) :: gr_ai, gr_aj
   real(wp) :: gr_aiai, gr_ajaj, gr_aiaj
   real(wp) :: grr_ai, grr_aj

   ! f(g)=keps/g derivatives
   real(wp) :: f1, f2, f3

   ! scalar K partials (r, ai, aj)
   real(wp) :: Kr, Krr, Krrr
   real(wp) :: Kai, Kaj
   real(wp) :: Kaiai, Kajaj, Kaiaj
   real(wp) :: Kr_ai, Kr_aj
   real(wp) :: Krr_ai, Krr_aj
   real(wp) :: Kr_aiai, Kr_ajaj, Kr_aiaj
   real(wp) :: Kaiaiai, Kajajaj, Kaiaiaj, Kaiajaj

   ! radial tensors (v-derivatives holding radii fixed)
   real(wp) :: e(3)
   real(wp) :: I3(3,3)
   real(wp) :: A, B, Bp, Ap
   real(wp) :: Kvv_ai(3,3), Kvv_aj(3,3)
   real(wp) :: Kvvv(3,3,3)
   real(wp) :: Kv_ai(3), Kv_aj(3)
   real(wp) :: Kv_aiai(3), Kv_ajaj(3), Kv_aiaj(3)

   ! coordinate-to-radius scalars
   real(wp) :: dai_k, dai_l, dai_m
   real(wp) :: daj_k, daj_l, daj_m
   real(wp) :: d2ai_kl, d2ai_km, d2ai_lm
   real(wp) :: d2aj_kl, d2aj_km, d2aj_lm
   real(wp) :: d3ai_klm, d3aj_klm

   real(wp) :: term
   real(wp) :: self_f1, self_f2, self_f3

   I3 = 0.0_wp
   I3(1,1)=1.0_wp; I3(2,2)=1.0_wp; I3(3,3)=1.0_wp

   d3Kdr3_ij(:, :, :, :, :, :) = 0.0_wp

   c = zetaP16o16   ! = zeta/16

   ! ==========================================================
   ! Diagonal self term: K_ii = keps / a_i  (exactly as *_full)
   ! ==========================================================
   if (i == j) then
      ai = brad(i)

      self_f1 = -self%keps / (ai*ai)
      self_f2 =  2.0_wp * self%keps / (ai*ai*ai)
      self_f3 = -6.0_wp * self%keps / (ai**4)

      do k = 1, nat
         do l = 1, nat
            do m = 1, nat
               do alpha = 1,3
                  dai_k = brdr(alpha,k,i)
                  do beta = 1,3
                     dai_l   = brdr(beta,l,i)
                     d2ai_kl = brdr2(alpha,k,beta,l,i)
                     do gamma = 1,3
                        dai_m    = brdr(gamma,m,i)
                        d2ai_km  = brdr2(alpha,k,gamma,m,i)
                        d2ai_lm  = brdr2(beta,l,gamma,m,i)
                        d3ai_klm = brdr3(alpha,k,beta,l,gamma,m,i)

                        term = 0.0_wp
                        term = term + self_f3 * (dai_k*dai_l*dai_m)
                        term = term + self_f2 * ( d2ai_kl*dai_m + d2ai_km*dai_l + d2ai_lm*dai_k )
                        term = term + self_f1 * d3ai_klm

                        d3Kdr3_ij(alpha,k,beta,l,gamma,m) = d3Kdr3_ij(alpha,k,beta,l,gamma,m) + term
                     end do
                  end do
               end do
            end do
         end do
      end do
      return
   end if

   ! ==========================================================
   ! Off-diagonal: match *_full semantics (computed for ip>jp and mirrored)
   ! ==========================================================
   if (i > j) then
      ip = i
      jp = j
   else
      ip = j
      jp = i
   end if

   ai = brad(ip)
   aj = brad(jp)

   v(:) = xyz(:, ip) - xyz(:, jp)
   r2   = dot_product(v, v)
   r    = sqrt(r2)
   if (r <= tiny_r) return

   invr  = 1.0_wp / r
   invr2 = invr * invr
   e(:)  = v(:) * invr

   u = sqrt(ai*aj)
   t = u + c*r

   a1  = u / t
   a16 = a1*a1
   a16 = a16*a16
   a16 = a16*a16
   a16 = a16*a16
   a17 = a1 * a16

   g     = r + u * a16
   invg  = 1.0_wp / g
   invg2 = invg*invg
   invg3 = invg2*invg
   invg4 = invg2*invg2

   ! --- g derivatives wrt r (radii fixed) ---
   gr   = 1.0_wp - zetaP16 * a17
   grr  = 17.0_wp * zetaP16 * c * a17 / t
   grrr = -4896.0_wp * c*c*c * a17 / (t*t)

   ! --- g derivatives wrt u (r fixed) ---
   gu   = a16 * (u + 17.0_wp*c*r) / t
   guu  = 272.0_wp * c*c * r2 * a16 / (u * t*t)
   guuu = 816.0_wp * c*c * r2 * (u**14) * (5.0_wp*c*r - u) / (t**19)

   ! --- mixed derivatives needed for Kr_ai etc ---
   gr_u  = -17.0_wp * zetaP16 * c * r * a16 / (t*t)
   gr_uu = -zetaP16 * (34.0_wp * c * r * (u**15) * (8.0_wp*c*r - u) / (t**19))
   grr_u = 17.0_wp * zetaP16 * c * (u**16) * (17.0_wp*c*r - u) / (t**19)

   ! --- u derivatives wrt ai,aj ---
   u_ai   = u / (2.0_wp*ai)
   u_aj   = u / (2.0_wp*aj)
   u_aiai = -u / (4.0_wp*ai*ai)
   u_ajaj = -u / (4.0_wp*aj*aj)
   u_aiaj =  u / (4.0_wp*ai*aj)

   u_aiaiai =  3.0_wp*u / (8.0_wp*ai**3)
   u_ajajaj =  3.0_wp*u / (8.0_wp*aj**3)
   u_aiaiaj = -u / (8.0_wp*ai*ai*aj)
   u_aiajaj = -u / (8.0_wp*ai*aj*aj)

   ! --- g radii derivatives (r fixed) ---
   g_ai   = gu * u_ai
   g_aj   = gu * u_aj
   g_aiai = guu*u_ai*u_ai + gu*u_aiai
   g_ajaj = guu*u_aj*u_aj + gu*u_ajaj
   g_aiaj = guu*u_ai*u_aj + gu*u_aiaj

   g_aiaiai = guuu*u_ai**3 + 3.0_wp*guu*u_ai*u_aiai + gu*u_aiaiai
   g_ajajaj = guuu*u_aj**3 + 3.0_wp*guu*u_aj*u_ajaj + gu*u_ajajaj
   g_aiaiaj = guuu*u_ai*u_ai*u_aj + guu*u_aiai*u_aj + 2.0_wp*guu*u_ai*u_aiaj + gu*u_aiaiaj
   g_aiajaj = guuu*u_aj*u_aj*u_ai + guu*u_ajaj*u_ai + 2.0_wp*guu*u_aj*u_aiaj + gu*u_aiajaj

   ! --- gr radii derivatives (r fixed) ---
   gr_ai   = gr_u * u_ai
   gr_aj   = gr_u * u_aj
   gr_aiai = gr_uu*u_ai*u_ai + gr_u*u_aiai
   gr_ajaj = gr_uu*u_aj*u_aj + gr_u*u_ajaj
   gr_aiaj = gr_uu*u_ai*u_aj + gr_u*u_aiaj

   grr_ai = grr_u * u_ai
   grr_aj = grr_u * u_aj

   ! --- f(g)=keps/g derivatives ---
   f1 = -self%keps * invg2
   f2 =  2.0_wp * self%keps * invg3
   f3 = -6.0_wp * self%keps * invg4

   ! --- scalar K derivatives ---
   Kr   = f1 * gr
   Krr  = f2 * gr*gr + f1 * grr
   Krrr = f3 * gr*gr*gr + 3.0_wp*f2*gr*grr + f1*grrr

   Kai = f1 * g_ai
   Kaj = f1 * g_aj

   Kaiai = f2*g_ai*g_ai + f1*g_aiai
   Kajaj = f2*g_aj*g_aj + f1*g_ajaj
   Kaiaj = f2*g_ai*g_aj + f1*g_aiaj

   Kr_ai = f2*g_ai*gr + f1*gr_ai
   Kr_aj = f2*g_aj*gr + f1*gr_aj

   Krr_ai = f3*g_ai*gr*gr + 2.0_wp*f2*gr*gr_ai + f2*g_ai*grr + f1*grr_ai
   Krr_aj = f3*g_aj*gr*gr + 2.0_wp*f2*gr*gr_aj + f2*g_aj*grr + f1*grr_aj

   Kr_aiai = f3*gr*(g_ai*g_ai) + f2*gr*g_aiai + 2.0_wp*f2*g_ai*gr_ai + f1*gr_aiai
   Kr_ajaj = f3*gr*(g_aj*g_aj) + f2*gr*g_ajaj + 2.0_wp*f2*g_aj*gr_aj + f1*gr_ajaj
   Kr_aiaj = f3*gr*(g_ai*g_aj) + f2*gr*g_aiaj + f2*(g_ai*gr_aj + g_aj*gr_ai) + f1*gr_aiaj

   Kaiaiai = f3*g_ai**3 + 3.0_wp*f2*g_ai*g_aiai + f1*g_aiaiai
   Kajajaj = f3*g_aj**3 + 3.0_wp*f2*g_aj*g_ajaj + f1*g_ajajaj
   Kaiaiaj = f3*(g_ai*g_ai*g_aj) + f2*(g_aiai*g_aj + 2.0_wp*g_ai*g_aiaj) + f1*g_aiaiaj
   Kaiajaj = f3*(g_aj*g_aj*g_ai) + f2*(g_ajaj*g_ai + 2.0_wp*g_aj*g_aiaj) + f1*g_aiajaj

   ! --- build radial tensors in v-space (radii fixed) ---
   A  = Krr - Kr*invr
   B  = Kr*invr
   Bp = (Krr*r - Kr) * invr2
   Ap = Krrr - Bp

   ! Kvv_ai / Kvv_aj
   do alpha = 1,3
      do beta = 1,3
         Kvv_ai(alpha,beta) = (Krr_ai - Kr_ai*invr) * e(alpha)*e(beta) + (Kr_ai*invr) * I3(alpha,beta)
         Kvv_aj(alpha,beta) = (Krr_aj - Kr_aj*invr) * e(alpha)*e(beta) + (Kr_aj*invr) * I3(alpha,beta)
      end do
   end do

   ! Kv_ai etc
   Kv_ai(:)   = Kr_ai   * e(:)
   Kv_aj(:)   = Kr_aj   * e(:)
   Kv_aiai(:) = Kr_aiai * e(:)
   Kv_ajaj(:) = Kr_ajaj * e(:)
   Kv_aiaj(:) = Kr_aiaj * e(:)

   ! Kvvv
   do alpha = 1,3
      do beta = 1,3
         do gamma = 1,3
            Kvvv(alpha,beta,gamma) = Ap * e(alpha)*e(beta)*e(gamma) &
               + (A*invr) * ( I3(alpha,gamma)*e(beta) + I3(beta,gamma)*e(alpha) - 2.0_wp*e(alpha)*e(beta)*e(gamma) ) &
               + Bp * I3(alpha,beta) * e(gamma)
         end do
      end do
   end do

   ! ==========================================================
   ! Assemble coordinate third derivative slab for (ip,jp)
   ! ==========================================================
   do k = 1, nat
      delk = 0; if (k==ip) delk=delk+1; if (k==jp) delk=delk-1
      do l = 1, nat
         dell = 0; if (l==ip) dell=dell+1; if (l==jp) dell=dell-1
         do m = 1, nat
            delm = 0; if (m==ip) delm=delm+1; if (m==jp) delm=delm-1

            do alpha = 1,3
               dai_k = brdr(alpha,k,ip); daj_k = brdr(alpha,k,jp)
               do beta = 1,3
                  dai_l = brdr(beta,l,ip); daj_l = brdr(beta,l,jp)
                  do gamma = 1,3
                     dai_m = brdr(gamma,m,ip); daj_m = brdr(gamma,m,jp)

                     d2ai_kl = brdr2(alpha,k,beta,l,ip)
                     d2ai_km = brdr2(alpha,k,gamma,m,ip)
                     d2ai_lm = brdr2(beta,l,gamma,m,ip)

                     d2aj_kl = brdr2(alpha,k,beta,l,jp)
                     d2aj_km = brdr2(alpha,k,gamma,m,jp)
                     d2aj_lm = brdr2(beta,l,gamma,m,jp)

                     d3ai_klm = brdr3(alpha,k,beta,l,gamma,m,ip)
                     d3aj_klm = brdr3(alpha,k,beta,l,gamma,m,jp)

                     term = 0.0_wp

                     ! (1) vvv term
                     if (delk/=0 .and. dell/=0 .and. delm/=0) then
                        term = term + real(delk*dell*delm,wp) * Kvvv(alpha,beta,gamma)
                     end if

                     ! (2) vv–a terms
                     if (delk/=0 .and. dell/=0) then
                        term = term + real(delk*dell,wp) * ( Kvv_ai(alpha,beta)*dai_m + Kvv_aj(alpha,beta)*daj_m )
                     end if
                     if (delk/=0 .and. delm/=0) then
                        term = term + real(delk*delm,wp) * ( Kvv_ai(alpha,gamma)*dai_l + Kvv_aj(alpha,gamma)*daj_l )
                     end if
                     if (dell/=0 .and. delm/=0) then
                        term = term + real(dell*delm,wp) * ( Kvv_ai(beta,gamma)*dai_k + Kvv_aj(beta,gamma)*daj_k )
                     end if

                     ! (3) v–aa / v–bb / v–ab terms
                     if (delk/=0) then
                        term = term + real(delk,wp) * ( Kv_aiai(alpha)*(dai_l*dai_m) + Kv_ajaj(alpha)*(daj_l*daj_m) &
                           + Kv_aiaj(alpha)*(dai_l*daj_m + daj_l*dai_m) )
                     end if
                     if (dell/=0) then
                        term = term + real(dell,wp) * ( Kv_aiai(beta)*(dai_k*dai_m) + Kv_ajaj(beta)*(daj_k*daj_m) &
                           + Kv_aiaj(beta)*(dai_k*daj_m + daj_k*dai_m) )
                     end if
                     if (delm/=0) then
                        term = term + real(delm,wp) * ( Kv_aiai(gamma)*(dai_k*dai_l) + Kv_ajaj(gamma)*(daj_k*daj_l) &
                           + Kv_aiaj(gamma)*(dai_k*daj_l + daj_k*dai_l) )
                     end if

                     ! (4) radii-only cubic terms
                     term = term + Kaiaiai*(dai_k*dai_l*dai_m) + Kajajaj*(daj_k*daj_l*daj_m)
                     term = term + Kaiaiaj*(dai_k*dai_l*daj_m + dai_k*daj_l*dai_m + daj_k*dai_l*dai_m)
                     term = term + Kaiajaj*(daj_k*daj_l*dai_m + daj_k*dai_l*daj_m + dai_k*daj_l*daj_m)

                     ! (5) K_pq terms with brdr2 (three pairings)
                     if (delm/=0) then
                        term = term + real(delm,wp) * ( d2ai_kl * Kv_ai(gamma) + d2aj_kl * Kv_aj(gamma) )
                     end if
                     term = term + d2ai_kl * (Kaiai*dai_m + Kaiaj*daj_m) + d2aj_kl * (Kaiaj*dai_m + Kajaj*daj_m)

                     if (dell/=0) then
                        term = term + real(dell,wp) * ( d2ai_km * Kv_ai(beta) + d2aj_km * Kv_aj(beta) )
                     end if
                     term = term + d2ai_km * (Kaiai*dai_l + Kaiaj*daj_l) + d2aj_km * (Kaiaj*dai_l + Kajaj*daj_l)

                     if (delk/=0) then
                        term = term + real(delk,wp) * ( d2ai_lm * Kv_ai(alpha) + d2aj_lm * Kv_aj(alpha) )
                     end if
                     term = term + d2ai_lm * (Kaiai*dai_k + Kaiaj*daj_k) + d2aj_lm * (Kaiaj*dai_k + Kajaj*daj_k)

                     ! (6) K_p * brdr3 term
                     term = term + Kai*d3ai_klm + Kaj*d3aj_klm

                     d3Kdr3_ij(alpha,k,beta,l,gamma,m) = d3Kdr3_ij(alpha,k,beta,l,gamma,m) + term
                  end do
               end do
            end do
         end do
      end do
   end do

end subroutine compute_p16_d3kdr3_ij



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


subroutine compute_still_dkdr_ij(self, nat, xyz, brad, brdr, i, j, dKdr_ij)
   class(still_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                 ! (3,nat)
   real(wp), intent(in) :: brad(:)                   ! (nat)
   real(wp), contiguous, intent(in) :: brdr(:, :, :) ! (3,nat,nat)
   integer, intent(in) :: i, j
   real(wp), contiguous, intent(out) :: dKdr_ij(:, :) ! (3,nat) = dK_ij / dr_k

   integer :: k
   real(wp), parameter :: a4 = 0.25_wp
   real(wp) :: rvec(3), r2
   real(wp) :: A, d, E, f2, invf, invf3
   real(wp) :: pref_pos
   real(wp) :: dK_dai, dK_daj

   dKdr_ij(:, :) = 0.0_wp

   if (i == j) then
      ! Diagonal: K_ii = keps/a_i  =>  dK_ii/dr_k = (-keps/a_i^2) * da_i/dr_k
      do k = 1, nat
         dKdr_ij(:, k) = dKdr_ij(:, k) + (-self%keps / (brad(i)*brad(i))) * brdr(:, k, i)
      end do
      return
   end if

   ! Off-diagonal element (i,j)
   rvec(:) = xyz(:, i) - xyz(:, j)
   r2      = dot_product(rvec, rvec)

   A = brad(i) * brad(j)
   d = a4 * r2 / A
   E = exp(-d)

   f2    = r2 + A * E
   invf  = 1.0_wp / sqrt(f2)
   invf3 = invf * invf * invf

   ! Explicit coordinate dependence (radii held fixed)
   pref_pos = -(1.0_wp - a4*E) * (self%keps * invf3)

   dKdr_ij(:, i) = dKdr_ij(:, i) + pref_pos * rvec(:)  ! k=i
   dKdr_ij(:, j) = dKdr_ij(:, j) - pref_pos * rvec(:)  ! k=j

   ! Chain rule via Born radii
   dK_dai = -0.5_wp * self%keps * E * (1.0_wp + d) * brad(j) * invf3
   dK_daj = -0.5_wp * self%keps * E * (1.0_wp + d) * brad(i) * invf3

   do k = 1, nat
      dKdr_ij(:, k) = dKdr_ij(:, k) &
         + dK_dai * brdr(:, k, i) &
         + dK_daj * brdr(:, k, j)
   end do
end subroutine compute_still_dkdr_ij


subroutine compute_still_d2kdr2_ij(self, nat, xyz, brad, brdr, brdr2, i, j, d2Kdr2_ij)
   class(still_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                        ! (3,nat)
   real(wp), intent(in) :: brad(:)                          ! (nat)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)        ! (3,nat,nat)
   real(wp), contiguous, intent(in) :: brdr2(:, :, :, :, :) ! (3,nat,3,nat,nat)
   integer, intent(in) :: i, j
   real(wp), contiguous, intent(out) :: d2Kdr2_ij(:, :, :, :) ! (3,nat,3,nat)

   integer :: k, l
   integer :: delk, dell
   real(wp), parameter :: a4 = 0.25_wp
   real(wp) :: v(3), r2
   real(wp) :: ai, aj, A, d, E, P, S
   real(wp) :: invf, invf3, invf5
   real(wp) :: C, dCdr2, dC_dai, dC_daj
   real(wp) :: Qi, Qj
   real(wp) :: dK_dai, dK_daj
   real(wp) :: d2K_dai2, d2K_daj2, d2K_daida_j
   real(wp) :: I3(3,3), vv(3,3), Hvv(3,3), M(3,3)
   real(wp) :: dk_i(3), dk_j(3), dl_i(3), dl_j(3)
   real(wp) :: coef1, coef2


   ! Identity matrix
   I3 = 0.0_wp
   I3(1,1) = 1.0_wp; I3(2,2) = 1.0_wp; I3(3,3) = 1.0_wp

   d2Kdr2_ij(:, :, :, :) = 0.0_wp

   ! -------------------------
   ! Diagonal self terms: K_ii = keps / a_i
   ! -------------------------
   if (i == j) then
      ai = brad(i)

      coef1 = -self%keps / (ai*ai)
      coef2 =  2.0_wp * self%keps / (ai*ai*ai)

      do k = 1, nat
         dk_i = brdr(:, k, i)
         do l = 1, nat
            dl_i = brdr(:, l, i)

            M = coef2 * (spread(dk_i,2,3)*spread(dl_i,1,3)) + coef1 * brdr2(:, k, :, l, i)
            d2Kdr2_ij(:, k, :, l) = d2Kdr2_ij(:, k, :, l) + M
         end do
      end do
      return
   end if

   ! -------------------------
   ! Off-diagonal element (i,j)
   ! -------------------------
   ai = brad(i)
   aj = brad(j)

   v(:) = xyz(:, i) - xyz(:, j)
   r2   = dot_product(v, v)

   A = ai * aj
   d = a4 * r2 / A
   E = exp(-d)
   P = 1.0_wp - a4 * E

   S     = r2 + A * E
   invf  = 1.0_wp / sqrt(S)
   invf3 = invf*invf*invf
   invf5 = invf3*invf*invf

   Qi = aj * E * (1.0_wp + d)
   Qj = ai * E * (1.0_wp + d)

   dK_dai = -0.5_wp * self%keps * Qi * invf3
   dK_daj = -0.5_wp * self%keps * Qj * invf3

   C = - self%keps * P * invf3

   dCdr2 = -self%keps * ( (a4*a4 * E / A) * invf3 - 1.5_wp * (P*P) * invf5 )

   vv  = spread(v, 2, 3) * spread(v, 1, 3)
   Hvv = C * I3 + (2.0_wp * dCdr2) * vv

   dC_dai = self%keps * (a4 * E * d / ai) * invf3 + 1.5_wp * self%keps * P * Qi * invf5
   dC_daj = self%keps * (a4 * E * d / aj) * invf3 + 1.5_wp * self%keps * P * Qj * invf5

   d2K_dai2 = -0.5_wp * self%keps * (aj * E * d*d / ai) * invf3 + 0.75_wp * self%keps * (Qi*Qi) * invf5
   d2K_daj2 = -0.5_wp * self%keps * (ai * E * d*d / aj) * invf3 + 0.75_wp * self%keps * (Qj*Qj) * invf5
   d2K_daida_j = -0.5_wp * self%keps * (E * (1.0_wp + d + d*d)) * invf3 + 0.75_wp * self%keps * (Qi*Qj) * invf5

   do k = 1, nat
      dk_i = brdr(:, k, i)
      dk_j = brdr(:, k, j)

      delk = 0
      if (k == i) delk = delk + 1
      if (k == j) delk = delk - 1

      do l = 1, nat
         dl_i = brdr(:, l, i)
         dl_j = brdr(:, l, j)

         dell = 0
         if (l == i) dell = dell + 1
         if (l == j) dell = dell - 1

         M = 0.0_wp

         if (delk /= 0 .and. dell /= 0) then
            M = M + real(delk*dell, wp) * Hvv
         end if

         if (delk /= 0) then
            M = M + real(delk, wp) * ( dC_dai * (spread(v,2,3)*spread(dl_i,1,3)) &
                                     + dC_daj * (spread(v,2,3)*spread(dl_j,1,3)) )
         end if

         if (dell /= 0) then
            M = M + real(dell, wp) * ( dC_dai * (spread(dk_i,2,3)*spread(v,1,3)) &
                                     + dC_daj * (spread(dk_j,2,3)*spread(v,1,3)) )
         end if

         M = M + d2K_dai2     * (spread(dk_i,2,3)*spread(dl_i,1,3))
         M = M + d2K_daj2     * (spread(dk_j,2,3)*spread(dl_j,1,3))
         M = M + d2K_daida_j  * ( (spread(dk_i,2,3)*spread(dl_j,1,3)) &
                                + (spread(dk_j,2,3)*spread(dl_i,1,3)) )

         M = M + dK_dai * brdr2(:, k, :, l, i) + dK_daj * brdr2(:, k, :, l, j)

         d2Kdr2_ij(:, k, :, l) = d2Kdr2_ij(:, k, :, l) + M
      end do
   end do

end subroutine compute_still_d2kdr2_ij


subroutine compute_still_d3kdr3_ij(self, nat, xyz, brad, brdr, brdr2, brdr3, i, j, d3Kdr3_ij)
   class(still_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                              ! (3,nat)
   real(wp), intent(in) :: brad(:)                                ! (nat)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)              ! (3,nat,nat)
   real(wp), contiguous, intent(in) :: brdr2(:, :, :, :, :)       ! (3,nat,3,nat,nat)
   real(wp), contiguous, intent(in) :: brdr3(:, :, :, :, :, :, :) ! (3,nat,3,nat,3,nat,nat)
   integer, intent(in) :: i, j
   real(wp), contiguous, intent(out) :: d3Kdr3_ij(:, :, :, :, :, :) ! (3,nat,3,nat,3,nat)

   integer :: ip, jp
   integer :: k, l, m
   integer :: alpha, beta, gamma
   integer :: delk, dell, delm
   real(wp), parameter :: a4 = 0.25_wp

   real(wp) :: v(3), r2
   real(wp) :: ai, aj, A, d, E, P, S
   real(wp) :: invf, invf3, invf5, invf7

   ! --- S-derivatives wrt (r2, ai, aj) ---
   real(wp) :: Sr, Sa, Sb
   real(wp) :: Srr, Srrr
   real(wp) :: Sra, Srb
   real(wp) :: Srra, Srrb
   real(wp) :: Saa2, Sbb2, Sab
   real(wp) :: Sraa, Srbb, Srab
   real(wp) :: Saa3, Sbb3, Saa2b, Sa2bb

   ! --- F(S)=keps*S^{-1/2} derivatives ---
   real(wp) :: F1, F2, F3

   ! --- K scalar partials needed ---
   real(wp) :: Kr, Krr, Krrr
   real(wp) :: Kai, Kaj
   real(wp) :: Kaiai, Kajaj, Kaiaj
   real(wp) :: Kr_ai, Kr_aj
   real(wp) :: Krr_ai, Krr_aj
   real(wp) :: Kr_aiai, Kr_ajaj, Kr_aiaj
   real(wp) :: Kaiaiai, Kajajaj, Kaiaiaj, Kaiajaj

   ! --- tensor building blocks ---
   real(wp) :: Kvvv(3,3,3)
   real(wp) :: Kvv_ai(3,3), Kvv_aj(3,3)
   real(wp) :: Kv_ai(3), Kv_aj(3)
   real(wp) :: Kv_aiai(3), Kv_ajaj(3), Kv_aiaj(3)
   real(wp) :: Kv(3)
   real(wp) :: I3(3,3)

   ! coordinate-to-radius derivatives (scalars)
   real(wp) :: dai_k, dai_l, dai_m
   real(wp) :: daj_k, daj_l, daj_m

   real(wp) :: d2ai_kl, d2ai_km, d2ai_lm
   real(wp) :: d2aj_kl, d2aj_km, d2aj_lm

   real(wp) :: d3ai_klm
   real(wp) :: d3aj_klm

   real(wp) :: term
   real(wp) :: f11, f22, f3_self

   I3 = 0.0_wp
   I3(1,1)=1.0_wp; I3(2,2)=1.0_wp; I3(3,3)=1.0_wp

   d3Kdr3_ij(:, :, :, :, :, :) = 0.0_wp

   ! -------------------------
   ! Diagonal self term: K_ii = keps / a_i   (identical to *_full readout)
   ! -------------------------
   if (i == j) then
      ai = brad(i)

      f11     = -self%keps / (ai*ai)
      f22     =  2.0_wp * self%keps / (ai*ai*ai)
      f3_self = -6.0_wp * self%keps / (ai**4)

      do k = 1, nat
         do l = 1, nat
            do m = 1, nat
               do alpha = 1,3
                  dai_k = brdr(alpha,k,i)
                  do beta = 1,3
                     dai_l   = brdr(beta,l,i)
                     d2ai_kl = brdr2(alpha,k,beta,l,i)
                     do gamma = 1,3
                        dai_m    = brdr(gamma,m,i)
                        d2ai_km  = brdr2(alpha,k,gamma,m,i)
                        d2ai_lm  = brdr2(beta,l,gamma,m,i)
                        d3ai_klm = brdr3(alpha,k,beta,l,gamma,m,i)

                        term = 0.0_wp
                        term = term + f3_self * (dai_k*dai_l*dai_m)
                        term = term + f22 * ( d2ai_kl*dai_m + d2ai_km*dai_l + d2ai_lm*dai_k )
                        term = term + f11 * d3ai_klm

                        d3Kdr3_ij(alpha,k,beta,l,gamma,m) = d3Kdr3_ij(alpha,k,beta,l,gamma,m) + term
                     end do
                  end do
               end do
            end do
         end do
      end do
      return
   end if

   ! -------------------------
   ! Off-diagonal: MUST match *_full semantics:
   ! *_full computes only for ip>jp and then mirrors into (jp,ip),
   ! so the value stored at (i,j) is always the ip>jp computation.
   ! -------------------------
   if (i > j) then
      ip = i
      jp = j
   else
      ip = j
      jp = i
   end if

   ai = brad(ip)
   aj = brad(jp)

   v(:) = xyz(:, ip) - xyz(:, jp)
   r2   = dot_product(v, v)

   A = ai * aj
   d = a4 * r2 / A
   E = exp(-d)

   S = r2 + A * E
   invf  = 1.0_wp / sqrt(S)
   invf3 = invf*invf*invf
   invf5 = invf3*invf*invf
   invf7 = invf5*invf*invf

   F1 = -0.5_wp   * self%keps * invf3
   F2 =  0.75_wp  * self%keps * invf5
   F3 = -1.875_wp * self%keps * invf7

   P  = 1.0_wp - a4*E
   Sr = P
   Srr  = (a4*a4 / A) * E
   Srrr = -(a4*a4*a4 / (A*A)) * E

   Sa = aj * E * (1.0_wp + d)
   Sb = ai * E * (1.0_wp + d)

   Saa2 = aj * E * (d*d) / ai
   Sbb2 = ai * E * (d*d) / aj
   Sab  = E * (1.0_wp + d + d*d)

   Sra = -(a4/ai) * E * d
   Srb = -(a4/aj) * E * d

   Srra = -(a4*a4/(ai*A)) * E * (1.0_wp - d)
   Srrb = -(a4*a4/(aj*A)) * E * (1.0_wp - d)

   Sraa = -a4 * E * d * (d - 2.0_wp) / (ai*ai)
   Srbb = -a4 * E * d * (d - 2.0_wp) / (aj*aj)
   Srab = -(a4/A) * E * d * (d - 1.0_wp)

   Saa3   = aj * E * d*d * (d - 3.0_wp) / (ai*ai)
   Sbb3   = ai * E * d*d * (d - 3.0_wp) / (aj*aj)
   Saa2b  = E * d*d * (d - 1.0_wp) / ai
   Sa2bb  = E * d*d * (d - 1.0_wp) / aj

   Kr  = F1 * Sr
   Kai = F1 * Sa
   Kaj = F1 * Sb

   Krr    = F2 * Sr*Sr + F1 * Srr
   Kaiai  = F2 * Sa*Sa + F1 * Saa2
   Kajaj  = F2 * Sb*Sb + F1 * Sbb2
   Kaiaj  = F2 * Sa*Sb + F1 * Sab

   Kr_ai  = F2 * Sr*Sa + F1 * Sra
   Kr_aj  = F2 * Sr*Sb + F1 * Srb

   Krrr = F3 * Sr*Sr*Sr + 3.0_wp*F2*Srr*Sr + F1*Srrr

   Krr_ai = F3 * Sr*Sr*Sa + F2*(Srr*Sa + 2.0_wp*Sr*Sra) + F1*Srra
   Krr_aj = F3 * Sr*Sr*Sb + F2*(Srr*Sb + 2.0_wp*Sr*Srb) + F1*Srrb

   Kr_aiai = F3 * Sr*Sa*Sa + F2*(Saa2*Sr + 2.0_wp*Sra*Sa) + F1*Sraa
   Kr_ajaj = F3 * Sr*Sb*Sb + F2*(Sbb2*Sr + 2.0_wp*Srb*Sb) + F1*Srbb
   Kr_aiaj = F3 * Sr*Sa*Sb + F2*(Sab*Sr + Sra*Sb + Srb*Sa) + F1*Srab

   Kaiaiai = F3 * Sa*Sa*Sa + 3.0_wp*F2*Saa2*Sa + F1*Saa3
   Kajajaj = F3 * Sb*Sb*Sb + 3.0_wp*F2*Sbb2*Sb + F1*Sbb3
   Kaiaiaj = F3 * Sa*Sa*Sb + F2*(Saa2*Sb + 2.0_wp*Sab*Sa) + F1*Saa2b
   Kaiajaj = F3 * Sa*Sb*Sb + F2*(Sbb2*Sa + 2.0_wp*Sab*Sb) + F1*Sa2bb

   Kv(:) = 2.0_wp * Kr * v(:)

   Kv_ai(:)   = 2.0_wp * Kr_ai   * v(:)
   Kv_aj(:)   = 2.0_wp * Kr_aj   * v(:)
   Kv_aiai(:) = 2.0_wp * Kr_aiai * v(:)
   Kv_ajaj(:) = 2.0_wp * Kr_ajaj * v(:)
   Kv_aiaj(:) = 2.0_wp * Kr_aiaj * v(:)

   Kvv_ai(:,:) = 2.0_wp*Kr_ai * I3(:,:) + 4.0_wp*Krr_ai * (spread(v,2,3)*spread(v,1,3))
   Kvv_aj(:,:) = 2.0_wp*Kr_aj * I3(:,:) + 4.0_wp*Krr_aj * (spread(v,2,3)*spread(v,1,3))

   Kvvv(:,:,:) = 0.0_wp
   do alpha = 1,3
      do beta = 1,3
         do gamma = 1,3
            Kvvv(alpha,beta,gamma) = 8.0_wp*Krrr * v(alpha)*v(beta)*v(gamma) &
               + 4.0_wp*Krr * ( I3(alpha,beta)*v(gamma) + I3(alpha,gamma)*v(beta) + I3(beta,gamma)*v(alpha) )
         end do
      end do
   end do

   ! Assemble coordinate third derivative slab for (ip,jp) — which equals the
   ! value stored by *_full at BOTH (ip,jp) and (jp,ip).
   do k = 1, nat
      delk = 0; if (k==ip) delk=delk+1; if (k==jp) delk=delk-1
      do l = 1, nat
         dell = 0; if (l==ip) dell=dell+1; if (l==jp) dell=dell-1
         do m = 1, nat
            delm = 0; if (m==ip) delm=delm+1; if (m==jp) delm=delm-1

            do alpha = 1,3
               dai_k = brdr(alpha,k,ip)
               daj_k = brdr(alpha,k,jp)
               do beta = 1,3
                  dai_l = brdr(beta,l,ip)
                  daj_l = brdr(beta,l,jp)
                  do gamma = 1,3
                     dai_m = brdr(gamma,m,ip)
                     daj_m = brdr(gamma,m,jp)

                     d2ai_kl = brdr2(alpha,k,beta,l,ip)
                     d2ai_km = brdr2(alpha,k,gamma,m,ip)
                     d2ai_lm = brdr2(beta,l,gamma,m,ip)

                     d2aj_kl = brdr2(alpha,k,beta,l,jp)
                     d2aj_km = brdr2(alpha,k,gamma,m,jp)
                     d2aj_lm = brdr2(beta,l,gamma,m,jp)

                     d3ai_klm = brdr3(alpha,k,beta,l,gamma,m,ip)
                     d3aj_klm = brdr3(alpha,k,beta,l,gamma,m,jp)

                     term = 0.0_wp

                     if (delk/=0 .and. dell/=0 .and. delm/=0) then
                        term = term + real(delk*dell*delm,wp) * Kvvv(alpha,beta,gamma)
                     end if

                     if (delk/=0 .and. dell/=0) then
                        term = term + real(delk*dell,wp) * Kvv_ai(alpha,beta) * dai_m
                        term = term + real(delk*dell,wp) * Kvv_aj(alpha,beta) * daj_m
                     end if
                     if (delk/=0 .and. delm/=0) then
                        term = term + real(delk*delm,wp) * Kvv_ai(alpha,gamma) * dai_l
                        term = term + real(delk*delm,wp) * Kvv_aj(alpha,gamma) * daj_l
                     end if
                     if (dell/=0 .and. delm/=0) then
                        term = term + real(dell*delm,wp) * Kvv_ai(beta,gamma) * dai_k
                        term = term + real(dell*delm,wp) * Kvv_aj(beta,gamma) * daj_k
                     end if

                     if (delk/=0) then
                        term = term + real(delk,wp) * Kv_aiai(alpha) * (dai_l*dai_m)
                        term = term + real(delk,wp) * Kv_ajaj(alpha) * (daj_l*daj_m)
                        term = term + real(delk,wp) * Kv_aiaj(alpha) * (dai_l*daj_m + daj_l*dai_m)
                     end if
                     if (dell/=0) then
                        term = term + real(dell,wp) * Kv_aiai(beta) * (dai_k*dai_m)
                        term = term + real(dell,wp) * Kv_ajaj(beta) * (daj_k*daj_m)
                        term = term + real(dell,wp) * Kv_aiaj(beta) * (dai_k*daj_m + daj_k*dai_m)
                     end if
                     if (delm/=0) then
                        term = term + real(delm,wp) * Kv_aiai(gamma) * (dai_k*dai_l)
                        term = term + real(delm,wp) * Kv_ajaj(gamma) * (daj_k*daj_l)
                        term = term + real(delm,wp) * Kv_aiaj(gamma) * (dai_k*daj_l + daj_k*dai_l)

                     end if

                     term = term + Kaiaiai * (dai_k*dai_l*dai_m)
                     term = term + Kajajaj * (daj_k*daj_l*daj_m)

                     term = term + Kaiaiaj * (dai_k*dai_l*daj_m + dai_k*daj_l*dai_m + daj_k*dai_l*dai_m)
                     term = term + Kaiajaj * (daj_k*daj_l*dai_m + daj_k*dai_l*daj_m + dai_k*daj_l*daj_m)

                     if (delm/=0) then
                        term = term + real(delm,wp) * d2ai_kl * Kv_ai(gamma)
                        term = term + real(delm,wp) * d2aj_kl * Kv_aj(gamma)
                     end if
                     term = term + d2ai_kl * ( Kaiai * dai_m + Kaiaj * daj_m )
                     term = term + d2aj_kl * ( Kaiaj * dai_m + Kajaj * daj_m )

                     if (dell/=0) then
                        term = term + real(dell,wp) * d2ai_km * Kv_ai(beta)
                        term = term + real(dell,wp) * d2aj_km * Kv_aj(beta)
                     end if
                     term = term + d2ai_km * ( Kaiai * dai_l + Kaiaj * daj_l )
                     term = term + d2aj_km * ( Kaiaj * dai_l + Kajaj * daj_l )

                     if (delk/=0) then
                        term = term + real(delk,wp) * d2ai_lm * Kv_ai(alpha)
                        term = term + real(delk,wp) * d2aj_lm * Kv_aj(alpha)
                     end if
                     term = term + d2ai_lm * ( Kaiai * dai_k + Kaiaj * daj_k )
                     term = term + d2aj_lm * ( Kaiaj * dai_k + Kajaj * daj_k )

                     term = term + Kai * d3ai_klm + Kaj * d3aj_klm

                     d3Kdr3_ij(alpha,k,beta,l,gamma,m) = d3Kdr3_ij(alpha,k,beta,l,gamma,m) + term
                  end do
               end do
            end do
         end do
      end do
   end do

end subroutine compute_still_d3kdr3_ij




end module tblite_solvation_kernel
