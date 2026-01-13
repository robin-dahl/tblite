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
   public :: compute_kernel_dkdr, compute_kernel_d2kdr2, compute_kernel_d3kdr3 ! convenience dispatcher by enum

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
      procedure(compute_kernel_d2kdr2_interface), deferred :: compute_kernel_d2kdr2
      procedure(compute_kernel_d3kdr3_interface), deferred :: compute_kernel_d3kdr3
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

      !> Full element-wise 2nd derivative tensor of the actual kernel matrix:
      subroutine compute_kernel_d2kdr2_interface(self, nat, xyz, brad, brdr, brdr2, d2Kdr2)
         import :: kernel_type, wp
         class(kernel_type), intent(in) :: self
         integer, intent(in) :: nat
         real(wp), intent(in) :: xyz(:, :)
         real(wp), intent(in) :: brad(:)
         real(wp), contiguous, intent(in) :: brdr(:, :, :)
         real(wp), contiguous, intent(in) :: brdr2(:, :, :, :, :)
         real(wp), contiguous, intent(out) :: d2Kdr2(:, :, :, :, :, :)
      end subroutine compute_kernel_d2kdr2_interface

      !> Full element-wise 3rd derivative tensor of the actual kernel matrix:
      subroutine compute_kernel_d3kdr3_interface(self, nat, xyz, brad, brdr, brdr2, brdr3, d3Kdr3)
         import :: kernel_type, wp
         class(kernel_type), intent(in) :: self
         integer, intent(in) :: nat
         real(wp), intent(in) :: xyz(:, :)
         real(wp), intent(in) :: brad(:)
         real(wp), contiguous, intent(in) :: brdr(:, :, :)
         real(wp), contiguous, intent(in) :: brdr2(:, :, :, :, :)
         real(wp), contiguous, intent(in) :: brdr3(:, :, :, :, :, :, :)
         real(wp), contiguous, intent(out) :: d3Kdr3(:, :, :, :, :, :, :, :)
      end subroutine compute_kernel_d3kdr3_interface
   end interface

   type, extends(kernel_type) :: still_kernel
   contains
      procedure :: add_kernel_mat   => add_still_mat
      procedure :: add_kernel_deriv => add_still_deriv
      procedure :: compute_kernel_dkdr => compute_still_dkdr_full
      procedure :: compute_kernel_d2kdr2 => compute_still_d2kdr2_full
      procedure :: compute_kernel_d3kdr3 => compute_still_d3kdr3_full
   end type still_kernel

   type, extends(kernel_type) :: p16_kernel
   contains
      procedure :: add_kernel_mat   => add_p16_mat
      procedure :: add_kernel_deriv => add_p16_deriv
      procedure :: compute_kernel_dkdr => compute_p16_dkdr_full
      procedure :: compute_kernel_d2kdr2 => compute_p16_d2kdr2_full
      procedure :: compute_kernel_d3kdr3 => compute_p16_d3kdr3_full
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

!> Convenience dispatcher (switches by kernel enum, returns full derivative tensor)
subroutine compute_kernel_d2kdr2(kernel_id, keps, nat, xyz, brad, brdr, brdr2, d2Kdr2)
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: brad(:)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)
   real(wp), contiguous, intent(in) :: brdr2(:, :, :, :, :)
   real(wp), contiguous, intent(out) :: d2Kdr2(:, :, :, :, :, :)

   class(kernel_type), allocatable :: kernel

   kernel = new_kernel(kernel_id, keps)
   call kernel%compute_kernel_d2kdr2(nat, xyz, brad, brdr, brdr2, d2Kdr2)
end subroutine compute_kernel_d2kdr2

!> Convenience dispatcher (switches by kernel enum, returns full derivative tensor)
subroutine compute_kernel_d3kdr3(kernel_id, keps, nat, xyz, brad, brdr, brdr2, brdr3, d3Kdr3)
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: brad(:)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)
   real(wp), contiguous, intent(in) :: brdr2(:, :, :, :, :)
   real(wp), contiguous, intent(in) :: brdr3(:, :, :, :, :, :, :)
   real(wp), contiguous, intent(out) :: d3Kdr3(:, :, :, :, :, :, :, :)

   class(kernel_type), allocatable :: kernel

   kernel = new_kernel(kernel_id, keps)
   call kernel%compute_kernel_d3kdr3(nat, xyz, brad, brdr, brdr2, brdr3, d3Kdr3)
end subroutine compute_kernel_d3kdr3

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

!> Full d2Kdr2 for the actual P16 kernel matrix (including diagonal self term)
subroutine compute_p16_d2kdr2_full(self, nat, xyz, brad, brdr, brdr2, d2Kdr2)
   class(p16_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                         ! (3,nat)
   real(wp), intent(in) :: brad(:)                           ! (nat)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)         ! (3,nat,nat)
   real(wp), contiguous, intent(in) :: brdr2(:, :, :, :, :)  ! (3,nat,3,nat,nat)
   real(wp), contiguous, intent(out) :: d2Kdr2(:, :, :, :, :, :) ! (3,nat,3,nat,nat,nat)

   integer :: i, j, k, l
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

   d2Kdr2(:, :, :, :, :, :) = 0.0_wp

   cp16 = zetaP16o16 

   ! -------------------------
   ! Off-diagonal: i > j, then mirror
   ! -------------------------
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

         u = sqrt(ai * aj)
         t = u + cp16 * r

         ! a1 = u/t; a16 = a1^16; a17 = a1^17
         a1  = u / t
         a16 = a1 * a1
         a16 = a16 * a16
         a16 = a16 * a16
         a16 = a16 * a16
         a17 = a1 * a16

         g    = r + u * a16
         invg = 1.0_wp / g
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

         ! Hessian wrt v (radii held fixed), isotropic formula:
         ! Hvv = (Kr/r) I + (Krr - Kr/r) (v⊗v)/r^2
         vv  = spread(v,2,3) * spread(v,1,3)
         Hvv = C * I3 + (Krr - C) * (vv * invr2)

         ! u-derivatives for radii chain terms (r fixed)
         ! g(u) = r + u^17 / t^16 = r + u*a16
         ! g_u  = u^16 (u + 17 c r) / t^17 = a16*(u + 17 c r)/t
         gu  = a16 * (u + 17.0_wp*cp16*r) / t
         ! g_uu = 272 c^2 r^2 u^15 / t^18 = 272 c^2 r^2 * a16 / (u*t^2)
         guu = 272.0_wp * cp16*cp16 * r2 * a16 / (u * t*t)

         ! u derivatives
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
         ! C = (K_r)/r, so dC/dai = (1/r) d(K_r)/dai
         ! K_r = -keps * g_r / g^2
         ! g_r(u) = 1 - zeta * u^17 / t^17  => g_r_u = -17 zeta c r * u^16 / t^18 = -17 zeta c r * a16 / t^2
         gr_u  = -17.0_wp * zetaP16 * cp16 * r * a16 / (t*t)
         gr_ai = gr_u * u_ai
         gr_aj = gr_u * u_aj

         Kr_ai = -self%keps * ( gr_ai * invg2 - 2.0_wp * gr * g_ai * invg3 )
         Kr_aj = -self%keps * ( gr_aj * invg2 - 2.0_wp * gr * g_aj * invg3 )

         dC_dai = Kr_ai * invr
         dC_daj = Kr_aj * invr

         ! Assemble full coordinate Hessian blocks for all (k,l)
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

               ! (2) Mixed v–a parts (changes in C through radii)
               if (delk /= 0) then
                  M = M + real(delk, wp) * ( dC_dai * (spread(v,2,3)*spread(dl_i,1,3)) &
                                           + dC_daj * (spread(v,2,3)*spread(dl_j,1,3)) )
               end if

               if (dell /= 0) then
                  M = M + real(dell, wp) * ( dC_dai * (spread(dk_i,2,3)*spread(v,1,3)) &
                                           + dC_daj * (spread(dk_j,2,3)*spread(v,1,3)) )
               end if

               ! (3) a–a parts (via brdr)
               M = M + d2K_dai2    * (spread(dk_i,2,3)*spread(dl_i,1,3))
               M = M + d2K_daj2    * (spread(dk_j,2,3)*spread(dl_j,1,3))
               M = M + d2K_daida_j * ( (spread(dk_i,2,3)*spread(dl_j,1,3)) &
                                     + (spread(dk_j,2,3)*spread(dl_i,1,3)) )

               ! (4) brdr2 terms
               M = M + dK_dai * brdr2(:, k, :, l, i) + dK_daj * brdr2(:, k, :, l, j)

               d2Kdr2(:, k, :, l, i, j) = d2Kdr2(:, k, :, l, i, j) + M
               d2Kdr2(:, k, :, l, j, i) = d2Kdr2(:, k, :, l, j, i) + M
            end do
         end do

      end do
   end do

   ! -------------------------
   ! Diagonal self terms: K_ii = keps / a_i   (same form as Still)
   ! -------------------------
   do i = 1, nat
      ai = brad(i)

      coef1 = -self%keps / (ai*ai)
      coef2 =  2.0_wp * self%keps / (ai*ai*ai)

      do k = 1, nat
         dk_i = brdr(:, k, i)
         do l = 1, nat
            dl_i = brdr(:, l, i)

            M = coef2 * (spread(dk_i,2,3)*spread(dl_i,1,3)) + coef1 * brdr2(:, k, :, l, i)

            d2Kdr2(:, k, :, l, i, i) = d2Kdr2(:, k, :, l, i, i) + M
         end do
      end do
   end do

end subroutine compute_p16_d2kdr2_full


!> Full d3Kdr3 for the actual P16 kernel matrix (including diagonal self term)
subroutine compute_p16_d3kdr3_full(self, nat, xyz, brad, brdr, brdr2, brdr3, d3Kdr3)
   class(p16_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                         ! (3,nat)
   real(wp), intent(in) :: brad(:)                           ! (nat)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)         ! (3,nat,nat)
   real(wp), contiguous, intent(in) :: brdr2(:, :, :, :, :)  ! (3,nat,3,nat,nat)
   real(wp), contiguous, intent(in) :: brdr3(:, :, :, :, :, :, :)  ! (3, nat, 3,nat,3,nat,nat)
   real(wp), contiguous, intent(out) :: d3Kdr3(:, :, :, :, :, :, :, :) ! (3, nat, 3,nat,3,nat,nat,nat)


   print *, "compute_p16_d3kdr3_full: Not yet implemented!"


end subroutine compute_p16_d3kdr3_full


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


! !> Full d2Kdr2 for the actual Still kernel matrix (including diagonal self term)
subroutine compute_still_d2kdr2_full(self, nat, xyz, brad, brdr, brdr2, d2Kdr2)
   class(still_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                      ! (3,nat)
   real(wp), intent(in) :: brad(:)                        ! (nat)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)      ! (3,nat,nat)
   real(wp), contiguous, intent(in) :: brdr2(:, :, :, :, :) ! (3,nat,3,nat,nat)
   real(wp), contiguous, intent(out) :: d2Kdr2(:, :, :, :, :, :) ! (3,nat,3,nat,nat,nat)

   integer :: i, j, k, l
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

   d2Kdr2(:, :, :, :, :, :) = 0.0_wp

   ! -------------------------
   ! Off-diagonal: i > j, then mirror (j,i)
   ! -------------------------
   do i = 1, nat
      ai = brad(i)
      do j = 1, i-1
         aj = brad(j)

         v(:) = xyz(:, i) - xyz(:, j)
         r2   = dot_product(v, v)

         A = ai * aj
         d = a4 * r2 / A
         E = exp(-d)
         P = 1.0_wp - a4 * E                     ! = 1 - (1/4)exp(-d)

         S    = r2 + A * E                       ! f^2
         invf = 1.0_wp / sqrt(S)
         invf3 = invf*invf*invf                  ! S^(-3/2)
         invf5 = invf3*invf*invf                 ! S^(-5/2)

         ! Useful “Born” blocks (also equal to ∂S/∂a_i and ∂S/∂a_j)
         Qi = aj * E * (1.0_wp + d)              ! = ∂S/∂a_i
         Qj = ai * E * (1.0_wp + d)              ! = ∂S/∂a_j

         ! First partials wrt radii (off-diagonal kernel)
         dK_dai = -0.5_wp * self%keps * Qi * invf3
         dK_daj = -0.5_wp * self%keps * Qj * invf3

         ! Gradient coefficient wrt v: ∂K/∂v = C * v
         C = - self%keps * P * invf3

         ! dC/dr2 (holding radii fixed)
         ! dC/dr2 = -keps * [ (a4^2 E / A) S^(-3/2) - (3/2) P^2 S^(-5/2) ]
         dCdr2 = -self%keps * ( (a4*a4 * E / A) * invf3 - 1.5_wp * (P*P) * invf5 )

         ! Hessian wrt v (holding radii fixed): Hvv = C*I + 2*dC/dr2 * (v⊗v)
         vv = spread(v, 2, 3) * spread(v, 1, 3)
         Hvv = C * I3 + (2.0_wp * dCdr2) * vv

         ! Mixed partials: ∂²K/(∂v ∂a_i) = v * (∂C/∂a_i)
         ! ∂C/∂a_i = keps*(a4*E*d/ai)*S^(-3/2) + (3/2)*keps*P*(∂S/∂a_i)*S^(-5/2)
         dC_dai = self%keps * (a4 * E * d / ai) * invf3 + 1.5_wp * self%keps * P * Qi * invf5
         dC_daj = self%keps * (a4 * E * d / aj) * invf3 + 1.5_wp * self%keps * P * Qj * invf5

         ! Second partials wrt radii
         ! d²K/da_i²
         d2K_dai2 = -0.5_wp * self%keps * (aj * E * d*d / ai) * invf3 + 0.75_wp * self%keps * (Qi*Qi) * invf5
         ! d²K/da_j²
         d2K_daj2 = -0.5_wp * self%keps * (ai * E * d*d / aj) * invf3 + 0.75_wp * self%keps * (Qj*Qj) * invf5
         ! d²K/(da_i da_j)
         d2K_daida_j = -0.5_wp * self%keps * (E * (1.0_wp + d + d*d)) * invf3 + 0.75_wp * self%keps * (Qi*Qj) * invf5

         ! Now assemble full coordinate Hessian blocks for all (k,l)
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

               ! (1) Explicit vv part (only if k,l hit i/j via deltas)
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

               ! (3) a–a parts (via brdr)
               M = M + d2K_dai2     * (spread(dk_i,2,3)*spread(dl_i,1,3))
               M = M + d2K_daj2     * (spread(dk_j,2,3)*spread(dl_j,1,3))
               M = M + d2K_daida_j  * ( (spread(dk_i,2,3)*spread(dl_j,1,3)) &
                                      + (spread(dk_j,2,3)*spread(dl_i,1,3)) )

               ! (4) second-radius-derivative parts (via brdr2)
               M = M + dK_dai * brdr2(:, k, :, l, i) + dK_daj * brdr2(:, k, :, l, j)

               d2Kdr2(:, k, :, l, i, j) = d2Kdr2(:, k, :, l, i, j) + M
               d2Kdr2(:, k, :, l, j, i) = d2Kdr2(:, k, :, l, j, i) + M  ! mirror symmetry K_ij = K_ji
            end do
         end do
      end do
   end do

   ! -------------------------
   ! Diagonal self terms: K_ii = keps / a_i
   ! -------------------------
   do i = 1, nat
      ai = brad(i)

      ! K = keps * a_i^{-1}
      ! ∂K/∂r_k = (-keps/a_i^2) * ∂a_i/∂r_k
      ! ∂²K/∂r_k∂r_l =
      !    (-keps/a_i^2) * ∂²a_i/∂r_k∂r_l
      !  + (2 keps/a_i^3) * (∂a_i/∂r_k) ⊗ (∂a_i/∂r_l)
      coef1 = -self%keps / (ai*ai)
      coef2 =  2.0_wp * self%keps / (ai*ai*ai)

      do k = 1, nat
         dk_i = brdr(:, k, i)
         do l = 1, nat
            dl_i = brdr(:, l, i)

            M = coef2 * (spread(dk_i,2,3)*spread(dl_i,1,3)) + coef1 * brdr2(:, k, :, l, i)

            d2Kdr2(:, k, :, l, i, i) = d2Kdr2(:, k, :, l, i, i) + M
         end do
      end do
   end do

end subroutine compute_still_d2kdr2_full


 !> Full d3Kdr3 for the actual Still kernel matrix (including diagonal self term)
subroutine compute_still_d3kdr3_full(self, nat, xyz, brad, brdr, brdr2, brdr3, d3Kdr3)
   class(still_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                      ! (3,nat)
   real(wp), intent(in) :: brad(:)                        ! (nat)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)      ! (3,nat,nat)
   real(wp), contiguous, intent(in) :: brdr2(:, :, :, :, :) ! (3,nat,3,nat,nat)
   real(wp), contiguous, intent(in) :: brdr3(:, :, :, :, :, :, :) ! (3,nat,3,nat,3,nat,nat)
   real(wp), contiguous, intent(out) :: d3Kdr3(:, :, :, :, :, :, :, :) ! (3,nat,3,nat,3,nat,nat,nat)

   integer :: i, j, k, l, m
   integer :: alpha, beta, gamma
   integer :: delk, dell, delm
   real(wp), parameter :: a4 = 0.25_wp

   real(wp) :: v(3), r2
   real(wp) :: ai, aj, A, d, E, P, S
   real(wp) :: invf, invf3, invf5, invf7

   ! --- S-derivatives wrt (r2, ai, aj) ---
   real(wp) :: Sr, Saa, Sbb, Sa, Sb
   real(wp) :: Srr, Sra, Srb, Saa2, Sbb2, Sab
   real(wp) :: Srrr, Srra, Srrb, Sraa, Srab, Srbb
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

   ! helpers
   real(wp) :: term
   real(wp) :: f11, f22, f3_self

   I3 = 0.0_wp
   I3(1,1)=1.0_wp; I3(2,2)=1.0_wp; I3(3,3)=1.0_wp

   d3Kdr3(:, :, :, :, :, :, :, :) = 0.0_wp

   ! =========================
   ! Off-diagonal i>j, mirror
   ! =========================
   do i = 1, nat
      ai = brad(i)
      do j = 1, i-1
         aj = brad(j)

         v(:) = xyz(:, i) - xyz(:, j)
         r2   = dot_product(v, v)

         A = ai * aj
         d = a4 * r2 / A
         E = exp(-d)

         S = r2 + A * E
         invf  = 1.0_wp / sqrt(S)
         invf3 = invf*invf*invf
         invf5 = invf3*invf*invf
         invf7 = invf5*invf*invf

         ! ---- F(S) derivatives ----
         F1 = -0.5_wp * self%keps * invf3
         F2 =  0.75_wp * self%keps * invf5
         F3 = -1.875_wp * self%keps * invf7   ! -15/8

         ! ---- S derivatives wrt r2 (call it "r") and radii a=ai, b=aj ----
         ! Using: Sr = 1 - a4 E = P
         P  = 1.0_wp - a4*E
         Sr = P
         Srr  = (a4*a4 / A) * E
         Srrr = -(a4*a4*a4 / (A*A)) * E

         ! Sa, Sb (your Qi, Qj)
         Sa = aj * E * (1.0_wp + d)
         Sb = ai * E * (1.0_wp + d)

         ! second wrt radii
         Saa2 = aj * E * (d*d) / ai
         Sbb2 = ai * E * (d*d) / aj
         Sab  = E * (1.0_wp + d + d*d)

         ! mixed r–a, r–b
         Sra = -(a4/ai) * E * d
         Srb = -(a4/aj) * E * d

         ! mixed rr–a, rr–b
         Srra = -(a4*a4/(ai*A)) * E * (1.0_wp - d)
         Srrb = -(a4*a4/(aj*A)) * E * (1.0_wp - d)

         ! mixed r–aa, r–bb, r–ab
         Sraa = -a4 * E * d * (d - 2.0_wp) / (ai*ai)
         Srbb = -a4 * E * d * (d - 2.0_wp) / (aj*aj)
         Srab = -(a4/A) * E * d * (d - 1.0_wp)

         ! third wrt radii
         Saa3   = aj * E * d*d * (d - 3.0_wp) / (ai*ai)
         Sbb3   = ai * E * d*d * (d - 3.0_wp) / (aj*aj)
         Saa2b  = E * d*d * (d - 1.0_wp) / ai
         Sa2bb  = E * d*d * (d - 1.0_wp) / aj

         ! ---- K scalar partials via chain rule (K = F(S)) ----
         ! first
         Kr  = F1 * Sr
         Kai = F1 * Sa
         Kaj = F1 * Sb

         ! second
         Krr    = F2 * Sr*Sr + F1 * Srr
         Kaiai = F2 * Sa*Sa + F1 * Saa2
         Kajaj  = F2 * Sb*Sb + F1 * Sbb2
         Kaiaj  = F2 * Sa*Sb + F1 * Sab

         Kr_ai  = F2 * Sr*Sa + F1 * Sra
         Kr_aj  = F2 * Sr*Sb + F1 * Srb

         ! third: (r,r,r)
         Krrr = F3 * Sr*Sr*Sr + 3.0_wp*F2*Srr*Sr + F1*Srrr

         ! third: (r,r,a) and (r,r,b)
         Krr_ai = F3 * Sr*Sr*Sa + F2*(Srr*Sa + 2.0_wp*Sr*Sra) + F1*Srra
         Krr_aj = F3 * Sr*Sr*Sb + F2*(Srr*Sb + 2.0_wp*Sr*Srb) + F1*Srrb

         ! third: (r,a,a), (r,b,b), (r,a,b)
         Kr_aiai = F3 * Sr*Sa*Sa + F2*(Saa2*Sr + 2.0_wp*Sra*Sa) + F1*Sraa
         Kr_ajaj = F3 * Sr*Sb*Sb + F2*(Sbb2*Sr + 2.0_wp*Srb*Sb) + F1*Srbb
         Kr_aiaj = F3 * Sr*Sa*Sb + F2*(Sab*Sr + Sra*Sb + Srb*Sa) + F1*Srab

         ! third: radii-only
         Kaiaiai = F3 * Sa*Sa*Sa + 3.0_wp*F2*Saa2*Sa + F1*Saa3
         Kajajaj = F3 * Sb*Sb*Sb + 3.0_wp*F2*Sbb2*Sb + F1*Sbb3
         Kaiaiaj = F3 * Sa*Sa*Sb + F2*(Saa2*Sb + 2.0_wp*Sab*Sa) + F1*Saa2b
         Kaiajaj = F3 * Sa*Sb*Sb + F2*(Sbb2*Sa + 2.0_wp*Sab*Sb) + F1*Sa2bb

         ! ---- build v-tensors from r2-derivatives ----
         ! Kv = ∂K/∂v = 2*Kr * v
         Kv(:) = 2.0_wp * Kr * v(:)

         ! K_v_ai = 2*K_{r,ai} * v etc
         Kv_ai(:)   = 2.0_wp * Kr_ai  * v(:)
         Kv_aj(:)   = 2.0_wp * Kr_aj  * v(:)
         Kv_aiai(:) = 2.0_wp * Kr_aiai * v(:)
         Kv_ajaj(:) = 2.0_wp * Kr_ajaj * v(:)
         Kv_aiaj(:) = 2.0_wp * Kr_aiaj * v(:)

         ! K_vv_ai = 2*K_{r,ai}*I + 4*K_{rr,ai}*(v⊗v)
         Kvv_ai(:,:) = 2.0_wp*Kr_ai * I3(:,:) + 4.0_wp*Krr_ai * (spread(v,2,3)*spread(v,1,3))
         Kvv_aj(:,:) = 2.0_wp*Kr_aj * I3(:,:) + 4.0_wp*Krr_aj * (spread(v,2,3)*spread(v,1,3))

         ! K_vvv tensor:
         !   8*Krrr * v⊗v⊗v  + 4*Krr * sym( I⊗v )
         Kvvv(:,:,:) = 0.0_wp
         do alpha = 1,3
            do beta = 1,3
               do gamma = 1,3
                  Kvvv(alpha,beta,gamma) = 8.0_wp*Krrr * v(alpha)*v(beta)*v(gamma) &
                     + 4.0_wp*Krr * ( I3(alpha,beta)*v(gamma) + I3(alpha,gamma)*v(beta) + I3(beta,gamma)*v(alpha) )
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
                     dai_k = brdr(alpha,k,i)
                     daj_k = brdr(alpha,k,j)
                     do beta = 1,3
                        dai_l = brdr(beta,l,i)
                        daj_l = brdr(beta,l,j)
                        do gamma = 1,3
                           dai_m = brdr(gamma,m,i)
                           daj_m = brdr(gamma,m,j)

                           d2ai_kl = brdr2(alpha,k,beta,l,i)
                           d2ai_km = brdr2(alpha,k,gamma,m,i)
                           d2ai_lm = brdr2(beta,l,gamma,m,i)

                           d2aj_kl = brdr2(alpha,k,beta,l,j)
                           d2aj_km = brdr2(alpha,k,gamma,m,j)
                           d2aj_lm = brdr2(beta,l,gamma,m,j)

                           d3ai_klm = brdr3(alpha,k,beta,l,gamma,m,i)
                           d3aj_klm = brdr3(alpha,k,beta,l,gamma,m,j)

                           term = 0.0_wp

                           ! ---------- Term 1: K_pqr y_p,x y_q,y y_r,z ----------
                           ! vvv
                           if (delk/=0 .and. dell/=0 .and. delm/=0) then
                              term = term + real(delk*dell*delm,wp) * Kvvv(alpha,beta,gamma)
                           end if

                           ! vv-ai / vv-aj (three placements)
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

                           ! v-aa, v-bb, v-ab (three placements of which coord gives v)
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

                           ! aaa / bbb / mixed radii-only
                           term = term + Kaiaiai * (dai_k*dai_l*dai_m)
                           term = term + Kajajaj * (daj_k*daj_l*daj_m)

                           term = term + Kaiaiaj * (dai_k*dai_l*daj_m + dai_k*daj_l*dai_m + daj_k*dai_l*dai_m)
                           term = term + Kaiajaj * (daj_k*daj_l*dai_m + daj_k*dai_l*daj_m + dai_k*daj_l*daj_m)

                           ! ---------- Term 2: K_pq ( y_p,xy y_q,z + perms ) ----------
                           ! (xy)=(k,l), z=(m)
                           ! p = ai/aj second-derivative; q = v/ai/aj first-derivative
                           if (delm/=0) then
                              term = term + real(delm,wp) * d2ai_kl * Kv_ai(gamma)
                              term = term + real(delm,wp) * d2aj_kl * Kv_aj(gamma)
                           end if
                           term = term + d2ai_kl * ( Kaiai * dai_m + Kaiaj * daj_m )
                           term = term + d2aj_kl * ( Kaiaj * dai_m + Kajaj * daj_m )

                           ! (xy)=(k,m), z=(l)
                           if (dell/=0) then
                              term = term + real(dell,wp) * d2ai_km * Kv_ai(beta)
                              term = term + real(dell,wp) * d2aj_km * Kv_aj(beta)
                           end if
                           term = term + d2ai_km * ( Kaiai * dai_l + Kaiaj * daj_l )
                           term = term + d2aj_km * ( Kaiaj * dai_l + Kajaj * daj_l )

                           ! (xy)=(l,m), z=(k)
                           if (delk/=0) then
                              term = term + real(delk,wp) * d2ai_lm * Kv_ai(alpha)
                              term = term + real(delk,wp) * d2aj_lm * Kv_aj(alpha)
                           end if
                           term = term + d2ai_lm * ( Kaiai * dai_k + Kaiaj * daj_k )
                           term = term + d2aj_lm * ( Kaiaj * dai_k + Kajaj * daj_k )

                           ! ---------- Term 3: K_p y_p,xyz ----------
                           term = term + Kai * d3ai_klm + Kaj * d3aj_klm

                           ! write to tensor
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

      f11 = -self%keps / (ai*ai)                ! f'(a)
      f22 =  2.0_wp * self%keps / (ai*ai*ai)    ! f''(a)
      f3_self = -6.0_wp * self%keps / (ai**4)  ! f'''(a)

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
                        term = term + f3_self * (dai_k*dai_l*dai_m)
                        term = term + f22 * ( d2ai_kl*dai_m + d2ai_km*dai_l + d2ai_lm*dai_k )
                        term = term + f11 * d3ai_klm

                        d3Kdr3(alpha,k,beta,l,gamma,m,i,i) = d3Kdr3(alpha,k,beta,l,gamma,m,i,i) + term
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do

end subroutine compute_still_d3kdr3_full



end module tblite_solvation_kernel
