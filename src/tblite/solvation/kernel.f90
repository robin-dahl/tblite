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

!> Kernel types with derivatives up to fourth order as well as a general setup of the A-matrices
module tblite_solvation_kernel
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   use tblite_blas, only : gemv
   implicit none
   private

   public :: kernel_type, new_kernel
   public :: still_kernel, p16_kernel
   public :: kernel_enum, kernel_enum_type


   type :: kernel_enum_type
      integer :: still = 1
      integer :: p16 = 2
   end type kernel_enum_type

   type(kernel_enum_type), parameter :: kernel_enum = kernel_enum_type()

   ! Abstract base class for kernel types
   type, abstract :: kernel_type
      real(wp) :: keps  
   contains
      procedure(add_kernel_mat_interface), deferred :: add_kernel_mat
      procedure(add_kernel_deriv_interface), deferred :: add_kernel_deriv
      procedure(compute_kernel_deriv_interface), deferred :: compute_kernel_deriv
   end type kernel_type


   abstract interface
      !> Add kernel contributions to interaction matrix
      subroutine add_kernel_mat_interface(self, nat, xyz, brad, amat)
         import :: kernel_type, wp
         !> Kernel instance
         class(kernel_type), intent(in) :: self
         !> Number of atoms
         integer, intent(in) :: nat
         !> Cartesian coordinates
         real(wp), intent(in) :: xyz(:, :)
         !> Born radii
         real(wp), intent(in) :: brad(:)
         !> Interaction matrix
         real(wp), intent(inout) :: amat(:, :)
      end subroutine add_kernel_mat_interface

      !> Add kernel derivative contributions to energy and gradient
      subroutine add_kernel_deriv_interface(self, nat, xyz, qat, brad, brdr, &
            & energy, gradient)
         import :: kernel_type, wp
         !> Kernel instance
         class(kernel_type), intent(in) :: self
         !> Number of atoms
         integer, intent(in) :: nat
         !> Cartesian coordinates
         real(wp), intent(in) :: xyz(:, :)
         !> Atomic partial charges
         real(wp), intent(in) :: qat(:)
         !> Born radii
         real(wp), intent(in) :: brad(:)
         !> Derivative of Born radii w.r.t. cartesian coordinates
         real(wp), contiguous, intent(in) :: brdr(:, :, :)
         !> Total Born solvation energy
         real(wp), intent(out) :: energy
         !> Derivatives of Born solvation energy
         real(wp), contiguous, intent(inout) :: gradient(:, :)
      end subroutine add_kernel_deriv_interface

      !> Compute kernel derivatives only (without charge multiplication)
      subroutine compute_kernel_deriv_interface(self, nat, xyz, brad, brdr, &
            & kernel_grad_spatial, kernel_grad_born)
         import :: kernel_type, wp
         !> Kernel instance
         class(kernel_type), intent(in) :: self
         !> Number of atoms
         integer, intent(in) :: nat
         !> Cartesian coordinates
         real(wp), intent(in) :: xyz(:, :)
         !> Born radii
         real(wp), intent(in) :: brad(:)
         !> Derivative of Born radii w.r.t. cartesian coordinates
         real(wp), contiguous, intent(in) :: brdr(:, :, :)
         !> Spatial kernel gradient (3, nat, nat)
         real(wp), contiguous, intent(out) :: kernel_grad_spatial(:, :, :)
         !> Born radii kernel gradient (nat, nat)
         real(wp), contiguous, intent(out) :: kernel_grad_born(:, :)
      end subroutine compute_kernel_deriv_interface
   end interface

   type, extends(kernel_type) :: still_kernel
      contains
         procedure :: add_kernel_mat => add_still_mat
         procedure :: add_kernel_deriv => add_still_deriv
         procedure :: compute_kernel_deriv => compute_still_deriv
   end type still_kernel

   type, extends(kernel_type) :: p16_kernel
      contains
         procedure :: add_kernel_mat => add_p16_mat
         procedure :: add_kernel_deriv => add_p16_deriv
         procedure :: compute_kernel_deriv => compute_p16_deriv
   end type p16_kernel


real(wp), parameter :: zetaP16 = 1.028_wp
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




subroutine add_p16_mat(self, nat, xyz, brad, Amat)
   !> Kernel instance
   class(p16_kernel), intent(in) :: self 
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Born radii
   real(wp), intent(in) :: brad(:)
   !> Interaction matrix
   real(wp), intent(inout) :: Amat(:, :)

   integer :: iat, jat
   real(wp) :: r1, ab, arg, fgb, dfgb, bp, vec(3)

   ! omp parallel do default(none) shared(Amat, ntpair, ppind, ddpair, brad, keps) &
   ! omp private(kk, iat, jat, r1, ab, arg, fgb, dfgb)
   do iat = 1, nat
      do jat = 1, iat - 1
         vec(:) = xyz(:, iat) - xyz(:, jat)
         r1 = norm2(vec)

         ab = sqrt(brad(iat) * brad(jat))
         arg = ab / (ab + zetaP16o16*r1) ! ab / (1 + ζR/(16·ab))
         arg = arg * arg ! ab / (1 + ζR/(16·ab))²
         arg = arg * arg ! ab / (1 + ζR/(16·ab))⁴
         arg = arg * arg ! ab / (1 + ζR/(16·ab))⁸
         arg = arg * arg ! ab / (1 + ζR/(16·ab))¹⁶
         fgb = r1 + ab*arg
         dfgb = 1.0_wp / fgb

         Amat(iat, jat) = self%keps*dfgb + Amat(iat, jat)
         Amat(jat, iat) = self%keps*dfgb + Amat(jat, iat)
      enddo
      ! self-energy part
      bp = 1.0_wp/brad(iat)
      Amat(iat, iat) = Amat(iat, iat) + self%keps*bp
   enddo

end subroutine add_p16_mat


subroutine add_p16_deriv(self, nat, xyz, qat, &
      & brad, brdr, energy, gradient)
   !> Kernel instance
   class(p16_kernel), intent(in) :: self 
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Born radii
   real(wp), intent(in) :: brad(:)
   !> Derivative of Born radii w.r.t. cartesian coordinates
   real(wp), contiguous, intent(in) :: brdr(:, :, :)
   !> Total Born solvation energy
   real(wp), intent(out) :: energy
   !> Deriatives of Born solvation energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)

   integer :: iat, jat
   real(wp) :: vec(3), r2, r1, ab, arg1, arg16, qq, fgb, dfgb, dfgb2, egb
   real(wp) :: dEdbri, dEdbrj, dG(3), ap, bp, dS(3, 3)
   real(wp), allocatable :: dEdbr(:)

   allocate(dEdbr(nat), source = 0.0_wp )

   egb = 0._wp
   dEdbr(:) = 0._wp

   ! GB energy and gradient
   ! omp parallel do default(none) reduction(+:egb, gradient, dEdbr) &
   ! omp private(iat, jat, vec, r1, r2, ab, arg1, arg16, fgb, dfgb, dfgb2, ap, &
   ! omp& bp, qq, dEdbri, dEdbrj, dG, dS) &
   ! omp shared(keps, qat, ntpair, ddpair, ppind, brad)
   do iat = 1, nat
      do jat = 1, iat - 1
         vec(:) = xyz(:, iat) - xyz(:, jat)
         r1 = norm2(vec)
         r2 = r1*r1

         qq = qat(iat)*qat(jat)

         ab = sqrt(brad(iat) * brad(jat))
         arg1 = ab / (ab + zetaP16o16*r1) ! 1 / (1 + ζR/(16·ab))
         arg16 = arg1 * arg1 ! 1 / (1 + ζR/(16·ab))²
         arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))⁴
         arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))⁸
         arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))¹⁶

         fgb = r1 + ab*arg16
         dfgb = 1.0_wp / fgb
         dfgb2 = dfgb * dfgb

         egb = egb + qq*self%keps*dfgb

         ! (1 - ζ/(1 + Rζ/(16 ab))^17)/(R + ab/(1 + Rζ/(16 ab))¹⁶)²
         ap = (1.0_wp - zetaP16 * arg1 * arg16) * dfgb2
         dG(:) = ap * vec * self%keps / r1 * qq
         gradient(:, iat) = gradient(:, iat) - dG
         gradient(:, jat) = gradient(:, jat) + dG

         ! -(Rζ/(2·ab²·(1 + Rζ/(16·ab))¹⁷) + 1/(2·ab·(1 + Rζ/(16·ab))¹⁶))/(R + ab/(1 + Rζ/(16·ab))¹⁶)²
         bp = -0.5_wp*(r1 * zetaP16 / ab * arg1 + 1.0_wp) / ab * arg16 * dfgb2
         dEdbri = brad(jat) * bp * self%keps * qq
         dEdbrj = brad(iat) * bp * self%keps * qq
         dEdbr(iat) = dEdbr(iat) + dEdbri
         dEdbr(jat) = dEdbr(jat) + dEdbrj

      end do

      ! self-energy part
      bp = 1._wp/brad(iat)
      qq = qat(iat)*bp
      egb = egb + 0.5_wp*qat(iat)*qq*self%keps
      dEdbri = -0.5_wp*self%keps*qq*bp
      dEdbr(iat) = dEdbr(iat) + dEdbri*qat(iat)
      !gradient = gradient + brdr(:, :, i) * dEdbri*qat(i)
   enddo

   ! contract with the Born radii derivatives
   call gemv(brdr, dEdbr, gradient, beta=1.0_wp)

   energy = egb

end subroutine add_p16_deriv


!> Compute P16 kernel derivatives only (without charge multiplication)
subroutine compute_p16_deriv(self, nat, xyz, brad, brdr, &
      & kernel_grad_spatial, kernel_grad_born)
   !> Kernel instance
   class(p16_kernel), intent(in) :: self 
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Born radii
   real(wp), intent(in) :: brad(:)
   !> Derivative of Born radii w.r.t. cartesian coordinates
   real(wp), contiguous, intent(in) :: brdr(:, :, :)
   !> Spatial kernel gradient (3, nat, nat)
   real(wp), contiguous, intent(out) :: kernel_grad_spatial(:, :, :)
   !> Born radii kernel gradient (nat, nat)
   real(wp), contiguous, intent(out) :: kernel_grad_born(:, :)

   integer :: iat, jat
   real(wp) :: vec(3), r1, ab, arg1, arg16, ap, bp
   real(wp), allocatable :: dKdbr(:)

   allocate(dKdbr(nat), source = 0.0_wp)
   
   kernel_grad_spatial(:, :, :) = 0.0_wp
   kernel_grad_born(:, :) = 0.0_wp

   ! Compute kernel derivatives (without charges)
   do iat = 1, nat
      do jat = 1, iat - 1
         vec(:) = xyz(:, iat) - xyz(:, jat)
         r1 = norm2(vec)

         ab = sqrt(brad(iat) * brad(jat))
         arg1 = ab / (ab + zetaP16o16*r1)
         arg16 = arg1 * arg1
         arg16 = arg16 * arg16
         arg16 = arg16 * arg16
         arg16 = arg16 * arg16

         ! Spatial derivative (similar to add_p16_deriv but without charges)
         ap = (1.0_wp - zetaP16 * arg1 * arg16) / (r1 + ab*arg16)**2 * self%keps
         
         kernel_grad_spatial(:, iat, jat) = ap * vec / r1
         kernel_grad_spatial(:, jat, iat) = -ap * vec / r1

         ! Born radii derivative
         bp = -0.5_wp*(r1 * zetaP16 / ab * arg1 + 1.0_wp) / ab * arg16 / (r1 + ab*arg16)**2 * self%keps
         
         kernel_grad_born(iat, jat) = bp * brad(jat)
         kernel_grad_born(jat, iat) = bp * brad(iat)
         
         dKdbr(iat) = dKdbr(iat) + bp * brad(jat)
         dKdbr(jat) = dKdbr(jat) + bp * brad(iat)
      enddo

      ! Self-energy kernel derivative
      bp = self%keps / brad(iat)
      dKdbr(iat) = dKdbr(iat) - bp / brad(iat)
   enddo

   ! Add contribution from Born radii position dependence
   do iat = 1, nat
      do jat = 1, nat
         kernel_grad_spatial(:, jat, iat) = kernel_grad_spatial(:, jat, iat) &
            & + brdr(:, jat, iat) * dKdbr(iat)
      enddo
   enddo

end subroutine compute_p16_deriv


pure subroutine add_still_mat(self, nat, xyz, brad, Amat)
   !> Kernel instance
   class(still_kernel), intent(in) :: self 
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Born radii
   real(wp), intent(in) :: brad(:)
   !> Interaction matrix
   real(wp), intent(inout) :: Amat(:, :)

   integer  :: i, j
   real(wp), parameter :: a13=1.0_wp/3.0_wp
   real(wp), parameter :: a4=0.25_wp
   real(wp) :: aa, vec(3), r1, r2, bp
   real(wp) :: dd, expd, fgb2, dfgb

   do i = 1, nat
      do j = 1, i - 1
         vec(:) = xyz(:, i) - xyz(:, j)
         r1 = norm2(vec)
         r2 = r1*r1

         aa = brad(i)*brad(j)
         dd = a4*r2/aa
         expd = exp(-dd)
         fgb2 = r2+aa*expd
         dfgb = 1.0_wp/sqrt(fgb2)
         Amat(i, j) = self%keps*dfgb + Amat(i, j)
         Amat(j, i) = self%keps*dfgb + Amat(j, i)
      enddo

      ! self-energy part
      bp = 1._wp/brad(i)
      Amat(i, i) = Amat(i, i) + self%keps*bp
   enddo

end subroutine add_still_mat


subroutine add_still_deriv(self, nat, xyz, qat, &
      & brad, brdr, energy, gradient)
   !> Kernel instance
   class(still_kernel), intent(in) :: self 
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Born radii
   real(wp), intent(in) :: brad(:)
   !> Derivative of Born radii w.r.t. cartesian coordinates
   real(wp), contiguous, intent(in) :: brdr(:, :, :)
   !> Total Born solvation energy
   real(wp), intent(out) :: energy
   !> Deriatives of Born solvation energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)

   integer :: i, j
   real(wp), parameter :: a13=1._wp/3._wp
   real(wp), parameter :: a4=0.25_wp
   real(wp) :: aa, r2, fgb2
   real(wp) :: qq, dd, expd, dfgb, dfgb2, dfgb3, egb, ap, bp
   real(wp) :: grddbi, grddbj
   real(wp) :: dr(3), r1, vec(3)
   real(wp), allocatable :: grddb(:)

   allocate(grddb(nat), source = 0.0_wp )

   egb = 0._wp
   grddb(:) = 0._wp

   ! GB energy and gradient

   ! compute energy and fgb direct and radii derivatives
   do i = 1, nat
      do j = 1, i - 1
         vec(:) = xyz(:, i) - xyz(:, j)
         r1 = norm2(vec)
         r2 = r1*r1

         ! dielectric scaling of the charges
         qq = qat(i)*qat(j)
         aa = brad(i)*brad(j)
         dd = a4*r2/aa
         expd = exp(-dd)
         fgb2 = r2+aa*expd
         dfgb2 = 1._wp/fgb2
         dfgb = sqrt(dfgb2)
         dfgb3 = dfgb2*dfgb*self%keps

         egb = egb + qq*self%keps*dfgb

         ap = (1._wp-a4*expd)*dfgb3
         dr = ap*vec
         gradient(:, i) = gradient(:, i) - dr*qq
         gradient(:, j) = gradient(:, j) + dr*qq

         bp = -0.5_wp*expd*(1._wp+dd)*dfgb3
         grddbi = brad(j)*bp
         grddbj = brad(i)*bp
         grddb(i) = grddb(i) + grddbi*qq
         grddb(j) = grddb(j) + grddbj*qq

      enddo

      ! self-energy part
      bp = 1._wp/brad(i)
      qq = qat(i)*bp
      egb = egb + 0.5_wp*qat(i)*qq*self%keps
      grddbi = -0.5_wp*self%keps*qq*bp
      grddb(i) = grddb(i) + grddbi*qat(i)
   enddo

   ! contract with the Born radii derivatives
   call gemv(brdr, grddb, gradient, beta=1.0_wp)

   energy = egb

end subroutine add_still_deriv


!> Compute kernel derivatives only (without charge multiplication)
subroutine compute_still_deriv(self, nat, xyz, brad, brdr, &
      & kernel_grad_spatial, kernel_grad_born)
   !> Kernel instance
   class(still_kernel), intent(in) :: self 
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Born radii
   real(wp), intent(in) :: brad(:)
   !> Derivative of Born radii w.r.t. cartesian coordinates
   real(wp), contiguous, intent(in) :: brdr(:, :, :)
   !> Spatial kernel gradient (3, nat, nat)
   real(wp), contiguous, intent(out) :: kernel_grad_spatial(:, :, :)
   !> Born radii kernel gradient (nat, nat)
   real(wp), contiguous, intent(out) :: kernel_grad_born(:, :)

   integer :: i, j
   real(wp), parameter :: a4=0.25_wp
   real(wp) :: aa, r2, fgb2
   real(wp) :: dd, expd, dfgb, dfgb2, dfgb3, ap, bp
   real(wp) :: r1, vec(3)
   real(wp), allocatable :: dKdbr(:)

   allocate(dKdbr(nat), source = 0.0_wp)
   
   kernel_grad_spatial(:, :, :) = 0.0_wp
   kernel_grad_born(:, :) = 0.0_wp

   ! Compute kernel derivatives (without charges)
   do i = 1, nat
      do j = 1, i - 1
         vec(:) = xyz(:, i) - xyz(:, j)
         r1 = norm2(vec)
         r2 = r1*r1

         aa = brad(i)*brad(j)
         dd = a4*r2/aa
         expd = exp(-dd)
         fgb2 = r2+aa*expd
         dfgb2 = 1._wp/fgb2
         dfgb = sqrt(dfgb2)
         dfgb3 = dfgb2*dfgb*self%keps

         ! Spatial derivative of kernel: ∂(κ/f_GB)/∂r_ij
         ap = (1._wp-a4*expd)*dfgb3
         
         ! Store directional derivative: ap * vec
         kernel_grad_spatial(:, i, j) = ap * vec
         kernel_grad_spatial(:, j, i) = -ap * vec

         ! Born radii derivative: ∂(κ/f_GB)/∂a_ij
         bp = -0.5_wp*expd*(1._wp+dd)*dfgb3
         
         ! Store kernel derivatives w.r.t. Born radii
         kernel_grad_born(i, j) = bp * brad(j)  ! ∂/∂r_{B,i}
         kernel_grad_born(j, i) = bp * brad(i)  ! ∂/∂r_{B,j}
         
         ! Accumulate for Born radii chain rule
         dKdbr(i) = dKdbr(i) + bp * brad(j)
         dKdbr(j) = dKdbr(j) + bp * brad(i)
      enddo

      ! Self-energy kernel derivative
      bp = self%keps / brad(i)
      dKdbr(i) = dKdbr(i) - bp / brad(i)
   enddo

   ! Add contribution from Born radii position dependence
   do i = 1, nat
      do j = 1, nat
         kernel_grad_spatial(:, j, i) = kernel_grad_spatial(:, j, i) &
            & + brdr(:, j, i) * dKdbr(i)
      enddo
   enddo

end subroutine compute_still_deriv



end module tblite_solvation_kernel