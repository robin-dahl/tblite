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
   public :: still_kernel, p16_kernel, coulomb_kernel
   public :: kernel_enum, kernel_enum_type
   public :: compute_kernel_d3Kdr3_ij 
   public :: compute_kernel_dKdr, compute_kernel_d2Kdr2, compute_kernel_d3Kdr3, compute_kernel_d4Kdr4
   public :: compute_coulomb_dKdr

   type :: kernel_enum_type
      integer :: still = 1
      integer :: p16   = 2
      integer :: coulomb = 3
   end type kernel_enum_type

   type(kernel_enum_type), parameter :: kernel_enum = kernel_enum_type()

   ! Abstract base class for kernel types
   type, abstract :: kernel_type
      real(wp) :: keps
   contains
      procedure(add_kernel_mat_interface),           deferred :: add_kernel_mat
      procedure(add_kernel_deriv_interface),         deferred :: add_kernel_deriv
      procedure(compute_kernel_dKdr_interface),      deferred :: compute_kernel_dKdr
      procedure(compute_kernel_d2Kdr2_interface),    deferred :: compute_kernel_d2Kdr2
      procedure(compute_kernel_d3Kdr3_ij_interface), deferred :: compute_kernel_d3Kdr3_ij
      procedure(compute_kernel_d3Kdr3_interface),    deferred :: compute_kernel_d3Kdr3
      procedure(compute_kernel_d4Kdr4_interface),    deferred :: compute_kernel_d4Kdr4
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

      subroutine compute_kernel_d3Kdr3_ij_interface(self, nat, xyz, brad, brdr, brdr2, brdr3, i, j, d3Kdr3_ij)
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
      end subroutine compute_kernel_d3Kdr3_ij_interface

      !> Element-wise kernel derivative interfaces
      subroutine compute_kernel_dKdr_interface(self, nat, xyz, brad, i, j, k, alpha, dKdr_elem, brdr)
         import :: kernel_type, wp
         class(kernel_type), intent(in) :: self
         integer, intent(in) :: nat
         real(wp), intent(in) :: xyz(:, :)
         real(wp), intent(in) :: brad(:)
         integer, intent(in) :: i, j
         integer, intent(in) :: k, alpha
         real(wp), intent(out) :: dKdr_elem
         real(wp), contiguous, intent(in), optional :: brdr(:, :, :)
      end subroutine compute_kernel_dKdr_interface

      subroutine compute_kernel_d2Kdr2_interface(self, nat, xyz, brad, i, j, k, alpha, l, beta, d2K_elem, brdr, brdr2)
         import :: kernel_type, wp
         class(kernel_type), intent(in) :: self
         integer, intent(in) :: nat
         real(wp), intent(in) :: xyz(:, :)
         real(wp), intent(in) :: brad(:)
         integer, intent(in) :: i, j
         integer, intent(in) :: k, alpha
         integer, intent(in) :: l, beta
         real(wp), intent(out) :: d2K_elem
         real(wp), contiguous, intent(in), optional :: brdr(:, :, :)
         real(wp), contiguous, intent(in), optional :: brdr2(:, :, :, :, :)
      end subroutine compute_kernel_d2Kdr2_interface

      subroutine compute_kernel_d3Kdr3_interface(self, nat, xyz, brad, i, j, k, alpha, l, beta, m, gamma, d3K_elem, &
            & brdr, brdr2, brdr3)
         import :: kernel_type, wp
         class(kernel_type), intent(in) :: self
         integer, intent(in) :: nat
         real(wp), intent(in) :: xyz(:, :)
         real(wp), intent(in) :: brad(:)
         integer, intent(in) :: i, j
         integer, intent(in) :: k, alpha
         integer, intent(in) :: l, beta
         integer, intent(in) :: m, gamma
         real(wp), intent(out) :: d3K_elem
         real(wp), contiguous, intent(in), optional :: brdr(:, :, :)
         real(wp), contiguous, intent(in), optional :: brdr2(:, :, :, :, :)
         real(wp), contiguous, intent(in), optional :: brdr3(:, :, :, :, :, :, :)
      end subroutine compute_kernel_d3Kdr3_interface

      subroutine compute_kernel_d4Kdr4_interface(self, nat, xyz, brad, i, j, k, alpha, l, beta, m, gamma, n, delta, d4K_elem, &
            & brdr, brdr2, brdr3, brdr4)
         import :: kernel_type, wp
         class(kernel_type), intent(in) :: self
         integer, intent(in) :: nat
         real(wp), intent(in) :: xyz(:, :)
         real(wp), intent(in) :: brad(:)
         integer, intent(in) :: i, j
         integer, intent(in) :: k, alpha, l, beta, m, gamma, n, delta
         real(wp), intent(out) :: d4K_elem
         real(wp), contiguous, intent(in), optional :: brdr(:, :, :)
         real(wp), contiguous, intent(in), optional :: brdr2(:, :, :, :, :)
         real(wp), contiguous, intent(in), optional :: brdr3(:, :, :, :, :, :, :)
         real(wp), contiguous, intent(in), optional :: brdr4(:, :, :, :, :, :, :, :, :)
      end subroutine compute_kernel_d4Kdr4_interface
      
   end interface

   type, extends(kernel_type) :: still_kernel
   contains
      procedure :: add_kernel_mat   => add_still_mat
      procedure :: add_kernel_deriv => add_still_deriv
      procedure :: compute_kernel_dKdr => compute_still_dKdr
      procedure :: compute_kernel_d2Kdr2 => compute_still_d2Kdr2
      procedure :: compute_kernel_d3Kdr3_ij => compute_still_d3Kdr3_ij
      procedure :: compute_kernel_d3Kdr3 => compute_still_d3Kdr3
      procedure :: compute_kernel_d4Kdr4 => compute_still_d4Kdr4
   end type still_kernel

   type, extends(kernel_type) :: p16_kernel
   contains
      procedure :: add_kernel_mat   => add_p16_mat
      procedure :: add_kernel_deriv => add_p16_deriv
      procedure :: compute_kernel_dKdr => compute_p16_dKdr
      procedure :: compute_kernel_d2Kdr2 => compute_p16_d2Kdr2
      procedure :: compute_kernel_d3Kdr3 => compute_p16_d3Kdr3
      procedure :: compute_kernel_d3Kdr3_ij => compute_p16_d3Kdr3_ij
      procedure :: compute_kernel_d4Kdr4 => compute_p16_d4Kdr4
   end type p16_kernel

   type, extends(kernel_type) :: coulomb_kernel
   contains
      procedure :: add_kernel_mat => add_coulomb_mat
      procedure :: add_kernel_deriv => add_coulomb_deriv
      procedure :: compute_kernel_dKdr => compute_coulomb_dKdr
      procedure :: compute_kernel_d2Kdr2 => compute_coulomb_d2Kdr2
      procedure :: compute_kernel_d3Kdr3_ij => compute_coulomb_d3Kdr3_ij
      procedure :: compute_kernel_d3Kdr3 => compute_coulomb_d3Kdr3
      procedure :: compute_kernel_d4Kdr4 => compute_coulomb_d4Kdr4
   end type

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
   case(kernel_enum%coulomb)
      allocate(coulomb_kernel :: kernel)
   case default
      allocate(p16_kernel :: kernel)
   end select

   kernel%keps = keps
end function new_kernel


subroutine compute_kernel_d3Kdr3_ij(kernel_id, keps, nat, xyz, brad, brdr, brdr2, brdr3, i, j, d3Kdr3_ij)
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
   call kernel%compute_kernel_d3Kdr3_ij(nat, xyz, brad, brdr, brdr2, brdr3, i, j, d3Kdr3_ij)
end subroutine compute_kernel_d3Kdr3_ij


!> Element-wise convenience dispatcher (switches by kernel enum, returns single derivative element)
subroutine compute_kernel_dKdr(kernel_id, keps, nat, xyz, brad, i, j, k, alpha, dKdr_elem, brdr)
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: brad(:)
   integer, intent(in) :: i, j
   integer, intent(in) :: k, alpha
   real(wp), intent(out) :: dKdr_elem
   real(wp), contiguous, intent(in), optional :: brdr(:, :, :)

   class(kernel_type), allocatable :: kernel

   kernel = new_kernel(kernel_id, keps)
   if (present(brdr)) then
      call kernel%compute_kernel_dKdr(nat, xyz, brad, i, j, k, alpha, dKdr_elem, brdr=brdr)
   else 
      call kernel%compute_kernel_dKdr(nat, xyz, brad, i, j, k, alpha, dKdr_elem)
   end if 

end subroutine compute_kernel_dKdr

!> Element-wise convenience dispatcher (switches by kernel enum, returns single derivative element)
subroutine compute_kernel_d2Kdr2(kernel_id, keps, nat, xyz, brad, i, j, k, alpha, l, beta, d2K_elem, &
      & brdr, brdr2)
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: brad(:)
   integer, intent(in) :: i, j
   integer, intent(in) :: k, alpha
   integer, intent(in) :: l, beta
   real(wp), intent(out) :: d2K_elem
   real(wp), contiguous, intent(in), optional :: brdr(:, :, :)
   real(wp), contiguous, intent(in), optional :: brdr2(:, :, :, :, :)

   class(kernel_type), allocatable :: kernel

   kernel = new_kernel(kernel_id, keps)
   if (present(brdr) .and. present(brdr2)) then
      call kernel%compute_kernel_d2Kdr2(nat, xyz, brad, i, j, k, alpha, l, beta, d2K_elem, &
         & brdr=brdr, brdr2=brdr2)
   else 
      call kernel%compute_kernel_d2Kdr2(nat, xyz, brad, i, j, k, alpha, l, beta, d2K_elem)
   end if
end subroutine compute_kernel_d2Kdr2

!> Element-wise convenience dispatcher (switches by kernel enum, returns single derivative element)
subroutine compute_kernel_d3Kdr3(kernel_id, keps, nat, xyz, brad, i, j, k, alpha, l, beta, m, gamma, d3K_elem, &
      & brdr, brdr2, brdr3)
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: brad(:)
   integer, intent(in) :: i, j
   integer, intent(in) :: k, alpha
   integer, intent(in) :: l, beta
   integer, intent(in) :: m, gamma
   real(wp), intent(out) :: d3K_elem
   real(wp), contiguous, intent(in), optional :: brdr(:, :, :)
   real(wp), contiguous, intent(in), optional :: brdr2(:, :, :, :, :)
   real(wp), contiguous, intent(in), optional :: brdr3(:, :, :, :, :, :, :)

   class(kernel_type), allocatable :: kernel

   kernel = new_kernel(kernel_id, keps)
   if (present(brdr) .and. present(brdr2) .and. present(brdr3)) then
      call kernel%compute_kernel_d3Kdr3(nat, xyz, brad, i, j, k, alpha, l, beta, m, gamma, d3K_elem, &
         & brdr=brdr, brdr2=brdr2, brdr3=brdr3)
   else
      call kernel%compute_kernel_d3Kdr3(nat, xyz, brad, i, j, k, alpha, l, beta, m, gamma, d3K_elem)
   end if
end subroutine compute_kernel_d3Kdr3

!> Element-wise convenience dispatcher (switches by kernel enum, returns single derivative element)
subroutine compute_kernel_d4Kdr4(kernel_id, keps, nat, xyz, brad, i, j, k, alpha, l, beta, m, gamma, n, delta, d4K_elem, &
      & brdr, brdr2, brdr3, brdr4)
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: keps
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: brad(:)
   integer, intent(in) :: i, j
   integer, intent(in) :: k, alpha, l, beta, m, gamma, n, delta
   real(wp), intent(out) :: d4K_elem
   real(wp), contiguous, intent(in), optional :: brdr(:, :, :)
   real(wp), contiguous, intent(in), optional :: brdr2(:, :, :, :, :)
   real(wp), contiguous, intent(in), optional :: brdr3(:, :, :, :, :, :, :)
   real(wp), contiguous, intent(in), optional :: brdr4(:, :, :, :, :, :, :, :, :)

   class(kernel_type), allocatable :: kernel

   kernel = new_kernel(kernel_id, keps)
   if (present(brdr) .and. present(brdr2) .and. present(brdr3) .and. present(brdr4)) then
      call kernel%compute_kernel_d4Kdr4(nat, xyz, brad, i, j, k, alpha, l, beta, m, gamma, n, delta, d4K_elem, &
         & brdr=brdr, brdr2=brdr2, brdr3=brdr3, brdr4=brdr4)
   else
      call kernel%compute_kernel_d4Kdr4(nat, xyz, brad, i, j, k, alpha, l, beta, m, gamma, n, delta, d4K_elem)
   end if
end subroutine compute_kernel_d4Kdr4




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Explicit kernel implementations with derivatives below
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

subroutine compute_still_dKdr(self, nat, xyz, brad, i, j, k, alpha, dKdr_elem, brdr)
   ! Element-wise first derivative for the Still kernel:
   !
   !   dKdr_elem = d K_ij / d r_{k,alpha}
   !
   ! Philosophy:
   ! - compute only the requested derivative entry (no (3,nat) storage)
   ! - keep current brdr input (computed as full tensor)
   ! - BUT: read the needed brdr element(s) into local scalars at the top, so later you can
   !   swap those reads to calls like brdr_elem(k,alpha,i) without touching the algebra.

   class(still_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                 ! (3,nat)
   real(wp), intent(in) :: brad(:)                   ! (nat)
   integer, intent(in) :: i, j
   integer, intent(in) :: k, alpha                   ! derivative site and component
   real(wp), intent(out) :: dKdr_elem
   real(wp), contiguous, intent(in), optional :: brdr(:, :, :) ! (3,nat,nat) 

   real(wp), parameter :: a4 = 0.25_wp
   real(wp), parameter :: tiny = 1.0e-30_wp

   real(wp) :: rvec(3), r2
   real(wp) :: A, d, E, f2, invf, invf3
   real(wp) :: pref_pos
   real(wp) :: dK_dai, dK_daj

   ! --- Cache the only brdr entries we might need for THIS derivative entry ---
   ! Once brdr becomes element-wise, you only replace these assignments.
   real(wp) :: brdr_ki_alpha, brdr_kj_alpha

   dKdr_elem = 0.0_wp

   ! Safety: alpha must be 1..3
   if (alpha < 1 .or. alpha > 3) return

   ! Cache brdr component(s) needed for chain rule term for this (k,alpha)
   ! (These are harmless even if i==j; still well-defined.)
   if (present(brdr)) then
      brdr_ki_alpha = brdr(alpha, k, i)
      brdr_kj_alpha = brdr(alpha, k, j)
   else
      brdr_ki_alpha = 0.0_wp
      brdr_kj_alpha = 0.0_wp
   end if

   if (i == j) then
      ! Diagonal: K_ii = keps / a_i
      ! dK_ii / d r_{k,alpha} = (-keps/a_i^2) * (d a_i / d r_{k,alpha})
      dKdr_elem = (-self%keps / (brad(i)*brad(i))) * brdr_ki_alpha
      return
   end if

   ! Off-diagonal element (i,j)
   rvec(:) = xyz(:, i) - xyz(:, j)
   r2      = dot_product(rvec, rvec)

   A = brad(i) * brad(j)
   if (abs(A) <= tiny) return   ! avoid divide-by-zero if radii are pathological

   d = a4 * r2 / A
   E = exp(-d)

   f2    = r2 + A * E
   if (f2 <= tiny) return       ! avoid sqrt issues
   invf  = 1.0_wp / sqrt(f2)
   invf3 = invf * invf * invf

   ! -------- Explicit coordinate dependence (radii held fixed) --------
   ! pref_pos = -(1 - a4*E) * (keps * invf^3)
   pref_pos = -(1.0_wp - a4*E) * (self%keps * invf3)

   ! Only k==i or k==j contributes via explicit coordinate dependence
   if (k == i) then
      dKdr_elem = dKdr_elem + pref_pos * rvec(alpha)
   else if (k == j) then
      dKdr_elem = dKdr_elem - pref_pos * rvec(alpha)
   end if

   ! -------- Chain rule via Born radii --------
   dK_dai = -0.5_wp * self%keps * E * (1.0_wp + d) * brad(j) * invf3
   dK_daj = -0.5_wp * self%keps * E * (1.0_wp + d) * brad(i) * invf3

   ! For this (k,alpha), only da_i/dr_{k,alpha} and da_j/dr_{k,alpha} are needed:
   dKdr_elem = dKdr_elem + dK_dai * brdr_ki_alpha + dK_daj * brdr_kj_alpha

end subroutine compute_still_dKdr

subroutine compute_still_d2Kdr2(self, nat, xyz, brad, i, j, k, alpha, l, beta, d2K_elem, &
      & brdr, brdr2)
   ! Element-wise 2nd derivative for Still kernel:
   !   d2K_elem = d^2 K_ij / ( d r_{k,alpha} d r_{l,beta} )
   !
   ! Same philosophy as before:
   ! - compute only one entry, no (3,nat,3,nat) storage
   ! - cache ONLY the needed Born-radius derivative entries (brdr, brdr2) as scalars
   !   so you can later replace them with element-wise providers easily.

   class(still_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                         ! (3,nat)
   real(wp), intent(in) :: brad(:)                           ! (nat)
   integer, intent(in) :: i, j
   integer, intent(in) :: k, alpha
   integer, intent(in) :: l, beta
   real(wp), intent(out) :: d2K_elem
   real(wp), contiguous, intent(in), optional :: brdr(:, :, :)         ! (3,nat,nat)
   real(wp), contiguous, intent(in), optional :: brdr2(:, :, :, :, :)  ! (3,nat,3,nat,nat)

   integer :: delk, dell
   real(wp), parameter :: a4 = 0.25_wp
   real(wp), parameter :: tiny = 1.0e-30_wp

   real(wp) :: v(3), r2
   real(wp) :: ai, aj, A, d, E, P, S
   real(wp) :: invf, invf3, invf5
   real(wp) :: C, dCdr2, dC_dai, dC_daj
   real(wp) :: Qi, Qj
   real(wp) :: dK_dai, dK_daj
   real(wp) :: d2K_dai2, d2K_daj2, d2K_daida_j
   real(wp) :: Iab, vv_ab, Hvv_ab
   real(wp) :: coef1, coef2

   ! --- Cache ONLY the needed Born-radius derivative entries for this (k,alpha,l,beta) ---
   ! brdr(alpha,k,atom), brdr(beta,l,atom), and brdr2(alpha,k,beta,l,atom)
   real(wp) :: dai_k, daj_k, dai_l, daj_l
   real(wp) :: d2ai_kl, d2aj_kl

   d2K_elem = 0.0_wp

   if (alpha < 1 .or. alpha > 3) return
   if (beta  < 1 .or. beta  > 3) return

   ! Cache the Born-radius derivative scalars (easy to swap later)
   if (present(brdr) .and. present(brdr2)) then
      dai_k   = brdr(alpha, k, i)
      daj_k   = brdr(alpha, k, j)
      dai_l   = brdr(beta , l, i)
      daj_l   = brdr(beta , l, j)
      d2ai_kl = brdr2(alpha, k, beta, l, i)
      d2aj_kl = brdr2(alpha, k, beta, l, j)
   else
      dai_k   = 0.0_wp
      daj_k   = 0.0_wp
      dai_l   = 0.0_wp
      daj_l   = 0.0_wp
      d2ai_kl = 0.0_wp
      d2aj_kl = 0.0_wp
   end if

   ! -------------------------
   ! Diagonal self term: K_ii = keps / a_i
   ! -------------------------
   if (i == j) then
      ai = brad(i)
      if (abs(ai) <= tiny) return

      coef1 = -self%keps / (ai*ai)
      coef2 =  2.0_wp * self%keps / (ai*ai*ai)

      ! d2K = f'' * (da_i)_k (da_i)_l + f' * d2a_i_kl
      d2K_elem = coef2 * (dai_k * dai_l) + coef1 * d2ai_kl
      return
   end if

   ! -------------------------
   ! Off-diagonal element (i,j)
   ! -------------------------
   ai = brad(i)
   aj = brad(j)
   if (abs(ai) <= tiny .or. abs(aj) <= tiny) return

   v(:) = xyz(:, i) - xyz(:, j)
   r2   = dot_product(v, v)

   A = ai * aj
   if (abs(A) <= tiny) return

   d = a4 * r2 / A
   E = exp(-d)
   P = 1.0_wp - a4 * E

   S     = r2 + A * E
   if (S <= tiny) return
   invf  = 1.0_wp / sqrt(S)
   invf3 = invf*invf*invf
   invf5 = invf3*invf*invf

   Qi = aj * E * (1.0_wp + d)
   Qj = ai * E * (1.0_wp + d)

   dK_dai = -0.5_wp * self%keps * Qi * invf3
   dK_daj = -0.5_wp * self%keps * Qj * invf3

   C = - self%keps * P * invf3

   dCdr2 = -self%keps * ( (a4*a4 * E / A) * invf3 - 1.5_wp * (P*P) * invf5 )

   Iab   = 0.0_wp
   if (alpha == beta) Iab = 1.0_wp
   vv_ab = v(alpha) * v(beta)

   Hvv_ab = C * Iab + (2.0_wp * dCdr2) * vv_ab

   dC_dai = self%keps * (a4 * E * d / ai) * invf3 + 1.5_wp * self%keps * P * Qi * invf5
   dC_daj = self%keps * (a4 * E * d / aj) * invf3 + 1.5_wp * self%keps * P * Qj * invf5

   d2K_dai2     = -0.5_wp * self%keps * (aj * E * d*d / ai) * invf3 + 0.75_wp * self%keps * (Qi*Qi) * invf5
   d2K_daj2     = -0.5_wp * self%keps * (ai * E * d*d / aj) * invf3 + 0.75_wp * self%keps * (Qj*Qj) * invf5
   d2K_daida_j  = -0.5_wp * self%keps * (E * (1.0_wp + d + d*d)) * invf3 + 0.75_wp * self%keps * (Qi*Qj) * invf5

   delk = 0; if (k == i) delk = delk + 1; if (k == j) delk = delk - 1
   dell = 0; if (l == i) dell = dell + 1; if (l == j) dell = dell - 1

   ! Now assemble the scalar entry that your original M(alpha,beta) contributed:

   ! (1) explicit coordinate Hessian block (only if k,l are i/j)
   if (delk /= 0 .and. dell /= 0) then
      d2K_elem = d2K_elem + real(delk*dell, wp) * Hvv_ab
   end if

   ! (2) mixed explicit/Radius chain pieces (rank-1 outer products)
   if (delk /= 0) then
      d2K_elem = d2K_elem + real(delk, wp) * ( dC_dai * v(alpha) * dai_l + dC_daj * v(alpha) * daj_l )
   end if

   if (dell /= 0) then
      d2K_elem = d2K_elem + real(dell, wp) * ( dC_dai * dai_k * v(beta) + dC_daj * daj_k * v(beta) )
   end if

   ! (3) pure radius chain terms (outer products)
   d2K_elem = d2K_elem + d2K_dai2    * (dai_k * dai_l)
   d2K_elem = d2K_elem + d2K_daj2    * (daj_k * daj_l)
   d2K_elem = d2K_elem + d2K_daida_j * (dai_k * daj_l + daj_k * dai_l)

   ! (4) radius Hessian contributions
   d2K_elem = d2K_elem + dK_dai * d2ai_kl + dK_daj * d2aj_kl

end subroutine compute_still_d2Kdr2

subroutine compute_still_d3Kdr3_ij(self, nat, xyz, brad, brdr, brdr2, brdr3, i, j, d3Kdr3_ij)
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

end subroutine compute_still_d3Kdr3_ij

subroutine compute_still_d3Kdr3(self, nat, xyz, brad, i, j, k, alpha, l, beta, m, gamma, d3K_elem, &
      brdr, brdr2, brdr3)
   ! Element-wise 3rd derivative for Still kernel:
   !   d3K_elem = d^3 K_ij / ( d r_{k,alpha} d r_{l,beta} d r_{m,gamma} )
   !
   ! Matches your *_full semantics by always evaluating in ordered (ip>jp) pair
   ! for the off-diagonal case.
   !
   ! Uses the exact same algebra as your tensor routine, but:
   ! - evaluates only the requested (alpha,beta,gamma; k,l,m) entry
   ! - caches ONLY the needed brdr/brdr2/brdr3 entries as scalars at the top of the code
   !   (so later swapping to element-wise radii derivative providers is localized).

   class(still_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                               ! (3,nat)
   real(wp), intent(in) :: brad(:)                                 ! (nat)
   integer, intent(in) :: i, j
   integer, intent(in) :: k, alpha
   integer, intent(in) :: l, beta
   integer, intent(in) :: m, gamma
   real(wp), intent(out) :: d3K_elem
   real(wp), contiguous, intent(in), optional :: brdr(:, :, :)               ! (3,nat,nat)
   real(wp), contiguous, intent(in), optional :: brdr2(:, :, :, :, :)        ! (3,nat,3,nat,nat)
   real(wp), contiguous, intent(in), optional :: brdr3(:, :, :, :, :, :, :)  ! (3,nat,3,nat,3,nat,nat)

   integer :: ip, jp
   integer :: delk, dell, delm
   real(wp), parameter :: a4 = 0.25_wp
   real(wp), parameter :: tiny = 1.0e-30_wp

   real(wp) :: v(3), r2
   real(wp) :: ai, aj, A, d, E, P, S
   real(wp) :: invf, invf3, invf5, invf7

   ! --- S-derivatives wrt (r2, ai, aj) ---
   real(wp) :: Sr, Srr, Srrr
   real(wp) :: Sa, Sb
   real(wp) :: Sra, Srb
   real(wp) :: Srra, Srrb
   real(wp) :: Saa2, Sbb2, Sab
   real(wp) :: Sraa, Srbb, Srab
   real(wp) :: Saa3, Sbb3, Saa2b, Sa2bb

   ! --- F(S)=keps*S^{-1/2} derivatives ---
   real(wp) :: F1, F2, F3

   ! --- K scalar partials needed ---
   real(wp) :: Kr
   real(wp) :: Kai, Kaj
   real(wp) :: Kaiai, Kajaj, Kaiaj
   real(wp) :: Kr_ai, Kr_aj
   real(wp) :: Krr, Krrr
   real(wp) :: Krr_ai, Krr_aj
   real(wp) :: Kr_aiai, Kr_ajaj, Kr_aiaj
   real(wp) :: Kaiaiai, Kajajaj, Kaiaiaj, Kaiajaj

   ! building blocks at the requested indices only
   real(wp) :: Iab, Iag, Ibg
   real(wp) :: Kvvv_abg
   real(wp) :: Kvv_ai_ab, Kvv_aj_ab
   real(wp) :: Kvv_ai_ag, Kvv_aj_ag
   real(wp) :: Kvv_ai_bg, Kvv_aj_bg
   real(wp) :: Kv_aiai_a, Kv_ajaj_a, Kv_aiaj_a
   real(wp) :: Kv_aiai_b, Kv_ajaj_b, Kv_aiaj_b
   real(wp) :: Kv_aiai_g, Kv_ajaj_g, Kv_aiaj_g
   real(wp) :: Kv_ai_a, Kv_aj_a
   real(wp) :: Kv_ai_b, Kv_aj_b
   real(wp) :: Kv_ai_g, Kv_aj_g

   ! --- Cache ONLY needed Born-radius derivative entries (scalars) ---
   real(wp) :: dai_k, daj_k, dai_l, daj_l, dai_m, daj_m
   real(wp) :: d2ai_kl, d2ai_km, d2ai_lm
   real(wp) :: d2aj_kl, d2aj_km, d2aj_lm
   real(wp) :: d3ai_klm, d3aj_klm

   real(wp) :: term

   d3K_elem = 0.0_wp
   if (alpha < 1 .or. alpha > 3) return
   if (beta  < 1 .or. beta  > 3) return
   if (gamma < 1 .or. gamma > 3) return

   ! -------------------------
   ! Diagonal: K_ii = keps/a_i
   ! -------------------------
   if (i == j) then
      ai = brad(i)
      if (abs(ai) <= tiny) return

      ! f'(a)= -keps/a^2, f''(a)= 2 keps/a^3, f'''(a)= -6 keps/a^4
      ! Element-wise chain rule:
      ! d3K = f''' a1 a1 a1 + f'' (a2 a1 + a2 a1 + a2 a1) + f' a3

      ! cache only what's needed for this (k,l,m; alpha,beta,gamma)
      if (present(brdr) .and. present(brdr2) .and. present(brdr3)) then
         dai_k   = brdr(alpha,k,i)
         dai_l   = brdr(beta ,l,i)
         dai_m   = brdr(gamma,m,i)

         d2ai_kl = brdr2(alpha,k,beta ,l,i)
         d2ai_km = brdr2(alpha,k,gamma,m,i)
         d2ai_lm = brdr2(beta ,l,gamma,m,i)

         d3ai_klm = brdr3(alpha,k,beta,l,gamma,m,i)
      else
         dai_k   = 0.0_wp
         dai_l   = 0.0_wp
         dai_m   = 0.0_wp

         d2ai_kl = 0.0_wp
         d2ai_km = 0.0_wp
         d2ai_lm = 0.0_wp

         d3ai_klm = 0.0_wp
      end if

      term = 0.0_wp
      term = term + (-6.0_wp*self%keps/(ai**4)) * (dai_k*dai_l*dai_m)
      term = term + ( 2.0_wp*self%keps/(ai**3)) * ( d2ai_kl*dai_m + d2ai_km*dai_l + d2ai_lm*dai_k )
      term = term + (-self%keps/(ai*ai))        * d3ai_klm

      d3K_elem = term
      return
   end if

   ! -------------------------
   ! Off-diagonal: match *_full semantics via ordered (ip>jp)
   ! -------------------------
   if (i > j) then
      ip = i; jp = j
   else
      ip = j; jp = i
   end if

   ai = brad(ip)
   aj = brad(jp)
   if (abs(ai) <= tiny .or. abs(aj) <= tiny) return

   v(:) = xyz(:, ip) - xyz(:, jp)
   r2   = dot_product(v, v)

   A = ai * aj
   if (abs(A) <= tiny) return

   d = a4 * r2 / A
   E = exp(-d)

   S = r2 + A * E
   if (S <= tiny) return
   invf  = 1.0_wp / sqrt(S)
   invf3 = invf*invf*invf
   invf5 = invf3*invf*invf
   invf7 = invf5*invf*invf

   F1 = -0.5_wp   * self%keps * invf3
   F2 =  0.75_wp  * self%keps * invf5
   F3 = -1.875_wp * self%keps * invf7

   P   = 1.0_wp - a4*E
   Sr  = P
   Srr  = (a4*a4 / A) * E
   Srrr = -(a4*a4*a4 / (A*A)) * E

   Sa = aj * E * (1.0_wp + d)
   Sb = ai * E * (1.0_wp + d)

   Saa2 = aj * E * (d*d) / ai
   Sbb2 = ai * E * (d*d) / aj
   Sab  = E * (1.0_wp + d + d*d)

   Sra  = -(a4/ai) * E * d
   Srb  = -(a4/aj) * E * d

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

   Krr    = F2*Sr*Sr + F1*Srr
   Kaiai  = F2*Sa*Sa + F1*Saa2
   Kajaj  = F2*Sb*Sb + F1*Sbb2
   Kaiaj  = F2*Sa*Sb + F1*Sab

   Kr_ai  = F2*Sr*Sa + F1*Sra
   Kr_aj  = F2*Sr*Sb + F1*Srb

   Krrr   = F3*Sr*Sr*Sr + 3.0_wp*F2*Srr*Sr + F1*Srrr

   Krr_ai = F3*Sr*Sr*Sa + F2*(Srr*Sa + 2.0_wp*Sr*Sra) + F1*Srra
   Krr_aj = F3*Sr*Sr*Sb + F2*(Srr*Sb + 2.0_wp*Sr*Srb) + F1*Srrb

   Kr_aiai = F3*Sr*Sa*Sa + F2*(Saa2*Sr + 2.0_wp*Sra*Sa) + F1*Sraa
   Kr_ajaj = F3*Sr*Sb*Sb + F2*(Sbb2*Sr + 2.0_wp*Srb*Sb) + F1*Srbb
   Kr_aiaj = F3*Sr*Sa*Sb + F2*(Sab*Sr + Sra*Sb + Srb*Sa) + F1*Srab

   Kaiaiai = F3*Sa*Sa*Sa + 3.0_wp*F2*Saa2*Sa + F1*Saa3
   Kajajaj = F3*Sb*Sb*Sb + 3.0_wp*F2*Sbb2*Sb + F1*Sbb3
   Kaiaiaj = F3*Sa*Sa*Sb + F2*(Saa2*Sb + 2.0_wp*Sab*Sa) + F1*Saa2b
   Kaiajaj = F3*Sa*Sb*Sb + F2*(Sbb2*Sa + 2.0_wp*Sab*Sb) + F1*Sa2bb

   ! del-factors (coordinate dependence through v)
   delk = 0; if (k==ip) delk=delk+1; if (k==jp) delk=delk-1
   dell = 0; if (l==ip) dell=dell+1; if (l==jp) dell=dell-1
   delm = 0; if (m==ip) delm=delm+1; if (m==jp) delm=delm-1

   ! Identity deltas at requested indices
   Iab = 0.0_wp; if (alpha == beta ) Iab = 1.0_wp
   Iag = 0.0_wp; if (alpha == gamma) Iag = 1.0_wp
   Ibg = 0.0_wp; if (beta  == gamma) Ibg = 1.0_wp

   ! Kvvv(alpha,beta,gamma)
   Kvvv_abg = 8.0_wp*Krrr * v(alpha)*v(beta)*v(gamma) &
            + 4.0_wp*Krr  * ( Iab*v(gamma) + Iag*v(beta) + Ibg*v(alpha) )

   ! Kvv_ai(alpha,beta) etc. (only the needed entries)
   Kvv_ai_ab = 2.0_wp*Kr_ai*Iab + 4.0_wp*Krr_ai * v(alpha)*v(beta)
   Kvv_aj_ab = 2.0_wp*Kr_aj*Iab + 4.0_wp*Krr_aj * v(alpha)*v(beta)

   Kvv_ai_ag = 2.0_wp*Kr_ai*Iag + 4.0_wp*Krr_ai * v(alpha)*v(gamma)
   Kvv_aj_ag = 2.0_wp*Kr_aj*Iag + 4.0_wp*Krr_aj * v(alpha)*v(gamma)

   Kvv_ai_bg = 2.0_wp*Kr_ai*Ibg + 4.0_wp*Krr_ai * v(beta)*v(gamma)
   Kvv_aj_bg = 2.0_wp*Kr_aj*Ibg + 4.0_wp*Krr_aj * v(beta)*v(gamma)

   ! Kv_ai(alpha) etc.
   Kv_ai_a = 2.0_wp*Kr_ai * v(alpha)
   Kv_aj_a = 2.0_wp*Kr_aj * v(alpha)

   Kv_ai_b = 2.0_wp*Kr_ai * v(beta)
   Kv_aj_b = 2.0_wp*Kr_aj * v(beta)

   Kv_ai_g = 2.0_wp*Kr_ai * v(gamma)
   Kv_aj_g = 2.0_wp*Kr_aj * v(gamma)

   ! Kv_aiai(alpha) etc.
   Kv_aiai_a = 2.0_wp*Kr_aiai * v(alpha)
   Kv_ajaj_a = 2.0_wp*Kr_ajaj * v(alpha)
   Kv_aiaj_a = 2.0_wp*Kr_aiaj * v(alpha)

   Kv_aiai_b = 2.0_wp*Kr_aiai * v(beta)
   Kv_ajaj_b = 2.0_wp*Kr_ajaj * v(beta)
   Kv_aiaj_b = 2.0_wp*Kr_aiaj * v(beta)

   Kv_aiai_g = 2.0_wp*Kr_aiai * v(gamma)
   Kv_ajaj_g = 2.0_wp*Kr_ajaj * v(gamma)
   Kv_aiaj_g = 2.0_wp*Kr_aiaj * v(gamma)

   ! Cache only required Born-radius derivative entries for ip/jp at (k,l,m; alpha,beta,gamma)
   if (present(brdr) .and. present(brdr2) .and. present(brdr3)) then
      dai_k = brdr(alpha,k,ip);  daj_k = brdr(alpha,k,jp)
      dai_l = brdr(beta ,l,ip);  daj_l = brdr(beta ,l,jp)
      dai_m = brdr(gamma,m,ip);  daj_m = brdr(gamma,m,jp)

      d2ai_kl = brdr2(alpha,k,beta ,l,ip)
      d2ai_km = brdr2(alpha,k,gamma,m,ip)
      d2ai_lm = brdr2(beta ,l,gamma,m,ip)

      d2aj_kl = brdr2(alpha,k,beta ,l,jp)
      d2aj_km = brdr2(alpha,k,gamma,m,jp)
      d2aj_lm = brdr2(beta ,l,gamma,m,jp)

      d3ai_klm = brdr3(alpha,k,beta,l,gamma,m,ip)
      d3aj_klm = brdr3(alpha,k,beta,l,gamma,m,jp)
   else
      dai_k = 0.0_wp;  daj_k = 0.0_wp
      dai_l = 0.0_wp;  daj_l = 0.0_wp
      dai_m = 0.0_wp;  daj_m = 0.0_wp

      d2ai_kl = 0.0_wp
      d2ai_km = 0.0_wp
      d2ai_lm = 0.0_wp

      d2aj_kl = 0.0_wp
      d2aj_km = 0.0_wp
      d2aj_lm = 0.0_wp

      d3ai_klm = 0.0_wp
      d3aj_klm = 0.0_wp
   end if

   term = 0.0_wp

   ! --- exact same assembly as your tensor code, but scalarized at the requested indices ---

   if (delk/=0 .and. dell/=0 .and. delm/=0) then
      term = term + real(delk*dell*delm,wp) * Kvvv_abg
   end if

   if (delk/=0 .and. dell/=0) then
      term = term + real(delk*dell,wp) * Kvv_ai_ab * dai_m
      term = term + real(delk*dell,wp) * Kvv_aj_ab * daj_m
   end if
   if (delk/=0 .and. delm/=0) then
      term = term + real(delk*delm,wp) * Kvv_ai_ag * dai_l
      term = term + real(delk*delm,wp) * Kvv_aj_ag * daj_l
   end if
   if (dell/=0 .and. delm/=0) then
      term = term + real(dell*delm,wp) * Kvv_ai_bg * dai_k
      term = term + real(dell*delm,wp) * Kvv_aj_bg * daj_k
   end if

   if (delk/=0) then
      term = term + real(delk,wp) * Kv_aiai_a * (dai_l*dai_m)
      term = term + real(delk,wp) * Kv_ajaj_a * (daj_l*daj_m)
      term = term + real(delk,wp) * Kv_aiaj_a * (dai_l*daj_m + daj_l*dai_m)
   end if
   if (dell/=0) then
      term = term + real(dell,wp) * Kv_aiai_b * (dai_k*dai_m)
      term = term + real(dell,wp) * Kv_ajaj_b * (daj_k*daj_m)
      term = term + real(dell,wp) * Kv_aiaj_b * (dai_k*daj_m + daj_k*dai_m)
   end if
   if (delm/=0) then
      term = term + real(delm,wp) * Kv_aiai_g * (dai_k*dai_l)
      term = term + real(delm,wp) * Kv_ajaj_g * (daj_k*daj_l)
      term = term + real(delm,wp) * Kv_aiaj_g * (dai_k*daj_l + daj_k*dai_l)
   end if

   term = term + Kaiaiai * (dai_k*dai_l*dai_m)
   term = term + Kajajaj * (daj_k*daj_l*daj_m)

   term = term + Kaiaiaj * (dai_k*dai_l*daj_m + dai_k*daj_l*dai_m + daj_k*dai_l*dai_m)
   term = term + Kaiajaj * (daj_k*daj_l*dai_m + daj_k*dai_l*daj_m + dai_k*daj_l*daj_m)

   if (delm/=0) then
      term = term + real(delm,wp) * d2ai_kl * Kv_ai_g
      term = term + real(delm,wp) * d2aj_kl * Kv_aj_g
   end if
   term = term + d2ai_kl * ( Kaiai * dai_m + Kaiaj * daj_m )
   term = term + d2aj_kl * ( Kaiaj * dai_m + Kajaj * daj_m )

   if (dell/=0) then
      term = term + real(dell,wp) * d2ai_km * Kv_ai_b
      term = term + real(dell,wp) * d2aj_km * Kv_aj_b
   end if
   term = term + d2ai_km * ( Kaiai * dai_l + Kaiaj * daj_l )
   term = term + d2aj_km * ( Kaiaj * dai_l + Kajaj * daj_l )

   if (delk/=0) then
      term = term + real(delk,wp) * d2ai_lm * Kv_ai_a
      term = term + real(delk,wp) * d2aj_lm * Kv_aj_a
   end if
   term = term + d2ai_lm * ( Kaiai * dai_k + Kaiaj * daj_k )
   term = term + d2aj_lm * ( Kaiaj * dai_k + Kajaj * daj_k )

   term = term + Kai * d3ai_klm + Kaj * d3aj_klm

   d3K_elem = term
end subroutine compute_still_d3Kdr3

subroutine compute_still_d4Kdr4(self, nat, xyz, brad, i, j, k, alpha, l, beta, m, gamma, n, delta, d4K_elem, &
      & brdr, brdr2, brdr3, brdr4)
   ! Element-wise 4th derivative for the Still kernel:
   !   d4K_elem = d^4 K_ij / ( d r_{k,alpha} d r_{l,beta} d r_{m,gamma} d r_{n,delta} )
   !
   ! Same algebra as compute_still_d4Kdr4_ij, but returns only one entry.
   ! Also caches the needed Born-radii derivative tensor entries into scalars
   ! at the beginning of the relevant index nesting, so later you can swap those
   ! assignments to element-wise brdr/brdr2/brdr3/brdr4 providers with minimal edits.

   class(still_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                                      ! (3,nat)
   real(wp), intent(in) :: brad(:)                                        ! (nat)
   integer, intent(in) :: i, j
   integer, intent(in) :: k, alpha, l, beta, m, gamma, n, delta
   real(wp), intent(out) :: d4K_elem
   real(wp), contiguous, intent(in), optional :: brdr(:, :, :)                      ! (3,nat,nat)
   real(wp), contiguous, intent(in), optional :: brdr2(:, :, :, :, :)               ! (3,nat,3,nat,nat)
   real(wp), contiguous, intent(in), optional :: brdr3(:, :, :, :, :, :, :)         ! (3,nat,3,nat,3,nat,nat)
   real(wp), contiguous, intent(in), optional :: brdr4(:, :, :, :, :, :, :, :, :)   ! (3,nat,3,nat,3,nat,3,nat,nat)

   integer :: ip, jp
   integer :: delk, dell, delm, deln
   real(wp), parameter :: a4 = 0.25_wp
   real(wp), parameter :: tiny = 1.0e-30_wp

   real(wp) :: v(3), r2
   real(wp) :: ai, aj, A0, B0, w0, E0, S0
   real(wp) :: invf, invf3, invf5, invf7, invf9
   real(wp) :: C1, C2, C3, C4
   real(wp) :: I3(3,3)

   ! ---- diagonal self-term constants ----
   real(wp) :: f1_self, f2_self, f3_self, f4_self

   ! ---- cached Born radii derivatives (ONLY the needed entries) ----
   real(wp) :: dai_k, dai_l, dai_m, dai_n
   real(wp) :: daj_k, daj_l, daj_m, daj_n

   real(wp) :: d2ai_kl, d2ai_km, d2ai_kn, d2ai_lm, d2ai_ln, d2ai_mn
   real(wp) :: d2aj_kl, d2aj_km, d2aj_kn, d2aj_lm, d2aj_ln, d2aj_mn

   real(wp) :: d3ai_klm, d3ai_kln, d3ai_kmn, d3ai_lmn
   real(wp) :: d3aj_klm, d3aj_kln, d3aj_kmn, d3aj_lmn

   real(wp) :: d4ai_klmn, d4aj_klmn

   ! ---- r2 derivatives (singles/pairs only; triples/quads are zero) ----
   real(wp) :: r_k, r_l, r_m, r_n
   real(wp) :: r_kl, r_km, r_kn, r_lm, r_ln, r_mn

   ! ---- A = ai*aj derivatives ----
   real(wp) :: A_k, A_l, A_m, A_n
   real(wp) :: A_kl, A_km, A_kn, A_lm, A_ln, A_mn
   real(wp) :: A_klm, A_kln, A_kmn, A_lmn
   real(wp) :: A_klmn

   ! ---- B = 1/A derivatives ----
   real(wp) :: B_k, B_l, B_m, B_n
   real(wp) :: B_kl, B_km, B_kn, B_lm, B_ln, B_mn
   real(wp) :: B_klm, B_kln, B_kmn, B_lmn
   real(wp) :: B_klmn

   ! ---- w = -a4*r2*B derivatives ----
   real(wp) :: w_k, w_l, w_m, w_n
   real(wp) :: w_kl, w_km, w_kn, w_lm, w_ln, w_mn
   real(wp) :: w_klm, w_kln, w_kmn, w_lmn
   real(wp) :: w_klmn

   ! ---- E = exp(w) derivatives ----
   real(wp) :: E_k, E_l, E_m, E_n
   real(wp) :: E_kl, E_km, E_kn, E_lm, E_ln, E_mn
   real(wp) :: E_klm, E_kln, E_kmn, E_lmn
   real(wp) :: E_klmn

   ! ---- S = r2 + A*E derivatives (singles..quad) ----
   real(wp) :: S_k, S_l, S_m, S_n
   real(wp) :: S_kl, S_km, S_kn, S_lm, S_ln, S_mn
   real(wp) :: S_klm, S_kln, S_kmn, S_lmn
   real(wp) :: S_klmn

   real(wp) :: term, sum1, sum2, sum3

   d4K_elem = 0.0_wp

   ! Safety: components must be 1..3
   if (alpha < 1 .or. alpha > 3) return
   if (beta  < 1 .or. beta  > 3) return
   if (gamma < 1 .or. gamma > 3) return
   if (delta < 1 .or. delta > 3) return

   ! Identity
   I3 = 0.0_wp
   I3(1,1)=1.0_wp; I3(2,2)=1.0_wp; I3(3,3)=1.0_wp

   ! -------------------------
   ! Diagonal self term: K_ii = keps / a_i
   ! -------------------------
   if (i == j) then
      ai = brad(i)
      if (abs(ai) <= tiny) return

      f1_self = -self%keps / (ai*ai)
      f2_self =  2.0_wp * self%keps / (ai**3)
      f3_self = -6.0_wp * self%keps / (ai**4)
      f4_self =  24.0_wp * self%keps / (ai**5)

      ! Cache only the needed radii derivatives for this entry
      if (present(brdr) .and. present(brdr2) .and. present(brdr3) .and. present(brdr4)) then
         dai_k     = brdr(alpha, k, i)
         dai_l     = brdr(beta , l, i)
         dai_m     = brdr(gamma, m, i)
         dai_n     = brdr(delta, n, i)

         d2ai_kl   = brdr2(alpha, k, beta , l, i)
         d2ai_km   = brdr2(alpha, k, gamma, m, i)
         d2ai_kn   = brdr2(alpha, k, delta, n, i)
         d2ai_lm   = brdr2(beta , l, gamma, m, i)
         d2ai_ln   = brdr2(beta , l, delta, n, i)
         d2ai_mn   = brdr2(gamma, m, delta, n, i)

         d3ai_klm  = brdr3(alpha, k, beta , l, gamma, m, i)
         d3ai_kln  = brdr3(alpha, k, beta , l, delta, n, i)
         d3ai_kmn  = brdr3(alpha, k, gamma, m, delta, n, i)
         d3ai_lmn  = brdr3(beta , l, gamma, m, delta, n, i)

         d4ai_klmn = brdr4(alpha, k, beta, l, gamma, m, delta, n, i)
      else
         dai_k     = 0.0_wp
         dai_l     = 0.0_wp
         dai_m     = 0.0_wp
         dai_n     = 0.0_wp

         d2ai_kl   = 0.0_wp
         d2ai_km   = 0.0_wp
         d2ai_kn   = 0.0_wp
         d2ai_lm   = 0.0_wp
         d2ai_ln   = 0.0_wp
         d2ai_mn   = 0.0_wp

         d3ai_klm  = 0.0_wp
         d3ai_kln  = 0.0_wp
         d3ai_kmn  = 0.0_wp
         d3ai_lmn  = 0.0_wp

         d4ai_klmn = 0.0_wp
      end if

      term = 0.0_wp

      ! f'''' * (a1 a1 a1 a1)
      term = term + f4_self * (dai_k*dai_l*dai_m*dai_n)

      ! f''' * sum_{(2,1,1)} a2 * a1 * a1  (6 terms)
      term = term + f3_self * ( d2ai_kl*dai_m*dai_n + d2ai_km*dai_l*dai_n + d2ai_kn*dai_l*dai_m &
                              + d2ai_lm*dai_k*dai_n + d2ai_ln*dai_k*dai_m + d2ai_mn*dai_k*dai_l )

      ! f'' * [ sum_{(3,1)} a3*a1 (4 terms) + sum_{(2,2)} a2*a2 (3 pairings) ]
      term = term + f2_self * ( d3ai_klm*dai_n + d3ai_kln*dai_m + d3ai_kmn*dai_l + d3ai_lmn*dai_k &
                              + d2ai_kl*d2ai_mn + d2ai_km*d2ai_ln + d2ai_kn*d2ai_lm )

      ! f' * a4
      term = term + f1_self * d4ai_klmn

      d4K_elem = term
      return
   end if

   ! -------------------------
   ! Off-diagonal semantics (match *_full):
   ! compute always in the ip>jp ordering, return that value for (i,j).
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
   if (abs(ai) <= tiny .or. abs(aj) <= tiny) return

   v(:) = xyz(:, ip) - xyz(:, jp)
   r2   = dot_product(v, v)

   A0 = ai * aj
   B0 = 1.0_wp / A0

   w0 = -a4 * r2 * B0
   E0 = exp(w0)

   S0 = r2 + A0 * E0
   if (S0 <= tiny) return

   invf  = 1.0_wp / sqrt(S0)
   invf3 = invf*invf*invf
   invf5 = invf3*invf*invf
   invf7 = invf5*invf*invf
   invf9 = invf7*invf*invf

   ! Coefficients for K = keps * S^{-1/2}
   C1 = -0.5_wp    * self%keps * invf3
   C2 =  0.75_wp   * self%keps * invf5
   C3 = -1.875_wp  * self%keps * invf7     ! -15/8
   C4 =  6.5625_wp * self%keps * invf9     ! 105/16

   ! Coordinate del-factors (for v = r_ip - r_jp)
   delk = 0; if (k==ip) delk=delk+1; if (k==jp) delk=delk-1
   dell = 0; if (l==ip) dell=dell+1; if (l==jp) dell=dell-1
   delm = 0; if (m==ip) delm=delm+1; if (m==jp) delm=delm-1
   deln = 0; if (n==ip) deln=deln+1; if (n==jp) deln=deln-1

   ! Cache ONLY the needed radii derivatives for this entry (ip and jp)
   if (present(brdr) .and. present(brdr2) .and. present(brdr3) .and. present(brdr4)) then
      dai_k     = brdr(alpha, k, ip);   daj_k     = brdr(alpha, k, jp)
      dai_l     = brdr(beta , l, ip);   daj_l     = brdr(beta , l, jp)
      dai_m     = brdr(gamma, m, ip);   daj_m     = brdr(gamma, m, jp)
      dai_n     = brdr(delta, n, ip);   daj_n     = brdr(delta, n, jp)

      d2ai_kl   = brdr2(alpha, k, beta , l, ip);   d2aj_kl   = brdr2(alpha, k, beta , l, jp)
      d2ai_km   = brdr2(alpha, k, gamma, m, ip);   d2aj_km   = brdr2(alpha, k, gamma, m, jp)
      d2ai_kn   = brdr2(alpha, k, delta, n, ip);   d2aj_kn   = brdr2(alpha, k, delta, n, jp)
      d2ai_lm   = brdr2(beta , l, gamma, m, ip);   d2aj_lm   = brdr2(beta , l, gamma, m, jp)
      d2ai_ln   = brdr2(beta , l, delta, n, ip);   d2aj_ln   = brdr2(beta , l, delta, n, jp)
      d2ai_mn   = brdr2(gamma, m, delta, n, ip);   d2aj_mn   = brdr2(gamma, m, delta, n, jp)

      d3ai_klm  = brdr3(alpha, k, beta , l, gamma, m, ip);   d3aj_klm  = brdr3(alpha, k, beta , l, gamma, m, jp)
      d3ai_kln  = brdr3(alpha, k, beta , l, delta, n, ip);   d3aj_kln  = brdr3(alpha, k, beta , l, delta, n, jp)
      d3ai_kmn  = brdr3(alpha, k, gamma, m, delta, n, ip);   d3aj_kmn  = brdr3(alpha, k, gamma, m, delta, n, jp)
      d3ai_lmn  = brdr3(beta , l, gamma, m, delta, n, ip);   d3aj_lmn  = brdr3(beta , l, gamma, m, delta, n, jp)

      d4ai_klmn = brdr4(alpha, k, beta, l, gamma, m, delta, n, ip)
      d4aj_klmn = brdr4(alpha, k, beta, l, gamma, m, delta, n, jp)

   else
      dai_k     = 0.0_wp;   daj_k     = 0.0_wp
      dai_l     = 0.0_wp;   daj_l     = 0.0_wp
      dai_m     = 0.0_wp;   daj_m     = 0.0_wp
      dai_n     = 0.0_wp;   daj_n     = 0.0_wp

      d2ai_kl   = 0.0_wp;   d2aj_kl   = 0.0_wp
      d2ai_km   = 0.0_wp;   d2aj_km   = 0.0_wp
      d2ai_kn   = 0.0_wp;   d2aj_kn   = 0.0_wp
      d2ai_lm   = 0.0_wp;   d2aj_lm   = 0.0_wp
      d2ai_ln   = 0.0_wp;   d2aj_ln   = 0.0_wp
      d2ai_mn   = 0.0_wp;   d2aj_mn   = 0.0_wp

      d3ai_klm  = 0.0_wp;   d3aj_klm  = 0.0_wp
      d3ai_kln  = 0.0_wp;   d3aj_kln  = 0.0_wp
      d3ai_kmn  = 0.0_wp;   d3aj_kmn  = 0.0_wp
      d3ai_lmn  = 0.0_wp;   d3aj_lmn  = 0.0_wp

      d4ai_klmn = 0.0_wp
      d4aj_klmn = 0.0_wp
   end if


   ! r2 first derivatives:
   r_k = 2.0_wp * real(delk,wp) * v(alpha)
   r_l = 2.0_wp * real(dell,wp) * v(beta)
   r_m = 2.0_wp * real(delm,wp) * v(gamma)
   r_n = 2.0_wp * real(deln,wp) * v(delta)

   ! r2 second derivatives:
   r_kl = 2.0_wp * real(delk*dell,wp) * I3(alpha,beta)
   r_km = 2.0_wp * real(delk*delm,wp) * I3(alpha,gamma)
   r_kn = 2.0_wp * real(delk*deln,wp) * I3(alpha,delta)
   r_lm = 2.0_wp * real(dell*delm,wp) * I3(beta,gamma)
   r_ln = 2.0_wp * real(dell*deln,wp) * I3(beta,delta)
   r_mn = 2.0_wp * real(delm*deln,wp) * I3(gamma,delta)

   ! =========================================================
   ! A = ai*aj derivatives (only those needed for this entry)
   ! =========================================================
   A_k = aj*dai_k + ai*daj_k
   A_l = aj*dai_l + ai*daj_l
   A_m = aj*dai_m + ai*daj_m
   A_n = aj*dai_n + ai*daj_n

   A_kl = aj*d2ai_kl + ai*d2aj_kl + (dai_k*daj_l + dai_l*daj_k)
   A_km = aj*d2ai_km + ai*d2aj_km + (dai_k*daj_m + dai_m*daj_k)
   A_kn = aj*d2ai_kn + ai*d2aj_kn + (dai_k*daj_n + dai_n*daj_k)
   A_lm = aj*d2ai_lm + ai*d2aj_lm + (dai_l*daj_m + dai_m*daj_l)
   A_ln = aj*d2ai_ln + ai*d2aj_ln + (dai_l*daj_n + dai_n*daj_l)
   A_mn = aj*d2ai_mn + ai*d2aj_mn + (dai_m*daj_n + dai_n*daj_m)

   A_klm = aj*d3ai_klm + ai*d3aj_klm &
         + (d2ai_kl*daj_m + d2ai_km*daj_l + d2ai_lm*daj_k) &
         + (dai_k*d2aj_lm + dai_l*d2aj_km + dai_m*d2aj_kl)

   A_kln = aj*d3ai_kln + ai*d3aj_kln &
         + (d2ai_kl*daj_n + d2ai_kn*daj_l + d2ai_ln*daj_k) &
         + (dai_k*d2aj_ln + dai_l*d2aj_kn + dai_n*d2aj_kl)

   A_kmn = aj*d3ai_kmn + ai*d3aj_kmn &
         + (d2ai_km*daj_n + d2ai_kn*daj_m + d2ai_mn*daj_k) &
         + (dai_k*d2aj_mn + dai_m*d2aj_kn + dai_n*d2aj_km)

   A_lmn = aj*d3ai_lmn + ai*d3aj_lmn &
         + (d2ai_lm*daj_n + d2ai_ln*daj_m + d2ai_mn*daj_l) &
         + (dai_l*d2aj_mn + dai_m*d2aj_ln + dai_n*d2aj_lm)

   A_klmn = aj*d4ai_klmn + ai*d4aj_klmn &
          + (d3ai_klm*daj_n + d3ai_kln*daj_m + d3ai_kmn*daj_l + d3ai_lmn*daj_k) &
          + (dai_k*d3aj_lmn + dai_l*d3aj_kmn + dai_m*d3aj_kln + dai_n*d3aj_klm) &
          + (d2ai_kl*d2aj_mn + d2ai_km*d2aj_ln + d2ai_kn*d2aj_lm &
           + d2ai_lm*d2aj_kn + d2ai_ln*d2aj_km + d2ai_mn*d2aj_kl)

   ! =========================================================
   ! B = 1/A derivatives (closed form up to 4th)
   ! =========================================================
   B_k = -(B0*B0) * A_k
   B_l = -(B0*B0) * A_l
   B_m = -(B0*B0) * A_m
   B_n = -(B0*B0) * A_n

   B_kl = 2.0_wp*(B0**3)*A_k*A_l - (B0*B0)*A_kl
   B_km = 2.0_wp*(B0**3)*A_k*A_m - (B0*B0)*A_km
   B_kn = 2.0_wp*(B0**3)*A_k*A_n - (B0*B0)*A_kn
   B_lm = 2.0_wp*(B0**3)*A_l*A_m - (B0*B0)*A_lm
   B_ln = 2.0_wp*(B0**3)*A_l*A_n - (B0*B0)*A_ln
   B_mn = 2.0_wp*(B0**3)*A_m*A_n - (B0*B0)*A_mn

   B_klm = -6.0_wp*(B0**4)*A_k*A_l*A_m &
         + 2.0_wp*(B0**3)*(A_kl*A_m + A_km*A_l + A_lm*A_k) &
         - (B0*B0)*A_klm

   B_kln = -6.0_wp*(B0**4)*A_k*A_l*A_n &
         + 2.0_wp*(B0**3)*(A_kl*A_n + A_kn*A_l + A_ln*A_k) &
         - (B0*B0)*A_kln

   B_kmn = -6.0_wp*(B0**4)*A_k*A_m*A_n &
         + 2.0_wp*(B0**3)*(A_km*A_n + A_kn*A_m + A_mn*A_k) &
         - (B0*B0)*A_kmn

   B_lmn = -6.0_wp*(B0**4)*A_l*A_m*A_n &
         + 2.0_wp*(B0**3)*(A_lm*A_n + A_ln*A_m + A_mn*A_l) &
         - (B0*B0)*A_lmn

   B_klmn = 24.0_wp*(B0**5)*A_k*A_l*A_m*A_n &
          - 6.0_wp*(B0**4)*( A_kl*A_m*A_n + A_km*A_l*A_n + A_kn*A_l*A_m &
                           + A_lm*A_k*A_n + A_ln*A_k*A_m + A_mn*A_k*A_l ) &
          + 2.0_wp*(B0**3)*( A_klm*A_n + A_kln*A_m + A_kmn*A_l + A_lmn*A_k &
                           + A_kl*A_mn + A_km*A_ln + A_kn*A_lm ) &
          - (B0*B0)*A_klmn

   ! =========================================================
   ! w = -a4 * r2 * B derivatives (r2 has only 1st/2nd derivatives)
   ! =========================================================
   w_k = -a4 * ( r_k*B0 + r2*B_k )
   w_l = -a4 * ( r_l*B0 + r2*B_l )
   w_m = -a4 * ( r_m*B0 + r2*B_m )
   w_n = -a4 * ( r_n*B0 + r2*B_n )

   w_kl = -a4 * ( r_kl*B0 + r_k*B_l + r_l*B_k + r2*B_kl )
   w_km = -a4 * ( r_km*B0 + r_k*B_m + r_m*B_k + r2*B_km )
   w_kn = -a4 * ( r_kn*B0 + r_k*B_n + r_n*B_k + r2*B_kn )
   w_lm = -a4 * ( r_lm*B0 + r_l*B_m + r_m*B_l + r2*B_lm )
   w_ln = -a4 * ( r_ln*B0 + r_l*B_n + r_n*B_l + r2*B_ln )
   w_mn = -a4 * ( r_mn*B0 + r_m*B_n + r_n*B_m + r2*B_mn )

   w_klm = -a4 * ( r_kl*B_m + r_km*B_l + r_lm*B_k &
                 + r_k*B_lm + r_l*B_km + r_m*B_kl &
                 + r2*B_klm )

   w_kln = -a4 * ( r_kl*B_n + r_kn*B_l + r_ln*B_k &
                 + r_k*B_ln + r_l*B_kn + r_n*B_kl &
                 + r2*B_kln )

   w_kmn = -a4 * ( r_km*B_n + r_kn*B_m + r_mn*B_k &
                 + r_k*B_mn + r_m*B_kn + r_n*B_km &
                 + r2*B_kmn )

   w_lmn = -a4 * ( r_lm*B_n + r_ln*B_m + r_mn*B_l &
                 + r_l*B_mn + r_m*B_ln + r_n*B_lm &
                 + r2*B_lmn )

   w_klmn = -a4 * ( r_kl*B_mn + r_km*B_ln + r_kn*B_lm &
                  + r_lm*B_kn + r_ln*B_km + r_mn*B_kl &
                  + r_k*B_lmn + r_l*B_kmn + r_m*B_kln + r_n*B_klm &
                  + r2*B_klmn )

   ! =========================================================
   ! E = exp(w) derivatives (Bell polynomial up to 4th)
   ! =========================================================
   E_k = E0 * w_k
   E_l = E0 * w_l
   E_m = E0 * w_m
   E_n = E0 * w_n

   E_kl = E0 * ( w_kl + w_k*w_l )
   E_km = E0 * ( w_km + w_k*w_m )
   E_kn = E0 * ( w_kn + w_k*w_n )
   E_lm = E0 * ( w_lm + w_l*w_m )
   E_ln = E0 * ( w_ln + w_l*w_n )
   E_mn = E0 * ( w_mn + w_m*w_n )

   E_klm = E0 * ( w_klm + w_kl*w_m + w_km*w_l + w_lm*w_k + w_k*w_l*w_m )
   E_kln = E0 * ( w_kln + w_kl*w_n + w_kn*w_l + w_ln*w_k + w_k*w_l*w_n )
   E_kmn = E0 * ( w_kmn + w_km*w_n + w_kn*w_m + w_mn*w_k + w_k*w_m*w_n )
   E_lmn = E0 * ( w_lmn + w_lm*w_n + w_ln*w_m + w_mn*w_l + w_l*w_m*w_n )

   E_klmn = E0 * ( w_klmn &
           + (w_klm*w_n + w_kln*w_m + w_kmn*w_l + w_lmn*w_k) &
           + (w_kl*w_mn + w_km*w_ln + w_kn*w_lm) &
           + (w_kl*w_m*w_n + w_km*w_l*w_n + w_kn*w_l*w_m &
            + w_lm*w_k*w_n + w_ln*w_k*w_m + w_mn*w_k*w_l) &
           + (w_k*w_l*w_m*w_n) )

   ! =========================================================
   ! S = r2 + A*E derivatives up to 4th
   ! =========================================================
   S_k = r_k + A_k*E0 + A0*E_k
   S_l = r_l + A_l*E0 + A0*E_l
   S_m = r_m + A_m*E0 + A0*E_m
   S_n = r_n + A_n*E0 + A0*E_n

   S_kl = r_kl + A_kl*E0 + A_k*E_l + A_l*E_k + A0*E_kl
   S_km = r_km + A_km*E0 + A_k*E_m + A_m*E_k + A0*E_km
   S_kn = r_kn + A_kn*E0 + A_k*E_n + A_n*E_k + A0*E_kn
   S_lm = r_lm + A_lm*E0 + A_l*E_m + A_m*E_l + A0*E_lm
   S_ln = r_ln + A_ln*E0 + A_l*E_n + A_n*E_l + A0*E_ln
   S_mn = r_mn + A_mn*E0 + A_m*E_n + A_n*E_m + A0*E_mn

   S_klm = A_klm*E0 &
         + (A_kl*E_m + A_km*E_l + A_lm*E_k) &
         + (A_k*E_lm + A_l*E_km + A_m*E_kl) &
         + A0*E_klm

   S_kln = A_kln*E0 &
         + (A_kl*E_n + A_kn*E_l + A_ln*E_k) &
         + (A_k*E_ln + A_l*E_kn + A_n*E_kl) &
         + A0*E_kln

   S_kmn = A_kmn*E0 &
         + (A_km*E_n + A_kn*E_m + A_mn*E_k) &
         + (A_k*E_mn + A_m*E_kn + A_n*E_km) &
         + A0*E_kmn

   S_lmn = A_lmn*E0 &
         + (A_lm*E_n + A_ln*E_m + A_mn*E_l) &
         + (A_l*E_mn + A_m*E_ln + A_n*E_lm) &
         + A0*E_lmn

   S_klmn = A_klmn*E0 &
          + (A_klm*E_n + A_kln*E_m + A_kmn*E_l + A_lmn*E_k) &
          + (A_kl*E_mn + A_km*E_ln + A_kn*E_lm + A_lm*E_kn + A_ln*E_km + A_mn*E_kl) &
          + (A_k*E_lmn + A_l*E_kmn + A_m*E_kln + A_n*E_klm) &
          + A0*E_klmn

   ! =========================================================
   ! 4th derivative of K = keps * S^{-1/2} via partitions
   ! =========================================================
   sum1 = S_klm*S_n + S_kln*S_m + S_kmn*S_l + S_lmn*S_k
   sum2 = S_kl*S_mn + S_km*S_ln + S_kn*S_lm
   sum3 = S_kl*S_m*S_n + S_km*S_l*S_n + S_kn*S_l*S_m &
        + S_lm*S_k*S_n + S_ln*S_k*S_m + S_mn*S_k*S_l

   term = 0.0_wp
   term = term + C1 * S_klmn
   term = term + C2 * (sum1 + sum2)
   term = term + C3 * sum3
   term = term + C4 * (S_k*S_l*S_m*S_n)

   d4K_elem = term
end subroutine compute_still_d4Kdr4



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

subroutine compute_p16_dKdr(self, nat, xyz, brad, i, j, k, alpha, dKdr_elem, brdr)
   class(p16_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                 ! (3,nat)
   real(wp), intent(in) :: brad(:)                   ! (nat)
   integer, intent(in) :: i, j
   integer, intent(in) :: k, alpha
   real(wp), intent(out) :: dKdr_elem
   real(wp), contiguous, intent(in), optional :: brdr(:, :, :) ! (3,nat,nat)

   integer :: ip, jp, delk
   real(wp) :: rvec(3), r, r2
   real(wp) :: ai, aj, ab
   real(wp) :: a1, a16
   real(wp) :: fgb, invfgb2
   real(wp) :: coef_pos
   real(wp) :: bp, dK_dai, dK_daj
   real(wp), parameter :: tiny_r = 1.0e-14_wp
   real(wp), parameter :: tiny = 1.0e-30_wp

   dKdr_elem = 0.0_wp
   if (alpha < 1 .or. alpha > 3) return

   ! Diagonal: K_ii = keps/a_i
   if (i == j) then
      ai = brad(i)
      if (abs(ai) <= tiny) return
      if (present(brdr)) then
         dKdr_elem = (-self%keps / (ai*ai)) * brdr(alpha, k, i)
      else 
         dKdr_elem = (-self%keps / (ai*ai)) * 0.0_wp 
      end if
      return
   end if

   ! Off-diagonal semantics (match your *_full): compute using ip>jp ordering
   if (i > j) then
      ip = i
      jp = j
   else
      ip = j
      jp = i
   end if

   ai = brad(ip)
   aj = brad(jp)
   if (abs(ai) <= tiny .or. abs(aj) <= tiny) return

   rvec(:) = xyz(:, ip) - xyz(:, jp)
   r2      = dot_product(rvec, rvec)
   r       = sqrt(r2)
   if (r <= tiny_r) return

   ab = sqrt(ai * aj)

   a1  = ab / (ab + zetaP16o16 * r)
   a16 = a1*a1
   a16 = a16*a16
   a16 = a16*a16
   a16 = a16*a16

   fgb     = r + ab * a16
   invfgb2 = 1.0_wp / (fgb * fgb)

   ! (A) Explicit coordinate part (radii held fixed): only k==ip or k==jp contributes
   delk = 0
   if (k == ip) delk = delk + 1
   if (k == jp) delk = delk - 1

   if (delk /= 0) then
      coef_pos = -self%keps * (1.0_wp - zetaP16 * a1 * a16) * invfgb2 / r
      dKdr_elem = dKdr_elem + real(delk, wp) * coef_pos * rvec(alpha)
   end if

   ! (B) Born radii partials for chain rule
   bp     = -0.5_wp * ( (r * zetaP16 / ab) * a1 + 1.0_wp ) / ab * a16 * invfgb2
   dK_dai = self%keps * aj * bp
   dK_daj = self%keps * ai * bp

   if (present(brdr)) then
      dKdr_elem = dKdr_elem + dK_dai * brdr(alpha, k, ip) + dK_daj * brdr(alpha, k, jp)
   end if
end subroutine compute_p16_dKdr

subroutine compute_p16_d2Kdr2(self, nat, xyz, brad, i, j, k, alpha, l, beta, d2K_elem, brdr, brdr2)
   class(p16_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                         ! (3,nat)
   real(wp), intent(in) :: brad(:)                           ! (nat)
   integer, intent(in) :: i, j
   integer, intent(in) :: k, alpha, l, beta
   real(wp), intent(out) :: d2K_elem
   real(wp), contiguous, intent(in), optional :: brdr(:, :, :)         ! (3,nat,nat)
   real(wp), contiguous, intent(in), optional :: brdr2(:, :, :, :, :)  ! (3,nat,3,nat,nat)

   integer :: ip, jp
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
   real(wp) :: I3ab, Hvv_ab
   real(wp) :: dk_i_a, dk_j_a, dl_i_b, dl_j_b
   real(wp) :: coef1, coef2
   real(wp), parameter :: tiny_r = 1.0e-14_wp
   real(wp), parameter :: tiny = 1.0e-30_wp

   d2K_elem = 0.0_wp
   if (alpha < 1 .or. alpha > 3) return
   if (beta  < 1 .or. beta  > 3) return

   I3ab = 0.0_wp
   if (alpha == beta) I3ab = 1.0_wp

   cp16 = zetaP16o16

   ! -------------------------
   ! Diagonal self terms: K_ii = keps / a_i
   ! -------------------------
   if (i == j) then
      ai = brad(i)
      if (abs(ai) <= tiny) return

      coef1 = -self%keps / (ai*ai)
      coef2 =  2.0_wp * self%keps / (ai*ai*ai)

      if (present(brdr)) then
         dk_i_a = brdr(alpha, k, i)
         dl_i_b = brdr(beta , l, i)
      else
         dk_i_a = 0.0_wp
         dl_i_b = 0.0_wp
      end if

      if (present(brdr2)) then
         d2K_elem = coef2 * dk_i_a * dl_i_b + coef1 * brdr2(alpha, k, beta, l, i)
      end if
      return
   end if

   ! -------------------------
   ! Off-diagonal: match *_full semantics (ip>jp)
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
   if (abs(ai) <= tiny .or. abs(aj) <= tiny) return

   v(:) = xyz(:, ip) - xyz(:, jp)
   r2   = dot_product(v, v)
   r    = sqrt(r2)
   if (r <= tiny_r) return

   invr  = 1.0_wp / r
   invr2 = invr * invr

   u = sqrt(ai * aj)
   t = u + cp16 * r

   a1  = u / t
   a16 = a1*a1
   a16 = a16*a16
   a16 = a16*a16
   a16 = a16*a16
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

   ! Hessian wrt v (radii held fixed): Hvv(alpha,beta)
   Hvv_ab = C * I3ab + (Krr - C) * (v(alpha) * v(beta) * invr2)

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
   d2K_dai2    = -self%keps * ( g_aiai * invg2 - 2.0_wp * (g_ai*g_ai) * invg3 )
   d2K_daj2    = -self%keps * ( g_ajaj * invg2 - 2.0_wp * (g_aj*g_aj) * invg3 )
   d2K_daida_j = -self%keps * ( g_aiaj * invg2 - 2.0_wp * (g_ai*g_aj) * invg3 )

   ! Need ∂C/∂a_i, ∂C/∂a_j for mixed v–a terms:
   gr_u  = -17.0_wp * zetaP16 * cp16 * r * a16 / (t*t)
   gr_ai = gr_u * u_ai
   gr_aj = gr_u * u_aj

   Kr_ai = -self%keps * ( gr_ai * invg2 - 2.0_wp * gr * g_ai * invg3 )
   Kr_aj = -self%keps * ( gr_aj * invg2 - 2.0_wp * gr * g_aj * invg3 )

   dC_dai = Kr_ai * invr
   dC_daj = Kr_aj * invr

   ! Coordinate deltas (for v = r_ip - r_jp)
   delk = 0
   if (k == ip) delk = delk + 1
   if (k == jp) delk = delk - 1

   dell = 0
   if (l == ip) dell = dell + 1
   if (l == jp) dell = dell - 1

   ! Cache only needed brdr components for this element
   if (present(brdr)) then
      dk_i_a = brdr(alpha, k, ip)
      dk_j_a = brdr(alpha, k, jp)
      dl_i_b = brdr(beta , l, ip)
      dl_j_b = brdr(beta , l, jp)
   else
      dk_i_a = 0.0_wp
      dk_j_a = 0.0_wp
      dl_i_b = 0.0_wp
      dl_j_b = 0.0_wp
   end if

   ! Assemble scalar entry (alpha,k ; beta,l)
   ! (1) Explicit vv part
   if (delk /= 0 .and. dell /= 0) then
      d2K_elem = d2K_elem + real(delk*dell, wp) * Hvv_ab
   end if

   ! (2) Mixed v–a parts
   if (delk /= 0) then
      d2K_elem = d2K_elem + real(delk, wp) * ( dC_dai * v(alpha) * dl_i_b &
                                            + dC_daj * v(alpha) * dl_j_b )
   end if

   if (dell /= 0) then
      d2K_elem = d2K_elem + real(dell, wp) * ( dC_dai * dk_i_a * v(beta) &
                                            + dC_daj * dk_j_a * v(beta) )
   end if

   ! (3) a–a parts
   d2K_elem = d2K_elem + d2K_dai2    * (dk_i_a * dl_i_b)
   d2K_elem = d2K_elem + d2K_daj2    * (dk_j_a * dl_j_b)
   d2K_elem = d2K_elem + d2K_daida_j * ( dk_i_a * dl_j_b + dk_j_a * dl_i_b )

   ! (4) brdr2 terms
   if (present(brdr2)) then
      d2K_elem = d2K_elem + dK_dai * brdr2(alpha, k, beta, l, ip) &
         + dK_daj * brdr2(alpha, k, beta, l, jp)
   end if

end subroutine compute_p16_d2Kdr2

subroutine compute_p16_d3Kdr3_ij(self, nat, xyz, brad, brdr, brdr2, brdr3, i, j, d3Kdr3_ij)
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

end subroutine compute_p16_d3Kdr3_ij

subroutine compute_p16_d3Kdr3(self, nat, xyz, brad, i, j, k, alpha, l, beta, m, gamma, d3K_elem, &
      brdr, brdr2, brdr3)
   ! Element-wise 3rd derivative for P16 kernel:
   !   d3K_elem = d^3 K_ij / ( d r_{k,alpha} d r_{l,beta} d r_{m,gamma} )
   !
   ! Same philosophy as Still:
   ! - compute only one entry, no (3,nat,3,nat,3,nat) storage
   ! - cache ONLY the needed Born-radius derivative entries (brdr, brdr2, brdr3) as scalars
   !   so you can later replace them with element-wise providers easily.

   class(p16_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                               ! (3,nat)
   real(wp), intent(in) :: brad(:)                                 ! (nat)
   integer, intent(in) :: i, j
   integer, intent(in) :: k, alpha
   integer, intent(in) :: l, beta
   integer, intent(in) :: m, gamma
   real(wp), intent(out) :: d3K_elem
   real(wp), contiguous, intent(in), optional :: brdr(:, :, :)               ! (3,nat,nat)
   real(wp), contiguous, intent(in), optional :: brdr2(:, :, :, :, :)        ! (3,nat,3,nat,nat)
   real(wp), contiguous, intent(in), optional :: brdr3(:, :, :, :, :, :, :)  ! (3,nat,3,nat,3,nat,nat)

   integer :: ip, jp
   integer :: delk, dell, delm

   real(wp), parameter :: tiny_r = 1.0e-14_wp
   real(wp), parameter :: tiny = 1.0e-30_wp
   real(wp) :: c

   real(wp) :: v(3), r2, r, invr, invr2
   real(wp) :: ai, aj, u, t, cp16
   real(wp) :: a1, a16, a17
   real(wp) :: g, invg, invg2, invg3, invg4

   ! --- g derivatives wrt r (radii fixed) ---
   real(wp) :: gr, grr, grrr

   ! --- g derivatives wrt u (r fixed) ---
   real(wp) :: gu, guu, guuu

   ! --- mixed g_r,u etc (only need gr_u, gr_uu, grr_u) ---
   real(wp) :: gr_u, gr_uu, grr_u

   ! --- u-derivatives wrt ai,aj ---
   real(wp) :: u_ai, u_aj
   real(wp) :: u_aiai, u_ajaj, u_aiaj
   real(wp) :: u_aiaiai, u_ajajaj
   real(wp) :: u_aiaiaj, u_aiajaj

   ! --- convert g-derivatives to ai/aj derivatives ---
   real(wp) :: g_ai, g_aj
   real(wp) :: g_aiai, g_ajaj, g_aiaj
   real(wp) :: g_aiaiai, g_ajajaj, g_aiaiaj, g_aiajaj

   real(wp) :: gr_ai, gr_aj
   real(wp) :: gr_aiai, gr_ajaj, gr_aiaj
   real(wp) :: grr_ai, grr_aj

   ! --- f(g)=keps/g derivatives ---
   real(wp) :: f1, f2, f3

   ! --- scalar K partials (r, ai, aj) ---
   real(wp) :: Kr, Krr, Krrr
   real(wp) :: Kai, Kaj
   real(wp) :: Kaiai, Kajaj, Kaiaj
   real(wp) :: Kr_ai, Kr_aj
   real(wp) :: Krr_ai, Krr_aj
   real(wp) :: Kr_aiai, Kr_ajaj, Kr_aiaj
   real(wp) :: Kaiaiai, Kajajaj, Kaiaiaj, Kaiajaj

   ! --- radial tensors (v-derivatives holding radii fixed) ---
   real(wp) :: e(3)
   real(wp) :: I3(3,3)
   real(wp) :: A, B, Bp, Ap
   real(wp) :: Kvv_ai_ab, Kvv_aj_ab
   real(wp) :: Kvv_ai_ag, Kvv_aj_ag
   real(wp) :: Kvv_ai_bg, Kvv_aj_bg
   real(wp) :: Kvvv_abg
   real(wp) :: Kv_ai_a, Kv_aj_a
   real(wp) :: Kv_ai_b, Kv_aj_b
   real(wp) :: Kv_ai_g, Kv_aj_g
   real(wp) :: Kv_aiai_a, Kv_ajaj_a, Kv_aiaj_a
   real(wp) :: Kv_aiai_b, Kv_ajaj_b, Kv_aiaj_b
   real(wp) :: Kv_aiai_g, Kv_ajaj_g, Kv_aiaj_g

   ! --- coordinate-to-radius scalars ---
   real(wp) :: dai_k, dai_l, dai_m
   real(wp) :: daj_k, daj_l, daj_m
   real(wp) :: d2ai_kl, d2ai_km, d2ai_lm
   real(wp) :: d2aj_kl, d2aj_km, d2aj_lm
   real(wp) :: d3ai_klm, d3aj_klm

   real(wp) :: term
   real(wp) :: self_f1, self_f2, self_f3

   I3 = 0.0_wp
   I3(1,1)=1.0_wp; I3(2,2)=1.0_wp; I3(3,3)=1.0_wp

   d3K_elem = 0.0_wp

   if (alpha < 1 .or. alpha > 3) return
   if (beta  < 1 .or. beta  > 3) return
   if (gamma < 1 .or. gamma > 3) return

   c = zetaP16o16

   ! ==========================================================
   ! Diagonal self term: K_ii = keps / a_i
   ! ==========================================================
   if (i == j) then
      ai = brad(i)
      if (abs(ai) <= tiny) return

      self_f1 = -self%keps / (ai*ai)
      self_f2 =  2.0_wp * self%keps / (ai*ai*ai)
      self_f3 = -6.0_wp * self%keps / (ai**4)

      ! Cache only the needed radii derivatives for this entry
      if (present(brdr) .and. present(brdr2) .and. present(brdr3)) then
         dai_k    = brdr(alpha,k,i)
         dai_l    = brdr(beta,l,i)
         dai_m    = brdr(gamma,m,i)
         d2ai_kl  = brdr2(alpha,k,beta,l,i)
         d2ai_km  = brdr2(alpha,k,gamma,m,i)
         d2ai_lm  = brdr2(beta,l,gamma,m,i)
         d3ai_klm = brdr3(alpha,k,beta,l,gamma,m,i)
      else
         dai_k    = 0.0_wp
         dai_l    = 0.0_wp
         dai_m    = 0.0_wp
         d2ai_kl  = 0.0_wp
         d2ai_km  = 0.0_wp
         d2ai_lm  = 0.0_wp
         d3ai_klm = 0.0_wp
      end if

      term = 0.0_wp
      term = term + self_f3 * (dai_k*dai_l*dai_m)
      term = term + self_f2 * ( d2ai_kl*dai_m + d2ai_km*dai_l + d2ai_lm*dai_k )
      term = term + self_f1 * d3ai_klm

      d3K_elem = term
      return
   end if

   ! ==========================================================
   ! Off-diagonal: match *_full semantics (computed for ip>jp)
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
   if (abs(ai) <= tiny .or. abs(aj) <= tiny) return

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

   ! Kvv_ai / Kvv_aj at needed indices
   Kvv_ai_ab = (Krr_ai - Kr_ai*invr) * e(alpha)*e(beta) + (Kr_ai*invr) * I3(alpha,beta)
   Kvv_aj_ab = (Krr_aj - Kr_aj*invr) * e(alpha)*e(beta) + (Kr_aj*invr) * I3(alpha,beta)

   Kvv_ai_ag = (Krr_ai - Kr_ai*invr) * e(alpha)*e(gamma) + (Kr_ai*invr) * I3(alpha,gamma)
   Kvv_aj_ag = (Krr_aj - Kr_aj*invr) * e(alpha)*e(gamma) + (Kr_aj*invr) * I3(alpha,gamma)

   Kvv_ai_bg = (Krr_ai - Kr_ai*invr) * e(beta)*e(gamma) + (Kr_ai*invr) * I3(beta,gamma)
   Kvv_aj_bg = (Krr_aj - Kr_aj*invr) * e(beta)*e(gamma) + (Kr_aj*invr) * I3(beta,gamma)

   ! Kv_ai etc at needed indices
   Kv_ai_a = Kr_ai * e(alpha)
   Kv_aj_a = Kr_aj * e(alpha)

   Kv_ai_b = Kr_ai * e(beta)
   Kv_aj_b = Kr_aj * e(beta)

   Kv_ai_g = Kr_ai * e(gamma)
   Kv_aj_g = Kr_aj * e(gamma)

   Kv_aiai_a = Kr_aiai * e(alpha)
   Kv_ajaj_a = Kr_ajaj * e(alpha)
   Kv_aiaj_a = Kr_aiaj * e(alpha)

   Kv_aiai_b = Kr_aiai * e(beta)
   Kv_ajaj_b = Kr_ajaj * e(beta)
   Kv_aiaj_b = Kr_aiaj * e(beta)

   Kv_aiai_g = Kr_aiai * e(gamma)
   Kv_ajaj_g = Kr_ajaj * e(gamma)
   Kv_aiaj_g = Kr_aiaj * e(gamma)

   ! Kvvv at the needed index triple
   Kvvv_abg = Ap * e(alpha)*e(beta)*e(gamma) &
      + (A*invr) * ( I3(alpha,gamma)*e(beta) + I3(beta,gamma)*e(alpha) - 2.0_wp*e(alpha)*e(beta)*e(gamma) ) &
      + Bp * I3(alpha,beta) * e(gamma)

   ! ==========================================================
   ! Coordinate del-factors (for v = r_ip - r_jp)
   ! ==========================================================
   delk = 0; if (k==ip) delk=delk+1; if (k==jp) delk=delk-1
   dell = 0; if (l==ip) dell=dell+1; if (l==jp) dell=dell-1
   delm = 0; if (m==ip) delm=delm+1; if (m==jp) delm=delm-1

   ! Cache ONLY the needed radii derivatives for this entry (ip and jp)
   if (present(brdr) .and. present(brdr2) .and. present(brdr3)) then
      dai_k = brdr(alpha,k,ip); daj_k = brdr(alpha,k,jp)
      dai_l = brdr(beta ,l,ip); daj_l = brdr(beta ,l,jp)
      dai_m = brdr(gamma,m,ip); daj_m = brdr(gamma,m,jp)

      d2ai_kl = brdr2(alpha,k,beta,l,ip)
      d2ai_km = brdr2(alpha,k,gamma,m,ip)
      d2ai_lm = brdr2(beta,l,gamma,m,ip)

      d2aj_kl = brdr2(alpha,k,beta,l,jp)
      d2aj_km = brdr2(alpha,k,gamma,m,jp)
      d2aj_lm = brdr2(beta,l,gamma,m,jp)

      d3ai_klm = brdr3(alpha,k,beta,l,gamma,m,ip)
      d3aj_klm = brdr3(alpha,k,beta,l,gamma,m,jp)
   else
      dai_k = 0.0_wp; daj_k = 0.0_wp
      dai_l = 0.0_wp; daj_l = 0.0_wp
      dai_m = 0.0_wp; daj_m = 0.0_wp

      d2ai_kl = 0.0_wp; d2aj_kl = 0.0_wp
      d2ai_km = 0.0_wp; d2aj_km = 0.0_wp
      d2ai_lm = 0.0_wp; d2aj_lm = 0.0_wp

      d3ai_klm = 0.0_wp
      d3aj_klm = 0.0_wp
   end if

   ! ==========================================================
   ! Assemble the scalar entry (exact same algebra as *_ij but scalarized)
   ! ==========================================================
   term = 0.0_wp

   ! (1) vvv term
   if (delk/=0 .and. dell/=0 .and. delm/=0) then
      term = term + real(delk*dell*delm,wp) * Kvvv_abg
   end if

   ! (2) vv–a terms
   if (delk/=0 .and. dell/=0) then
      term = term + real(delk*dell,wp) * ( Kvv_ai_ab*dai_m + Kvv_aj_ab*daj_m )
   end if
   if (delk/=0 .and. delm/=0) then
      term = term + real(delk*delm,wp) * ( Kvv_ai_ag*dai_l + Kvv_aj_ag*daj_l )
   end if
   if (dell/=0 .and. delm/=0) then
      term = term + real(dell*delm,wp) * ( Kvv_ai_bg*dai_k + Kvv_aj_bg*daj_k )
   end if

   ! (3) v–aa / v–bb / v–ab terms
   if (delk/=0) then
      term = term + real(delk,wp) * ( Kv_aiai_a*(dai_l*dai_m) + Kv_ajaj_a*(daj_l*daj_m) &
         + Kv_aiaj_a*(dai_l*daj_m + daj_l*dai_m) )
   end if
   if (dell/=0) then
      term = term + real(dell,wp) * ( Kv_aiai_b*(dai_k*dai_m) + Kv_ajaj_b*(daj_k*daj_m) &
         + Kv_aiaj_b*(dai_k*daj_m + daj_k*dai_m) )
   end if
   if (delm/=0) then
      term = term + real(delm,wp) * ( Kv_aiai_g*(dai_k*dai_l) + Kv_ajaj_g*(daj_k*daj_l) &
         + Kv_aiaj_g*(dai_k*daj_l + daj_k*dai_l) )
   end if

   ! (4) radii-only cubic terms
   term = term + Kaiaiai*(dai_k*dai_l*dai_m) + Kajajaj*(daj_k*daj_l*daj_m)
   term = term + Kaiaiaj*(dai_k*dai_l*daj_m + dai_k*daj_l*dai_m + daj_k*dai_l*dai_m)
   term = term + Kaiajaj*(daj_k*daj_l*dai_m + daj_k*dai_l*daj_m + dai_k*daj_l*daj_m)

   ! (5) K_pq terms with brdr2 (three pairings)
   if (delm/=0) then
      term = term + real(delm,wp) * ( d2ai_kl * Kv_ai_g + d2aj_kl * Kv_aj_g )
   end if
   term = term + d2ai_kl * (Kaiai*dai_m + Kaiaj*daj_m) + d2aj_kl * (Kaiaj*dai_m + Kajaj*daj_m)

   if (dell/=0) then
      term = term + real(dell,wp) * ( d2ai_km * Kv_ai_b + d2aj_km * Kv_aj_b )
   end if
   term = term + d2ai_km * (Kaiai*dai_l + Kaiaj*daj_l) + d2aj_km * (Kaiaj*dai_l + Kajaj*daj_l)

   if (delk/=0) then
      term = term + real(delk,wp) * ( d2ai_lm * Kv_ai_a + d2aj_lm * Kv_aj_a )
   end if
   term = term + d2ai_lm * (Kaiai*dai_k + Kaiaj*daj_k) + d2aj_lm * (Kaiaj*dai_k + Kajaj*daj_k)

   ! (6) K_p * brdr3 term
   term = term + Kai*d3ai_klm + Kaj*d3aj_klm

   d3K_elem = term
end subroutine compute_p16_d3Kdr3

subroutine compute_p16_d4Kdr4(self, nat, xyz, brad, i, j, k, alpha, l, beta, m, gamma, n, delta, d4K_elem, &
      & brdr, brdr2, brdr3, brdr4)
   ! Element-wise 4th derivative for the P16 kernel:
   !   d4K_elem = d^4 K_ij / ( d r_{k,alpha} d r_{l,beta} d r_{m,gamma} d r_{n,delta} )
   !
   ! Same philosophy as Still:
   ! - compute only one entry, no (3,nat,3,nat,3,nat,3,nat) storage
   ! - cache ONLY the needed Born-radius derivative entries as scalars
   ! - for off-diagonal, match your *_full semantics: compute in ip>jp ordering

   class(p16_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                                      ! (3,nat)
   real(wp), intent(in) :: brad(:)                                        ! (nat)
   integer, intent(in) :: i, j
   integer, intent(in) :: k, alpha, l, beta, m, gamma, n, delta
   real(wp), intent(out) :: d4K_elem
   real(wp), contiguous, intent(in), optional :: brdr(:, :, :)                      ! (3,nat,nat)
   real(wp), contiguous, intent(in), optional :: brdr2(:, :, :, :, :)               ! (3,nat,3,nat,nat)
   real(wp), contiguous, intent(in), optional :: brdr3(:, :, :, :, :, :, :)         ! (3,nat,3,nat,3,nat,nat)
   real(wp), contiguous, intent(in), optional :: brdr4(:, :, :, :, :, :, :, :, :)   ! (3,nat,3,nat,3,nat,3,nat,nat)

   integer :: ip, jp
   integer :: delk, dell, delm, deln
   real(wp), parameter :: tiny_r = 1.0e-14_wp
   real(wp), parameter :: tiny = 1.0e-30_wp
   real(wp) :: c

   real(wp) :: v(3), r2, r

   real(wp) :: ai, aj
   real(wp) :: A0, u0, t0, invt0, z0, y0, g0, invg0
   real(wp) :: z2, z4, z8, z16, z15, z14, z13, z12
   real(wp) :: y1, y2, y3, y4

   ! self-term chain coeffs for keps/ai
   real(wp) :: self_f1, self_f2, self_f3, self_f4

   real(wp) :: I3(3,3)

   ! ===== coordinate ai/aj derivatives for this index quadruple =====
   real(wp) :: dai_k, dai_l, dai_m, dai_n
   real(wp) :: daj_k, daj_l, daj_m, daj_n

   real(wp) :: d2ai_kl, d2ai_km, d2ai_kn, d2ai_lm, d2ai_ln, d2ai_mn
   real(wp) :: d2aj_kl, d2aj_km, d2aj_kn, d2aj_lm, d2aj_ln, d2aj_mn

   real(wp) :: d3ai_klm, d3ai_kln, d3ai_kmn, d3ai_lmn
   real(wp) :: d3aj_klm, d3aj_kln, d3aj_kmn, d3aj_lmn

   real(wp) :: d4ai_klmn, d4aj_klmn

   ! ===== r2 derivatives (only 1st/2nd nonzero) =====
   real(wp) :: s_a, s_b, s_c, s_d
   real(wp) :: s_ab, s_ac, s_ad, s_bc, s_bd, s_cd

   ! ===== r derivatives up to 4th =====
   real(wp) :: r_a, r_b, r_c, r_d
   real(wp) :: r_ab, r_ac, r_ad, r_bc, r_bd, r_cd
   real(wp) :: r_abc, r_abd, r_acd, r_bcd
   real(wp) :: r_abcd

   real(wp) :: R1, R2c, R3c, R4c   ! coefficients for r = sqrt(s)

   ! ===== A = ai*aj derivatives up to 4th =====
   real(wp) :: A_a, A_b, A_c, A_d
   real(wp) :: A_ab, A_ac, A_ad, A_bc, A_bd, A_cd
   real(wp) :: A_abc, A_abd, A_acd, A_bcd
   real(wp) :: A_abcd

   ! ===== u = sqrt(A) derivatives up to 4th =====
   real(wp) :: u_a, u_b, u_c, u_d
   real(wp) :: u_ab, u_ac, u_ad, u_bc, u_bd, u_cd
   real(wp) :: u_abc, u_abd, u_acd, u_bcd
   real(wp) :: u_abcd

   real(wp) :: U1, U2c, U3c, U4c   ! coefficients for u = sqrt(A)

   ! ===== t = u + c*r and invt = 1/t derivatives up to 4th =====
   real(wp) :: t_a, t_b, t_c, t_d
   real(wp) :: t_ab, t_ac, t_ad, t_bc, t_bd, t_cd
   real(wp) :: t_abc, t_abd, t_acd, t_bcd
   real(wp) :: t_abcd

   real(wp) :: invt_a, invt_b, invt_c, invt_d
   real(wp) :: invt_ab, invt_ac, invt_ad, invt_bc, invt_bd, invt_cd
   real(wp) :: invt_abc, invt_abd, invt_acd, invt_bcd
   real(wp) :: invt_abcd

   ! ===== z = u * invt derivatives up to 4th =====
   real(wp) :: z_a, z_b, z_c, z_d
   real(wp) :: z_ab, z_ac, z_ad, z_bc, z_bd, z_cd
   real(wp) :: z_abc, z_abd, z_acd, z_bcd
   real(wp) :: z_abcd

   ! ===== y = z^16 derivatives up to 4th =====
   real(wp) :: y_a, y_b, y_c, y_d
   real(wp) :: y_ab, y_ac, y_ad, y_bc, y_bd, y_cd
   real(wp) :: y_abc, y_abd, y_acd, y_bcd
   real(wp) :: y_abcd

   ! ===== p = u*y derivatives up to 4th =====
   real(wp) :: p_a, p_b, p_c, p_d
   real(wp) :: p_ab, p_ac, p_ad, p_bc, p_bd, p_cd
   real(wp) :: p_abc, p_abd, p_acd, p_bcd
   real(wp) :: p_abcd

   ! ===== g = r + p derivatives =====
   real(wp) :: g_a, g_b, g_c, g_d
   real(wp) :: g_ab, g_ac, g_ad, g_bc, g_bd, g_cd
   real(wp) :: g_abc, g_abd, g_acd, g_bcd
   real(wp) :: g_abcd

   ! ===== invg = 1/g derivatives (up to 4th) =====
   real(wp) :: invg_abcd

   real(wp) :: term
   real(wp) :: sum1, sum2, sum3

   d4K_elem = 0.0_wp

   ! Safety: components must be 1..3
   if (alpha < 1 .or. alpha > 3) return
   if (beta  < 1 .or. beta  > 3) return
   if (gamma < 1 .or. gamma > 3) return
   if (delta < 1 .or. delta > 3) return

   I3 = 0.0_wp
   I3(1,1)=1.0_wp; I3(2,2)=1.0_wp; I3(3,3)=1.0_wp

   c = zetaP16o16   ! = zeta/16

   ! ==========================================================
   ! Diagonal self term: K_ii = keps / a_i
   ! ==========================================================
   if (i == j) then
      ai = brad(i)
      if (abs(ai) <= tiny) return

      self_f1 = -self%keps / (ai*ai)
      self_f2 =  2.0_wp * self%keps / (ai**3)
      self_f3 = -6.0_wp * self%keps / (ai**4)
      self_f4 =  24.0_wp * self%keps / (ai**5)

      ! Cache only the needed radii derivatives for this entry
      if (present(brdr) .and. present(brdr2) .and. present(brdr3) .and. present(brdr4)) then
         dai_k     = brdr(alpha, k, i)
         dai_l     = brdr(beta , l, i)
         dai_m     = brdr(gamma, m, i)
         dai_n     = brdr(delta, n, i)

         d2ai_kl   = brdr2(alpha, k, beta , l, i)
         d2ai_km   = brdr2(alpha, k, gamma, m, i)
         d2ai_kn   = brdr2(alpha, k, delta, n, i)
         d2ai_lm   = brdr2(beta , l, gamma, m, i)
         d2ai_ln   = brdr2(beta , l, delta, n, i)
         d2ai_mn   = brdr2(gamma, m, delta, n, i)

         d3ai_klm  = brdr3(alpha, k, beta , l, gamma, m, i)
         d3ai_kln  = brdr3(alpha, k, beta , l, delta, n, i)
         d3ai_kmn  = brdr3(alpha, k, gamma, m, delta, n, i)
         d3ai_lmn  = brdr3(beta , l, gamma, m, delta, n, i)

         d4ai_klmn = brdr4(alpha, k, beta, l, gamma, m, delta, n, i)
      else
         dai_k     = 0.0_wp
         dai_l     = 0.0_wp
         dai_m     = 0.0_wp
         dai_n     = 0.0_wp

         d2ai_kl   = 0.0_wp
         d2ai_km   = 0.0_wp
         d2ai_kn   = 0.0_wp
         d2ai_lm   = 0.0_wp
         d2ai_ln   = 0.0_wp
         d2ai_mn   = 0.0_wp

         d3ai_klm  = 0.0_wp
         d3ai_kln  = 0.0_wp
         d3ai_kmn  = 0.0_wp
         d3ai_lmn  = 0.0_wp

         d4ai_klmn = 0.0_wp
      end if

      term = 0.0_wp

      ! f'''' * (a1 a1 a1 a1)
      term = term + self_f4 * (dai_k*dai_l*dai_m*dai_n)

      ! f''' * sum_{(2,1,1)}  (6 terms)
      term = term + self_f3 * ( d2ai_kl*dai_m*dai_n + d2ai_km*dai_l*dai_n + d2ai_kn*dai_l*dai_m &
                              + d2ai_lm*dai_k*dai_n + d2ai_ln*dai_k*dai_m + d2ai_mn*dai_k*dai_l )

      ! f'' * [ sum_{(3,1)} (4 terms) + sum_{(2,2)} (3 terms) ]
      term = term + self_f2 * ( d3ai_klm*dai_n + d3ai_kln*dai_m + d3ai_kmn*dai_l + d3ai_lmn*dai_k &
                              + d2ai_kl*d2ai_mn + d2ai_km*d2ai_ln + d2ai_kn*d2ai_lm )

      ! f' * a4
      term = term + self_f1 * d4ai_klmn

      d4K_elem = term
      return
   end if

   ! ==========================================================
   ! Off-diagonal: match *_full semantics (computed for ip>jp)
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
   if (abs(ai) <= tiny .or. abs(aj) <= tiny) return

   v(:) = xyz(:, ip) - xyz(:, jp)
   r2   = dot_product(v, v)
   r    = sqrt(r2)
   if (r <= tiny_r) return

   A0 = ai*aj
   u0 = sqrt(A0)
   t0 = u0 + c*r
   invt0 = 1.0_wp / t0

   z0 = u0 * invt0

   ! z^16 etc by repeated squaring
   z2  = z0*z0
   z4  = z2*z2
   z8  = z4*z4
   z16 = z8*z8
   y0  = z16

   ! powers needed for y-derivative coefficients
   z15 = z16 / z0
   z14 = z15 / z0
   z13 = z14 / z0
   z12 = z13 / z0

   y1 = 16.0_wp * z15
   y2 = 240.0_wp * z14
   y3 = 3360.0_wp * z13
   y4 = 43680.0_wp * z12

   g0    = r + u0*y0
   invg0 = 1.0_wp / g0

   ! coefficients for r = sqrt(s) where s=r2
   R1  = 0.5_wp / r
   R2c = -0.25_wp / (r*r*r)
   R3c = 0.375_wp / (r**5)      ! 3/8
   R4c = -0.9375_wp / (r**7)    ! -15/16

   ! coefficients for u = sqrt(A)
   U1  = 0.5_wp / u0
   U2c = -0.25_wp / (u0**3)
   U3c = 0.375_wp / (u0**5)
   U4c = -0.9375_wp / (u0**7)

   ! Coordinate del-factors (for v = r_ip - r_jp)
   delk = 0; if (k==ip) delk=delk+1; if (k==jp) delk=delk-1
   dell = 0; if (l==ip) dell=dell+1; if (l==jp) dell=dell-1
   delm = 0; if (m==ip) delm=delm+1; if (m==jp) delm=delm-1
   deln = 0; if (n==ip) deln=deln+1; if (n==jp) deln=deln-1

   ! Cache ONLY the needed radii derivatives for this entry (ip and jp)
   if (present(brdr) .and. present(brdr2) .and. present(brdr3) .and. present(brdr4)) then
      dai_k     = brdr(alpha, k, ip);   daj_k     = brdr(alpha, k, jp)
      dai_l     = brdr(beta , l, ip);   daj_l     = brdr(beta , l, jp)
      dai_m     = brdr(gamma, m, ip);   daj_m     = brdr(gamma, m, jp)
      dai_n     = brdr(delta, n, ip);   daj_n     = brdr(delta, n, jp)

      d2ai_kl   = brdr2(alpha, k, beta , l, ip);   d2aj_kl   = brdr2(alpha, k, beta , l, jp)
      d2ai_km   = brdr2(alpha, k, gamma, m, ip);   d2aj_km   = brdr2(alpha, k, gamma, m, jp)
      d2ai_kn   = brdr2(alpha, k, delta, n, ip);   d2aj_kn   = brdr2(alpha, k, delta, n, jp)
      d2ai_lm   = brdr2(beta , l, gamma, m, ip);   d2aj_lm   = brdr2(beta , l, gamma, m, jp)
      d2ai_ln   = brdr2(beta , l, delta, n, ip);   d2aj_ln   = brdr2(beta , l, delta, n, jp)
      d2ai_mn   = brdr2(gamma, m, delta, n, ip);   d2aj_mn   = brdr2(gamma, m, delta, n, jp)

      d3ai_klm  = brdr3(alpha, k, beta , l, gamma, m, ip);   d3aj_klm  = brdr3(alpha, k, beta , l, gamma, m, jp)
      d3ai_kln  = brdr3(alpha, k, beta , l, delta, n, ip);   d3aj_kln  = brdr3(alpha, k, beta , l, delta, n, jp)
      d3ai_kmn  = brdr3(alpha, k, gamma, m, delta, n, ip);   d3aj_kmn  = brdr3(alpha, k, gamma, m, delta, n, jp)
      d3ai_lmn  = brdr3(beta , l, gamma, m, delta, n, ip);   d3aj_lmn  = brdr3(beta , l, gamma, m, delta, n, jp)

      d4ai_klmn = brdr4(alpha, k, beta, l, gamma, m, delta, n, ip)
      d4aj_klmn = brdr4(alpha, k, beta, l, gamma, m, delta, n, jp)

   else
      dai_k     = 0.0_wp;   daj_k     = 0.0_wp
      dai_l     = 0.0_wp;   daj_l     = 0.0_wp
      dai_m     = 0.0_wp;   daj_m     = 0.0_wp
      dai_n     = 0.0_wp;   daj_n     = 0.0_wp

      d2ai_kl   = 0.0_wp;   d2aj_kl   = 0.0_wp
      d2ai_km   = 0.0_wp;   d2aj_km   = 0.0_wp
      d2ai_kn   = 0.0_wp;   d2aj_kn   = 0.0_wp
      d2ai_lm   = 0.0_wp;   d2aj_lm   = 0.0_wp
      d2ai_ln   = 0.0_wp;   d2aj_ln   = 0.0_wp
      d2ai_mn   = 0.0_wp;   d2aj_mn   = 0.0_wp

      d3ai_klm  = 0.0_wp;   d3aj_klm  = 0.0_wp
      d3ai_kln  = 0.0_wp;   d3aj_kln  = 0.0_wp
      d3ai_kmn  = 0.0_wp;   d3aj_kmn  = 0.0_wp
      d3ai_lmn  = 0.0_wp;   d3aj_lmn  = 0.0_wp

      d4ai_klmn = 0.0_wp
      d4aj_klmn = 0.0_wp
   end if

   ! r2 first derivatives:
   s_a = 2.0_wp * real(delk,wp) * v(alpha)
   s_b = 2.0_wp * real(dell,wp) * v(beta)
   s_c = 2.0_wp * real(delm,wp) * v(gamma)
   s_d = 2.0_wp * real(deln,wp) * v(delta)

   ! r2 second derivatives:
   s_ab = 2.0_wp * real(delk*dell,wp) * I3(alpha,beta)
   s_ac = 2.0_wp * real(delk*delm,wp) * I3(alpha,gamma)
   s_ad = 2.0_wp * real(delk*deln,wp) * I3(alpha,delta)
   s_bc = 2.0_wp * real(dell*delm,wp) * I3(beta,gamma)
   s_bd = 2.0_wp * real(dell*deln,wp) * I3(beta,delta)
   s_cd = 2.0_wp * real(delm*deln,wp) * I3(gamma,delta)

   ! =========================================================
   ! r = sqrt(s) derivatives (s has only 1st/2nd)
   ! =========================================================
   r_a = R1 * s_a
   r_b = R1 * s_b
   r_c = R1 * s_c
   r_d = R1 * s_d

   r_ab = R1*s_ab + R2c*s_a*s_b
   r_ac = R1*s_ac + R2c*s_a*s_c
   r_ad = R1*s_ad + R2c*s_a*s_d
   r_bc = R1*s_bc + R2c*s_b*s_c
   r_bd = R1*s_bd + R2c*s_b*s_d
   r_cd = R1*s_cd + R2c*s_c*s_d

   r_abc = R2c*(s_ab*s_c + s_ac*s_b + s_bc*s_a) + R3c*s_a*s_b*s_c
   r_abd = R2c*(s_ab*s_d + s_ad*s_b + s_bd*s_a) + R3c*s_a*s_b*s_d
   r_acd = R2c*(s_ac*s_d + s_ad*s_c + s_cd*s_a) + R3c*s_a*s_c*s_d
   r_bcd = R2c*(s_bc*s_d + s_bd*s_c + s_cd*s_b) + R3c*s_b*s_c*s_d

   r_abcd = R2c*(s_ab*s_cd + s_ac*s_bd + s_ad*s_bc) &
          + R3c*( s_ab*s_c*s_d + s_ac*s_b*s_d + s_ad*s_b*s_c &
                + s_bc*s_a*s_d + s_bd*s_a*s_c + s_cd*s_a*s_b ) &
          + R4c*s_a*s_b*s_c*s_d

   ! =========================================================
   ! A = ai*aj derivatives (same structure as Still d4)
   ! =========================================================
   A_a = aj*dai_k + ai*daj_k
   A_b = aj*dai_l + ai*daj_l
   A_c = aj*dai_m + ai*daj_m
   A_d = aj*dai_n + ai*daj_n

   A_ab = aj*d2ai_kl + ai*d2aj_kl + (dai_k*daj_l + dai_l*daj_k)
   A_ac = aj*d2ai_km + ai*d2aj_km + (dai_k*daj_m + dai_m*daj_k)
   A_ad = aj*d2ai_kn + ai*d2aj_kn + (dai_k*daj_n + dai_n*daj_k)
   A_bc = aj*d2ai_lm + ai*d2aj_lm + (dai_l*daj_m + dai_m*daj_l)
   A_bd = aj*d2ai_ln + ai*d2aj_ln + (dai_l*daj_n + dai_n*daj_l)
   A_cd = aj*d2ai_mn + ai*d2aj_mn + (dai_m*daj_n + dai_n*daj_m)

   A_abc = aj*d3ai_klm + ai*d3aj_klm &
         + (d2ai_kl*daj_m + d2ai_km*daj_l + d2ai_lm*daj_k) &
         + (dai_k*d2aj_lm + dai_l*d2aj_km + dai_m*d2aj_kl)

   A_abd = aj*d3ai_kln + ai*d3aj_kln &
         + (d2ai_kl*daj_n + d2ai_kn*daj_l + d2ai_ln*daj_k) &
         + (dai_k*d2aj_ln + dai_l*d2aj_kn + dai_n*d2aj_kl)

   A_acd = aj*d3ai_kmn + ai*d3aj_kmn &
         + (d2ai_km*daj_n + d2ai_kn*daj_m + d2ai_mn*daj_k) &
         + (dai_k*d2aj_mn + dai_m*d2aj_kn + dai_n*d2aj_km)

   A_bcd = aj*d3ai_lmn + ai*d3aj_lmn &
         + (d2ai_lm*daj_n + d2ai_ln*daj_m + d2ai_mn*daj_l) &
         + (dai_l*d2aj_mn + dai_m*d2aj_ln + dai_n*d2aj_lm)

   A_abcd = aj*d4ai_klmn + ai*d4aj_klmn &
          + (d3ai_klm*daj_n + d3ai_kln*daj_m + d3ai_kmn*daj_l + d3ai_lmn*daj_k) &
          + (dai_k*d3aj_lmn + dai_l*d3aj_kmn + dai_m*d3aj_kln + dai_n*d3aj_klm) &
          + (d2ai_kl*d2aj_mn + d2ai_km*d2aj_ln + d2ai_kn*d2aj_lm &
           + d2ai_lm*d2aj_kn + d2ai_ln*d2aj_km + d2ai_mn*d2aj_kl)

   ! =========================================================
   ! u = sqrt(A) derivatives (composition)
   ! =========================================================
   u_a = U1 * A_a
   u_b = U1 * A_b
   u_c = U1 * A_c
   u_d = U1 * A_d

   u_ab = U1*A_ab + U2c*A_a*A_b
   u_ac = U1*A_ac + U2c*A_a*A_c
   u_ad = U1*A_ad + U2c*A_a*A_d
   u_bc = U1*A_bc + U2c*A_b*A_c
   u_bd = U1*A_bd + U2c*A_b*A_d
   u_cd = U1*A_cd + U2c*A_c*A_d

   u_abc = U1*A_abc + U2c*(A_ab*A_c + A_ac*A_b + A_bc*A_a) + U3c*A_a*A_b*A_c
   u_abd = U1*A_abd + U2c*(A_ab*A_d + A_ad*A_b + A_bd*A_a) + U3c*A_a*A_b*A_d
   u_acd = U1*A_acd + U2c*(A_ac*A_d + A_ad*A_c + A_cd*A_a) + U3c*A_a*A_c*A_d
   u_bcd = U1*A_bcd + U2c*(A_bc*A_d + A_bd*A_c + A_cd*A_b) + U3c*A_b*A_c*A_d

   sum1 = A_abc*A_d + A_abd*A_c + A_acd*A_b + A_bcd*A_a
   sum2 = A_ab*A_cd + A_ac*A_bd + A_ad*A_bc
   sum3 = A_ab*A_c*A_d + A_ac*A_b*A_d + A_ad*A_b*A_c + A_bc*A_a*A_d + A_bd*A_a*A_c + A_cd*A_a*A_b
   u_abcd = U1*A_abcd + U2c*(sum1 + sum2) + U3c*sum3 + U4c*A_a*A_b*A_c*A_d

   ! =========================================================
   ! t = u + c*r derivatives
   ! =========================================================
   t_a = u_a + c*r_a
   t_b = u_b + c*r_b
   t_c = u_c + c*r_c
   t_d = u_d + c*r_d

   t_ab = u_ab + c*r_ab
   t_ac = u_ac + c*r_ac
   t_ad = u_ad + c*r_ad
   t_bc = u_bc + c*r_bc
   t_bd = u_bd + c*r_bd
   t_cd = u_cd + c*r_cd

   t_abc = u_abc + c*r_abc
   t_abd = u_abd + c*r_abd
   t_acd = u_acd + c*r_acd
   t_bcd = u_bcd + c*r_bcd

   t_abcd = u_abcd + c*r_abcd

   ! =========================================================
   ! invt = 1/t derivatives (closed form)
   ! =========================================================
   invt_a = -(invt0*invt0) * t_a
   invt_b = -(invt0*invt0) * t_b
   invt_c = -(invt0*invt0) * t_c
   invt_d = -(invt0*invt0) * t_d

   invt_ab = 2.0_wp*(invt0**3)*t_a*t_b - (invt0*invt0)*t_ab
   invt_ac = 2.0_wp*(invt0**3)*t_a*t_c - (invt0*invt0)*t_ac
   invt_ad = 2.0_wp*(invt0**3)*t_a*t_d - (invt0*invt0)*t_ad
   invt_bc = 2.0_wp*(invt0**3)*t_b*t_c - (invt0*invt0)*t_bc
   invt_bd = 2.0_wp*(invt0**3)*t_b*t_d - (invt0*invt0)*t_bd
   invt_cd = 2.0_wp*(invt0**3)*t_c*t_d - (invt0*invt0)*t_cd

   invt_abc = -6.0_wp*(invt0**4)*t_a*t_b*t_c &
            + 2.0_wp*(invt0**3)*(t_ab*t_c + t_ac*t_b + t_bc*t_a) &
            - (invt0*invt0)*t_abc

   invt_abd = -6.0_wp*(invt0**4)*t_a*t_b*t_d &
            + 2.0_wp*(invt0**3)*(t_ab*t_d + t_ad*t_b + t_bd*t_a) &
            - (invt0*invt0)*t_abd

   invt_acd = -6.0_wp*(invt0**4)*t_a*t_c*t_d &
            + 2.0_wp*(invt0**3)*(t_ac*t_d + t_ad*t_c + t_cd*t_a) &
            - (invt0*invt0)*t_acd

   invt_bcd = -6.0_wp*(invt0**4)*t_b*t_c*t_d &
            + 2.0_wp*(invt0**3)*(t_bc*t_d + t_bd*t_c + t_cd*t_b) &
            - (invt0*invt0)*t_bcd

   sum1 = t_ab*t_c*t_d + t_ac*t_b*t_d + t_ad*t_b*t_c + t_bc*t_a*t_d + t_bd*t_a*t_c + t_cd*t_a*t_b
   sum2 = t_abc*t_d + t_abd*t_c + t_acd*t_b + t_bcd*t_a
   sum3 = t_ab*t_cd + t_ac*t_bd + t_ad*t_bc
   invt_abcd = 24.0_wp*(invt0**5)*t_a*t_b*t_c*t_d &
             - 6.0_wp*(invt0**4)*sum1 &
             + 2.0_wp*(invt0**3)*(sum2 + sum3) &
             - (invt0*invt0)*t_abcd

   ! =========================================================
   ! z = u * invt derivatives (product rule)
   ! =========================================================
   z_a = u_a*invt0 + u0*invt_a
   z_b = u_b*invt0 + u0*invt_b
   z_c = u_c*invt0 + u0*invt_c
   z_d = u_d*invt0 + u0*invt_d

   z_ab = u_ab*invt0 + u_a*invt_b + u_b*invt_a + u0*invt_ab
   z_ac = u_ac*invt0 + u_a*invt_c + u_c*invt_a + u0*invt_ac
   z_ad = u_ad*invt0 + u_a*invt_d + u_d*invt_a + u0*invt_ad
   z_bc = u_bc*invt0 + u_b*invt_c + u_c*invt_b + u0*invt_bc
   z_bd = u_bd*invt0 + u_b*invt_d + u_d*invt_b + u0*invt_bd
   z_cd = u_cd*invt0 + u_c*invt_d + u_d*invt_c + u0*invt_cd

   z_abc = u_abc*invt0 + (u_ab*invt_c + u_ac*invt_b + u_bc*invt_a) &
         + (u_a*invt_bc + u_b*invt_ac + u_c*invt_ab) + u0*invt_abc

   z_abd = u_abd*invt0 + (u_ab*invt_d + u_ad*invt_b + u_bd*invt_a) &
         + (u_a*invt_bd + u_b*invt_ad + u_d*invt_ab) + u0*invt_abd

   z_acd = u_acd*invt0 + (u_ac*invt_d + u_ad*invt_c + u_cd*invt_a) &
         + (u_a*invt_cd + u_c*invt_ad + u_d*invt_ac) + u0*invt_acd

   z_bcd = u_bcd*invt0 + (u_bc*invt_d + u_bd*invt_c + u_cd*invt_b) &
         + (u_b*invt_cd + u_c*invt_bd + u_d*invt_bc) + u0*invt_bcd

   z_abcd = u_abcd*invt0 &
          + (u_abc*invt_d + u_abd*invt_c + u_acd*invt_b + u_bcd*invt_a) &
          + (u_ab*invt_cd + u_ac*invt_bd + u_ad*invt_bc + u_bc*invt_ad + u_bd*invt_ac + u_cd*invt_ab) &
          + (u_a*invt_bcd + u_b*invt_acd + u_c*invt_abd + u_d*invt_abc) &
          + u0*invt_abcd

   ! =========================================================
   ! y = z^16 derivatives (composition)
   ! =========================================================
   y_a = y1 * z_a
   y_b = y1 * z_b
   y_c = y1 * z_c
   y_d = y1 * z_d

   y_ab = y1*z_ab + y2*z_a*z_b
   y_ac = y1*z_ac + y2*z_a*z_c
   y_ad = y1*z_ad + y2*z_a*z_d
   y_bc = y1*z_bc + y2*z_b*z_c
   y_bd = y1*z_bd + y2*z_b*z_d
   y_cd = y1*z_cd + y2*z_c*z_d

   y_abc = y1*z_abc + y2*(z_ab*z_c + z_ac*z_b + z_bc*z_a) + y3*z_a*z_b*z_c
   y_abd = y1*z_abd + y2*(z_ab*z_d + z_ad*z_b + z_bd*z_a) + y3*z_a*z_b*z_d
   y_acd = y1*z_acd + y2*(z_ac*z_d + z_ad*z_c + z_cd*z_a) + y3*z_a*z_c*z_d
   y_bcd = y1*z_bcd + y2*(z_bc*z_d + z_bd*z_c + z_cd*z_b) + y3*z_b*z_c*z_d

   sum1 = z_abc*z_d + z_abd*z_c + z_acd*z_b + z_bcd*z_a
   sum2 = z_ab*z_cd + z_ac*z_bd + z_ad*z_bc
   sum3 = z_ab*z_c*z_d + z_ac*z_b*z_d + z_ad*z_b*z_c + z_bc*z_a*z_d + z_bd*z_a*z_c + z_cd*z_a*z_b
   y_abcd = y1*z_abcd + y2*(sum1 + sum2) + y3*sum3 + y4*z_a*z_b*z_c*z_d

   ! =========================================================
   ! p = u*y derivatives (product rule)
   ! =========================================================
   p_a = u_a*y0 + u0*y_a
   p_b = u_b*y0 + u0*y_b
   p_c = u_c*y0 + u0*y_c
   p_d = u_d*y0 + u0*y_d

   p_ab = u_ab*y0 + u_a*y_b + u_b*y_a + u0*y_ab
   p_ac = u_ac*y0 + u_a*y_c + u_c*y_a + u0*y_ac
   p_ad = u_ad*y0 + u_a*y_d + u_d*y_a + u0*y_ad
   p_bc = u_bc*y0 + u_b*y_c + u_c*y_b + u0*y_bc
   p_bd = u_bd*y0 + u_b*y_d + u_d*y_b + u0*y_bd
   p_cd = u_cd*y0 + u_c*y_d + u_d*y_c + u0*y_cd

   p_abc = u_abc*y0 + (u_ab*y_c + u_ac*y_b + u_bc*y_a) &
         + (u_a*y_bc + u_b*y_ac + u_c*y_ab) + u0*y_abc

   p_abd = u_abd*y0 + (u_ab*y_d + u_ad*y_b + u_bd*y_a) &
         + (u_a*y_bd + u_b*y_ad + u_d*y_ab) + u0*y_abd

   p_acd = u_acd*y0 + (u_ac*y_d + u_ad*y_c + u_cd*y_a) &
         + (u_a*y_cd + u_c*y_ad + u_d*y_ac) + u0*y_acd

   p_bcd = u_bcd*y0 + (u_bc*y_d + u_bd*y_c + u_cd*y_b) &
         + (u_b*y_cd + u_c*y_bd + u_d*y_bc) + u0*y_bcd

   p_abcd = u_abcd*y0 &
          + (u_abc*y_d + u_abd*y_c + u_acd*y_b + u_bcd*y_a) &
          + (u_ab*y_cd + u_ac*y_bd + u_ad*y_bc + u_bc*y_ad + u_bd*y_ac + u_cd*y_ab) &
          + (u_a*y_bcd + u_b*y_acd + u_c*y_abd + u_d*y_abc) &
          + u0*y_abcd

   ! =========================================================
   ! g = r + p derivatives
   ! =========================================================
   g_a = r_a + p_a
   g_b = r_b + p_b
   g_c = r_c + p_c
   g_d = r_d + p_d

   g_ab = r_ab + p_ab
   g_ac = r_ac + p_ac
   g_ad = r_ad + p_ad
   g_bc = r_bc + p_bc
   g_bd = r_bd + p_bd
   g_cd = r_cd + p_cd

   g_abc = r_abc + p_abc
   g_abd = r_abd + p_abd
   g_acd = r_acd + p_acd
   g_bcd = r_bcd + p_bcd

   g_abcd = r_abcd + p_abcd

   ! =========================================================
   ! invg_abcd = d^4(1/g) using closed-form inverse identities
   ! =========================================================
   sum1 = g_ab*g_c*g_d + g_ac*g_b*g_d + g_ad*g_b*g_c + g_bc*g_a*g_d + g_bd*g_a*g_c + g_cd*g_a*g_b
   sum2 = g_abc*g_d + g_abd*g_c + g_acd*g_b + g_bcd*g_a
   sum3 = g_ab*g_cd + g_ac*g_bd + g_ad*g_bc

   invg_abcd =  24.0_wp*(invg0**5)*g_a*g_b*g_c*g_d &
             -  6.0_wp*(invg0**4)*sum1 &
             +  2.0_wp*(invg0**3)*(sum2 + sum3) &
             - (invg0*invg0)*g_abcd

   term = self%keps * invg_abcd

   d4K_elem = term
end subroutine compute_p16_d4Kdr4


!==============================================================================
! Coulomb kernel
!==============================================================================

pure subroutine add_coulomb_mat(self, nat, xyz, brad, Amat)
   class(coulomb_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: brad(:)
   real(wp), intent(inout) :: Amat(:, :)

   integer  :: i, j
   real(wp) :: vec(3), r1, invr

   ! Classic Coulomb interaction kernel: 1 / R
   do i = 1, nat
      do j = 1, i - 1
         vec(:) = xyz(:, i) - xyz(:, j)
         r1 = norm2(vec)

         ! NOTE: assumes no coincident atoms (r1 > 0)
         invr = 1.0_wp / r1

         Amat(i, j) = Amat(i, j) + self%keps*invr
         Amat(j, i) = Amat(j, i) + self%keps*invr
      enddo

      ! No Coulomb self-term; keep diagonal unchanged (and keep brad arg for interface compatibility)
      ! Amat(i, i) = Amat(i, i)
   enddo
end subroutine add_coulomb_mat

subroutine add_coulomb_deriv(self, nat, xyz, qat, brad, brdr, energy, gradient)
   class(coulomb_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: qat(:)
   real(wp), intent(in) :: brad(:)
   real(wp), contiguous, intent(in) :: brdr(:, :, :)
   real(wp), intent(out) :: energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)

   integer :: i, j
   real(wp) :: vec(3), r1, r2, invr, invr3
   real(wp) :: qq
   real(wp) :: dr(3)
   real(wp) :: e_coul
   real(wp), allocatable :: grddb(:)


   ! Keep for interface compatibility (unused for pure Coulomb kernel)
   allocate(grddb(nat), source = 0.0_wp)

   e_coul = 0.0_wp
   grddb(:) = 0.0_wp

   do i = 1, nat
      do j = 1, i - 1
         vec(:) = xyz(:, i) - xyz(:, j)
         r1 = norm2(vec)
         r2 = r1*r1

         ! NOTE: assumes no coincident atoms (r1 > 0)
         invr  = 1.0_wp / r1
         invr3 = invr / r2    ! = 1 / r^3

         qq = qat(i)*qat(j)

         ! Energy contribution: keps * q_i q_j / r_ij
         e_coul = e_coul + self%keps * qq * invr

         ! d/dr (1/r) = - r_vec / r^3
         dr = self%keps * invr3 * vec

         ! Gradient on coordinates (matches sign convention used in your Still routine)
         gradient(:, i) = gradient(:, i) - dr*qq
         gradient(:, j) = gradient(:, j) + dr*qq
      enddo

      ! No self energy term for classical Coulomb (and no born radius dependence)
   enddo

   ! Keep call for interface compatibility; grddb is zero so this is a no-op.
   call gemv(brdr, grddb, gradient, beta=1.0_wp)

   energy = e_coul
end subroutine add_coulomb_deriv

subroutine compute_coulomb_dKdr(self, nat, xyz, brad, i, j, k, alpha, dKdr_elem, brdr)
   class(coulomb_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                  ! (3,nat)
   real(wp), intent(in) :: brad(:)                    ! (nat)  (unused; interface compat)
   integer, intent(in) :: i, j
   integer, intent(in) :: k          ! atom index for derivative
   integer, intent(in) :: alpha      ! 1..3 cartesian component
   real(wp), intent(out) :: dKdr_elem
   real(wp), contiguous, intent(in), optional :: brdr(:, :, :)  ! (3,nat,nat) (unused; interface compat)

   real(wp) :: rvec(3), r2, r1, invr3
   real(wp) :: pref_pos

   dKdr_elem = 0.0_wp

   ! K_ii = 0 for classic Coulomb kernel => derivative is zero
   if (i == j) return

   ! Only atoms i or j contribute to dK_ij / dr_k
   if (k /= i .and. k /= j) return

   ! Safety: alpha must be 1..3
   if (alpha < 1 .or. alpha > 3) return

   rvec(:) = xyz(:, i) - xyz(:, j)
   r2      = dot_product(rvec, rvec)
   r1      = sqrt(r2)

   ! assumes r1 > 0 (no coincident atoms)
   invr3   = 1.0_wp / (r1*r2)     ! 1/r^3

   ! pref_pos = -self%keps * invr3
   pref_pos = -invr3

   if (k == i) then
      ! dK_ij / dr_i = pref_pos * (r_i - r_j)
      dKdr_elem = pref_pos * rvec(alpha)
   else
      ! dK_ij / dr_j = -pref_pos * (r_i - r_j)
      dKdr_elem = -pref_pos * rvec(alpha)
   end if

end subroutine compute_coulomb_dKdr

subroutine compute_coulomb_d2Kdr2(self, nat, xyz, brad, i, j, k, alpha, l, beta, d2K_elem, brdr, brdr2)
   ! Element-wise second derivative for the Coulomb kernel:
   !   d2K_elem = d^2 K_ij / ( d r_{k,alpha} d r_{l,beta} )
   !
   ! For K_ij = 1/|r_i - r_j| (or keps/|...| if you include self%keps),
   ! with v = r_i - r_j, the Hessian w.r.t. v is:
   !   d^2(1/r)/dv_a dv_b = (3 v_a v_b - r^2 delta_ab)/r^5
   !
   ! Then chain rule:
   !   d/d r_k = delk * d/dv   with delk = +1 if k=i, -1 if k=j, 0 otherwise
   !   => d^2/d r_k d r_l = delk*dell * d^2/dv^2

   class(coulomb_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                         ! (3,nat)
   real(wp), intent(in) :: brad(:)                           ! (nat) (unused; interface compat)
   integer, intent(in) :: i, j
   integer, intent(in) :: k, alpha
   integer, intent(in) :: l, beta
   real(wp), intent(out) :: d2K_elem
   real(wp), contiguous, intent(in), optional :: brdr(:, :, :)         ! (3,nat,nat) (unused; interface compat)
   real(wp), contiguous, intent(in), optional :: brdr2(:, :, :, :, :)  ! (3,nat,3,nat,nat) (unused; interface compat)

   integer :: delk, dell
   real(wp) :: v(3), r2, r1, r5
   real(wp) :: H_ab

   d2K_elem = 0.0_wp

   ! K_ii = 0 => all derivatives zero
   if (i == j) return

   ! Safety: components must be 1..3
   if (alpha < 1 .or. alpha > 3) return
   if (beta  < 1 .or. beta  > 3) return

   ! Only k,l in {i,j} can contribute
   delk = 0
   if (k == i) delk = delk + 1
   if (k == j) delk = delk - 1
   if (delk == 0) return

   dell = 0
   if (l == i) dell = dell + 1
   if (l == j) dell = dell - 1
   if (dell == 0) return

   v(:) = xyz(:, i) - xyz(:, j)
   r2   = dot_product(v, v)
   r1   = sqrt(r2)

   ! assumes r1 > 0
   r5 = r2*r2*r1

   ! H_ab = (3 v_a v_b - r^2 delta_ab) / r^5
   H_ab = (3.0_wp * v(alpha) * v(beta)) / r5
   if (alpha == beta) H_ab = H_ab - (r2 / r5)

   ! If your kernel is actually K_ij = self%keps / r, multiply by self%keps:
   ! H_ab = self%keps * H_ab

   d2K_elem = real(delk*dell, wp) * H_ab
end subroutine compute_coulomb_d2Kdr2

subroutine compute_coulomb_d3Kdr3_ij(self, nat, xyz, brad, brdr, brdr2, brdr3, i, j, d3Kdr3_ij)
   class(coulomb_kernel), intent(in) :: self
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
   real(wp) :: v(3), r2, r1, r5, r7
   real(wp) :: I3(3,3)
   real(wp) :: T(3,3,3)
   real(wp) :: fac, invr5, invr7
   real(wp) :: keps

   keps = 1.0_wp

   ! Identity matrix
   I3 = 0.0_wp
   I3(1,1)=1.0_wp; I3(2,2)=1.0_wp; I3(3,3)=1.0_wp

   d3Kdr3_ij(:, :, :, :, :, :) = 0.0_wp

   ! Diagonal: classic Coulomb kernel has no self term => all derivatives are zero
   if (i == j) then
      return
   end if

   ! Match *_full semantics: always build from ordered (ip>jp) pair
   if (i > j) then
      ip = i
      jp = j
   else
      ip = j
      jp = i
   end if

   ! Off-diagonal element (ip,jp): K = keps / r, r = |v|
   v(:) = xyz(:, ip) - xyz(:, jp)
   r2   = dot_product(v, v)
   r1   = sqrt(r2)

   ! NOTE: assumes no coincident atoms (r1 > 0)
   r5 = r2*r2*r1          ! r^5
   r7 = r5*r2             ! r^7

   invr5 = 1.0_wp / r5
   invr7 = 1.0_wp / r7

   ! Third derivative w.r.t. v-components:
   ! d^3(1/r)/dv_a dv_b dv_c =
   !    3 ( δ_ab v_c + δ_ac v_b + δ_bc v_a ) / r^5  -  15 v_a v_b v_c / r^7
   do alpha = 1,3
      do beta = 1,3
         do gamma = 1,3
            T(alpha,beta,gamma) = keps * ( &
               3.0_wp * ( I3(alpha,beta)*v(gamma) + I3(alpha,gamma)*v(beta) + I3(beta,gamma)*v(alpha) ) * invr5 &
               - 15.0_wp * v(alpha)*v(beta)*v(gamma) * invr7 )
         end do
      end do
   end do

   ! Map v-derivatives to coordinate derivatives using delk*dell*delm factors
   do k = 1, nat
      delk = 0; if (k==ip) delk = delk + 1; if (k==jp) delk = delk - 1
      do l = 1, nat
         dell = 0; if (l==ip) dell = dell + 1; if (l==jp) dell = dell - 1
         do m = 1, nat
            delm = 0; if (m==ip) delm = delm + 1; if (m==jp) delm = delm - 1

            if (delk /= 0 .and. dell /= 0 .and. delm /= 0) then
               fac = real(delk*dell*delm, wp)
               do alpha = 1,3
                  do beta = 1,3
                     do gamma = 1,3
                        d3Kdr3_ij(alpha,k,beta,l,gamma,m) = d3Kdr3_ij(alpha,k,beta,l,gamma,m) + fac*T(alpha,beta,gamma)
                     end do
                  end do
               end do
            end if

         end do
      end do
   end do

   ! No Born-radius chain rule terms for pure Coulomb kernel:
   ! brad, brdr, brdr2, brdr3 are present only for interface compatibility.
end subroutine compute_coulomb_d3Kdr3_ij

subroutine compute_coulomb_d3Kdr3(self, nat, xyz, brad, i, j, k, alpha, l, beta, m, gamma, d3K_elem, brdr, brdr2, brdr3)
   ! Element-wise third derivative for the Coulomb kernel:
   !   d3K_elem = d^3 K_ij / ( d r_{k,alpha} d r_{l,beta} d r_{m,gamma} )
   !
   ! Uses v = r_ip - r_jp with (ip,jp) being the ordered pair (max(i,j), min(i,j))
   ! to match your *_full semantics.
   !
   ! For K = keps / r, r = |v|:
   ! d^3(1/r)/dv_a dv_b dv_c =
   !    3 ( δ_ab v_c + δ_ac v_b + δ_bc v_a ) / r^5  -  15 v_a v_b v_c / r^7
   !
   ! Then map to coordinate derivatives with del factors:
   !   d/d r_k = delk * d/dv, delk = +1 if k=ip, -1 if k=jp, 0 otherwise
   !   => d^3/d r_k d r_l d r_m = delk*dell*delm * d^3/dv^3

   class(coulomb_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                               ! (3,nat)
   real(wp), intent(in) :: brad(:)                                 ! (nat) (unused; interface compat)
   integer, intent(in) :: i, j
   integer, intent(in) :: k, alpha
   integer, intent(in) :: l, beta
   integer, intent(in) :: m, gamma
   real(wp), intent(out) :: d3K_elem
   real(wp), contiguous, intent(in), optional :: brdr(:, :, :)               ! (3,nat,nat) (unused; interface compat)
   real(wp), contiguous, intent(in), optional :: brdr2(:, :, :, :, :)        ! (unused; interface compat)
   real(wp), contiguous, intent(in), optional :: brdr3(:, :, :, :, :, :, :)  ! (unused; interface compat)

   integer :: ip, jp
   integer :: delk, dell, delm
   real(wp) :: v(3), r2, r1, r5, r7
   real(wp) :: invr5, invr7
   real(wp) :: T_abg

   d3K_elem = 0.0_wp

   ! K_ii = 0 => all derivatives zero
   if (i == j) return

   ! Safety: components must be 1..3
   if (alpha < 1 .or. alpha > 3) return
   if (beta  < 1 .or. beta  > 3) return
   if (gamma < 1 .or. gamma > 3) return

   ! Match your ordered-pair semantics
   if (i > j) then
      ip = i
      jp = j
   else
      ip = j
      jp = i
   end if

   ! Only k,l,m in {ip,jp} can contribute
   delk = 0; if (k == ip) delk = delk + 1; if (k == jp) delk = delk - 1
   if (delk == 0) return

   dell = 0; if (l == ip) dell = dell + 1; if (l == jp) dell = dell - 1
   if (dell == 0) return

   delm = 0; if (m == ip) delm = delm + 1; if (m == jp) delm = delm - 1
   if (delm == 0) return

   v(:) = xyz(:, ip) - xyz(:, jp)
   r2   = dot_product(v, v)
   r1   = sqrt(r2)

   ! assumes r1 > 0
   r5 = r2*r2*r1
   r7 = r5*r2

   invr5 = 1.0_wp / r5
   invr7 = 1.0_wp / r7

   ! Element-wise third derivative w.r.t. v-components:
   ! T_abg = keps * [ 3(δ_ab v_g + δ_ag v_b + δ_bg v_a)/r^5 - 15 v_a v_b v_g / r^7 ]
   T_abg = 0.0_wp
   if (alpha == beta)  T_abg = T_abg + 3.0_wp * v(gamma) * invr5
   if (alpha == gamma) T_abg = T_abg + 3.0_wp * v(beta)  * invr5
   if (beta  == gamma) T_abg = T_abg + 3.0_wp * v(alpha) * invr5

   T_abg = T_abg - 15.0_wp * v(alpha)*v(beta)*v(gamma) * invr7

   ! If your kernel is actually K_ij = self%keps / r, multiply by self%keps:
   ! T_abg = self%keps * T_abg

   d3K_elem = real(delk*dell*delm, wp) * T_abg
end subroutine compute_coulomb_d3Kdr3

subroutine compute_coulomb_d4Kdr4(self, nat, xyz, brad, i, j, k, alpha, l, beta, m, gamma, n, delta, d4K_elem, &
      & brdr, brdr2, brdr3, brdr4)
   ! Element-wise 4th derivative for the Coulomb kernel:
   !   d4K_elem = d^4 K_ij / ( d r_{k,alpha} d r_{l,beta} d r_{m,gamma} d r_{n,delta} )
   !
   ! Matches your *_full semantics by always evaluating with ordered (ip>jp) pair.
   !
   ! For K = keps / r with v = r_ip - r_jp:
   !
   ! ∂_a∂_b∂_c∂_d (1/r) =
   !   (105 v_a v_b v_c v_d
   !    - 15 r^2 * s1
   !    +  3 r^4 * s2) / r^9
   !
   ! where
   !   s1 = sym(δ_ab v_c v_d)  (6 terms)
   !   s2 = sym(δ_ab δ_cd)     (3 terms)
   !
   ! Map to coordinate derivatives via del factors:
   !   d/d r_k = delk * d/dv, delk = +1 if k=ip, -1 if k=jp, 0 otherwise
   !   => d^4/d r_k d r_l d r_m d r_n = delk*dell*delm*deln * d^4/dv^4

   class(coulomb_kernel), intent(in) :: self
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)                                      ! (3,nat)
   real(wp), intent(in) :: brad(:)                                        ! (nat) (unused; interface compat)
   integer, intent(in) :: i, j
   integer, intent(in) :: k, alpha
   integer, intent(in) :: l, beta
   integer, intent(in) :: m, gamma
   integer, intent(in) :: n, delta
   real(wp), intent(out) :: d4K_elem
   real(wp), contiguous, intent(in), optional :: brdr(:, :, :)                      ! (unused)
   real(wp), contiguous, intent(in), optional :: brdr2(:, :, :, :, :)               ! (unused)
   real(wp), contiguous, intent(in), optional :: brdr3(:, :, :, :, :, :, :)         ! (unused)
   real(wp), contiguous, intent(in), optional :: brdr4(:, :, :, :, :, :, :, :, :)   ! (unused)

   integer :: ip, jp
   integer :: delk, dell, delm, deln
   real(wp) :: v(3), r2, r1, r4, r9
   real(wp) :: s1, s2
   real(wp) :: T4_abgd
   real(wp), parameter :: tiny_r = 1.0e-14_wp

   d4K_elem = 0.0_wp

   ! K_ii = 0 => all derivatives zero
   if (i == j) return

   ! Safety: components must be 1..3
   if (alpha < 1 .or. alpha > 3) return
   if (beta  < 1 .or. beta  > 3) return
   if (gamma < 1 .or. gamma > 3) return
   if (delta < 1 .or. delta > 3) return

   ! Match your ordered-pair semantics
   if (i > j) then
      ip = i
      jp = j
   else
      ip = j
      jp = i
   end if

   ! Only k,l,m,n in {ip,jp} can contribute
   delk = 0; if (k == ip) delk = delk + 1; if (k == jp) delk = delk - 1
   if (delk == 0) return

   dell = 0; if (l == ip) dell = dell + 1; if (l == jp) dell = dell - 1
   if (dell == 0) return

   delm = 0; if (m == ip) delm = delm + 1; if (m == jp) delm = delm - 1
   if (delm == 0) return

   deln = 0; if (n == ip) deln = deln + 1; if (n == jp) deln = deln - 1
   if (deln == 0) return

   v(:) = xyz(:, ip) - xyz(:, jp)
   r2   = dot_product(v, v)
   r1   = sqrt(r2)
   if (r1 <= tiny_r) return

   r4 = r2*r2
   r9 = r4*r4*r1   ! r^9

   ! s1 = sym(δ_ab v_c v_d) with indices (alpha,beta,gamma,delta)
   s1 = 0.0_wp
   if (alpha == beta)  s1 = s1 + v(gamma)*v(delta)
   if (alpha == gamma) s1 = s1 + v(beta) *v(delta)
   if (alpha == delta) s1 = s1 + v(beta) *v(gamma)
   if (beta  == gamma) s1 = s1 + v(alpha)*v(delta)
   if (beta  == delta) s1 = s1 + v(alpha)*v(gamma)
   if (gamma == delta) s1 = s1 + v(alpha)*v(beta)

   ! s2 = sym(δ_ab δ_cd) with indices (alpha,beta,gamma,delta)
   s2 = 0.0_wp
   if (alpha == beta  .and. gamma == delta) s2 = s2 + 1.0_wp
   if (alpha == gamma .and. beta  == delta) s2 = s2 + 1.0_wp
   if (alpha == delta .and. beta  == gamma) s2 = s2 + 1.0_wp

   T4_abgd = ( 105.0_wp * v(alpha)*v(beta)*v(gamma)*v(delta) &
             -  15.0_wp * r2 * s1 &
             +   3.0_wp * r4 * s2 ) / r9

   ! If your kernel is actually K_ij = self%keps / r, multiply by self%keps:
   ! T4_abgd = self%keps * T4_abgd

   d4K_elem = real(delk*dell*delm*deln, wp) * T4_abgd
end subroutine compute_coulomb_d4Kdr4



end module tblite_solvation_kernel
