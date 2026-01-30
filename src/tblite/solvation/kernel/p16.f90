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

!> @file tblite/solvation/kernel/p16.f90
!> Provides the P16 kernel implementation for ALPB solvation.

module tblite_solvation_kernel_p16
   use mctc_env, only : wp
   use tblite_blas, only : gemv
   use tblite_solvation_kernel_type, only : kernel_type
   implicit none
   private

   public :: p16_kernel 

   type, extends(kernel_type) :: p16_kernel
   contains
      procedure :: add_kernel_mat       => add_p16_mat
      procedure :: add_kernel_deriv     => add_p16_deriv
      procedure :: kernel_d1_pair       => p16_d1_pair
      procedure :: kernel_d1_pair_dborn => p16_d1_pair_dborn
      procedure :: kernel_d2_pair       => p16_d2_pair
      procedure :: kernel_d2_pair_dborn => p16_d2_pair_dborn
      procedure :: kernel_d3_pair       => p16_d3_pair
      procedure :: kernel_d3_pair_dborn => p16_d3_pair_dborn
      procedure :: kernel_d4_pair       => p16_d4_pair
      procedure :: kernel_d4_pair_dborn => p16_d4_pair_dborn
      procedure :: kernel_d5_pair       => p16_d5_pair
   end type p16_kernel

   real(wp), parameter :: zetaP16    = 1.028_wp
   real(wp), parameter :: zetaP16o16 = zetaP16 / 16.0_wp

contains

subroutine add_p16_mat(self, nat, xyz, brad, Amat)
   !> Instance of P16 kernel
   class(p16_kernel), intent(in) :: self
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Born radii
   real(wp), intent(in) :: brad(:)
   !> Charge-charge interaction matrix
   real(wp), intent(inout) :: Amat(:, :)

   integer :: iat, jat
   real(wp) :: r1, ab, arg, fgb, dfgb, bp, vec(3)

   Amat = 0.0_wp

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
   !> Instance of P16 kernel
   class(p16_kernel), intent(in) :: self
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Born radii
   real(wp), intent(in) :: brad(:)
   !> Born radii derivatives
   real(wp), contiguous, intent(in) :: brdr(:, :, :)
   !> Solvation energy
   real(wp), intent(out) :: energy
   !> Molecular gradient
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


! subroutine add_p16_deriv(self, nat, xyz, qat, dpat, qpat, brad, brdr, energy, gradient)
!   class(p16_kernel), intent(in) :: self
!   integer, intent(in) :: nat
!   real(wp), intent(in) :: xyz(:,:)
!   real(wp), intent(in) :: qat(:)
!   real(wp), intent(in) :: dpat(:,:)
!   real(wp), intent(in) :: qpat(:,:)
!   real(wp), intent(in) :: brad(:)
!   real(wp), contiguous, intent(in) :: brdr(:,:,:)
!   real(wp), intent(out) :: energy
!   real(wp), contiguous, intent(inout) :: gradient(:,:)

  

! end subroutine add_p16_deriv


  subroutine p16_d1_pair(self, ra, rb, bornA, bornB, d1)
    ! d1(i) = ∂K/∂xa_i
    !> Instance of P16 kernel
    class(p16_kernel), intent(in) :: self
    !> Cartesian coordinates of atom A
    real(wp), intent(in)  :: ra(3)
    !> Cartesian coordinates of atom B
    real(wp), intent(in)  :: rb(3)
    !> Born radius of atom A
    real(wp), intent(in)  :: bornA
    !> Born radius of atom B
    real(wp), intent(in)  :: bornB
    !> First derivative of P16 kernel
    real(wp), intent(out) :: d1(3)

    real(wp) :: rv(3), r, k1, k2, k3, k4

    rv = ra - rb
    r = sqrt(dot_product(rv, rv))
    if (r == 0.0_wp) then
      d1 = 0.0_wp
      return
    end if

    call p16_radial_k_derivs(r, bornA, bornB, k1, k2, k3, k4)
    d1 = self%keps * (k1 / r) * rv
  end subroutine p16_d1_pair

  subroutine p16_d1_pair_dborn(self, ra, rb, bornA, bornB, d1_bA, d1_bB)
  !> Instance of P16 kernel
  class(p16_kernel), intent(in) :: self
  !> Cartesian coordinates of atom A
  real(wp), intent(in)  :: ra(3)
  !> Cartesian coordinates of atom B
  real(wp), intent(in)  :: rb(3)
  !> Born radius of atom A
  real(wp), intent(in)  :: bornA
  !> Born radius of atom B
  real(wp), intent(in)  :: bornB
  !> Derivative wrt Born radius of atom A
  real(wp), intent(out) :: d1_bA(3)
  !> Derivative wrt Born radius of atom B
  real(wp), intent(out) :: d1_bB(3)

  real(wp) :: rv(3), r, u(3), prod, ab, facA, facB
  real(wp) :: k1,k2,k3,k4, dk1ab,dk2ab,dk3ab,dk4ab

  rv = ra - rb
  r  = sqrt(dot_product(rv, rv))
  if (r == 0.0_wp) then
    d1_bA = 0.0_wp; d1_bB = 0.0_wp
    return
  end if
  u = rv / r

  prod = bornA*bornB
  if (prod <= 0.0_wp) then
    d1_bA = 0.0_wp; d1_bB = 0.0_wp
    return
  end if
  ab   = sqrt(prod)
  facA = 0.5_wp * bornB / ab
  facB = 0.5_wp * bornA / ab

  ! FIX: use keyword arguments for optional dk*ab outputs
  call p16_radial_k_derivs(r, bornA, bornB, k1, k2, k3, k4, &
                           dk1ab=dk1ab, dk2ab=dk2ab, dk3ab=dk3ab, dk4ab=dk4ab)

  d1_bA = self%keps * (dk1ab*facA) * u
  d1_bB = self%keps * (dk1ab*facB) * u
end subroutine p16_d1_pair_dborn




  subroutine p16_d2_pair(self, ra, rb, bornA, bornB, d2)
    ! d2(i,j) = ∂²K/∂xa_i∂xa_j
    !> Instance of P16 kernel
    class(p16_kernel), intent(in) :: self
    !> Cartesian coordinates of atom A
    real(wp), intent(in)  :: ra(3)
    !> Cartesian coordinates of atom B
    real(wp), intent(in)  :: rb(3)
    !> Born radius of atom A
    real(wp), intent(in)  :: bornA
    !> Born radius of atom B
    real(wp), intent(in)  :: bornB
    !> Second derivative of P16 kernel
    real(wp), intent(out) :: d2(3,3)

    real(wp) :: rv(3), r, u(3), k1, k2, k3, k4
    real(wp) :: a, b
    integer :: i, j

    rv = ra - rb
    r = sqrt(dot_product(rv, rv))
    if (r == 0.0_wp) then
      d2 = 0.0_wp
      return
    end if
    u = rv / r

    call p16_radial_k_derivs(r, bornA, bornB, k1, k2, k3, k4)

    ! Hessian for radial K(r):
    ! d2 = (k2 - k1/r) u⊗u + (k1/r) I
    a = k1 / r
    b = k2 - a

    do i = 1,3
      do j = 1,3
        d2(i,j) = self%keps * b * u(i)*u(j)
        if (i == j) d2(i,j) = d2(i,j) + self%keps * a
      end do
    end do
  end subroutine p16_d2_pair

subroutine p16_d2_pair_dborn(self, ra, rb, bornA, bornB, d2_bA, d2_bB)
  !> Instance of P16 kernel
  class(p16_kernel), intent(in) :: self
  !> Cartesian coordinates of atom A
  real(wp), intent(in)  :: ra(3)
  !> Cartesian coordinates of atom B
  real(wp), intent(in)  :: rb(3)
  !> Born radius of atom A
  real(wp), intent(in)  :: bornA
  !> Born radius of atom B
  real(wp), intent(in)  :: bornB
  !> Derivative wrt Born radius of atom A
  real(wp), intent(out) :: d2_bA(3,3)
  !> Derivative wrt Born radius of atom B
  real(wp), intent(out) :: d2_bB(3,3)

  real(wp) :: rv(3), r, u(3), prod, ab, facA, facB
  real(wp) :: k1,k2,k3,k4, dk1ab,dk2ab,dk3ab,dk4ab
  real(wp) :: daA, dbA, daB, dbB
  integer :: i,j

  rv = ra - rb
  r  = sqrt(dot_product(rv, rv))
  if (r == 0.0_wp) then
    d2_bA = 0.0_wp; d2_bB = 0.0_wp
    return
  end if
  u = rv / r

  prod = bornA*bornB
  if (prod <= 0.0_wp) then
    d2_bA = 0.0_wp; d2_bB = 0.0_wp
    return
  end if
  ab   = sqrt(prod)
  facA = 0.5_wp * bornB / ab
  facB = 0.5_wp * bornA / ab

  ! FIX: use keyword arguments for optional dk*ab outputs
  call p16_radial_k_derivs(r, bornA, bornB, k1, k2, k3, k4, &
                           dk1ab=dk1ab, dk2ab=dk2ab, dk3ab=dk3ab, dk4ab=dk4ab)

  ! a = k1/r ; b = k2 - a
  daA = (dk1ab*facA) / r
  daB = (dk1ab*facB) / r
  dbA = (dk2ab*facA) - daA
  dbB = (dk2ab*facB) - daB

  do i=1,3
    do j=1,3
      d2_bA(i,j) = self%keps * dbA * u(i)*u(j)
      d2_bB(i,j) = self%keps * dbB * u(i)*u(j)
      if (i==j) then
        d2_bA(i,j) = d2_bA(i,j) + self%keps * daA
        d2_bB(i,j) = d2_bB(i,j) + self%keps * daB
      end if
    end do
  end do
end subroutine p16_d2_pair_dborn

  subroutine p16_d3_pair(self, ra, rb, bornA, bornB, d3)
    ! d3(i,j,k) = ∂³K/∂xa_i∂xa_j∂xa_k
    !> Instance of P16 kernel
    class(p16_kernel), intent(in) :: self
    !> Cartesian coordinates of atom A
    real(wp), intent(in)  :: ra(3)
    !> Cartesian coordinates of atom B
    real(wp), intent(in)  :: rb(3)
    !> Born radius of atom A
    real(wp), intent(in)  :: bornA
    !> Born radius of atom B
    real(wp), intent(in)  :: bornB
    !> Third derivative of P16 kernel
    real(wp), intent(out) :: d3(3,3,3)

    real(wp) :: rv(3), r, u(3), k1, k2, k3, k4
    real(wp) :: b, c, d
    integer :: i, j, k

    rv = ra - rb
    r = sqrt(dot_product(rv, rv))
    if (r == 0.0_wp) then
      d3 = 0.0_wp
      return
    end if
    u = rv / r

    call p16_radial_k_derivs(r, bornA, bornB, k1, k2, k3, k4)

    ! d3 = c(δij u_k + δik u_j + δjk u_i) + d u_i u_j u_k
    b = k2 - k1/r
    c = b / r
    d = k3 - 3.0_wp*b/r

    d3 = 0.0_wp
    do i = 1,3
      do j = 1,3
        do k = 1,3
          d3(i,j,k) = self%keps * d * u(i)*u(j)*u(k)
          if (i == j) d3(i,j,k) = d3(i,j,k) + self%keps * c*u(k)
          if (i == k) d3(i,j,k) = d3(i,j,k) + self%keps * c*u(j)
          if (j == k) d3(i,j,k) = d3(i,j,k) + self%keps * c*u(i)
        end do
      end do
    end do
  end subroutine p16_d3_pair


subroutine p16_d3_pair_dborn(self, ra, rb, bornA, bornB, d3_bA, d3_bB)
  !> Instance of P16 kernel
  class(p16_kernel), intent(in) :: self
  !> Cartesian coordinates of atom A
  real(wp), intent(in)  :: ra(3)
  !> Cartesian coordinates of atom B
  real(wp), intent(in)  :: rb(3)
  !> Born radius of atom A
  real(wp), intent(in)  :: bornA
  !> Born radius of atom B
  real(wp), intent(in)  :: bornB
  !> Derivative wrt Born radius of atom A
  real(wp), intent(out) :: d3_bA(3,3,3)
  !> Derivative wrt Born radius of atom B
  real(wp), intent(out) :: d3_bB(3,3,3)

  real(wp) :: rv(3), r, u(3), prod, ab, facA, facB
  real(wp) :: k1,k2,k3,k4, dk1ab,dk2ab,dk3ab,dk4ab
  real(wp) :: dbA, dbB, dcA, dcB, ddA, ddB
  real(wp) :: bA, bB
  integer :: i,j,k

  rv = ra - rb
  r  = sqrt(dot_product(rv, rv))
  if (r == 0.0_wp) then
    d3_bA = 0.0_wp; d3_bB = 0.0_wp
    return
  end if
  u = rv / r

  prod = bornA*bornB
  if (prod <= 0.0_wp) then
    d3_bA = 0.0_wp; d3_bB = 0.0_wp
    return
  end if
  ab   = sqrt(prod)
  facA = 0.5_wp * bornB / ab
  facB = 0.5_wp * bornA / ab

  ! FIX: use keyword arguments for optional dk*ab outputs
  call p16_radial_k_derivs(r, bornA, bornB, k1, k2, k3, k4, &
                           dk1ab=dk1ab, dk2ab=dk2ab, dk3ab=dk3ab, dk4ab=dk4ab)

  bA  = k2 - k1/r
  bB  = bA   ! same b value; derivatives differ below

  dbA = (dk2ab*facA) - (dk1ab*facA)/r
  dbB = (dk2ab*facB) - (dk1ab*facB)/r

  ! c = b/r ; d = k3 - 3 b/r
  dcA = dbA / r
  dcB = dbB / r
  ddA = (dk3ab*facA) - 3.0_wp*dbA/r
  ddB = (dk3ab*facB) - 3.0_wp*dbB/r

  d3_bA = 0.0_wp
  d3_bB = 0.0_wp
  do i=1,3
    do j=1,3
      do k=1,3
        d3_bA(i,j,k) = self%keps * ddA * u(i)*u(j)*u(k)
        d3_bB(i,j,k) = self%keps * ddB * u(i)*u(j)*u(k)
        if (i==j) then
          d3_bA(i,j,k) = d3_bA(i,j,k) + self%keps * dcA*u(k)
          d3_bB(i,j,k) = d3_bB(i,j,k) + self%keps * dcB*u(k)
        end if
        if (i==k) then
          d3_bA(i,j,k) = d3_bA(i,j,k) + self%keps * dcA*u(j)
          d3_bB(i,j,k) = d3_bB(i,j,k) + self%keps * dcB*u(j)
        end if
        if (j==k) then
          d3_bA(i,j,k) = d3_bA(i,j,k) + self%keps * dcA*u(i)
          d3_bB(i,j,k) = d3_bB(i,j,k) + self%keps * dcB*u(i)
        end if
      end do
    end do
  end do
end subroutine p16_d3_pair_dborn
  


  subroutine p16_d4_pair(self, ra, rb, bornA, bornB, d4)
    ! d4(i,j,k,l) = ∂⁴K/∂xa_i∂xa_j∂xa_k∂xa_l
    !> Instance of P16 kernel
    class(p16_kernel), intent(in) :: self
    !> Cartesian coordinates of atom A
    real(wp), intent(in)  :: ra(3)
    !> Cartesian coordinates of atom B
    real(wp), intent(in)  :: rb(3)
    !> Born radius of atom A
    real(wp), intent(in)  :: bornA
    !> Born radius of atom B
    real(wp), intent(in)  :: bornB
    !> Fourth derivative of P16 kernel
    real(wp), intent(out) :: d4(3,3,3,3)

    real(wp) :: rv(3), r, u(3), k1, k2, k3, k4
    real(wp) :: b, e, f, g
    integer :: i, j, k, l
    logical :: dij, dik, dil, djk, djl, dkl

    rv = ra - rb
    r = sqrt(dot_product(rv, rv))
    if (r == 0.0_wp) then
      d4 = 0.0_wp
      return
    end if
    u = rv / r

    call p16_radial_k_derivs(r, bornA, bornB, k1, k2, k3, k4)

    b = k2 - k1/r
    e = b / (r*r)
    f = (k3 - 3.0_wp*b/r) / r
    g = k4 - 6.0_wp*k3/r + 15.0_wp*k2/(r*r) - 15.0_wp*k1/(r**3)

    d4 = 0.0_wp
    do i = 1,3
      do j = 1,3
        dij = (i==j)
        do k = 1,3
          dik = (i==k); djk = (j==k)
          do l = 1,3
            dil = (i==l); djl = (j==l); dkl = (k==l)

            d4(i,j,k,l) = self%keps * g*u(i)*u(j)*u(k)*u(l)

            if (dij) d4(i,j,k,l) = d4(i,j,k,l) + self%keps * f*u(k)*u(l)
            if (dik) d4(i,j,k,l) = d4(i,j,k,l) + self%keps * f*u(j)*u(l)
            if (dil) d4(i,j,k,l) = d4(i,j,k,l) + self%keps * f*u(j)*u(k)
            if (djk) d4(i,j,k,l) = d4(i,j,k,l) + self%keps * f*u(i)*u(l)
            if (djl) d4(i,j,k,l) = d4(i,j,k,l) + self%keps * f*u(i)*u(k)
            if (dkl) d4(i,j,k,l) = d4(i,j,k,l) + self%keps * f*u(i)*u(j)

            if (dij .and. dkl) d4(i,j,k,l) = d4(i,j,k,l) + self%keps * e
            if (dik .and. djl) d4(i,j,k,l) = d4(i,j,k,l) + self%keps * e
            if (dil .and. djk) d4(i,j,k,l) = d4(i,j,k,l) + self%keps * e
          end do
        end do
      end do
    end do
  end subroutine p16_d4_pair


 subroutine p16_d4_pair_dborn(self, ra, rb, bornA, bornB, d4_bA, d4_bB)
  !> Instance of P16 kernel
  class(p16_kernel), intent(in) :: self
  !> Cartesian coordinates of atom A
  real(wp), intent(in)  :: ra(3)
  !> Cartesian coordinates of atom B
  real(wp), intent(in)  :: rb(3)
  !> Born radius of atom A
  real(wp), intent(in)  :: bornA
  !> Born radius of atom B
  real(wp), intent(in)  :: bornB
  !> Derivative wrt Born radius of atom A
  real(wp), intent(out) :: d4_bA(3,3,3,3)
  !> Derivative wrt Born radius of atom B
  real(wp), intent(out) :: d4_bB(3,3,3,3)

  real(wp) :: rv(3), r, u(3), prod, ab, facA, facB
  real(wp) :: k1,k2,k3,k4, dk1ab,dk2ab,dk3ab,dk4ab
  real(wp) :: dbA, dbB, deA, deB, dfA, dfB, dgA, dgB
  integer :: i,j,k,l
  logical :: dij, dik, dil, djk, djl, dkl

  rv = ra - rb
  r  = sqrt(dot_product(rv, rv))
  if (r == 0.0_wp) then
    d4_bA = 0.0_wp; d4_bB = 0.0_wp
    return
  end if
  u = rv / r

  prod = bornA*bornB
  if (prod <= 0.0_wp) then
    d4_bA = 0.0_wp; d4_bB = 0.0_wp
    return
  end if
  ab   = sqrt(prod)
  facA = 0.5_wp * bornB / ab
  facB = 0.5_wp * bornA / ab

  ! FIX: use keyword arguments for optional dk*ab outputs
  call p16_radial_k_derivs(r, bornA, bornB, k1, k2, k3, k4, &
                           dk1ab=dk1ab, dk2ab=dk2ab, dk3ab=dk3ab, dk4ab=dk4ab)

  ! b = k2 - k1/r
  dbA = (dk2ab*facA) - (dk1ab*facA)/r
  dbB = (dk2ab*facB) - (dk1ab*facB)/r

  ! e = b/r^2
  deA = dbA / (r*r)
  deB = dbB / (r*r)

  ! f = (k3 - 3b/r)/r
  dfA = ((dk3ab*facA) - 3.0_wp*dbA/r) / r
  dfB = ((dk3ab*facB) - 3.0_wp*dbB/r) / r

  ! g = k4 - 6 k3/r +15 k2/r^2 -15 k1/r^3
  dgA = (dk4ab*facA) - 6.0_wp*(dk3ab*facA)/r + 15.0_wp*(dk2ab*facA)/(r*r) - 15.0_wp*(dk1ab*facA)/(r**3)
  dgB = (dk4ab*facB) - 6.0_wp*(dk3ab*facB)/r + 15.0_wp*(dk2ab*facB)/(r*r) - 15.0_wp*(dk1ab*facB)/(r**3)

  d4_bA = 0.0_wp
  d4_bB = 0.0_wp
  do i=1,3
    do j=1,3
      dij = (i==j)
      do k=1,3
        dik = (i==k); djk = (j==k)
        do l=1,3
          dil = (i==l); djl = (j==l); dkl = (k==l)

          d4_bA(i,j,k,l) = self%keps * dgA * u(i)*u(j)*u(k)*u(l)
          d4_bB(i,j,k,l) = self%keps * dgB * u(i)*u(j)*u(k)*u(l)

          if (dij) then
            d4_bA(i,j,k,l) = d4_bA(i,j,k,l) + self%keps * dfA * u(k)*u(l)
            d4_bB(i,j,k,l) = d4_bB(i,j,k,l) + self%keps * dfB * u(k)*u(l)
          end if
          if (dik) then
            d4_bA(i,j,k,l) = d4_bA(i,j,k,l) + self%keps * dfA * u(j)*u(l)
            d4_bB(i,j,k,l) = d4_bB(i,j,k,l) + self%keps * dfB * u(j)*u(l)
          end if
          if (dil) then
            d4_bA(i,j,k,l) = d4_bA(i,j,k,l) + self%keps * dfA * u(j)*u(k)
            d4_bB(i,j,k,l) = d4_bB(i,j,k,l) + self%keps * dfB * u(j)*u(k)
          end if
          if (djk) then
            d4_bA(i,j,k,l) = d4_bA(i,j,k,l) + self%keps * dfA * u(i)*u(l)
            d4_bB(i,j,k,l) = d4_bB(i,j,k,l) + self%keps * dfB * u(i)*u(l)
          end if
          if (djl) then
            d4_bA(i,j,k,l) = d4_bA(i,j,k,l) + self%keps * dfA * u(i)*u(k)
            d4_bB(i,j,k,l) = d4_bB(i,j,k,l) + self%keps * dfB * u(i)*u(k)
          end if
          if (dkl) then
            d4_bA(i,j,k,l) = d4_bA(i,j,k,l) + self%keps * dfA * u(i)*u(j)
            d4_bB(i,j,k,l) = d4_bB(i,j,k,l) + self%keps * dfB * u(i)*u(j)
          end if

          if (dij .and. dkl) then
            d4_bA(i,j,k,l) = d4_bA(i,j,k,l) + self%keps * deA
            d4_bB(i,j,k,l) = d4_bB(i,j,k,l) + self%keps * deB
          end if
          if (dik .and. djl) then
            d4_bA(i,j,k,l) = d4_bA(i,j,k,l) + self%keps * deA
            d4_bB(i,j,k,l) = d4_bB(i,j,k,l) + self%keps * deB
          end if
          if (dil .and. djk) then
            d4_bA(i,j,k,l) = d4_bA(i,j,k,l) + self%keps * deA
            d4_bB(i,j,k,l) = d4_bB(i,j,k,l) + self%keps * deB
          end if
        end do
      end do
    end do
  end do
end subroutine p16_d4_pair_dborn


    subroutine p16_d5_pair(self, ra, rb, bornA, bornB, d5)
    ! d5(i,j,k,l,m) = ∂⁵K/∂xa_i∂xa_j∂xa_k∂xa_l∂xa_m  (radial K(r))
    !> Instance of P16 kernel
    class(p16_kernel), intent(in) :: self
    !> Cartesian coordinates of atom A
    real(wp), intent(in)  :: ra(3)
    !> Cartesian coordinates of atom B
    real(wp), intent(in)  :: rb(3)
    !> Born radius of atom A
    real(wp), intent(in)  :: bornA
    !> Born radius of atom B
    real(wp), intent(in)  :: bornB
    !> Fifth derivative of P16 kernel
    real(wp), intent(out) :: d5(3,3,3,3,3)

    real(wp) :: rv(3), r, u(3)
    real(wp) :: k1, k2, k3, k4, k5
    real(wp) :: A, B, C
    integer  :: i, j, k, l, m
    logical  :: dij, dik, dil, dim, djk, djl, djm, dkl, dkm, dlm

    rv = ra - rb
    r  = sqrt(dot_product(rv, rv))
    if (r == 0.0_wp) then
      d5 = 0.0_wp
      return
    end if
    u = rv / r

    call p16_radial_k_derivs(r, bornA, bornB, k1, k2, k3, k4, k5)

    ! Isotropic decomposition:
    ! d5 = A uuuuu + B Sym(δ uuu) + C Sym(δδ u)
    C = (r*r*k3 - 3.0_wp*r*k2 + 3.0_wp*k1) / (r**4)
    B = (-15.0_wp*k1 + 15.0_wp*k2*r - 6.0_wp*k3*r*r + k4*r**3) / (r**4)
    A = k5 - 10.0_wp*k4/r + 45.0_wp*k3/(r*r) - 105.0_wp*k2/(r**3) + 105.0_wp*k1/(r**4)

    d5 = 0.0_wp

    do i = 1,3
      do j = 1,3
        dij = (i==j)
        do k = 1,3
          dik = (i==k); djk = (j==k)
          do l = 1,3
            dil = (i==l); djl = (j==l); dkl = (k==l)
            do m = 1,3
              dim = (i==m); djm = (j==m); dkm = (k==m); dlm = (l==m)

              ! A term
              d5(i,j,k,l,m) = self%keps * A * u(i)*u(j)*u(k)*u(l)*u(m)

              ! B terms: 10 permutations
              if (dij) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * B * u(k)*u(l)*u(m)
              if (dik) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * B * u(j)*u(l)*u(m)
              if (dil) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * B * u(j)*u(k)*u(m)
              if (dim) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * B * u(j)*u(k)*u(l)

              if (djk) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * B * u(i)*u(l)*u(m)
              if (djl) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * B * u(i)*u(k)*u(m)
              if (djm) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * B * u(i)*u(k)*u(l)

              if (dkl) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * B * u(i)*u(j)*u(m)
              if (dkm) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * B * u(i)*u(j)*u(l)
              if (dlm) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * B * u(i)*u(j)*u(k)

              ! C terms: 15 permutations
              if (dij .and. dkl) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * C * u(m)
              if (dij .and. dkm) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * C * u(l)
              if (dij .and. dlm) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * C * u(k)

              if (dik .and. djl) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * C * u(m)
              if (dik .and. djm) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * C * u(l)
              if (dik .and. dlm) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * C * u(j)

              if (dil .and. djk) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * C * u(m)
              if (dil .and. djm) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * C * u(k)
              if (dil .and. dkm) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * C * u(j)

              if (dim .and. djk) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * C * u(l)
              if (dim .and. djl) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * C * u(k)
              if (dim .and. dkl) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * C * u(j)

              if (djk .and. dlm) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * C * u(i)
              if (djl .and. dkm) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * C * u(i)
              if (djm .and. dkl) d5(i,j,k,l,m) = d5(i,j,k,l,m) + self%keps * C * u(i)

            end do
          end do
        end do
      end do
    end do

  end subroutine p16_d5_pair


!> Internal helper routine to compute radial derivatives of P16 kernel
pure subroutine p16_radial_k_derivs(r, bornA, bornB, k1, k2, k3, k4, k5, dk1ab, dk2ab, dk3ab, dk4ab)
  real(wp), intent(in)  :: r, bornA, bornB
  real(wp), intent(out) :: k1, k2, k3, k4
  real(wp), intent(out), optional :: k5                    ! d5K/dr5
  real(wp), intent(out), optional :: dk1ab, dk2ab, dk3ab, dk4ab  ! d(ki)/d(ab)

  real(wp) :: prod, ab, c, q, arg16
  real(wp) :: acoef, a1, a2, a3
  real(wp) :: y1, y2, y3, y4
  real(wp) :: f0, f1, f2, f3, f4
  real(wp) :: invf, invf2, invf3, invf4, invf5

  ! Only for k5
  real(wp) :: a4, y5, f5, invf6

  ! for derivatives wrt ab
  real(wp) :: darg16, dacoef, da1, da2, da3
  real(wp) :: dy1, dy2, dy3, dy4
  real(wp) :: df0, df1, df2, df3, df4
  real(wp) :: dinvf2, dinvf3, dinvf4, dinvf5
  real(wp) :: n2, n3, n4, dn2, dn3, dn4

  prod = bornA * bornB
  if (prod <= 0.0_wp .or. r == 0.0_wp) then
    k1 = 0.0_wp; k2 = 0.0_wp; k3 = 0.0_wp; k4 = 0.0_wp
    if (present(k5)) k5 = 0.0_wp
    if (present(dk1ab)) dk1ab = 0.0_wp
    if (present(dk2ab)) dk2ab = 0.0_wp
    if (present(dk3ab)) dk3ab = 0.0_wp
    if (present(dk4ab)) dk4ab = 0.0_wp
    return
  end if

  ab = sqrt(prod)
  c  = zetaP16o16
  q  = ab + c*r

  ! arg16 = (ab/q)^16
  arg16 = (ab/q)
  arg16 = arg16*arg16
  arg16 = arg16*arg16
  arg16 = arg16*arg16
  arg16 = arg16*arg16

  acoef = -16.0_wp * c / q
  a1    =  16.0_wp * c*c / (q*q)
  a2    = -32.0_wp * c**3 / (q**3)
  a3    =  96.0_wp * c**4 / (q**4)

  y1 = acoef * arg16
  y2 = (a1 + acoef*acoef) * arg16
  y3 = (a2 + 3.0_wp*acoef*a1 + acoef**3) * arg16
  y4 = (a3 + 4.0_wp*acoef*a2 + 3.0_wp*a1*a1 + 6.0_wp*acoef*acoef*a1 + acoef**4) * arg16

  f0 = r + ab*arg16
  f1 = 1.0_wp + ab*y1
  f2 = ab*y2
  f3 = ab*y3
  f4 = ab*y4

  invf  = 1.0_wp / f0
  invf2 = invf*invf
  invf3 = invf2*invf
  invf4 = invf3*invf
  invf5 = invf4*invf

  k1 = -f1 * invf2
  k2 = (2.0_wp*f1*f1 - f0*f2) * invf3
  k3 = (-6.0_wp*f1**3 + 6.0_wp*f0*f1*f2 - f0*f0*f3) * invf4
  k4 = ( 24.0_wp*f1**4 - 36.0_wp*f0*f1*f1*f2 + 6.0_wp*f0*f0*f2*f2 &
       + 8.0_wp*f0*f0*f1*f3 - f0**3*f4 ) * invf5

  if (present(k5)) then
    a4 = -384.0_wp * c**5 / (q**5)
    y5 = (a4 + 5.0_wp*acoef*a3 + 10.0_wp*a1*a2 + 10.0_wp*acoef*acoef*a2 &
        + 15.0_wp*acoef*a1*a1 + 10.0_wp*acoef**3*a1 + acoef**5) * arg16
    f5   = ab*y5
    invf6 = invf5*invf

    k5 = ( -120.0_wp*f1**5 + 240.0_wp*f0*f1**3*f2 - 90.0_wp*f0*f0*f1*f2*f2 &
         - 60.0_wp*f0*f0*f1*f1*f3 + 20.0_wp*f0**3*f2*f3 + 10.0_wp*f0**3*f1*f4 &
         - f0**4*f5 ) * invf6
  end if

  if (present(dk1ab) .or. present(dk2ab) .or. present(dk3ab) .or. present(dk4ab)) then
    darg16 = arg16 * 16.0_wp * (1.0_wp/ab - 1.0_wp/q)

    dacoef =  16.0_wp * c / (q*q)
    da1    = -32.0_wp * c*c / (q**3)
    da2    =  96.0_wp * c**3 / (q**4)
    da3    = -384.0_wp * c**4 / (q**5)

    dy1 = dacoef*arg16 + acoef*darg16
    dy2 = (da1 + 2.0_wp*acoef*dacoef)*arg16 + (a1 + acoef*acoef)*darg16
    dy3 = (da2 + 3.0_wp*(dacoef*a1 + acoef*da1) + 3.0_wp*acoef*acoef*dacoef)*arg16 &
        + (a2 + 3.0_wp*acoef*a1 + acoef**3)*darg16
    dy4 = (da3 + 4.0_wp*(dacoef*a2 + acoef*da2) + 6.0_wp*a1*da1 &
        + 6.0_wp*(2.0_wp*acoef*dacoef*a1 + acoef*acoef*da1) + 4.0_wp*acoef**3*dacoef)*arg16 &
        + (a3 + 4.0_wp*acoef*a2 + 3.0_wp*a1*a1 + 6.0_wp*acoef*acoef*a1 + acoef**4)*darg16

    df0 = arg16 + ab*darg16
    df1 = y1 + ab*dy1
    df2 = y2 + ab*dy2
    df3 = y3 + ab*dy3
    df4 = y4 + ab*dy4

    dinvf2 = -2.0_wp * invf3 * df0
    dinvf3 = -3.0_wp * invf4 * df0
    dinvf4 = -4.0_wp * (invf4*invf) * df0
    dinvf5 = -5.0_wp * (invf5*invf) * df0

    if (present(dk1ab)) dk1ab = -(df1*invf2 + f1*dinvf2)

    if (present(dk2ab)) then
      n2  = 2.0_wp*f1*f1 - f0*f2
      dn2 = 4.0_wp*f1*df1 - df0*f2 - f0*df2
      dk2ab = dn2*invf3 + n2*dinvf3
    end if

    if (present(dk3ab)) then
      n3  = -6.0_wp*f1**3 + 6.0_wp*f0*f1*f2 - f0*f0*f3
      dn3 = -18.0_wp*f1*f1*df1 &
          + 6.0_wp*(df0*f1*f2 + f0*df1*f2 + f0*f1*df2) &
          - (2.0_wp*f0*df0*f3 + f0*f0*df3)
      dk3ab = dn3*invf4 + n3*dinvf4
    end if

    if (present(dk4ab)) then
      n4  = 24.0_wp*f1**4 - 36.0_wp*f0*f1*f1*f2 + 6.0_wp*f0*f0*f2*f2 &
          + 8.0_wp*f0*f0*f1*f3 - f0**3*f4
      dn4 = 96.0_wp*f1**3*df1 &
          - 36.0_wp*(df0*f1*f1*f2 + f0*2.0_wp*f1*df1*f2 + f0*f1*f1*df2) &
          + 6.0_wp*(2.0_wp*f0*df0*f2*f2 + f0*f0*2.0_wp*f2*df2) &
          + 8.0_wp*(2.0_wp*f0*df0*f1*f3 + f0*f0*df1*f3 + f0*f0*f1*df3) &
          - (3.0_wp*f0*f0*df0*f4 + f0**3*df4)
      dk4ab = dn4*invf5 + n4*dinvf5
    end if
  end if
end subroutine p16_radial_k_derivs





end module tblite_solvation_kernel_p16
