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
      procedure :: add_kernel_mat => add_still_mat
      procedure :: add_kernel_deriv => add_still_deriv
      procedure :: add_kernel_deriv_multipole_contributions => add_still_deriv_multipole_contributions
      procedure :: kernel_d1_pair => still_d1_pair
      procedure :: kernel_d1_pair_dborn => still_d1_pair_dborn
      procedure :: kernel_d2_pair => still_d2_pair
      procedure :: kernel_d2_pair_dborn => still_d2_pair_dborn
      procedure :: kernel_d3_pair => still_d3_pair
      procedure :: kernel_d3_pair_dborn => still_d3_pair_dborn
      procedure :: kernel_d4_pair => still_d4_pair
      procedure :: kernel_d4_pair_dborn => still_d4_pair_dborn
      procedure :: kernel_d5_pair => still_d5_pair
   end type still_kernel

contains

   pure subroutine add_still_mat(self, nat, xyz, brad, Amat)
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
   end subroutine add_still_mat

   subroutine add_still_deriv(self, nat, xyz, qat, brad, brdr, energy, gradient)
      !> Instance of Still kernel
      class(still_kernel), intent(in) :: self
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

      integer :: i, j
      real(wp), parameter :: a4 = 0.25_wp
      real(wp) :: aa, r2, fgb2
      real(wp) :: qq, dd, expd, dfgb, dfgb2, dfgb3, egb, ap, bp
      real(wp) :: grddbi, grddbj
      real(wp) :: dr(3), r1, vec(3)
      real(wp), allocatable :: grddb(:)

      allocate (grddb(nat), source=0.0_wp)

      egb = 0.0_wp
      grddb(:) = 0.0_wp

      do i = 1, nat
         do j = 1, i-1
            vec(:) = xyz(:, i)-xyz(:, j)
            r1 = norm2(vec)
            r2 = r1*r1

            qq = qat(i)*qat(j)
            aa = brad(i)*brad(j)
            dd = a4*r2/aa
            expd = exp(-dd)
            fgb2 = r2+aa*expd
            dfgb2 = 1.0_wp/fgb2
            dfgb = sqrt(dfgb2)
            dfgb3 = dfgb2*dfgb*self%keps

            egb = egb+qq*self%keps*dfgb

            ! Frozen radii:
            ap = (1.0_wp-a4*expd)*dfgb3
            dr = ap*vec
            gradient(:, i) = gradient(:, i)-dr*qq
            gradient(:, j) = gradient(:, j)+dr*qq

            ! Born radii dependence:
            bp = -0.5_wp*expd*(1.0_wp+dd)*dfgb3
            grddbi = brad(j)*bp
            grddbj = brad(i)*bp
            grddb(i) = grddb(i)+grddbi*qq
            grddb(j) = grddb(j)+grddbj*qq
         end do

         bp = 1.0_wp/brad(i)
         qq = qat(i)*bp
         egb = egb+0.5_wp*qat(i)*qq*self%keps
         grddbi = -0.5_wp*self%keps*qq*bp
         grddb(i) = grddb(i)+grddbi*qat(i)
      end do

      call gemv(brdr, grddb, gradient, beta=1.0_wp)
      energy = egb
   end subroutine add_still_deriv

   subroutine add_still_deriv_multipole_contributions(self, nat, xyz, q_at, mu_at, q_at2, brad, brdr, gradient)
      use mctc_env, only: wp
      use tblite_blas, only: gemv
      implicit none

      !> Instance of Still kernel
      class(still_kernel), intent(in) :: self
      !> Number of atoms
      integer, intent(in) :: nat
      !> Cartesian coordinates (3,nat)
      real(wp), intent(in) :: xyz(:, :)
      !> Atomic partial charges (nat)
      real(wp), intent(in) :: q_at(:)
      !> Atomic dipole moments (3,nat)
      real(wp), intent(in) :: mu_at(:, :)
      !> Atomic quadrupole moments (6,nat), packed as (xx,xy,yy,xz,yz,zz), NOT doubled
      real(wp), intent(in) :: q_at2(:, :)
      !> Born radii (nat)
      real(wp), intent(in) :: brad(:)
      !> Born radii derivatives (3,nat,nat)
      real(wp), contiguous, intent(in) :: brdr(:, :, :)
      !> Nuclear gradient (3,nat)
      real(wp), contiguous, intent(inout) :: gradient(:, :)

      !> Loop indices for atoms
      integer :: iat, jat
      !> Loop indices for Cartesian directions and multipole components
      integer :: ic, ipk, iqk
      !> Auxiliary indices for tensor contractions
      integer :: ia1, ia2, ig1, ig2

      !> Distance vector between atoms
      real(wp) :: vec(3)
      !> Distance between atoms and its powers
      real(wp) :: rij, invr, invr2, invr3
      !> Unit vector and its derivatives
      real(wp) :: uvec(3), duvec(3, 3)
      !> Derivative of inverse distance
      real(wp) :: dinvr(3)

      !> Kernel derivatives up to 5th order
      real(wp) :: d1(3), d2(3, 3), d3(3, 3, 3), d4(3, 3, 3, 3), d5(3, 3, 3, 3, 3)
      !> First-order kernel derivatives with respect to Born radii
      real(wp) :: d1_br(3), d1_bs(3)
      !> Second-order kernel derivatives with respect to Born radii
      real(wp) :: d2_br(3, 3), d2_bs(3, 3)
      !> Third-order kernel derivatives with respect to Born radii
      real(wp) :: d3_br(3, 3, 3), d3_bs(3, 3, 3)
      !> Fourth-order kernel derivatives with respect to Born radii
      real(wp) :: d4_br(3, 3, 3, 3), d4_bs(3, 3, 3, 3)

      !> Charge on source atom
      real(wp) :: qsrc
      !> Dipole moments on response and source atoms
      real(wp) :: mresp(3), msrc(3)
      !> Quadrupole moments on response and source atoms (packed format)
      real(wp) :: qresp6(6), qsrc6(6)
      !> Born radii for response and source atoms
      real(wp) :: brad_resp, brad_src

      !> Generalized Born function and kernel prefactor
      real(wp) :: gpar, kpp, coef
      !> Derivatives of gpar, kpp, and coef with respect to coordinates
      real(wp) :: dgpar(3), dkpp(3), dcoef(3)
      !> Second derivatives and temporary storage
      real(wp) :: d2u(3), tmp3(3)

      !> Generalized Born function derivatives with respect to Born radii
      real(wp) :: gpar_br, gpar_bs, kpp_br, kpp_bs, coef_br, coef_bs
      !> Fifth-order contraction and its derivatives
      real(wp) :: s5, ds5(3), s5_br, s5_bs

      !> Temporary storage for Born radii derivatives
      real(wp), allocatable :: grddb(:)

      !> Temporary variables for Born radii derivative calculations
      real(wp) :: t, uu12, dtc, sym2, dwabge, wabge, brad_contrib

      !> Index arrays for unpacking quadrupole tensor components
      integer, parameter :: pa(6) = [1, 1, 2, 1, 2, 3]
      integer, parameter :: pb(6) = [1, 2, 2, 3, 3, 3]
      !> Factor array for symmetric tensor components
      integer, parameter :: pf(6) = [1, 2, 1, 2, 2, 1]

      allocate (grddb(nat), source=0.0_wp)

      ! Loop over ordered pairs: response = iat, source = jat
      do iat = 1, nat
         mresp(:) = mu_at(:, iat)
         qresp6(:) = q_at2(:, iat)
         brad_resp = brad(iat)

         do jat = 1, nat
            if (jat == iat) cycle

            qsrc = q_at(jat)
            msrc(:) = mu_at(:, jat)
            qsrc6(:) = q_at2(:, jat)
            brad_src = brad(jat)

            ! vec = xyz(:, jat) - xyz(:, iat) points from response to source
            vec = xyz(:, jat)-xyz(:, iat)
            rij = sqrt(dot_product(vec, vec))
            if (rij == 0.0_wp) cycle

            invr = 1.0_wp/rij
            invr2 = invr*invr
            invr3 = invr2*invr

            uvec = vec*invr

            ! Derivatives w.r.t. response position xyz(:, iat):
            ! du/dR_resp = -(I - u u^T) / r
            duvec = 0.0_wp
            do ic = 1, 3
               duvec(1, ic) = -(merge(1.0_wp, 0.0_wp, 1 == ic)-uvec(1)*uvec(ic))*invr
               duvec(2, ic) = -(merge(1.0_wp, 0.0_wp, 2 == ic)-uvec(2)*uvec(ic))*invr
               duvec(3, ic) = -(merge(1.0_wp, 0.0_wp, 3 == ic)-uvec(3)*uvec(ic))*invr
            end do

            ! d(1/r)/dR_resp = +u / r^2
            dinvr = invr2*uvec

            ! Kernel derivatives: response center = xyz(:, iat)
            call self%kernel_d1_pair(xyz(:, iat), xyz(:, jat), brad_resp, brad_src, d1)
            call self%kernel_d2_pair(xyz(:, iat), xyz(:, jat), brad_resp, brad_src, d2)
            call self%kernel_d3_pair(xyz(:, iat), xyz(:, jat), brad_resp, brad_src, d3)
            call self%kernel_d4_pair(xyz(:, iat), xyz(:, jat), brad_resp, brad_src, d4)
            call self%kernel_d5_pair(xyz(:, iat), xyz(:, jat), brad_resp, brad_src, d5)

            ! Derivatives w.r.t. Born radii
            call self%kernel_d1_pair_dborn(xyz(:, iat), xyz(:, jat), brad_resp, brad_src, d1_br, d1_bs)
            call self%kernel_d2_pair_dborn(xyz(:, iat), xyz(:, jat), brad_resp, brad_src, d2_br, d2_bs)
            call self%kernel_d3_pair_dborn(xyz(:, iat), xyz(:, jat), brad_resp, brad_src, d3_br, d3_bs)
            call self%kernel_d4_pair_dborn(xyz(:, iat), xyz(:, jat), brad_resp, brad_src, d4_br, d4_bs)
            call self%kernel_d3_pair(xyz(:, iat), xyz(:, jat), brad_resp, brad_src, d3)
            call self%kernel_d4_pair(xyz(:, iat), xyz(:, jat), brad_resp, brad_src, d4)
            call self%kernel_d5_pair(xyz(:, iat), xyz(:, jat), brad_resp, brad_src, d5)

            ! Derivatives w.r.t. Born radii
            call self%kernel_d1_pair_dborn(xyz(:, iat), xyz(:, jat), brad_resp, brad_src, d1_br, d1_bs)
            call self%kernel_d2_pair_dborn(xyz(:, iat), xyz(:, jat), brad_resp, brad_src, d2_br, d2_bs)
            call self%kernel_d3_pair_dborn(xyz(:, iat), xyz(:, jat), brad_resp, brad_src, d3_br, d3_bs)
            call self%kernel_d4_pair_dborn(xyz(:, iat), xyz(:, jat), brad_resp, brad_src, d4_br, d4_bs)

            ! Construct coefficient for SQ/QQ terms
            gpar = dot_product(d1, uvec)
            d2u = matmul(d2, uvec)
            kpp = dot_product(uvec, d2u)
            coef = (kpp+gpar*invr)/3.0_wp

            do ic = 1, 3
               dgpar(ic) = dot_product(d2(:, ic), uvec)+dot_product(d1, duvec(:, ic))
               tmp3 = matmul(d3(:, :, ic), uvec)
               dkpp(ic) = 2.0_wp*dot_product(duvec(:, ic), d2u)+dot_product(uvec, tmp3)
               dcoef(ic) = (dkpp(ic)+dgpar(ic)*invr+gpar*dinvr(ic))/3.0_wp
            end do

            ! Derivatives of coefficient w.r.t. Born radii
            gpar_br = dot_product(d1_br, uvec)
            gpar_bs = dot_product(d1_bs, uvec)
            kpp_br = dot_product(uvec, matmul(d2_br, uvec))
            kpp_bs = dot_product(uvec, matmul(d2_bs, uvec))
            coef_br = (kpp_br+gpar_br*invr)/3.0_wp
            coef_bs = (kpp_bs+gpar_bs*invr)/3.0_wp

            ! ============================================================
            ! SD: E += qsrc * mresp · d1
            ! d/dR_resp uses d2; Born chain uses d1_dborn
            ! ============================================================
            gradient(:, iat) = gradient(:, iat)+qsrc*matmul(transpose(d2), mresp)
            gradient(:, jat) = gradient(:, jat)-qsrc*matmul(transpose(d2), mresp)

            grddb(iat) = grddb(iat)+qsrc*dot_product(mresp, d1_br)
            grddb(jat) = grddb(jat)+qsrc*dot_product(mresp, d1_bs)

            ! ============================================================
            ! DD: E += 0.5 * mresp^T * (-d2) * msrc
            ! d/dR_resp uses -d3; Born chain uses -d2_dborn
            ! ============================================================
            do ic = 1, 3
               t = -0.5_wp*dot_product(mresp, matmul(d3(:, :, ic), msrc))
               gradient(ic, iat) = gradient(ic, iat)+t
               gradient(ic, jat) = gradient(ic, jat)-t
            end do

            grddb(iat) = grddb(iat)-0.5_wp*dot_product(mresp, matmul(d2_br, msrc))
            grddb(jat) = grddb(jat)-0.5_wp*dot_product(mresp, matmul(d2_bs, msrc))

            ! ============================================================
            ! DQ: E += mresp · [-(1/3) pf * d3] · qsrc6
            ! d/dR_resp uses d4; Born chain uses d3_dborn
            ! ============================================================
            do ic = 1, 3
               t = 0.0_wp
               do ipk = 1, 6
                  ia1 = pa(ipk)
                  ia2 = pb(ipk)
                  t = t+real(pf(ipk), wp)*qsrc6(ipk)*dot_product(mresp, d4(:, ia1, ia2, ic))
               end do
               t = -(1.0_wp/3.0_wp)*t
               gradient(ic, iat) = gradient(ic, iat)+t
               gradient(ic, jat) = gradient(ic, jat)-t
            end do

            brad_contrib = 0.0_wp
            do ipk = 1, 6
               ia1 = pa(ipk)
               ia2 = pb(ipk)
               brad_contrib = brad_contrib+real(pf(ipk), wp)*qsrc6(ipk)*dot_product(mresp, d3_br(:, ia1, ia2))
            end do
            grddb(iat) = grddb(iat)-(1.0_wp/3.0_wp)*brad_contrib

            brad_contrib = 0.0_wp
            do ipk = 1, 6
               ia1 = pa(ipk)
               ia2 = pb(ipk)
               brad_contrib = brad_contrib+real(pf(ipk), wp)*qsrc6(ipk)*dot_product(mresp, d3_bs(:, ia1, ia2))
            end do
            grddb(jat) = grddb(jat)-(1.0_wp/3.0_wp)*brad_contrib

            ! ============================================================
            ! SQ: E += qsrc * qresp6 · tc, tc_p = pf(p)*coef*u_a u_b
            ! d/dR_resp uses dcoef and du; Born chain uses coef_dborn
            ! ============================================================
            do ic = 1, 3
               t = 0.0_wp
               do ipk = 1, 6
                  ia1 = pa(ipk)
                  ia2 = pb(ipk)
                  uu12 = uvec(ia1)*uvec(ia2)
                  dtc = real(pf(ipk), wp)*(dcoef(ic)*uu12+coef*(duvec(ia1, ic)*uvec(ia2)+uvec(ia1)*duvec(ia2, ic)))
                  t = t+qresp6(ipk)*dtc
               end do
               gradient(ic, iat) = gradient(ic, iat)+qsrc*t
               gradient(ic, jat) = gradient(ic, jat)-qsrc*t
            end do

            brad_contrib = 0.0_wp
            do ipk = 1, 6
               ia1 = pa(ipk)
               ia2 = pb(ipk)
               uu12 = uvec(ia1)*uvec(ia2)
               brad_contrib = brad_contrib+qresp6(ipk)*(real(pf(ipk), wp)*coef_br*uu12)
            end do
            grddb(iat) = grddb(iat)+qsrc*brad_contrib

            brad_contrib = 0.0_wp
            do ipk = 1, 6
               ia1 = pa(ipk)
               ia2 = pb(ipk)
               uu12 = uvec(ia1)*uvec(ia2)
               brad_contrib = brad_contrib+qresp6(ipk)*(real(pf(ipk), wp)*coef_bs*uu12)
            end do
            grddb(jat) = grddb(jat)+qsrc*brad_contrib

            ! ============================================================
            ! QQ: E += 0.5 * qresp6^T * [pfpf*((1/3)d4 + 0.5*sym2*s5)] * qsrc6
            ! s5 = coef / r^2
            ! d/dR_resp uses d5 and ds5; Born chain uses d4_dborn and s5_dborn
            ! ============================================================
            s5 = coef*invr2
            do ic = 1, 3
               ds5(ic) = dcoef(ic)*invr2+coef*(2.0_wp*invr3*uvec(ic))
            end do
            s5_br = coef_br*invr2
            s5_bs = coef_bs*invr2

            do ic = 1, 3
               t = 0.0_wp
               do ipk = 1, 6
                  ia1 = pa(ipk)
                  ia2 = pb(ipk)
                  do iqk = 1, 6
                     ig1 = pa(iqk)
                     ig2 = pb(iqk)

                     sym2 = 0.0_wp
                     if (ia1 == ia2 .and. ig1 == ig2) sym2 = sym2+1.0_wp
                     if (ia1 == ig1 .and. ia2 == ig2) sym2 = sym2+1.0_wp
                     if (ia1 == ig2 .and. ia2 == ig1) sym2 = sym2+1.0_wp

                     dwabge = (1.0_wp/3.0_wp)*d5(ia1, ia2, ig1, ig2, ic)+0.5_wp*sym2*ds5(ic)

                     t = t+qresp6(ipk)*(real(pf(ipk)*pf(iqk), wp)*dwabge)*qsrc6(iqk)
                  end do
               end do
               t = 0.5_wp*t
               gradient(ic, iat) = gradient(ic, iat)+t
               gradient(ic, jat) = gradient(ic, jat)-t
            end do

            brad_contrib = 0.0_wp
            do ipk = 1, 6
               ia1 = pa(ipk)
               ia2 = pb(ipk)
               do iqk = 1, 6
                  ig1 = pa(iqk)
                  ig2 = pb(iqk)

                  sym2 = 0.0_wp
                  if (ia1 == ia2 .and. ig1 == ig2) sym2 = sym2+1.0_wp
                  if (ia1 == ig1 .and. ia2 == ig2) sym2 = sym2+1.0_wp
                  if (ia1 == ig2 .and. ia2 == ig1) sym2 = sym2+1.0_wp

                  wabge = (1.0_wp/3.0_wp)*d4_br(ia1, ia2, ig1, ig2)+0.5_wp*sym2*s5_br
                  brad_contrib = brad_contrib+qresp6(ipk)*(real(pf(ipk)*pf(iqk), wp)*wabge)*qsrc6(iqk)
               end do
            end do
            grddb(iat) = grddb(iat)+0.5_wp*brad_contrib

            brad_contrib = 0.0_wp
            do ipk = 1, 6
               ia1 = pa(ipk)
               ia2 = pb(ipk)
               do iqk = 1, 6
                  ig1 = pa(iqk)
                  ig2 = pb(iqk)

                  sym2 = 0.0_wp
                  if (ia1 == ia2 .and. ig1 == ig2) sym2 = sym2+1.0_wp
                  if (ia1 == ig1 .and. ia2 == ig2) sym2 = sym2+1.0_wp
                  if (ia1 == ig2 .and. ia2 == ig1) sym2 = sym2+1.0_wp

                  wabge = (1.0_wp/3.0_wp)*d4_bs(ia1, ia2, ig1, ig2)+0.5_wp*sym2*s5_bs
                  brad_contrib = brad_contrib+qresp6(ipk)*(real(pf(ipk)*pf(iqk), wp)*wabge)*qsrc6(iqk)
               end do
            end do
            grddb(jat) = grddb(jat)+0.5_wp*brad_contrib

         end do
      end do

      ! Born chain rule: dE/dR += sum_k (dE/db_k) (db_k/dR)
      call gemv(brdr, grddb, gradient, beta=1.0_wp)

      deallocate (grddb)

   end subroutine add_still_deriv_multipole_contributions

   subroutine still_d1_pair(self, rA, rB, bornA, bornB, d1)
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
   end subroutine still_d1_pair

   subroutine still_d1_pair_dborn(self, rA, rB, bornA, bornB, d1_bA, d1_bB)
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
   end subroutine still_d1_pair_dborn

   subroutine still_d2_pair(self, rA, rB, bornA, bornB, d2)
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
   end subroutine still_d2_pair

   subroutine still_d2_pair_dborn(self, rA, rB, bornA, bornB, d2_bA, d2_bB)
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
   end subroutine still_d2_pair_dborn

   subroutine still_d3_pair(self, rA, rB, bornA, bornB, d3)
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
   end subroutine still_d3_pair

   subroutine still_d3_pair_dborn(self, rA, rB, bornA, bornB, d3_bA, d3_bB)
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
   end subroutine still_d3_pair_dborn

   subroutine still_d4_pair(self, rA, rB, bornA, bornB, d4)
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
   end subroutine still_d4_pair

   subroutine still_d4_pair_dborn(self, rA, rB, bornA, bornB, d4_bA, d4_bB)
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
   end subroutine still_d4_pair_dborn

   subroutine still_d5_pair(self, rA, rB, bornA, bornB, d5)
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
   end subroutine still_d5_pair

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
