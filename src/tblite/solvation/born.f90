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

!> @file tblite/solvation/born.f90
!> Provides a Born radii integrator

!> Integrator for Born radii based on the Onufriev-Bashford-Case model
module tblite_solvation_born
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_convert, only : aatoau
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   implicit none
   private

   public :: new_born_integrator

   ! -------- 2nd-order scalar AD type (w.r.t. r) --------
   type :: ad2
   real(wp) :: v
   real(wp) :: d1
   real(wp) :: d2
end type ad2

type :: ad3
   real(wp) :: v
   real(wp) :: d1
   real(wp) :: d2
   real(wp) :: d3
end type ad3
type :: ad4
   real(wp) :: v
   real(wp) :: d1
   real(wp) :: d2
   real(wp) :: d3
   real(wp) :: d4
end type ad4

interface operator(+)
   module procedure add_ad2, add_ad3, add_ad4
end interface
interface operator(-)
   module procedure sub_ad2, neg_ad2, sub_ad3, neg_ad3, sub_ad4, neg_ad4
end interface
interface operator(*)
   module procedure mul_ad2, mul_ad3, mul_ad4
end interface
interface operator(/)
   module procedure div_ad2, div_ad3, div_ad4
end interface



   !> Implementation of GBOBC integrator
   type, public :: born_integrator
      !> van der Waals radii of the particles
      real(wp), allocatable :: vdwr(:)
      !> pair descreening approximation radii
      real(wp), allocatable :: rho(:)
      !> offset van der Waals radii
      real(wp), allocatable :: svdw(:)
      !> cut-off radius for the Born radius NN list
      real(wp) :: lrcut
      !> Scaling factor for Born radii
      real(wp) :: born_scale
      !> Volume polynome correction, default parameters correspond to GBOBCII
      real(wp) :: obc(3)
   contains
      !> Calculate Born radii for a given geometry
      procedure :: get_rad
   end type born_integrator

   real(wp), parameter :: lrcut_default = 35.0_wp * aatoau
   real(wp), parameter :: born_scale_default = 1.0_wp
   real(wp), parameter :: born_offset_default = 0.0_wp
   real(wp), parameter :: descreening_default = 0.8_wp
   real(wp), parameter :: obc_default(3) = [1.0_wp, 0.8_wp, 4.85_wp]


contains

!> Create new Born radii integrator
subroutine new_born_integrator(self, mol, vdwrad, descreening, born_scale, born_offset, &
      & obc, rcutoff)
   !> Instance of the Born integrator
   type(born_integrator), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Van-der-Waals Radii
   real(wp), intent(in) :: vdwRad(:)
   !> Dielectric descreening parameter
   real(wp), intent(in), optional :: descreening(:)
   !> Scaling factor for Born radii
   real(wp), intent(in), optional :: born_scale
   !> Offset parameter for Born radii integration
   real(wp), intent(in), optional :: born_offset
   !> GBOBC integrator parameters
   real(wp), intent(in), optional :: obc(3)
   !> Real-space cutoff for Born radii integration
   real(wp), intent(in), optional :: rCutoff

   self%lrcut = lrcut_default
   if (present(rCutoff)) then
      self%lrcut = rCutoff
   end if

   self%born_scale = born_scale_default
   if (present(born_scale)) then
      self%born_scale = born_scale
   end if

   self%obc = obc_default
   if (present(obc)) then
      self%obc = obc
   end if

   self%vdwr = vdwRad(mol%id)

   if (present(descreening)) then
      self%rho = self%vdwr * descreening(mol%id)
   else
      self%rho = self%vdwr * descreening_default
   end if

   if (present(born_offset)) then
      self%svdw = self%vdwr - born_offset
   else
      self%svdw = self%vdwr - born_offset_default
   end if
end subroutine new_born_integrator

!> Calculate Born radii
subroutine get_rad(self, mol, rad, draddr, dradd2r, dradd3r, dradd4r)
   class(born_integrator), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: rad(:)
   real(wp), intent(out), optional :: draddr(:, :, :)
   real(wp), intent(out), optional :: dradd2r(:, :, :, :, :)
   real(wp), intent(out), optional :: dradd3r(:, :, :, :, :, :, :)
   real(wp), intent(out), optional :: dradd4r(:, :, :, :, :, :, :, :, :)

   type(adjacency_list) :: list
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), allocatable :: brdr(:, :, :)
   real(wp), allocatable :: brddr(:, :, :, :, :)
   real(wp), allocatable :: brd3dr(:, :, :, :, :, :, :)
   real(wp), allocatable :: brd4dr(:, :, :, :, :, :, :, :, :)

   call new_adjacency_list(list, mol, trans, self%lrcut)

   allocate(brdr(3, mol%nat, mol%nat))
   if (present(dradd2r) .or. present(dradd3r) .or. present(dradd4r)) &
      allocate(brddr(3, mol%nat, 3, mol%nat, mol%nat))
   if (present(dradd3r) .or. present(dradd4r)) &
      allocate(brd3dr(3, mol%nat, 3, mol%nat, 3, mol%nat, mol%nat))
   if (present(dradd4r)) &
      allocate(brd4dr(3, mol%nat, 3, mol%nat, 3, mol%nat, 3, mol%nat, mol%nat))

   if (present(dradd4r)) then
      call compute_bornr(mol%nat, mol%xyz, list, self%vdwr, self%rho, self%svdw, self%born_scale, self%obc, &
         & rad, brdr, brddr, brd3dr, brd4dr)
   else if (present(dradd3r)) then
      call compute_bornr(mol%nat, mol%xyz, list, self%vdwr, self%rho, self%svdw, self%born_scale, self%obc, &
         & rad, brdr, brddr, brd3dr)
   else if (present(dradd2r)) then
      call compute_bornr(mol%nat, mol%xyz, list, self%vdwr, self%rho, self%svdw, self%born_scale, self%obc, &
         & rad, brdr, brddr)
   else
      call compute_bornr(mol%nat, mol%xyz, list, self%vdwr, self%rho, self%svdw, self%born_scale, self%obc, &
         & rad, brdr)
   end if

   if (present(draddr))  draddr(:, :, :) = brdr
   if (present(dradd2r)) dradd2r(:, :, :, :, :) = brddr
   if (present(dradd3r)) dradd3r(:, :, :, :, :, :, :) = brd3dr
   if (present(dradd4r)) dradd4r(:, :, :, :, :, :, :, :, :) = brd4dr
end subroutine get_rad




subroutine compute_bornr(nat, xyz, list, vdwr, rho, svdw, c1, obc, brad, brdr, brddr, brd3dr, brd4dr)
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   type(adjacency_list), intent(in) :: list
   real(wp), intent(in) :: vdwr(:), rho(:), svdw(:), c1, obc(3)
   real(wp), intent(out) :: brad(:)
   real(wp), intent(out) :: brdr(:, :, :)
   real(wp), intent(out), optional :: brddr(:, :, :, :, :)
   real(wp), intent(out), optional :: brd3dr(:, :, :, :, :, :, :)
   real(wp), intent(out), optional :: brd4dr(:, :, :, :, :, :, :, :, :)

   integer :: iat
   real(wp) :: R, fp, fpp, fppp, fpppp

   if (present(brd4dr)) then
      call compute_psi(nat, xyz, list, vdwr, rho, brad, brdr, brddr, brd3dr, brd4dr)
      call compute_bornr_d4(nat, vdwr, svdw, c1, obc, brad, brdr, brddr, brd3dr, brd4dr)
   else if (present(brd3dr)) then
      call compute_psi(nat, xyz, list, vdwr, rho, brad, brdr, brddr, brd3dr)
      call compute_bornr_d3(nat, vdwr, svdw, c1, obc, brad, brdr, brddr, brd3dr)
   else if (present(brddr)) then
      call compute_psi(nat, xyz, list, vdwr, rho, brad, brdr, brddr)
      call compute_bornr_d2(nat, vdwr, svdw, c1, obc, brad, brdr, brddr)
   else
      call compute_psi(nat, xyz, list, vdwr, rho, brad, brdr)
      do iat = 1, nat
         call obc_map_d1d2d3d4(brad(iat), svdw(iat), vdwr(iat), c1, obc, R, fp, fpp, fppp, fpppp)
         brad(iat) = R
         brdr(:, :, iat) = fp * brdr(:, :, iat)
      end do
   end if
end subroutine compute_bornr



subroutine compute_bornr_d2(nat, vdwr, svdw, c1, obc, brad, brdr, brddr)
   integer, intent(in) :: nat
   real(wp), intent(in) :: vdwr(:), svdw(:), c1, obc(3)
   real(wp), intent(inout) :: brad(:)
   real(wp), intent(inout) :: brdr(:, :, :)
   real(wp), intent(inout) :: brddr(:, :, :, :, :)

   integer :: iat, k, l, a, b
   real(wp) :: R, fp, fpp
   real(wp) :: gpsi(3, nat)

   do iat = 1, nat
      gpsi(:, :) = brdr(:, :, iat)

      call obc_map_d1d2(brad(iat), svdw(iat), vdwr(iat), c1, obc, R, fp, fpp)

      brad(iat) = R
      brdr(:, :, iat) = fp * brdr(:, :, iat)

      do k = 1, nat
         do l = 1, nat
            do a = 1, 3
               do b = 1, 3
                  brddr(a, k, b, l, iat) = fp * brddr(a, k, b, l, iat) + fpp * gpsi(a, k) * gpsi(b, l)
               end do
            end do
         end do
      end do
   end do
end subroutine compute_bornr_d2

pure subroutine obc_map_d1d2(psi, svdwi, vdwri, c1, obc, R, fp, fpp)
   real(wp), intent(in) :: psi, svdwi, vdwri, c1, obc(3)
   real(wp), intent(out) :: R, fp, fpp

   real(wp) :: s, B, u, up, upp
   real(wp) :: t, ch, sech2
   real(wp) :: s1, v1, A
   real(wp) :: u_psi, u_psipsi

   s  = 0.5_wp * svdwi
   B  = s * psi

   ! u(B) = a B - b B^2 + c B^3
   u   = B * (obc(1) + B * (obc(3) * B - obc(2)))
   up  = obc(1) - 2.0_wp * obc(2) * B + 3.0_wp * obc(3) * B * B
   upp = -2.0_wp * obc(2) + 6.0_wp * obc(3) * B

   t  = tanh(u)
   ch = cosh(u)
   sech2 = 1.0_wp / (ch * ch)

   s1 = 1.0_wp / svdwi
   v1 = 1.0_wp / vdwri

   A = s1 - v1 * t

   R = c1 / A

   u_psi    = up  * s
   u_psipsi = upp * s * s

   fp  = c1 * v1 * sech2 * u_psi / (A * A)

   fpp = c1 * v1 * sech2 / (A * A) * (u_psipsi - 2.0_wp * t * u_psi * u_psi) &
       + 2.0_wp * c1 * v1 * v1 * (sech2 * sech2) * (u_psi * u_psi) / (A * A * A)
end subroutine obc_map_d1d2


subroutine compute_bornr_d3(nat, vdwr, svdw, c1, obc, brad, brdr, brddr, brd3dr)
   integer, intent(in) :: nat
   real(wp), intent(in) :: vdwr(:), svdw(:), c1, obc(3)
   real(wp), intent(inout) :: brad(:)
   real(wp), intent(inout) :: brdr(:, :, :)
   real(wp), intent(inout) :: brddr(:, :, :, :, :)
   real(wp), intent(inout) :: brd3dr(:, :, :, :, :, :, :)

   integer :: iat, k, l, m, a, b, c
   real(wp) :: R, fp, fpp, fppp
   real(wp) :: g1, g2, g3
   real(wp) :: h12, h13, h23

   do iat = 1, nat
      call obc_map_d1d2d3(brad(iat), svdw(iat), vdwr(iat), c1, obc, R, fp, fpp, fppp)

      ! Third derivative first (needs old brdr/brddr)
      do k = 1, nat
         do l = 1, nat
            do m = 1, nat
               do a = 1, 3
                  g1 = brdr(a, k, iat)
                  do b = 1, 3
                     g2 = brdr(b, l, iat)
                     h12 = brddr(a, k, b, l, iat)
                     do c = 1, 3
                        g3 = brdr(c, m, iat)
                        h13 = brddr(a, k, c, m, iat)
                        h23 = brddr(b, l, c, m, iat)

                        brd3dr(a, k, b, l, c, m, iat) = &
                           fp   * brd3dr(a, k, b, l, c, m, iat) &
                         + fpp  * (h12*g3 + h13*g2 + h23*g1) &
                         + fppp * (g1*g2*g3)
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! Hessian update
      do k = 1, nat
         do l = 1, nat
            do a = 1, 3
               do b = 1, 3
                  brddr(a, k, b, l, iat) = fp * brddr(a, k, b, l, iat) + fpp * brdr(a, k, iat) * brdr(b, l, iat)
               end do
            end do
         end do
      end do

      ! Gradient update
      brdr(:, :, iat) = fp * brdr(:, :, iat)
      brad(iat) = R
   end do
end subroutine compute_bornr_d3


pure subroutine obc_map_d1d2d3(psi, svdwi, vdwri, c1, obc, R, fp, fpp, fppp)
   real(wp), intent(in) :: psi, svdwi, vdwri, c1, obc(3)
   real(wp), intent(out) :: R, fp, fpp, fppp

   real(wp) :: s, B, u, up, upp, uppp
   real(wp) :: t, ch, sech2
   real(wp) :: s1, v1, A
   real(wp) :: u1, u2, u3
   real(wp) :: t1, t2, t3
   real(wp) :: A1, A2, A3
   real(wp) :: A2inv, A3inv, A4inv

   s  = 0.5_wp * svdwi
   B  = s * psi

   ! u(B) = aB - bB^2 + cB^3
   u    = B * (obc(1) + B * (obc(3) * B - obc(2)))
   up   = obc(1) - 2.0_wp * obc(2) * B + 3.0_wp * obc(3) * B * B          ! du/dB
   upp  = -2.0_wp * obc(2) + 6.0_wp * obc(3) * B                           ! d2u/dB2
   uppp = 6.0_wp * obc(3)                                                  ! d3u/dB3

   t  = tanh(u)
   ch = cosh(u)
   sech2 = 1.0_wp / (ch * ch)   ! sech^2(u)

   s1 = 1.0_wp / svdwi
   v1 = 1.0_wp / vdwri

   u1 = up   * s
   u2 = upp  * s * s
   u3 = uppp * s * s * s

   t1 = sech2 * u1
   t2 = sech2 * u2 - 2.0_wp * sech2 * t * (u1*u1)
   t3 = sech2 * ( u3 - 6.0_wp * t * u1 * u2 + (4.0_wp*t*t - 2.0_wp*sech2) * (u1*u1*u1) )

   A  = s1 - v1 * t
   A1 = -v1 * t1
   A2 = -v1 * t2
   A3 = -v1 * t3

   R = c1 / A

   A2inv = 1.0_wp / (A*A)
   A3inv = A2inv / A
   A4inv = A3inv / A

   fp   = -c1 * A1 * A2inv
   fpp  =  2.0_wp*c1*(A1*A1)*A3inv - c1*A2*A2inv
   fppp = -6.0_wp*c1*(A1*A1*A1)*A4inv + 6.0_wp*c1*(A1*A2)*A3inv - c1*A3*A2inv
end subroutine obc_map_d1d2d3


subroutine compute_bornr_d4(nat, vdwr, svdw, c1, obc, brad, brdr, brddr, brd3dr, brd4dr)
   integer, intent(in) :: nat
   real(wp), intent(in) :: vdwr(:), svdw(:), c1, obc(3)
   real(wp), intent(inout) :: brad(:)
   real(wp), intent(inout) :: brdr(:, :, :)
   real(wp), intent(inout) :: brddr(:, :, :, :, :)
   real(wp), intent(inout) :: brd3dr(:, :, :, :, :, :, :)
   real(wp), intent(inout) :: brd4dr(:, :, :, :, :, :, :, :, :)

   integer :: iat, k,l,m,n, a,b,c,d
   real(wp) :: R, fp, fpp, fppp, fpppp

   real(wp) :: g1,g2,g3,g4
   real(wp) :: H12,H13,H14,H23,H24,H34
   real(wp) :: T123,T124,T134,T234

   do iat = 1, nat
      call obc_map_d1d2d3d4(brad(iat), svdw(iat), vdwr(iat), c1, obc, R, fp, fpp, fppp, fpppp)

      ! ---- 4th derivative update first (needs old g/H/T/Q) ----
      do k = 1, nat
         do l = 1, nat
            do m = 1, nat
               do n = 1, nat
                  do a = 1, 3
                     g1 = brdr(a, k, iat)
                     do b = 1, 3
                        g2  = brdr(b, l, iat)
                        H12 = brddr(a, k, b, l, iat)

                        do c = 1, 3
                           g3  = brdr(c, m, iat)
                           H13 = brddr(a, k, c, m, iat)
                           H23 = brddr(b, l, c, m, iat)
                           T123 = brd3dr(a, k, b, l, c, m, iat)

                           do d = 1, 3
                              g4  = brdr(d, n, iat)
                              H14 = brddr(a, k, d, n, iat)
                              H24 = brddr(b, l, d, n, iat)
                              H34 = brddr(c, m, d, n, iat)

                              T124 = brd3dr(a, k, b, l, d, n, iat)
                              T134 = brd3dr(a, k, c, m, d, n, iat)
                              T234 = brd3dr(b, l, c, m, d, n, iat)

                              brd4dr(a, k, b, l, c, m, d, n, iat) = &
                                 fp    * brd4dr(a, k, b, l, c, m, d, n, iat) &
                               + fpp   * ( T123*g4 + T124*g3 + T134*g2 + T234*g1 &
                                           + H12*H34 + H13*H24 + H14*H23 ) &
                               + fppp  * ( H12*g3*g4 + H13*g2*g4 + H14*g2*g3 &
                                           + H23*g1*g4 + H24*g1*g3 + H34*g1*g2 ) &
                               + fpppp * ( g1*g2*g3*g4 )
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! ---- 3rd derivative update (same as your d3, but uses fp/fpp/fppp) ----
      do k = 1, nat
         do l = 1, nat
            do m = 1, nat
               do a = 1, 3
                  g1 = brdr(a, k, iat)
                  do b = 1, 3
                     g2  = brdr(b, l, iat)
                     H12 = brddr(a, k, b, l, iat)
                     do c = 1, 3
                        g3  = brdr(c, m, iat)
                        H13 = brddr(a, k, c, m, iat)
                        H23 = brddr(b, l, c, m, iat)

                        brd3dr(a, k, b, l, c, m, iat) = &
                           fp   * brd3dr(a, k, b, l, c, m, iat) &
                         + fpp  * (H12*g3 + H13*g2 + H23*g1) &
                         + fppp * (g1*g2*g3)
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! ---- Hessian update ----
      do k = 1, nat
         do l = 1, nat
            do a = 1, 3
               do b = 1, 3
                  brddr(a, k, b, l, iat) = fp * brddr(a, k, b, l, iat) &
                                        + fpp * brdr(a, k, iat) * brdr(b, l, iat)
               end do
            end do
         end do
      end do

      ! ---- Gradient + value update ----
      brdr(:, :, iat) = fp * brdr(:, :, iat)
      brad(iat) = R
   end do
end subroutine compute_bornr_d4


pure subroutine obc_map_d1d2d3d4(psi, svdwi, vdwri, c1, obc, R, fp, fpp, fppp, fpppp)
   real(wp), intent(in) :: psi, svdwi, vdwri, c1, obc(3)
   real(wp), intent(out) :: R, fp, fpp, fppp, fpppp

   real(wp) :: s, B, u, up, upp, uppp
   real(wp) :: t, ch, sech2
   real(wp) :: s1, v1, A
   real(wp) :: u1, u2, u3
   real(wp) :: t1, t2, t3, t4
   real(wp) :: A1, A2, A3, A4
   real(wp) :: A5inv

   s  = 0.5_wp * svdwi
   B  = s * psi

   ! u(B) = aB - bB^2 + cB^3
   u    = B * (obc(1) + B * (obc(3) * B - obc(2)))
   up   = obc(1) - 2.0_wp * obc(2) * B + 3.0_wp * obc(3) * B * B
   upp  = -2.0_wp * obc(2) + 6.0_wp * obc(3) * B
   uppp = 6.0_wp * obc(3)

   t  = tanh(u)
   ch = cosh(u)
   sech2 = 1.0_wp / (ch * ch)

   s1 = 1.0_wp / svdwi
   v1 = 1.0_wp / vdwri

   ! u derivatives wrt psi
   u1 = up   * s
   u2 = upp  * s * s
   u3 = uppp * s * s * s
   ! u4 = 0 for cubic u(B)

   ! tanh(u) derivatives wrt psi
   t1 = sech2 * u1
   t2 = sech2 * u2 - 2.0_wp * sech2 * t * (u1*u1)
   t3 = sech2 * ( u3 - 6.0_wp * t * u1 * u2 + (6.0_wp*t*t - 2.0_wp) * (u1*u1*u1) )
   t4 = sech2 * ( -8.0_wp * t * u1 * u3 - 6.0_wp * t * (u2*u2) &
                  + 12.0_wp * (2.0_wp - 3.0_wp*sech2) * (u1*u1*u2) &
                  + 8.0_wp  * t * (3.0_wp*sech2 - 1.0_wp) * (u1*u1*u1*u1) )

   A  = s1 - v1 * t
   A1 = -v1 * t1
   A2 = -v1 * t2
   A3 = -v1 * t3
   A4 = -v1 * t4

   R = c1 / A

   ! Use closed-form 4th derivative of c1/A in terms of A, A1..A4
   A5inv = 1.0_wp / (A*A*A*A*A)

   fp    = -A1*c1/(A*A)

   fpp   = c1 * (-A*A2 + 2.0_wp*A1*A1) / (A*A*A)

   fppp  = c1 * (-A*A*A3 + 6.0_wp*A*A1*A2 - 6.0_wp*A1*A1*A1) / (A*A*A*A)

   fpppp = c1 * ( -A*A*A*A4 + A*A*(8.0_wp*A1*A3 + 6.0_wp*A2*A2) &
                  - 36.0_wp*A*A1*A1*A2 + 24.0_wp*A1*A1*A1*A1 ) * A5inv
end subroutine obc_map_d1d2d3d4




pure subroutine compute_psi(nat, xyz, list, vdwr, rho, psi, dpsidr, d2psidr2, d3psidr3, d4psidr4)
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   type(adjacency_list), intent(in) :: list
   real(wp), intent(in) :: vdwr(:), rho(:)
   real(wp), intent(out) :: psi(:)
   real(wp), intent(out) :: dpsidr(:, :, :)
   real(wp), intent(out), optional :: d2psidr2(:, :, :, :, :)
   real(wp), intent(out), optional :: d3psidr3(:, :, :, :, :, :, :)
   real(wp), intent(out), optional :: d4psidr4(:, :, :, :, :, :, :, :, :)

   if (present(d4psidr4)) then
      call compute_psi_d4(nat, xyz, list, vdwr, rho, psi, dpsidr, d2psidr2, d3psidr3, d4psidr4)
   else if (present(d3psidr3)) then
      call compute_psi_d3(nat, xyz, list, vdwr, rho, psi, dpsidr, d2psidr2, d3psidr3)
   else if (present(d2psidr2)) then
      call compute_psi_d2(nat, xyz, list, vdwr, rho, psi, dpsidr, d2psidr2)
   else
      call compute_psi_d1(nat, xyz, list, vdwr, rho, psi, dpsidr)
   end if
end subroutine compute_psi




pure subroutine compute_psi_d1(nat, xyz, list, vdwr, rho, psi, dpsidr)
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Neighbourlist
   type(adjacency_list), intent(in) :: list
   !> Van-der-Waals radii
   real(wp), intent(in) :: vdwr(:)
   !> Descreened van-der-Waals radii
   real(wp), intent(in) :: rho(:)
   !> Integrated value of Psi
   real(wp), intent(out) :: psi(:)
   !> Derivative of Psi w.r.t. cartesian coordinates
   real(wp), intent(out) :: dpsidr(:, :, :)
   real(wp), allocatable :: dpsitr(:, :)

   integer  :: iat, jat, img, inl
   real(wp) :: vec(3), r, rhoi, rhoj
   real(wp) :: gi, gj, ap, am, lnab, rhab, ab, dgi, dgj
   real(wp) :: drjj(3)
   real(wp) :: rh1, rhr1, r24, r1, aprh1, r12
   real(wp) :: rvdwi, rvdwj
   logical :: ijov, jiov

   allocate(dpsitr(3, nat))
   psi(:) = 0.0_wp
   dpsidr(:, :, :) = 0.0_wp
   dpsitr(:, :) = 0.0_wp

   do iat = 1, nat
      inl = list%inl(iat)
      do img = 1, list%nnl(iat)
         jat = list%nlat(inl+img)

         vec(:) = xyz(:, iat) - xyz(:, jat)
         r = norm2(vec)

         rhoi = rho(iat)
         rhoj = rho(jat)
         rvdwi = vdwr(iat)
         rvdwj = vdwr(jat)

         ijov = r < (rvdwi+rhoj)
         jiov = r < (rhoi+rvdwj)

         if (.not.(ijov .or. jiov)) then
            ! nonoverlaping spheres
            if(abs(rhoi-rhoj) < 1.e-8_wp) then
               ! equal reduced radii
               r1 = 1.0_wp/r
               ap = r+rhoj
               am = r-rhoj
               ab = ap*am
               rhab = rhoj/ab
               lnab = 0.5_wp*log(am/ap)*r1
               gi = rhab+lnab
               dgi = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
               ! accumulate psi
               psi(iat) = psi(iat)+gi
               psi(jat) = psi(jat)+gi
               ! accumulate psi gradient
               drjj(:) = dgi*vec(:)
               dpsitr(:, iat) = dpsitr(:, iat)+drjj(:)
               dpsidr(:, jat, iat) = dpsidr(:, jat, iat)-drjj(:)
               dpsitr(:, jat) = dpsitr(:, jat)-drjj(:)
               dpsidr(:, iat, jat) = dpsidr(:, iat, jat)+drjj(:)
            else
               ! unequal reduced radii
               ! ij contribution
               r1 = 1.0_wp/r
               ap = r+rhoj
               am = r-rhoj
               ab = ap*am
               rhab = rhoj/ab
               lnab = 0.5_wp*log(am/ap)*r1
               gi = rhab+lnab
               dgi = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
               ! ji contribution
               ap = r+rhoi
               am = r-rhoi
               ab = ap*am
               rhab = rhoi/ab
               lnab = 0.5_wp*log(am/ap)*r1
               gj = rhab+lnab
               dgj = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
               ! accumulate psi
               psi(iat) = psi(iat)+gi
               psi(jat) = psi(jat)+gj
               ! accumulate psi gradient
               drjj(:) = dgi*vec(:)
               dpsitr(:, iat) = dpsitr(:, iat)+drjj(:)
               dpsidr(:, jat, iat) = dpsidr(:, jat, iat)-drjj(:)

               drjj(:) = dgj*vec(:)
               dpsitr(:, jat) = dpsitr(:, jat)-drjj(:)
               dpsidr(:, iat, jat) = dpsidr(:, iat, jat)+drjj(:)
            end if

         else if (.not.ijov .and. jiov) then

            ! ij contribution
            r1 = 1.0_wp/r
            ap = r+rhoj
            am = r-rhoj
            ab = ap*am
            rhab = rhoj/ab
            lnab = 0.5_wp*log(am/ap)*r1
            gi = rhab+lnab
            dgi = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
            ! accumulate psi
            psi(iat) = psi(iat)+gi
            ! accumulate psi gradient
            drjj(:) = dgi*vec(:)
            dpsitr(:, iat) = dpsitr(:, iat)+drjj(:)
            dpsidr(:, jat, iat) = dpsidr(:, jat, iat)-drjj(:)

            if((r+rhoi) > rvdwj) then
               ! ji contribution
               r1 = 1.0_wp/r
               r12 = 0.5_wp*r1
               r24 = r12*r12

               ap = r+rhoi
               am = r-rhoi
               rh1 = 1.0_wp/rvdwj
               rhr1 = 1.0_wp/ap
               aprh1 = ap*rh1
               lnab = log(aprh1)

               gj = rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

               dgj = rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
                  &         rhoi*r24*(rhr1-rh1*aprh1)+ &
                  &         r12*(r1*lnab-rhr1)
               dgj = dgj*r1
               ! accumulate psi
               psi(jat) = psi(jat)+gj
               ! accumulate psi gradient
               drjj(:) = dgj*vec(:)
               dpsitr(:, jat) = dpsitr(:, jat)-drjj(:)
               dpsidr(:, iat, jat) = dpsidr(:, iat, jat)+drjj(:)
            end if

         else if (ijov .and. .not.jiov) then

            if((r+rhoj) > rvdwi) then
               ! ij contribution
               r1 = 1.0_wp/r
               r12 = 0.5_wp*r1
               r24 = r12*r12

               ap = r+rhoj
               am = r-rhoj
               rh1 = 1.0_wp/rvdwi
               rhr1 = 1.0_wp/ap
               aprh1 = ap*rh1
               lnab = log(aprh1)

               gi = rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

               dgi = rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
                  &         rhoj*r24*(rhr1-rh1*aprh1)+ &
                  &         r12*(r1*lnab-rhr1)
               dgi = dgi*r1
               ! accumulate psi
               psi(iat) = psi(iat)+gi
               ! accumulate psi gradient
               drjj(:) = dgi*vec(:)
               dpsitr(:, iat) = dpsitr(:, iat)+drjj(:)
               dpsidr(:, jat, iat) = dpsidr(:, jat, iat)-drjj(:)
            end if

            ! ji contribution
            ap = r+rhoi
            am = r-rhoi
            ab = ap*am
            rhab = rhoi/ab
            lnab = 0.5_wp*log(am/ap)*r1
            gj = rhab+lnab
            dgj = -2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
            ! accumulate psi
            psi(jat) = psi(jat)+gj
            ! accumulate psi gradient
            drjj(:) = dgj*vec(:)
            dpsitr(:, jat) = dpsitr(:, jat)-drjj(:)
            dpsidr(:, iat, jat) = dpsidr(:, iat, jat)+drjj(:)

         else if (ijov .and. jiov) then
            ! overlaping spheres
            if((r+rhoj) > rvdwi) then
               ! ij contribution
               r1 = 1.0_wp/r
               r12 = 0.5_wp*r1
               r24 = r12*r12

               ap = r+rhoj
               am = r-rhoj
               rh1 = 1.0_wp/rvdwi
               rhr1 = 1.0_wp/ap
               aprh1 = ap*rh1
               lnab = log(aprh1)

               gi = rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

               dgi = rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
                  &         rhoj*r24*(rhr1-rh1*aprh1)+ &
                  &         r12*(r1*lnab-rhr1)
               dgi = dgi*r1
               ! accumulate psi
               psi(iat) = psi(iat)+gi
               ! accumulate psi gradient
               drjj(:) = dgi*vec(:)
               dpsitr(:, iat) = dpsitr(:, iat)+drjj(:)
               dpsidr(:, jat, iat) = dpsidr(:, jat, iat)-drjj(:)
            end if

            if((r+rhoi) > rvdwj) then
               ! ji contribution
               r1 = 1.0_wp/r
               r12 = 0.5_wp*r1
               r24 = r12*r12

               ap = r+rhoi
               am = r-rhoi
               rh1 = 1.0_wp/rvdwj
               rhr1 = 1.0_wp/ap
               aprh1 = ap*rh1
               lnab = log(aprh1)

               gj = rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

               dgj = rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
                  &         rhoi*r24*(rhr1-rh1*aprh1)+ &
                  &         r12*(r1*lnab-rhr1)
               dgj = dgj*r1
               ! accumulate psi
               psi(jat) = psi(jat)+gj
               ! accumulate psi gradient
               drjj(:) = dgj*vec(:)
               dpsitr(:, jat) = dpsitr(:, jat)-drjj(:)
               dpsidr(:, iat, jat) = dpsidr(:, iat, jat)+drjj(:)
            end if

         end if

      end do
   end do

   ! save one-center terms
   do iat = 1, nat
      dpsidr(:, iat, iat) = dpsitr(:, iat)
   end do
end subroutine compute_psi_d1


pure subroutine compute_psi_d2(nat, xyz, list, vdwr, rho, psi, dpsidr, d2psidr2)
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   type(adjacency_list), intent(in) :: list
   real(wp), intent(in) :: vdwr(:)
   real(wp), intent(in) :: rho(:)
   real(wp), intent(out) :: psi(:)
   real(wp), intent(out) :: dpsidr(:, :, :)
   real(wp), intent(out) :: d2psidr2(:, :, :, :, :)

   real(wp), allocatable :: dpsitr(:, :)
   real(wp), allocatable :: d2psitr(:, :, :)
   integer  :: iat, jat, img, inl
   real(wp) :: vec(3), r, rhoi, rhoj, rvdwi, rvdwj
   logical :: ijov, jiov

   type(ad2) :: rad
   type(ad2) :: gi, gj

   allocate(dpsitr(3, nat))
   allocate(d2psitr(3, 3, nat))

   psi(:) = 0.0_wp
   dpsidr(:, :, :) = 0.0_wp
   d2psidr2(:, :, :, :, :) = 0.0_wp
   dpsitr(:, :) = 0.0_wp
   d2psitr(:, :, :) = 0.0_wp

   do iat = 1, nat
      inl = list%inl(iat)
      do img = 1, list%nnl(iat)
         jat = list%nlat(inl+img)

         vec(:) = xyz(:, iat) - xyz(:, jat)
         r = norm2(vec)

         rhoi = rho(iat)
         rhoj = rho(jat)
         rvdwi = vdwr(iat)
         rvdwj = vdwr(jat)

         ijov = r < (rvdwi + rhoj)
         jiov = r < (rhoi + rvdwj)

         if (.not.(ijov .or. jiov)) then
            ! nonoverlapping
            rad = ad2_var(r)

            if (abs(rhoi - rhoj) < 1.e-8_wp) then
               gi = g_nonoverlap(rad, rhoj)
               call accum_pair(iat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr)
               call accum_pair(jat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr)
            else
               gi = g_nonoverlap(rad, rhoj)
               gj = g_nonoverlap(rad, rhoi)
               call accum_pair(iat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr)
               call accum_pair(jat, iat, jat, vec, r, gj, psi, dpsidr, dpsitr, d2psidr2, d2psitr)
            end if

         else if (.not.ijov .and. jiov) then
            ! i gets nonoverlap with rhoj
            rad = ad2_var(r)
            gi = g_nonoverlap(rad, rhoj)
            call accum_pair(iat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr)

            ! j may get overlap term
            if ((r + rhoi) > rvdwj) then
               gj = g_overlap(rad, rhoi, rvdwj)
               call accum_pair(jat, iat, jat, vec, r, gj, psi, dpsidr, dpsitr, d2psidr2, d2psitr)
            end if

         else if (ijov .and. .not.jiov) then
            rad = ad2_var(r)

            if ((r + rhoj) > rvdwi) then
               gi = g_overlap(rad, rhoj, rvdwi)
               call accum_pair(iat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr)
            end if

            gj = g_nonoverlap(rad, rhoi)
            call accum_pair(jat, iat, jat, vec, r, gj, psi, dpsidr, dpsitr, d2psidr2, d2psitr)

         else
            ! ijov .and. jiov
            rad = ad2_var(r)

            if ((r + rhoj) > rvdwi) then
               gi = g_overlap(rad, rhoj, rvdwi)
               call accum_pair(iat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr)
            end if

            if ((r + rhoi) > rvdwj) then
               gj = g_overlap(rad, rhoi, rvdwj)
               call accum_pair(jat, iat, jat, vec, r, gj, psi, dpsidr, dpsitr, d2psidr2, d2psitr)
            end if
         end if

      end do
   end do

   ! save one-center (diagonal wrt moved-atom == owner)
   do iat = 1, nat
      dpsidr(:, iat, iat) = dpsitr(:, iat)
      d2psidr2(:, iat, :, iat, iat) = d2psitr(:, :, iat)
   end do
end subroutine compute_psi_d2



pure subroutine compute_psi_d3(nat, xyz, list, vdwr, rho, psi, dpsidr, d2psidr2, d3psidr3)
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   type(adjacency_list), intent(in) :: list
   real(wp), intent(in) :: vdwr(:), rho(:)

   real(wp), intent(out) :: psi(:)
   real(wp), intent(out) :: dpsidr(:, :, :)
   real(wp), intent(out) :: d2psidr2(:, :, :, :, :)
   real(wp), intent(out) :: d3psidr3(:, :, :, :, :, :, :)

   real(wp), allocatable :: dpsitr(:, :)
   real(wp), allocatable :: d2psitr(:, :, :)
   real(wp), allocatable :: d3psitr(:, :, :, :)

   integer  :: iat, jat, img, inl
   real(wp) :: vec(3), r, rhoi, rhoj, rvdwi, rvdwj
   logical :: ijov, jiov

   type(ad3) :: rad
   type(ad3) :: gi, gj

   allocate(dpsitr(3, nat))
   allocate(d2psitr(3, 3, nat))
   allocate(d3psitr(3, 3, 3, nat))

   psi(:) = 0.0_wp
   dpsidr(:, :, :) = 0.0_wp
   d2psidr2(:, :, :, :, :) = 0.0_wp
   d3psidr3(:, :, :, :, :, :, :) = 0.0_wp

   dpsitr(:, :) = 0.0_wp
   d2psitr(:, :, :) = 0.0_wp
   d3psitr(:, :, :, :) = 0.0_wp

   do iat = 1, nat
      inl = list%inl(iat)
      do img = 1, list%nnl(iat)
         jat = list%nlat(inl+img)

         vec(:) = xyz(:, iat) - xyz(:, jat)
         r = norm2(vec)

         rhoi = rho(iat)
         rhoj = rho(jat)
         rvdwi = vdwr(iat)
         rvdwj = vdwr(jat)

         ijov = r < (rvdwi + rhoj)
         jiov = r < (rhoi + rvdwj)

         rad = ad3_var(r)

         if (.not.(ijov .or. jiov)) then
            ! nonoverlapping spheres
            if (abs(rhoi - rhoj) < 1.e-8_wp) then
               gi = g_nonoverlap3(rad, rhoj)
               call accum_pair3(iat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr)
               call accum_pair3(jat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr)
            else
               gi = g_nonoverlap3(rad, rhoj)
               gj = g_nonoverlap3(rad, rhoi)
               call accum_pair3(iat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr)
               call accum_pair3(jat, iat, jat, vec, r, gj, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr)
            end if

         else if (.not.ijov .and. jiov) then
            ! i gets nonoverlap with rhoj
            gi = g_nonoverlap3(rad, rhoj)
            call accum_pair3(iat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr)

            ! j may get overlap term
            if ((r + rhoi) > rvdwj) then
               gj = g_overlap3(rad, rhoi, rvdwj)
               call accum_pair3(jat, iat, jat, vec, r, gj, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr)
            end if

         else if (ijov .and. .not.jiov) then
            if ((r + rhoj) > rvdwi) then
               gi = g_overlap3(rad, rhoj, rvdwi)
               call accum_pair3(iat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr)
            end if

            gj = g_nonoverlap3(rad, rhoi)
            call accum_pair3(jat, iat, jat, vec, r, gj, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr)

         else
            ! ijov .and. jiov (overlap for both, if allowed)
            if ((r + rhoj) > rvdwi) then
               gi = g_overlap3(rad, rhoj, rvdwi)
               call accum_pair3(iat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr)
            end if

            if ((r + rhoi) > rvdwj) then
               gj = g_overlap3(rad, rhoi, rvdwj)
               call accum_pair3(jat, iat, jat, vec, r, gj, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr)
            end if
         end if
      end do
   end do

   ! Save one-center terms (fully diagonal blocks)
   do iat = 1, nat
      dpsidr(:, iat, iat) = dpsitr(:, iat)
      d2psidr2(:, iat, :, iat, iat) = d2psitr(:, :, iat)
      d3psidr3(:, iat, :, iat, :, iat, iat) = d3psitr(:, :, :, iat)
   end do
end subroutine compute_psi_d3



pure subroutine compute_psi_d4(nat, xyz, list, vdwr, rho, psi, dpsidr, d2psidr2, d3psidr3, d4psidr4)
   integer, intent(in) :: nat
   real(wp), intent(in) :: xyz(:, :)
   type(adjacency_list), intent(in) :: list
   real(wp), intent(in) :: vdwr(:), rho(:)

   real(wp), intent(out) :: psi(:)
   real(wp), intent(out) :: dpsidr(:, :, :)
   real(wp), intent(out) :: d2psidr2(:, :, :, :, :)
   real(wp), intent(out) :: d3psidr3(:, :, :, :, :, :, :)
   real(wp), intent(out) :: d4psidr4(:, :, :, :, :, :, :, :, :)

   real(wp), allocatable :: dpsitr(:, :)
   real(wp), allocatable :: d2psitr(:, :, :)
   real(wp), allocatable :: d3psitr(:, :, :, :)
   real(wp), allocatable :: d4psitr(:, :, :, :, :)

   integer  :: iat, jat, img, inl
   real(wp) :: vec(3), r, rhoi, rhoj, rvdwi, rvdwj
   logical :: ijov, jiov

   type(ad4) :: rad
   type(ad4) :: gi, gj

   allocate(dpsitr(3, nat))
   allocate(d2psitr(3, 3, nat))
   allocate(d3psitr(3, 3, 3, nat))
   allocate(d4psitr(3, 3, 3, 3, nat))

   psi(:) = 0.0_wp
   dpsidr(:, :, :) = 0.0_wp
   d2psidr2(:, :, :, :, :) = 0.0_wp
   d3psidr3(:, :, :, :, :, :, :) = 0.0_wp
   d4psidr4(:, :, :, :, :, :, :, :, :) = 0.0_wp

   dpsitr(:, :) = 0.0_wp
   d2psitr(:, :, :) = 0.0_wp
   d3psitr(:, :, :, :) = 0.0_wp
   d4psitr(:, :, :, :, :) = 0.0_wp

   do iat = 1, nat
      inl = list%inl(iat)
      do img = 1, list%nnl(iat)
         jat = list%nlat(inl+img)

         vec(:) = xyz(:, iat) - xyz(:, jat)
         r = norm2(vec)

         rhoi = rho(iat)
         rhoj = rho(jat)
         rvdwi = vdwr(iat)
         rvdwj = vdwr(jat)

         ijov = r < (rvdwi + rhoj)
         jiov = r < (rhoi + rvdwj)

         rad = ad4_var(r)

         if (.not.(ijov .or. jiov)) then
            if (abs(rhoi - rhoj) < 1.e-8_wp) then
               gi = g_nonoverlap4(rad, rhoj)
               call accum_pair4(iat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr, d4psidr4, d4psitr)
               call accum_pair4(jat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr, d4psidr4, d4psitr)
            else
               gi = g_nonoverlap4(rad, rhoj)
               gj = g_nonoverlap4(rad, rhoi)
               call accum_pair4(iat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr, d4psidr4, d4psitr)
               call accum_pair4(jat, iat, jat, vec, r, gj, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr, d4psidr4, d4psitr)
            end if

         else if (.not.ijov .and. jiov) then
            gi = g_nonoverlap4(rad, rhoj)
            call accum_pair4(iat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr, d4psidr4, d4psitr)

            if ((r + rhoi) > rvdwj) then
               gj = g_overlap4(rad, rhoi, rvdwj)
               call accum_pair4(jat, iat, jat, vec, r, gj, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr, d4psidr4, d4psitr)
            end if

         else if (ijov .and. .not.jiov) then
            if ((r + rhoj) > rvdwi) then
               gi = g_overlap4(rad, rhoj, rvdwi)
               call accum_pair4(iat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr, d4psidr4, d4psitr)
            end if

            gj = g_nonoverlap4(rad, rhoi)
            call accum_pair4(jat, iat, jat, vec, r, gj, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr, d4psidr4, d4psitr)

         else
            if ((r + rhoj) > rvdwi) then
               gi = g_overlap4(rad, rhoj, rvdwi)
               call accum_pair4(iat, iat, jat, vec, r, gi, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr, d4psidr4, d4psitr)
            end if

            if ((r + rhoi) > rvdwj) then
               gj = g_overlap4(rad, rhoi, rvdwj)
               call accum_pair4(jat, iat, jat, vec, r, gj, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr, d4psidr4, d4psitr)
            end if
         end if
      end do
   end do

   do iat = 1, nat
      dpsidr(:, iat, iat) = dpsitr(:, iat)
      d2psidr2(:, iat, :, iat, iat) = d2psitr(:, :, iat)
      d3psidr3(:, iat, :, iat, :, iat, iat) = d3psitr(:, :, :, iat)
      d4psidr4(:, iat, :, iat, :, iat, :, iat, iat) = d4psitr(:, :, :, :, iat)
   end do
end subroutine compute_psi_d4







! -------- 2nd-order scalar AD type (w.r.t. r) --------


pure elemental function ad2_c(x) result(a)
   real(wp), intent(in) :: x
   type(ad2) :: a
   a%v = x; a%d1 = 0.0_wp; a%d2 = 0.0_wp
end function ad2_c

pure elemental function ad2_var(x) result(a)
   real(wp), intent(in) :: x
   type(ad2) :: a
   a%v = x; a%d1 = 1.0_wp; a%d2 = 0.0_wp
end function ad2_var

pure elemental function add_ad2(a,b) result(c)
   type(ad2), intent(in) :: a,b
   type(ad2) :: c
   c%v = a%v + b%v
   c%d1 = a%d1 + b%d1
   c%d2 = a%d2 + b%d2
end function

pure elemental function sub_ad2(a,b) result(c)
   type(ad2), intent(in) :: a,b
   type(ad2) :: c
   c%v = a%v - b%v
   c%d1 = a%d1 - b%d1
   c%d2 = a%d2 - b%d2
end function

pure elemental function neg_ad2(a) result(c)
   type(ad2), intent(in) :: a
   type(ad2) :: c
   c%v = -a%v
   c%d1 = -a%d1
   c%d2 = -a%d2
end function

pure elemental function mul_ad2(a,b) result(c)
   type(ad2), intent(in) :: a,b
   type(ad2) :: c
   c%v  = a%v*b%v
   c%d1 = a%d1*b%v + a%v*b%d1
   c%d2 = a%d2*b%v + 2.0_wp*a%d1*b%d1 + a%v*b%d2
end function

pure elemental function div_ad2(a,b) result(c)
   type(ad2), intent(in) :: a,b
   type(ad2) :: c
   real(wp) :: bv2, bv3
   bv2 = b%v*b%v
   bv3 = bv2*b%v
   c%v  = a%v / b%v
   c%d1 = (a%d1*b%v - a%v*b%d1) / bv2
   c%d2 = (a%d2*b%v - a%v*b%d2) / bv2 - 2.0_wp*(a%d1*b%v - a%v*b%d1)*b%d1 / bv3
end function

pure elemental function log_ad2(a) result(c)
   type(ad2), intent(in) :: a
   type(ad2) :: c
   real(wp) :: inv
   inv = 1.0_wp/a%v
   c%v  = log(a%v)
   c%d1 = a%d1*inv
   c%d2 = (a%d2*a%v - a%d1*a%d1) * (inv*inv)
end function


pure elemental function ad3_c(x) result(a)
   real(wp), intent(in) :: x
   type(ad3) :: a
   a%v = x; a%d1 = 0.0_wp; a%d2 = 0.0_wp; a%d3 = 0.0_wp
end function ad3_c

pure elemental function ad3_var(x) result(a)
   real(wp), intent(in) :: x
   type(ad3) :: a
   a%v = x; a%d1 = 1.0_wp; a%d2 = 0.0_wp; a%d3 = 0.0_wp
end function ad3_var

pure elemental function add_ad3(a,b) result(c)
   type(ad3), intent(in) :: a,b
   type(ad3) :: c
   c%v = a%v + b%v
   c%d1 = a%d1 + b%d1
   c%d2 = a%d2 + b%d2
   c%d3 = a%d3 + b%d3
end function add_ad3

pure elemental function sub_ad3(a,b) result(c)
   type(ad3), intent(in) :: a,b
   type(ad3) :: c
   c%v = a%v - b%v
   c%d1 = a%d1 - b%d1
   c%d2 = a%d2 - b%d2
   c%d3 = a%d3 - b%d3
end function sub_ad3

pure elemental function neg_ad3(a) result(c)
   type(ad3), intent(in) :: a
   type(ad3) :: c
   c%v = -a%v
   c%d1 = -a%d1
   c%d2 = -a%d2
   c%d3 = -a%d3
end function neg_ad3

pure elemental function mul_ad3(a,b) result(c)
   type(ad3), intent(in) :: a,b
   type(ad3) :: c
   c%v  = a%v*b%v
   c%d1 = a%d1*b%v + a%v*b%d1
   c%d2 = a%d2*b%v + 2.0_wp*a%d1*b%d1 + a%v*b%d2
   c%d3 = a%d3*b%v + 3.0_wp*a%d2*b%d1 + 3.0_wp*a%d1*b%d2 + a%v*b%d3
end function mul_ad3

pure elemental function inv_ad3(b) result(r)
   type(ad3), intent(in) :: b
   type(ad3) :: r
   real(wp) :: b0, b1, b2, b3, b02, b03, b04
   b0 = b%v; b1 = b%d1; b2 = b%d2; b3 = b%d3
   b02 = b0*b0
   b03 = b02*b0
   b04 = b03*b0

   r%v  = 1.0_wp/b0
   r%d1 = -b1/b02
   r%d2 = (2.0_wp*b1*b1 - b0*b2)/b03
   r%d3 = (-6.0_wp*b1*b1*b1 + 6.0_wp*b0*b1*b2 - b0*b0*b3)/b04
end function inv_ad3

pure elemental function div_ad3(a,b) result(c)
   type(ad3), intent(in) :: a,b
   type(ad3) :: c
   c = a * inv_ad3(b)
end function div_ad3

pure elemental function ad3_log(a) result(c)
   type(ad3), intent(in) :: a
   type(ad3) :: c
   real(wp) :: x, x2, x3
   x  = a%v
   x2 = x*x
   x3 = x2*x
   c%v  = log(x)                          ! intrinsic log(real)
   c%d1 = a%d1/x
   c%d2 = (a%d2*x - a%d1*a%d1)/x2
   c%d3 = (a%d3*x2 - 3.0_wp*a%d2*x*a%d1 + 2.0_wp*a%d1*a%d1*a%d1)/x3
end function ad3_log


pure elemental function ad4_c(x) result(a)
   real(wp), intent(in) :: x
   type(ad4) :: a
   a%v = x; a%d1 = 0.0_wp; a%d2 = 0.0_wp; a%d3 = 0.0_wp; a%d4 = 0.0_wp
end function ad4_c

pure elemental function ad4_var(x) result(a)
   real(wp), intent(in) :: x
   type(ad4) :: a
   a%v = x; a%d1 = 1.0_wp; a%d2 = 0.0_wp; a%d3 = 0.0_wp; a%d4 = 0.0_wp
end function ad4_var

pure elemental function add_ad4(a,b) result(c)
   type(ad4), intent(in) :: a,b
   type(ad4) :: c
   c%v  = a%v  + b%v
   c%d1 = a%d1 + b%d1
   c%d2 = a%d2 + b%d2
   c%d3 = a%d3 + b%d3
   c%d4 = a%d4 + b%d4
end function add_ad4

pure elemental function sub_ad4(a,b) result(c)
   type(ad4), intent(in) :: a,b
   type(ad4) :: c
   c%v  = a%v  - b%v
   c%d1 = a%d1 - b%d1
   c%d2 = a%d2 - b%d2
   c%d3 = a%d3 - b%d3
   c%d4 = a%d4 - b%d4
end function sub_ad4

pure elemental function neg_ad4(a) result(c)
   type(ad4), intent(in) :: a
   type(ad4) :: c
   c%v  = -a%v
   c%d1 = -a%d1
   c%d2 = -a%d2
   c%d3 = -a%d3
   c%d4 = -a%d4
end function neg_ad4

pure elemental function mul_ad4(a,b) result(c)
   type(ad4), intent(in) :: a,b
   type(ad4) :: c
   c%v  = a%v*b%v
   c%d1 = a%d1*b%v + a%v*b%d1
   c%d2 = a%d2*b%v + 2.0_wp*a%d1*b%d1 + a%v*b%d2
   c%d3 = a%d3*b%v + 3.0_wp*a%d2*b%d1 + 3.0_wp*a%d1*b%d2 + a%v*b%d3
   c%d4 = a%d4*b%v + 4.0_wp*a%d3*b%d1 + 6.0_wp*a%d2*b%d2 + 4.0_wp*a%d1*b%d3 + a%v*b%d4
end function mul_ad4

pure elemental function inv_ad4(b) result(r)
   type(ad4), intent(in) :: b
   type(ad4) :: r
   real(wp) :: b0,b1,b2,b3,b4
   real(wp) :: b02,b03,b04,b05

   b0 = b%v;  b1 = b%d1; b2 = b%d2; b3 = b%d3; b4 = b%d4
   b02 = b0*b0
   b03 = b02*b0
   b04 = b03*b0
   b05 = b04*b0

   r%v  = 1.0_wp/b0
   r%d1 = -b1/b02
   r%d2 = (2.0_wp*b1*b1 - b0*b2)/b03
   r%d3 = (-6.0_wp*b1*b1*b1 + 6.0_wp*b0*b1*b2 - b0*b0*b3)/b04
   r%d4 = (-b0*b0*b0*b4 + b0*b0*(8.0_wp*b1*b3 + 6.0_wp*b2*b2) &
           - 36.0_wp*b0*b1*b1*b2 + 24.0_wp*b1*b1*b1*b1) / b05
end function inv_ad4

pure elemental function div_ad4(a,b) result(c)
   type(ad4), intent(in) :: a,b
   type(ad4) :: c
   c = a * inv_ad4(b)
end function div_ad4

pure elemental function ad4_log(a) result(c)
   type(ad4), intent(in) :: a
   type(ad4) :: c
   real(wp) :: x, x2, x3, x4
   x  = a%v
   x2 = x*x
   x3 = x2*x
   x4 = x3*x

   c%v  = log(x)
   c%d1 = a%d1/x
   c%d2 = (a%d2*x - a%d1*a%d1)/x2
   c%d3 = (a%d3*x2 - 3.0_wp*a%d2*x*a%d1 + 2.0_wp*a%d1*a%d1*a%d1)/x3
   c%d4 = (a%d4*x3 - x2*(4.0_wp*a%d1*a%d3 + 3.0_wp*a%d2*a%d2) &
           + 12.0_wp*x*a%d1*a%d1*a%d2 - 6.0_wp*a%d1*a%d1*a%d1*a%d1) / x4
end function ad4_log




! -------- Pair contribution functions (same formulas, AD-enabled) --------

pure elemental function g_nonoverlap(r, rho_s) result(g)
   type(ad2), intent(in) :: r
   real(wp), intent(in) :: rho_s
   type(ad2) :: g
   type(ad2) :: ap, am, ab, r1, lnab, rhab

   ap = r + ad2_c(rho_s)
   am = r - ad2_c(rho_s)
   ab = ap * am
   rhab = ad2_c(rho_s) / ab
   r1 = ad2_c(1.0_wp) / r
   lnab = ad2_c(0.5_wp) * log_ad2(am / ap) * r1
   g = rhab + lnab
end function g_nonoverlap

pure elemental function g_overlap(r, rho_s, rvdw_t) result(g)
   type(ad2), intent(in) :: r
   real(wp), intent(in) :: rho_s, rvdw_t
   type(ad2) :: g
   type(ad2) :: ap, am, r1, r12, rhr1, aprh1, lnab
   real(wp) :: rh1

   rh1 = 1.0_wp/rvdw_t
   r1  = ad2_c(1.0_wp) / r
   r12 = ad2_c(0.5_wp) * r1

   ap = r + ad2_c(rho_s)
   am = r - ad2_c(rho_s)

   rhr1 = ad2_c(1.0_wp) / ap
   aprh1 = ap * ad2_c(rh1)
   lnab = log_ad2(aprh1)

   g = ad2_c(rh1) - rhr1 + r12 * ( ad2_c(0.5_wp) * am * (rhr1 - ad2_c(rh1)*aprh1) - lnab )
end function g_overlap


pure elemental function g_nonoverlap3(r, rho_s) result(g)
   type(ad3), intent(in) :: r
   real(wp), intent(in) :: rho_s
   type(ad3) :: g
   type(ad3) :: ap, am, ab, r1, lnab, rhab

   ap = r + ad3_c(rho_s)
   am = r - ad3_c(rho_s)
   ab = ap * am
   rhab = ad3_c(rho_s) / ab
   r1 = ad3_c(1.0_wp) / r
   lnab = ad3_c(0.5_wp) * ad3_log(am / ap) * r1
   g = rhab + lnab
end function g_nonoverlap3

pure elemental function g_overlap3(r, rho_s, rvdw_t) result(g)
   type(ad3), intent(in) :: r
   real(wp), intent(in) :: rho_s, rvdw_t
   type(ad3) :: g
   type(ad3) :: ap, am, r1, r12, rhr1, aprh1, lnab
   real(wp) :: rh1

   rh1 = 1.0_wp/rvdw_t
   r1  = ad3_c(1.0_wp) / r
   r12 = ad3_c(0.5_wp) * r1

   ap = r + ad3_c(rho_s)
   am = r - ad3_c(rho_s)

   rhr1 = ad3_c(1.0_wp) / ap
   aprh1 = ap * ad3_c(rh1)
   lnab = ad3_log(aprh1)

   g = ad3_c(rh1) - rhr1 + r12 * ( ad3_c(0.5_wp) * am * (rhr1 - ad3_c(rh1)*aprh1) - lnab )
end function g_overlap3

pure elemental function g_nonoverlap4(r, rho_s) result(g)
   type(ad4), intent(in) :: r
   real(wp), intent(in) :: rho_s
   type(ad4) :: g
   type(ad4) :: ap, am, ab, r1, lnab, rhab

   ap   = r + ad4_c(rho_s)
   am   = r - ad4_c(rho_s)
   ab   = ap * am
   rhab = ad4_c(rho_s) / ab
   r1   = ad4_c(1.0_wp) / r
   lnab = ad4_c(0.5_wp) * ad4_log(am / ap) * r1
   g    = rhab + lnab
end function g_nonoverlap4

pure elemental function g_overlap4(r, rho_s, rvdw_t) result(g)
   type(ad4), intent(in) :: r
   real(wp), intent(in) :: rho_s, rvdw_t
   type(ad4) :: g
   type(ad4) :: ap, am, r1, r12, rhr1, aprh1, lnab
   real(wp) :: rh1

   rh1 = 1.0_wp/rvdw_t
   r1  = ad4_c(1.0_wp) / r
   r12 = ad4_c(0.5_wp) * r1

   ap = r + ad4_c(rho_s)
   am = r - ad4_c(rho_s)

   rhr1  = ad4_c(1.0_wp) / ap
   aprh1 = ap * ad4_c(rh1)
   lnab  = ad4_log(aprh1)

   g = ad4_c(rh1) - rhr1 + r12 * ( ad4_c(0.5_wp) * am * (rhr1 - ad4_c(rh1)*aprh1) - lnab )
end function g_overlap4



! -------- Accumulation of one scalar g(r) into psi(owner), grad, Hess --------
pure subroutine accum_pair(owner, p, q, vec, r, g, psi, dpsidr, dpsitr, d2psidr2, d2psitr)
   integer, intent(in) :: owner, p, q
   real(wp), intent(in) :: vec(3), r
   type(ad2), intent(in) :: g
   real(wp), intent(inout) :: psi(:)
   real(wp), intent(inout) :: dpsidr(:, :, :)
   real(wp), intent(inout) :: dpsitr(:, :)
   real(wp), intent(inout) :: d2psidr2(:, :, :, :, :)
   real(wp), intent(inout) :: d2psitr(:, :, :)

   real(wp) :: gp, gpp
   real(wp) :: g1, dd
   real(wp) :: dr(3)
   real(wp) :: Hv(3,3)
   integer :: a,b

   ! scalar derivatives wrt r
   gp  = g%d1
   gpp = g%d2

   ! gradient scalar factor: g'(r)/r
   g1 = gp / r

   ! Hessian scalar factor for outer term: (g'' r - g') / r^3
   dd = (gpp*r - gp) / (r*r*r)

   psi(owner) = psi(owner) + g%v

   dr(:) = g1 * vec(:)

   ! grad wrt p: +dr ; wrt q: -dr
   if (p == owner) then
      dpsitr(:, owner) = dpsitr(:, owner) + dr(:)
   else
      dpsidr(:, p, owner) = dpsidr(:, p, owner) + dr(:)
   end if

   if (q == owner) then
      dpsitr(:, owner) = dpsitr(:, owner) - dr(:)
   else
      dpsidr(:, q, owner) = dpsidr(:, q, owner) - dr(:)
   end if

   ! Build Hv = g1 I + dd * (vec vec^T)
   Hv(:, :) = 0.0_wp
   do a = 1,3
      Hv(a,a) = g1
   end do
   do a = 1,3
      do b = 1,3
         Hv(a,b) = Hv(a,b) + dd * vec(a) * vec(b)
      end do
   end do

   call add_block(owner, p, p, +1.0_wp, Hv, d2psidr2, d2psitr)
   call add_block(owner, p, q, -1.0_wp, Hv, d2psidr2, d2psitr)
   call add_block(owner, q, p, -1.0_wp, Hv, d2psidr2, d2psitr)
   call add_block(owner, q, q, +1.0_wp, Hv, d2psidr2, d2psitr)
end subroutine accum_pair

pure subroutine add_block(owner, aidx, bidx, sgn, Hv, d2psidr2, d2psitr)
   integer, intent(in) :: owner, aidx, bidx
   real(wp), intent(in) :: sgn
   real(wp), intent(in) :: Hv(3,3)
   real(wp), intent(inout) :: d2psidr2(:, :, :, :, :)
   real(wp), intent(inout) :: d2psitr(:, :, :)

   if (aidx == owner .and. bidx == owner) then
      d2psitr(:, :, owner) = d2psitr(:, :, owner) + sgn * Hv(:, :)
   else
      d2psidr2(:, aidx, :, bidx, owner) = d2psidr2(:, aidx, :, bidx, owner) + sgn * Hv(:, :)
   end if
end subroutine add_block



pure subroutine accum_pair3(owner, p, q, vec, r, g, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr)
   integer, intent(in) :: owner, p, q
   real(wp), intent(in) :: vec(3), r
   type(ad3), intent(in) :: g
   real(wp), intent(inout) :: psi(:)
   real(wp), intent(inout) :: dpsidr(:, :, :)
   real(wp), intent(inout) :: dpsitr(:, :)
   real(wp), intent(inout) :: d2psidr2(:, :, :, :, :)
   real(wp), intent(inout) :: d2psitr(:, :, :)
   real(wp), intent(inout) :: d3psidr3(:, :, :, :, :, :, :)
   real(wp), intent(inout) :: d3psitr(:, :, :, :)

   real(wp) :: gp, gpp, gppp
   real(wp) :: g1, dd, beta
   real(wp) :: dr(3)
   real(wp) :: Hv(3,3)
   real(wp) :: Tv(3,3,3)
   integer :: a,b,c
   integer :: i1,i2,i3
   integer :: idx(2)
   real(wp) :: sgn(2)
   real(wp) :: s

   gp   = g%d1
   gpp  = g%d2
   gppp = g%d3

   ! gradient factor: g'(r)/r
   g1 = gp / r

   ! Hessian outer factor: (g'' r - g') / r^3
   dd = (gpp*r - gp) / (r*r*r)

   ! Third-derivative outer factor:
   ! beta = g'''/r^3 - 3 g''/r^4 + 3 g'/r^5
   beta = gppp/(r*r*r) - 3.0_wp*gpp/(r*r*r*r) + 3.0_wp*gp/(r*r*r*r*r)

   psi(owner) = psi(owner) + g%v

   ! ---- Gradient accumulation ----
   dr(:) = g1 * vec(:)

   if (p == owner) then
      dpsitr(:, owner) = dpsitr(:, owner) + dr(:)
   else
      dpsidr(:, p, owner) = dpsidr(:, p, owner) + dr(:)
   end if

   if (q == owner) then
      dpsitr(:, owner) = dpsitr(:, owner) - dr(:)
   else
      dpsidr(:, q, owner) = dpsidr(:, q, owner) - dr(:)
   end if

   ! ---- Hessian (same as your AD2 version) ----
   Hv(:, :) = 0.0_wp
   do a = 1,3
      Hv(a,a) = g1
   end do
   do a = 1,3
      do b = 1,3
         Hv(a,b) = Hv(a,b) + dd * vec(a) * vec(b)
      end do
   end do

   call add_block2(owner, p, p, +1.0_wp, Hv, d2psidr2, d2psitr)
   call add_block2(owner, p, q, -1.0_wp, Hv, d2psidr2, d2psitr)
   call add_block2(owner, q, p, -1.0_wp, Hv, d2psidr2, d2psitr)
   call add_block2(owner, q, q, +1.0_wp, Hv, d2psidr2, d2psitr)

   ! ---- Third derivative tensor w.r.t. v = r_p - r_q ----
   ! T_abc = dd*(δ_ab v_c + δ_ac v_b + δ_bc v_a) + beta*v_a v_b v_c
   Tv(:, :, :) = 0.0_wp
   do a = 1,3
      do b = 1,3
         do c = 1,3
            Tv(a,b,c) = Tv(a,b,c) + beta * vec(a)*vec(b)*vec(c)
            if (a == b) Tv(a,b,c) = Tv(a,b,c) + dd * vec(c)
            if (a == c) Tv(a,b,c) = Tv(a,b,c) + dd * vec(b)
            if (b == c) Tv(a,b,c) = Tv(a,b,c) + dd * vec(a)
         end do
      end do
   end do

   ! sign bookkeeping: derivative wrt p -> +, wrt q -> -
   idx(1) = p; sgn(1) = +1.0_wp
   idx(2) = q; sgn(2) = -1.0_wp

   do i1 = 1,2
      do i2 = 1,2
         do i3 = 1,2
            s = sgn(i1)*sgn(i2)*sgn(i3)
            call add_block3(owner, idx(i1), idx(i2), idx(i3), s, Tv, d3psidr3, d3psitr)
         end do
      end do
   end do
end subroutine accum_pair3


pure subroutine add_block2(owner, aidx, bidx, sgn, Hv, d2psidr2, d2psitr)
   integer, intent(in) :: owner, aidx, bidx
   real(wp), intent(in) :: sgn
   real(wp), intent(in) :: Hv(3,3)
   real(wp), intent(inout) :: d2psidr2(:, :, :, :, :)
   real(wp), intent(inout) :: d2psitr(:, :, :)

   if (aidx == owner .and. bidx == owner) then
      d2psitr(:, :, owner) = d2psitr(:, :, owner) + sgn * Hv(:, :)
   else
      d2psidr2(:, aidx, :, bidx, owner) = d2psidr2(:, aidx, :, bidx, owner) + sgn * Hv(:, :)
   end if
end subroutine add_block2


pure subroutine add_block3(owner, aidx, bidx, cidx, sgn, Tv, d3psidr3, d3psitr)
   integer, intent(in) :: owner, aidx, bidx, cidx
   real(wp), intent(in) :: sgn
   real(wp), intent(in) :: Tv(3,3,3)
   real(wp), intent(inout) :: d3psidr3(:, :, :, :, :, :, :)
   real(wp), intent(inout) :: d3psitr(:, :, :, :)

   if (aidx == owner .and. bidx == owner .and. cidx == owner) then
      d3psitr(:, :, :, owner) = d3psitr(:, :, :, owner) + sgn * Tv(:, :, :)
   else
      d3psidr3(:, aidx, :, bidx, :, cidx, owner) = d3psidr3(:, aidx, :, bidx, :, cidx, owner) + sgn * Tv(:, :, :)
   end if
end subroutine add_block3


pure subroutine add_block4(owner, aidx, bidx, cidx, didx, sgn, Qv, d4psidr4, d4psitr)
   integer, intent(in) :: owner, aidx, bidx, cidx, didx
   real(wp), intent(in) :: sgn
   real(wp), intent(in) :: Qv(3,3,3,3)
   real(wp), intent(inout) :: d4psidr4(:, :, :, :, :, :, :, :, :)
   real(wp), intent(inout) :: d4psitr(:, :, :, :, :)

   if (aidx == owner .and. bidx == owner .and. cidx == owner .and. didx == owner) then
      d4psitr(:, :, :, :, owner) = d4psitr(:, :, :, :, owner) + sgn * Qv(:, :, :, :)
   else
      d4psidr4(:, aidx, :, bidx, :, cidx, :, didx, owner) = &
         d4psidr4(:, aidx, :, bidx, :, cidx, :, didx, owner) + sgn * Qv(:, :, :, :)
   end if
end subroutine add_block4


pure subroutine accum_pair4(owner, p, q, vec, r, g, psi, dpsidr, dpsitr, d2psidr2, d2psitr, d3psidr3, d3psitr, d4psidr4, d4psitr)
   integer, intent(in) :: owner, p, q
   real(wp), intent(in) :: vec(3), r
   type(ad4), intent(in) :: g

   real(wp), intent(inout) :: psi(:)
   real(wp), intent(inout) :: dpsidr(:, :, :)
   real(wp), intent(inout) :: dpsitr(:, :)
   real(wp), intent(inout) :: d2psidr2(:, :, :, :, :)
   real(wp), intent(inout) :: d2psitr(:, :, :)
   real(wp), intent(inout) :: d3psidr3(:, :, :, :, :, :, :)
   real(wp), intent(inout) :: d3psitr(:, :, :, :)
   real(wp), intent(inout) :: d4psidr4(:, :, :, :, :, :, :, :, :)
   real(wp), intent(inout) :: d4psitr(:, :, :, :, :)

   real(wp) :: gp, gpp, gppp, gpppp
   real(wp) :: g1, dd, beta, gamma
   real(wp) :: dr(3)
   real(wp) :: Hv(3,3)
   real(wp) :: Tv(3,3,3)
   real(wp) :: Qv(3,3,3,3)
   integer :: a,b,c,d
   integer :: i1,i2,i3,i4
   integer :: idx(2)
   real(wp) :: sgn(2)
   real(wp) :: s

   gp    = g%d1
   gpp   = g%d2
   gppp  = g%d3
   gpppp = g%d4

   g1 = gp / r
   dd = (gpp*r - gp) / (r*r*r)
   beta  = gppp/(r*r*r) - 3.0_wp*gpp/(r*r*r*r) + 3.0_wp*gp/(r*r*r*r*r)
   gamma = gpppp/(r*r*r*r) - 6.0_wp*gppp/(r*r*r*r*r) + 15.0_wp*gpp/(r*r*r*r*r*r) - 15.0_wp*gp/(r*r*r*r*r*r*r)

   psi(owner) = psi(owner) + g%v

   ! ---- Gradient ----
   dr(:) = g1 * vec(:)

   if (p == owner) then
      dpsitr(:, owner) = dpsitr(:, owner) + dr(:)
   else
      dpsidr(:, p, owner) = dpsidr(:, p, owner) + dr(:)
   end if
   if (q == owner) then
      dpsitr(:, owner) = dpsitr(:, owner) - dr(:)
   else
      dpsidr(:, q, owner) = dpsidr(:, q, owner) - dr(:)
   end if

   ! ---- Hessian ----
   Hv(:, :) = 0.0_wp
   do a = 1,3
      Hv(a,a) = g1
   end do
   do a = 1,3
      do b = 1,3
         Hv(a,b) = Hv(a,b) + dd * vec(a) * vec(b)
      end do
   end do
   call add_block2(owner, p, p, +1.0_wp, Hv, d2psidr2, d2psitr)
   call add_block2(owner, p, q, -1.0_wp, Hv, d2psidr2, d2psitr)
   call add_block2(owner, q, p, -1.0_wp, Hv, d2psidr2, d2psitr)
   call add_block2(owner, q, q, +1.0_wp, Hv, d2psidr2, d2psitr)

   ! ---- 3rd tensor ----
   Tv(:, :, :) = 0.0_wp
   do a = 1,3
      do b = 1,3
         do c = 1,3
            Tv(a,b,c) = Tv(a,b,c) + beta * vec(a)*vec(b)*vec(c)
            if (a == b) Tv(a,b,c) = Tv(a,b,c) + dd * vec(c)
            if (a == c) Tv(a,b,c) = Tv(a,b,c) + dd * vec(b)
            if (b == c) Tv(a,b,c) = Tv(a,b,c) + dd * vec(a)
         end do
      end do
   end do
   call add_block3(owner, p, p, p, +1.0_wp, Tv, d3psidr3, d3psitr)
   call add_block3(owner, p, p, q, -1.0_wp, Tv, d3psidr3, d3psitr)
   call add_block3(owner, p, q, p, -1.0_wp, Tv, d3psidr3, d3psitr)
   call add_block3(owner, q, p, p, -1.0_wp, Tv, d3psidr3, d3psitr)
   call add_block3(owner, p, q, q, +1.0_wp, Tv, d3psidr3, d3psitr)
   call add_block3(owner, q, p, q, +1.0_wp, Tv, d3psidr3, d3psitr)
   call add_block3(owner, q, q, p, +1.0_wp, Tv, d3psidr3, d3psitr)
   call add_block3(owner, q, q, q, -1.0_wp, Tv, d3psidr3, d3psitr)

   ! ---- 4th tensor (isotropic radial formula) ----
   ! Q_abcd = dd*(δ_ab δ_cd + δ_ac δ_bd + δ_bc δ_ad)
   !        + beta*(δ_ab v_c v_d + δ_ac v_b v_d + δ_bc v_a v_d + δ_ad v_b v_c + δ_bd v_a v_c + δ_cd v_a v_b)
   !        + gamma*v_a v_b v_c v_d
   Qv(:, :, :, :) = 0.0_wp
   do a = 1,3
      do b = 1,3
         do c = 1,3
            do d = 1,3
               Qv(a,b,c,d) = Qv(a,b,c,d) + gamma * vec(a)*vec(b)*vec(c)*vec(d)

               if (a == b .and. c == d) Qv(a,b,c,d) = Qv(a,b,c,d) + dd
               if (a == c .and. b == d) Qv(a,b,c,d) = Qv(a,b,c,d) + dd
               if (b == c .and. a == d) Qv(a,b,c,d) = Qv(a,b,c,d) + dd

               if (a == b) Qv(a,b,c,d) = Qv(a,b,c,d) + beta * vec(c)*vec(d)
               if (a == c) Qv(a,b,c,d) = Qv(a,b,c,d) + beta * vec(b)*vec(d)
               if (b == c) Qv(a,b,c,d) = Qv(a,b,c,d) + beta * vec(a)*vec(d)
               if (a == d) Qv(a,b,c,d) = Qv(a,b,c,d) + beta * vec(b)*vec(c)
               if (b == d) Qv(a,b,c,d) = Qv(a,b,c,d) + beta * vec(a)*vec(c)
               if (c == d) Qv(a,b,c,d) = Qv(a,b,c,d) + beta * vec(a)*vec(b)
            end do
         end do
      end do
   end do

   idx(1) = p; sgn(1) = +1.0_wp
   idx(2) = q; sgn(2) = -1.0_wp

   do i1 = 1,2
      do i2 = 1,2
         do i3 = 1,2
            do i4 = 1,2
               s = sgn(i1)*sgn(i2)*sgn(i3)*sgn(i4)
               call add_block4(owner, idx(i1), idx(i2), idx(i3), idx(i4), s, Qv, d4psidr4, d4psitr)
            end do
         end do
      end do
   end do
end subroutine accum_pair4





end module tblite_solvation_born
