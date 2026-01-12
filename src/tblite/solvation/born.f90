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


interface operator(+); module procedure add_ad2; end interface
interface operator(-); module procedure sub_ad2, neg_ad2; end interface
interface operator(*); module procedure mul_ad2; end interface
interface operator(/); module procedure div_ad2; end interface

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

   interface log
   module procedure log_ad2
   end interface

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
subroutine get_rad(self, mol, rad, draddr, dradd2r)
   !> Instance of the Born integrator
   class(born_integrator), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Born radii
   real(wp), intent(out) :: rad(:)
   !> Derivative of Born radii w.r.t. cartesian displacements
   real(wp), intent(out), optional :: draddr(:, :, :)
   !> Second derivative of Born radii w.r.t. cartesian displacements
   !> Layout: (3, nat, 3, nat, nat_owner)
   real(wp), intent(out), optional :: dradd2r(:, :, :, :, :)

   type(adjacency_list) :: list
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), allocatable :: brdr(:, :, :)
   real(wp), allocatable :: brddr(:, :, :, :, :)
   

   call new_adjacency_list(list, mol, trans, self%lrcut)

   allocate(brdr(3, mol%nat, mol%nat))
   if (present(dradd2r)) allocate(brddr(3, mol%nat, 3, mol%nat, mol%nat))

   if (present(dradd2r)) then
      call compute_bornr(mol%nat, mol%xyz, list, &
         & self%vdwr, self%rho, self%svdw, self%born_scale, self%obc, rad, brdr, brddr)
   else
      call compute_bornr(mol%nat, mol%xyz, list, &
         & self%vdwr, self%rho, self%svdw, self%born_scale, self%obc, rad, brdr)
   end if

   if (present(draddr)) then
      draddr(:, :, :) = brdr
   end if
   if (present(dradd2r)) then
      dradd2r(:, :, :, :, :) = brddr
   end if
end subroutine get_rad


subroutine compute_bornr(nat, xyz, list, vdwr, rho, svdw, c1, obc, &
      & brad, brdr, brddr)
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
   !> van-der-Waals radii with offset
   real(wp), intent(in) :: svdw(:)
   !> Scaling factor for the Born radii
   real(wp), intent(in) :: c1
   !> Volume polynome correction
   real(wp), intent(in) :: obc(3)
   !> Born radii
   real(wp), intent(out) :: brad(:)
   !> Derivative of Born radii w.r.t. cartesian coordinates
   real(wp), intent(out) :: brdr(:, :, :)
   !> Second derivative of Born radii w.r.t. cartesian coordinates
   real(wp), intent(out), optional :: brddr(:, :, :, :, :)

   integer :: iat
   real(wp) :: br, dpsi, svdwi, vdwri, s1, v1, s2, arg, arg2
   real(wp) :: th, ch
   real(wp) :: R, fp, fpp

   call compute_psi(nat, xyz, list, vdwr, rho, brad, brdr)

   do iat = 1, nat

      br = brad(iat)

      svdwi = svdw(iat)
      vdwri = vdwr(iat)
      s1 = 1.0_wp/svdwi
      v1 = 1.0_wp/vdwri
      s2 = 0.5_wp*svdwi

      br = br*s2

      arg2 = br*(obc(3)*br-obc(2))
      arg = br*(obc(1)+arg2)
      arg2 = 2.0_wp*arg2+obc(1)+obc(3)*br*br

      th = tanh(arg)
      ch = cosh(arg)

      br = 1.0_wp/(s1-v1*th)
      ! Include GBMV2-like scaling
      br = c1*br

      dpsi = ch*(s1-v1*th)
      dpsi = s2*v1*arg2/(dpsi*dpsi)
      dpsi = c1*dpsi

      brad(iat) = br
      brdr(:, :, iat) = brdr(:, :, iat) * dpsi

   end do

   ! Can probably be optimized later
   if (present(brddr)) then
      call compute_psi(nat, xyz, list, vdwr, rho, brad, brdr, brddr)
      call compute_bornr_d2(nat, vdwr, svdw, c1, obc, brad, brdr, brddr)
   else
      call compute_psi(nat, xyz, list, vdwr, rho, brad, brdr)
      do iat = 1, nat
         call obc_map_d1d2(brad(iat), svdw(iat), vdwr(iat), c1, obc, R, fp, fpp)
         brad(iat) = R
         brdr(:, :, iat) = brdr(:, :, iat) * fp
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


pure subroutine compute_psi(nat, xyz, list, vdwr, rho, psi, dpsidr, d2psidr2)
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

   real(wp), intent(out), optional :: d2psidr2(:, :, :, :, :)

    ! Probably can be optimized later
   if (present(d2psidr2)) then
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
   lnab = ad2_c(0.5_wp) * log(am / ap) * r1
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
   lnab = log(aprh1)

   g = ad2_c(rh1) - rhr1 + r12 * ( ad2_c(0.5_wp) * am * (rhr1 - ad2_c(rh1)*aprh1) - lnab )
end function g_overlap

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



end module tblite_solvation_born
