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

!> @file tblite/solvation/radii_scaling.f90
!> Provides radii scaling models for implicit solvation cavities.

!> Radii scaling models for implicit solvation.
module tblite_solvation_radii_scaling
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use multicharge, only : get_eeq_charges
   use mctc_ncoord, only : new_ncoord, ncoord_type, cn_count
   use mctc_ncoord_type, only : ncoord_type
   use tblite_solvation_data_draco, only : get_alpha, min_rad, &
      & eeq_prefac_water_cpcm, eeq_expo_water_cpcm, eeq_k_water_cpcm, &
      & eeq_prefac_other_solvents_cpcm, eeq_expo_other_solvents_cpcm, &
      & eeq_k_other_solvents_cpcm, eeq_o_shift_other_solvents_cpcm, &
      & eeq_prefac_water_cosmo, eeq_expo_water_cosmo, eeq_k_water_cosmo, &
      & eeq_prefac_other_solvents_cosmo, eeq_expo_other_solvents_cosmo, &
      & eeq_k_other_solvents_cosmo, eeq_o_shift_other_solvents_cosmo, &
      & eeq_prefac_water_smd, eeq_expo_water_smd, eeq_k_water_smd, &
      & eeq_prefac_other_solvents_smd, eeq_expo_other_solvents_smd, &
      & eeq_k_other_solvents_smd, eeq_o_shift_other_solvents_smd
   implicit none
   private

   public :: radii_scaling_type, draco_radii_scaling, new_draco_radii_scaling
   public :: draco

   !> Abstract base class for solvation cavity radii scaling models.
   type, public, abstract :: radii_scaling_type
   contains
      procedure(scale_radii), deferred :: scale
   end type radii_scaling_type

   abstract interface
      subroutine scale_radii(self, mol, radii_in, radii_out, q, solvent, cn, &
            & atoms_to_scale, dqdr, dcndr, drdr)
         import :: radii_scaling_type, structure_type, wp
         class(radii_scaling_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         real(wp), intent(in) :: radii_in(:)
         real(wp), intent(out) :: radii_out(:)
         real(wp), intent(in) :: q(:)
         character(len=*), intent(in) :: solvent
         real(wp), intent(in), optional :: cn(:)
         integer, intent(in), optional :: atoms_to_scale(:)
         real(wp), intent(in), optional, contiguous :: dqdr(:, :, :)
         real(wp), intent(in), optional, contiguous :: dcndr(:, :, :)
         real(wp), intent(out), optional, contiguous :: drdr(:, :, :)
      end subroutine scale_radii
   end interface

   !> DRACO radii scaling parameters.
   type, public, extends(radii_scaling_type) :: draco_radii_scaling
      !> Radii family the parameters are applied to.
      character(len=:), allocatable :: radtype
      !> Prevent radii from becoming smaller than the DRACO fit domain.
      logical :: damp_small_rad = .false.
      !> Scaling parameters.
      real(wp) :: prefac(94) = 0.0_wp
      real(wp) :: expo(94) = 0.0_wp
      real(wp) :: o_shift(94) = 0.0_wp
      real(wp) :: k1(94) = 0.0_wp
   contains
      procedure :: scale => scale_draco_radii
   end type draco_radii_scaling

contains

!> Create a DRACO radii scaling model with bundled tblite parameters.
subroutine new_draco_radii_scaling(self, radtype, solvent, damp_small_rad)
   !> DRACO radii scaling model.
   type(draco_radii_scaling), intent(out) :: self
   !> Radii family, one of cpcm, cosmo, or smd.
   character(len=*), intent(in) :: radtype
   !> Solvent for parameter selection.
   character(len=*), intent(in) :: solvent
   !> Prevent radii from becoming smaller than the DRACO fit domain.
   logical, intent(in), optional :: damp_small_rad

   self%radtype = radtype
   if (present(damp_small_rad)) self%damp_small_rad = damp_small_rad

   select case(radtype)
   case("cpcm")
      select case(solvent)
      case("water")
         self%prefac = eeq_prefac_water_cpcm
         self%expo = eeq_expo_water_cpcm
         self%k1 = eeq_k_water_cpcm
         self%o_shift = 0.0_wp
      case default
         self%prefac = eeq_prefac_other_solvents_cpcm
         self%expo = eeq_expo_other_solvents_cpcm
         self%k1 = eeq_k_other_solvents_cpcm
         self%o_shift = eeq_o_shift_other_solvents_cpcm
      end select
   case("cosmo")
      select case(solvent)
      case("water")
         self%prefac = eeq_prefac_water_cosmo
         self%expo = eeq_expo_water_cosmo
         self%k1 = eeq_k_water_cosmo
         self%o_shift = 0.0_wp
      case default
         self%prefac = eeq_prefac_other_solvents_cosmo
         self%expo = eeq_expo_other_solvents_cosmo
         self%k1 = eeq_k_other_solvents_cosmo
         self%o_shift = eeq_o_shift_other_solvents_cosmo
      end select
   case("smd")
      select case(solvent)
      case("water")
         self%prefac = eeq_prefac_water_smd
         self%expo = eeq_expo_water_smd
         self%k1 = eeq_k_water_smd
         self%o_shift = 0.0_wp
      case default
         self%prefac = eeq_prefac_other_solvents_smd
         self%expo = eeq_expo_other_solvents_smd
         self%k1 = eeq_k_other_solvents_smd
         self%o_shift = eeq_o_shift_other_solvents_smd
      end select
   end select
end subroutine new_draco_radii_scaling

!> Convenience interface for applying DRACO directly from solvation models.
subroutine draco(mol, radii_in, radii_out, solvent, radtype, q, cn, atoms_to_scale, &
      & damp_small_rad, dqdr, dcndr, drdr)
   !> Molecular structure data.
   type(structure_type), intent(in) :: mol
   !> Input radii in bohr.
   real(wp), intent(in) :: radii_in(:)
   !> Scaled radii in bohr.
   real(wp), intent(out) :: radii_out(:)
   !> Solvent for parameter selection.
   character(len=*), intent(in) :: solvent
   !> Radii family, one of cpcm, cosmo, or smd.
   character(len=*), intent(in), optional :: radtype
   !> Atomic charges. If absent, EEQ charges are used.
   real(wp), intent(in), optional :: q(:)
   !> Coordination numbers. If absent, GFN coordination numbers are used.
   real(wp), intent(in), optional :: cn(:)
   !> Atomic numbers whose radii should be scaled.
   integer, intent(in), optional :: atoms_to_scale(:)
   !> Prevent radii from becoming smaller than the DRACO fit domain.
   logical, intent(in), optional :: damp_small_rad
   !> Charge derivatives w.r.t. Cartesian coordinates.
   real(wp), intent(in), optional, contiguous :: dqdr(:, :, :)
   !> Coordination-number derivatives w.r.t. Cartesian coordinates.
   real(wp), intent(in), optional, contiguous :: dcndr(:, :, :)
   !> Radii derivatives w.r.t. Cartesian coordinates.
   real(wp), intent(out), optional, contiguous :: drdr(:, :, :)

   type(draco_radii_scaling) :: scale
   character(len=:), allocatable :: local_radtype
   real(wp), allocatable :: qscratch(:), cnscratch(:)
   real(wp), allocatable :: dqdr_scratch(:, :, :), dcndr_scratch(:, :, :)
   class(ncoord_type), allocatable :: ncoord

   type(error_type), allocatable :: error

   local_radtype = "cosmo"
   if (present(radtype)) local_radtype = radtype

   call new_draco_radii_scaling(scale, local_radtype, solvent, damp_small_rad)

   if (present(q)) then
      qscratch = q
   else
      allocate(qscratch(mol%nat))
      if (present(drdr)) then
         call get_eeq_charges(mol, error, qscratch, dqdr=dqdr_scratch)
      else
         call get_eeq_charges(mol, error, qscratch)
      end if
   end if

   if (present(cn)) then
      cnscratch = cn
   else
      allocate(cnscratch(mol%nat))
      call new_ncoord(ncoord, mol, cn_count%erf_en, error)
      if (present(drdr)) then
         allocate(dcndr_scratch(3, mol%nat, mol%nat))
         call ncoord%get_cn(mol, cnscratch, dcndr_scratch)
      else
         call ncoord%get_cn(mol, cnscratch)
      end if
   end if

   if (present(drdr)) then
      if (present(dqdr) .and. present(dcndr)) then
         call scale%scale(mol, radii_in, radii_out, qscratch, solvent, cnscratch, &
            & atoms_to_scale, dqdr, dcndr, drdr)
      else if (allocated(dqdr_scratch) .and. allocated(dcndr_scratch)) then
         call scale%scale(mol, radii_in, radii_out, qscratch, solvent, cnscratch, &
            & atoms_to_scale, dqdr_scratch, dcndr_scratch, drdr)
      else
         call scale%scale(mol, radii_in, radii_out, qscratch, solvent, cnscratch, &
            & atoms_to_scale)
         drdr = 0.0_wp
      end if
   else
      call scale%scale(mol, radii_in, radii_out, qscratch, solvent, cnscratch, &
         & atoms_to_scale)
   end if
end subroutine draco

!> Apply DRACO radii scaling.
subroutine scale_draco_radii(self, mol, radii_in, radii_out, q, solvent, cn, &
      & atoms_to_scale, dqdr, dcndr, drdr)
   !> DRACO radii scaling model.
   class(draco_radii_scaling), intent(in) :: self
   !> Molecular structure data.
   type(structure_type), intent(in) :: mol
   !> Input radii in bohr.
   real(wp), intent(in) :: radii_in(:)
   !> Scaled radii in bohr.
   real(wp), intent(out) :: radii_out(:)
   !> Atomic charges.
   real(wp), intent(in) :: q(:)
   !> Solvent for parameter selection.
   character(len=*), intent(in) :: solvent
   !> Coordination numbers.
   real(wp), intent(in), optional :: cn(:)
   !> Atomic numbers whose radii should be scaled.
   integer, intent(in), optional :: atoms_to_scale(:)
   !> Charge derivatives w.r.t. Cartesian coordinates.
   real(wp), intent(in), optional, contiguous :: dqdr(:, :, :)
   !> Coordination-number derivatives w.r.t. Cartesian coordinates.
   real(wp), intent(in), optional, contiguous :: dcndr(:, :, :)
   !> Radii derivatives w.r.t. Cartesian coordinates.
   real(wp), intent(out), optional, contiguous :: drdr(:, :, :)

   real(wp) :: arg, cn_iat, outer, scale, shift
   integer :: iat, izp
   logical :: grad

   grad = present(dqdr) .and. present(dcndr) .and. present(drdr)
   if (grad) drdr = 0.0_wp

   do iat = 1, mol%nat
      izp = mol%num(mol%id(iat))
      if (scale_atom(izp, atoms_to_scale)) then
         cn_iat = 0.0_wp
         if (present(cn)) cn_iat = cn(iat)
         arg = q(iat) + self%k1(izp)*q(iat)*cn_iat - self%expo(izp)
         scale = erf(self%prefac(izp)*arg) + 1.0_wp
         if (self%damp_small_rad .and. scale < min_rad(self%radtype)) then
            scale = min_rad(self%radtype)
         end if
         shift = 0.0_wp
         if (izp == 8 .and. get_alpha(solvent) < 0.43_wp) then
            shift = self%o_shift(izp)*(0.43_wp - get_alpha(solvent))
         end if
         radii_out(iat) = (scale + shift) * radii_in(iat)
         if (grad) then
            outer = 2.0_wp * exp(-self%prefac(izp)**2 * arg**2) / sqrt(pi)
            drdr(:, :, iat) = self%prefac(izp) * &
               & (dqdr(:, :, iat) + self%k1(izp)*dqdr(:, :, iat)*cn_iat &
               & + q(iat)*self%k1(izp)*dcndr(:, :, iat)) * outer * radii_in(iat)
         end if
      else
         radii_out(iat) = radii_in(iat)
      end if
   end do
end subroutine scale_draco_radii

!> Check whether an element should be scaled.
pure function scale_atom(izp, atoms_to_scale) result(scale)
   !> Atomic number.
   integer, intent(in) :: izp
   !> Atomic numbers whose radii should be scaled.
   integer, intent(in), optional :: atoms_to_scale(:)
   !> Whether this atom should be scaled.
   logical :: scale

   if (present(atoms_to_scale)) then
      scale = any(atoms_to_scale == izp)
   else
      scale = izp > 0 .and. izp <= 94
   end if
end function scale_atom

end module tblite_solvation_radii_scaling

