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

!> @file tblite/solvation/draco/type.f90
!> Provides DRACO implementation

!> DRACO radii scaling
module tblite_solvation_draco_type
   use mctc_env, only: wp
   use mctc_io_convert, only: autoaa, aatoau
   use mctc_io_constants, only: pi
   use mctc_io, only: structure_type
   use mctc_io, only: toNumber => to_number
   use tblite_solvation_data_draco, only: get_alpha_cpcm, get_shift_cpcm, get_kcn_cpcm
   use tblite_solvation_data_draco, only: get_alpha_smd, get_shift_smd, get_kcn_smd
   use tblite_solvation_data_draco, only: get_alpha_cosmo, get_shift_cosmo, get_kcn_cosmo
   use tblite_solvation_data
   implicit none
   private
   public :: draco_type, new_draco

   !> Draco type
   type :: draco_type
      !> Default radii from the solvation model
      real(wp), allocatable :: defaultradii(:)

      !> Shift parameter for the error function (adaps scaled radii relative to the original values)
      real(wp), allocatable :: shift(:)

      !> Sensitivity of atomic radius toward change in enveronment (paramter a in error function)
      real(wp), allocatable :: alpha(:)

      !> CN dependence of the effective charge (differentiates between various bonding motifs)
      real(wp), allocatable :: kcn(:)

      !> Cutoff radius
      real(wp) :: minrad



   contains
      procedure :: radii_adjustment
      procedure :: get_fscale
   end type draco_type

contains
   subroutine new_draco(self, mol, radtype) ! Gets the default radii based on solvation model
      !> DRACO trype
      class(draco_type), intent(out) :: self

      !> Solvation model that defines the intital radii
      character(len=*), intent(in) :: radtype

      !> Molecuar structure type
      type(structure_type), intent(in) :: mol


      select case (radtype)
      case ('cpcm')
          self%defaultradii = get_vdw_rad_cpcm(mol%num)
          self%kcn = get_kcn_cpcm(mol%num)
          self%alpha = get_alpha_cpcm(mol%num)
          self%shift = get_shift_cpcm(mol%num)

          self%minrad = 0.4_wp
      case ('smd')
          self%defaultradii = get_vdw_rad_smd(mol%num)
          self%kcn = get_kcn_smd(mol%num)
          self%alpha = get_alpha_smd(mol%num)
          self%shift = get_shift_smd(mol%num)

          self%minrad = 0.5_wp
      case ('cosmo')
          self%defaultradii = get_vdw_rad_cosmo(mol%num)
          self%kcn = get_kcn_cosmo(mol%num)
          self%alpha = get_alpha_cosmo(mol%num)
          self%shift = get_shift_cosmo(mol%num)

          self%minrad = 0.4_wp
      end select
   end subroutine new_draco


   !> Subroutine that rescales the radii, for now implemented only for water as a solvent
   subroutine radii_adjustment(self, mol, q, cn, scaled_radii, dsraddr, dsraddL, dqdr, dcndr, dqdL, dcndL)
      !> DRACO class
      class(draco_type), intent(in) :: self

      !> Molecular structure data
      type(structure_type), intent(in) :: mol

      !> Atomic partial charges
      real(wp), intent(in) :: q(:)

      !> Coordination numbers
      real(wp), intent(in) :: cn(:)

      !> Radii, adjusted with scaling function
      real(wp), allocatable, intent(out) :: scaled_radii(:)

      !> Derivative of scaled radii w.r.t. atom positions
      real(wp), intent(out), optional :: dsraddr(:, :, :)

      !> Derivative of scaled radii w.r.t. lattice vectors
      real(wp), intent(out), optional :: dsraddL(:, :, :)

      !> Derivative of partial charges w.r.t. atom positions
      real(wp), intent(in), optional :: dqdr(:, :, :)

      !> Derivative of coordination numbers w.r.t. atom positions
      real(wp), intent(in), optional ::  dcndr(:, :, :)

      !> Derivative of partial charges w.r.t. lattice vectors
      real(wp), intent(in), optional :: dqdL(:, :, :)

      !> Derivative of coordination numbers w.r.t. lattice vectors
      real(wp), intent(in), optional ::  dcndL(:, :, :)
      
      
      ! Scaling function
      real(wp), allocatable :: fscale(:) 

      !> Derivative of scaling function w.r.t. atom positions
      real(wp), allocatable :: dfscaledr(:, :, :)

      !> Derivative of scaling function w.r.t. lattice vectors
      real(wp), allocatable :: dfscaledL(:, :, :)

      integer :: iat, izp
      logical :: grad

      allocate (fscale(mol%nat), scaled_radii(mol%nat))

      grad = present(dsraddr) .and. present(dqdr) .and. present(dcndr) &
         & .and. present(dsraddL) .and. present(dqdL) .and. present(dcndL)
      if (grad) then
         allocate(dfscaledr(3, mol%nat, mol%nat))
         allocate(dfscaledL(3, 3, mol%nat))
         call self%get_fscale(mol, q, cn, fscale, dqdr, dcndr, dqdL, dcndL, dfscaledr, dfscaledL)

         ! Scale the radii
         do iat = 1, mol%nat
             izp = mol%id(iat)
             scaled_radii(iat) = self%defaultradii(izp) * fscale(iat) * autoaa
             
             ! Get the derivative
             dsraddr(:, :, iat) = self%defaultradii(izp) * dfscaledr(:, :, iat) * autoaa
             dsraddL(:, :, iat) = self%defaultradii(izp) * dfscaledL(:, :, iat) * autoaa
         end do
      else
          call self%get_fscale(mol, q, cn, fscale)

          ! Scale the radii
          do iat = 1, mol%nat
              izp = mol%id(iat)
              scaled_radii(iat) = self%defaultradii(izp) * fscale(iat) * autoaa

              if (scaled_radii(iat) <= self%minrad) scaled_radii(iat) = self%minrad
              
          end do
      end if
   end subroutine radii_adjustment


   !> Subroutine that gets the scaling function for the radii and optionally its derivative
   subroutine get_fscale(self, mol, q, cn, fscale, dqdr, dcndr, dqdL, dcndL, dfscaledr, dfscaledL) 
      !> DRACO type
      class(draco_type), intent(in) :: self

      !> Molecular structure data
      type(structure_type), intent(in) :: mol

      !> Partial charges
      real(wp), intent(in) :: q(:)

      !> Coordination numbers
      real(wp), intent(in) :: cn(:)

      !> Scaling function
      real(wp), intent(out) :: fscale(:)
  
      !> Derivative of partial charges w.r.t. atom positions
      real(wp), intent(in), optional :: dqdr(:, :, :)

      !> Derivative of coordination numbers w.r.t. atom positions
      real(wp), intent(in), optional :: dcndr(:, :, :)

      !> Derivative of partial charges w.r.t. lattice vectors
      real(wp), intent(in), optional :: dqdL(:, :, :)

      !> Derivative of coordination numbers w.r.t. lattice vectors
      real(wp), intent(in), optional :: dcndL(:, :, :)

      !> Derivative of scaling function w.r.t. atom positions
      real(wp), intent(out), optional :: dfscaledr(:, :, :)

      !> Derivative of scaling function w.r.t. lattice vectors
      real(wp), intent(out), optional :: dfscaledL(:, :, :)

      !> Effecvtive partial charges
      real(wp) :: qeff

      !> Derivative of effective partial charges w.r.t. atom positions
      real(wp), allocatable :: dqeffdr(:, :, :)

      !> Derivative of effective partial charges w.r.t. lattice vectors
      real(wp), allocatable ::dqeffdL(:, :, :)

      !> Partial derivative of the scaling function w.r.t. effectrive partial charge
      real(wp) :: dfscaledqeff

      integer :: iat, izp, isp
      logical :: grad
      
      grad = present(dqdr) .and. present(dcndr) .and. present(dfscaledr) .and. &
         & present(dqdL) .and. present(dcndL) .and. present(dfscaledL)
      if (grad) then
         allocate(dqeffdr(3, mol%nat, mol%nat))
         allocate(dqeffdL(3, 3, mol%nat))
         do iat = 1, mol%nat
           izp = mol%id(iat)

           qeff = q(iat) + self%kcn(izp) * q(iat) * cn(iat)
           fscale(iat) = erf(self%alpha(izp) * (qeff - self%shift(izp))) + 1

           dfscaledqeff = 2/sqrt(pi) * self%alpha(izp) * exp(-self%alpha(izp)**2 * (qeff - self%shift(izp))**2)
           
           dqeffdr(:, :, iat) = (self%kcn(izp) * cn(iat) + 1) * dqdr(:, :, iat) + self%kcn(izp) * q(iat) * dcndr(:, :, iat)
           dfscaledr(:, :, iat) = dfscaledqeff * dqeffdr(:, :, iat)

           dqeffdL(:, :, iat) = (self%kcn(izp) * cn(iat) + 1) * dqdL(:, :, iat) + self%kcn(izp) * q(iat) * dcndL(:, :, iat)
           dfscaledL(:, :, iat) = dfscaledqeff * dqeffdL(:, :, iat)
         end do
      else
         do iat = 1, mol%nat
           izp = mol%id(iat)

           qeff = q(iat) + self%kcn(izp) * q(iat) * cn(iat)
           fscale(iat) = erf(self%alpha(izp) * (qeff - self%shift(izp))) + 1
         end do
      end if
   end subroutine get_fscale

end module tblite_solvation_draco_type
