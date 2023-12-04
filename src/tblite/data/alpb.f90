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

!> @file tblite/data/alpb.f90
!> Provides ALPB/GBSA parameters

!> ALPB/GBSA parameters
module tblite_data_alpb
   use mctc_env, only : error_type, fatal_error
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   use tblite_solvation_alpb, only: alpb_input
   use mctc_io, only : structure_type
   use mctc_io_convert, only : aatoau
   implicit none
   private

   public :: get_alpb_param

   type :: alpb_parameter
      real(wp) :: epsv = 0.0_wp
      real(wp) :: smass = 0.0_wp
      real(wp) :: rhos = 0.0_wp
      real(wp) :: c1 = 0.0_wp
      real(wp) :: rprobe = 0.0_wp
      real(wp) :: gshift = 0.0_wp
      real(wp) :: soset = 0.0_wp
      real(wp) :: alpha = 0.0_wp
      real(wp) :: gamscale(94) = 0.0_wp
      real(wp) :: sx(94) = 0.0_wp
      real(wp) :: tmp(94) = 0.0_wp
   end type alpb_parameter

   include 'alpb/param_gbsa_acetone.fh'
   include 'alpb/param_gbsa_acetonitrile.fh'
   include 'alpb/param_gbsa_benzene.fh'
   include 'alpb/param_gbsa_ch2cl2.fh'
   include 'alpb/param_gbsa_chcl3.fh'
   include 'alpb/param_gbsa_cs2.fh'
   include 'alpb/param_gbsa_dmso.fh'
   include 'alpb/param_gbsa_ether.fh'
   include 'alpb/param_gbsa_h2o.fh'
   include 'alpb/param_gbsa_methanol.fh'
   include 'alpb/param_gbsa_thf.fh'
   include 'alpb/param_gbsa_toluene.fh'
   include 'alpb/param_gbsa_dmf.fh'
   include 'alpb/param_gbsa_nhexan.fh'

   include 'alpb/param_alpb_acetone.fh'
   include 'alpb/param_alpb_acetonitrile.fh'
   include 'alpb/param_alpb_aniline.fh'
   include 'alpb/param_alpb_benzaldehyde.fh'
   include 'alpb/param_alpb_benzene.fh'
   include 'alpb/param_alpb_ch2cl2.fh'
   include 'alpb/param_alpb_chcl3.fh'
   include 'alpb/param_alpb_cs2.fh'
   include 'alpb/param_alpb_dioxane.fh'
   include 'alpb/param_alpb_dmf.fh'
   include 'alpb/param_alpb_dmso.fh'
   include 'alpb/param_alpb_ether.fh'
   include 'alpb/param_alpb_ethylacetate.fh'
   include 'alpb/param_alpb_furane.fh'
   include 'alpb/param_alpb_hexadecane.fh'
   include 'alpb/param_alpb_hexane.fh'
   include 'alpb/param_alpb_nitromethane.fh'
   include 'alpb/param_alpb_octanol.fh'
   include 'alpb/param_alpb_phenol.fh'
   include 'alpb/param_alpb_thf.fh'
   include 'alpb/param_alpb_toluene.fh'
   include 'alpb/param_alpb_water.fh'
   include 'alpb/param_alpb_woctanol.fh'
   include 'alpb/param_alpb_methanol.fh'
   include 'alpb/param_alpb_ethanol.fh'

contains


!> Get ALPB/GBSA parameters
subroutine get_alpb_param(input, mol, error)
   !> Input of ALPB
   type(alpb_input), intent(inout) :: input
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> internal parameters to be used
   type(alpb_parameter), allocatable :: param
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   select case(input%kernel)
   case(1)
      if (input%method == 'gfn2') then
         select case(input%solvent)
         case('acetone');      param = gfn2_acetone
         case('acetonitrile'); param = gfn2_acetonitrile
         case('benzene');      param = gfn2_benzene
         case('ch2cl2','dichlormethane'); param = gfn2_ch2cl2
         case('chcl3','chloroform');      param = gfn2_chcl3
         case('cs2');          param = gfn2_cs2
         case('dmso');         param = gfn2_dmso
         case('ether');        param = gfn2_ether
         case('h2o','water');  param = gfn2_h2o
         case('methanol');     param = gfn2_methanol
         case('thf');          param = gfn2_thf
         case('toluene');      param = gfn2_toluene
         case('dmf');          param = gfn2_dmf
         case('nhexan','n-hexan','nhexane','n-hexane','hexane');
            param = gfn2_nhexan
         end select
      else
         select case(input%solvent)
         case('acetone');      param = gfn1_acetone
         case('acetonitrile'); param = gfn1_acetonitrile
         case('benzene');      param = gfn1_benzene
         case('ch2cl2','dichlormethane'); param = gfn1_ch2cl2
         case('chcl3','chloroform');      param = gfn1_chcl3
         case('cs2');          param = gfn1_cs2
         case('dmso');         param = gfn1_dmso
         case('ether');        param = gfn1_ether
         case('h2o','water');  param = gfn1_h2o
         case('methanol');     param = gfn1_methanol
         case('thf');          param = gfn1_thf
         case('toluene');      param = gfn1_toluene
         end select
      end if
   case(2)
      if (input%method == 'gfn2') then
         select case(input%solvent)
         case('acetone');      param = gfn2_alpb_acetone
         case('acetonitrile'); param = gfn2_alpb_acetonitrile
         case('aniline');      param = gfn2_alpb_aniline
         case('benzaldehyde');      param = gfn2_alpb_benzaldehyde
         case('benzene');      param = gfn2_alpb_benzene
         case('dioxane');      param = gfn2_alpb_dioxane
         case('ethylacetate');      param = gfn2_alpb_ethylacetate
         case('furane');      param = gfn2_alpb_furane
         case('hexadecane');      param = gfn2_alpb_hexadecane
         case('nitromethane');      param = gfn2_alpb_nitromethane
         case('octanol');      param = gfn2_alpb_octanol
         case('woctanol');      param = gfn2_alpb_woctanol
         case('phenol');      param = gfn2_alpb_phenol 
         case('ch2cl2','dichlormethane'); param = gfn2_alpb_ch2cl2
         case('chcl3','chloroform');      param = gfn2_alpb_chcl3
         case('cs2');          param = gfn2_alpb_cs2
         case('dmso');         param = gfn2_alpb_dmso
         case('ether');        param = gfn2_alpb_ether
         case('h2o','water');  param = gfn2_alpb_water
         case('methanol');     param = gfn2_alpb_methanol 
         case('thf');          param = gfn2_alpb_thf
         case('toluene');      param = gfn2_alpb_toluene
         case('dmf');          param = gfn2_alpb_dmf
         case('ethanol');      param = gfn2_alpb_ethanol
         case('nhexan','n-hexan','nhexane','n-hexane','hexane');
            param = gfn2_alpb_hexane
         end select
      else
         select case(input%solvent)
         case('acetone');      param = gfn1_alpb_acetone
         case('acetonitrile'); param = gfn1_alpb_acetonitrile
         case('aniline');      param = gfn1_alpb_aniline
         case('benzaldehyde');      param = gfn1_alpb_benzaldehyde
         case('benzene');      param = gfn1_alpb_benzene
         case('dioxane');      param = gfn1_alpb_dioxane
         case('ethylacetate');      param = gfn1_alpb_ethylacetate
         case('furane');      param = gfn1_alpb_furane
         case('hexadecane');      param = gfn1_alpb_hexadecane
         case('nitromethane');      param = gfn1_alpb_nitromethane
         case('octanol');      param = gfn1_alpb_octanol
         case('woctanol');      param = gfn1_alpb_woctanol
         case('phenol');      param = gfn1_alpb_phenol 
         case('ch2cl2','dichlormethane'); param = gfn1_alpb_ch2cl2
         case('chcl3','chloroform');      param = gfn1_alpb_chcl3
         case('cs2');          param = gfn1_alpb_cs2
         case('dmso');         param = gfn1_alpb_dmso
         case('ether');        param = gfn1_alpb_ether
         case('h2o','water');  param = gfn1_alpb_water
         case('methanol');     param = gfn1_alpb_methanol
         case('ethanol');      param = gfn1_alpb_ethanol
         case('thf');          param = gfn1_alpb_thf
         case('toluene');      param = gfn1_alpb_toluene
         case('dmf');          param = gfn1_alpb_dmf
         case('nhexan','n-hexan','nhexane','n-hexane','hexane');
            param = gfn1_alpb_hexane
         end select
      end if
   end select

   if (.not.allocated(param)) then
      call fatal_error(error, "Unknown solvent")
   end if
 
   call load_alpb_param(input, mol, param)

end subroutine get_alpb_param

!> Get ALPB/GBSA parameters
subroutine load_alpb_param(input, mol, param)
   !> Input of ALPB
   type(alpb_input), intent(inout) :: input
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> internal parameters to be used
   type(alpb_parameter), intent(in) :: param

   !> born scale 
   input%born_scale = param%c1

   !> born offset  
   input%born_offset = param%soset * 0.1_wp * aatoau

   !> dielectric constant
   input%dielectric_const = param%epsv

   if (.not. allocated(input%descreening)) then
      allocate(input%descreening(mol%nid))
   end if

   input%descreening = param%sx(mol%num)

end subroutine load_alpb_param

end module tblite_data_alpb
