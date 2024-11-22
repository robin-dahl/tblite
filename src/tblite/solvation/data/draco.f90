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

!> @file tblite/solvation/data/draco.f90
!> Provides parameters for DRACO

!> DRACO radii scaling parameters
module tblite_solvation_data_draco
   use mctc_env, only: wp
   use mctc_io, only: to_number
   implicit none
   private
   public :: get_alpha_cpcm, get_shift_cpcm, get_kcn_cpcm
   public :: get_alpha_smd, get_shift_smd, get_kcn_smd
   public :: get_alpha_cosmo, get_shift_cosmo, get_kcn_cosmo

   !>  convert bohr (a.u.) to Ångström and back
   real(wp), parameter :: autoaa = 0.52917726_wp
   real(wp), parameter :: aatoau = 1.0_wp/autoaa

   !> CPCM ----
   !> Interface for getting the alpha parameter of the scaling function based on either the atomic number or symbol
   interface get_alpha_cpcm
      module procedure get_alpha_cpcm_num
      module procedure get_alpha_cpcm_sym
   end interface get_alpha_cpcm
   !> Interface for getting the shift parameter of the scaling function based on either the atomic number or symbol
   interface get_shift_cpcm
      module procedure get_shift_cpcm_num
      module procedure get_shift_cpcm_sym
   end interface get_shift_cpcm
   !> Interface for getting the k parameter of the scaling function based on either the atomic number or symbol
   interface get_kcn_cpcm
      module procedure get_kcn_cpcm_num
      module procedure get_kcn_cpcm_sym
   end interface get_kcn_cpcm

   !> SMD -----
   !> Interface for getting the alpha parameter of the scaling function based on either the atomic number or symbol
   interface get_alpha_smd
      module procedure get_alpha_smd_num
      module procedure get_alpha_smd_sym
   end interface get_alpha_smd
   !> Interface for getting the shift parameter of the scaling function based on either the atomic number or symbol
   interface get_shift_smd
      module procedure get_shift_smd_num
      module procedure get_shift_smd_sym
   end interface get_shift_smd
   !> Interface for getting the k parameter of the scaling function based on either the atomic number or symbol
   interface get_kcn_smd
      module procedure get_kcn_smd_num
      module procedure get_kcn_smd_sym
   end interface get_kcn_smd

   !> COSMO -----
   !> Interface for getting the alpha parameter of the scaling function based on either the atomic number or symbol
   interface get_alpha_cosmo
      module procedure get_alpha_cosmo_num
      module procedure get_alpha_cosmo_sym
   end interface get_alpha_cosmo
   !> Interface for getting the shift parameter of the scaling function based on either the atomic number or symbol
   interface get_shift_cosmo
      module procedure get_shift_cosmo_num
      module procedure get_shift_cosmo_sym
   end interface get_shift_cosmo
   !> Interface for getting the k parameter of the scaling function based on either the atomic number or symbol
   interface get_kcn_cosmo
      module procedure get_kcn_cosmo_num
      module procedure get_kcn_cosmo_sym
   end interface get_kcn_cosmo

   !> CPCM ####################################################################
   !> EEQ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !> Water
   real(wp), parameter :: eeq_prefac_water_cpcm(94) = [ &
   &-0.48580467_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.00467076_wp, 0.00731119_wp, 0.14219428_wp, &
   &1.51840280_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.51695360_wp, &
   &0.52575628_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.71772284_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]
   real(wp), parameter :: eeq_expo_water_cpcm(94) = [ &
   &0.73703909_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 31.23704066_wp, 5.18246260_wp, 3.11305778_wp, &
   &-0.75255631_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, -0.27512118_wp, &
   &-0.97118624_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, -0.18646034_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]
   real(wp), parameter :: eeq_k_water_cpcm(94) = [ &
   &2.22602441_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 43.69356793_wp, -4.77192548_wp, -2.41506096_wp, &
   &0.87149428_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.18863602_wp, &
   &2.52012598_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, -2.23857263_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]

   !> Other solvents
   real(wp), parameter :: eeq_prefac_other_solvents_cpcm(94) = [ &
   &0.04334600_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, -1.60595905_wp, -0.00709949_wp, 0.02601606_wp, &
   &14.20750071_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, -1.63941002_wp, &
   &0.44756913_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.19537750_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]
   real(wp), parameter :: eeq_expo_other_solvents_cpcm(94) = [ &
   &0.66291388_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, -0.07256963_wp, 44.96439531_wp, 7.40419900_wp, &
   &-0.32022526_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, -0.12559499_wp, &
   &-1.32277741_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, -1.78741237_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]
   real(wp), parameter :: eeq_k_other_solvents_cpcm(94) = [ &
   &-20.27800374_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, -0.22493811_wp, -12.62008638_wp, 8.46988686_wp, &
   &-0.01542781_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, -0.34163280_wp, &
   &1.41095131_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 4.81934081_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]
   real(wp), parameter :: eeq_o_shift_other_solvents_cpcm(94) = [ &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 1.37365379_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]

   !> SMD #####################################################################
   !> EEQ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !> Water
   real(wp), parameter :: eeq_prefac_water_smd(94) = [ &
   &-0.20299412_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.00116557_wp, 0.00112715_wp, 0.05390754_wp, &
   &0.62174879_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.20563583_wp, &
   &0.04342682_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.35845544_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]
   real(wp), parameter :: eeq_expo_water_smd(94) = [ &
   &0.35663583_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 53.34588775_wp, 121.08319293_wp, 2.14269861_wp, &
   &-0.43639762_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.15907634_wp, &
   &1.21454052_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.16155323_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]
   real(wp), parameter :: eeq_k_water_smd(94) = [ &
   &0.90110573_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 8.75241480_wp, -42.05520406_wp, -0.31820080_wp, &
   &0.80691234_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.62967096_wp, &
   &0.32705664_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, -3.18641147_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]

   !> Other solvents
   real(wp), parameter :: eeq_prefac_other_solvents_smd(94) = [ &
   &0.05447275_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, -1.77608024_wp, -0.00508844_wp, 0.01738052_wp, &
   &12.00126366_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, -1.37395461_wp, &
   &0.51956849_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.13813859_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]
   real(wp), parameter :: eeq_expo_other_solvents_smd(94) = [ &
   &-2.50038085_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, -0.05353459_wp, 43.30983074_wp, -4.67837598_wp, &
   &-0.31514792_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, -0.21519048_wp, &
   &-1.07041829_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, -0.31621463_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]
   real(wp), parameter :: eeq_k_other_solvents_smd(94) = [ &
   &-20.59413200_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, -0.24782848_wp, -12.57240511_wp, 18.63767256_wp, &
   &-0.01085468_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, -0.31346972_wp, &
   &1.43231694_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 5.81656950_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]
   real(wp), parameter :: eeq_o_shift_other_solvents_smd(94) = [ &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 1.50193776_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]

   !> COSMO ###################################################################
   !> EEQ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !> Water
   real(wp), parameter :: eeq_prefac_water_cosmo(94) = [ &
   &-0.30092734_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.00476617_wp, 0.00739481_wp, 0.11798678_wp, &
   &1.59282208_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.54293999_wp, &
   &0.53158958_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.71064467_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]
   real(wp), parameter :: eeq_expo_water_cosmo(94) = [ &
   &0.67947964_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 30.04851006_wp, 3.85013316_wp, 3.41416781_wp, &
   &-0.77593279_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, -0.29072110_wp, &
   &-1.01334593_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, -0.19747693_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]
   real(wp), parameter :: eeq_k_water_cosmo(94) = [ &
   &2.26579612_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 19.17030480_wp, -4.55761058_wp, -2.48750963_wp, &
   &0.96219818_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.21396610_wp, &
   &2.52804811_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, -2.13665118_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]

   !> Other solvents
   real(wp), parameter :: eeq_prefac_other_solvents_cosmo(94) = [ &
   &0.01913447_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, -1.00863834_wp, -0.03035434_wp, 0.01923511_wp, &
   &-0.06368306_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.17129345_wp, &
   &0.55742109_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.09589146_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]
   real(wp), parameter :: eeq_expo_other_solvents_cosmo(94) = [ &
   &4.19343256_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, -0.08382289_wp, 9.38436873_wp, 9.79438998_wp, &
   &0.05316351_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 1.62927856_wp, &
   &-0.46721507_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, -2.74780407_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]
   real(wp), parameter :: eeq_k_other_solvents_cosmo(94) = [ &
   &-37.09431303_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, -0.16521795_wp, -2.33403961_wp, 8.69270475_wp, &
   &-10.03287872_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.23999841_wp, &
   &-0.72711179_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.82036564_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]
   real(wp), parameter :: eeq_o_shift_other_solvents_cosmo(94) = [ &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 1.43381474_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
   &0.0000_wp, 0.0000_wp]

   !!!!!!!!!!!!!

contains

   !!! CPCM -----
   !> Get the alpha parameter of the error function for cpcm based on atomic number
   elemental function get_alpha_cpcm_num(number) result(alpha)
      integer, intent(in) :: number
      real(wp) :: alpha

      alpha = eeq_prefac_water_cpcm(number)

   end function get_alpha_cpcm_num

   !> Get the alpha parameter of the error function for cpcm based on atomic symbol
   elemental function get_alpha_cpcm_sym(symbol) result(alpha)
      !> Element symbol
      character(len=*), intent(in) :: symbol
      !> Alpha parameter
      real(wp) :: alpha

      alpha = eeq_prefac_water_cpcm(to_number(symbol))
   end function get_alpha_cpcm_sym

   !> Get the shift parameter of the error function for cpcm based on atomic number
   elemental function get_shift_cpcm_num(number) result(shift)
      integer, intent(in) :: number
      real(wp) :: shift

      shift = eeq_expo_water_cpcm(number)

   end function get_shift_cpcm_num

   !> Get the shift parameter of the error function for cpcm based on atomic symbol
   elemental function get_shift_cpcm_sym(symbol) result(shift)
      character(len=*), intent(in) :: symbol
      real(wp) :: shift

      shift = eeq_expo_water_cpcm(to_number(symbol))
   end function get_shift_cpcm_sym

   !> Get the k parameter of the error function for cpcm based on atomic number
   elemental function get_kcn_cpcm_num(number) result(kcn)
      integer, intent(in) :: number
      real(wp) :: kcn

      kcn = eeq_k_water_cpcm(number)

   end function get_kcn_cpcm_num

   !> Get the k parameter of the error function for cpcm based on atomic symbol
   elemental function get_kcn_cpcm_sym(symbol) result(kcn)
      character(len=*), intent(in) :: symbol
      real(wp) :: kcn

      kcn = eeq_k_water_cpcm(to_number(symbol))

   end function get_kcn_cpcm_sym

    !!! SMD -----
   !> Get the alpha parameter of the error function for smd based atomic number
   elemental function get_alpha_smd_num(number) result(alpha)
      integer, intent(in) :: number
      real(wp) :: alpha

      alpha = eeq_prefac_water_smd(number)

   end function get_alpha_smd_num

   !> Get the alpha parameter of the error function for smd based atomic symbol
   elemental function get_alpha_smd_sym(symbol) result(alpha)
      character(len=*), intent(in) :: symbol
      real(wp) :: alpha

      alpha = eeq_prefac_water_smd(to_number(symbol))

   end function get_alpha_smd_sym

   !> Get the shift parameter of the error function for smd based on atomic number
   elemental function get_shift_smd_num(number) result(shift)
      integer, intent(in) :: number
      real(wp) :: shift

      shift = eeq_expo_water_smd(number)

   end function get_shift_smd_num

   !> Get the shift parameter of the error function for smd based on atomic symbol
   elemental function get_shift_smd_sym(symbol) result(shift)
      character(len=*), intent(in) :: symbol
      real(wp) :: shift

      shift = eeq_expo_water_smd(to_number(symbol))

   end function get_shift_smd_sym

   !> Get the k parameter of the error function for smd based on atomic number
   elemental function get_kcn_smd_num(number) result(kcn)
      integer, intent(in) :: number
      real(wp) :: kcn

      kcn = eeq_k_water_smd(number)

   end function get_kcn_smd_num

   !> Get the k parameter of the error function for smd based on atomic symbol
   elemental function get_kcn_smd_sym(symbol) result(kcn)
      character(len=*), intent(in) :: symbol
      real(wp) :: kcn

      kcn = eeq_k_water_smd(to_number(symbol))

   end function get_kcn_smd_sym

    !!! COSMO -----
   !> Get the alpha parameter of the error function for cosmo based atomic number
   elemental function get_alpha_cosmo_num(number) result(alpha)
      integer, intent(in) :: number
      real(wp) :: alpha

      alpha = eeq_prefac_water_cosmo(number)

   end function get_alpha_cosmo_num

   !> Get the alpha parameter of the error function for cosmo based atomic symbol
   elemental function get_alpha_cosmo_sym(symbol) result(alpha)
      character(len=*), intent(in) :: symbol
      real(wp) :: alpha

      alpha = eeq_prefac_water_cosmo(to_number(symbol))

   end function get_alpha_cosmo_sym

   !> Get the shift parameter of the error function for cosmo based on atomic number
   elemental function get_shift_cosmo_num(number) result(shift)
      integer, intent(in) :: number
      real(wp) :: shift

      shift = eeq_expo_water_cosmo(number)

   end function get_shift_cosmo_num

   !> Get the shift parameter of the error function for cosmo based on atomic symbol
   elemental function get_shift_cosmo_sym(symbol) result(shift)
      character(len=*), intent(in) :: symbol
      real(wp) :: shift

      shift = eeq_expo_water_cosmo(to_number(symbol))

   end function get_shift_cosmo_sym

   !> Get the k parameter of the error function for cosmo based on atomic number
   elemental function get_kcn_cosmo_num(number) result(kcn)
      integer, intent(in) :: number
      real(wp) :: kcn

      kcn = eeq_k_water_cosmo(number)

   end function get_kcn_cosmo_num

   !> Get the k parameter of the error function for cosmo based on atomic symbol
   elemental function get_kcn_cosmo_sym(symbol) result(kcn)
      character(len=*), intent(in) :: symbol
      real(wp) :: kcn

      kcn = eeq_k_water_cosmo(to_number(symbol))

      !-999

   end function get_kcn_cosmo_sym

end module tblite_solvation_data_draco
