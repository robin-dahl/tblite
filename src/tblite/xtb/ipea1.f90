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

module tblite_xtb_ipea1
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_ortho, only : orthogonalize
   use tblite_basis_type, only : basis_type, new_basis, cgto_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_coulomb_charge, only : new_effective_coulomb, harmonic_average
   use tblite_coulomb_thirdorder, only : new_onsite_thirdorder
   use tblite_disp, only : d3_dispersion, new_d3_dispersion
   use tblite_ncoord, only : new_ncoord
   use tblite_param_paulingen, only : get_pauling_en
   use tblite_repulsion, only : new_repulsion
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_h0, only : new_hamiltonian
   use tblite_xtb_spec, only : tb_h0spec
   implicit none
   private

   public :: new_ipea1_calculator
   public :: ipea1_h0spec

   integer, parameter :: max_elem = 86
   integer, parameter :: max_shell = 3

   !> Use older eV to Eh conversion for consistency
   real(wp), parameter :: evtoau = 1.0_wp / 27.21138505_wp

   !> Exponents of repulsion term for IPEA1-xTB repulsion
   real(wp), parameter :: rep_alpha(max_elem) = [&
      & 2.133166_wp, 1.382907_wp, 0.630666_wp, 0.830517_wp, 1.066430_wp, &
      & 1.206815_wp, 1.622255_wp, 2.002624_wp, 2.538852_wp, 3.038727_wp, &
      & 0.678744_wp, 0.796201_wp, 0.929219_wp, 0.948165_wp, 1.047895_wp, &
      & 1.170894_wp, 1.404155_wp, 1.323756_wp, 0.581529_wp, 0.665588_wp, &
      & 0.760117_wp, 0.817434_wp, 1.067056_wp, 0.984779_wp, 1.034098_wp, &
      & 1.136974_wp, 1.202626_wp, 1.405045_wp, 1.224139_wp, 1.093312_wp, &
      & 0.982124_wp, 1.076416_wp, 1.227794_wp, 1.229699_wp, 1.153580_wp, &
      & 1.335287_wp, 0.554032_wp, 0.657904_wp, 0.760144_wp, 0.709735_wp, &
      & 0.909696_wp, 0.961038_wp, 1.028240_wp, 1.066144_wp, 1.131380_wp, &
      & 1.226066_wp, 1.029976_wp, 0.928849_wp, 0.822482_wp, 0.954713_wp, &
      & 0.976956_wp, 0.954405_wp, 0.949181_wp, 1.074785_wp, 0.579919_wp, &
      & 0.606485_wp, 1.311200_wp, 0.839861_wp, 0.847281_wp, 0.854701_wp, &
      & 0.862121_wp, 0.869541_wp, 0.876961_wp, 0.884381_wp, 0.891801_wp, &
      & 0.899221_wp, 0.906641_wp, 0.914061_wp, 0.921481_wp, 0.928901_wp, &
      & 0.936321_wp, 0.831440_wp, 1.001602_wp, 1.054268_wp, 1.295110_wp, &
      & 1.163798_wp, 1.280092_wp, 1.356651_wp, 1.464120_wp, 1.147709_wp, &
      & 0.959310_wp, 0.975700_wp, 0.952375_wp, 1.013118_wp, 0.964652_wp, &
      & 0.998641_wp]

   !> Effective nuclear charge for IPEA1-xTB repulsion
   real(wp), parameter :: rep_zeff(max_elem) = [&
      &  1.292559_wp,  0.440231_wp,  2.747587_wp,  4.557601_wp,  5.255954_wp, &
      &  4.427268_wp,  5.500005_wp,  5.178714_wp,  7.370899_wp,  9.102523_wp, &
      &  9.192606_wp, 12.176772_wp, 16.283595_wp, 16.898359_wp, 16.296596_wp, &
      & 15.831707_wp, 17.000000_wp, 17.153132_wp, 20.831436_wp, 19.840212_wp, &
      & 23.231297_wp, 15.054336_wp, 21.158314_wp, 22.335723_wp, 21.161245_wp, &
      & 24.885518_wp, 25.477503_wp, 17.151988_wp, 24.562118_wp, 35.552646_wp, &
      & 36.575935_wp, 39.096346_wp, 38.315669_wp, 34.959310_wp, 35.000000_wp, &
      & 36.000000_wp, 39.653032_wp, 38.924904_wp, 39.000000_wp, 37.288814_wp, &
      & 35.440124_wp, 38.163829_wp, 43.000000_wp, 44.492732_wp, 45.241537_wp, &
      & 41.939426_wp, 47.198422_wp, 49.016827_wp, 51.718417_wp, 54.503455_wp, &
      & 52.995653_wp, 50.537298_wp, 53.000000_wp, 52.500985_wp, 65.029838_wp, &
      & 46.532974_wp, 48.337542_wp, 30.638143_wp, 34.130718_wp, 37.623294_wp, &
      & 41.115870_wp, 44.608445_wp, 48.101021_wp, 51.593596_wp, 55.086172_wp, &
      & 58.578748_wp, 62.071323_wp, 65.563899_wp, 69.056474_wp, 72.549050_wp, &
      & 76.041625_wp, 55.350542_wp, 50.858491_wp, 74.000000_wp, 75.000000_wp, &
      & 76.000000_wp, 77.000000_wp, 78.000000_wp, 79.000000_wp, 80.000000_wp, &
      & 84.050696_wp,101.661715_wp, 83.000000_wp, 84.000000_wp, 85.000000_wp, &
      & 86.000000_wp]

   !> Atomic hardnesses used in second order electrostatics
   real(wp), parameter :: chemical_hardness(max_elem) = [&
      & 0.404839_wp, 1.441379_wp, 0.196176_wp, 0.315230_wp, 0.608113_wp, &
      & 0.639692_wp, 0.492599_wp, 0.717333_wp, 0.574104_wp, 0.612878_wp, &
      & 0.146475_wp, 0.332890_wp, 0.221658_wp, 0.438331_wp, 0.637801_wp, &
      & 0.455078_wp, 0.559471_wp, 0.529906_wp, 0.114358_wp, 0.134187_wp, &
      & 0.897973_wp, 1.678103_wp, 1.094925_wp, 0.317872_wp, 0.231875_wp, &
      & 0.900000_wp, 0.213854_wp, 0.257961_wp, 0.232867_wp, 0.384579_wp, &
      & 0.351929_wp, 0.878878_wp, 0.417457_wp, 0.546433_wp, 0.732530_wp, &
      & 0.820312_wp, 0.075735_wp, 0.122861_wp, 0.351290_wp, 0.176042_wp, &
      & 0.180389_wp, 0.392238_wp, 0.405474_wp, 0.305394_wp, 0.293973_wp, &
      & 0.311234_wp, 0.440109_wp, 0.332147_wp, 0.285832_wp, 0.372861_wp, &
      & 0.535840_wp, 0.283009_wp, 0.520472_wp, 0.935394_wp, 0.085110_wp, &
      & 0.137819_wp, 0.495969_wp, 0.350000_wp, 0.342306_wp, 0.334612_wp, &
      & 0.326917_wp, 0.319223_wp, 0.311529_wp, 0.303835_wp, 0.296140_wp, &
      & 0.288446_wp, 0.280752_wp, 0.273058_wp, 0.265364_wp, 0.257669_wp, &
      & 0.249975_wp, 0.280206_wp, 0.209307_wp, 0.298893_wp, 0.192431_wp, &
      & 1.511128_wp, 0.391625_wp, 0.306418_wp, 0.430317_wp, 0.335841_wp, &
      & 0.672340_wp, 1.000000_wp, 1.081470_wp, 1.091248_wp, 1.264162_wp, &
      & 0.798170_wp]

   !> Scaling factors for shell electrostatics
   real(wp), parameter :: shell_hardness(0:2, max_elem) = 1.0_wp + reshape([&
      & 0.0_wp, 0.0000000_wp, 0.0000000_wp,  0.0_wp, 0.0000000_wp, 0.0000000_wp, &
      & 0.0_wp,-0.2241991_wp, 0.0000000_wp,  0.0_wp,-0.0053663_wp, 0.0000000_wp, &
      & 0.0_wp,-0.2662342_wp, 0.0000000_wp,  0.0_wp,-0.0273283_wp, 0.0000000_wp, &
      & 0.0_wp, 0.0701245_wp, 0.0000000_wp,  0.0_wp,-0.2619316_wp, 0.0000000_wp, &
      & 0.0_wp, 0.0808738_wp, 0.0000000_wp,  0.0_wp,-0.3892542_wp, 0.0000000_wp, &
      & 0.0_wp,-0.4545627_wp, 0.0000000_wp,  0.0_wp, 0.1823758_wp, 0.0000000_wp, &
      & 0.0_wp, 0.0503564_wp, 0.0000000_wp,  0.0_wp,-0.5925834_wp, 0.0000000_wp, &
      & 0.0_wp,-0.3150116_wp, 0.0000000_wp,  0.0_wp,-0.1956114_wp, 0.0000000_wp, &
      & 0.0_wp,-0.3994898_wp, 0.0000000_wp,  0.0_wp,-0.1450000_wp, 0.0000000_wp, &
      & 0.0_wp,-0.5332978_wp, 0.0000000_wp,  0.0_wp, 1.1522018_wp, 0.0000000_wp, &
      & 0.0_wp,-0.2000000_wp,-0.6976887_wp,  0.0_wp,-0.1500000_wp,-0.7677638_wp, &
      & 0.0_wp,-0.2000000_wp,-0.5176731_wp,  0.0_wp,-0.2500000_wp, 0.6084446_wp, &
      & 0.0_wp,-0.2500000_wp, 1.2062125_wp,  0.0_wp,-0.2000000_wp,-0.3348701_wp, &
      & 0.0_wp,-0.2000000_wp, 0.3750103_wp,  0.0_wp,-0.2000000_wp, 0.2009132_wp, &
      & 0.0_wp, 0.0000000_wp,-0.1582039_wp,  0.0_wp,-0.1515628_wp, 0.0000000_wp, &
      & 0.0_wp,-0.3894500_wp, 0.1000000_wp,  0.0_wp,-0.6463392_wp,-0.1300000_wp, &
      & 0.0_wp,-0.0758566_wp,-0.1000000_wp,  0.0_wp,-0.4319966_wp,-0.2500000_wp, &
      & 0.0_wp,-0.1440020_wp,-0.1000000_wp,  0.0_wp,-0.3743296_wp, 0.0000000_wp, &
      & 0.0_wp,-0.5181667_wp, 0.0000000_wp,  0.0_wp,-0.8003590_wp, 0.0000000_wp, &
      & 0.0_wp,-0.0800000_wp,-0.4159186_wp,  0.0_wp,-0.2500000_wp, 0.3471956_wp, &
      & 0.0_wp,-0.2000000_wp, 0.9898269_wp,  0.0_wp,-0.2500000_wp,-0.2957569_wp, &
      & 0.0_wp,-0.2000000_wp, 0.2642680_wp,  0.0_wp,-0.1500000_wp, 0.1772831_wp, &
      & 0.0_wp,-0.2500000_wp, 0.3782936_wp,  0.0_wp,-0.2500000_wp, 0.2865563_wp, &
      & 0.0_wp,-0.1500000_wp,-0.1981564_wp,  0.0_wp,-0.3010874_wp, 0.0000000_wp, &
      & 0.0_wp,-0.3631014_wp, 0.0000000_wp,  0.0_wp,-0.3294996_wp,-0.1500000_wp, &
      & 0.0_wp,-0.1806602_wp,-0.2000000_wp,  0.0_wp,-0.1102393_wp,-0.2000000_wp, &
      & 0.0_wp,-0.0387987_wp,-0.1500000_wp,  0.0_wp,-0.3435282_wp,-0.1500000_wp, &
      & 0.0_wp,-0.7035550_wp, 0.0000000_wp,  0.0_wp,-0.8801363_wp, 0.0000000_wp, &
      & 0.0_wp,-0.1500000_wp,-0.6396752_wp,  0.0_wp,-0.1500000_wp,-0.5245538_wp, &
      & 0.0_wp,-0.1500000_wp,-0.5064761_wp,  0.0_wp,-0.1500000_wp,-0.4883984_wp, &
      & 0.0_wp,-0.1500000_wp,-0.4703207_wp,  0.0_wp,-0.1500000_wp,-0.4522429_wp, &
      & 0.0_wp,-0.1500000_wp,-0.4341652_wp,  0.0_wp,-0.1500000_wp,-0.4160875_wp, &
      & 0.0_wp,-0.1500000_wp,-0.3980098_wp,  0.0_wp,-0.1500000_wp,-0.3799321_wp, &
      & 0.0_wp,-0.1500000_wp,-0.3618544_wp,  0.0_wp,-0.1500000_wp,-0.3437767_wp, &
      & 0.0_wp,-0.1500000_wp,-0.3256989_wp,  0.0_wp,-0.1500000_wp,-0.3076212_wp, &
      & 0.0_wp,-0.1500000_wp,-0.2895435_wp,  0.0_wp,-0.1500000_wp, 0.0920885_wp, &
      & 0.0_wp,-0.1500000_wp, 0.5567170_wp,  0.0_wp,-0.1500000_wp, 0.2658721_wp, &
      & 0.0_wp,-0.2000000_wp, 0.1569632_wp,  0.0_wp,-0.1000000_wp,-0.1931486_wp, &
      & 0.0_wp,-0.1500000_wp, 0.7240972_wp,  0.0_wp,-0.2000000_wp, 0.0257168_wp, &
      & 0.0_wp,-0.1500000_wp,-0.3868172_wp,  0.0_wp,-0.4553637_wp, 0.0000000_wp, &
      & 0.0_wp,-0.7973634_wp, 0.0000000_wp,  0.0_wp,-0.8045028_wp, 0.0000000_wp, &
      & 0.0_wp,-0.6727954_wp, 0.0000000_wp,  0.0_wp,-0.3955914_wp, 0.0000000_wp, &
      & 0.0_wp,-0.3402676_wp, 0.0000000_wp,  0.0_wp,-0.2380762_wp, 0.0000000_wp],&
      & shape(shell_hardness))

   !> Third order Hubbard derivatives
   real(wp), parameter :: p_hubbard_derivs(max_elem) = 0.1_wp * [&
      & 1.500000_wp, 1.500000_wp, 1.846487_wp, 1.500000_wp, 0.500000_wp, &
      & 1.500000_wp,-1.000000_wp,-2.000000_wp, 2.000000_wp, 1.600000_wp, &
      & 1.200000_wp, 1.100000_wp, 1.200000_wp, 1.500000_wp, 1.616190_wp, &
      & 1.020159_wp, 2.395510_wp, 0.829312_wp, 0.732923_wp, 1.116963_wp, &
      &-1.784792_wp,-0.000358_wp, 0.800000_wp, 0.800000_wp, 0.300000_wp, &
      & 0.500000_wp, 0.300000_wp, 1.000000_wp, 0.920598_wp, 1.400000_wp, &
      & 1.400000_wp, 1.400000_wp, 1.300000_wp, 1.799283_wp,-0.500000_wp, &
      & 1.000000_wp, 1.500000_wp, 1.300000_wp, 1.400000_wp, 0.460210_wp, &
      & 0.212639_wp,-0.335700_wp, 0.500000_wp, 0.001205_wp, 0.622690_wp, &
      & 0.500000_wp, 0.021628_wp, 1.362587_wp, 1.063557_wp, 1.746331_wp, &
      &-0.678152_wp, 1.387384_wp,-0.500000_wp,-0.800000_wp, 1.500000_wp, &
      & 1.500000_wp, 1.500000_wp, 1.200000_wp, 1.200000_wp, 1.200000_wp, &
      & 1.200000_wp, 1.200000_wp, 1.200000_wp, 1.200000_wp, 1.200000_wp, &
      & 1.200000_wp, 1.200000_wp, 1.200000_wp, 1.200000_wp, 1.200000_wp, &
      & 1.200000_wp, 0.436362_wp,-0.294270_wp, 0.139734_wp, 0.300000_wp, &
      &-0.170295_wp,-1.130878_wp, 0.291755_wp, 2.232497_wp,-0.040460_wp, &
      & 0.992527_wp, 1.500000_wp, 1.179183_wp,-0.035874_wp,-0.860502_wp, &
      &-0.838429_wp]

   !> Number of shells
   integer, parameter :: nshell(max_elem) = [&
      & 2, 1, 2, 2, 3, 3, 3, 3, 3, 3, 2, 2, 3, 3, 3, 3, 3, 3, 2, 3, &
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, &
      & 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, &
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, &
      & 2, 2, 2, 3, 3, 3]

   !> Angular momentum of each shell
   integer, parameter :: ang_shell(max_shell, max_elem) = reshape([&
      & 0, 0, 0,  0, 0, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0, &
      & 0, 1, 0,  0, 1, 0,  0, 1, 2,  0, 1, 0,  0, 1, 0,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 0,  0, 1, 2,  2, 0, 1, &
      & 2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1, &
      & 2, 0, 1,  0, 1, 0,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 0,  0, 1, 2,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1, &
      & 2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  0, 1, 0,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 0,  0, 1, 2, &
      & 2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1, &
      & 2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1, &
      & 2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1,  2, 0, 1, &
      & 2, 0, 1,  2, 0, 1,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2], shape(ang_shell))

   !> Principal quantum number of each shell
   integer, parameter :: principal_quantum_number(max_shell, max_elem) = reshape([&
      & 1, 2, 0,  1, 0, 0,  2, 2, 0,  2, 2, 0,  2, 2, 3,  2, 2, 3,  2, 2, 3, &
      & 2, 2, 3,  2, 2, 3,  2, 2, 3,  2, 2, 0,  2, 2, 0,  3, 3, 3,  3, 3, 3, &
      & 3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  4, 4, 0,  4, 4, 3,  3, 4, 4, &
      & 3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4,  3, 4, 4, &
      & 3, 4, 4,  4, 4, 0,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4, &
      & 4, 4, 4,  5, 5, 0,  5, 5, 4,  4, 5, 5,  4, 5, 5,  4, 5, 5,  4, 5, 5, &
      & 4, 5, 5,  4, 5, 5,  4, 5, 5,  4, 5, 5,  4, 5, 5,  5, 5, 0,  5, 5, 5, &
      & 5, 5, 5,  5, 5, 5,  5, 5, 5,  5, 5, 5,  5, 5, 5,  5, 5, 0,  5, 5, 4, &
      & 5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6, &
      & 5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6, &
      & 5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6,  5, 6, 6, &
      & 5, 6, 6,  5, 6, 6,  6, 6, 0,  6, 6, 0,  6, 6, 0,  6, 6, 0,  6, 6, 5, &
      & 6, 6, 5,  6, 6, 5], shape(principal_quantum_number))

   !> Shell polynomials to scale Hamiltonian elements
   real(wp), parameter :: p_shpoly(0:2, max_elem) = 0.01_wp * reshape([&
      &  0.000000_wp,  0.000000_wp,  0.000000_wp,   8.084149_wp,  0.000000_wp,  0.000000_wp, &
      &-13.007335_wp,  6.802043_wp,  0.000000_wp, -22.404993_wp, -6.611030_wp,  0.000000_wp, &
      & -5.397599_wp, 14.708712_wp,  0.000000_wp,  -2.692247_wp,  8.177943_wp,  0.000000_wp, &
      &-10.560649_wp,  0.940101_wp,  0.000000_wp,  -9.907396_wp,  3.920125_wp,  0.000000_wp, &
      &-11.637428_wp, -9.003663_wp,  0.000000_wp,  -2.115896_wp,-15.124326_wp,  0.000000_wp, &
      &  7.549004_wp, 13.576090_wp,  0.000000_wp, -18.095421_wp, 21.451533_wp,  0.000000_wp, &
      &-21.085827_wp, 24.805127_wp, 26.405814_wp, -14.201582_wp, -3.893343_wp, 25.499221_wp, &
      &-20.056162_wp, -7.612324_wp, 27.896649_wp, -20.287694_wp, -7.124132_wp, 17.511019_wp, &
      & -9.341919_wp,-12.605796_wp, 14.905217_wp,  -0.082808_wp, -9.217948_wp, 12.204172_wp, &
      & 12.482844_wp, 22.323655_wp,  0.000000_wp, -11.421376_wp, 14.628284_wp, 10.129602_wp, &
      & 12.477046_wp, 46.956747_wp,-10.433716_wp,   8.253155_wp, 30.482831_wp,-18.082141_wp, &
      &-15.627155_wp, 21.787754_wp,-23.674752_wp,   7.351214_wp, 18.517229_wp,-16.744068_wp, &
      &  4.364955_wp,  9.242032_wp,-16.231235_wp, -10.476481_wp, 11.176576_wp,-22.084984_wp, &
      &  4.121503_wp, 14.917653_wp,-20.140234_wp,  10.160508_wp, 12.528539_wp,-20.892120_wp, &
      & -8.590501_wp,  7.141500_wp,-15.142399_wp, -20.097296_wp, -0.515231_wp,  0.000000_wp, &
      & -9.919005_wp, 10.598388_wp, 19.671655_wp, -18.139462_wp, -6.247789_wp, 19.611782_wp, &
      &-17.489686_wp,-10.836128_wp, 17.626604_wp, -17.147845_wp,-11.828189_wp, 15.353343_wp, &
      &-17.815502_wp,-14.058044_wp,  5.468245_wp, -25.437273_wp,-12.813227_wp, 10.440712_wp, &
      & -7.450752_wp, 16.670533_wp,  0.000000_wp,  -6.087125_wp,  2.115262_wp, 17.076466_wp, &
      & 10.950764_wp, 45.679760_wp,-28.061976_wp,  30.048862_wp, 24.685487_wp,-14.065000_wp, &
      & 15.379439_wp, 31.393684_wp,-25.998052_wp,   5.815301_wp, 20.432554_wp,-18.092116_wp, &
      & 24.977603_wp,  1.953838_wp,-23.231470_wp,  15.281981_wp,  1.340798_wp,-23.099524_wp, &
      & 10.450086_wp, 15.559547_wp,-23.540560_wp,  10.407257_wp, 16.980983_wp,-17.752834_wp, &
      &-12.856324_wp, -2.440401_wp, -1.478311_wp, -17.616123_wp,  6.247124_wp,  0.000000_wp, &
      &-10.488459_wp, 13.780348_wp,  5.584366_wp, -14.906472_wp, -2.497008_wp, 10.683419_wp, &
      &-15.500737_wp, -8.070432_wp,  5.919672_wp, -17.920966_wp,-11.924030_wp, 25.088164_wp, &
      &-21.954071_wp,-10.823970_wp, 12.522287_wp, -22.530281_wp,-16.667114_wp,  8.021956_wp, &
      & -1.460631_wp, 15.879494_wp,  0.000000_wp,  -5.468018_wp,  4.368854_wp, 14.328052_wp, &
      & -3.988102_wp, 40.847293_wp,-44.208463_wp,   6.148475_wp, 42.873822_wp,-36.440945_wp, &
      &  7.806576_wp, 42.846148_wp,-36.021673_wp,   9.464678_wp, 42.818474_wp,-35.602402_wp, &
      & 11.122779_wp, 42.790801_wp,-35.183130_wp,  12.780881_wp, 42.763127_wp,-34.763859_wp, &
      & 14.438982_wp, 42.735454_wp,-34.344587_wp,  16.097083_wp, 42.707780_wp,-33.925315_wp, &
      & 17.755185_wp, 42.680106_wp,-33.506044_wp,  19.413286_wp, 42.652433_wp,-33.086772_wp, &
      & 21.071387_wp, 42.624759_wp,-32.667501_wp,  22.729489_wp, 42.597085_wp,-32.248229_wp, &
      & 24.387590_wp, 42.569412_wp,-31.828957_wp,  26.045692_wp, 42.541738_wp,-31.409686_wp, &
      & 27.703793_wp, 42.514065_wp,-30.990414_wp,  15.014122_wp, 25.669742_wp,-20.141582_wp, &
      & 29.782424_wp, 38.084212_wp,-23.077812_wp,  35.195571_wp, 21.973526_wp,-12.013846_wp, &
      & 17.414335_wp, -0.067497_wp,-23.115824_wp,  28.734400_wp, -0.984767_wp,-15.758126_wp, &
      & 25.774929_wp,  1.135126_wp,-17.060525_wp,  38.415536_wp,  5.655587_wp,-17.070367_wp, &
      &-20.627998_wp,  0.482588_wp,-13.276408_wp, -11.990346_wp,  2.697576_wp,  0.000000_wp, &
      & -1.275250_wp, -1.152862_wp,  0.000000_wp,  -6.987542_wp, -6.568936_wp,  0.000000_wp, &
      &-28.955169_wp,-13.130855_wp,  0.000000_wp, -18.477865_wp,-14.037423_wp, 13.809093_wp, &
      &-21.965390_wp,-12.804436_wp, 16.836546_wp, -22.139701_wp,-20.539955_wp, 17.249637_wp],&
      & shape(p_shpoly))

   !> Reference occupation of the atom
   real(wp), parameter :: reference_occ(0:2, max_elem) = reshape([&
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp, &
      & 2.0_wp, 3.0_wp, 0.0_wp,  2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp, &
      & 2.0_wp, 6.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 2.0_wp,  2.0_wp, 0.0_wp, 3.0_wp,  2.0_wp, 0.0_wp, 4.0_wp, &
      & 2.0_wp, 0.0_wp, 5.0_wp,  2.0_wp, 0.0_wp, 6.0_wp,  2.0_wp, 0.0_wp, 7.0_wp, &
      & 2.0_wp, 0.0_wp, 8.0_wp,  2.0_wp, 0.0_wp, 9.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 2.0_wp,  2.0_wp, 0.0_wp, 3.0_wp,  2.0_wp, 0.0_wp, 4.0_wp, &
      & 2.0_wp, 0.0_wp, 5.0_wp,  2.0_wp, 0.0_wp, 6.0_wp,  2.0_wp, 0.0_wp, 7.0_wp, &
      & 2.0_wp, 0.0_wp, 8.0_wp,  2.0_wp, 0.0_wp, 9.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp, &
      & 2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 1.0_wp,  2.0_wp, 0.0_wp, 2.0_wp, &
      & 2.0_wp, 0.0_wp, 3.0_wp,  2.0_wp, 0.0_wp, 4.0_wp,  2.0_wp, 0.0_wp, 5.0_wp, &
      & 2.0_wp, 0.0_wp, 6.0_wp,  2.0_wp, 0.0_wp, 7.0_wp,  2.0_wp, 0.0_wp, 8.0_wp, &
      & 2.0_wp, 0.0_wp, 9.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 1.0_wp, 0.0_wp, &
      & 2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp,  2.0_wp, 4.0_wp, 0.0_wp, &
      & 2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp], shape(reference_occ))

   !> Exponent of the Slater function
   real(wp), parameter :: slater_exponent(max_shell, max_elem) = reshape([&
      & 1.306199_wp, 0.398421_wp, 0.000000_wp,  2.133698_wp, 0.000000_wp, 0.000000_wp, &
      & 0.675700_wp, 0.643466_wp, 0.000000_wp,  0.976724_wp, 1.325974_wp, 0.000000_wp, &
      & 1.426737_wp, 1.472693_wp, 1.157688_wp,  1.979793_wp, 1.850537_wp, 0.633095_wp, &
      & 2.126160_wp, 2.068425_wp, 0.671040_wp,  2.575264_wp, 2.159448_wp, 0.958679_wp, &
      & 2.304073_wp, 2.314354_wp, 0.570243_wp,  3.200000_wp, 2.294365_wp, 2.684436_wp, &
      & 0.566702_wp, 0.391147_wp, 0.000000_wp,  0.947394_wp, 0.538324_wp, 0.000000_wp, &
      & 1.497753_wp, 1.232966_wp, 0.606937_wp,  1.521960_wp, 1.609138_wp, 1.168971_wp, &
      & 2.002025_wp, 1.791833_wp, 1.327886_wp,  2.492075_wp, 1.997904_wp, 1.747541_wp, &
      & 2.559226_wp, 2.054846_wp, 1.762522_wp,  3.502323_wp, 2.287983_wp, 1.761181_wp, &
      & 0.841791_wp, 0.771618_wp, 0.000000_wp,  1.321845_wp, 0.734954_wp, 0.947032_wp, &
      & 1.510174_wp, 1.278989_wp, 0.873934_wp,  1.526221_wp, 1.746162_wp, 0.810143_wp, &
      & 1.706872_wp, 1.520194_wp, 1.100000_wp,  1.860646_wp, 1.214281_wp, 1.130000_wp, &
      & 1.735856_wp, 1.025035_wp, 1.270000_wp,  2.164423_wp, 1.738248_wp, 1.300000_wp, &
      & 2.205668_wp, 1.091436_wp, 1.350000_wp,  2.226169_wp, 1.315045_wp, 1.350000_wp, &
      & 2.386767_wp, 1.388075_wp, 1.350000_wp,  1.739691_wp, 1.291181_wp, 0.000000_wp, &
      & 2.025803_wp, 1.524895_wp, 0.524764_wp,  2.516954_wp, 1.671762_wp, 0.799301_wp, &
      & 2.478788_wp, 1.905779_wp, 1.176250_wp,  2.998140_wp, 2.073392_wp, 1.745074_wp, &
      & 2.886237_wp, 2.190987_wp, 1.789395_wp,  2.828105_wp, 1.965472_wp, 1.512609_wp, &
      & 0.809529_wp, 0.950253_wp, 0.000000_wp,  1.458742_wp, 0.730658_wp, 1.028147_wp, &
      & 2.300000_wp, 1.593058_wp, 1.170000_wp,  1.787871_wp, 1.411952_wp, 1.230000_wp, &
      & 1.552976_wp, 1.475705_wp, 1.200000_wp,  1.765483_wp, 1.920748_wp, 1.220000_wp, &
      & 2.120497_wp, 1.789115_wp, 1.250000_wp,  2.352683_wp, 1.883645_wp, 1.370000_wp, &
      & 2.436353_wp, 2.000000_wp, 1.470000_wp,  2.380741_wp, 1.984859_wp, 1.550000_wp, &
      & 2.644658_wp, 2.056123_wp, 1.620000_wp,  1.961223_wp, 1.409748_wp, 0.000000_wp, &
      & 2.217568_wp, 1.536904_wp, 0.880728_wp,  2.443152_wp, 1.782901_wp, 0.936038_wp, &
      & 2.687886_wp, 2.003834_wp, 1.078859_wp,  2.873897_wp, 2.169037_wp, 1.511922_wp, &
      & 3.117622_wp, 2.248195_wp, 1.831809_wp,  3.128524_wp, 2.316580_wp, 1.888452_wp, &
      & 0.779877_wp, 0.810404_wp, 0.000000_wp,  1.387083_wp, 0.532658_wp, 0.853415_wp, &
      & 3.000000_wp, 1.492677_wp, 1.350000_wp,  3.000000_wp, 1.553483_wp, 1.380859_wp, &
      & 2.992307_wp, 1.578839_wp, 1.385620_wp,  2.984614_wp, 1.604196_wp, 1.390381_wp, &
      & 2.976922_wp, 1.629552_wp, 1.395142_wp,  2.969229_wp, 1.654909_wp, 1.399903_wp, &
      & 2.961536_wp, 1.680265_wp, 1.404664_wp,  2.953843_wp, 1.705622_wp, 1.409425_wp, &
      & 2.946150_wp, 1.730979_wp, 1.414186_wp,  2.938457_wp, 1.756335_wp, 1.418947_wp, &
      & 2.930765_wp, 1.781692_wp, 1.423708_wp,  2.923072_wp, 1.807048_wp, 1.428469_wp, &
      & 2.915379_wp, 1.832405_wp, 1.433230_wp,  2.907686_wp, 1.857761_wp, 1.437991_wp, &
      & 2.899993_wp, 1.883118_wp, 1.442752_wp,  2.433439_wp, 1.918710_wp, 1.450000_wp, &
      & 1.747032_wp, 1.980803_wp, 1.400000_wp,  2.151802_wp, 1.685637_wp, 1.400000_wp, &
      & 2.127769_wp, 1.886445_wp, 1.450000_wp,  2.607119_wp, 1.940513_wp, 1.650000_wp, &
      & 2.636374_wp, 2.145197_wp, 1.650000_wp,  2.622072_wp, 2.341428_wp, 1.650000_wp, &
      & 3.093746_wp, 2.050985_wp, 1.750000_wp,  2.012351_wp, 1.773458_wp, 0.000000_wp, &
      & 2.607880_wp, 1.680493_wp, 0.000000_wp,  2.817165_wp, 1.919947_wp, 0.000000_wp, &
      & 2.998330_wp, 2.249513_wp, 0.000000_wp,  3.150662_wp, 2.382063_wp, 1.241625_wp, &
      & 3.516922_wp, 2.392024_wp, 1.380239_wp,  3.520683_wp, 2.535389_wp, 1.418875_wp],&
      & shape(slater_exponent))

   !> Atomic level information
   real(wp), parameter :: p_selfenergy(max_shell, max_elem) = evtoau * reshape([&
      &-10.521941_wp, -3.560151_wp,  0.000000_wp, -22.121015_wp,  0.000000_wp,  0.000000_wp, &
      & -6.778577_wp, -3.316533_wp,  0.000000_wp,  -9.155548_wp, -6.934385_wp,  0.000000_wp, &
      &-12.479504_wp, -7.044782_wp, -4.032136_wp, -15.968100_wp, -9.777107_wp, -3.366990_wp, &
      &-24.040776_wp,-14.426520_wp, -5.882278_wp, -26.695120_wp,-16.014457_wp, -6.622058_wp, &
      &-27.565998_wp,-17.439238_wp, -6.120500_wp, -31.167487_wp,-18.268975_wp,  1.487984_wp, &
      & -5.425941_wp, -2.970587_wp,  0.000000_wp, -10.585393_wp, -2.983087_wp,  0.000000_wp, &
      &-12.916245_wp, -3.441043_wp, -1.751415_wp, -14.506128_wp, -7.557337_wp, -2.508113_wp, &
      &-16.954206_wp, -9.866726_wp, -0.997152_wp, -24.488821_wp,-12.136640_wp, -1.796070_wp, &
      &-24.452163_wp,-13.465451_wp, -1.836065_wp, -31.395427_wp,-17.412901_wp, -1.119399_wp, &
      & -5.815562_wp, -3.747255_wp,  0.000000_wp,  -7.979180_wp, -2.517008_wp, -2.752355_wp, &
      & -6.601089_wp, -3.541803_wp, -1.687209_wp,  -6.171156_wp, -4.647280_wp, -0.606001_wp, &
      & -6.314200_wp,-12.834159_wp,  0.035792_wp,  -8.095169_wp, -6.662737_wp, -2.647661_wp, &
      & -6.756964_wp, -6.497660_wp, -4.391631_wp,  -8.609628_wp, -6.896217_wp, -6.443170_wp, &
      & -8.957293_wp, -5.330498_wp, -3.914720_wp,  -9.491654_wp, -6.891146_wp, -3.002778_wp, &
      &-11.157415_wp, -8.659820_wp, -4.486136_wp,  -9.538781_wp, -4.393779_wp,  0.000000_wp, &
      &-14.013394_wp, -5.616872_wp, -2.437453_wp, -15.795287_wp,-10.050005_wp, -4.221097_wp, &
      &-21.397516_wp, -9.716689_wp, -2.252045_wp, -23.000000_wp,-11.353754_wp, -0.978922_wp, &
      &-19.875752_wp,-12.818655_wp, -3.348113_wp, -20.280017_wp,-15.200155_wp, -4.253986_wp, &
      & -7.616948_wp, -4.369842_wp,  0.000000_wp,  -6.840171_wp, -3.338573_wp, -1.715680_wp, &
      & -5.731066_wp, -8.748292_wp, -0.838555_wp,  -5.698303_wp, -2.340849_wp, -5.079938_wp, &
      & -7.621949_wp, -6.455094_wp, -2.522398_wp,  -7.338804_wp, -7.762642_wp, -3.971551_wp, &
      & -8.690050_wp, -5.089073_wp, -4.878724_wp, -10.960165_wp, -6.304229_wp, -5.569969_wp, &
      &-11.935915_wp, -4.883179_wp, -4.427854_wp, -11.057352_wp, -5.667970_wp, -3.149700_wp, &
      & -9.967352_wp, -5.974518_wp, -3.478715_wp,  -9.271670_wp, -3.857452_wp,  0.000000_wp, &
      &-14.755266_wp, -4.611036_wp, -4.455688_wp, -22.241945_wp, -7.801265_wp, -3.306963_wp, &
      &-19.959912_wp, -7.703762_wp, -2.025267_wp, -27.798530_wp,-10.496142_wp, -1.080464_wp, &
      &-23.832631_wp,-11.604442_wp, -2.025327_wp, -21.969064_wp,-11.870978_wp, -2.697796_wp, &
      & -6.341379_wp, -3.944275_wp,  0.000000_wp,  -6.452630_wp, -3.975353_wp, -2.305768_wp, &
      & -5.872226_wp, -6.500000_wp, -0.727921_wp,  -5.032003_wp, -6.275363_wp,  0.291196_wp, &
      & -4.944984_wp, -6.271128_wp,  0.241817_wp,  -4.857964_wp, -6.266893_wp,  0.192438_wp, &
      & -4.770945_wp, -6.262657_wp,  0.143059_wp,  -4.683925_wp, -6.258422_wp,  0.093680_wp, &
      & -4.596906_wp, -6.254187_wp,  0.044301_wp,  -4.509886_wp, -6.249952_wp, -0.005078_wp, &
      & -4.422867_wp, -6.245716_wp, -0.054457_wp,  -4.335848_wp, -6.241481_wp, -0.103836_wp, &
      & -4.248828_wp, -6.237246_wp, -0.153215_wp,  -4.161809_wp, -6.233011_wp, -0.202593_wp, &
      & -4.074789_wp, -6.228775_wp, -0.251972_wp,  -3.987770_wp, -6.224540_wp, -0.301351_wp, &
      & -3.900750_wp, -6.220305_wp, -0.350730_wp,  -3.817190_wp, -4.131659_wp, -2.759205_wp, &
      & -9.232014_wp, -8.600553_wp, -1.266176_wp,  -7.943440_wp, -2.878936_wp, -2.606135_wp, &
      & -7.620328_wp, -8.366866_wp, -7.118201_wp, -10.102898_wp, -3.655133_wp, -7.060522_wp, &
      &-11.130049_wp, -5.401042_wp, -6.032947_wp, -11.750154_wp, -5.493306_wp, -4.896782_wp, &
      &-10.306258_wp, -7.467067_wp, -3.211351_wp, -11.634473_wp, -4.320883_wp,  0.000000_wp, &
      &-10.454936_wp, -5.031237_wp,  0.000000_wp,  -8.264605_wp, -6.658984_wp,  0.000000_wp, &
      &-18.199529_wp, -8.408056_wp,  0.000000_wp, -23.908422_wp, -8.889548_wp, -0.921251_wp, &
      &-21.752193_wp,-10.031093_wp, -0.852571_wp, -18.381647_wp,-10.236606_wp, -0.973687_wp],&
      & shape(p_selfenergy))

   integer, parameter :: ipea1_kinds(max_elem) = [&
      &  1,                                                 1, &! H-He
      &  0, 0,                               0, 1, 1, 1, 1, 1, &! Li-Ne
      &  0, 0,                               0, 1, 1, 1, 1, 1, &! Na-Ar
      &  0, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, &! K-Kr
      &  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, &! Rb-Xe
      &  0, 0, &! Cs/Ba
      &        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &!La-Lu
      &        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]  ! Lu-Rn

   !> Specification of the IPEA1-xTB effective Hamiltonian
   type, extends(tb_h0spec) :: ipea1_h0spec
      real(wp) :: kshell(0:2, 0:2)
      real(wp), allocatable :: kpair(:, :)
      logical, allocatable :: valence(:, :)
   contains
      !> Generator for the self energy / atomic levels of the Hamiltonian
      procedure :: get_selfenergy
      !> Generator for the coordination number dependent shift of the self energy
      procedure :: get_cnshift
      !> Generator for the enhancement factor to for scaling Hamiltonian elements
      procedure :: get_hscale
      !> Generator for the polynomial parameters for the distant dependent scaling
      procedure :: get_shpoly
      !> Generator for the reference occupation numbers of the atoms
      procedure :: get_reference_occ
   end type ipea1_h0spec

   interface ipea1_h0spec
      module procedure :: new_ipea1_h0spec
   end interface ipea1_h0spec

contains


subroutine new_ipea1_calculator(calc, mol)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(out) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   call add_basis(calc, mol)
   call add_ncoord(calc, mol)
   call add_hamiltonian(calc, mol)
   call add_repulsion(calc, mol)
   call add_dispersion(calc, mol)
   call add_coulomb(calc, mol)

end subroutine new_ipea1_calculator

subroutine add_basis(calc, mol)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   integer :: isp, izp, ish, stat, ng, il
   integer, allocatable :: nsh_id(:)
   integer :: ang_idx(0:2), ortho(max_shell)
   type(cgto_type), allocatable :: cgto(:, :)

   nsh_id = nshell(mol%num)
   allocate(cgto(maxval(nsh_id), mol%nid))
   do isp = 1, mol%nid
      ang_idx = 0
      ortho = 0
      izp = mol%num(isp)
      do ish = 1, nsh_id(isp)
         il = ang_shell(ish, izp)
         ng = number_of_primitives(ish, izp, ang_idx(il) == 0)
         if (ang_idx(il) > 0) then
            ortho(ish) = ang_idx(il)
         else
            ang_idx(il) = ish
         end if
         call slater_to_gauss(ng, principal_quantum_number(ish, izp), il, &
            & slater_exponent(ish, izp), cgto(ish, isp), .true., stat)
      end do

      do ish = 1, nsh_id(isp)
         if (ortho(ish) > 0) then
            call orthogonalize(cgto(ortho(ish), isp), cgto(ish, isp))
         end if
      end do
   end do

   call new_basis(calc%bas, mol, nsh_id, cgto, 1.0_wp)

end subroutine add_basis

pure function number_of_primitives(ish, izp, valence) result(nprim)
   integer, intent(in) :: ish
   integer, intent(in) :: izp
   logical, intent(in) :: valence
   integer :: nprim

   nprim = 0
   if (izp <= 2) then
      select case(ang_shell(ish, izp))
      case(0)
         nprim = merge(4, 3, valence)
      case(1:)
         nprim = 3
      end select
   else
      select case(ang_shell(ish, izp))
      case(0)
         nprim = merge(6, 3, principal_quantum_number(ish, izp) > 5 .or. valence)
      case(1)
         nprim = 6
      case(2:)
         nprim = 4
      end select
   end if

end function number_of_primitives

subroutine add_hamiltonian(calc, mol)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   call new_hamiltonian(calc%h0, mol, calc%bas, new_ipea1_h0spec(mol))
end subroutine add_hamiltonian

subroutine add_dispersion(calc, mol)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   real(wp), parameter :: s6 = 1.0_wp, s8 = 2.4_wp, a1 = 0.63_wp, a2 = 5.0_wp, s9 = 0.0_wp
   type(d3_dispersion), allocatable :: tmp

   allocate(tmp)
   call new_d3_dispersion(tmp, mol, s6=s6, s8=s8, a1=a1, a2=a2, s9=s9)
   call move_alloc(tmp, calc%dispersion)
end subroutine add_dispersion

subroutine add_ncoord(calc, mol)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   call new_ncoord(calc%ncoord, mol, cn_type="exp")
end subroutine add_ncoord

subroutine add_repulsion(calc, mol)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   real(wp), allocatable :: alpha(:), zeff(:)

   allocate(calc%repulsion)
   alpha = rep_alpha(mol%num)
   zeff = rep_zeff(mol%num)
   call new_repulsion(calc%repulsion, mol, alpha, zeff, 1.5_wp, 1.5_wp, 1.0_wp)
end subroutine add_repulsion

subroutine add_coulomb(calc, mol)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   real(wp), allocatable :: hardness(:, :), hubbard_derivs(:, :)

   allocate(calc%coulomb)
   allocate(calc%coulomb%es2)
   call get_shell_hardness(mol, calc%bas, hardness)
   call new_effective_coulomb(calc%coulomb%es2, mol, 1.0_wp, hardness, harmonic_average, &
      & calc%bas%nsh_id)

   allocate(calc%coulomb%es3)
   hubbard_derivs = spread(p_hubbard_derivs(mol%num), 1, 1)
   call new_onsite_thirdorder(calc%coulomb%es3, mol, hubbard_derivs)

end subroutine add_coulomb

subroutine get_shell_hardness(mol, bas, hardness)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Shell resolved hardness parameters
   real(wp), allocatable, intent(out) :: hardness(:, :)

   integer :: isp, izp, ish, il

   allocate(hardness(maxval(bas%nsh_id), mol%nid))
   hardness(:, :) = 0.0_wp
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         il = bas%cgto(ish, isp)%ang
         hardness(ish, isp) = chemical_hardness(izp) * shell_hardness(il, izp)
      end do
   end do
end subroutine get_shell_hardness


pure function new_ipea1_h0spec(mol) result(self)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Instance of the Hamiltonian specification
   type(ipea1_h0spec) :: self

   real(wp), parameter :: kshell(0:2) = [2.15_wp, 2.30_wp, 2.0_wp]
   integer :: isp, jsp, il, jl, izp, ish
   integer :: ang_idx(0:2)

   allocate(self%kpair(mol%nid, mol%nid))
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         self%kpair(jsp, isp) = get_pair_param(mol%num(jsp), mol%num(isp))
      end do
   end do

   do il = 0, 2
      do jl = 0, 2
         self%kshell(jl, il) = 0.5_wp * (kshell(jl) + kshell(il))
      end do
   end do
   self%kshell(0, 1) = 2.25_wp
   self%kshell(1, 0) = 2.25_wp

   allocate(self%valence(3, mol%nid))
   do isp = 1, mol%nid
      ang_idx = 0
      izp = mol%num(isp)
      do ish = 1, nshell(izp)
         il = ang_shell(ish, izp)
         self%valence(ish, isp) = ang_idx(il) == 0
         if (self%valence(ish, isp)) ang_idx(il) = ish
      end do
   end do

end function new_ipea1_h0spec

pure function get_pair_param(jzp, izp) result(kpair)
   integer, intent(in) :: izp, jzp
   real(wp) :: kpair
   integer :: itr, jtr
   real(wp), parameter :: kp(3) = [1.1_wp, 1.2_wp, 1.2_wp]

   itr = get_dblock_row(izp)
   jtr = get_dblock_row(jzp)
   if (itr > 0 .and. jtr > 0) then
      kpair = 0.5_wp * (kp(itr) + kp(jtr))
   else
      kpair = 1.0_wp
   end if
end function get_pair_param

elemental function get_dblock_row(zp) result(tr)
   integer, intent(in) :: zp
   integer :: tr

   if (zp > 20 .and. zp < 30) then
      tr = 1
   else if (zp > 38 .and. zp < 48) then
      tr = 2
   else if (zp > 56 .and. zp < 80) then
      tr = 3
   else
      tr = 0
   end if

end function get_dblock_row

!> Generator for the enhancement factor to for scaling Hamiltonian elements
subroutine get_hscale(self, mol, bas, hscale)
   !> Instance of the Hamiltonian specification
   class(ipea1_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Scaling parameters for the Hamiltonian elements
   real(wp), intent(out) :: hscale(:, :, :, :)

   real(wp), parameter :: enscale = -7.0e-3_wp
   real(wp), parameter :: kdiff = 2.25_wp
   integer :: isp, jsp, izp, jzp, ish, jsh, il, jl
   real(wp) :: den, enp, km

   hscale(:, :, :, :) = 0.0_wp

   do isp = 1, mol%nid
      izp = mol%num(isp)
      do jsp = 1, mol%nid
         jzp = mol%num(jsp)
         den = (get_pauling_en(izp) - get_pauling_en(jzp))**2
         do ish = 1, bas%nsh_id(isp)
            il = bas%cgto(ish, isp)%ang
            do jsh = 1, bas%nsh_id(jsp)
               jl = bas%cgto(jsh, jsp)%ang
               if (self%valence(ish, isp) .and. self%valence(jsh, jsp)) then
                  enp = 1.0_wp + enscale * den
                  km = self%kpair(jsp, isp) * self%kshell(jl, il) * enp
               else if (self%valence(ish, isp) .and. .not.self%valence(jsh, jsp)) then
                  km = 0.5_wp * (self%kshell(il, il) + kdiff)
               else if (.not.self%valence(ish, isp) .and. self%valence(jsh, jsp)) then
                  km = 0.5_wp * (self%kshell(jl, jl) + kdiff)
               else
                  km = kdiff
               end if
               hscale(jsh, ish, jsp, isp) = km
            end do
         end do
      end do
   end do
end subroutine get_hscale


!> Generator for the self energy / atomic levels of the Hamiltonian
subroutine get_selfenergy(self, mol, bas, selfenergy)
   !> Instance of the Hamiltonian specification
   class(ipea1_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Self energy / atomic levels
   real(wp), intent(out) :: selfenergy(:, :)

   integer :: isp, izp, ish

   selfenergy(:, :) = 0.0_wp

   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         selfenergy(ish, isp) = p_selfenergy(ish, izp)
      end do
   end do
end subroutine get_selfenergy


!> Generator of the coordination number dependent shift of the self energy
subroutine get_cnshift(self, mol, bas, kcn)
   !> Instance of the Hamiltonian specification
   class(ipea1_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Coordination number dependent shift
   real(wp), intent(out) :: kcn(:, :)

   integer :: isp, izp, ish, il, ik

   real(wp), parameter :: cnshell(2, 0:2) = reshape(&
      &[1.0_wp, 1.0_wp,-0.3_wp,-0.3_wp,-0.5_wp, 0.5_wp], &
      & shape(cnshell)) * 0.01_wp

   kcn(:, :) = 0.0_wp
   do isp = 1, mol%nid
      izp = mol%num(isp)
      ik = ipea1_kinds(izp)
      if (ik > 0) then
         do ish = 1, bas%nsh_id(isp)
            il = bas%cgto(ish, isp)%ang
            kcn(ish, isp) = - p_selfenergy(ish, izp) * cnshell(ik, il)
         end do
      end if
   end do
end subroutine get_cnshift


!> Generator for the polynomial parameters for the distant dependent scaling
subroutine get_shpoly(self, mol, bas, shpoly)
   !> Instance of the Hamiltonian specification
   class(ipea1_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Polynomial parameters for distant dependent scaleing
   real(wp), intent(out) :: shpoly(:, :)

   integer :: isp, izp, ish

   shpoly(:, :) = 0.0_wp

   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         shpoly(ish, isp) = p_shpoly(bas%cgto(ish, isp)%ang, izp)
      end do
   end do
end subroutine get_shpoly


subroutine get_reference_occ(self, mol, bas, refocc)
   !> Instance of the Hamiltonian specification
   class(ipea1_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Reference occupation numbers
   real(wp), intent(out) :: refocc(:, :)

   integer :: isp, izp, ish

   refocc(:, :) = 0.0_wp

   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         if (self%valence(ish, isp)) then
            refocc(ish, isp) = reference_occ(bas%cgto(ish, isp)%ang, izp)
         else
            refocc(ish, isp) = 0.0_wp
         end if
      end do
   end do
end subroutine get_reference_occ


subroutine get_hubbard_derivs(mol, bas, hubbard_derivs)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Shell resolved Hubbard derivatives
   real(wp), allocatable, intent(out) :: hubbard_derivs(:, :)

   real(wp), parameter :: shell_hubbard_derivs(0:2) = [1.0_wp, 0.5_wp, 0.25_wp]

   integer :: isp, izp, ish, il

   allocate(hubbard_derivs(maxval(bas%nsh_id), mol%nid))
   hubbard_derivs(:, :) = 0.0_wp
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, bas%nsh_id(isp)
         il = bas%cgto(ish, isp)%ang
         hubbard_derivs(ish, isp) = p_hubbard_derivs(izp) * shell_hubbard_derivs(il)
      end do
   end do
end subroutine get_hubbard_derivs

end module tblite_xtb_ipea1