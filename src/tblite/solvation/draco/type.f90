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

        
        contains
        procedure :: radii_adjustment
        procedure :: get_fscale
    end type draco_type

contains

    subroutine new_draco(self, mol, radtype) ! Gets the default radii based on solvation model

        class(draco_type), intent(out) :: self

        character(len=*), intent(in) :: radtype

        type(structure_type), intent(in) :: mol
        
       
        select case (radtype)
        case('cpcm')
            self%defaultradii = get_vdw_rad_cpcm(mol%num)
            self%kcn = get_kcn_cpcm(mol%num)
            self%alpha = get_alpha_cpcm(mol%num)
            self%shift = get_shift_cpcm(mol%num)

        case ('smd')
            self%defaultradii = get_vdw_rad_smd(mol%num)
            self%kcn = get_kcn_smd(mol%num)
            self%alpha = get_alpha_smd(mol%num)
            self%shift = get_shift_smd(mol%num)

        case ('cosmo')
            self%defaultradii = get_vdw_rad_cosmo(mol%num)
            self%kcn = get_kcn_cosmo(mol%num)
            self%alpha = get_alpha_cosmo(mol%num)
            self%shift = get_shift_cosmo(mol%num)
        end select

    end subroutine new_draco
    
    !> Subroutine that rescales the radii, for now implemented only for water as a solvent
    subroutine radii_adjustment(self, mol, q, cn, scaled_radii) 
        
        class(draco_type), intent(in) :: self
        
        type(structure_type), intent(in) :: mol
        
        real(wp), intent(in), dimension(:) :: q

        real(wp), intent(in), dimension(:) :: cn
        
        real(wp), dimension(mol%nat), intent(out) :: scaled_radii 

        real(wp), allocatable :: fscale(:)
        integer :: iat, izp

        !>  convert bohr (a.u.) to Ångström and back
        real(wp), parameter :: autoaa = 0.52917726_wp
  
        allocate(fscale(mol%nat))

        call self%get_fscale(mol, q, cn, fscale)

        !if(present())

        ! Scale the radii
        do iat = 1, mol%nat
            izp = mol%id(iat)
            scaled_radii(iat) = self%defaultradii(izp) * fscale(iat) * autoaa 
        end do
 
    end subroutine radii_adjustment

    !> Subroutine that gets the scalin function for the radii based on effective partial charges and coordination numbers
    subroutine get_fscale(self, mol, q, cn, fscale) 
        
        class(draco_type), intent(in) :: self

        type(structure_type), intent(in) :: mol
        !> Partial charges
        real(wp), intent(in) :: q(:)
        !> Coordination numbers
        real(wp), intent(in) :: cn(:)
        !> Scaling function      
        real(wp), intent(out) :: fscale(:) 

        ! Effecvtive partial charges
        real(wp), allocatable :: qeff(:)
        ! General and specific indices for elements
        integer :: iat, izp, isp

        allocate(qeff(mol%nat))

        ! Scale the radii
        do iat = 1, mol%nat
            izp = mol%id(iat)
            
            qeff(iat) = q(iat) + self%kcn(izp) * q(iat) * cn(iat)

            fscale(iat) = erf( self%alpha(izp) * (qeff(iat) - self%shift(izp)) ) + 1

        end do

    end subroutine get_fscale


    ! subroutine get_dfscale(self, mol, q, cn, fscale)

    

    ! end subroutine get_dfscale 
    



end module tblite_solvation_draco_type