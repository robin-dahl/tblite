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

module test_solvation_draco
   use mctc_env, only: wp
   use mctc_io_convert, only: autoaa, aatoau
   use mctc_env_testing, only: new_unittest, unittest_type, error_type, check, &
     & test_failed
   use mctc_io, only: structure_type
   use mstore, only: get_structure
   use tblite_container, only: container_cache
   use tblite_scf_potential, only: potential_type
   use tblite_solvation_alpb
   use tblite_solvation_born
   use tblite_solvation_data, only: solvent_data, get_solvent_data, &
      & get_vdw_rad_cosmo, get_vdw_rad_bondi, get_vdw_rad_d3
   use tblite_solvation_draco, only: new_draco, draco_type
   use tblite_solvation_data_alpb, only: get_alpb_param
   use tblite_wavefunction, only: wavefunction_type, new_wavefunction, &
      &  eeq_guess
   use tblite_xtb_calculator, only: xtb_calculator
   use tblite_xtb_gfn2, only: new_gfn2_calculator
   use tblite_ncoord, only: new_ncoord, ncoord_type
   use tblite_disp_d4, only: get_eeq_charges
   implicit none
   private
   public :: collect_solvation_draco

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains
   !> Collect all exported unit tests
   subroutine collect_solvation_draco(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("cpcm-radii-adjust", test_r_cpcm_adjust), &
                  new_unittest("smd-radii-adjust", test_r_smd_adjust), &
                  new_unittest("cosmo-radii-adjust", test_r_cosmo_adjust), &
                  new_unittest("scaled_radii_gradient", test_gradient_cpcm) &
                  ]
   end subroutine collect_solvation_draco

   !> Test the dynamic radii adjustment
   subroutine test_r_adjust(error, mol, draco, ref)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> DRACO class
      type(draco_type), intent(in) :: draco

      !> Molecular structure data
      type(structure_type), intent(in) :: mol

      !> Reference energy
      real(wp), intent(in) :: ref(:)

      !> Coordination numbers
      real(wp), allocatable :: cn(:)

      !> Partial charges
      real(wp), allocatable :: qat(:)

      !> Radii adjusted with scaling function
      real(wp), allocatable :: scaled_radii(:)

      
      class(ncoord_type), allocatable :: ncoord

      allocate (cn(mol%nat), qat(mol%nat), scaled_radii(mol%nat))

      ! Get eeq charges and coordiantion numbers
      call get_eeq_charges(mol, qat)
      call new_ncoord(ncoord, mol, cn_type="gfn")
      call ncoord%get_cn(mol, cn)

      ! Use DRACO to sclae the radii
      call draco%radii_adjustment(mol, qat, cn, scaled_radii)

      print '(3es21.13)', scaled_radii

      if (any(abs(scaled_radii - ref) > thr2)) then
          call test_failed(error, "Scaled radii do not match reference")
          print '(a)', 'Scaled radii:'
          print '(3es21.13)', scaled_radii
          print '(a)', 'Difference to ref:'
          print '(3es21.13)', scaled_radii - ref
      end if
   end subroutine test_r_adjust

   !> Test the analytical gradient of the scaled radii w.r.t. the atom positions
   subroutine test_g(error, mol, draco) 
      !> Error handling

      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type), intent(inout) :: mol
      
      !> DRACO class
      type(draco_type), intent(in) :: draco

      real(wp), parameter :: step = 1.0e-4_wp
      real(wp), allocatable :: dsraddr(:, :, :), dsraddL(:, :, :), numgr(:, :, :), numgL(:, :, :)
      real(wp), parameter :: unity(3, 3)  = reshape(&
         & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
      real(wp), allocatable :: dqdr(:, :, :), dcndr(:, :, :), dqdL(:, :, :), dcndL(:, :, :)
      real(wp), allocatable :: qat(:), cn(:), scaled_radii(:)
      real(wp), allocatable :: rr(:), rl(:), lr(:), ll(:), eps(:, :)
      
      integer :: ii, ic
      real(wp), allocatable :: dpat(:, :)
      class(ncoord_type), allocatable :: ncoord
      call new_ncoord(ncoord, mol, cn_type="gfn")

      allocate(cn(mol%nat), qat(mol%nat), scaled_radii(mol%nat), rr(mol%nat), rl(mol%nat))
      allocate(numgr(3, mol%nat, mol%nat))
      allocate(dsraddr(3, mol%nat, mol%nat), dsraddL(3, 3, mol%nat), dpat(3, mol%nat))
      do ii = 1, mol%nat
         do ic = 1, 3
            rr = 0.0_wp
            rl = 0.0_wp

            ! Deflection to the right
            mol%xyz(ic, ii) = mol%xyz(ic, ii) + step
            call get_eeq_charges(mol, qat)
            call ncoord%get_cn(mol, cn)
            call draco%radii_adjustment(mol, qat, cn, rr)

            ! Deflection to the left
            mol%xyz(ic, ii) = mol%xyz(ic, ii) - 2*step
            call get_eeq_charges(mol, qat)
            call ncoord%get_cn(mol, cn)
            call draco%radii_adjustment(mol, qat, cn, rl)

            ! Back to initial structure
            mol%xyz(ic, ii) = mol%xyz(ic, ii) + step 

            ! Calcualtion of numerical gradient
            numgr(ic, ii, :) = 0.5_wp*(rr - rl)/step
         end do
      end do
      
      allocate(numgL(3, 3, mol%nat), eps(3, 3))
      eps(:, :) = unity
      do ic = 1, 3
         do ii = 1, 3
            eps(ii, ic) = eps(ii, ic) + step
            mol%xyz(:, :) = matmul(eps, mol%xyz)
            call get_eeq_charges(mol, qat)
            call ncoord%get_cn(mol, cn)
            if (any(mol%periodic)) mol%lattice(:, :) = matmul(eps, mol%lattice)
            call draco%radii_adjustment(mol, qat, cn, lr)

            eps(ii, ic) = eps(ii, ic) - 2*step
            mol%xyz(:, :) = matmul(eps, mol%xyz)
            call get_eeq_charges(mol, qat)
            call ncoord%get_cn(mol, cn)
            if (any(mol%periodic)) mol%lattice(:, :) = matmul(eps, mol%lattice)
            call draco%radii_adjustment(mol, qat, cn, ll)

            numgL(ii, ic, :) = 0.5_wp*(lr - ll)/step

            eps(ii, ic) = eps(ii, ic) + step
            mol%xyz(:, :) = matmul(eps, mol%xyz)
            if (any(mol%periodic)) mol%lattice(:, :) = matmul(eps, mol%lattice)
         end do
      end do
      !print '(3es20.13)', numgL
        
      allocate (dqdr(3, mol%nat, mol%nat), dcndr(3, mol%nat, mol%nat))
      allocate (dqdL(3, 3, mol%nat), dcndL(3, 3, mol%nat))

      
      call get_eeq_charges(mol, qat, dqdr=dqdr, dqdL=dqdL)
      call ncoord%get_cn(mol, cn, dcndr=dcndr, dcndL=dcndL)
      ! Get the analytical gradients
      call draco%radii_adjustment(mol, qat, cn, scaled_radii, dsraddr, dsraddL, dqdr, dcndr, dqdL, dcndL)

      if (any(abs(dsraddr - numgr) > thr2)) then
         call test_failed(error, "Gradient w.r.t r does not match")
         print '(a)', 'dsraddr'
         print '(3es20.13)', dsraddr
         print '(a)', "numgr"
         print '(3es20.13)', numgr
         print '(a)', "dsraddr - numgr"
         print '(3es20.13)', dsraddr - numgr
      end if

      ! if (any(abs(dsraddL - numgL) > thr2)) then
      !    call test_failed(error, "Gradient w.r.t L does not match")
      !    print '(a)', 'dsraddL:'
      !    print '(3es20.13)', dsraddL
      !    print '(a)', "numgL:"
      !    print '(3es20.13)', numgL
      !    print '(a)', "dsraddL - numgL:"
      !    print '(3es20.13)', dsraddL - numgL
      ! end if
   end subroutine test_g

   subroutine test_r_cpcm_adjust(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type) :: mol

      !> DRACO class
      type(draco_type) :: draco

      !> Reference
      real(wp), parameter :: ref(16) = [&
      &1.5668219785499E+00_wp, 2.3040000000000E+00_wp, 1.6865825254874E+00_wp,&
      &2.6652696084507E+00_wp, 2.3040000000000E+00_wp, 1.5269659439116E+00_wp,&
      &1.5532212184893E+00_wp, 2.5200000000000E+00_wp, 1.6225088275790E+00_wp,&
      &1.5309489956972E+00_wp, 1.4345096465537E+00_wp, 2.3040000000000E+00_wp,&
      &2.5200000000000E+00_wp, 1.3089371798953E+00_wp, 1.5819517127436E+00_wp,&
      &2.3040000000000E+00_wp]

      call get_structure(mol, "MB16-43", "04")
      call new_draco(draco, mol, 'cpcm')
      call test_r_adjust(error, mol, draco, ref)
   end subroutine test_r_cpcm_adjust

   subroutine test_r_smd_adjust(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type) :: mol

      !> DRACO class
      type(draco_type) :: draco

      !> Reference
      real(wp), parameter :: ref(16) = [&
      &1.2390428427580E+00_wp, 1.9200000000000E+00_wp, 1.2615521116284E+00_wp, &
      &1.7560521433416E+00_wp, 1.9200000000000E+00_wp, 1.2304816656895E+00_wp, &
      &1.2315634651943E+00_wp, 2.4700000000000E+00_wp, 1.2463742450074E+00_wp, &
      &1.2248897626319E+00_wp, 1.7072333614339E+00_wp, 1.8400000000000E+00_wp, &
      &2.4700000000000E+00_wp, 1.3118205039369E+00_wp, 1.2381437933040E+00_wp, &
      &1.9200000000000E+00_wp]

      call get_structure(mol, "MB16-43", "04")
      call new_draco(draco, mol, 'smd')
      call test_r_adjust(error, mol, draco, ref)
   end subroutine test_r_smd_adjust

   subroutine test_r_cosmo_adjust(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type) :: mol

      !> DRACO class
      type(draco_type) :: draco

      !> Reference
      real(wp), parameter :: ref(16) = [&
      &1.4239147358894E+00_wp, 2.0480000000000E+00_wp, 1.5005239643615E+00_wp, &
      &2.6388987366261E+00_wp, 2.0480000000000E+00_wp, 1.3988305709959E+00_wp, &
      &1.4155564697294E+00_wp, 2.2000000000000E+00_wp, 1.4594599743493E+00_wp, &
      &1.4016403136306E+00_wp, 1.5569074248726E+00_wp, 2.1530000000000E+00_wp, &
      &2.2000000000000E+00_wp, 1.2595826301673E+00_wp, 1.4336763841499E+00_wp, &
      &2.0480000000000E+00_wp]

      call get_structure(mol, "MB16-43", "04")
      call new_draco(draco, mol, 'cosmo')
      call test_r_adjust(error, mol, draco, ref)
   end subroutine test_r_cosmo_adjust

   subroutine test_gradient_cpcm(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type) :: mol

      !> DRACO class
      type(draco_type) :: draco

      call get_structure(mol, "MB16-43", "04")
      call new_draco(draco, mol, 'cpcm')
      call test_g(error, mol, draco)
   end subroutine test_gradient_cpcm
end module test_solvation_draco
