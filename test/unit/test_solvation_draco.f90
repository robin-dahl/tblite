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
    use mctc_env, only : wp
    use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
       & test_failed
    use mctc_io, only : structure_type
    use mstore, only : get_structure
    use tblite_container, only : container_cache
    use tblite_scf_potential, only : potential_type
    use tblite_solvation_alpb
    use tblite_solvation_born
    use tblite_solvation_data, only : solvent_data, get_solvent_data, &
       & get_vdw_rad_cosmo, get_vdw_rad_bondi, get_vdw_rad_d3
    use tblite_solvation_draco, only: new_draco, draco_type
    use tblite_solvation_data_alpb, only: get_alpb_param
    use tblite_wavefunction, only : wavefunction_type, new_wavefunction, &
       &  eeq_guess
    use tblite_xtb_calculator, only : xtb_calculator
    use tblite_xtb_gfn2, only : new_gfn2_calculator
    use tblite_ncoord, only: new_ncoord
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
       new_unittest("cosmo-radii-adjust", test_r_cosmo_adjust) &
       ]
 
 end subroutine collect_solvation_draco
 
 
 subroutine test_r_adjust(error, mol, draco, ref)
 
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    type(draco_type), intent(in) :: draco
 
    !> Molecular structure data
    type(structure_type), intent(in) :: mol

    !> Reference energy
    real(wp), intent(in) :: ref(:)

    type(xtb_calculator) :: calc
    type(wavefunction_type) :: wfn
    
    real(wp), allocatable :: cn(:)
    real(wp), allocatable :: scaled_radii(:)

    call new_gfn2_calculator(calc, mol, error) !dummy
    call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 0.0_wp)
    call eeq_guess(mol, calc, wfn)

    allocate(cn(mol%nat), scaled_radii(mol%nat))

    call new_ncoord(calc%ncoord, mol, "gfn")
    call calc%ncoord%get_cn(mol, cn)
    
    call draco%radii_adjustment(mol, wfn%qat(:,1), cn, scaled_radii)

 
    if (any(abs(scaled_radii - ref) > thr2)) then
       call test_failed(error, "Scaled radii do not match reference")
       print '(a)', 'Scaled radii:'
       print '(3es21.13)', scaled_radii
       print '(a)', 'Difference to ref:'
       print '(3es21.13)', scaled_radii-ref
    end if
 end subroutine test_r_adjust
 
 
 subroutine test_g(error, mol, input, qat, method)
 
    !> Error handling
    type(error_type), allocatable, intent(out) :: error
 
    !> Molecular structure data
    type(structure_type), intent(inout) :: mol
 
    !> Solvation model input
    type(alpb_input), intent(in) :: input
 
    !> Atomic partial charges
    real(wp), intent(in) :: qat(:)
 
    !> Method for parameter selection
    character(len=*), intent(in), optional :: method
 
    type(alpb_solvation) :: solv
    type(alpb_input), allocatable :: scratch_input
    type(wavefunction_type) :: wfn
    type(potential_type) :: pot
    type(container_cache) :: cache
    real(wp), parameter :: step = 1.0e-4_wp
    real(wp), allocatable :: gradient(:, :), numg(:, :)
    real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3)
    integer :: ii, ic
 
    wfn%qat = reshape(qat, [size(qat), 1])
    allocate(pot%vat(size(qat, 1), 1))
 
    scratch_input = input
 
    if (allocated(input%solvent) .and. present(method)) then   
       call get_alpb_param(scratch_input, mol, method, error)
       if(allocated(error)) then
          call test_failed(error, "No ALPB/GBSA parameters found for the method/solvent")
          return
       end if
    end if
 
    solv = alpb_solvation(mol, scratch_input, method)
 
    allocate(numg(3, mol%nat), gradient(3, mol%nat))
    do ii = 1, mol%nat
       do ic = 1, 3
          er = 0.0_wp
          el = 0.0_wp
          mol%xyz(ic, ii) = mol%xyz(ic, ii) + step
          call solv%update(mol, cache)
          call solv%get_potential(mol, cache, wfn, pot)
          call solv%get_energy(mol, cache, wfn, er)
 
          mol%xyz(ic, ii) = mol%xyz(ic, ii) - 2*step
          call solv%update(mol, cache)
          call solv%get_potential(mol, cache, wfn, pot)
          call solv%get_energy(mol, cache, wfn, el)
 
          mol%xyz(ic, ii) = mol%xyz(ic, ii) + step
          numg(ic, ii) = 0.5_wp*(sum(er) - sum(el))/step
       end do
    end do
 
    energy = 0.0_wp
    gradient(:, :) = 0.0_wp
 
    call solv%update(mol, cache)
    call solv%get_potential(mol, cache, wfn, pot)
    call solv%get_energy(mol, cache, wfn, energy)
    call solv%get_gradient(mol, cache, wfn, gradient, sigma)
 
    if (any(abs(gradient - numg) > thr2)) then
       call test_failed(error, "Gradient does not match")
       print '(3es20.13)', gradient
       print '(a)', "---"
       print '(3es20.13)', numg
       print '(a)', "---"
       print '(3es20.13)', gradient - numg
    end if
 end subroutine test_g
 
 
 
 subroutine test_r_cpcm_adjust(error)
 
    !> Error handling
    type(error_type), allocatable, intent(out) :: error
 
    type(structure_type) :: mol
    
    type(draco_type) :: draco

    real(wp), parameter :: ref(16) = [&
    &1.5668221239151E+00_wp, 2.3040002137584E+00_wp, 1.6865826819636E+00_wp, &
    &2.6652698557266E+00_wp, 2.3040002137584E+00_wp, 1.5269660855791E+00_wp, &
    &1.5532213625926E+00_wp, 2.5200002337982E+00_wp, 1.6225089781106E+00_wp, &
    &1.5309491377342E+00_wp, 1.4345097796433E+00_wp, 2.3040002137584E+00_wp, &
    &2.5200002337982E+00_wp, 1.3089373013346E+00_wp, 1.5819518595124E+00_wp, &
    &2.3040002137584E+00_wp ]

 
    call get_structure(mol, "MB16-43", "04")

    call new_draco(draco, mol, 'cpcm')

    call test_r_adjust(error, mol, draco, ref) 
 
 end subroutine test_r_cpcm_adjust
 

 subroutine test_r_smd_adjust(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    type(structure_type) :: mol
    
    type(draco_type) :: draco

    real(wp), parameter :: ref(16) = [&
    &1.2390429577128E+00_wp, 1.9200001781320E+00_wp, 1.2615522286715E+00_wp,&
    &1.7560523062630E+00_wp, 1.9200001781320E+00_wp, 1.2304817798500E+00_wp,&
    &1.2315635794551E+00_wp, 2.4700002291593E+00_wp, 1.2463743606423E+00_wp,&
    &1.2248898762735E+00_wp, 1.7072335198260E+00_wp, 1.8400001707098E+00_wp,&
    &2.4700002291593E+00_wp, 1.3118206256438E+00_wp, 1.2381439081753E+00_wp,&
    &1.9200001781320E+00_wp ]

    call get_structure(mol, "MB16-43", "04")

    call new_draco(draco, mol, 'smd')

    call test_r_adjust(error, mol, draco, ref)

 end subroutine test_r_smd_adjust


 subroutine test_r_cosmo_adjust(error)

    !> Error handling
    type(error_type), allocatable, intent(out) :: error

    type(structure_type) :: mol

    type(draco_type) :: draco

    real(wp), parameter :: ref(16) = [&
    &1.4239148679961E+00_wp, 2.0480001900074E+00_wp, 1.5005241035757E+00_wp, &
    &2.6388989814554E+00_wp, 2.0480001900074E+00_wp, 1.3988307007753E+00_wp, &
    &1.4155566010606E+00_wp, 2.2000002041095E+00_wp, 1.4594601097537E+00_wp, &
    &1.4016404436706E+00_wp, 1.5569075693179E+00_wp, 2.1530001997490E+00_wp, &
    &2.2000002041095E+00_wp, 1.2595827470277E+00_wp, 1.4336765171621E+00_wp, &
    2.0480001900074E+00_wp ]
    

    call get_structure(mol, "MB16-43", "04")

    call new_draco(draco, mol, 'cosmo')

    call test_r_adjust(error, mol, draco, ref)

end subroutine test_r_cosmo_adjust

 
 end module test_solvation_draco
 