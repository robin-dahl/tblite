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

module test_solvation_ddx
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_basis_type, only : basis_type
   use tblite_container, only : container_cache
   use tblite_context_type, only : context_type
   use tblite_features, only : get_tblite_feature
   use tblite_integral_type, only : integral_type
   use tblite_scf_potential, only : potential_type
   use tblite_solvation_ddx, only : ddx_solvation, ddx_solvation_model, ddx_input, new_ddx
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, eeq_guess
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint

   implicit none
   private

   public :: collect_solvation_ddx

   real(wp), parameter :: acc = 0.01_wp
   real(wp), parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = 10*sqrt(epsilon(1.0_wp))
   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp


contains


!> Collect all exported unit tests
subroutine collect_solvation_ddx(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   if (get_tblite_feature("ddx")) then
      testsuite = [ &
         new_unittest("energy-mol-cosmo", test_e_cosmo_m01), &
         new_unittest("energy-mol-cpcm", test_e_cpcm_m01), &
         new_unittest("energy-mol-pcm", test_e_pcm_m01), &
         new_unittest("gradient-mol-num-cosmo", test_g_cosmo_m02), &
         new_unittest("gradient-mol-num-cosmo-dipole", test_g_cosmo_dipole_m02), &
         new_unittest("gradient-mol-num-cosmo-quadrupole", test_g_cosmo_quadrupole_m02), &
         new_unittest("gradient-mol-num-cpcm", test_g_cpcm_m02), &
         new_unittest("gradient-mol-num-pcm", test_g_pcm_m02), &
         new_unittest("potential-mol-cosmo", test_p_cosmo_m03), &
         new_unittest("potential-mol-cosmo-dipole", test_p_cosmo_dipole_m03), &
         new_unittest("potential-mol-cosmo-quadrupole", test_p_cosmo_quadrupole_m03), &
         new_unittest("potential-mol-cpcm", test_p_cpcm_m03), &
         new_unittest("potential-mol-pcm", test_p_pcm_m03) &
         ]
   else
      testsuite = [new_unittest("ddx-disabled", test_ddx_disabled)]
   end if

end subroutine collect_solvation_ddx

subroutine test_ddx_disabled(error)
   type(error_type), allocatable, intent(out) :: error

   call check(error, .not.get_tblite_feature("ddx"), "ddX feature is unexpectedly enabled")
end subroutine test_ddx_disabled


subroutine test_e(error, model, mol, qat, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Solvation model (COSMO=100, CPCM=101, PCM=200)
   integer, intent(in) :: model

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Reference energy
   real(wp), intent(in) :: ref

   type(ddx_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: eps = 80.0_wp
   integer, parameter :: nang = 302
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp) :: energy(mol%nat)

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1), source=0.0_wp)
   energy = 0.0_wp

   call new_ddx(solv, mol, ddx_input(ddx_model=model, dielectric_const=eps, nang=nang), error)
   if (allocated(error)) return

   call solv%update(mol, cache)
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)

   if (abs(sum(energy) - ref) > thr) then
      call test_failed(error, "Energy does not match reference")
      print '(a)', 'Energy:'
      print '(3es20.13)', sum(energy)
      print '(a)', "---"
      print '(a)', 'Reference:'
      print '(3es20.13)', ref
      print '(a)', "---"
      print '(a)', 'Difference:'
      print '(3es20.13)', sum(energy) - ref
   end if
end subroutine test_e


subroutine test_g(error, model, mol, qat, dpat, qpat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Solvation model (COSMO=100, CPCM=101, PCM=200)
   integer, intent(in) :: model

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Atom-resolved dipole moments
   real(wp), intent(in), optional :: dpat(:, :)

   !> Atom-resolved quadrupole moments
   real(wp), intent(in), optional :: qpat(:, :)

   type(ddx_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: eps = 80.0_wp
   integer, parameter :: nang = 302
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp), allocatable :: gradient(:, :), numg(:, :)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3)
   integer :: ii, ic

   wfn%qat = reshape(qat, [size(qat), 1])
   if (present(dpat)) then
      wfn%dpat = reshape(dpat, [3, mol%nat, 1])
   end if
   if (present(qpat)) then
      wfn%qpat = reshape(qpat, [6, mol%nat, 1])
   end if
   allocate(pot%vat(size(qat, 1), 1))
   if (present(dpat)) then
      allocate(pot%vdp(3, mol%nat, 1))
   end if
   if (present(qpat)) then
      allocate(pot%vqp(6, mol%nat, 1))
   end if

   call new_ddx(solv, mol, ddx_input(ddx_model=model, dielectric_const=eps, nang=nang, &
      & use_dipoles=present(dpat), use_quadrupoles=present(qpat)), error)
   if (allocated(error)) return

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

   if (any(abs(gradient - numg) > thr)) then
      call test_failed(error, "Gradient does not match")
      print '(3es20.13)', gradient
      print '(a)', "---"
      print '(3es20.13)', numg
      print '(a)', "---"
      print '(3es20.13)', gradient - numg
   end if
end subroutine test_g

subroutine test_p(error, model, mol, qat, dpat, qpat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Solvation model (COSMO=100, CPCM=101, PCM=200)
   integer, intent(in) :: model

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Atom-resolved dipole moments
   real(wp), intent(in), optional :: dpat(:, :)

   !> Atom-resolved quadrupole moments
   real(wp), intent(in), optional :: qpat(:, :)

   type(ddx_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache), allocatable :: cache
   real(wp), parameter :: eps = 80.0_wp
   integer, parameter :: nang = 302
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr = 1e+3_wp*sqrt(epsilon(1.0_wp))
   real(wp), allocatable :: vat(:), vdp(:, :), vqp(:, :)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat)
   integer :: ii, ic

   wfn%qat = reshape(qat, [size(qat), 1])
   if (present(dpat)) then
      wfn%dpat = reshape(dpat, [3, mol%nat, 1])
   end if
   if (present(qpat)) then
      wfn%qpat = reshape(qpat, [6, mol%nat, 1])
   end if
   allocate(pot%vat(size(qat, 1), 1))
   if (present(dpat)) then
      allocate(pot%vdp(3, mol%nat, 1))
   end if
   if (present(qpat)) then
      allocate(pot%vqp(6, mol%nat, 1))
   end if

   call new_ddx(solv, mol, ddx_input(ddx_model=model, dielectric_const=eps, nang=nang, &
      & use_dipoles=present(dpat), use_quadrupoles=present(qpat)), error)
   if (allocated(error)) return

   allocate(cache)
   call solv%update(mol, cache)

   allocate(vat(mol%nat))
   if (present(dpat)) then
      allocate(vdp(3, mol%nat))
   end if
   if (present(qpat)) then
      allocate(vqp(6, mol%nat))
   end if

   do ii = 1, mol%nat
      er = 0.0_wp
      el = 0.0_wp
      wfn%qat(ii, 1) = wfn%qat(ii, 1) + step
      call solv%get_potential(mol, cache, wfn, pot)
      call solv%get_energy(mol, cache, wfn, er)

      wfn%qat(ii, 1) = wfn%qat(ii, 1) - 2*step
      call solv%get_potential(mol, cache, wfn, pot)
      call solv%get_energy(mol, cache, wfn, el)

      wfn%qat(ii, 1) = wfn%qat(ii, 1) + step
      vat(ii) = 0.5_wp*(sum(er) - sum(el))/step
   end do

   if (present(dpat)) then
      do ii = 1, mol%nat
         do ic = 1, 3
            er = 0.0_wp
            el = 0.0_wp
            wfn%dpat(ic, ii, 1) = wfn%dpat(ic, ii, 1) + step
            call solv%get_potential(mol, cache, wfn, pot)
            call solv%get_energy(mol, cache, wfn, er)

            wfn%dpat(ic, ii, 1) = wfn%dpat(ic, ii, 1) - 2*step
            call solv%get_potential(mol, cache, wfn, pot)
            call solv%get_energy(mol, cache, wfn, el)

            wfn%dpat(ic, ii, 1) = wfn%dpat(ic, ii, 1) + step
            vdp(ic, ii) = 0.5_wp*(sum(er) - sum(el))/step
         end do
      end do
   end if

   if (present(qpat)) then
      do ii = 1, mol%nat
         do ic = 1, 6
            er = 0.0_wp
            el = 0.0_wp
            wfn%qpat(ic, ii, 1) = wfn%qpat(ic, ii, 1) + step
            call solv%get_potential(mol, cache, wfn, pot)
            call solv%get_energy(mol, cache, wfn, er)

            wfn%qpat(ic, ii, 1) = wfn%qpat(ic, ii, 1) - 2*step
            call solv%get_potential(mol, cache, wfn, pot)
            call solv%get_energy(mol, cache, wfn, el)

            wfn%qpat(ic, ii, 1) = wfn%qpat(ic, ii, 1) + step
            vqp(ic, ii) = 0.5_wp*(sum(er) - sum(el))/step
         end do
      end do
   end if

   energy = 0.0_wp
   pot%vat(:, :) = 0.0_wp
   if (allocated(pot%vdp)) pot%vdp(:, :, :) = 0.0_wp
   if (allocated(pot%vqp)) pot%vqp(:, :, :) = 0.0_wp
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)

   if (any(abs([pot%vat] - vat) > thr)) then
      call test_failed(error, "Charge-dependent potential does not match")
      print '(a)', 'analytical'
      print '(3es20.13)', pot%vat
      print '(a)', "---"
      print '(a)', 'numerical'
      print '(3es20.13)', vat
      print '(a)', "---"
      print '(a)', 'diff'
      print '(3es20.13)', [pot%vat] - vat
      return
   end if

   if (present(dpat)) then
      if (any(abs(pot%vdp(:, :, 1) - vdp) > thr)) then
         call test_failed(error, "Dipole-dependent potential does not match")
         print '(a)', 'analytical'
         print '(3es20.13)', pot%vdp(:, :, 1)
         print '(a)', "---"
         print '(a)', 'numerical'
         print '(3es20.13)', vdp
         print '(a)', "---"
         print '(a)', 'diff'
         print '(3es20.13)', pot%vdp(:, :, 1) - vdp
         return
      end if
   end if

   if (present(qpat)) then
      if (any(abs(pot%vqp(:, :, 1) - vqp) > thr)) then
         call test_failed(error, "Quadrupole-dependent potential does not match")
         print '(a)', 'analytical'
         print '(3es20.13)', pot%vqp(:, :, 1)
         print '(a)', "---"
         print '(a)', 'numerical'
         print '(3es20.13)', vqp
         print '(a)', "---"
         print '(a)', 'diff'
         print '(3es20.13)', pot%vqp(:, :, 1) - vqp
         return
      end if
   end if
end subroutine test_p


subroutine test_e_cosmo_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.73347900345264E-1_wp, 1.07626888948184E-1_wp,-3.66999593831010E-1_wp,&
      & 4.92833325937897E-2_wp,-1.83332156197733E-1_wp, 2.33302086605469E-1_wp,&
      & 6.61837152062315E-2_wp,-5.43944165050002E-1_wp,-2.70264356583716E-1_wp,&
      & 2.66618968841682E-1_wp, 2.62725033202480E-1_wp,-7.15315510172571E-2_wp,&
      &-3.73300777019193E-1_wp, 3.84585237785621E-2_wp,-5.05851088366940E-1_wp,&
      & 5.17677238544189E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   ! COSMO radii reference:
   ! call test_e(error, ddx_solvation_model%cosmo, mol, qat, -3.4697720884118800E-2_wp)
   call test_e(error, ddx_solvation_model%cosmo, mol, qat, -1.9188724249112E-2_wp)

end subroutine test_e_cosmo_m01

subroutine test_e_cpcm_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.73347900345264E-1_wp, 1.07626888948184E-1_wp,-3.66999593831010E-1_wp,&
      & 4.92833325937897E-2_wp,-1.83332156197733E-1_wp, 2.33302086605469E-1_wp,&
      & 6.61837152062315E-2_wp,-5.43944165050002E-1_wp,-2.70264356583716E-1_wp,&
      & 2.66618968841682E-1_wp, 2.62725033202480E-1_wp,-7.15315510172571E-2_wp,&
      &-3.73300777019193E-1_wp, 3.84585237785621E-2_wp,-5.05851088366940E-1_wp,&
      & 5.17677238544189E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   ! COSMO radii reference:
   ! call test_e(error, ddx_solvation_model%cpcm, mol, qat, -3.4914581639644553E-002_wp)
   call test_e(error, ddx_solvation_model%cpcm, mol, qat, -1.9308653775669E-2_wp)

end subroutine test_e_cpcm_m01

subroutine test_e_pcm_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.73347900345264E-1_wp, 1.07626888948184E-1_wp,-3.66999593831010E-1_wp,&
      & 4.92833325937897E-2_wp,-1.83332156197733E-1_wp, 2.33302086605469E-1_wp,&
      & 6.61837152062315E-2_wp,-5.43944165050002E-1_wp,-2.70264356583716E-1_wp,&
      & 2.66618968841682E-1_wp, 2.62725033202480E-1_wp,-7.15315510172571E-2_wp,&
      &-3.73300777019193E-1_wp, 3.84585237785621E-2_wp,-5.05851088366940E-1_wp,&
      & 5.17677238544189E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   ! COSMO radii reference:
   ! call test_e(error, ddx_solvation_model%pcm, mol, qat, -3.3624259293951506E-2_wp)
   call test_e(error, ddx_solvation_model%pcm, mol, qat, -1.8492346683564E-2_wp)

end subroutine test_e_pcm_m01

subroutine test_g_cosmo_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.50000000000000E-1_wp,-2.50000000000000E-1_wp, 5.00000000000000E-1_wp,&
      &-5.00000000000000E-1_wp]

   call get_structure(mol, "Heavy28", "bih3")
   call test_g(error, ddx_solvation_model%cosmo, mol, qat)

end subroutine test_g_cosmo_m02

subroutine test_g_cosmo_dipole_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.50000000000000E-1_wp,-2.50000000000000E-1_wp, 5.00000000000000E-1_wp,&
      &-5.00000000000000E-1_wp]
   real(wp), parameter :: dpat(3, 4) = reshape([&
      & 1.00000000000000E-1_wp,-2.00000000000000E-2_wp, 3.00000000000000E-2_wp,&
      &-4.00000000000000E-2_wp, 8.00000000000000E-2_wp,-1.00000000000000E-2_wp,&
      & 6.00000000000000E-2_wp, 2.00000000000000E-2_wp,-7.00000000000000E-2_wp,&
      &-3.00000000000000E-2_wp,-5.00000000000000E-2_wp, 4.00000000000000E-2_wp], &
      & [3, 4])

   call get_structure(mol, "Heavy28", "bih3")
   call test_g(error, ddx_solvation_model%cosmo, mol, qat, dpat=dpat)

end subroutine test_g_cosmo_dipole_m02

subroutine test_g_cosmo_quadrupole_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.50000000000000E-1_wp,-2.50000000000000E-1_wp, 5.00000000000000E-1_wp,&
      &-5.00000000000000E-1_wp]
   real(wp), parameter :: dpat(3, 4) = reshape([&
      & 1.00000000000000E-1_wp,-2.00000000000000E-2_wp, 3.00000000000000E-2_wp,&
      &-4.00000000000000E-2_wp, 8.00000000000000E-2_wp,-1.00000000000000E-2_wp,&
      & 6.00000000000000E-2_wp, 2.00000000000000E-2_wp,-7.00000000000000E-2_wp,&
      &-3.00000000000000E-2_wp,-5.00000000000000E-2_wp, 4.00000000000000E-2_wp], &
      & [3, 4])
   real(wp), parameter :: qpat(6, 4) = reshape([&
      & 2.00000000000000E-2_wp, 1.00000000000000E-2_wp,-3.00000000000000E-2_wp,&
      & 4.00000000000000E-3_wp,-2.00000000000000E-3_wp, 1.00000000000000E-2_wp,&
      &-1.00000000000000E-2_wp, 3.00000000000000E-3_wp, 2.00000000000000E-2_wp,&
      &-5.00000000000000E-3_wp, 7.00000000000000E-3_wp,-2.00000000000000E-2_wp,&
      & 3.00000000000000E-2_wp,-4.00000000000000E-3_wp, 1.00000000000000E-2_wp,&
      & 8.00000000000000E-3_wp,-6.00000000000000E-3_wp,-4.00000000000000E-2_wp,&
      &-2.00000000000000E-2_wp, 5.00000000000000E-3_wp,-1.00000000000000E-2_wp,&
      &-7.00000000000000E-3_wp, 2.00000000000000E-3_wp, 3.00000000000000E-2_wp], &
      & [6, 4])

   call get_structure(mol, "Heavy28", "bih3")
   call test_g(error, ddx_solvation_model%cosmo, mol, qat, dpat=dpat, qpat=qpat)

end subroutine test_g_cosmo_quadrupole_m02

subroutine test_g_cpcm_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.50000000000000E-1_wp,-2.50000000000000E-1_wp, 5.00000000000000E-1_wp,&
      &-5.00000000000000E-1_wp]

   call get_structure(mol, "Heavy28", "bih3")
   call test_g(error, ddx_solvation_model%cpcm, mol, qat)

end subroutine test_g_cpcm_m02

subroutine test_g_pcm_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.50000000000000E-1_wp,-2.50000000000000E-1_wp, 5.00000000000000E-1_wp,&
      &-5.00000000000000E-1_wp]

   call get_structure(mol, "Heavy28", "bih3")
   call test_g(error, ddx_solvation_model%pcm, mol, qat)

end subroutine test_g_pcm_m02

subroutine test_p_cosmo_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.50000000000000E-1_wp,-2.50000000000000E-1_wp, 5.00000000000000E-1_wp,&
      &-5.00000000000000E-1_wp]

   call get_structure(mol, "Heavy28", "bih3")
   call test_p(error, ddx_solvation_model%cosmo, mol, qat)

end subroutine test_p_cosmo_m03

subroutine test_p_cosmo_dipole_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.50000000000000E-1_wp,-2.50000000000000E-1_wp, 5.00000000000000E-1_wp,&
      &-5.00000000000000E-1_wp]
   real(wp), parameter :: dpat(3, 4) = reshape([&
      & 1.00000000000000E-1_wp,-2.00000000000000E-2_wp, 3.00000000000000E-2_wp,&
      &-4.00000000000000E-2_wp, 8.00000000000000E-2_wp,-1.00000000000000E-2_wp,&
      & 6.00000000000000E-2_wp, 2.00000000000000E-2_wp,-7.00000000000000E-2_wp,&
      &-3.00000000000000E-2_wp,-5.00000000000000E-2_wp, 4.00000000000000E-2_wp], &
      & [3, 4])

   call get_structure(mol, "Heavy28", "bih3")
   call test_p(error, ddx_solvation_model%cosmo, mol, qat, dpat=dpat)

end subroutine test_p_cosmo_dipole_m03

subroutine test_p_cosmo_quadrupole_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.50000000000000E-1_wp,-2.50000000000000E-1_wp, 5.00000000000000E-1_wp,&
      &-5.00000000000000E-1_wp]
   real(wp), parameter :: dpat(3, 4) = reshape([&
      & 1.00000000000000E-1_wp,-2.00000000000000E-2_wp, 3.00000000000000E-2_wp,&
      &-4.00000000000000E-2_wp, 8.00000000000000E-2_wp,-1.00000000000000E-2_wp,&
      & 6.00000000000000E-2_wp, 2.00000000000000E-2_wp,-7.00000000000000E-2_wp,&
      &-3.00000000000000E-2_wp,-5.00000000000000E-2_wp, 4.00000000000000E-2_wp], &
      & [3, 4])
   real(wp), parameter :: qpat(6, 4) = reshape([&
      & 2.00000000000000E-2_wp, 1.00000000000000E-2_wp,-3.00000000000000E-2_wp,&
      & 4.00000000000000E-3_wp,-2.00000000000000E-3_wp, 1.00000000000000E-2_wp,&
      &-1.00000000000000E-2_wp, 3.00000000000000E-3_wp, 2.00000000000000E-2_wp,&
      &-5.00000000000000E-3_wp, 7.00000000000000E-3_wp,-2.00000000000000E-2_wp,&
      & 3.00000000000000E-2_wp,-4.00000000000000E-3_wp, 1.00000000000000E-2_wp,&
      & 8.00000000000000E-3_wp,-6.00000000000000E-3_wp,-4.00000000000000E-2_wp,&
      &-2.00000000000000E-2_wp, 5.00000000000000E-3_wp,-1.00000000000000E-2_wp,&
      &-7.00000000000000E-3_wp, 2.00000000000000E-3_wp, 3.00000000000000E-2_wp], &
      & [6, 4])

   call get_structure(mol, "Heavy28", "bih3")
   call test_p(error, ddx_solvation_model%cosmo, mol, qat, dpat=dpat, qpat=qpat)

end subroutine test_p_cosmo_quadrupole_m03

subroutine test_p_cpcm_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.50000000000000E-1_wp,-2.50000000000000E-1_wp, 5.00000000000000E-1_wp,&
      &-5.00000000000000E-1_wp]

   call get_structure(mol, "Heavy28", "bih3")
   call test_p(error, ddx_solvation_model%cpcm, mol, qat)

end subroutine test_p_cpcm_m03

subroutine test_p_pcm_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.50000000000000E-1_wp,-2.50000000000000E-1_wp, 5.00000000000000E-1_wp,&
      &-5.00000000000000E-1_wp]

   call get_structure(mol, "Heavy28", "bih3")
   call test_p(error, ddx_solvation_model%pcm, mol, qat)

end subroutine test_p_pcm_m03


end module test_solvation_ddx
