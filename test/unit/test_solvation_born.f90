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

module test_solvation_born
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
   use tblite_solvation_data_alpb, only: get_alpb_param
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: collect_solvation_born


   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains


!> Collect all exported unit tests
subroutine collect_solvation_born(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("born-1", test_mb01), &
      new_unittest("born-2", test_mb02), &
      new_unittest("born-3", test_mb03), &
      new_unittest("energy-p16", test_e_p16), &
      new_unittest("energy-alpb-gfn1", test_e_alpb_gfn1_all_solvents), &
      new_unittest("energy-alpb-gfn2", test_e_alpb_gfn2_all_solvents), &
      new_unittest("energy-still", test_e_still), &
      new_unittest("energy-gbsa-gfn1", test_e_gbsa_gfn1_all_solvents), &
      new_unittest("energy-gbsa-gfn2", test_e_gbsa_gfn2_all_solvents), &
      new_unittest("energy-charged-p16", test_e_charged_p16), &
      new_unittest("energy-charged-alpb-gfn1", test_e_charged_alpb_gfn1), &
      new_unittest("energy-charged-alpb-gfn2", test_e_charged_alpb_gfn2), &
      new_unittest("energy-charged-still", test_e_charged_still), &
      new_unittest("energy-charged-gbsa-gfn1", test_e_charged_gbsa_gfn1), &
      new_unittest("energy-charged-gbsa-gfn2", test_e_charged_gbsa_gfn2), &
      new_unittest("gradient-p16", test_g_p16), &
      new_unittest("gradient-alpb", test_g_alpb), &
      new_unittest("gradient-still", test_g_still), &
      new_unittest("gradient-gbsa", test_g_gbsa), &
      new_unittest("potential-p16", test_p_p16), &
      new_unittest("potential-alpb", test_p_alpb), &
      new_unittest("potential-still", test_p_still), &
      new_unittest("potential-gbsa", test_p_gbsa), &
      new_unittest("unsupported-solvent", test_unsupported_solvent, should_fail=.true.) &
      ]

end subroutine collect_solvation_born


subroutine test_numg(error, gbobc, mol)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Born radii integrator
   type(born_integrator), intent(inout) :: gbobc
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   integer :: iat, ic
   real(wp), allocatable :: rad(:), sr(:), sl(:), draddr(:, :, :), numg(:, :, :), drad2r(:,:,:,:,:)
   real(wp), parameter :: step = 1.0e-5_wp

   allocate(rad(mol%nat), sr(mol%nat), sl(mol%nat))
   allocate(draddr(3, mol%nat, mol%nat), numg(3, mol%nat, mol%nat))

   call gbobc%get_rad(mol, rad, draddr)

   do iat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call gbobc%get_rad(mol, sr)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call gbobc%get_rad(mol, sl)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step

         numg(ic, iat, :) = 0.5_wp * (sr - sl) / step
      end do
   end do

   if (any(abs(numg - draddr) > thr2)) then
      call test_failed(error, "Born radii derivative does not match finite difference solution")
   end if
end subroutine test_numg

subroutine test_mb01(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   real(wp), allocatable :: rvdw(:), rad(:), draddr(:, :, :)
   real(wp), parameter :: ref(16) = [&
      & 4.07331798531438E+0_wp, 2.56451650429167E+0_wp, 2.97676954448882E+0_wp, &
      & 2.88833112759592E+0_wp, 2.59127476011008E+0_wp, 2.63750279425510E+0_wp, &
      & 3.56149571036025E+0_wp, 2.89090958281373E+0_wp, 4.53815592283277E+0_wp, &
      & 2.46342847720303E+0_wp, 2.72707461251522E+0_wp, 3.52933932564532E+0_wp, &
      & 3.66934919146868E+0_wp, 3.66697827876019E+0_wp, 3.37809284756764E+0_wp, &
      & 4.57932013403544E+0_wp]

   call get_structure(mol, "MB16-43", "01")

   allocate(rad(mol%nat), draddr(3, mol%nat, mol%nat))
   rvdw = get_vdw_rad_d3(mol%num)

   call new_born_integrator(gbobc, mol, rvdw)
   call gbobc%get_rad(mol, rad, draddr)

   if (any(abs(rad - ref) > thr2)) then
      call test_failed(error, "Born radii area values do not match")
      print '(es20.14e1)', rad
      return
   end if

   call test_numg(error, gbobc, mol)

end subroutine test_mb01

subroutine test_mb02(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   real(wp), allocatable :: rvdw(:), rad(:), draddr(:, :, :)
   real(wp), parameter :: ref(16) = [&
      & 3.43501192886252E+0_wp, 5.02404916999371E+0_wp, 4.72253865039692E+0_wp, &
      & 3.52096661104217E+0_wp, 4.76023330956437E+0_wp, 3.46119195863261E+0_wp, &
      & 3.17361370475619E+0_wp, 2.90775065608382E+0_wp, 4.94595287355805E+0_wp, &
      & 3.32592657749444E+0_wp, 4.54348353409109E+0_wp, 4.30924297002105E+0_wp, &
      & 3.47420343563851E+0_wp, 2.82302349370343E+0_wp, 6.67552064739394E+0_wp, &
      & 4.23898159491675E+0_wp]

   call get_structure(mol, "MB16-43", "02")

   allocate(rad(mol%nat), draddr(3, mol%nat, mol%nat))
   rvdw = get_vdw_rad_bondi(mol%num)

   call new_born_integrator(gbobc, mol, rvdw)
   call gbobc%get_rad(mol, rad, draddr)

   if (any(abs(rad - ref) > thr2)) then
      call test_failed(error, "Born radii area values do not match")
      print '(es20.14e1)', rad
      return
   end if

   call test_numg(error, gbobc, mol)

end subroutine test_mb02

subroutine test_mb03(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   real(wp), allocatable :: rvdw(:), rad(:), draddr(:, :, :)
   real(wp), parameter :: ref(16) = [&
      & 4.94764986698701E+0_wp, 3.95122747438791E+0_wp, 4.57075556245289E+0_wp, &
      & 5.46368225070994E+0_wp, 8.24269261139398E+0_wp, 5.68405471762112E+0_wp, &
      & 5.51002325309604E+0_wp, 4.75597020148093E+0_wp, 4.21190195089894E+0_wp, &
      & 4.32836770082885E+0_wp, 4.12684869499911E+0_wp, 5.15171226623248E+0_wp, &
      & 4.83223856996055E+0_wp, 3.02638025720185E+0_wp, 4.05683426506167E+0_wp, &
      & 4.63569783992096E+0_wp]

   call get_structure(mol, "MB16-43", "03")

   allocate(rad(mol%nat), draddr(3, mol%nat, mol%nat))
   rvdw = get_vdw_rad_cosmo(mol%num)

   call new_born_integrator(gbobc, mol, rvdw)
   call gbobc%get_rad(mol, rad, draddr)

   if (any(abs(rad - ref) > thr2)) then
      call test_failed(error, "Born radii area values do not match")
      print '(es20.14e1)', rad
      return
   end if

   call test_numg(error, gbobc, mol)

end subroutine test_mb03


subroutine test_e(error, mol, input, qat, dpat, qpat, ref, method)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solvation model input
   type(alpb_input), intent(in) :: input

   !> Atomic partial charges for this structure
   real(wp), contiguous, intent(in) :: qat(:)

   !> Atomic dipole moments for this structure
   real(wp), contiguous, intent(in) :: dpat(:, :)

   !> Atomic quadrupole moments for this structure
   real(wp), contiguous, intent(in) :: qpat(:, :)

   !> Reference energy
   real(wp), intent(in) :: ref

   !> Method for parameter selection
   character(len=*), intent(in), optional :: method

   type(alpb_solvation) :: solv
   type(alpb_input), allocatable :: scratch_input
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp) :: energy(mol%nat)

   energy = 0.0_wp
   wfn%qat = reshape(qat, [size(qat), 1])
   wfn%dpat = reshape(dpat, [shape(dpat), 1])
   wfn%qpat = reshape(qpat, [shape(qpat), 1])
   allocate(pot%vat(size(qat, 1), 1))
   allocate(pot%vdp(3, size(qat, 1), 1))
   allocate(pot%vqp(6, size(qat, 1), 1))

   scratch_input = input


   if (allocated(input%solvent) .and. present(method)) then 
      call get_alpb_param(scratch_input, mol, method, error)
      if(allocated(error)) then
         call test_failed(error, "No ALPB/GBSA parameters found for the method/solvent")
         return
      end if
   end if

   solv = alpb_solvation(mol, scratch_input, method)

   call solv%update(mol, cache)
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)

   if (abs(sum(energy) - ref) > thr) then
      call test_failed(error, "Energy does not match reference")
      print *, abs(sum(energy) - ref)
   end if
end subroutine test_e


subroutine test_g(error, mol, input, qat, dpat, qpat, method)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Solvation model input
   type(alpb_input), intent(in) :: input

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   real(wp), intent(in) :: dpat(:)
   real(wp), intent(in) :: qpat(:)

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
   wfn%dpat = reshape(dpat, [3, size(qat, 1), 1])
   wfn%qpat = reshape(qpat, [6, size(qat, 1), 1])
   allocate(pot%vat(size(qat, 1), 1))
   allocate(pot%vdp(3, size(qat, 1), 1))
   allocate(pot%vqp(6, size(qat, 1), 1))

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
      if (input%do_multipoles) then
         call test_failed(error, "Multipole gradient does not match")
      else
         call test_failed(error, "Gradient does not match")
      end if
      print '(3es20.13)', gradient
      print '(a)', "---"
      print '(3es20.13)', numg
      print '(a)', "---"
      print '(3es20.13)', gradient - numg
   end if
end subroutine test_g


subroutine test_p(error, mol, input, qat, dpat, qpat, method)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Solvation model input
   type(alpb_input), intent(in) :: input

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   real(wp), intent(in) :: dpat(:)
   real(wp), intent(in) :: qpat(:)

   !> Method for parameter selection
   character(len=*), optional, intent(in) :: method

   type(alpb_solvation) :: solv
   type(alpb_input), allocatable :: scratch_input
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), allocatable :: vat(:)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat)
   integer :: ii

   wfn%qat = reshape(qat, [size(qat), 1])
   wfn%dpat = reshape(dpat, [3, size(qat, 1), 1])
   wfn%qpat = reshape(qpat, [6, size(qat, 1), 1])
   allocate(pot%vat(size(qat, 1), 1))
   allocate(pot%vdp(3, size(qat, 1), 1))
   allocate(pot%vqp(6, size(qat, 1), 1))

   scratch_input = input

   if (allocated(input%solvent) .and. present(method)) then   
      call get_alpb_param(scratch_input, mol, method, error)
      if(allocated(error)) then
         call test_failed(error, "No ALPB/GBSA parameters found for the method/solvent")
         return
      end if
   end if
   
   solv = alpb_solvation(mol, scratch_input, method)

   call solv%update(mol, cache)

   allocate(vat(mol%nat))
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

   energy = 0.0_wp
   pot%vat(:, :) = 0.0_wp
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)

   if (any(abs([pot%vat] - vat) > thr2)) then
      if (input%do_multipoles) then
         call test_failed(error, "Multipole potential does not match")
      else
         call test_failed(error, "Potential does not match")
      end if
      print '(3es20.13)', pot%vat
      print '(a)', "---"
      print '(3es20.13)', vat
      print '(a)', "---"
      print '(3es20.13)', [pot%vat] - vat
   end if
end subroutine test_p


subroutine test_e_p16(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      &-8.99891426876454E-2_wp, 9.42167151049675E-2_wp,-1.49390480695042E-1_wp, &
      &-2.99114077043291E-1_wp, 4.85528544756154E-1_wp,-6.83158233623334E-2_wp, &
      & 1.49989659279676E-2_wp, 2.79361802448324E-1_wp,-1.24070454865924E-1_wp, &
      &-9.36742767297382E-2_wp,-2.19062009156706E-1_wp, 2.14544824574956E-1_wp, &
      & 3.06160853875672E-1_wp,-3.86107183037218E-1_wp,-1.51404848137219E-3_wp, &
      & 3.64257893712230E-2_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      &-3.64731729820407E-2_wp,-2.13686275288388E-2_wp, 6.68806426144540E-2_wp, &
      & 4.04240419990876E-2_wp,-1.49843414114636E-1_wp,-5.58836609439025E-3_wp, &
      &-6.05242232810174E-2_wp, 9.15234565804751E-2_wp,-7.40096212350601E-2_wp, &
      &-2.00610828069480E-2_wp, 8.79914756075611E-2_wp, 1.17474526864424E-1_wp, &
      &-1.23227254106242E-1_wp, 5.17071359855638E-2_wp,-3.45193834725206E-1_wp, &
      & 8.38659891802097E-2_wp,-1.85766281220651E-3_wp, 2.16613946001174E-3_wp, &
      & 7.28212028877821E-2_wp,-1.11656708078832E-1_wp, 4.44987301144173E-2_wp, &
      & 5.21199968878805E-1_wp, 2.42009354648746E-1_wp,-2.59373648270648E-1_wp, &
      & 2.42771833481858E-3_wp, 1.33481324585323E-1_wp, 4.70393507351052E-2_wp, &
      & 6.59447581952794E-2_wp,-4.77190059636476E-2_wp, 1.41459747150361E-1_wp, &
      &-5.74635892974854E-3_wp,-9.43452985224502E-2_wp,-4.14652054192348E-2_wp, &
      & 1.18721461568371E-1_wp,-8.07370651006255E-2_wp, 1.47715903211519E-1_wp, &
      & 1.80051594470901E-1_wp, 2.31074137154501E-1_wp, 4.61275280422502E-1_wp, &
      & 7.00280160506826E-2_wp,-6.56712228336652E-2_wp, 5.83187156889666E-2_wp, &
      &-7.73514384205769E-2_wp,-6.11970211951996E-3_wp,-1.18441889940799E-1_wp, &
      &-2.63615278011550E-2_wp,-8.98662961570448E-2_wp, 9.70285404014370E-2_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      & 2.07660448410257E-2_wp, 4.77562367062922E-2_wp,-6.38136462827382E-2_wp, &
      &-6.08076485231580E-2_wp,-2.89106396726373E-3_wp, 4.30476014417125E-2_wp, &
      & 5.32310785129656E-2_wp, 5.56201199842301E-2_wp,-1.94185137737050E-1_wp, &
      &-1.69435262877121E-2_wp, 4.06941497321145E-2_wp, 1.40954059224083E-1_wp, &
      &-1.77841513693408E-2_wp,-2.75394740661722E-2_wp, 8.23579661010067E-3_wp, &
      & 3.34177535391165E-2_wp,-5.71071859403073E-3_wp, 9.54835475924012E-3_wp, &
      & 2.42936845781245E-2_wp, 5.48198008679115E-3_wp, 9.79238907789204E-4_wp, &
      & 5.02629146832239E-3_wp,-3.26133215380267E-2_wp,-2.52729234859127E-2_wp, &
      & 6.13790503829809E-1_wp, 5.32657270579180E-1_wp, 1.54187752434240E-1_wp, &
      &-5.16490059386104E-1_wp, 5.67205701727938E-2_wp,-7.67978256264049E-1_wp, &
      & 7.83646726355242E-2_wp, 7.93496900696582E-3_wp,-6.64025879132266E-2_wp, &
      &-2.39938346396998E-2_wp, 2.27466901433374E-2_wp,-1.19620847222975E-2_wp, &
      &-2.05612814844119E-2_wp,-4.95270994769136E-2_wp, 3.87946677272715E-2_wp, &
      & 2.09792192994460E-2_wp, 1.57023828153490E-2_wp,-1.82333862428595E-2_wp, &
      & 2.06990008614840E-1_wp,-9.40931225186795E-1_wp,-9.81260233409964E-1_wp, &
      & 1.81683320051775E-1_wp,-2.29434219507657E-1_wp, 7.74270224795126E-1_wp, &
      &-4.49851593643217E-2_wp,-1.43120792142963E-2_wp, 8.76827127392763E-2_wp, &
      & 7.95606845308971E-5_wp, 1.49470897135840E-2_wp,-4.26975533749547E-2_wp, &
      &-3.98095736465717E-2_wp,-4.34657048996160E-2_wp, 1.60460563069101E-2_wp, &
      & 1.10714431698586E-2_wp, 2.00213367980184E-2_wp, 2.37635173396616E-2_wp, &
      & 2.31027463622808E-2_wp,-3.26295375143408E-2_wp,-1.33550726444072E-1_wp, &
      &-1.53667922610297E-1_wp, 3.82022829991620E-2_wp, 1.10447980081791E-1_wp, &
      &-4.50996897167389E-1_wp,-1.02330031032349E+0_wp, 8.40160246681131E-2_wp, &
      &-3.56277305140736E-1_wp,-6.25766606417323E-1_wp, 3.66980872499278E-1_wp, &
      & 2.08713890794846E-1_wp, 8.09936390516330E-1_wp,-1.01293604136016E+0_wp, &
      & 3.35301352962707E-1_wp, 3.38610231741285E-1_wp, 8.04222150565318E-1_wp, &
      & 2.87096625886200E-2_wp, 3.36901296815917E-2_wp,-2.68793522207332E-3_wp, &
      &-2.05780563238218E-2_wp,-6.21642900499840E-3_wp,-2.60217273665459E-2_wp, &
      &-1.50087159351410E-2_wp,-3.69461704433247E-3_wp,-1.82860627486028E-2_wp, &
      & 4.40210284335156E-2_wp, 2.60057075458535E-2_wp, 3.32947786837439E-2_wp, &
      & 4.27111264045819E-1_wp, 4.85129382065812E-2_wp,-6.56980004019461E-1_wp, &
      &-1.30593579650094E-1_wp,-2.07993009790954E-1_wp, 2.29868739973643E-1_wp],&
      & shape(qpat))
   real(wp), parameter :: feps = 80.0_wp

   call get_structure(mol, "MB16-43", "04")

   ! Without multipoles
   call test_e(error, mol, alpb_input(feps, kernel=born_kernel%p16, alpb=.true., do_multipoles=.false.), &
      & qat, dpat, qpat, -7.2619981784854421E-003_wp) ! cosmo radii

   ! With multipoles
   call test_e(error, mol, alpb_input(feps, kernel=born_kernel%p16, alpb=.true., do_multipoles=.true.), &
      & qat, dpat, qpat, -3.5672123597123508E-003_wp) ! cosmo radii

end subroutine test_e_p16

subroutine test_e_alpb_gfn1_all_solvents(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   real(wp), parameter :: qat(*) = [&
      &-8.99890404486076E-2_wp, 9.42168087556583E-2_wp,-1.49387631509499E-1_wp, &
      &-2.99114121895542E-1_wp, 4.85527734875224E-1_wp,-6.83156326406137E-2_wp, &
      & 1.50011293889337E-2_wp, 2.79368544459465E-1_wp,-1.24072452878322E-1_wp, &
      &-9.36760994051244E-2_wp,-2.19062123031622E-1_wp, 2.14538817685587E-1_wp, &
      & 3.06156726072831E-1_wp,-3.86105514712244E-1_wp,-1.51265171389388E-3_wp, &
      & 3.64255069977693E-2_wp]
   real(wp), parameter :: dpat(3, 16) = 0.0_wp
   real(wp), parameter :: qpat(6, 16) = 0.0_wp 

   type(alpb_input) :: input
   integer, parameter :: nsolvents = 25
   character(len=*), parameter :: solvents(*) = [character(len=nsolvents):: &
      & "acetone", "acetonitrile", "aniline", "benzaldehyde", "benzene", &
      & "ch2cl2", "chcl3", "cs2", "dioxane", "dmf", "dmso", "ethanol", &
      & "ether", "ethylacetate", "furane", "hexadecane", "hexane", &
      & "nitromethane", "methanol", "octanol", "phenol", "thf", "toluene", &
      & "water", "woctanol"] 
   real(wp), parameter :: refs(*) = [&
      &-9.52291351801132E-4_wp,-2.66014621365225E-3_wp,-2.81206079506230E-3_wp, &
      &-2.16891247905536E-3_wp,-2.07471669939404E-3_wp,-2.66920593654818E-3_wp, &
      &-4.02007164995396E-3_wp,-1.97291661714106E-3_wp,-3.01602678705293E-3_wp, &
      &-2.12324260213683E-3_wp,-2.36246169433000E-3_wp,-1.50930796674902E-3_wp, &
      &-1.31316585664219E-3_wp,-2.21518979088728E-3_wp,-2.49094454727624E-3_wp, &
      &-1.87322698372293E-3_wp,-1.62359321564825E-3_wp,-2.78662564996861E-3_wp, &
      &-2.45908702328905E-3_wp,-1.48161394314782E-3_wp,-3.81441431379714E-3_wp, &
      &-1.82160716037726E-3_wp,-2.14756160288507E-3_wp,-2.03939585800164E-3_wp, &
      &-1.39250235603086E-3_wp]
   integer :: i

   call get_structure(mol, "MB16-43", "04")

   ! Check GFN1/ALPB for all available solvents
   do i = 1, nsolvents
      solvent = get_solvent_data(solvents(i))
      input = alpb_input(solvent%eps, solvent=solvent%solvent, &
         & kernel=born_kernel%p16, alpb=.true., do_multipoles=.false.)
      call test_e(error, mol, input, qat, dpat, qpat, refs(i), method='gfn1')
      if(allocated(error)) return
   end do 

end subroutine test_e_alpb_gfn1_all_solvents

subroutine test_e_alpb_gfn2_all_solvents(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   real(wp), parameter :: qat(*) = [&
      &-8.99890404486076E-2_wp, 9.42168087556583E-2_wp,-1.49387631509499E-1_wp, &
      &-2.99114121895542E-1_wp, 4.85527734875224E-1_wp,-6.83156326406137E-2_wp, &
      & 1.50011293889337E-2_wp, 2.79368544459465E-1_wp,-1.24072452878322E-1_wp, &
      &-9.36760994051244E-2_wp,-2.19062123031622E-1_wp, 2.14538817685587E-1_wp, &
      & 3.06156726072831E-1_wp,-3.86105514712244E-1_wp,-1.51265171389388E-3_wp, &
      & 3.64255069977693E-2_wp]
   real(wp), parameter :: dpat(3, 16) = 0.0_wp
   real(wp), parameter :: qpat(6, 16) = 0.0_wp
   type(alpb_input) :: input
   integer, parameter :: nsolvents = 25
   character(len=*), parameter :: solvents(*) = [character(len=nsolvents):: &
      & "acetone", "acetonitrile", "aniline", "benzaldehyde", "benzene", &
      & "ch2cl2", "chcl3", "cs2", "dioxane", "dmf", "dmso", "ethanol", &
      & "ether", "ethylacetate", "furane", "hexadecane", "hexane", &
      & "nitromethane", "methanol", "octanol", "phenol", "thf", "toluene", &
      & "water", "woctanol"] 
   real(wp), parameter :: refs(*) = [&
      &-2.57648699865429E-3_wp,-3.28866717415013E-3_wp,-3.88900453989673E-3_wp, &
      &-2.04759081304808E-3_wp,-2.69053080068258E-3_wp,-3.95180642311958E-3_wp, &
      &-5.19631219958033E-3_wp,-3.03037971986796E-3_wp,-4.39796647846814E-3_wp, &
      &-2.60629751634379E-3_wp,-2.96301910862054E-3_wp,-3.47370352219562E-3_wp, &
      &-2.48602622170971E-3_wp,-2.56338439803488E-3_wp,-4.22408900629525E-3_wp, &
      &-3.45855011247743E-3_wp,-4.34096195676291E-3_wp,-2.71441664243754E-3_wp, &
      &-3.63211525626785E-3_wp,-2.83492998452246E-3_wp,-7.29246427793122E-3_wp, &
      &-2.59130816165545E-3_wp,-2.22897169805095E-3_wp,-4.70172426317049E-3_wp, &
      &-2.89954818985142E-3_wp]
   integer :: i

   call get_structure(mol, "MB16-43", "04")

   ! Check GFN2/ALPB for all available solvents
   do i = 1, nsolvents
      solvent = get_solvent_data(solvents(i))
      input = alpb_input(solvent%eps, solvent=solvent%solvent, &
         & kernel=born_kernel%p16, alpb=.true., do_multipoles=.false.)
      call test_e(error, mol, input, qat, dpat, qpat, refs(i), method='gfn2') 
      if(allocated(error)) return
   end do 

end subroutine test_e_alpb_gfn2_all_solvents


subroutine test_e_still(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      & 2.29216160742194E-1_wp, 4.71025340158282E-2_wp,-3.52853843840331E-2_wp, &
      & 3.09109494571158E-2_wp, 3.28899921234171E-1_wp,-2.01752645265381E-1_wp, &
      & 5.46543133415123E-2_wp,-1.09690925298550E-1_wp,-3.47346448335524E-1_wp, &
      & 1.84561311715549E-1_wp,-2.07552642508510E-1_wp, 4.67143385624953E-1_wp, &
      &-1.84214559781700E-2_wp,-1.05015851128102E-1_wp, 6.52406633071011E-2_wp, &
      &-3.82663886540142E-1_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      &-1.36866725567795E-1_wp, 8.98162734808879E-2_wp, 1.93992067964322E-1_wp, &
      &-7.94862211844545E-2_wp, 1.72156341648282E-1_wp, 1.09913940611171E-1_wp, &
      & 4.98989463722096E-3_wp, 5.29170302861464E-2_wp,-2.76614889872771E-2_wp, &
      & 9.88532719472691E-2_wp,-1.40034980230081E-2_wp, 5.34276261339990E-3_wp, &
      &-1.43509540880904E-1_wp,-1.53177661284083E-1_wp,-1.44819803333425E-1_wp, &
      &-2.04156814828549E-2_wp,-1.36595645633772E-1_wp, 9.21268600746920E-2_wp, &
      & 4.82467726581279E-2_wp,-6.83638402837403E-2_wp, 2.25756519899868E-2_wp, &
      &-1.07827737715754E-1_wp,-1.57029392143043E-2_wp,-4.52817748323062E-3_wp, &
      &-5.75565163095502E-2_wp, 8.04433988598797E-2_wp,-6.18156660422769E-2_wp, &
      &-5.40875924342744E-2_wp, 8.84669124477697E-2_wp, 5.79993467185712E-2_wp, &
      & 1.17647171474739E-1_wp, 1.47697439492523E-2_wp, 1.01604279159870E-1_wp, &
      &-5.04764399493363E-1_wp,-9.83938684357523E-3_wp,-3.12571551890408E-1_wp, &
      & 3.55860816794651E-2_wp,-1.31041235123739E-1_wp,-8.75524850448995E-2_wp, &
      & 7.29800329344290E-2_wp,-9.87358811335129E-2_wp,-9.28538506341673E-2_wp, &
      &-8.73682094598224E-2_wp,-2.05679212716268E-1_wp,-1.33550747139953E-1_wp, &
      & 8.30672094549445E-2_wp, 6.79087546790233E-2_wp, 7.61092730005571E-3_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      &-3.21923529804590E-2_wp, 1.01247911778500E-1_wp,-1.60153050884742E-1_wp, &
      & 2.09688971314230E-1_wp,-5.53870278646362E-1_wp, 1.92345403865201E-1_wp, &
      & 1.42745048693680E-1_wp, 3.96650238849656E-1_wp,-1.85556055670220E-2_wp, &
      & 2.09624026315852E-1_wp,-1.00919824206596E-1_wp,-1.24189443126656E-1_wp, &
      &-4.76658863743303E-4_wp, 2.90471575684891E-2_wp, 1.07961071436757E-1_wp, &
      &-9.53915574675821E-2_wp, 3.36258209533697E-2_wp,-1.07484412573014E-1_wp, &
      & 1.27092159973006E-1_wp, 1.21655068472347E-2_wp,-5.40941977072504E-2_wp, &
      &-1.59997403617732E-2_wp,-3.07740637834542E-6_wp,-7.29979622657559E-2_wp, &
      & 5.53303994778954E-1_wp,-2.99239292851686E-1_wp,-2.51102054122008E-1_wp, &
      &-6.26827633889481E-1_wp,-8.42332212885064E-1_wp,-3.02201940656946E-1_wp, &
      & 2.19857432174074E-1_wp, 2.27516201446880E-1_wp,-2.10413952758986E-1_wp, &
      &-4.07388945048188E-1_wp,-1.26641550036282E-1_wp,-9.44347941508889E-3_wp, &
      & 2.93948297678871E-2_wp,-3.43251555808793E-2_wp,-2.11422152270078E-2_wp, &
      & 5.41372742690424E-2_wp,-3.39231131165861E-2_wp,-8.25261454087947E-3_wp, &
      & 2.37070130746082E-1_wp, 3.85988938198882E-3_wp,-1.72140997530014E-1_wp, &
      &-1.84923389243440E-1_wp,-2.83535528712180E-3_wp,-6.49291332160717E-2_wp, &
      &-5.91277580612433E-3_wp, 2.78962136409584E-2_wp,-4.44146350092125E-2_wp, &
      & 4.21395120616379E-2_wp, 1.77901916517195E-2_wp, 5.03274108153378E-2_wp, &
      &-1.69793065817916E-2_wp,-4.35671755332885E-2_wp,-8.26169504201363E-3_wp, &
      &-3.08726223882655E-2_wp, 7.19218791920170E-2_wp, 2.52410016238050E-2_wp, &
      &-1.16489039702432E-1_wp, 1.21709469870565E-1_wp, 1.61694380019930E-1_wp, &
      &-1.78938261587284E-1_wp,-2.40115227627245E-1_wp,-4.52053403174976E-2_wp, &
      & 9.96286547963425E-1_wp,-1.19749514806617E-1_wp,-1.35617319657923E+0_wp, &
      & 8.55455083538778E-1_wp,-9.90077412493467E-1_wp, 3.59886648615801E-1_wp, &
      & 5.22916600071566E-2_wp,-3.91715496559764E-2_wp, 6.83298122566691E-3_wp, &
      &-2.53231481131226E-2_wp,-5.89822395445022E-2_wp,-5.91246412328236E-2_wp, &
      & 4.75850655614168E-3_wp,-1.05070993090875E-2_wp,-9.80745234500818E-3_wp, &
      &-1.41037515552590E-2_wp, 3.45290550320865E-2_wp, 5.04894578886597E-3_wp, &
      & 3.84056205679418E-1_wp, 4.90928417189913E-1_wp,-4.20009072583480E-1_wp, &
      & 1.27378622840086E-1_wp, 2.70322396743841E-1_wp, 3.59528669040627E-2_wp, &
      & 1.02687867732048E-1_wp,-9.41266714173517E-2_wp,-1.57143312769043E-2_wp, &
      &-1.16797766346470E-1_wp,-2.25346145494165E-1_wp,-8.69735364551442E-2_wp],&
      & shape(qpat))

   real(wp), parameter :: feps = 80.0_wp

   call get_structure(mol, "MB16-43", "05")
   
   ! Without multipoles
   call test_e(error, mol, alpb_input(feps, kernel=born_kernel%still, alpb=.false., do_multipoles=.false.), &
      & qat, dpat, qpat, -5.8173811610460532E-003_wp) ! cosmo radii

   ! With multipoles
   call test_e(error, mol, alpb_input(feps, kernel=born_kernel%still, alpb=.false., do_multipoles=.true.), &
      & qat, dpat, qpat, -3.7859501826777529E-003_wp) ! cosmo radii
             
end subroutine test_e_still

subroutine test_e_gbsa_gfn1_all_solvents(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   real(wp), parameter :: qat(*) = [&
      & 2.29208822115185E-1_wp, 4.70816658242009E-2_wp,-3.52834459119718E-2_wp, &
      & 3.09012847396269E-2_wp, 3.28898468854920E-1_wp,-2.01747405019535E-1_wp, &
      & 5.46554362391008E-2_wp,-1.09681283064574E-1_wp,-3.47340505091849E-1_wp, &
      & 1.84567817865267E-1_wp,-2.07552337277185E-1_wp, 4.67140380802351E-1_wp, &
      &-1.84261319200178E-2_wp,-1.05015595324833E-1_wp, 6.52511545312054E-2_wp, &
      &-3.82658324740237E-1_wp]
   real(wp), parameter :: dpat(3, 16) = 0.0_wp
   real(wp), parameter :: qpat(6, 16) = 0.0_wp 
   type(alpb_input) :: input
   integer, parameter :: nsolvents = 12
   character(len=*), parameter :: solvents(*) = [character(len=nsolvents):: &
      & "acetone", "acetonitrile", "benzene", "ch2cl2", "chcl3", "cs2", &
      & "dmso", "ether", "methanol", "thf", "toluene", "water"] 
   real(wp), parameter :: refs(*) = [&
      &-7.27084470362122E-3_wp,-7.24038921758113E-3_wp,-5.71082502695492E-3_wp, &
      &-7.51692928907610E-3_wp,-7.51692928907610E-3_wp,-8.59114683994077E-3_wp, &
      &-7.28115528178506E-3_wp,-5.27553261667375E-3_wp,-6.50347034570797E-3_wp, &
      &-5.38952465735327E-3_wp,-5.71082502695492E-3_wp,-7.46977962503038E-3_wp]
   integer :: i

   call get_structure(mol, "MB16-43", "05")

   ! Check GFN1/GBSA for all available solvents
   do i = 1, nsolvents
      solvent = get_solvent_data(solvents(i))
      input = alpb_input(solvent%eps, solvent=solvent%solvent, &
         & kernel=born_kernel%still, alpb=.false., do_multipoles=.false.)
      call test_e(error, mol, input, qat, dpat, qpat, refs(i), method='gfn1') 
      if(allocated(error)) return
   end do 

end subroutine test_e_gbsa_gfn1_all_solvents

subroutine test_e_gbsa_gfn2_all_solvents(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   real(wp), parameter :: qat(*) = [&
      & 2.29208822115185E-1_wp, 4.70816658242009E-2_wp,-3.52834459119718E-2_wp, &
      & 3.09012847396269E-2_wp, 3.28898468854920E-1_wp,-2.01747405019535E-1_wp, &
      & 5.46554362391008E-2_wp,-1.09681283064574E-1_wp,-3.47340505091849E-1_wp, &
      & 1.84567817865267E-1_wp,-2.07552337277185E-1_wp, 4.67140380802351E-1_wp, &
      &-1.84261319200178E-2_wp,-1.05015595324833E-1_wp, 6.52511545312054E-2_wp, &
      &-3.82658324740237E-1_wp]
   real(wp), parameter :: dpat(3, 16) = 0.0_wp
   real(wp), parameter :: qpat(6, 16) = 0.0_wp 
   type(alpb_input) :: input
   integer, parameter :: nsolvents = 14
   character(len=*), parameter :: solvents(*) = [character(len=nsolvents):: &
      & "acetone", "acetonitrile", "benzene", "ch2cl2", "chcl3", "cs2", &
      & "dmf", "dmso", "ether", "hexane", "methanol", "thf", "toluene", &
      & "water"]
   real(wp), parameter :: refs(*) = [&
      &-2.82214546333666E-3_wp,-3.95653118976006E-3_wp,-3.08316535309036E-3_wp, &
      &-3.93390014572867E-3_wp,-5.40357618800033E-3_wp,-3.69550929789389E-3_wp, &
      &-4.08719342920095E-3_wp,-3.34945147131324E-3_wp,-3.22980115949444E-3_wp, &
      &-3.38991798342756E-3_wp,-3.88172249120498E-3_wp,-3.65698696265639E-3_wp, &
      &-3.08316535309036E-3_wp,-5.08662848390225E-3_wp]
   integer :: i

   call get_structure(mol, "MB16-43", "05")

   ! Check GFN2/GBSA for all available solvents
   do i = 1, nsolvents
      solvent = get_solvent_data(solvents(i))
      input = alpb_input(solvent%eps, solvent=solvent%solvent, &
         & kernel=born_kernel%still, alpb=.false., do_multipoles=.false.)
      call test_e(error, mol, input, qat, dpat, qpat, refs(i), method='gfn2') 
      if(allocated(error)) return
   end do 

end subroutine test_e_gbsa_gfn2_all_solvents


subroutine test_e_charged_p16(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(alpb_input) :: input
   real(wp), parameter :: qat(*) = [&
      & 2.91737167991831E-1_wp,-4.44921281646692E-1_wp, 7.03886898899541E-2_wp, &
      & 4.32458130745259E-2_wp, 3.90393832017139E-2_wp, 9.56766032864156E-2_wp, &
      & 6.01866392558221E-2_wp,-3.31751380042108E-1_wp, 1.56707781479786E-1_wp, &
      & 8.89443731456522E-2_wp,-6.03526873354187E-2_wp, 3.83058968439732E-2_wp, &
      & 9.51289165540499E-2_wp,-9.20860756307965E-2_wp, 7.12122392663975E-2_wp, &
      & 2.76742553014572E-1_wp,-4.38083248713717E-1_wp,-1.76798298145349E-1_wp, &
      & 2.14382387079637E-1_wp, 3.25215078814299E-1_wp,-4.30846325150764E-1_wp, &
      & 1.12908046155543E-1_wp, 6.61648434849400E-2_wp, 9.89171412174742E-2_wp, &
      & 5.80323177221752E-2_wp,-4.43388471029542E-1_wp, 3.38493755177292E-1_wp, &
      &-3.52969753592199E-1_wp, 4.59552645375000E-1_wp,-5.33442813398395E-1_wp, &
      &-5.23460418980519E-1_wp,-3.14158345734236E-1_wp, 8.63777453496141E-2_wp, &
      & 5.51896920278464E-2_wp, 5.60421076086023E-2_wp, 1.02875052304244E-1_wp, &
      & 5.43513415521492E-2_wp,-3.13791816210054E-1_wp, 1.68862983166011E-1_wp, &
      & 8.93759201291419E-2_wp,-8.90121909290432E-2_wp, 4.15150867426933E-2_wp, &
      & 1.22221651251480E-1_wp,-8.26447904844349E-2_wp, 1.00154406589010E-1_wp, &
      & 2.77572579099978E-1_wp,-4.29147369583175E-1_wp,-1.78581481555413E-1_wp, &
      & 2.09487890121871E-1_wp, 3.17109649407645E-1_wp,-4.62948099575570E-1_wp, &
      & 9.81620738022878E-2_wp, 5.14984224707990E-2_wp, 9.63222020737258E-2_wp, &
      & 3.80443799704811E-2_wp,-4.41189291092377E-1_wp, 3.13549746324888E-1_wp, &
      &-4.41335746902051E-1_wp, 3.01219329594079E-1_wp]
   real(wp), parameter :: dpat(3, 59) = 0.0_wp
   real(wp), parameter :: qpat(6, 59) = 0.0_wp 
   real(wp), parameter :: feps = 80.0_wp


   call get_structure(mol, "UPU23", "0a")
   input = alpb_input(feps, kernel=born_kernel%still, alpb=.false., do_multipoles=.false.)
   call test_e(error, mol, input, qat, dpat, qpat, -6.2623428747454107E-2_wp) ! cosmo radii

end subroutine test_e_charged_p16

subroutine test_e_charged_alpb_gfn1(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(alpb_input) :: input
   real(wp), parameter :: qat(*) = [&
      & 2.91737167991831E-1_wp,-4.44921281646692E-1_wp, 7.03886898899541E-2_wp, &
      & 4.32458130745259E-2_wp, 3.90393832017139E-2_wp, 9.56766032864156E-2_wp, &
      & 6.01866392558221E-2_wp,-3.31751380042108E-1_wp, 1.56707781479786E-1_wp, &
      & 8.89443731456522E-2_wp,-6.03526873354187E-2_wp, 3.83058968439732E-2_wp, &
      & 9.51289165540499E-2_wp,-9.20860756307965E-2_wp, 7.12122392663975E-2_wp, &
      & 2.76742553014572E-1_wp,-4.38083248713717E-1_wp,-1.76798298145349E-1_wp, &
      & 2.14382387079637E-1_wp, 3.25215078814299E-1_wp,-4.30846325150764E-1_wp, &
      & 1.12908046155543E-1_wp, 6.61648434849400E-2_wp, 9.89171412174742E-2_wp, &
      & 5.80323177221752E-2_wp,-4.43388471029542E-1_wp, 3.38493755177292E-1_wp, &
      &-3.52969753592199E-1_wp, 4.59552645375000E-1_wp,-5.33442813398395E-1_wp, &
      &-5.23460418980519E-1_wp,-3.14158345734236E-1_wp, 8.63777453496141E-2_wp, &
      & 5.51896920278464E-2_wp, 5.60421076086023E-2_wp, 1.02875052304244E-1_wp, &
      & 5.43513415521492E-2_wp,-3.13791816210054E-1_wp, 1.68862983166011E-1_wp, &
      & 8.93759201291419E-2_wp,-8.90121909290432E-2_wp, 4.15150867426933E-2_wp, &
      & 1.22221651251480E-1_wp,-8.26447904844349E-2_wp, 1.00154406589010E-1_wp, &
      & 2.77572579099978E-1_wp,-4.29147369583175E-1_wp,-1.78581481555413E-1_wp, &
      & 2.09487890121871E-1_wp, 3.17109649407645E-1_wp,-4.62948099575570E-1_wp, &
      & 9.81620738022878E-2_wp, 5.14984224707990E-2_wp, 9.63222020737258E-2_wp, &
      & 3.80443799704811E-2_wp,-4.41189291092377E-1_wp, 3.13549746324888E-1_wp, &
      &-4.41335746902051E-1_wp, 3.01219329594079E-1_wp]
   real(wp), parameter :: dpat(6, 59) = 0.0_wp
   real(wp), parameter :: qpat(6, 59) = 0.0_wp 

   call get_structure(mol, "UPU23", "0a")
   solvent = get_solvent_data("water")
   input = alpb_input(solvent%eps, solvent=solvent%solvent, &
      & kernel=born_kernel%p16, alpb=.true., do_multipoles=.false.)
   call test_e(error, mol, input, qat, dpat, qpat, -9.7339246821001216E-002_wp, method='gfn1')
   
end subroutine test_e_charged_alpb_gfn1

subroutine test_e_charged_alpb_gfn2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(alpb_input) :: input
   real(wp), parameter :: qat(*) = [&
      & 2.91737167991831E-1_wp,-4.44921281646692E-1_wp, 7.03886898899541E-2_wp, &
      & 4.32458130745259E-2_wp, 3.90393832017139E-2_wp, 9.56766032864156E-2_wp, &
      & 6.01866392558221E-2_wp,-3.31751380042108E-1_wp, 1.56707781479786E-1_wp, &
      & 8.89443731456522E-2_wp,-6.03526873354187E-2_wp, 3.83058968439732E-2_wp, &
      & 9.51289165540499E-2_wp,-9.20860756307965E-2_wp, 7.12122392663975E-2_wp, &
      & 2.76742553014572E-1_wp,-4.38083248713717E-1_wp,-1.76798298145349E-1_wp, &
      & 2.14382387079637E-1_wp, 3.25215078814299E-1_wp,-4.30846325150764E-1_wp, &
      & 1.12908046155543E-1_wp, 6.61648434849400E-2_wp, 9.89171412174742E-2_wp, &
      & 5.80323177221752E-2_wp,-4.43388471029542E-1_wp, 3.38493755177292E-1_wp, &
      &-3.52969753592199E-1_wp, 4.59552645375000E-1_wp,-5.33442813398395E-1_wp, &
      &-5.23460418980519E-1_wp,-3.14158345734236E-1_wp, 8.63777453496141E-2_wp, &
      & 5.51896920278464E-2_wp, 5.60421076086023E-2_wp, 1.02875052304244E-1_wp, &
      & 5.43513415521492E-2_wp,-3.13791816210054E-1_wp, 1.68862983166011E-1_wp, &
      & 8.93759201291419E-2_wp,-8.90121909290432E-2_wp, 4.15150867426933E-2_wp, &
      & 1.22221651251480E-1_wp,-8.26447904844349E-2_wp, 1.00154406589010E-1_wp, &
      & 2.77572579099978E-1_wp,-4.29147369583175E-1_wp,-1.78581481555413E-1_wp, &
      & 2.09487890121871E-1_wp, 3.17109649407645E-1_wp,-4.62948099575570E-1_wp, &
      & 9.81620738022878E-2_wp, 5.14984224707990E-2_wp, 9.63222020737258E-2_wp, &
      & 3.80443799704811E-2_wp,-4.41189291092377E-1_wp, 3.13549746324888E-1_wp, &
      &-4.41335746902051E-1_wp, 3.01219329594079E-1_wp]
   real(wp), parameter :: dpat(6, 59) = 0.0_wp
   real(wp), parameter :: qpat(6, 59) = 0.0_wp

   call get_structure(mol, "UPU23", "0a")
   solvent = get_solvent_data("water")
   input = alpb_input(solvent%eps, solvent=solvent%solvent, &
      & kernel=born_kernel%p16, alpb=.true., do_multipoles=.false.)
   call test_e(error, mol, input, qat, dpat, qpat, -0.10736560684364888_wp, method='gfn2')
   
end subroutine test_e_charged_alpb_gfn2

subroutine test_e_charged_still(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(alpb_input) :: input
   real(wp), parameter :: qat(*) = [&
      & 2.91737167991831E-1_wp,-4.44921281646692E-1_wp, 7.03886898899541E-2_wp, &
      & 4.32458130745259E-2_wp, 3.90393832017139E-2_wp, 9.56766032864156E-2_wp, &
      & 6.01866392558221E-2_wp,-3.31751380042108E-1_wp, 1.56707781479786E-1_wp, &
      & 8.89443731456522E-2_wp,-6.03526873354187E-2_wp, 3.83058968439732E-2_wp, &
      & 9.51289165540499E-2_wp,-9.20860756307965E-2_wp, 7.12122392663975E-2_wp, &
      & 2.76742553014572E-1_wp,-4.38083248713717E-1_wp,-1.76798298145349E-1_wp, &
      & 2.14382387079637E-1_wp, 3.25215078814299E-1_wp,-4.30846325150764E-1_wp, &
      & 1.12908046155543E-1_wp, 6.61648434849400E-2_wp, 9.89171412174742E-2_wp, &
      & 5.80323177221752E-2_wp,-4.43388471029542E-1_wp, 3.38493755177292E-1_wp, &
      &-3.52969753592199E-1_wp, 4.59552645375000E-1_wp,-5.33442813398395E-1_wp, &
      &-5.23460418980519E-1_wp,-3.14158345734236E-1_wp, 8.63777453496141E-2_wp, &
      & 5.51896920278464E-2_wp, 5.60421076086023E-2_wp, 1.02875052304244E-1_wp, &
      & 5.43513415521492E-2_wp,-3.13791816210054E-1_wp, 1.68862983166011E-1_wp, &
      & 8.93759201291419E-2_wp,-8.90121909290432E-2_wp, 4.15150867426933E-2_wp, &
      & 1.22221651251480E-1_wp,-8.26447904844349E-2_wp, 1.00154406589010E-1_wp, &
      & 2.77572579099978E-1_wp,-4.29147369583175E-1_wp,-1.78581481555413E-1_wp, &
      & 2.09487890121871E-1_wp, 3.17109649407645E-1_wp,-4.62948099575570E-1_wp, &
      & 9.81620738022878E-2_wp, 5.14984224707990E-2_wp, 9.63222020737258E-2_wp, &
      & 3.80443799704811E-2_wp,-4.41189291092377E-1_wp, 3.13549746324888E-1_wp, &
      &-4.41335746902051E-1_wp, 3.01219329594079E-1_wp]
   real(wp), parameter :: dpat(6, 59) = 0.0_wp
   real(wp), parameter :: qpat(6, 59) = 0.0_wp
   real(wp), parameter :: feps = 80.0_wp

   call get_structure(mol, "UPU23", "0a")
   input = alpb_input(feps, kernel=born_kernel%still, alpb=.false., do_multipoles=.false.)
   call test_e(error, mol, input, qat, dpat, qpat, -6.2623428747454107E-2_wp) ! cosmo radii

end subroutine test_e_charged_still

subroutine test_e_charged_gbsa_gfn1(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(alpb_input) :: input
   real(wp), parameter :: qat(*) = [&
      & 2.91737167991831E-1_wp,-4.44921281646692E-1_wp, 7.03886898899541E-2_wp, &
      & 4.32458130745259E-2_wp, 3.90393832017139E-2_wp, 9.56766032864156E-2_wp, &
      & 6.01866392558221E-2_wp,-3.31751380042108E-1_wp, 1.56707781479786E-1_wp, &
      & 8.89443731456522E-2_wp,-6.03526873354187E-2_wp, 3.83058968439732E-2_wp, &
      & 9.51289165540499E-2_wp,-9.20860756307965E-2_wp, 7.12122392663975E-2_wp, &
      & 2.76742553014572E-1_wp,-4.38083248713717E-1_wp,-1.76798298145349E-1_wp, &
      & 2.14382387079637E-1_wp, 3.25215078814299E-1_wp,-4.30846325150764E-1_wp, &
      & 1.12908046155543E-1_wp, 6.61648434849400E-2_wp, 9.89171412174742E-2_wp, &
      & 5.80323177221752E-2_wp,-4.43388471029542E-1_wp, 3.38493755177292E-1_wp, &
      &-3.52969753592199E-1_wp, 4.59552645375000E-1_wp,-5.33442813398395E-1_wp, &
      &-5.23460418980519E-1_wp,-3.14158345734236E-1_wp, 8.63777453496141E-2_wp, &
      & 5.51896920278464E-2_wp, 5.60421076086023E-2_wp, 1.02875052304244E-1_wp, &
      & 5.43513415521492E-2_wp,-3.13791816210054E-1_wp, 1.68862983166011E-1_wp, &
      & 8.93759201291419E-2_wp,-8.90121909290432E-2_wp, 4.15150867426933E-2_wp, &
      & 1.22221651251480E-1_wp,-8.26447904844349E-2_wp, 1.00154406589010E-1_wp, &
      & 2.77572579099978E-1_wp,-4.29147369583175E-1_wp,-1.78581481555413E-1_wp, &
      & 2.09487890121871E-1_wp, 3.17109649407645E-1_wp,-4.62948099575570E-1_wp, &
      & 9.81620738022878E-2_wp, 5.14984224707990E-2_wp, 9.63222020737258E-2_wp, &
      & 3.80443799704811E-2_wp,-4.41189291092377E-1_wp, 3.13549746324888E-1_wp, &
      &-4.41335746902051E-1_wp, 3.01219329594079E-1_wp]
   real(wp), parameter :: dpat(6, 59) = 0.0_wp
   real(wp), parameter :: qpat(6, 59) = 0.0_wp

   call get_structure(mol, "UPU23", "0a")
   solvent = get_solvent_data("water")
   input = alpb_input(solvent%eps, solvent=solvent%solvent, &
      & kernel=born_kernel%still, alpb=.false., do_multipoles=.false.)
   call test_e(error, mol, input, qat, dpat, qpat, -0.11225040798405941_wp, method='gfn1')
   
end subroutine test_e_charged_gbsa_gfn1

subroutine test_e_charged_gbsa_gfn2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(alpb_input) :: input
   real(wp), parameter :: qat(*) = [&
      & 2.91737167991831E-1_wp,-4.44921281646692E-1_wp, 7.03886898899541E-2_wp, &
      & 4.32458130745259E-2_wp, 3.90393832017139E-2_wp, 9.56766032864156E-2_wp, &
      & 6.01866392558221E-2_wp,-3.31751380042108E-1_wp, 1.56707781479786E-1_wp, &
      & 8.89443731456522E-2_wp,-6.03526873354187E-2_wp, 3.83058968439732E-2_wp, &
      & 9.51289165540499E-2_wp,-9.20860756307965E-2_wp, 7.12122392663975E-2_wp, &
      & 2.76742553014572E-1_wp,-4.38083248713717E-1_wp,-1.76798298145349E-1_wp, &
      & 2.14382387079637E-1_wp, 3.25215078814299E-1_wp,-4.30846325150764E-1_wp, &
      & 1.12908046155543E-1_wp, 6.61648434849400E-2_wp, 9.89171412174742E-2_wp, &
      & 5.80323177221752E-2_wp,-4.43388471029542E-1_wp, 3.38493755177292E-1_wp, &
      &-3.52969753592199E-1_wp, 4.59552645375000E-1_wp,-5.33442813398395E-1_wp, &
      &-5.23460418980519E-1_wp,-3.14158345734236E-1_wp, 8.63777453496141E-2_wp, &
      & 5.51896920278464E-2_wp, 5.60421076086023E-2_wp, 1.02875052304244E-1_wp, &
      & 5.43513415521492E-2_wp,-3.13791816210054E-1_wp, 1.68862983166011E-1_wp, &
      & 8.93759201291419E-2_wp,-8.90121909290432E-2_wp, 4.15150867426933E-2_wp, &
      & 1.22221651251480E-1_wp,-8.26447904844349E-2_wp, 1.00154406589010E-1_wp, &
      & 2.77572579099978E-1_wp,-4.29147369583175E-1_wp,-1.78581481555413E-1_wp, &
      & 2.09487890121871E-1_wp, 3.17109649407645E-1_wp,-4.62948099575570E-1_wp, &
      & 9.81620738022878E-2_wp, 5.14984224707990E-2_wp, 9.63222020737258E-2_wp, &
      & 3.80443799704811E-2_wp,-4.41189291092377E-1_wp, 3.13549746324888E-1_wp, &
      &-4.41335746902051E-1_wp, 3.01219329594079E-1_wp]
   real(wp), parameter :: dpat(6, 59) = 0.0_wp
   real(wp), parameter :: qpat(6, 59) = 0.0_wp

   call get_structure(mol, "UPU23", "0a")
   solvent = get_solvent_data("water")
   input = alpb_input(solvent%eps, solvent=solvent%solvent, &
      & kernel=born_kernel%still, alpb=.false., do_multipoles=.false.)
   call test_e(error, mol, input, qat, dpat, qpat, -9.5967790364628852E-002_wp, method='gfn2')
   
end subroutine test_e_charged_gbsa_gfn2


subroutine test_g_p16(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(alpb_input) :: input
   real(wp), parameter :: qat(*) = [&
      &-3.28160099939119E-1_wp, 3.63839789415764E-1_wp,-9.39678438468329E-1_wp,&
      &-5.67353131753718E-1_wp, 3.91549236321241E-1_wp,-7.25527696006913E-1_wp,&
      & 2.81658498913997E-1_wp, 6.37609035711367E-1_wp,-3.32341795239006E-1_wp,&
      & 1.57147288398166E-1_wp,-2.15641708745235E-1_wp, 6.23132575019337E-1_wp,&
      & 6.86303747651841E-1_wp,-5.41860284045708E-1_wp, 2.80264616247842E-1_wp,&
      & 2.29058416549171E-1_wp]

   real(wp), parameter :: dpat(*) = [&
      &-1.02866383946914E-1_wp,-6.59701808627942E-2_wp, 1.67698043002308E-1_wp,&
      &-6.23739039114127E-2_wp,-1.97051591839775E-1_wp,-1.22265382825162E-1_wp,&
      &-1.64257857747057E-2_wp, 2.27490812781041E-2_wp,-1.84302388618745E-2_wp,&
      &-8.56553766483229E-4_wp, 2.74842733013248E-3_wp, 2.55949906672251E-3_wp,&
      &-1.61115426160772E-1_wp, 1.66196888766802E-1_wp,-1.20546783549511E-1_wp,&
      & 9.88324537231105E-2_wp, 2.50428455703373E-3_wp,-2.88635579792292E-3_wp,&
      & 8.34445890794251E-2_wp,-1.34421662914771E-1_wp, 3.99666225357878E-2_wp,&
      &-4.83529798793855E-2_wp, 5.74113854317782E-3_wp, 5.87539955830004E-2_wp,&
      &-2.97630886764780E-3_wp, 2.19742923198916E-1_wp, 6.43178710881168E-2_wp,&
      & 9.28242957167868E-2_wp,-6.85043369046399E-2_wp, 2.00259137721060E-1_wp,&
      & 7.09776971104250E-5_wp, 2.03464908299606E-2_wp, 1.45200228907336E-2_wp,&
      &-2.15470436895541E-2_wp,-1.67828530829147E-2_wp,-3.81704914987720E-2_wp,&
      &-6.92951149767597E-2_wp, 9.85295808358629E-2_wp,-4.67219521543993E-3_wp,&
      &-1.41470268241085E-2_wp, 1.89740357101272E-2_wp,-6.41196409236768E-3_wp,&
      &-8.82996388817125E-2_wp,-1.65518095011688E-2_wp,-1.35247354235483E-1_wp,&
      &-2.61342733605818E-1_wp,-6.17705247109865E-2_wp, 4.23216561956390E-1_wp]

    real(wp), parameter :: qpat(*) = [&
    & 4.51357475089842E-2_wp,-4.17318742479749E-2_wp,-1.36136078209980E-2_wp,&
      & 9.12863796761614E-2_wp, 7.51195530423392E-2_wp,-3.15221396879885E-2_wp,&
      &-9.80963767772843E-2_wp, 1.45280277246739E-1_wp,-5.11292962085315E-1_wp,&
      &-1.16185111757740E-2_wp,-2.61668573633159E-1_wp, 6.09389338862585E-1_wp,&
      & 6.52037683652264E-3_wp, 3.69171710073827E-2_wp,-1.39436488066209E-2_wp,&
      &-2.83004471610422E-2_wp, 3.62439775801960E-2_wp, 7.42327197009803E-3_wp,&
      &-3.86296110845069E-2_wp,-1.21112119157970E-2_wp, 9.97808370410505E-3_wp,&
      &-1.91278825077846E-2_wp, 6.47092359099180E-2_wp, 2.86515273804006E-2_wp,&
      & 8.56200033892328E-1_wp, 9.74235458476180E-1_wp, 1.18086922059802E-1_wp,&
      &-8.71949362470118E-1_wp, 2.97351048258474E-1_wp,-9.74286955952133E-1_wp,&
      &-2.09879239083559E-1_wp, 9.31236036478151E-2_wp, 5.30742905878537E-2_wp,&
      & 2.06227339580677E-2_wp,-7.56884351806246E-3_wp, 1.56804948495703E-1_wp,&
      & 2.93163202989695E-2_wp, 1.59225627857236E-1_wp,-1.42600503718415E-1_wp,&
      &-3.26169646892855E-2_wp, 6.57594862333770E-2_wp, 1.13284183419448E-1_wp,&
      &-3.36105260187027E-1_wp,-6.53815224929778E-1_wp, 5.72904099232633E-2_wp,&
      & 3.91795155531810E-1_wp, 1.49074708451067E-1_wp, 2.78814850263761E-1_wp,&
      & 1.39849278576512E-1_wp, 1.62052454092380E-2_wp,-2.52552662853771E-1_wp,&
      & 4.12583681109575E-3_wp,-1.03140151821042E-1_wp, 1.12703384277264E-1_wp,&
      & 1.14113595393239E-1_wp, 7.38176190606366E-2_wp, 1.63093603236892E-1_wp,&
      &-2.24757223472586E-1_wp, 1.63465059762335E-1_wp,-2.77207198630140E-1_wp,&
      & 4.49629181418368E-1_wp, 4.25710033896662E-1_wp,-4.29670969707868E-1_wp,&
      &-8.92661030885367E-1_wp,-1.88058091049599E-1_wp,-1.99582117104920E-2_wp,&
      & 7.59715738599860E-1_wp, 8.23127259363381E-1_wp,-1.17603646025185E+0_wp,&
      &-1.44601626881170E+0_wp, 6.22081607059793E-1_wp, 4.16320721651909E-1_wp,&
      &-3.31955764455255E-1_wp, 1.15087423515139E+0_wp,-2.25645212472335E-1_wp,&
      & 1.18140494652646E+0_wp,-7.02380712786817E-1_wp, 5.57600976927600E-1_wp,&
      & 5.67085367397220E-2_wp, 1.13737943078387E-1_wp,-1.02296461455588E-2_wp,&
      &-8.06936077116478E-2_wp, 5.38804500826868E-2_wp,-4.64788905941666E-2_wp,&
      & 2.11496816821796E-2_wp,-3.21916493532053E-2_wp, 1.22390310502235E-1_wp,&
      &-1.62481131000237E-1_wp,-4.62447548254617E-2_wp,-1.43539992184418E-1_wp,&
      & 4.32313899320035E-1_wp,-2.95562991660850E-1_wp,-1.10104884940963E+0_wp,&
      & 3.42643149336475E-1_wp,-6.71823237278679E-1_wp, 6.68734950089545E-1_wp]

   real(wp), parameter :: feps = 80.0_wp

   call get_structure(mol, "MB16-43", "06")
   ! Without multipoles
   input = alpb_input(feps, kernel=born_kernel%p16, alpb = .true., do_multipoles=.false.)
   call test_g(error, mol, input, qat, dpat, qpat, method='gfn2')

   ! With multipoles
   input = alpb_input(feps, kernel=born_kernel%p16, alpb = .true., do_multipoles=.true.)
   call test_g(error, mol, input, qat, dpat, qpat, method='gfn2')

end subroutine test_g_p16

subroutine test_g_alpb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(alpb_input) :: input
   real(wp), parameter :: qat(*) = [&
      &-3.28160099939119E-1_wp, 3.63839789415764E-1_wp,-9.39678438468329E-1_wp,&
      &-5.67353131753718E-1_wp, 3.91549236321241E-1_wp,-7.25527696006913E-1_wp,&
      & 2.81658498913997E-1_wp, 6.37609035711367E-1_wp,-3.32341795239006E-1_wp,&
      & 1.57147288398166E-1_wp,-2.15641708745235E-1_wp, 6.23132575019337E-1_wp,&
      & 6.86303747651841E-1_wp,-5.41860284045708E-1_wp, 2.80264616247842E-1_wp,&
      & 2.29058416549171E-1_wp]

   real(wp), parameter :: dpat(*) = [&
      &-1.02866383946914E-1_wp,-6.59701808627942E-2_wp, 1.67698043002308E-1_wp,&
      &-6.23739039114127E-2_wp,-1.97051591839775E-1_wp,-1.22265382825162E-1_wp,&
      &-1.64257857747057E-2_wp, 2.27490812781041E-2_wp,-1.84302388618745E-2_wp,&
      &-8.56553766483229E-4_wp, 2.74842733013248E-3_wp, 2.55949906672251E-3_wp,&
      &-1.61115426160772E-1_wp, 1.66196888766802E-1_wp,-1.20546783549511E-1_wp,&
      & 9.88324537231105E-2_wp, 2.50428455703373E-3_wp,-2.88635579792292E-3_wp,&
      & 8.34445890794251E-2_wp,-1.34421662914771E-1_wp, 3.99666225357878E-2_wp,&
      &-4.83529798793855E-2_wp, 5.74113854317782E-3_wp, 5.87539955830004E-2_wp,&
      &-2.97630886764780E-3_wp, 2.19742923198916E-1_wp, 6.43178710881168E-2_wp,&
      & 9.28242957167868E-2_wp,-6.85043369046399E-2_wp, 2.00259137721060E-1_wp,&
      & 7.09776971104250E-5_wp, 2.03464908299606E-2_wp, 1.45200228907336E-2_wp,&
      &-2.15470436895541E-2_wp,-1.67828530829147E-2_wp,-3.81704914987720E-2_wp,&
      &-6.92951149767597E-2_wp, 9.85295808358629E-2_wp,-4.67219521543993E-3_wp,&
      &-1.41470268241085E-2_wp, 1.89740357101272E-2_wp,-6.41196409236768E-3_wp,&
      &-8.82996388817125E-2_wp,-1.65518095011688E-2_wp,-1.35247354235483E-1_wp,&
      &-2.61342733605818E-1_wp,-6.17705247109865E-2_wp, 4.23216561956390E-1_wp]

    real(wp), parameter :: qpat(*) = [&
    & 4.51357475089842E-2_wp,-4.17318742479749E-2_wp,-1.36136078209980E-2_wp,&
      & 9.12863796761614E-2_wp, 7.51195530423392E-2_wp,-3.15221396879885E-2_wp,&
      &-9.80963767772843E-2_wp, 1.45280277246739E-1_wp,-5.11292962085315E-1_wp,&
      &-1.16185111757740E-2_wp,-2.61668573633159E-1_wp, 6.09389338862585E-1_wp,&
      & 6.52037683652264E-3_wp, 3.69171710073827E-2_wp,-1.39436488066209E-2_wp,&
      &-2.83004471610422E-2_wp, 3.62439775801960E-2_wp, 7.42327197009803E-3_wp,&
      &-3.86296110845069E-2_wp,-1.21112119157970E-2_wp, 9.97808370410505E-3_wp,&
      &-1.91278825077846E-2_wp, 6.47092359099180E-2_wp, 2.86515273804006E-2_wp,&
      & 8.56200033892328E-1_wp, 9.74235458476180E-1_wp, 1.18086922059802E-1_wp,&
      &-8.71949362470118E-1_wp, 2.97351048258474E-1_wp,-9.74286955952133E-1_wp,&
      &-2.09879239083559E-1_wp, 9.31236036478151E-2_wp, 5.30742905878537E-2_wp,&
      & 2.06227339580677E-2_wp,-7.56884351806246E-3_wp, 1.56804948495703E-1_wp,&
      & 2.93163202989695E-2_wp, 1.59225627857236E-1_wp,-1.42600503718415E-1_wp,&
      &-3.26169646892855E-2_wp, 6.57594862333770E-2_wp, 1.13284183419448E-1_wp,&
      &-3.36105260187027E-1_wp,-6.53815224929778E-1_wp, 5.72904099232633E-2_wp,&
      & 3.91795155531810E-1_wp, 1.49074708451067E-1_wp, 2.78814850263761E-1_wp,&
      & 1.39849278576512E-1_wp, 1.62052454092380E-2_wp,-2.52552662853771E-1_wp,&
      & 4.12583681109575E-3_wp,-1.03140151821042E-1_wp, 1.12703384277264E-1_wp,&
      & 1.14113595393239E-1_wp, 7.38176190606366E-2_wp, 1.63093603236892E-1_wp,&
      &-2.24757223472586E-1_wp, 1.63465059762335E-1_wp,-2.77207198630140E-1_wp,&
      & 4.49629181418368E-1_wp, 4.25710033896662E-1_wp,-4.29670969707868E-1_wp,&
      &-8.92661030885367E-1_wp,-1.88058091049599E-1_wp,-1.99582117104920E-2_wp,&
      & 7.59715738599860E-1_wp, 8.23127259363381E-1_wp,-1.17603646025185E+0_wp,&
      &-1.44601626881170E+0_wp, 6.22081607059793E-1_wp, 4.16320721651909E-1_wp,&
      &-3.31955764455255E-1_wp, 1.15087423515139E+0_wp,-2.25645212472335E-1_wp,&
      & 1.18140494652646E+0_wp,-7.02380712786817E-1_wp, 5.57600976927600E-1_wp,&
      & 5.67085367397220E-2_wp, 1.13737943078387E-1_wp,-1.02296461455588E-2_wp,&
      &-8.06936077116478E-2_wp, 5.38804500826868E-2_wp,-4.64788905941666E-2_wp,&
      & 2.11496816821796E-2_wp,-3.21916493532053E-2_wp, 1.22390310502235E-1_wp,&
      &-1.62481131000237E-1_wp,-4.62447548254617E-2_wp,-1.43539992184418E-1_wp,&
      & 4.32313899320035E-1_wp,-2.95562991660850E-1_wp,-1.10104884940963E+0_wp,&
      & 3.42643149336475E-1_wp,-6.71823237278679E-1_wp, 6.68734950089545E-1_wp]

   call get_structure(mol, "MB16-43", "06")
   solvent = get_solvent_data("water")

   ! Without multipoles
   input = alpb_input(solvent%eps, solvent=solvent%solvent, &
      & kernel=born_kernel%p16, alpb=.true., do_multipoles=.false.)
   call test_g(error, mol, input, qat, dpat, qpat)

   ! With multipoles
   input = alpb_input(solvent%eps, solvent=solvent%solvent, &
      & kernel=born_kernel%p16, alpb=.true., do_multipoles=.true.)
   call test_g(error, mol, input, qat, dpat, qpat)

end subroutine test_g_alpb


subroutine test_g_still(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(alpb_input) :: input
   real(wp), parameter :: qat(*) = [&
      &-3.28160099939119E-1_wp, 3.63839789415764E-1_wp,-9.39678438468329E-1_wp,&
      &-5.67353131753718E-1_wp, 3.91549236321241E-1_wp,-7.25527696006913E-1_wp,&
      & 2.81658498913997E-1_wp, 6.37609035711367E-1_wp,-3.32341795239006E-1_wp,&
      & 1.57147288398166E-1_wp,-2.15641708745235E-1_wp, 6.23132575019337E-1_wp,&
      & 6.86303747651841E-1_wp,-5.41860284045708E-1_wp, 2.80264616247842E-1_wp,&
      & 2.29058416549171E-1_wp]

   real(wp), parameter :: dpat(*) = [&
      &-1.02866383946914E-1_wp,-6.59701808627942E-2_wp, 1.67698043002308E-1_wp,&
      &-6.23739039114127E-2_wp,-1.97051591839775E-1_wp,-1.22265382825162E-1_wp,&
      &-1.64257857747057E-2_wp, 2.27490812781041E-2_wp,-1.84302388618745E-2_wp,&
      &-8.56553766483229E-4_wp, 2.74842733013248E-3_wp, 2.55949906672251E-3_wp,&
      &-1.61115426160772E-1_wp, 1.66196888766802E-1_wp,-1.20546783549511E-1_wp,&
      & 9.88324537231105E-2_wp, 2.50428455703373E-3_wp,-2.88635579792292E-3_wp,&
      & 8.34445890794251E-2_wp,-1.34421662914771E-1_wp, 3.99666225357878E-2_wp,&
      &-4.83529798793855E-2_wp, 5.74113854317782E-3_wp, 5.87539955830004E-2_wp,&
      &-2.97630886764780E-3_wp, 2.19742923198916E-1_wp, 6.43178710881168E-2_wp,&
      & 9.28242957167868E-2_wp,-6.85043369046399E-2_wp, 2.00259137721060E-1_wp,&
      & 7.09776971104250E-5_wp, 2.03464908299606E-2_wp, 1.45200228907336E-2_wp,&
      &-2.15470436895541E-2_wp,-1.67828530829147E-2_wp,-3.81704914987720E-2_wp,&
      &-6.92951149767597E-2_wp, 9.85295808358629E-2_wp,-4.67219521543993E-3_wp,&
      &-1.41470268241085E-2_wp, 1.89740357101272E-2_wp,-6.41196409236768E-3_wp,&
      &-8.82996388817125E-2_wp,-1.65518095011688E-2_wp,-1.35247354235483E-1_wp,&
      &-2.61342733605818E-1_wp,-6.17705247109865E-2_wp, 4.23216561956390E-1_wp]

    real(wp), parameter :: qpat(*) = [&
    & 4.51357475089842E-2_wp,-4.17318742479749E-2_wp,-1.36136078209980E-2_wp,&
      & 9.12863796761614E-2_wp, 7.51195530423392E-2_wp,-3.15221396879885E-2_wp,&
      &-9.80963767772843E-2_wp, 1.45280277246739E-1_wp,-5.11292962085315E-1_wp,&
      &-1.16185111757740E-2_wp,-2.61668573633159E-1_wp, 6.09389338862585E-1_wp,&
      & 6.52037683652264E-3_wp, 3.69171710073827E-2_wp,-1.39436488066209E-2_wp,&
      &-2.83004471610422E-2_wp, 3.62439775801960E-2_wp, 7.42327197009803E-3_wp,&
      &-3.86296110845069E-2_wp,-1.21112119157970E-2_wp, 9.97808370410505E-3_wp,&
      &-1.91278825077846E-2_wp, 6.47092359099180E-2_wp, 2.86515273804006E-2_wp,&
      & 8.56200033892328E-1_wp, 9.74235458476180E-1_wp, 1.18086922059802E-1_wp,&
      &-8.71949362470118E-1_wp, 2.97351048258474E-1_wp,-9.74286955952133E-1_wp,&
      &-2.09879239083559E-1_wp, 9.31236036478151E-2_wp, 5.30742905878537E-2_wp,&
      & 2.06227339580677E-2_wp,-7.56884351806246E-3_wp, 1.56804948495703E-1_wp,&
      & 2.93163202989695E-2_wp, 1.59225627857236E-1_wp,-1.42600503718415E-1_wp,&
      &-3.26169646892855E-2_wp, 6.57594862333770E-2_wp, 1.13284183419448E-1_wp,&
      &-3.36105260187027E-1_wp,-6.53815224929778E-1_wp, 5.72904099232633E-2_wp,&
      & 3.91795155531810E-1_wp, 1.49074708451067E-1_wp, 2.78814850263761E-1_wp,&
      & 1.39849278576512E-1_wp, 1.62052454092380E-2_wp,-2.52552662853771E-1_wp,&
      & 4.12583681109575E-3_wp,-1.03140151821042E-1_wp, 1.12703384277264E-1_wp,&
      & 1.14113595393239E-1_wp, 7.38176190606366E-2_wp, 1.63093603236892E-1_wp,&
      &-2.24757223472586E-1_wp, 1.63465059762335E-1_wp,-2.77207198630140E-1_wp,&
      & 4.49629181418368E-1_wp, 4.25710033896662E-1_wp,-4.29670969707868E-1_wp,&
      &-8.92661030885367E-1_wp,-1.88058091049599E-1_wp,-1.99582117104920E-2_wp,&
      & 7.59715738599860E-1_wp, 8.23127259363381E-1_wp,-1.17603646025185E+0_wp,&
      &-1.44601626881170E+0_wp, 6.22081607059793E-1_wp, 4.16320721651909E-1_wp,&
      &-3.31955764455255E-1_wp, 1.15087423515139E+0_wp,-2.25645212472335E-1_wp,&
      & 1.18140494652646E+0_wp,-7.02380712786817E-1_wp, 5.57600976927600E-1_wp,&
      & 5.67085367397220E-2_wp, 1.13737943078387E-1_wp,-1.02296461455588E-2_wp,&
      &-8.06936077116478E-2_wp, 5.38804500826868E-2_wp,-4.64788905941666E-2_wp,&
      & 2.11496816821796E-2_wp,-3.21916493532053E-2_wp, 1.22390310502235E-1_wp,&
      &-1.62481131000237E-1_wp,-4.62447548254617E-2_wp,-1.43539992184418E-1_wp,&
      & 4.32313899320035E-1_wp,-2.95562991660850E-1_wp,-1.10104884940963E+0_wp,&
      & 3.42643149336475E-1_wp,-6.71823237278679E-1_wp, 6.68734950089545E-1_wp]

   real(wp), parameter :: feps = 80.0_wp

   call get_structure(mol, "MB16-43", "07")

   ! Without multipoles
   input = alpb_input(feps, kernel=born_kernel%still, alpb=.true., do_multipoles=.false.)
   call test_g(error, mol, input, qat, dpat, qpat)

   ! With multipoles
   input = alpb_input(feps, kernel=born_kernel%still, alpb=.true., do_multipoles=.true.)
   call test_g(error, mol, input, qat, dpat, qpat)

end subroutine test_g_still

subroutine test_g_gbsa(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(alpb_input) :: input
   real(wp), parameter :: qat(*) = [&
      &-3.28160099939119E-1_wp, 3.63839789415764E-1_wp,-9.39678438468329E-1_wp,&
      &-5.67353131753718E-1_wp, 3.91549236321241E-1_wp,-7.25527696006913E-1_wp,&
      & 2.81658498913997E-1_wp, 6.37609035711367E-1_wp,-3.32341795239006E-1_wp,&
      & 1.57147288398166E-1_wp,-2.15641708745235E-1_wp, 6.23132575019337E-1_wp,&
      & 6.86303747651841E-1_wp,-5.41860284045708E-1_wp, 2.80264616247842E-1_wp,&
      & 2.29058416549171E-1_wp]

   real(wp), parameter :: dpat(*) = [&
      &-1.02866383946914E-1_wp,-6.59701808627942E-2_wp, 1.67698043002308E-1_wp,&
      &-6.23739039114127E-2_wp,-1.97051591839775E-1_wp,-1.22265382825162E-1_wp,&
      &-1.64257857747057E-2_wp, 2.27490812781041E-2_wp,-1.84302388618745E-2_wp,&
      &-8.56553766483229E-4_wp, 2.74842733013248E-3_wp, 2.55949906672251E-3_wp,&
      &-1.61115426160772E-1_wp, 1.66196888766802E-1_wp,-1.20546783549511E-1_wp,&
      & 9.88324537231105E-2_wp, 2.50428455703373E-3_wp,-2.88635579792292E-3_wp,&
      & 8.34445890794251E-2_wp,-1.34421662914771E-1_wp, 3.99666225357878E-2_wp,&
      &-4.83529798793855E-2_wp, 5.74113854317782E-3_wp, 5.87539955830004E-2_wp,&
      &-2.97630886764780E-3_wp, 2.19742923198916E-1_wp, 6.43178710881168E-2_wp,&
      & 9.28242957167868E-2_wp,-6.85043369046399E-2_wp, 2.00259137721060E-1_wp,&
      & 7.09776971104250E-5_wp, 2.03464908299606E-2_wp, 1.45200228907336E-2_wp,&
      &-2.15470436895541E-2_wp,-1.67828530829147E-2_wp,-3.81704914987720E-2_wp,&
      &-6.92951149767597E-2_wp, 9.85295808358629E-2_wp,-4.67219521543993E-3_wp,&
      &-1.41470268241085E-2_wp, 1.89740357101272E-2_wp,-6.41196409236768E-3_wp,&
      &-8.82996388817125E-2_wp,-1.65518095011688E-2_wp,-1.35247354235483E-1_wp,&
      &-2.61342733605818E-1_wp,-6.17705247109865E-2_wp, 4.23216561956390E-1_wp]

    real(wp), parameter :: qpat(*) = [&
    & 4.51357475089842E-2_wp,-4.17318742479749E-2_wp,-1.36136078209980E-2_wp,&
      & 9.12863796761614E-2_wp, 7.51195530423392E-2_wp,-3.15221396879885E-2_wp,&
      &-9.80963767772843E-2_wp, 1.45280277246739E-1_wp,-5.11292962085315E-1_wp,&
      &-1.16185111757740E-2_wp,-2.61668573633159E-1_wp, 6.09389338862585E-1_wp,&
      & 6.52037683652264E-3_wp, 3.69171710073827E-2_wp,-1.39436488066209E-2_wp,&
      &-2.83004471610422E-2_wp, 3.62439775801960E-2_wp, 7.42327197009803E-3_wp,&
      &-3.86296110845069E-2_wp,-1.21112119157970E-2_wp, 9.97808370410505E-3_wp,&
      &-1.91278825077846E-2_wp, 6.47092359099180E-2_wp, 2.86515273804006E-2_wp,&
      & 8.56200033892328E-1_wp, 9.74235458476180E-1_wp, 1.18086922059802E-1_wp,&
      &-8.71949362470118E-1_wp, 2.97351048258474E-1_wp,-9.74286955952133E-1_wp,&
      &-2.09879239083559E-1_wp, 9.31236036478151E-2_wp, 5.30742905878537E-2_wp,&
      & 2.06227339580677E-2_wp,-7.56884351806246E-3_wp, 1.56804948495703E-1_wp,&
      & 2.93163202989695E-2_wp, 1.59225627857236E-1_wp,-1.42600503718415E-1_wp,&
      &-3.26169646892855E-2_wp, 6.57594862333770E-2_wp, 1.13284183419448E-1_wp,&
      &-3.36105260187027E-1_wp,-6.53815224929778E-1_wp, 5.72904099232633E-2_wp,&
      & 3.91795155531810E-1_wp, 1.49074708451067E-1_wp, 2.78814850263761E-1_wp,&
      & 1.39849278576512E-1_wp, 1.62052454092380E-2_wp,-2.52552662853771E-1_wp,&
      & 4.12583681109575E-3_wp,-1.03140151821042E-1_wp, 1.12703384277264E-1_wp,&
      & 1.14113595393239E-1_wp, 7.38176190606366E-2_wp, 1.63093603236892E-1_wp,&
      &-2.24757223472586E-1_wp, 1.63465059762335E-1_wp,-2.77207198630140E-1_wp,&
      & 4.49629181418368E-1_wp, 4.25710033896662E-1_wp,-4.29670969707868E-1_wp,&
      &-8.92661030885367E-1_wp,-1.88058091049599E-1_wp,-1.99582117104920E-2_wp,&
      & 7.59715738599860E-1_wp, 8.23127259363381E-1_wp,-1.17603646025185E+0_wp,&
      &-1.44601626881170E+0_wp, 6.22081607059793E-1_wp, 4.16320721651909E-1_wp,&
      &-3.31955764455255E-1_wp, 1.15087423515139E+0_wp,-2.25645212472335E-1_wp,&
      & 1.18140494652646E+0_wp,-7.02380712786817E-1_wp, 5.57600976927600E-1_wp,&
      & 5.67085367397220E-2_wp, 1.13737943078387E-1_wp,-1.02296461455588E-2_wp,&
      &-8.06936077116478E-2_wp, 5.38804500826868E-2_wp,-4.64788905941666E-2_wp,&
      & 2.11496816821796E-2_wp,-3.21916493532053E-2_wp, 1.22390310502235E-1_wp,&
      &-1.62481131000237E-1_wp,-4.62447548254617E-2_wp,-1.43539992184418E-1_wp,&
      & 4.32313899320035E-1_wp,-2.95562991660850E-1_wp,-1.10104884940963E+0_wp,&
      & 3.42643149336475E-1_wp,-6.71823237278679E-1_wp, 6.68734950089545E-1_wp]

      

   call get_structure(mol, "MB16-43", "07")
   solvent = get_solvent_data("water")

   ! Without multipoles
   input = alpb_input(solvent%eps, solvent=solvent%solvent, &
      & kernel=born_kernel%still, alpb=.false., do_multipoles=.false.)
   call test_g(error, mol, input, qat, dpat, qpat, method='gfn2')

   ! With multipoles
   input = alpb_input(solvent%eps, solvent=solvent%solvent, &
      & kernel=born_kernel%still, alpb=.false., do_multipoles=.true.)
   call test_g(error, mol, input, qat, dpat, qpat, method='gfn2')

end subroutine test_g_gbsa


subroutine test_p_p16(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(alpb_input) :: input
   real(wp), parameter :: qat(*) = [&
      &-3.28160099939119E-1_wp, 3.63839789415764E-1_wp,-9.39678438468329E-1_wp,&
      &-5.67353131753718E-1_wp, 3.91549236321241E-1_wp,-7.25527696006913E-1_wp,&
      & 2.81658498913997E-1_wp, 6.37609035711367E-1_wp,-3.32341795239006E-1_wp,&
      & 1.57147288398166E-1_wp,-2.15641708745235E-1_wp, 6.23132575019337E-1_wp,&
      & 6.86303747651841E-1_wp,-5.41860284045708E-1_wp, 2.80264616247842E-1_wp,&
      & 2.29058416549171E-1_wp]

   real(wp), parameter :: dpat(*) = [&
      &-1.02866383946914E-1_wp,-6.59701808627942E-2_wp, 1.67698043002308E-1_wp,&
      &-6.23739039114127E-2_wp,-1.97051591839775E-1_wp,-1.22265382825162E-1_wp,&
      &-1.64257857747057E-2_wp, 2.27490812781041E-2_wp,-1.84302388618745E-2_wp,&
      &-8.56553766483229E-4_wp, 2.74842733013248E-3_wp, 2.55949906672251E-3_wp,&
      &-1.61115426160772E-1_wp, 1.66196888766802E-1_wp,-1.20546783549511E-1_wp,&
      & 9.88324537231105E-2_wp, 2.50428455703373E-3_wp,-2.88635579792292E-3_wp,&
      & 8.34445890794251E-2_wp,-1.34421662914771E-1_wp, 3.99666225357878E-2_wp,&
      &-4.83529798793855E-2_wp, 5.74113854317782E-3_wp, 5.87539955830004E-2_wp,&
      &-2.97630886764780E-3_wp, 2.19742923198916E-1_wp, 6.43178710881168E-2_wp,&
      & 9.28242957167868E-2_wp,-6.85043369046399E-2_wp, 2.00259137721060E-1_wp,&
      & 7.09776971104250E-5_wp, 2.03464908299606E-2_wp, 1.45200228907336E-2_wp,&
      &-2.15470436895541E-2_wp,-1.67828530829147E-2_wp,-3.81704914987720E-2_wp,&
      &-6.92951149767597E-2_wp, 9.85295808358629E-2_wp,-4.67219521543993E-3_wp,&
      &-1.41470268241085E-2_wp, 1.89740357101272E-2_wp,-6.41196409236768E-3_wp,&
      &-8.82996388817125E-2_wp,-1.65518095011688E-2_wp,-1.35247354235483E-1_wp,&
      &-2.61342733605818E-1_wp,-6.17705247109865E-2_wp, 4.23216561956390E-1_wp]

    real(wp), parameter :: qpat(*) = [&
    & 4.51357475089842E-2_wp,-4.17318742479749E-2_wp,-1.36136078209980E-2_wp,&
      & 9.12863796761614E-2_wp, 7.51195530423392E-2_wp,-3.15221396879885E-2_wp,&
      &-9.80963767772843E-2_wp, 1.45280277246739E-1_wp,-5.11292962085315E-1_wp,&
      &-1.16185111757740E-2_wp,-2.61668573633159E-1_wp, 6.09389338862585E-1_wp,&
      & 6.52037683652264E-3_wp, 3.69171710073827E-2_wp,-1.39436488066209E-2_wp,&
      &-2.83004471610422E-2_wp, 3.62439775801960E-2_wp, 7.42327197009803E-3_wp,&
      &-3.86296110845069E-2_wp,-1.21112119157970E-2_wp, 9.97808370410505E-3_wp,&
      &-1.91278825077846E-2_wp, 6.47092359099180E-2_wp, 2.86515273804006E-2_wp,&
      & 8.56200033892328E-1_wp, 9.74235458476180E-1_wp, 1.18086922059802E-1_wp,&
      &-8.71949362470118E-1_wp, 2.97351048258474E-1_wp,-9.74286955952133E-1_wp,&
      &-2.09879239083559E-1_wp, 9.31236036478151E-2_wp, 5.30742905878537E-2_wp,&
      & 2.06227339580677E-2_wp,-7.56884351806246E-3_wp, 1.56804948495703E-1_wp,&
      & 2.93163202989695E-2_wp, 1.59225627857236E-1_wp,-1.42600503718415E-1_wp,&
      &-3.26169646892855E-2_wp, 6.57594862333770E-2_wp, 1.13284183419448E-1_wp,&
      &-3.36105260187027E-1_wp,-6.53815224929778E-1_wp, 5.72904099232633E-2_wp,&
      & 3.91795155531810E-1_wp, 1.49074708451067E-1_wp, 2.78814850263761E-1_wp,&
      & 1.39849278576512E-1_wp, 1.62052454092380E-2_wp,-2.52552662853771E-1_wp,&
      & 4.12583681109575E-3_wp,-1.03140151821042E-1_wp, 1.12703384277264E-1_wp,&
      & 1.14113595393239E-1_wp, 7.38176190606366E-2_wp, 1.63093603236892E-1_wp,&
      &-2.24757223472586E-1_wp, 1.63465059762335E-1_wp,-2.77207198630140E-1_wp,&
      & 4.49629181418368E-1_wp, 4.25710033896662E-1_wp,-4.29670969707868E-1_wp,&
      &-8.92661030885367E-1_wp,-1.88058091049599E-1_wp,-1.99582117104920E-2_wp,&
      & 7.59715738599860E-1_wp, 8.23127259363381E-1_wp,-1.17603646025185E+0_wp,&
      &-1.44601626881170E+0_wp, 6.22081607059793E-1_wp, 4.16320721651909E-1_wp,&
      &-3.31955764455255E-1_wp, 1.15087423515139E+0_wp,-2.25645212472335E-1_wp,&
      & 1.18140494652646E+0_wp,-7.02380712786817E-1_wp, 5.57600976927600E-1_wp,&
      & 5.67085367397220E-2_wp, 1.13737943078387E-1_wp,-1.02296461455588E-2_wp,&
      &-8.06936077116478E-2_wp, 5.38804500826868E-2_wp,-4.64788905941666E-2_wp,&
      & 2.11496816821796E-2_wp,-3.21916493532053E-2_wp, 1.22390310502235E-1_wp,&
      &-1.62481131000237E-1_wp,-4.62447548254617E-2_wp,-1.43539992184418E-1_wp,&
      & 4.32313899320035E-1_wp,-2.95562991660850E-1_wp,-1.10104884940963E+0_wp,&
      & 3.42643149336475E-1_wp,-6.71823237278679E-1_wp, 6.68734950089545E-1_wp]

   real(wp), parameter :: feps = 80.0_wp

   call get_structure(mol, "MB16-43", "08")
   ! Test without multipoles
   input = alpb_input(feps, kernel=born_kernel%p16, alpb=.true., do_multipoles=.false.)
   call test_p(error, mol, input, qat, dpat, qpat)

   ! Test with multipoles
   input = alpb_input(feps, kernel=born_kernel%p16, alpb=.true., do_multipoles=.true.)
   call test_p(error, mol, input, qat, dpat, qpat)

end subroutine test_p_p16

subroutine test_p_alpb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(alpb_input) :: input

   real(wp), parameter :: qat(*) = [&
      &-3.28160099939119E-1_wp, 3.63839789415764E-1_wp,-9.39678438468329E-1_wp,&
      &-5.67353131753718E-1_wp, 3.91549236321241E-1_wp,-7.25527696006913E-1_wp,&
      & 2.81658498913997E-1_wp, 6.37609035711367E-1_wp,-3.32341795239006E-1_wp,&
      & 1.57147288398166E-1_wp,-2.15641708745235E-1_wp, 6.23132575019337E-1_wp,&
      & 6.86303747651841E-1_wp,-5.41860284045708E-1_wp, 2.80264616247842E-1_wp,&
      & 2.29058416549171E-1_wp]

   real(wp), parameter :: dpat(*) = [&
      &-1.02866383946914E-1_wp,-6.59701808627942E-2_wp, 1.67698043002308E-1_wp,&
      &-6.23739039114127E-2_wp,-1.97051591839775E-1_wp,-1.22265382825162E-1_wp,&
      &-1.64257857747057E-2_wp, 2.27490812781041E-2_wp,-1.84302388618745E-2_wp,&
      &-8.56553766483229E-4_wp, 2.74842733013248E-3_wp, 2.55949906672251E-3_wp,&
      &-1.61115426160772E-1_wp, 1.66196888766802E-1_wp,-1.20546783549511E-1_wp,&
      & 9.88324537231105E-2_wp, 2.50428455703373E-3_wp,-2.88635579792292E-3_wp,&
      & 8.34445890794251E-2_wp,-1.34421662914771E-1_wp, 3.99666225357878E-2_wp,&
      &-4.83529798793855E-2_wp, 5.74113854317782E-3_wp, 5.87539955830004E-2_wp,&
      &-2.97630886764780E-3_wp, 2.19742923198916E-1_wp, 6.43178710881168E-2_wp,&
      & 9.28242957167868E-2_wp,-6.85043369046399E-2_wp, 2.00259137721060E-1_wp,&
      & 7.09776971104250E-5_wp, 2.03464908299606E-2_wp, 1.45200228907336E-2_wp,&
      &-2.15470436895541E-2_wp,-1.67828530829147E-2_wp,-3.81704914987720E-2_wp,&
      &-6.92951149767597E-2_wp, 9.85295808358629E-2_wp,-4.67219521543993E-3_wp,&
      &-1.41470268241085E-2_wp, 1.89740357101272E-2_wp,-6.41196409236768E-3_wp,&
      &-8.82996388817125E-2_wp,-1.65518095011688E-2_wp,-1.35247354235483E-1_wp,&
      &-2.61342733605818E-1_wp,-6.17705247109865E-2_wp, 4.23216561956390E-1_wp]

    real(wp), parameter :: qpat(*) = [&
    & 4.51357475089842E-2_wp,-4.17318742479749E-2_wp,-1.36136078209980E-2_wp,&
      & 9.12863796761614E-2_wp, 7.51195530423392E-2_wp,-3.15221396879885E-2_wp,&
      &-9.80963767772843E-2_wp, 1.45280277246739E-1_wp,-5.11292962085315E-1_wp,&
      &-1.16185111757740E-2_wp,-2.61668573633159E-1_wp, 6.09389338862585E-1_wp,&
      & 6.52037683652264E-3_wp, 3.69171710073827E-2_wp,-1.39436488066209E-2_wp,&
      &-2.83004471610422E-2_wp, 3.62439775801960E-2_wp, 7.42327197009803E-3_wp,&
      &-3.86296110845069E-2_wp,-1.21112119157970E-2_wp, 9.97808370410505E-3_wp,&
      &-1.91278825077846E-2_wp, 6.47092359099180E-2_wp, 2.86515273804006E-2_wp,&
      & 8.56200033892328E-1_wp, 9.74235458476180E-1_wp, 1.18086922059802E-1_wp,&
      &-8.71949362470118E-1_wp, 2.97351048258474E-1_wp,-9.74286955952133E-1_wp,&
      &-2.09879239083559E-1_wp, 9.31236036478151E-2_wp, 5.30742905878537E-2_wp,&
      & 2.06227339580677E-2_wp,-7.56884351806246E-3_wp, 1.56804948495703E-1_wp,&
      & 2.93163202989695E-2_wp, 1.59225627857236E-1_wp,-1.42600503718415E-1_wp,&
      &-3.26169646892855E-2_wp, 6.57594862333770E-2_wp, 1.13284183419448E-1_wp,&
      &-3.36105260187027E-1_wp,-6.53815224929778E-1_wp, 5.72904099232633E-2_wp,&
      & 3.91795155531810E-1_wp, 1.49074708451067E-1_wp, 2.78814850263761E-1_wp,&
      & 1.39849278576512E-1_wp, 1.62052454092380E-2_wp,-2.52552662853771E-1_wp,&
      & 4.12583681109575E-3_wp,-1.03140151821042E-1_wp, 1.12703384277264E-1_wp,&
      & 1.14113595393239E-1_wp, 7.38176190606366E-2_wp, 1.63093603236892E-1_wp,&
      &-2.24757223472586E-1_wp, 1.63465059762335E-1_wp,-2.77207198630140E-1_wp,&
      & 4.49629181418368E-1_wp, 4.25710033896662E-1_wp,-4.29670969707868E-1_wp,&
      &-8.92661030885367E-1_wp,-1.88058091049599E-1_wp,-1.99582117104920E-2_wp,&
      & 7.59715738599860E-1_wp, 8.23127259363381E-1_wp,-1.17603646025185E+0_wp,&
      &-1.44601626881170E+0_wp, 6.22081607059793E-1_wp, 4.16320721651909E-1_wp,&
      &-3.31955764455255E-1_wp, 1.15087423515139E+0_wp,-2.25645212472335E-1_wp,&
      & 1.18140494652646E+0_wp,-7.02380712786817E-1_wp, 5.57600976927600E-1_wp,&
      & 5.67085367397220E-2_wp, 1.13737943078387E-1_wp,-1.02296461455588E-2_wp,&
      &-8.06936077116478E-2_wp, 5.38804500826868E-2_wp,-4.64788905941666E-2_wp,&
      & 2.11496816821796E-2_wp,-3.21916493532053E-2_wp, 1.22390310502235E-1_wp,&
      &-1.62481131000237E-1_wp,-4.62447548254617E-2_wp,-1.43539992184418E-1_wp,&
      & 4.32313899320035E-1_wp,-2.95562991660850E-1_wp,-1.10104884940963E+0_wp,&
      & 3.42643149336475E-1_wp,-6.71823237278679E-1_wp, 6.68734950089545E-1_wp]

   call get_structure(mol, "MB16-43", "08")
   solvent = get_solvent_data("water")

   ! Test without multipoles
   input = alpb_input(solvent%eps, solvent=solvent%solvent, &
      & kernel=born_kernel%p16, alpb=.true., do_multipoles=.false.)
   call test_p(error, mol, input, qat, dpat, qpat, method='gfn2')

   ! Test with multipoles
   input = alpb_input(solvent%eps, solvent=solvent%solvent, &
      & kernel=born_kernel%p16, alpb=.true., do_multipoles=.true.)
   call test_p(error, mol, input, qat, dpat, qpat, method='gfn2')

end subroutine test_p_alpb


subroutine test_p_still(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(alpb_input) :: input

   real(wp), parameter :: qat(*) = [&
      &-3.28160099939119E-1_wp, 3.63839789415764E-1_wp,-9.39678438468329E-1_wp,&
      &-5.67353131753718E-1_wp, 3.91549236321241E-1_wp,-7.25527696006913E-1_wp,&
      & 2.81658498913997E-1_wp, 6.37609035711367E-1_wp,-3.32341795239006E-1_wp,&
      & 1.57147288398166E-1_wp,-2.15641708745235E-1_wp, 6.23132575019337E-1_wp,&
      & 6.86303747651841E-1_wp,-5.41860284045708E-1_wp, 2.80264616247842E-1_wp,&
      & 2.29058416549171E-1_wp]

   real(wp), parameter :: dpat(*) = [&
      &-1.02866383946914E-1_wp,-6.59701808627942E-2_wp, 1.67698043002308E-1_wp,&
      &-6.23739039114127E-2_wp,-1.97051591839775E-1_wp,-1.22265382825162E-1_wp,&
      &-1.64257857747057E-2_wp, 2.27490812781041E-2_wp,-1.84302388618745E-2_wp,&
      &-8.56553766483229E-4_wp, 2.74842733013248E-3_wp, 2.55949906672251E-3_wp,&
      &-1.61115426160772E-1_wp, 1.66196888766802E-1_wp,-1.20546783549511E-1_wp,&
      & 9.88324537231105E-2_wp, 2.50428455703373E-3_wp,-2.88635579792292E-3_wp,&
      & 8.34445890794251E-2_wp,-1.34421662914771E-1_wp, 3.99666225357878E-2_wp,&
      &-4.83529798793855E-2_wp, 5.74113854317782E-3_wp, 5.87539955830004E-2_wp,&
      &-2.97630886764780E-3_wp, 2.19742923198916E-1_wp, 6.43178710881168E-2_wp,&
      & 9.28242957167868E-2_wp,-6.85043369046399E-2_wp, 2.00259137721060E-1_wp,&
      & 7.09776971104250E-5_wp, 2.03464908299606E-2_wp, 1.45200228907336E-2_wp,&
      &-2.15470436895541E-2_wp,-1.67828530829147E-2_wp,-3.81704914987720E-2_wp,&
      &-6.92951149767597E-2_wp, 9.85295808358629E-2_wp,-4.67219521543993E-3_wp,&
      &-1.41470268241085E-2_wp, 1.89740357101272E-2_wp,-6.41196409236768E-3_wp,&
      &-8.82996388817125E-2_wp,-1.65518095011688E-2_wp,-1.35247354235483E-1_wp,&
      &-2.61342733605818E-1_wp,-6.17705247109865E-2_wp, 4.23216561956390E-1_wp]

    real(wp), parameter :: qpat(*) = [&
    & 4.51357475089842E-2_wp,-4.17318742479749E-2_wp,-1.36136078209980E-2_wp,&
      & 9.12863796761614E-2_wp, 7.51195530423392E-2_wp,-3.15221396879885E-2_wp,&
      &-9.80963767772843E-2_wp, 1.45280277246739E-1_wp,-5.11292962085315E-1_wp,&
      &-1.16185111757740E-2_wp,-2.61668573633159E-1_wp, 6.09389338862585E-1_wp,&
      & 6.52037683652264E-3_wp, 3.69171710073827E-2_wp,-1.39436488066209E-2_wp,&
      &-2.83004471610422E-2_wp, 3.62439775801960E-2_wp, 7.42327197009803E-3_wp,&
      &-3.86296110845069E-2_wp,-1.21112119157970E-2_wp, 9.97808370410505E-3_wp,&
      &-1.91278825077846E-2_wp, 6.47092359099180E-2_wp, 2.86515273804006E-2_wp,&
      & 8.56200033892328E-1_wp, 9.74235458476180E-1_wp, 1.18086922059802E-1_wp,&
      &-8.71949362470118E-1_wp, 2.97351048258474E-1_wp,-9.74286955952133E-1_wp,&
      &-2.09879239083559E-1_wp, 9.31236036478151E-2_wp, 5.30742905878537E-2_wp,&
      & 2.06227339580677E-2_wp,-7.56884351806246E-3_wp, 1.56804948495703E-1_wp,&
      & 2.93163202989695E-2_wp, 1.59225627857236E-1_wp,-1.42600503718415E-1_wp,&
      &-3.26169646892855E-2_wp, 6.57594862333770E-2_wp, 1.13284183419448E-1_wp,&
      &-3.36105260187027E-1_wp,-6.53815224929778E-1_wp, 5.72904099232633E-2_wp,&
      & 3.91795155531810E-1_wp, 1.49074708451067E-1_wp, 2.78814850263761E-1_wp,&
      & 1.39849278576512E-1_wp, 1.62052454092380E-2_wp,-2.52552662853771E-1_wp,&
      & 4.12583681109575E-3_wp,-1.03140151821042E-1_wp, 1.12703384277264E-1_wp,&
      & 1.14113595393239E-1_wp, 7.38176190606366E-2_wp, 1.63093603236892E-1_wp,&
      &-2.24757223472586E-1_wp, 1.63465059762335E-1_wp,-2.77207198630140E-1_wp,&
      & 4.49629181418368E-1_wp, 4.25710033896662E-1_wp,-4.29670969707868E-1_wp,&
      &-8.92661030885367E-1_wp,-1.88058091049599E-1_wp,-1.99582117104920E-2_wp,&
      & 7.59715738599860E-1_wp, 8.23127259363381E-1_wp,-1.17603646025185E+0_wp,&
      &-1.44601626881170E+0_wp, 6.22081607059793E-1_wp, 4.16320721651909E-1_wp,&
      &-3.31955764455255E-1_wp, 1.15087423515139E+0_wp,-2.25645212472335E-1_wp,&
      & 1.18140494652646E+0_wp,-7.02380712786817E-1_wp, 5.57600976927600E-1_wp,&
      & 5.67085367397220E-2_wp, 1.13737943078387E-1_wp,-1.02296461455588E-2_wp,&
      &-8.06936077116478E-2_wp, 5.38804500826868E-2_wp,-4.64788905941666E-2_wp,&
      & 2.11496816821796E-2_wp,-3.21916493532053E-2_wp, 1.22390310502235E-1_wp,&
      &-1.62481131000237E-1_wp,-4.62447548254617E-2_wp,-1.43539992184418E-1_wp,&
      & 4.32313899320035E-1_wp,-2.95562991660850E-1_wp,-1.10104884940963E+0_wp,&
      & 3.42643149336475E-1_wp,-6.71823237278679E-1_wp, 6.68734950089545E-1_wp]

   real(wp), parameter :: feps = 80.0_wp

   call get_structure(mol, "MB16-43", "09")

   ! Test without multipoles
   input = alpb_input(feps, kernel=born_kernel%still, alpb=.true., do_multipoles=.false.)
   call test_p(error, mol, input, qat, dpat, qpat)

   ! Test with multipoles
   input = alpb_input(feps, kernel=born_kernel%still, alpb=.true., do_multipoles=.true.)
   call test_p(error, mol, input, qat, dpat, qpat)

end subroutine test_p_still

subroutine test_p_gbsa(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   type(alpb_input) :: input
    real(wp), parameter :: qat(*) = [&
      &-3.28160099939119E-1_wp, 3.63839789415764E-1_wp,-9.39678438468329E-1_wp,&
      &-5.67353131753718E-1_wp, 3.91549236321241E-1_wp,-7.25527696006913E-1_wp,&
      & 2.81658498913997E-1_wp, 6.37609035711367E-1_wp,-3.32341795239006E-1_wp,&
      & 1.57147288398166E-1_wp,-2.15641708745235E-1_wp, 6.23132575019337E-1_wp,&
      & 6.86303747651841E-1_wp,-5.41860284045708E-1_wp, 2.80264616247842E-1_wp,&
      & 2.29058416549171E-1_wp]

   real(wp), parameter :: dpat(*) = [&
      &-1.02866383946914E-1_wp,-6.59701808627942E-2_wp, 1.67698043002308E-1_wp,&
      &-6.23739039114127E-2_wp,-1.97051591839775E-1_wp,-1.22265382825162E-1_wp,&
      &-1.64257857747057E-2_wp, 2.27490812781041E-2_wp,-1.84302388618745E-2_wp,&
      &-8.56553766483229E-4_wp, 2.74842733013248E-3_wp, 2.55949906672251E-3_wp,&
      &-1.61115426160772E-1_wp, 1.66196888766802E-1_wp,-1.20546783549511E-1_wp,&
      & 9.88324537231105E-2_wp, 2.50428455703373E-3_wp,-2.88635579792292E-3_wp,&
      & 8.34445890794251E-2_wp,-1.34421662914771E-1_wp, 3.99666225357878E-2_wp,&
      &-4.83529798793855E-2_wp, 5.74113854317782E-3_wp, 5.87539955830004E-2_wp,&
      &-2.97630886764780E-3_wp, 2.19742923198916E-1_wp, 6.43178710881168E-2_wp,&
      & 9.28242957167868E-2_wp,-6.85043369046399E-2_wp, 2.00259137721060E-1_wp,&
      & 7.09776971104250E-5_wp, 2.03464908299606E-2_wp, 1.45200228907336E-2_wp,&
      &-2.15470436895541E-2_wp,-1.67828530829147E-2_wp,-3.81704914987720E-2_wp,&
      &-6.92951149767597E-2_wp, 9.85295808358629E-2_wp,-4.67219521543993E-3_wp,&
      &-1.41470268241085E-2_wp, 1.89740357101272E-2_wp,-6.41196409236768E-3_wp,&
      &-8.82996388817125E-2_wp,-1.65518095011688E-2_wp,-1.35247354235483E-1_wp,&
      &-2.61342733605818E-1_wp,-6.17705247109865E-2_wp, 4.23216561956390E-1_wp]

    real(wp), parameter :: qpat(*) = [&
    & 4.51357475089842E-2_wp,-4.17318742479749E-2_wp,-1.36136078209980E-2_wp,&
      & 9.12863796761614E-2_wp, 7.51195530423392E-2_wp,-3.15221396879885E-2_wp,&
      &-9.80963767772843E-2_wp, 1.45280277246739E-1_wp,-5.11292962085315E-1_wp,&
      &-1.16185111757740E-2_wp,-2.61668573633159E-1_wp, 6.09389338862585E-1_wp,&
      & 6.52037683652264E-3_wp, 3.69171710073827E-2_wp,-1.39436488066209E-2_wp,&
      &-2.83004471610422E-2_wp, 3.62439775801960E-2_wp, 7.42327197009803E-3_wp,&
      &-3.86296110845069E-2_wp,-1.21112119157970E-2_wp, 9.97808370410505E-3_wp,&
      &-1.91278825077846E-2_wp, 6.47092359099180E-2_wp, 2.86515273804006E-2_wp,&
      & 8.56200033892328E-1_wp, 9.74235458476180E-1_wp, 1.18086922059802E-1_wp,&
      &-8.71949362470118E-1_wp, 2.97351048258474E-1_wp,-9.74286955952133E-1_wp,&
      &-2.09879239083559E-1_wp, 9.31236036478151E-2_wp, 5.30742905878537E-2_wp,&
      & 2.06227339580677E-2_wp,-7.56884351806246E-3_wp, 1.56804948495703E-1_wp,&
      & 2.93163202989695E-2_wp, 1.59225627857236E-1_wp,-1.42600503718415E-1_wp,&
      &-3.26169646892855E-2_wp, 6.57594862333770E-2_wp, 1.13284183419448E-1_wp,&
      &-3.36105260187027E-1_wp,-6.53815224929778E-1_wp, 5.72904099232633E-2_wp,&
      & 3.91795155531810E-1_wp, 1.49074708451067E-1_wp, 2.78814850263761E-1_wp,&
      & 1.39849278576512E-1_wp, 1.62052454092380E-2_wp,-2.52552662853771E-1_wp,&
      & 4.12583681109575E-3_wp,-1.03140151821042E-1_wp, 1.12703384277264E-1_wp,&
      & 1.14113595393239E-1_wp, 7.38176190606366E-2_wp, 1.63093603236892E-1_wp,&
      &-2.24757223472586E-1_wp, 1.63465059762335E-1_wp,-2.77207198630140E-1_wp,&
      & 4.49629181418368E-1_wp, 4.25710033896662E-1_wp,-4.29670969707868E-1_wp,&
      &-8.92661030885367E-1_wp,-1.88058091049599E-1_wp,-1.99582117104920E-2_wp,&
      & 7.59715738599860E-1_wp, 8.23127259363381E-1_wp,-1.17603646025185E+0_wp,&
      &-1.44601626881170E+0_wp, 6.22081607059793E-1_wp, 4.16320721651909E-1_wp,&
      &-3.31955764455255E-1_wp, 1.15087423515139E+0_wp,-2.25645212472335E-1_wp,&
      & 1.18140494652646E+0_wp,-7.02380712786817E-1_wp, 5.57600976927600E-1_wp,&
      & 5.67085367397220E-2_wp, 1.13737943078387E-1_wp,-1.02296461455588E-2_wp,&
      &-8.06936077116478E-2_wp, 5.38804500826868E-2_wp,-4.64788905941666E-2_wp,&
      & 2.11496816821796E-2_wp,-3.21916493532053E-2_wp, 1.22390310502235E-1_wp,&
      &-1.62481131000237E-1_wp,-4.62447548254617E-2_wp,-1.43539992184418E-1_wp,&
      & 4.32313899320035E-1_wp,-2.95562991660850E-1_wp,-1.10104884940963E+0_wp,&
      & 3.42643149336475E-1_wp,-6.71823237278679E-1_wp, 6.68734950089545E-1_wp]

   call get_structure(mol, "MB16-43", "09")
   solvent = get_solvent_data("water")

   ! Test GFN2/GBSA with multipoles off
   input = alpb_input(solvent%eps, solvent=solvent%solvent, &
      & kernel=born_kernel%still, alpb=.false., do_multipoles=.false.)
   call test_p(error, mol, input, qat, dpat, qpat, method='gfn2')

   ! Test with multipoles
   input = alpb_input(solvent%eps, solvent=solvent%solvent, &
      & kernel=born_kernel%still, alpb=.false., do_multipoles=.true.)
   call test_p(error, mol, input, qat, dpat, qpat, method='gfn2')

end subroutine test_p_gbsa


subroutine test_unsupported_solvent(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(solvent_data) :: solvent
   real(wp), parameter :: qat(*) = [&
      &-8.99890404486076E-2_wp, 9.42168087556583E-2_wp,-1.49387631509499E-1_wp, &
      &-2.99114121895542E-1_wp, 4.85527734875224E-1_wp,-6.83156326406137E-2_wp, &
      & 1.50011293889337E-2_wp, 2.79368544459465E-1_wp,-1.24072452878322E-1_wp, &
      &-9.36760994051244E-2_wp,-2.19062123031622E-1_wp, 2.14538817685587E-1_wp, &
      & 3.06156726072831E-1_wp,-3.86105514712244E-1_wp,-1.51265171389388E-3_wp, &
      & 3.64255069977693E-2_wp]
   type(alpb_input) :: input

   call get_structure(mol, "MB16-43", "04")

   ! Check GFN1/GBSA unavailable solvent
   solvent = get_solvent_data("aniline")
   input = alpb_input(solvent%eps, solvent=solvent%solvent, &
      & kernel=born_kernel%still, alpb=.false., do_multipoles=.false.)
   call get_alpb_param(input, mol, 'gfn1', error)

end subroutine test_unsupported_solvent


end module test_solvation_born
