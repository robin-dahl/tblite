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
   use tblite_container, only : container_cache
   use tblite_scf_potential, only : potential_type
   use tblite_solvation_ddx, only : ddx_solvation, ddx_solvation_model, ddx_input 
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, eeq_guess
   use tblite_basis_type, only : basis_type
   use tblite_integral_type, only : integral_type
   use tblite_context_type, only : context_type
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

   testsuite = [ &
      ! new_unittest("energy-mol-cosmo", test_e_cosmo_m01), &
      ! new_unittest("energy-mol-pcm", test_e_pcm_m01), &
      ! new_unittest("energy-mol-lpb", test_e_lpb_m01), &
      ! new_unittest("gradient-mol-num-cosmo", test_g_num_cosmo_m02), &
      ! new_unittest("gradient-mol-cosmo", test_g_cosmo_m02) &
      ! new_unittest("gradient-mol-num-pcm", test_g_pcm_m02), &
      ! new_unittest("gradient-mol-lpb", test_g_lpb_m02), &
      new_unittest("potential-mol-cosmo", test_p_cosmo_m03) &
      ! new_unittest("potential-mol-pcm", test_p_pcm_m03) &
      ! new_unittest("potential-mol-lpb", test_p_lpb_m03) &
      ]

end subroutine collect_solvation_ddx


subroutine test_e(error, model, mol, qat, ref, kappa)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Solvation model (COSMO=11, CPCM=12, PCM=2, LPB=3)
   integer, intent(in) :: model

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Reference energy
   real(wp), intent(in) :: ref

   !> Debye-Hückel screening parameter (only used in LPB)
   real(wp), optional, intent(in) :: kappa

   type(ddx_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: feps = 80.0_wp, rscale = 1.0_wp
   integer, parameter :: nang = 302
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp) :: energy(mol%nat)

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1), source=0.0_wp)
   energy = 0.0_wp

   if (present(kappa)) then
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale, kappa=kappa))
   else
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale))
   end if

   call solv%update(mol, cache)
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)

   if (abs(sum(energy) - ref) > thr) then
      call test_failed(error, "Energy does not match reference")
      print *, sum(energy)
   end if
end subroutine test_e


subroutine test_g_num(error, model, mol, qat, kappa)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Solvation model (COSMO=11, CPCM=12, PCM=2, LPB=3)
   integer, intent(in) :: model

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Debye-Hückel screening parameter (only used in LPB)
   real(wp), optional, intent(in) :: kappa

   type(ddx_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: feps = 80.0_wp, rscale = 1.0_wp
   integer, parameter :: nang = 302
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp), allocatable :: gradient(:, :), numg(:, :)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3)
   integer :: ii, ic

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))

   if (present(kappa)) then
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale, kappa=kappa))
   else
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale))
   end if

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
end subroutine test_g_num

subroutine test_g(error, model, mol, qat, ref, kappa)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Solvation model (COSMO=11, CPCM=12, PCM=2, LPB=3)
   integer, intent(in) :: model

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Reference gradient
   real(wp), intent(in) :: ref(:, :)

   !> Debye-Hückel screening parameter (only used in LPB)
   real(wp), optional, intent(in) :: kappa

   type(ddx_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: feps = 78.0_wp, rscale = 1.0_wp
   integer, parameter :: nang = 302
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp), allocatable :: gradient(:, :), numg(:, :)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3)
   integer :: ii, ic

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))


   if (present(kappa)) then
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale, kappa=kappa))
   else
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale))
   end if

   allocate(gradient(3, mol%nat))

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp

   call solv%update(mol, cache)
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)
   call solv%get_gradient(mol, cache, wfn, gradient, sigma)

   if (any(abs(gradient - ref) > thr)) then
      call test_failed(error, "Gradient does not match")
      print '(3es20.13)', gradient
      print '(a)', "---"
      print '(3es20.13)', ref
      print '(a)', "---"
      print '(3es20.13)', gradient - ref
   end if
end subroutine test_g


subroutine test_p(error, model, mol, qat, dpat, kappa)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Solvation model (COSMO=11, CPCM=12, PCM=2, LPB=3)
   integer, intent(in) :: model

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Atomic dipole
   real(wp), intent(in) :: dpat(:, :)

   !> Debye-Hückel screening parameter (only used in LPB)
   real(wp), optional, intent(in) :: kappa

   type(ddx_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache), allocatable :: cache
   real(wp), parameter :: feps = 80.0_wp, rscale = 1.0_wp
   integer, parameter :: nang = 302
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr = 1e+3_wp*sqrt(epsilon(1.0_wp))
   real(wp), allocatable :: vat(:)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat)
   integer :: ii

   wfn%qat = reshape(qat, [size(qat), 1])
   wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])
   allocate(pot%vat(size(qat, 1), 1), source=0.0_wp)
   allocate(pot%vdp(3, size(qat, 1), 1), source=0.0_wp)

   if (present(kappa)) then
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale, kappa=kappa))
   else
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale))
   end if

   allocate(cache)
   call solv%update(mol, cache)

   allocate(vat(mol%nat), source=0.0_wp)
   do ii = 1, mol%nat
      er = 0.0_wp
      el = 0.0_wp
      wfn%qat(ii, 1) = wfn%qat(ii, 1) + step
      wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])
      call solv%get_potential(mol, cache, wfn, pot)
      call solv%get_energy(mol, cache, wfn, er)

      wfn%qat(ii, 1) = wfn%qat(ii, 1) - 2*step
      wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])
      call solv%get_potential(mol, cache, wfn, pot)
      call solv%get_energy(mol, cache, wfn, el)

      wfn%qat(ii, 1) = wfn%qat(ii, 1) + step
      wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])
      vat(ii) = 0.5_wp*(sum(er) - sum(el))/step
   end do

   energy = 0.0_wp
   pot%vat(:, :) = 0.0_wp
   pot%vdp(:, :, :) = 0.0_wp
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)

   if (any(abs([pot%vat] - vat) > thr)) then
      call test_failed(error, "Potential does not match")
      print '(a)', 'analytical'
      print '(3es20.13)', pot%vat
      print '(a)', "---"
      print '(a)', 'numerical'
      print '(3es20.13)', vat
      print '(a)', "---"
      print '(a)', 'diff'
      print '(3es20.13)', [pot%vat] - vat
   end if
end subroutine test_p

subroutine test_dp(error, model, mol, qat, dpat, kappa)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Solvation model (COSMO=11, CPCM=12, PCM=2, LPB=3)
   integer, intent(in) :: model

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Atomic partial charges (kept constant here)
   real(wp), intent(in) :: qat(:)
   !> Atomic dipoles (varied here)
   real(wp), intent(in) :: dpat(:, :)

   !> Debye-Hückel screening parameter (only used in LPB)
   real(wp), optional, intent(in) :: kappa

   type(ddx_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache), allocatable :: cache
   real(wp), parameter :: feps = 80.0_wp, rscale = 1.0_wp
   integer, parameter :: nang = 302
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr  = 1e+3_wp*sqrt(epsilon(1.0_wp))
   real(wp) :: vdp(3, mol%nat)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat)
   integer :: ii, k

   ! Set baseline wfn to the *given* monopoles and dipoles
   wfn%qat  = reshape(qat, [size(qat), 1])
   wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])

   allocate(pot%vat(size(qat, 1), 1), source=0.0_wp)
   allocate(pot%vdp(3, size(qat, 1), 1), source=0.0_wp)

   if (present(kappa)) then
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale, kappa=kappa))
   else
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale))
   end if

   allocate(cache)
   call solv%update(mol, cache)


   !--- Numerical vdp via central differences wrt dipole components ---
   vdp = 0.0_wp
   do ii = 1, mol%nat
      do k = 1, 3
         er = 0.0_wp
         el = 0.0_wp

         ! +step on component k of atom ii, keep everything else fixed
         wfn%qat  = reshape(qat,  [size(qat), 1])                 ! keep monopoles constant
         wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])
         wfn%dpat(k, ii, 1) = wfn%dpat(k, ii, 1) + step
         call solv%get_energy(mol, cache, wfn, er)

         ! -step on the same component
         wfn%qat  = reshape(qat,  [size(qat), 1])
         wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])
         wfn%dpat(k, ii, 1) = wfn%dpat(k, ii, 1) - step
         call solv%get_energy(mol, cache, wfn, el)

         ! central finite difference
         vdp(k, ii) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   !--- Analytical potentials/energy at the baseline (unshifted) wfn ---
   energy       = 0.0_wp
   wfn%qat      = reshape(qat,  [size(qat), 1])
   wfn%dpat     = reshape(dpat, [3, size(dpat, 2), 1])
   pot%vat(:, :)      = 0.0_wp
   pot%vdp(:, :, :)   = 0.0_wp
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)

   write(*,*) pot%vdp
   write(*,*) vdp
   stop

   !--- Compare analytical vs numerical dipole potential ---
   if (any(abs(pot%vdp(:, :, 1) - vdp) > thr)) then
      call test_failed(error, "Dipole potential does not match")
      print '(a)', 'analytical (pot%vdp(:, :, 1))'
      print '(3es20.13)', pot%vdp(:, :, 1)
      print '(a)', "---"
      print '(a)', 'numerical (vdp)'
      print '(3es20.13)', vdp
      print '(a)', "---"
      print '(a)', 'diff (analytical - numerical)'
      print '(3es20.13)', pot%vdp(:, :, 1) - vdp
   end if

end subroutine test_dp
subroutine test_qp_new(error, model, mol, qat, dpat, qpat, kappa)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Solvation model (COSMO=11, CPCM=12, PCM=2, LPB=3)
   integer, intent(in) :: model

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Atomic partial charges (kept constant here)
   real(wp), intent(in) :: qat(:)
   !> Atomic dipoles (kept constant here)
   real(wp), intent(in) :: dpat(:, :)
   !> Atomic quadrupoles (6 components: xx, xy, yy, xz, yz, zz)
   real(wp), intent(in) :: qpat(:, :)

   !> Debye-Hückel screening parameter (only used in LPB)
   real(wp), optional, intent(in) :: kappa

   type(ddx_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache), allocatable :: cache
   real(wp), parameter :: feps = 80.0_wp, rscale = 1.0_wp
   integer, parameter :: nang = 302
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: thr  = 1e+3_wp*sqrt(epsilon(1.0_wp))

   real(wp) :: vqp(6, mol%nat)
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat)
   integer :: ii, k

   ! Set baseline wfn to the given multipoles
   wfn%qat  = reshape(qat,  [size(qat), 1])
   wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])
   wfn%qpat = reshape(qpat, [6, size(qpat, 2), 1])

   allocate(pot%vat(size(qat, 1), 1), source=0.0_wp)
   allocate(pot%vdp(3, size(qat, 1), 1), source=0.0_wp)
   allocate(pot%vqp(6, size(qat, 1), 1), source=0.0_wp)

   if (present(kappa)) then
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale, kappa=kappa))
   else
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale))
   end if

   allocate(cache)
   call solv%update(mol, cache)

   !--- Numerical vqp via central differences wrt quadrupole components ---
   vqp = 0.0_wp
   do ii = 1, mol%nat
      do k = 1, 6
         er = 0.0_wp
         el = 0.0_wp

         ! +step on component k of atom ii, keep everything else fixed
         wfn%qat  = reshape(qat,  [size(qat), 1])
         wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])
         wfn%qpat = reshape(qpat, [6, size(qpat, 2), 1])
         wfn%qpat(k, ii, 1) = wfn%qpat(k, ii, 1) + step
         call solv%get_energy(mol, cache, wfn, er)

         ! -step on the same component
         wfn%qat  = reshape(qat,  [size(qat), 1])
         wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])
         wfn%qpat = reshape(qpat, [6, size(qpat, 2), 1])
         wfn%qpat(k, ii, 1) = wfn%qpat(k, ii, 1) - step
         call solv%get_energy(mol, cache, wfn, el)

         ! central finite difference
         vqp(k, ii) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   !--- Analytical potentials/energy at the baseline (unshifted) wfn ---
   energy            = 0.0_wp
   wfn%qat           = reshape(qat,  [size(qat), 1])
   wfn%dpat          = reshape(dpat, [3, size(dpat, 2), 1])
   wfn%qpat          = reshape(qpat, [6, size(qpat, 2), 1])
   pot%vat(:, :)     = 0.0_wp
   pot%vdp(:, :, :)  = 0.0_wp
   pot%vqp(:, :, :)  = 0.0_wp
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)


   !--- Compare analytical vs numerical quadrupole potential derivative ---
   if (any(abs(pot%vqp(:, :, 1) - vqp) > thr)) then
      call test_failed(error, "Quadrupole potential does not match")
      print '(a)', 'analytical (pot%vqp(:, :, 1))'
      print '(6es20.13)', pot%vqp(:, :, 1)
      print '(a)', "---"
      print '(a)', 'numerical (vqp)'
      print '(6es20.13)', vqp
      print '(a)', "---"
      print '(a)', 'diff (analytical - numerical)'
      print '(6es20.13)', pot%vqp(:, :, 1) - vqp
   end if

end subroutine test_qp_new

subroutine test_qp_traceless(error, model, mol, qat, dpat, qpat, kappa)
   implicit none
   !----------------- interfaces -----------------
   type(error_type),        allocatable, intent(out) :: error
   integer,                             intent(in)  :: model
   type(structure_type),                 intent(in)  :: mol
   real(wp),                            intent(in)  :: qat(:)
   real(wp),                            intent(in)  :: dpat(:, :)
   ! Traceless Cartesian quadrupoles, 6 per atom: (xx, xy, yy, xz, yz, zz)
   real(wp),                            intent(in)  :: qpat(:, :)
   real(wp), optional,                  intent(in)  :: kappa

   !----------------- solver containers -----------------
   type(ddx_solvation)                          :: solv
   type(container_cache), allocatable           :: cache
   type(wavefunction_type)                      :: wfn
   type(potential_type)                         :: pot

   !----------------- sizes & numerics -----------------
   integer                                      :: nat, nqp
   real(wp)                                     :: step, thr
   real(wp), parameter                          :: feps  = 80.0_wp
   real(wp), parameter                          :: rscale= 1.0_wp
   integer,  parameter                          :: nang  = 302

   !----------------- working arrays -----------------
   ! allocatables that we size at runtime:
   real(wp), allocatable                        :: vqp(:,:)   ! (6, nat)
   real(wp), allocatable                        :: energy(:), er(:), el(:)

   ! fixed-size local vectors/matrices for one atom:
   real(wp)                                     :: D(6,6)     ! 6 traceless directions
   real(wp)                                     :: s(6)       ! directional derivs along D(:,j)
   real(wp)                                     :: q0(6), qp(6), qm(6)

   integer                                      :: ii, j

   !================= shape checks & allocs =================
   nqp = size(qpat, 1)
   nat = mol%nat
   if (nqp /= 6) then
      call test_failed(error, "Expected 6 traceless Cartesian components (xx,xy,yy,xz,yz,zz).")
      return
   end if

   allocate(vqp(6, nat))
   allocate(energy(nat), er(nat), el(nat))

   step = 1.0e-4_wp
   thr  = 1.0e+3_wp * sqrt(epsilon(1.0_wp))

   !================= define 6 traceless directions =================
   ! Each column j is the perturbation for component j, preserving trace.
   ! D_xx
   D(:,1) = (/ +1.0_wp, 0.0_wp, -0.5_wp, 0.0_wp, 0.0_wp, -0.5_wp /)
   ! D_xy
   D(:,2) = (/  0.0_wp, +1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,  0.0_wp /)
   ! D_yy
   D(:,3) = (/ -0.5_wp, 0.0_wp, +1.0_wp, 0.0_wp, 0.0_wp, -0.5_wp /)
   ! D_xz
   D(:,4) = (/  0.0_wp, 0.0_wp, 0.0_wp, +1.0_wp, 0.0_wp,  0.0_wp /)
   ! D_yz
   D(:,5) = (/  0.0_wp, 0.0_wp, 0.0_wp,  0.0_wp, +1.0_wp, 0.0_wp /)
   ! D_zz
   D(:,6) = (/ -0.5_wp, 0.0_wp, -0.5_wp, 0.0_wp, 0.0_wp, +1.0_wp /)

   !================= build solver & baseline WFN =================
   if (present(kappa)) then
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale, kappa=kappa))
   else
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale))
   end if
   allocate(cache)
   call solv%update(mol, cache)

   wfn%qat  = reshape(qat,  [size(qat), 1])
   wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])
   wfn%qpat = reshape(qpat, [nqp, size(qpat, 2), 1])

   ! analytical potential container (limit to 6 traceless comps)
   allocate(pot%vat(size(qat,1), 1), source=0.0_wp)
   allocate(pot%vqp(6, size(qat,1), 1), source=0.0_wp)

   !================= numerical derivatives =================
   vqp = 0.0_wp
   do ii = 1, nat
      q0 = qpat(:, ii)

      do j = 1, 6
         ! +h along D(:,j)
         qp = q0 + step * D(:, j)
         wfn%qat  = reshape(qat,  [size(qat), 1])
         wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])
         wfn%qpat = reshape(qpat, [nqp, size(qpat, 2), 1])
         wfn%qpat(:, ii, 1) = qp
         call solv%get_energy(mol, cache, wfn, er)

         ! -h along D(:,j)
         qm = q0 - step * D(:, j)
         wfn%qat  = reshape(qat,  [size(qat), 1])
         wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])
         wfn%qpat = reshape(qpat, [nqp, size(qpat, 2), 1])
         wfn%qpat(:, ii, 1) = qm
         call solv%get_energy(mol, cache, wfn, el)

         s(j) = 0.5_wp * (sum(er) - sum(el)) / step
      end do

      ! Map back to component gradient.
      ! Off-diagonals are direct. For the diagonals:
      !   s_xx = (3/2) g_xx,  s_yy = (3/2) g_yy,  s_zz = (3/2) g_zz  (trace fixed)
      vqp(1, ii) = (2.0_wp/3.0_wp) * s(1)   ! dE/d q_xx
      vqp(2, ii) =               s(2)       ! dE/d q_xy
      vqp(3, ii) = (2.0_wp/3.0_wp) * s(3)   ! dE/d q_yy
      vqp(4, ii) =               s(4)       ! dE/d q_xz
      vqp(5, ii) =               s(5)       ! dE/d q_yz
      vqp(6, ii) = (2.0_wp/3.0_wp) * s(6)   ! dE/d q_zz
   end do

   !================= analytical & comparison =================
   wfn%qat  = reshape(qat,  [size(qat), 1])
   wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])
   wfn%qpat = reshape(qpat, [nqp, size(qpat, 2), 1])

   pot%vat(:, :)    = 0.0_wp
   pot%vqp(:, :, :) = 0.0_wp
   call solv%get_potential(mol, cache, wfn, pot)
   call solv%get_energy(mol, cache, wfn, energy)

   if (any(abs(pot%vqp(1:6, :, 1) - vqp) > thr)) then
      call test_failed(error, "Quadrupole potential (6-comp traceless) does not match numerical derivative.")
      print '(a)', 'analytical (pot%vqp(1:6, :, 1))'
      print '(3es20.13)', pot%vqp(1:6, :, 1)
      print '(a)', 'numerical (vqp)'
      print '(3es20.13)', vqp
      print '(a)', 'diff (analytical - numerical)'
      print '(3es20.13)', pot%vqp(1:6, :, 1) - vqp
   end if
end subroutine test_qp_traceless


! subroutine test_qp(error, model, mol, qat, dpat, qpat, kappa)

!    !> Error handling
!    type(error_type), allocatable, intent(out) :: error

!    !> Solvation model (COSMO=11, CPCM=12, PCM=2, LPB=3)
!    integer, intent(in) :: model

!    !> Molecular structure data
!    type(structure_type), intent(in) :: mol

!    !> Atomic partial charges (kept constant here)
!    real(wp), intent(in) :: qat(:)
!    !> Atomic dipoles (kept constant here)
!    real(wp), intent(in) :: dpat(:, :)
!    !> Atomic quadrupoles (varied here)
!    real(wp), intent(in) :: qpat(:, :)

!    !> Debye-Hückel screening parameter (only used in LPB)
!    real(wp), optional, intent(in) :: kappa

!    type(ddx_solvation) :: solv
!    type(wavefunction_type) :: wfn
!    type(potential_type)   :: pot
!    type(container_cache), allocatable :: cache

!    real(wp), parameter :: feps = 80.0_wp, rscale = 1.0_wp
!    integer,  parameter :: nang = 302
!    real(wp) :: step = 1.0e-4_wp
!    real(wp), parameter :: thr  = 1e+3_wp*sqrt(epsilon(1.0_wp))

!    integer :: ii, l, nqp, nat
!    real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat)

!    ! Numerical potential wrt quadrupoles
!    real(wp) :: vqp(size(qpat, 1), mol%nat)

!    ! Number of quadrupole components per atom (usually 6 or 9)
!    nqp = size(qpat, 1)
!    nat = mol%nat


!    ! Set baseline wfn to the given monopoles, dipoles, and quadrupoles
!    wfn%qat  = reshape(qat,  [size(qat), 1])
!    wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])
!    wfn%qpat = reshape(qpat, [nqp, size(qpat, 2), 1])

   
!    ! Allocate only what we need the solver to fill
!    allocate(pot%vat(size(qat, 1), 1), source=0.0_wp)
!    ! NOTE: container field name for quadrupole potential may differ in your codebase.
!    ! If your type uses a different name, rename pot%vqp accordingly.
!    allocate(pot%vqp(nqp, size(qat, 1), 1), source=0.0_wp)

!    if (present(kappa)) then
!       solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale, kappa=kappa))
!    else
!       solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale))
!    end if

!    allocate(cache)
!    call solv%update(mol, cache)

!    !--- Numerical vqp via central differences wrt quadrupole components ---
!    vqp = 0.0_wp
!    do ii = 1, nat
!       do l = 1, nqp
!          er = 0.0_wp
!          el = 0.0_wp

!          ! +step on quadrupole component l of atom ii; keep monopoles and dipoles fixed
!          wfn%qat  = reshape(qat,  [size(qat), 1])
!          wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])
!          wfn%qpat = reshape(qpat, [nqp, size(qpat, 2), 1])

!          wfn%qat  = reshape(qat,  [size(qat), 1])
!          wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])
!          wfn%qpat = reshape(qpat, [nqp, size(qpat, 2), 1])
!          wfn%qpat(l, ii, 1) = wfn%qpat(l, ii, 1) + step
!          call solv%get_energy(mol, cache, wfn, er)

!          ! -step on the same component
!          wfn%qat  = reshape(qat,  [size(qat), 1])
!          wfn%dpat = reshape(dpat, [3, size(dpat, 2), 1])
!          wfn%qpat = reshape(qpat, [nqp, size(qpat, 2), 1])
!          wfn%qpat(l, ii, 1) = wfn%qpat(l, ii, 1) - step
!          call solv%get_energy(mol, cache, wfn, el)

!          ! central finite difference
!          vqp(l, ii) = 0.5_wp*(sum(er) - sum(el))/step
!       end do
!    end do

!    !--- Analytical potentials/energy at the baseline (unshifted) wfn ---
!    energy         = 0.0_wp
!    wfn%qat        = reshape(qat,  [size(qat), 1])
!    wfn%dpat       = reshape(dpat, [3, size(dpat, 2), 1])
!    wfn%qpat       = reshape(qpat, [nqp, size(qpat, 2), 1])
!    pot%vat(:, :)  = 0.0_wp
!    pot%vqp(:, :, :) = 0.0_wp
!    call solv%get_potential(mol, cache, wfn, pot)
!    call solv%get_energy(mol, cache, wfn, energy)

!    !--- Compare analytical vs numerical quadrupole potential ---
!    if (any(abs(pot%vqp(:, :, 1) - vqp) > thr)) then
!       call test_failed(error, "Quadrupole potential does not match")
!       print '(a)', 'analytical (pot%vqp(:, :, 1))'
!       print '(3es20.13)', pot%vqp(:, :, 1)
!       print '(a)', "---"
!       print '(a)', 'numerical (vqp)'
!       print '(3es20.13)', vqp
!       print '(a)', "---"
!       print '(a)', 'diff (analytical - numerical)'
!       print '(3es20.13)', pot%vqp(:, :, 1) - vqp
!    end if

! end subroutine test_qp


! subroutine test_qp(error, model, mol, qat, dpat, qpat, kappa)
!    use mctc_env,                 only : wp, error_type
!    use mctc_io,                  only : structure_type
!    use tblite_solvation_ddx,     only : ddx_solvation, ddx_input
!    use tblite_wavefunction_type, only : wavefunction_type
!    use tblite_scf_potential,     only : potential_type
!    use tblite_container_cache,   only : container_cache
!    implicit none

!    ! I/O
!    type(error_type), allocatable, intent(out) :: error
!    integer,               intent(in) :: model                 ! COSMO=11, CPCM=12, PCM=2, LPB=3
!    type(structure_type),  intent(in) :: mol
!    real(wp),              intent(in) :: qat(:)                ! (nat)
!    real(wp),              intent(in) :: dpat(:, :)            ! (3, nat)
!    ! qpat must be (6, nat) in order (xx, xy, yy, xz, yz, zz).
!    ! Stored diagonals are traceless as a triplet (xx+yy+zz ≈ 0) in many implementations.
!    real(wp),              intent(in) :: qpat(:, :)
!    real(wp),    optional, intent(in) :: kappa                 ! LPB only

!    ! Types
!    type(ddx_solvation)                :: solv
!    type(wavefunction_type)            :: wfn
!    type(potential_type)               :: pot
!    type(container_cache), allocatable :: cache

!    ! Params
!    integer,  parameter :: nqp = 6, nang = 302
!    integer,  parameter :: ixx=1, ixy=2, iyy=3, ixz=4, iyz=5, izz=6
!    real(wp), parameter :: feps = 80.0_wp, rscale = 1.0_wp
!    real(wp), parameter :: step = 1.0e-4_wp
!    real(wp), parameter :: thr_abs = 1.0e+3_wp*sqrt(epsilon(1.0_wp))
!    real(wp), parameter :: thr_rel = 1.0e-7_wp

!    ! Locals
!    integer :: nat, ii, l
!    real(wp) :: er(mol%nat), el(mol%nat)
!    real(wp), allocatable :: vnum(:,:), vana(:,:)
!    logical :: offdiag_times_two

!    ! --- sanity ---
!    nat = mol%nat
!    if (size(qpat,1) /= nqp .or. size(qpat,2) /= nat) then
!       call test_failed(error, "qpat must be size (6, nat) in order (xx,xy,yy,xz,yz,zz).")
!       return
!    end if
!    if (size(dpat,1) /= 3 .or. size(dpat,2) /= nat) then
!       call test_failed(error, "dpat must be (3, nat).")
!       return
!    end if
!    if (size(qat,1) /= nat) then
!       call test_failed(error, "qat must be length nat.")
!       return
!    end if

!    ! --- baseline wavefunction ---
!    wfn%qat  = reshape(qat,  [nat, 1])
!    wfn%dpat = reshape(dpat, [3,   nat, 1])
!    wfn%qpat = reshape(qpat, [6,   nat, 1])

!    allocate(pot%vat(nat, 1), source=0.0_wp)
!    allocate(pot%vqp(6, nat, 1), source=0.0_wp)

!    if (present(kappa)) then
!       solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale, kappa=kappa))
!    else
!       solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, rscale=rscale))
!    end if

!    allocate(cache)
!    call solv%update(mol, cache)

!    ! --- numerical derivatives along the STORED 6 directions
!    ! For diagonals, apply trace-preserving bumps so directions match solver's traceless diagonal basis.
!    allocate(vnum(6, nat), source=0.0_wp)

!    do ii = 1, nat
!       do l = 1, 6
!          ! +step
!          wfn%qat  = reshape(qat,  [nat, 1])
!          wfn%dpat = reshape(dpat, [3,   nat, 1])
!          wfn%qpat = reshape(qpat, [6,   nat, 1])
!          call bump_traceless_stored(wfn%qpat(:,ii,1), l, +step)
!          call solv%get_energy(mol, cache, wfn, er)

!          ! -step
!          wfn%qat  = reshape(qat,  [nat, 1])
!          wfn%dpat = reshape(dpat, [3,   nat, 1])
!          wfn%qpat = reshape(qpat, [6,   nat, 1])
!          call bump_traceless_stored(wfn%qpat(:,ii,1), l, -step)
!          call solv%get_energy(mol, cache, wfn, el)

!          vnum(l, ii) = 0.5_wp*(sum(er) - sum(el))/step
!       end do
!    end do

!    ! --- analytic gradient from solver, same storage (6, nat)
!    pot%vat(:, :)    = 0.0_wp
!    pot%vqp(:, :, :) = 0.0_wp
!    call solv%get_potential(mol, cache, wfn, pot)

!    allocate(vana(6, nat))
!    vana(:,:) = pot%vqp(:,:,1)

!    ! --- auto-detect off-diagonal ×2 packing and correct if needed
!    offdiag_times_two = detect_offdiag_times_two(vnum, vana, nat)
!    if (offdiag_times_two) then
!       vana(ixy,:) = 0.5_wp*vana(ixy,:)
!       vana(ixz,:) = 0.5_wp*vana(ixz,:)
!       vana(iyz,:) = 0.5_wp*vana(iyz,:)
!       write(*,'(a)') "note: detected off-diagonal ×2 packing; halving analytic off-diagonals for comparison."
!    end if

!    ! --- compare directly (component-wise, same (6,nat) storage)
!    call check_match_or_fail(error, vnum, vana, nat, thr_abs, thr_rel)

! contains

!    ! Apply a bump to the l-th stored component while keeping the tensor traceless
!    ! for the diagonal triplet. Off-diagonals are already traceless and bumped directly.
!    pure subroutine bump_traceless_stored(q6, lidx, h)
!       real(wp), intent(inout) :: q6(6)     ! (xx,xy,yy,xz,yz,zz)
!       integer,  intent(in)    :: lidx
!       real(wp), intent(in)    :: h
!       integer, parameter :: ixx=1, ixy=2, iyy=3, ixz=4, iyz=5, izz=6
!       select case (lidx)
!       case (ixx)
!          q6(ixx) = q6(ixx) + h
!          q6(iyy) = q6(iyy) - 0.5_wp*h
!          q6(izz) = q6(izz) - 0.5_wp*h
!       case (iyy)
!          q6(iyy) = q6(iyy) + h
!          q6(ixx) = q6(ixx) - 0.5_wp*h
!          q6(izz) = q6(izz) - 0.5_wp*h
!       case (izz)
!          q6(izz) = q6(izz) + h
!          q6(ixx) = q6(ixx) - 0.5_wp*h
!          q6(iyy) = q6(iyy) - 0.5_wp*h
!       case default
!          q6(lidx) = q6(lidx) + h   ! xy/xz/yz: symmetric off-diagonals; trace already zero
!       end select
!    end subroutine bump_traceless_stored

!    ! Heuristic: if off-diagonal analytic components look ~2× the FD,
!    ! assume packing uses 2*G_ij and request halving for comparison.
!    logical function detect_offdiag_times_two(vn, va, n)
!       real(wp), intent(in) :: vn(6,n), va(6,n)
!       integer,  intent(in) :: n
!       real(wp) :: ratios(3*n), tmp
!       integer :: k, i, j
!       k = 0
!       do i = 1, n
!          if (abs(va(ixy,i)) > 0.0_wp .and. abs(vn(ixy,i)) > 0.0_wp) then
!             k = k + 1; ratios(k) = abs(vn(ixy,i)/va(ixy,i))
!          end if
!          if (abs(va(ixz,i)) > 0.0_wp .and. abs(vn(ixz,i)) > 0.0_wp) then
!             k = k + 1; ratios(k) = abs(vn(ixz,i)/va(ixz,i))
!          end if
!          if (abs(va(iyz,i)) > 0.0_wp .and. abs(vn(iyz,i)) > 0.0_wp) then
!             k = k + 1; ratios(k) = abs(vn(iyz,i)/va(iyz,i))
!          end if
!       end do
!       if (k < 3) then
!          detect_offdiag_times_two = .false.
!          return
!       end if
!       ! crude median-of-k
!       do i = 1, k-1
!          do j = i+1, k
!             if (ratios(j) < ratios(i)) then
!                tmp = ratios(i); ratios(i) = ratios(j); ratios(j) = tmp
!             end if
!          end do
!       end do
!       tmp = ratios((k+1)/2)
!       detect_offdiag_times_two = (tmp > 0.45_wp .and. tmp < 0.55_wp)  ! near 0.5
!    end function detect_offdiag_times_two

!    subroutine check_match_or_fail(error, vnum6n, vana6n, n, rabs, rrel)
!       use mctc_env, only : wp, error_type
!       type(error_type), allocatable, intent(out) :: error
!       integer, intent(in) :: n
!       real(wp), intent(in) :: vnum6n(6,n), vana6n(6,n)
!       real(wp), intent(in) :: rabs, rrel
!       integer :: ii, l
!       real(wp) :: diff, denom
!       logical :: failed
!       failed = .false.

!       do ii = 1, n
!          do l = 1, 6
!             diff  = abs(vnum6n(l,ii) - vana6n(l,ii))
!             denom = max(1.0_wp, abs(vana6n(l,ii)))
!             if (diff > max(rabs, rrel*denom)) then
!                failed = .true.
!                exit
!             end if
!          end do
!          if (failed) exit
!       end do

!       if (failed) then
!          call test_failed(error, "Quadrupole potential mismatch (direct packed, traceless diagonal directions).")
!          do ii = 1, n
!             write(*,'(a,i4)') 'atom', ii
!             write(*,'(a,6es20.13)') '  vnum: ', vnum6n(:,ii)
!             write(*,'(a,6es20.13)') '  vana: ', vana6n(:,ii)
!             write(*,'(a,6es20.13)') '  diff: ', vnum6n(:,ii) - vana6n(:,ii)
!          end do
!       end if
!    end subroutine check_match_or_fail

! end subroutine test_qp




subroutine test_qp(error, model, mol, qat, dpat, qpat, kappa)
   use mctc_env,                 only : wp, error_type
   use mctc_io,                  only : structure_type
   use tblite_solvation_ddx,     only : ddx_solvation, ddx_input
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_scf_potential,     only : potential_type
   use tblite_container_cache,   only : container_cache
   implicit none

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Solvation model (COSMO=11, CPCM=12, PCM=2, LPB=3)
   integer, intent(in) :: model
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Atomic partial charges (kept constant)
   real(wp), intent(in) :: qat(:)
   !> Atomic dipoles (kept constant)
   real(wp), intent(in) :: dpat(:, :)
   !> Atomic quadrupoles in lower-tri order (xx,xy,yy,xz,yz,zz), traceless
   real(wp), intent(in) :: qpat(:, :)
   !> Debye–Hückel screening parameter (only used in LPB)
   real(wp), optional, intent(in) :: kappa

   type(ddx_solvation)                :: solv
   type(wavefunction_type)            :: wfn
   type(potential_type)               :: pot
   type(container_cache), allocatable :: cache

   integer,  parameter :: nqp = 6, nang = 302
   integer,  parameter :: ixx=1, ixy=2, iyy=3, ixz=4, iyz=5, izz=6
   real(wp), parameter :: feps = 80.0_wp
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp), parameter :: thr_abs = 1.0e+3_wp*sqrt(epsilon(1.0_wp))
   real(wp), parameter :: thr_rel = 1.0e-7_wp

   integer :: ii, l, nat
   real(wp) :: er(mol%nat), el(mol%nat)
   real(wp), allocatable :: vnum(:,:), vdir(:,:)    ! numerical FD and analytic directional

   ! --- sizes & sanity ---
   nat = mol%nat
   if (size(qpat,1) /= nqp) then
      call test_failed(error, "test_qp expects qpat with 6 components (xx,xy,yy,xz,yz,zz).")
      return
   end if
   if (size(qpat,2) /= nat) then
      call test_failed(error, "test_qp size mismatch: size(qpat,2) must equal mol%nat.")
      return
   end if

   ! --- baseline wfn ---
   wfn%qat  = reshape(qat,  [size(qat),         1])
   wfn%dpat = reshape(dpat, [3, size(dpat, 2),  1])
   wfn%qpat = reshape(qpat, [nqp, size(qpat,2), 1])

   allocate(pot%vat(size(qat,1), 1), source=0.0_wp)
   allocate(pot%vqp(nqp, size(qpat,2), 1), source=0.0_wp)

   if (present(kappa)) then
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, kappa=kappa))
   else
      solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang))
   end if

   allocate(cache)
   call solv%update(mol, cache)

   ! --- numerical FD constrained to traceless manifold (option 1) ---
   allocate(vnum(nqp, nat), source=0.0_wp)

   do ii = 1, nat
      do l = 1, nqp
         er = 0.0_wp; el = 0.0_wp

         ! +step
         wfn%qat  = reshape(qat,  [size(qat), 1])
         wfn%dpat = reshape(dpat, [3, size(dpat,2), 1])
         wfn%qpat = reshape(qpat, [nqp, size(qpat,2), 1])
         call apply_traceless_bump(wfn%qpat(:,ii,1), l, +step)
         call enforce_traceless(wfn%qpat(:,ii,1))
         call solv%get_energy(mol, cache, wfn, er)

         ! -step
         wfn%qat  = reshape(qat,  [size(qat), 1])
         wfn%dpat = reshape(dpat, [3, size(dpat,2), 1])
         wfn%qpat = reshape(qpat, [nqp, size(qpat,2), 1])
         call apply_traceless_bump(wfn%qpat(:,ii,1), l, -step)
         call enforce_traceless(wfn%qpat(:,ii,1))
         call solv%get_energy(mol, cache, wfn, el)

         vnum(l, ii) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   ! --- analytic dual at baseline (convert to STORED 6-vector dual) ---
   pot%vat(:, :)    = 0.0_wp
   pot%vqp(:, :, :) = 0.0_wp
   call solv%get_potential(mol, cache, wfn, pot)  ! raw symmetric dual

   allocate(vdir(nqp, nat), source=0.0_wp)

   do ii = 1, nat
      ! Start from raw symmetric dual g_raw = pot%vqp(:,ii,1)
      ! Map to STORED dual: off-diagonals ×2 so that  δE = g_stored · δq_stored
      ! (diagonals unchanged)
      vdir(ixx,ii) = 0.0_wp
      vdir(ixy,ii) = 0.0_wp
      vdir(iyy,ii) = 0.0_wp
      vdir(ixz,ii) = 0.0_wp
      vdir(iyz,ii) = 0.0_wp
      vdir(izz,ii) = 0.0_wp

      call directional_dual_from_raw(pot%vqp(:,ii,1), vdir(:,ii))
   end do

   ! Now vdir(:,ii) holds the analytic **directional** values:
   !   for l=xx: g·(1,0,-1/2,0,0,-1/2), etc.; for off-diags: plain component.
   ! Compare vnum (FD) vs vdir (analytic) directly.

call check_match_or_fail(error, vnum, vdir, nat, thr_abs, thr_rel)
contains

   pure subroutine apply_traceless_bump(q6, lidx, h)
      ! Apply a perturbation of size h to component lidx while preserving trace.
      real(wp), intent(inout) :: q6(6)     ! (xx,xy,yy,xz,yz,zz)
      integer, intent(in)     :: lidx
      real(wp), intent(in)    :: h
      select case (lidx)
      case (ixx)
         q6(ixx) = q6(ixx) + h
         q6(iyy) = q6(iyy) - 0.5_wp*h
         q6(izz) = q6(izz) - 0.5_wp*h
      case (iyy)
         q6(iyy) = q6(iyy) + h
         q6(ixx) = q6(ixx) - 0.5_wp*h
         q6(izz) = q6(izz) - 0.5_wp*h
      case (izz)
         q6(izz) = q6(izz) + h
         q6(ixx) = q6(ixx) - 0.5_wp*h
         q6(iyy) = q6(iyy) - 0.5_wp*h
      case default
         q6(lidx) = q6(lidx) + h   ! xy/xz/yz: already trace-free
      end select
   end subroutine apply_traceless_bump

   pure subroutine enforce_traceless(q6)
      real(wp), intent(inout) :: q6(6)
      real(wp) :: tthird
      tthird = (q6(ixx) + q6(iyy) + q6(izz)) / 3.0_wp
      q6(ixx) = q6(ixx) - tthird
      q6(iyy) = q6(iyy) - tthird
      q6(izz) = q6(izz) - tthird
   end subroutine enforce_traceless

  pure subroutine directional_dual_from_raw(graw6, gdir6)
   ! Use RAW symmetric dual directly for directional comparison.
   ! No off-diagonal scaling. Diagonal directions are the constrained
   ! traceless combos matching the FD bumps.
   real(wp), intent(in)  :: graw6(6)   ! (xx,xy,yy,xz,yz,zz) raw dual from solver
   real(wp), intent(out) :: gdir6(6)   ! directional values in "stored" index order
   integer, parameter :: ixx=1, ixy=2, iyy=3, ixz=4, iyz=5, izz=6

   ! Optionally remove isotropic part; it cancels in the combos anyway,
   ! but doing it keeps things numerically clean.
   real(wp) :: gxx, gyy, gzz, tthird
   gxx = graw6(ixx); gyy = graw6(iyy); gzz = graw6(izz)
   tthird = (gxx + gyy + gzz)/3.0_wp
   gxx = gxx - tthird
   gyy = gyy - tthird
   gzz = gzz - tthird

   ! Directional dots for the constrained diagonal directions
   gdir6(ixx) = gxx - 0.5_wp*(gyy + gzz)
   gdir6(iyy) = gyy - 0.5_wp*(gxx + gzz)
   gdir6(izz) = gzz - 0.5_wp*(gxx + gyy)

   ! Off-diagonals: use raw dual directly (no ×2)
   gdir6(ixy) = graw6(ixy)
   gdir6(ixz) = graw6(ixz)
   gdir6(iyz) = graw6(iyz)
end subroutine directional_dual_from_raw


   subroutine check_match_or_fail(error, vnum6n, vdir6n, nat, rabs, rrel)
   use mctc_env, only : wp, error_type
   implicit none
   type(error_type), allocatable, intent(out) :: error
   integer, intent(in) :: nat
   real(wp), intent(in) :: vnum6n(6, nat), vdir6n(6, nat)
   real(wp), intent(in) :: rabs, rrel

   integer :: ii, l
   real(wp) :: diff, denom
   logical :: failed
   failed = .false.

   do ii = 1, nat
      do l = 1, 6
         diff  = abs(vnum6n(l,ii) - vdir6n(l,ii))
         denom = max(1.0_wp, abs(vdir6n(l,ii)))
         if (diff > max(rabs, rrel*denom)) then
            failed = .true.
            exit
         end if
      end do
      if (failed) exit
   end do

   if (failed) then
      call test_failed(error, "Quadrupole potential mismatch (FD vs analytic along SAME constrained directions).")
      print '(a)', 'numerical vnum (directional, stored)'
      do ii = 1, nat
         print '(3es20.13)', vnum6n(:,ii)
      end do
      print '(a)', 'analytic vdir (directional, stored)'
      do ii = 1, nat
         print '(3es20.13)', vdir6n(:,ii)
      end do
      print '(a)', 'diff (numerical - analytic)'
      do ii = 1, nat
         print '(3es20.13)', vnum6n(:,ii) - vdir6n(:,ii)
      end do
   end if
end subroutine check_match_or_fail


end subroutine test_qp


! subroutine test_qp(error, model, mol, qat, dpat, qpat, kappa)
!    use mctc_env,              only : wp, error_type
!    use mctc_io,               only : structure_type
!    use tblite_solvation_ddx,  only : ddx_solvation, ddx_input
!    use tblite_wavefunction_type, only : wavefunction_type
!    use tblite_scf_potential,  only : potential_type
!    use tblite_container_cache,only : container_cache
!    implicit none

!    !> Error handling
!    type(error_type), allocatable, intent(out) :: error
!    !> Solvation model (COSMO=11, CPCM=12, PCM=2, LPB=3)
!    integer, intent(in) :: model
!    !> Molecular structure data
!    type(structure_type), intent(in) :: mol
!    !> Atomic partial charges (kept constant)
!    real(wp), intent(in) :: qat(:)
!    !> Atomic dipoles (kept constant)
!    real(wp), intent(in) :: dpat(:, :)
!    !> Atomic quadrupoles in Cartesian lower-triangle order (xx,xy,yy,xz,yz,zz), traceless
!    real(wp), intent(in) :: qpat(:, :)
!    !> Debye–Hückel screening parameter (only used in LPB)
!    real(wp), optional, intent(in) :: kappa

!    type(ddx_solvation)                  :: solv
!    type(wavefunction_type)              :: wfn
!    type(potential_type)                 :: pot
!    type(container_cache), allocatable   :: cache

!    integer,  parameter :: nqp = 6, nang = 302
!    real(wp), parameter :: feps = 80.0_wp
!    real(wp), parameter :: step = 1.0e-4_wp
!    real(wp), parameter :: thr_abs = 1.0e+3_wp*sqrt(epsilon(1.0_wp))
!    real(wp), parameter :: thr_rel = 1.0e-7_wp

!    integer :: ii, l, nat
!    real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat)

!    ! Work arrays
!    real(wp), allocatable :: vnum(:,:), vnum_raw(:,:), vnum_proj(:,:)
!    real(wp), allocatable :: vanal_cart(:,:), vanal_proj(:,:), diff(:,:)

!    ! --- sizes & sanity ---
!    nat = mol%nat
!    if (size(qpat,1) /= nqp) then
!       call test_failed(error, "test_qp expects qpat with 6 Cartesian components (xx,xy,yy,xz,yz,zz).")
!       return
!    end if
!    if (size(qpat,2) /= nat) then
!       call test_failed(error, "test_qp size mismatch: size(qpat,2) must equal mol%nat.")
!       return
!    end if

!    ! --- baseline wfn ---
!    wfn%qat  = reshape(qat,  [size(qat),         1])
!    wfn%dpat = reshape(dpat, [3, size(dpat, 2),  1])
!    wfn%qpat = reshape(qpat, [nqp, size(qpat,2), 1])

!    allocate(pot%vat(size(qat,1), 1), source=0.0_wp)
!    allocate(pot%vqp(nqp, size(qpat,2), 1), source=0.0_wp)

!    if (present(kappa)) then
!       solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang, kappa=kappa))
!    else
!       solv = ddx_solvation(mol, ddx_input(feps, model, nang=nang))
!    end if

!    allocate(cache)
!    call solv%update(mol, cache)

!    ! --- numerical FD wrt stored 6D components (xx,xy,yy,xz,yz,zz) ---
!    allocate(vnum(nqp, nat), source=0.0_wp)
!    do ii = 1, nat
!       do l = 1, nqp
!          er = 0.0_wp; el = 0.0_wp

!          ! +step
!          wfn%qat  = reshape(qat,  [size(qat), 1])
!          wfn%dpat = reshape(dpat, [3, size(dpat,2), 1])
!          wfn%qpat = reshape(qpat, [nqp, size(qpat,2), 1])
!          wfn%qpat(l, ii, 1) = wfn%qpat(l, ii, 1) + step
!          call solv%get_energy(mol, cache, wfn, er)

!          ! -step
!          wfn%qat  = reshape(qat,  [size(qat), 1])
!          wfn%dpat = reshape(dpat, [3, size(dpat,2), 1])
!          wfn%qpat = reshape(qpat, [nqp, size(qpat,2), 1])
!          wfn%qpat(l, ii, 1) = wfn%qpat(l, ii, 1) - step
!          call solv%get_energy(mol, cache, wfn, el)

!          vnum(l, ii) = 0.5_wp*(sum(er) - sum(el))/step
!       end do
!    end do

!    ! --- analytic potential at baseline (raw symmetric dual from solver) ---
!    energy = 0.0_wp
!    wfn%qat  = reshape(qat,  [size(qat), 1])
!    wfn%dpat = reshape(dpat, [3, size(dpat,2), 1])
!    wfn%qpat = reshape(qpat, [nqp, size(qpat,2), 1])
!    pot%vat(:, :)    = 0.0_wp
!    pot%vqp(:, :, :) = 0.0_wp
!    call solv%get_potential(mol, cache, wfn, pot)   ! pot%vqp(:,ii,1) is CONJUGATE TO RAW Q_ij (no ×2)

!    ! --- compare apples-to-apples ---
!    ! (A) Convert FD (stored) -> raw symmetric dual by halving off-diagonals once.
!    allocate(vnum_raw(nqp, nat)); vnum_raw = vnum
!    vnum_raw(2,:) = 0.5_wp * vnum_raw(2,:)   ! xy
!    vnum_raw(4,:) = 0.5_wp * vnum_raw(4,:)   ! xz
!    vnum_raw(5,:) = 0.5_wp * vnum_raw(5,:)   ! yz

!    ! (B) Project BOTH sides to the traceless dual space on the diagonals
!    allocate(vnum_proj(nqp, nat));  vnum_proj  = 0.0_wp
!    allocate(vanal_cart(nqp, nat)); vanal_cart = 0.0_wp
!    allocate(vanal_proj(nqp, nat)); vanal_proj = 0.0_wp
!    do ii = 1, nat
!       vanal_cart(:, ii) = pot%vqp(:, ii, 1)                 ! raw dual from solver
!       vnum_proj(:,  ii) = proj_traceless6( vnum_raw(:, ii) )
!       vanal_proj(:, ii) = proj_traceless6( vanal_cart(:, ii) )
!    end do

!    ! (C) Check thresholds
!    allocate(diff(nqp, nat)); diff = vnum_proj - vanal_proj
!    do ii = 1, nat
!       do l = 1, nqp
!          if (.not. pass_thresh(abs(diff(l,ii)), abs(vanal_proj(l,ii)), thr_abs, thr_rel)) then
!             call test_failed(error, "Quadrupole potential mismatch (compare in RAW+traceless dual; FD off-diags halved).")
!             print '(a)', 'numerical (FD→raw) vnum_proj'
!             print '(3es20.13)', vnum_proj
!             print '(a)', 'analytic (raw) projected vanal_proj'
!             print '(3es20.13)', vanal_proj
!             print '(a)', 'diff (numerical - analytic)'
!             print '(3es20.13)', diff
!             return
!          end if
!       end do
!    end do

! contains

!    pure function proj_traceless6(v) result(u)
!       ! Project a 6-vector (xx,xy,yy,xz,yz,zz) to the traceless dual on diagonals:
!       ! diagonals: u = (I - (1/3)11^T) * v_diag; off-diagonals pass through.
!       real(wp), intent(in) :: v(6)
!       real(wp) :: u(6)
!       real(wp) :: vxx, vyy, vzz, tthird
!       vxx = v(1); vyy = v(3); vzz = v(6)
!       tthird = (vxx + vyy + vzz) / 3.0_wp
!       u(1) = vxx - tthird
!       u(3) = vyy - tthird
!       u(6) = vzz - tthird
!       u(2) = v(2)
!       u(4) = v(4)
!       u(5) = v(5)
!    end function proj_traceless6

!    pure logical function pass_thresh(d, a, rabs, rrel) result(ok)
!       ! Accept if |diff| <= max(absolute_thresh, relative_thresh * max(1, |analytic|))
!       real(wp), intent(in) :: d, a, rabs, rrel
!       real(wp) :: denom
!       denom = max(1.0_wp, abs(a))
!       ok = (d <= max(rabs, rrel*denom))
!    end function pass_thresh

! end subroutine test_qp




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
   call test_e(error, ddx_solvation_model%cosmo, mol, qat, -3.4697720884118800E-2_wp)

end subroutine test_e_cosmo_m01

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
   call test_e(error, ddx_solvation_model%pcm, mol, qat, -3.3624259293951506E-2_wp)

end subroutine test_e_pcm_m01

subroutine test_e_lpb_m01(error)

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
   call test_e(error, ddx_solvation_model%lpb, mol, qat, -3.1747235888764228E-2_wp, kappa=0.5_wp)

end subroutine test_e_lpb_m01


subroutine test_g_num_cosmo_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.38394711236234E-2_wp,-1.68354976558608E-1_wp,-3.47642833746823E-1_wp,&
      &-7.05489267186003E-1_wp, 7.73548301641266E-1_wp, 2.30207581365386E-1_wp,&
      & 1.02748501676354E-1_wp, 9.47818107467040E-2_wp, 2.44260351729187E-2_wp,&
      & 2.34984927037408E-1_wp,-3.17839896393030E-1_wp, 6.67112994818879E-1_wp,&
      &-4.78119977010488E-1_wp, 6.57536027459275E-2_wp, 1.08259054549882E-1_wp,&
      &-3.58215329983396E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_g_num(error, ddx_solvation_model%cosmo, mol, qat)

end subroutine test_g_num_cosmo_m02

subroutine test_g_cosmo_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.38394711236234E-2_wp,-1.68354976558608E-1_wp,-3.47642833746823E-1_wp,&
      &-7.05489267186003E-1_wp, 7.73548301641266E-1_wp, 2.30207581365386E-1_wp,&
      & 1.02748501676354E-1_wp, 9.47818107467040E-2_wp, 2.44260351729187E-2_wp,&
      & 2.34984927037408E-1_wp,-3.17839896393030E-1_wp, 6.67112994818879E-1_wp,&
      &-4.78119977010488E-1_wp, 6.57536027459275E-2_wp, 1.08259054549882E-1_wp,&
      &-3.58215329983396E-1_wp]

   real(wp), parameter :: ref(3,16) = reshape([&
      & -8.72885014e-04_wp,  1.29722126e-03_wp,  4.82135743e-04_wp,  2.86976918e-03_wp, &
      &  4.91109654e-05_wp, -4.29260562e-04_wp, -5.67289804e-03_wp,  8.14714673e-04_wp, &
      & -4.12468023e-04_wp, -1.64554349e-04_wp,  1.83163634e-03_wp,  8.14645705e-04_wp, &
      & -2.92212009e-03_wp, -1.29547615e-03_wp,  1.77475155e-03_wp,  1.83559208e-03_wp, &
      & -2.82980413e-04_wp,  1.71557120e-03_wp,  2.73353837e-03_wp, -2.05688333e-03_wp, &
      &  6.38494126e-03_wp,  5.78537370e-04_wp, -2.37129787e-03_wp,  2.90753892e-04_wp, &
      &  5.08752936e-04_wp,  2.15648519e-04_wp, -1.25921100e-03_wp,  1.01160620e-03_wp, &
      & -4.55980356e-03_wp,  1.90794668e-03_wp, -5.63872868e-03_wp,  8.21396325e-04_wp, &
      &  5.02731375e-03_wp,  4.01930567e-03_wp, -1.74485028e-03_wp, -1.42778390e-02_wp, &
      & -6.01680881e-03_wp,  2.69775116e-04_wp,  1.42002404e-02_wp, -7.83354901e-04_wp, &
      & -1.47082304e-03_wp,  5.63906422e-04_wp,  6.02706612e-03_wp,  5.75784078e-04_wp, &
      & -7.23463064e-04_wp, -1.65844186e-04_wp, -5.66640254e-03_wp,  1.66094445e-04_wp &
      ], shape=[3,16])

   call get_structure(mol, "MB16-43", "02")
   call test_g(error, ddx_solvation_model%cosmo, mol, qat, ref)

end subroutine test_g_cosmo_m02

subroutine test_g_pcm_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.38394711236234E-2_wp,-1.68354976558608E-1_wp,-3.47642833746823E-1_wp,&
      &-7.05489267186003E-1_wp, 7.73548301641266E-1_wp, 2.30207581365386E-1_wp,&
      & 1.02748501676354E-1_wp, 9.47818107467040E-2_wp, 2.44260351729187E-2_wp,&
      & 2.34984927037408E-1_wp,-3.17839896393030E-1_wp, 6.67112994818879E-1_wp,&
      &-4.78119977010488E-1_wp, 6.57536027459275E-2_wp, 1.08259054549882E-1_wp,&
      &-3.58215329983396E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_g_num(error, ddx_solvation_model%pcm, mol, qat)

end subroutine test_g_pcm_m02

subroutine test_g_lpb_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.38394711236234E-2_wp,-1.68354976558608E-1_wp,-3.47642833746823E-1_wp,&
      &-7.05489267186003E-1_wp, 7.73548301641266E-1_wp, 2.30207581365386E-1_wp,&
      & 1.02748501676354E-1_wp, 9.47818107467040E-2_wp, 2.44260351729187E-2_wp,&
      & 2.34984927037408E-1_wp,-3.17839896393030E-1_wp, 6.67112994818879E-1_wp,&
      &-4.78119977010488E-1_wp, 6.57536027459275E-2_wp, 1.08259054549882E-1_wp,&
      &-3.58215329983396E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_g_num(error, ddx_solvation_model%lpb, mol, qat, kappa=0.5_wp)

end subroutine test_g_lpb_m02


subroutine test_p_cosmo_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-1.77788256288236E-1_wp,-8.22943267808161E-1_wp, 4.04578389873281E-2_wp,&
      & 5.79710531992282E-1_wp, 6.99601887637659E-1_wp, 6.84309612639107E-2_wp,&
      &-3.42971414989811E-1_wp, 4.64954031865410E-2_wp, 6.77012204116428E-2_wp,&
      & 8.49931225363225E-2_wp,-5.22285304699699E-1_wp,-2.92515001764712E-1_wp,&
      &-3.98375452377043E-1_wp, 2.09769668297792E-1_wp, 7.23140464830357E-1_wp,&
      & 3.65775987838250E-2_wp]

   real(wp), parameter :: dpat(3,16) = reshape([&
     3.08663678818639E-02_wp, -1.23672242032681E-02_wp, -7.80809552815601E-02_wp,&
     1.11378962331252E-02_wp,  2.49592863297026E-02_wp, -3.43860910100973E-02_wp,&
    -3.14360758213388E-02_wp, -1.02773078762668E-01_wp,  1.35869342433299E-01_wp,&
    -6.59221036649827E-02_wp, -4.98746708517068E-03_wp,  8.18392699352346E-03_wp,&
     2.94110714733907E-01_wp,  4.32992951768027E-01_wp,  5.73566729874146E-02_wp,&
     5.57444981972002E-02_wp,  8.64117606916050E-02_wp, -2.90348073828545E-02_wp,&
     6.06660962349574E-02_wp,  1.67983241027530E-02_wp, -3.13513724081013E-02_wp,&
     4.01524751946939E-02_wp, -1.68321289994974E-01_wp, -2.83766398426090E-03_wp,&
    -1.28822537767878E-01_wp,  1.18980146094674E-01_wp, -1.87396729496736E-02_wp,&
     4.05166776408743E-02_wp, -6.16446805339326E-02_wp, -1.46993806131036E-01_wp,&
    -2.11957220453610E-02_wp, -2.17404002169060E-03_wp, -1.13706491905540E-02_wp,&
     1.05822848372919E-01_wp, -1.32384216104661E-02_wp, -1.10132480312151E-01_wp,&
     1.98128712421592E-01_wp,  1.90132721681466E-01_wp,  1.04924269271602E-01_wp,&
    -3.02008458727527E-02_wp, -8.72261745788219E-02_wp,  1.30187020977396E-01_wp,&
     1.21967269288880E-01_wp,  3.36027641726082E-02_wp,  1.75601803507984E-02_wp,&
     8.14663591946045E-02_wp,  6.23935351174657E-02_wp,  3.46979182796507E-02_wp],&
     [3,16])

  real(wp), parameter :: qpat(6,16) = reshape([ &
     0.61568297965979824_wp,  0.15216743018229489_wp, -6.8858163770286529E-002_wp, &
     0.20034415597194447_wp,  7.0998034740450036E-002_wp, -0.54682481588951060_wp, &
    -0.18035001616235569_wp, -1.5290927213447330E-002_wp,  5.6748939159116274E-002_wp, &
    -0.12469658457423533_wp, -0.12284775097921814_wp,  0.12360107700323897_wp, &
     0.12288957868143997_wp, -4.6154720858250807E-002_wp, -3.2959513350267464E-002_wp, &
     5.3131785701090901E-002_wp,  0.19734030046031048_wp, -8.9930065331172454E-002_wp, &
     0.43195863777460369_wp,  0.26975537056921134_wp,  0.14248632607348222_wp, &
    -6.6425231603841756E-003_wp,  9.9391879939014979E-003_wp, -0.57444496384808597_wp, &
    -2.1930144731313788_wp, -0.92392459741636790_wp,  1.6825024012991843_wp, &
    -1.3996979526542677_wp, -1.4954678946689095_wp,  0.51051207183219416_wp, &
    -0.86133556568282454_wp,  4.2053815129143968E-002_wp,  0.34472931268682705_wp, &
    -0.10178254541038746_wp,  4.3590080464292080E-002_wp,  0.51660625299599761_wp, &
     0.60970669320980697_wp,  0.30110153863653299_wp, -0.49090396185915242_wp, &
     0.15688347235376510_wp,  0.12033699914653465_wp, -0.11880273135065400_wp, &
     0.12095062175035461_wp,  8.1530674251783619E-002_wp, -0.24841814373094451_wp, &
    -6.2667267017851099E-004_wp,  1.0482410400213037E-002_wp,  0.12746752198058992_wp, &
    -9.0581099943795915E-002_wp,  0.20067815983035059_wp, -4.1512872600256524E-002_wp, &
    -4.8667836801105271E-002_wp,  4.0939086164223308E-002_wp,  0.13209397254405233_wp, &
     0.10265014793892824_wp,  3.1101633555667963E-002_wp,  8.5063335086427563E-002_wp, &
     8.4166467055277791E-002_wp, -0.11946027934068097_wp, -0.18771348302535568_wp, &
     3.9847012397439041E-002_wp,  6.5551551551080349E-003_wp, -2.3253670772282020E-002_wp, &
     2.3163091024144163E-002_wp,  2.0074269396591028E-003_wp, -1.6593341625158353E-002_wp, &
    -4.2477411473465974E-002_wp, -8.9887638410775864E-002_wp,  0.13195734956967245_wp, &
     1.5153820167078599E-002_wp,  8.8852027548639259E-002_wp, -8.9479938096208050E-002_wp, &
     0.10879019619465635_wp,  0.39714952843489526_wp, -2.5492394806908569E-004_wp, &
     0.24165082351889752_wp,  0.38441548100362644_wp, -0.10853527224658738_wp, &
     8.6674579935310289E-002_wp, -2.9705533015954605E-002_wp,  9.7008797832631045E-003_wp, &
     4.4393275271057221E-002_wp,  0.12962247713961847_wp, -9.6375459718573392E-002_wp, &
    -9.6698478856955261E-003_wp, -0.39475198227878894_wp, -0.18088530469418387_wp, &
     0.22479716036464653_wp,  0.40311368079845267_wp,  0.19055515257987959_wp, &
    -9.9359049825878815E-002_wp, -0.10573490419971270_wp,  2.6836974999905311E-002_wp, &
    -2.6928537127614872E-002_wp, -3.2408785761527094E-002_wp,  7.2522074825973559E-002_wp &
   ], [6,16])

   call get_structure(mol, "MB16-43", "03")
   ! call test_dp(error, ddx_solvation_model%cosmo, mol, qat, dpat)
   call test_qp_traceless(error, ddx_solvation_model%cosmo, mol, qat, dpat, qpat)

end subroutine test_p_cosmo_m03

! subroutine test_p_pcm_m03(error)

!    !> Error handling
!    type(error_type), allocatable, intent(out) :: error

!    type(structure_type) :: mol
!    real(wp), parameter :: qat(*) = [&
!       &-1.77788256288236E-1_wp,-8.22943267808161E-1_wp, 4.04578389873281E-2_wp,&
!       & 5.79710531992282E-1_wp, 6.99601887637659E-1_wp, 6.84309612639107E-2_wp,&
!       &-3.42971414989811E-1_wp, 4.64954031865410E-2_wp, 6.77012204116428E-2_wp,&
!       & 8.49931225363225E-2_wp,-5.22285304699699E-1_wp,-2.92515001764712E-1_wp,&
!       &-3.98375452377043E-1_wp, 2.09769668297792E-1_wp, 7.23140464830357E-1_wp,&
!       & 3.65775987838250E-2_wp]

!    call get_structure(mol, "MB16-43", "03")
!    call test_p(error, ddx_solvation_model%pcm, mol, qat)

! end subroutine test_p_pcm_m03

! subroutine test_p_lpb_m03(error)

!    !> Error handling
!    type(error_type), allocatable, intent(out) :: error

!    type(structure_type) :: mol
!    real(wp), parameter :: qat(*) = [&
!       &-1.77788256288236E-1_wp,-8.22943267808161E-1_wp, 4.04578389873281E-2_wp,&
!       & 5.79710531992282E-1_wp, 6.99601887637659E-1_wp, 6.84309612639107E-2_wp,&
!       &-3.42971414989811E-1_wp, 4.64954031865410E-2_wp, 6.77012204116428E-2_wp,&
!       & 8.49931225363225E-2_wp,-5.22285304699699E-1_wp,-2.92515001764712E-1_wp,&
!       &-3.98375452377043E-1_wp, 2.09769668297792E-1_wp, 7.23140464830357E-1_wp,&
!       & 3.65775987838250E-2_wp]

!    call get_structure(mol, "MB16-43", "03")
!    call test_p(error, ddx_solvation_model%lpb, mol, qat, kappa=0.5_wp)

! end subroutine test_p_lpb_m03

end module test_solvation_ddx
