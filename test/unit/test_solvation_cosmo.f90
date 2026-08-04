! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later

module test_solvation_cosmo
   use mctc_env, only : wp
   use mctc_env_testing, only : error_type, new_unittest, test_failed, unittest_type
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_container, only : container_cache
   use tblite_scf_potential, only : potential_type
   use tblite_solvation_cosmo, only : cosmo_input, cosmo_solvation, cosmo_solvation_model, &
      & get_gaussian_multipole_kernels, new_cosmo
   use tblite_wavefunction, only : wavefunction_type
   implicit none
   private

   public :: collect_solvation_cosmo

contains

subroutine collect_solvation_cosmo(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-cosmo-monopole", &
         test_e_cosmo_monopole_m01), &
      new_unittest("energy-cosmo-dipole", &
         test_e_cosmo_dipole_m03), &
      new_unittest("energy-cosmo-quadrupole", &
         test_e_cosmo_quadrupole_m03), &
      new_unittest("gaussian-multipole-finite-difference", &
         test_gaussian_multipole_finite_difference), &
      new_unittest("gaussian-multipole-point-limit", &
         test_gaussian_multipole_point_limit) &
      ]
end subroutine collect_solvation_cosmo


subroutine test_e(error, model, mol, monopoles, dipoles, quadrupoles, ref, qat, dpat, qpat)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Solvation model (COSMO=1, CPCM=2)
   integer, intent(in) :: model

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Whether to include monopoles in the solute representation
   logical, intent(in) :: monopoles

   !> Whether to include dipoles in the solute representation
   logical, intent(in) :: dipoles

   !> Whether to include quadrupoles in the solute representation
   logical, intent(in) :: quadrupoles

   !> Reference energy
   real(wp), intent(in) :: ref

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Atomic dipole moments
   real(wp), intent(in), optional :: dpat(:,:)

   !> Atomic quadrupole moments
   real(wp), intent(in), optional :: qpat(:,:)
   

   type(cosmo_solvation) :: solv
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: eps = 80.0_wp
   integer, parameter :: nang = 302
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp) :: energy(mol%nat)

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1), source=0.0_wp)
   if (present(dpat)) then
      wfn%dpat = reshape(dpat, [3, mol%nat, 1])
      allocate(pot%vdp(3, mol%nat, 1), source=0.0_wp)
   end if
   if (present(qpat)) then
      wfn%qpat = reshape(qpat, [6, mol%nat, 1])
      allocate(pot%vqp(6, mol%nat, 1))
   end if
   energy = 0.0_wp

   call new_cosmo(solv, mol, cosmo_input(model=model, dielectric_const=eps, nang=nang, &
      & monopoles=monopoles, dipoles=dipoles, quadrupoles=quadrupoles), error)
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

!> Test COSMO solvation energy with point monopole representation of solute electron density 
subroutine test_e_cosmo_monopole_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   logical :: monopoles = .true.
   logical :: dipoles = .false.
   logical :: quadrupoles = .false.

   real(wp), parameter :: qat(*) = [&
      & 7.73347900345264E-1_wp, 1.07626888948184E-1_wp,-3.66999593831010E-1_wp,&
      & 4.92833325937897E-2_wp,-1.83332156197733E-1_wp, 2.33302086605469E-1_wp,&
      & 6.61837152062315E-2_wp,-5.43944165050002E-1_wp,-2.70264356583716E-1_wp,&
      & 2.66618968841682E-1_wp, 2.62725033202480E-1_wp,-7.15315510172571E-2_wp,&
      &-3.73300777019193E-1_wp, 3.84585237785621E-2_wp,-5.05851088366940E-1_wp,&
      & 5.17677238544189E-1_wp]

   call get_structure(mol, "MB16-43", "01")

   call test_e(error, cosmo_solvation_model%cosmo, mol, monopoles, dipoles, &
      & quadrupoles, -1.9661425306307E-02_wp, qat)

end subroutine test_e_cosmo_monopole_m01

!> Test COSMO solvation energy with point monopole+dipole representation of solute electron density
subroutine test_e_cosmo_dipole_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   logical :: monopoles = .true.
   logical :: dipoles = .true.
   logical :: quadrupoles = .false.

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

   call test_e(error, cosmo_solvation_model%cosmo, mol, monopoles, dipoles, &
      & quadrupoles, -3.7657586483210E-02_wp, qat, dpat=dpat)

end subroutine test_e_cosmo_dipole_m03

!> Test COSMO solvation energy with point monopole+dipole+quadrupole representation of solute electron density
subroutine test_e_cosmo_quadrupole_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   logical :: monopoles = .true.
   logical :: dipoles = .true.
   logical :: quadrupoles = .true.

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

   call test_e(error, cosmo_solvation_model%cosmo, mol, monopoles, dipoles, &
      & quadrupoles, -3.7794562071892E-02_wp, qat, dpat=dpat, qpat=qpat)

end subroutine test_e_cosmo_quadrupole_m03

subroutine test_gaussian_multipole_finite_difference(error)
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: xis(3) = [0.45_wp, 1.20_wp, 2.75_wp]
   real(wp), parameter :: vecs(3, 3) = reshape([ &
      1.30_wp, -0.70_wp, 2.10_wp, &
     -1.80_wp,  1.10_wp, 0.60_wp, &
      0.55_wp,  1.65_wp, -1.25_wp], [3, 3])
   real(wp), parameter :: qmat(3, 3) = reshape([ &
       0.70_wp,  0.20_wp, -0.10_wp, &
       0.20_wp, -0.40_wp,  0.30_wp, &
      -0.10_wp,  0.30_wp, -0.30_wp], [3, 3])
   real(wp), parameter :: qvec(6) = [ &
      qmat(1, 1), qmat(1, 2), qmat(2, 2), &
      qmat(1, 3), qmat(2, 3), qmat(3, 3)]
   real(wp), parameter :: h = 5.0e-4_wp
   real(wp), parameter :: dipole_tol = 2.0e-9_wp
   real(wp), parameter :: quadrupole_tol = 2.0e-7_wp

   real(wp) :: center(3), shifted(3), c0, c1(3), c2(6)
   real(wp) :: dnum(3), hnum(3, 3), f0, qnum, qimpl
   integer :: icase, ia, ib

   do icase = 1, size(xis)
      ! Put the surface Gaussian at the origin, so s_i-R_A=vecs(:,icase).
      center = -vecs(:, icase)
      call get_gaussian_multipole_kernels(vecs(:, icase), xis(icase), &
         c0, c1, c2)
      f0 = gaussian_monopole(center, xis(icase))

      do ia = 1, 3
         shifted = center
         shifted(ia) = center(ia) - 2.0_wp*h
         dnum(ia) = gaussian_monopole(shifted, xis(icase))
         shifted(ia) = center(ia) - h
         dnum(ia) = dnum(ia) - 8.0_wp*gaussian_monopole(shifted, xis(icase))
         shifted(ia) = center(ia) + h
         dnum(ia) = dnum(ia) + 8.0_wp*gaussian_monopole(shifted, xis(icase))
         shifted(ia) = center(ia) + 2.0_wp*h
         dnum(ia) = (dnum(ia) - gaussian_monopole(shifted, xis(icase))) &
            & /(12.0_wp*h)

         if (abs(c1(ia) - dnum(ia)) > dipole_tol) then
            call test_failed(error, &
               "Gaussian dipole kernel disagrees with finite differences")
            return
         end if

         shifted = center
         shifted(ia) = center(ia) - 2.0_wp*h
         hnum(ia, ia) = -gaussian_monopole(shifted, xis(icase))
         shifted(ia) = center(ia) - h
         hnum(ia, ia) = hnum(ia, ia) &
            & + 16.0_wp*gaussian_monopole(shifted, xis(icase))
         shifted(ia) = center(ia) + h
         hnum(ia, ia) = hnum(ia, ia) &
            & + 16.0_wp*gaussian_monopole(shifted, xis(icase))
         shifted(ia) = center(ia) + 2.0_wp*h
         hnum(ia, ia) = (hnum(ia, ia) &
            & - gaussian_monopole(shifted, xis(icase)) - 30.0_wp*f0) &
            & /(12.0_wp*h*h)
      end do

      do ia = 1, 3
         do ib = ia + 1, 3
            shifted = center
            shifted(ia) = center(ia) + h
            shifted(ib) = center(ib) + h
            hnum(ia, ib) = gaussian_monopole(shifted, xis(icase))
            shifted(ib) = center(ib) - h
            hnum(ia, ib) = hnum(ia, ib) &
               & - gaussian_monopole(shifted, xis(icase))
            shifted(ia) = center(ia) - h
            hnum(ia, ib) = hnum(ia, ib) &
               & + gaussian_monopole(shifted, xis(icase))
            shifted(ib) = center(ib) + h
            hnum(ia, ib) = hnum(ia, ib) &
               & - gaussian_monopole(shifted, xis(icase))
            hnum(ia, ib) = hnum(ia, ib)/(4.0_wp*h*h)
            hnum(ib, ia) = hnum(ia, ib)
         end do
      end do

      qnum = sum(qmat*hnum)/3.0_wp
      qimpl = dot_product(qvec, c2)
      if (abs(qimpl - qnum) > quadrupole_tol) then
         call test_failed(error, &
            "Gaussian quadrupole kernel disagrees with finite differences")
         return
      end if
   end do
end subroutine test_gaussian_multipole_finite_difference

subroutine test_gaussian_multipole_point_limit(error)
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: vec(3) = [1.25_wp, -0.85_wp, 1.70_wp]
   real(wp), parameter :: xi = 30.0_wp
   real(wp), parameter :: tol = 5.0e-13_wp
   real(wp) :: c0, c1(3), c2(6), r2, r1, ref2(6)

   call get_gaussian_multipole_kernels(vec, xi, c0, c1, c2)
   r2 = sum(vec**2)
   r1 = sqrt(r2)
   ref2 = [vec(1)*vec(1), 2.0_wp*vec(1)*vec(2), &
      & vec(2)*vec(2), 2.0_wp*vec(1)*vec(3), &
      & 2.0_wp*vec(2)*vec(3), vec(3)*vec(3)]/(r1*r2*r2)

   if (abs(c0 - 1.0_wp/r1) > tol &
      & .or. maxval(abs(c1 - vec/(r1*r2))) > tol &
      & .or. maxval(abs(c2 - ref2)) > tol) then
      call test_failed(error, &
         "Gaussian multipole kernels do not recover the point limit")
   end if
end subroutine test_gaussian_multipole_point_limit

pure function gaussian_monopole(center, xi) result(value)
   real(wp), intent(in) :: center(3)
   real(wp), intent(in) :: xi
   real(wp) :: value, distance

   distance = sqrt(sum(center**2))
   value = erf(xi*distance)/distance
end function gaussian_monopole

end module test_solvation_cosmo
