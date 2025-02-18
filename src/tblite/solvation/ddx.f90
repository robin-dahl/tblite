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

!> @file tblite/solvation/ddx.f90
!> Provides a polarizable continuum model

!> Implicit solvation model based on a polarizable dielectric continuum
module tblite_solvation_ddx
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use tblite_blas, only : dot, gemv
   use tblite_container_cache, only : container_cache
   use tblite_mesh_lebedev, only : grid_size, get_angular_grid, list_bisection
   use tblite_scf_info, only : scf_info, atom_resolved
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_solvation_data, only : get_vdw_rad_cosmo
   use tblite_solvation_type, only : solvation_type

   use omp_lib, only : omp_get_max_threads

   use iso_fortran_env, only : output_unit

   use tblite_disp_d4, only: get_eeq_charges


   use ddx, only: ddx_type, ddx_error_type, check_error, ddinit, ddx_state_type, allocate_state
   use ddx, only: ddrun
   use ddx, only: setup, fill_guess, solve, fill_guess_adjoint, solve_adjoint, solvation_force_terms
   use ddx_core, only: ddx_electrostatics_type
   use ddx_multipolar_solutes, only: multipole_electrostatics, multipole_force_terms, multipole_psi


   implicit none
   private

   public :: ddx_solvation, new_ddx, ddx_input, ddx_cache


   !> Input for ddX solvation
   type :: ddx_input
      !> Dielectric constant
      real(wp) :: dielectric_const
      !> ddx model
      integer :: ddx_model = 1
      !> Scaling of van-der-Waals radii
      real(wp) :: rscale = 1.0_wp
      !> Accuracy for iterative solver
      real(wp) :: conv = 1.0e-9_wp
      !> Regularization parameter
      real(wp) :: eta = 1.0_wp
      !> Number of grid points for each atom
      integer :: nang = grid_size(6)
      !> Maximum angular momentum of basis functions
      integer :: lmax = 1
      !> Van-der-Waals radii for all species
      real(wp), allocatable :: rvdw(:)
      !> Number of OMP threads
      integer :: nproc = 1
   end type ddx_input

   !> Definition of polarizable continuum model
   type, extends(solvation_type) :: ddx_solvation
      ! !> ddX instance
      type(ddx_input) :: ddx_input
      !> Dielectric function
      real(wp) :: keps
      !> Dielctric constant
      real(wp) :: dielectric_const
      !> Van-der-Waal radii for all atoms
      real(wp), allocatable :: rvdw(:)
   contains
      !> Update cache from container
      procedure :: update
      !> Return dependency on density
      procedure :: variable_info
      !> Get electric field energy
      procedure :: get_energy
      !> Get electric field potential
      procedure :: get_potential
      !> Get electric field gradient
      procedure :: get_gradient
   end type ddx_solvation

   !> Provide constructor for ddX solvation
   interface ddx_solvation
      module procedure :: create_ddx
   end interface ddX_solvation

   !> Restart data for ddX calculation
   type :: ddx_cache
      !> ddX instance
      type(ddx_type) :: ddx 
      !> ddX container with quantities common to all models
      type(ddx_state_type) :: ddx_state
      !> ddX container for the electrostatic properties  
      type(ddx_electrostatics_type) :: ddx_electrostatics
      !> ddX error handling
      type(ddx_error_type) :: ddx_error
      !> Interaction matrix with surface charges jmat(ncav, nat)
      real(wp), allocatable :: jmat(:, :)
      !> Solvation energy as returned by ddx
      real(wp) :: esolv
      !> ddx 
      real(wp), allocatable :: zeta(:)
      !> ddx multipole, (1, mol%nat)
      real(wp), allocatable :: multipoles(:, :)
      !> ddx forces (i.e. gradient of the solvation energy)
      real(wp), allocatable :: force(:, :)
   end type ddx_cache

   real(wp), parameter :: alpha_alpb = 0.571412_wp

contains

!> Create new electric field container
subroutine new_ddx(self, mol, input, error)
   !> Instance of the solvation model
   type(ddx_solvation), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input for ddX solvation
   type(ddx_input), intent(in) :: input
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, izp

   ! Set label
   if (input%ddx_model == 1) then
      self%label = "ddcpcm/ddcosmo solvation model"
   else if (input%ddx_model == 2) then
      self%label = "ddpcm solvation model"
   else if (input%ddx_model == 3) then
      self%label = "ddlpb solvation model"
   end if

   ! Get number of OMP threads
   self%ddx_input%nproc = omp_get_max_threads()

   ! Get radii for all atoms
   allocate(self%rvdw(mol%nat))
   self%rvdw(mol%nat) = 0.0_wp
   if (allocated(input%rvdw)) then
      self%rvdw(:) = input%rscale * input%rvdw(mol%id)
   else
      do iat = 1, mol%nat
         izp = mol%num(mol%id(iat))
         self%rvdw(iat) = input%rscale * get_vdw_rad_cosmo(izp)
      end do
   end if

   ! Calculate Epsilon and keps
   self%dielectric_const = input%dielectric_const
   ! self%keps = -(1.0_wp/self%dielectric_const - 1.0_wp) / (1.0_wp + alpha_alpb)
   self%keps = (self%dielectric_const - 1.0_wp) / (self%dielectric_const + 0.5_wp)



end subroutine new_ddx

!> Type constructor for ddX splvation
function create_ddx(mol, input) result(self)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input for ddX solvation
   type(ddx_input), intent(in) :: input
   !> Instance of the solvation model
   type(ddx_solvation) :: self
   !> Error handling
   type(error_type), allocatable :: error

   !> Create new instance of the solvation model
   call new_ddx(self, mol, input, error)
   if (allocated(error)) then
      call fatal_error(error)
   end if

end function create_ddx


!> Update cache from container
subroutine update(self, mol, cache)
   !> Instance of the solvation model
   class(ddx_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   type(ddx_cache), pointer :: ptr
   call taint(cache, ptr)

   !> Electrostatics at the cavity points
   ! Electric potential
   ptr%ddx_electrostatics%do_phi = .true.
   ! Electric field
   ptr%ddx_electrostatics%do_e = .true.
   ! Electric field gradient
   ptr%ddx_electrostatics%do_g = .true.


   ! initialize ddx
   call ddinit(self%ddx_input%ddx_model, mol%nat, mol%xyz, self%rvdw, self%dielectric_const, ptr%ddx, &
      & ptr%ddx_error, force=1, &
      & ngrid=self%ddx_input%nang, lmax=self%ddx_input%lmax, &
      & nproc=self%ddx_input%nproc , eta=self%ddx_input%eta)
   call check_error(ptr%ddx_error)

   ! Allocate forces
   if (allocated(ptr%force)) then
      deallocate(ptr%force)
   end if
   allocate(ptr%force(3, ptr%ddx%params%nsph))

   ! Initialize ddx_state
   call allocate_state(ptr%ddx%params, ptr%ddx%constants,  &
      ptr%ddx_state, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Allocate and get coulomb matrix
   if (allocated(ptr%jmat)) then
      deallocate(ptr%jmat)
   end if
   allocate(ptr%jmat(ptr%ddx%constants%ncav, mol%nat), source=0.0_wp)
   call get_coulomb_matrix(mol%xyz, ptr%ddx%constants%ccav, ptr%jmat)

   ! Allocate zeta
   ! allocate(ptr%zeta(ptr%ddx%constants%ncav))
   ! if (allocated(ptr%zeta)) then
   !    deallocate(ptr%zeta)
   ! end if
   ! allocate(ptr%zeta(ptr%ddx%params%nsph*ptr%ddx%params%ngrid))

   ! Allocate multipoles
   if (allocated(ptr%multipoles)) then
      deallocate(ptr%multipoles)
   end if
   allocate(ptr%multipoles(1, mol%nat))


   call fill_guess(ptr%ddx%params, ptr%ddx%constants, &
         & ptr%ddx%workspace, ptr%ddx_state, self%ddx_input%conv, ptr%ddx_error)
      call check_error(ptr%ddx_error)

      call fill_guess_adjoint(ptr%ddx%params, ptr%ddx%constants, &
         & ptr%ddx%workspace, ptr%ddx_state, self%ddx_input%conv, ptr%ddx_error)
      call check_error(ptr%ddx_error)


end subroutine update

!> Get electric field energy
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the solvation model
   class(ddx_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Solvation free energy
   real(wp), intent(inout) :: energies(:)
   real(wp), external :: ddot

   logical, save :: first_scf = .true.

   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   type(ddx_cache), pointer :: ptr

   call taint(cache, ptr)

   ! write(*,*) 'X energy beginning', ptr%ddx_state%xs
   ! write(*,*) 'Psi energy beginning', ptr%ddx_state%psi

   ptr%ddx_electrostatics%phi_cav = 0.0_wp
   ! multipoles_electrostatics already computes phi_cav, so no need to recompute here. However, the version match
   call get_phi(wfn%qat(:, 1), ptr%jmat, ptr%ddx_electrostatics%phi_cav)

   ! Calculate psi (representation of the solute density in spherical harmonics)
   ptr%ddx_state%psi = 0.0_wp
   call get_psi(wfn%qat(:, 1), ptr%ddx_state%psi)

   
   ptr%multipoles(1, :) = wfn%qat(:, 1) / sqrt(4.0_wp*pi)
   call multipole_electrostatics(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%multipoles, 0, ptr%ddx_electrostatics, ptr%ddx_error)

   call setup(ptr%ddx%params,ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, ptr%ddx_electrostatics, &
      & ptr%ddx_state%psi, ptr%ddx_error)
   call check_error(ptr%ddx_error)


   ! Calculate all solvation terms
   ! todo: this currently also always calculates all force terms
   ! call ddrun(ptr%ddx, ptr%ddx_state, ptr%ddx_electrostatics, ptr%ddx_state%psi, &
   !    self%ddx_input%conv, ptr%esolv, ptr%ddx_error, ptr%force) 
   ! call check_error(ptr%ddx_error)

   call solve(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, self%ddx_input%conv, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   !> Add solvation energy to total energy
   ! The factor of 0.5 is correct here
   energies(:) = energies +  self%keps * 0.5_wp * sum(ptr%ddx_state%xs * ptr%ddx_state%psi, 1) 
   ! energies(:) = energies +  self%keps * 0.5_wp * dot_product(ptr%ddx_state%xs(1,:), ptr%ddx_state%psi(1,:)) 

   ! write(*,*) 'X energy end', ptr%ddx_state%xs
   ! write(*,*) 'Psi energy end', ptr%ddx_state%psi
end subroutine get_energy

!> Get electric field potential
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the solvation model
   class(ddx_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   logical, save :: first_scf = .true.
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   type(ddx_cache), pointer :: ptr

   call taint(cache, ptr)

   ! write(*,*) 'X pot beginning', ptr%ddx_state%xs
   ! write(*,*) 'Psi pot beginning', ptr%ddx_state%psi

   ptr%ddx_state%psi = 0.0_wp
   call get_psi(wfn%qat(:, 1), ptr%ddx_state%psi)
   ptr%multipoles(1, :) = wfn%qat(:, 1) / sqrt(4.0_wp*pi)
   call multipole_electrostatics(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%multipoles, 0, ptr%ddx_electrostatics, ptr%ddx_error)

      call setup(ptr%ddx%params,ptr%ddx%constants, &
         & ptr%ddx%workspace, ptr%ddx_state, ptr%ddx_electrostatics, &
         & ptr%ddx_state%psi, ptr%ddx_error)
      call check_error(ptr%ddx_error)

   call solve(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, self%ddx_input%conv, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   call solve_adjoint(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, self%ddx_input%conv, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Contract with the Coulomb matrix
   call gemv(ptr%jmat, ptr%ddx_state%zeta, pot%vat(:, 1), alpha=-1.0_wp, beta=1.0_wp, trans='t') ! to give second term of 2-30

   pot%vat(:, 1) = 0.5_wp * self%keps * (pot%vat(:, 1) + sqrt(4.0_wp*pi) * ptr%ddx_state%xs(1, :)) ! see Eq. 2-30

   ! write(*,*) 'X pot end', ptr%ddx_state%xs
   ! write(*,*) 'Psi pot end', ptr%ddx_state%psi

end subroutine get_potential

!> Get electric field gradient
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the solvation model
   class(ddx_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the solvation free energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the solvation free energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   type(ddx_cache), pointer :: ptr
   call view(cache, ptr)

   call solvation_force_terms(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, ptr%ddx_electrostatics, ptr%force, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Given a multipolar distribution in real spherical harmonics and centered on the
   ! spheres, Compute the contributions to the forces stemming from electrostatic 
   ! interactions with multipolar distribution in real spherical harmonics and centered ! on the spheres
   call multipole_force_terms(ptr%ddx%params, ptr%ddx%constants, ptr%ddx%workspace, &
      ptr%ddx_state, 0, ptr%multipoles, ptr%force, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Calculate the gradient of the solvation energy
   ! We needed the factor 2 here because 0.5 was falsly included in keps
   ptr%force = self%keps * ptr%force 

   ! call write_2d_matrix(ptr%force, "force", unit=output_unit, step=5)

   ! Add the gradient of the solvation energy to the total gradient
   gradient =  gradient + ptr%force

   ! call test_potential

end subroutine get_gradient


!> Return dependency on density
pure function variable_info(self) result(info)
   !> Instance of the solvation model
   class(ddx_solvation), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=atom_resolved)
end function variable_info


subroutine taint(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(ddx_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(ddx_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

subroutine view(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(ddx_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(ddx_cache)
      ptr => target
   end select
end subroutine view

!> Evaluate the Coulomb interactions between the atomic sides (xyz) and the
!> surface elements of the cavity (ccav).
subroutine get_coulomb_matrix(xyz, ccav, jmat)
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: ccav(:, :)
   real(wp), intent(inout) :: jmat(:, :)

   integer :: ic, j
   real(wp) :: vec(3), d2, d

   jmat(:, :) = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(ccav, xyz, jmat) private(ic, j, vec, d2, d)
   do ic = 1, size(ccav, 2)
      do j = 1, size(xyz, 2)
         vec(:) = ccav(:, ic) - xyz(:, j)
         d2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         d = sqrt(d2)
         jmat(ic, j) = 1.0_wp / d
      end do
   end do

end subroutine get_coulomb_matrix

!> Routine to compute the psi vector
subroutine get_psi(charge, psi)
   real(wp), intent(in) :: charge(:)
   real(wp), intent(out) :: psi(:, :)

   integer :: iat
   real(wp), parameter :: fac = sqrt(4.0_wp*pi) 

   psi(:,:) = 0.0_wp

   do iat = 1, size(charge)
      psi(1, iat) = fac*charge(iat)
   end do

end subroutine get_psi

!> Routine to compute the potential vector
subroutine get_phi(charge, jmat, phi)
   real(wp), intent(in) :: charge(:)
   real(wp), intent(in) :: jmat(:, :)
   real(wp), intent(out) :: phi(:)

   phi(:) = 0.0_wp

   call gemv(jmat, charge, phi)

end subroutine get_phi


subroutine get_zeta(self, keps, zeta)
   type(ddx_cache), intent(in) :: self
   real(wp), intent(in) :: keps
   real(wp), intent(out) :: zeta(self%ddx%constants%ncav)

   integer :: iat, its, ii

   zeta = 0.0_wp

   ii = 0
   do iat = 1, self%ddx%params%nsph
      do its = 1, self%ddx%params%ngrid
         if (self%ddx%constants%ui(its, iat) > 0.0_wp) then
            ii = ii + 1
            zeta(ii) = self%ddx%constants%wgrid(its) &
               & * self%ddx%constants%ui(its, iat) &
               & * dot_product(self%ddx%constants%vgrid(:, its), self%ddx_state%s(:, iat))
         end if
      end do
   end do

 end subroutine get_zeta


subroutine write_solvation_file()

end subroutine write_solvation_file




!**************************************************************
!********************* Dev Routines ***************************
!**************************************************************

subroutine write_vector(vector, name, unit)
   implicit none
   real(wp),intent(in) :: vector(:)
   character(len=*),intent(in),optional :: name
   integer, intent(in),optional :: unit
   integer :: d
   integer :: i, j, k, l, istep, iunit

   d = size(vector, dim=1)

   if (present(unit)) then
       iunit = unit
   else
       iunit = output_unit
   end if

   if (present(name)) write(iunit,'(/,"vector:",1x,a)') name

   do j = 1, d
      write(iunit, '(i6)', advance='no') j
      write(iunit, '(1x,f15.8)', advance='no') vector(j)
      write(iunit, '(a)')
   end do

end subroutine write_vector


subroutine write_2d_matrix(matrix, name, unit, step)
   implicit none
   real(wp),intent(in) :: matrix(:, :)
   character(len=*),intent(in),optional :: name
   integer, intent(in),optional :: unit
   integer, intent(in),optional :: step
   integer :: d1, d2
   integer :: i, j, k, l, istep, iunit

   d1 = size(matrix, dim=1)
   d2 = size(matrix, dim=2)

   if (present(unit)) then
       iunit = unit
   else
       iunit = output_unit
   end if

   if (present(step)) then
       istep = step
   else
       istep = 5
   end if

   if (present(name)) write(iunit,'(/,"matrix:",1x,a)') name

   do i = 1, d2, istep
      l = min(i+istep-1,d2)
      write(iunit,'(6x)',advance='no')
      do k = i, l
         write(iunit,'(6x,i12,3x)',advance='no') k
      end do
      write(iunit,'(a)')
      do j = 1, d1
         write(iunit,'(i6)',advance='no') j
         do k = i, l
            write(iunit,'(1x,e20.8)',advance='no') matrix(j,k)
         end do
         write(iunit,'(a)')
      end do
   end do

end subroutine write_2d_matrix

end module tblite_solvation_ddx
