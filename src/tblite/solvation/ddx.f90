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
   use tblite_mesh_lebedev, only : grid_size
   use tblite_scf_info, only : scf_info, atom_resolved
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_solvation_data, only : get_vdw_rad_cosmo
   use tblite_solvation_type, only : solvation_type

   use omp_lib, only : omp_get_max_threads

   use ddx, only: ddx_type, ddx_error_type, check_error, ddinit, ddx_state_type, allocate_state, ddrun
   use ddx, only: setup, fill_guess, solve, fill_guess_adjoint, solve_adjoint, solvation_force_terms
   use ddx_core, only: ddx_electrostatics_type
   use ddx_multipolar_solutes, only: multipole_electrostatics, multipole_force_terms, multipole_psi

   use ddx_gradients, only: get_dgriddr, get_dvdr, get_dtdr, get_dsdr, get_duidr, get_dphidr, get_dgdr, &
      & contract_grad_l_dynrad

   ! ! for DRACO 
   use tblite_ncoord, only: new_ncoord, ncoord_type
   use tblite_solvation_draco, only: new_draco, draco_type
   use tblite_solvation_draco_type, only: draco_input, create_draco_input
   use tblite_solvation_data_draco, only: get_draco_param
   use tblite_disp_d4, only: get_eeq_charges


   implicit none
   private

   public :: ddx_solvation, new_ddx, ddx_input, ddx_cache, ddx_solvation_model


   !> Possible solvation models to be used within the dd framework
   type :: enum_ddx_solvation_model
      ! COSMO and CPCM temporary defined as 11 and 12 to get correct feps in the first step
      ! Labels are dumbed to 1 before passing the model input to the ddX routine,
      ! where both models are handeled the same way
      !> Conductor like screening model
      integer :: cosmo = 11
      !> Conductor-like polarizable continuum model
      integer :: cpcm = 12
      !> Polarizable continuum model
      integer :: pcm = 2
      !> Linearized Poisson-Boltzmann model
      integer :: lpb = 3
   end type enum_ddx_solvation_model

   !> Actual enumerator for the dd solvation models
   type(enum_ddx_solvation_model), parameter :: ddx_solvation_model = enum_ddx_solvation_model()


   !> Input for ddX solvation
   type :: ddx_input
      !> Dielectric constant
      real(wp) :: dielectric_const
      !> ddx model
      integer :: ddx_model = ddx_solvation_model%cosmo
      !> Scaling of van-der-Waals radii
      real(wp) :: rscale = 1.0_wp
      !> Accuracy for iterative solver
      real(wp) :: conv = 1.0e-10_wp
      !> Regularization parameter
      real(wp) :: eta = 0.1_wp
      !> Number of grid points for each atom
      integer :: nang = grid_size(8)
      !> Maximum angular momentum of basis functions
      integer :: lmax = 1
      !> Van-der-Waals radii for all species
      real(wp), allocatable :: rvdw(:)
      !> Number of OMP threads
      integer :: nproc = 1
      !> Debye-Hückel screening length (only used for LPB)
      real(wp) :: kappa = 0.0_wp
      !> Shift of the characteristic function 
      ! (default value depends on the model, for COSMO/CPCM it is -1)
      real(wp) :: shift = -1.0_wp
      !> Maximum number of iterations for the iterative solver
      integer :: max_iter = 100
      !> Number of extrapolation points for the Jacobi/DIIS solver
      integer :: jacobi_ndiis = 20
      !> Maximal degree of multipole spherical harmonics
      integer :: pm = 10
      !> Maximal degree of local spherical harmonics
      integer :: pl = 10
      !> Handling of the sparse matrices 
      integer :: incore = 0
      !> 1 to use FMM acceleration and 0 otherwise
      integer :: enable_fmm = 0
   end type ddx_input

   !> Definition of polarizable continuum model
   type, extends(solvation_type) :: ddx_solvation
      !> ddX instance
      type(ddx_input) :: ddx_input
      !> Dielectric function
      real(wp) :: feps
      !> Dielctric constant
      real(wp) :: dielectric_const
      !> Van-der-Waal radii for all atoms
      real(wp), allocatable :: rvdw(:)
      !> Derivative of the radii
      real(wp), allocatable :: drdr(:,:,:)
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
      !> ddX potential
      real(wp), allocatable :: ddx_pot(:)
      !> Solvation energy as returned by ddx
      real(wp) :: esolv
      !> ddx multipole, (1, mol%nat)
      real(wp), allocatable :: multipoles(:, :)
      !> ddx forces (i.e. gradient of the solvation energy)
      real(wp), allocatable :: force(:, :)
   end type ddx_cache

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

   !> DRACO type
   type(draco_type) :: draco
   type(draco_input) :: input_draco

   class(ncoord_type), allocatable :: ncoord
   real(wp) :: cn(mol%nat), q(mol%nat)
   real(wp), allocatable :: radii(:)
   type(wavefunction_type):: wfn
   real(wp), allocatable :: draddr(:,:,:), dcndr(:,:,:), dqdr(:,:,:)

   integer :: iat, izp, nlines
   real(wp) :: feps_param 
   character(len=100) :: line

   ! Set label
   if (input%ddx_model == ddx_solvation_model%cosmo) then
      self%label = "ddcosmo solvation model"
   else if (input%ddx_model == ddx_solvation_model%cpcm) then
      self%label = "ddcpcm solvation model"
   else if (input%ddx_model == ddx_solvation_model%pcm) then
      self%label = "ddpcm solvation model"
      else if (input%ddx_model == ddx_solvation_model%lpb) then
      self%label = "ddlpb solvation model"
   end if

   ! Set model 
   self%ddx_input%ddx_model = input%ddx_model

   ! Get number of OMP threads
   self%ddx_input%nproc = omp_get_max_threads()

   ! Get radii for all atoms
   allocate(self%rvdw(mol%nat), source=0.0_wp)
   if (allocated(input%rvdw)) then
      self%rvdw(:) = input%rscale * input%rvdw(mol%id)
   else
      do iat = 1, mol%nat
         izp = mol%num(mol%id(iat))
         self%rvdw(iat) = input%rscale * get_vdw_rad_cosmo(izp)
      end do
   end if

   open(unit=95, file='draco_radii.txt', status='old', form='formatted')
   nlines = 0
   do
      read(95, '(A)', iostat=iat) line
      if (iat /= 0) exit
      nlines = nlines + 1
   end do
   close(95)

   ! Allocate and read values
   open(unit=95, file='draco_radii.txt', status='old', form='formatted')
   do iat = 1, nlines
      read(95, *) self%rvdw(iat)
   end do
   close(95)

   ! Convert radii to Bohr
   self%rvdw = self%rvdw / 0.529177210903 ! + 1.5_wp

   ! Calculate epsilon and feps
   self%dielectric_const = input%dielectric_const
   if (input%ddx_model == ddx_solvation_model%cosmo) then
      feps_param = 0.5_wp
      self%feps = (self%dielectric_const - 1.0_wp) / (self%dielectric_const + feps_param)
   else if (input%ddx_model == ddx_solvation_model%cpcm) then
      feps_param = 0.0_wp
      self%feps = (self%dielectric_const - 1.0_wp) / (self%dielectric_const + feps_param)
   else 
      self%feps = 1.0_wp
   end if

   ! Get Debye-Hückel screening length (only used for LPB)
   self%ddx_input%kappa = input%kappa

   ! Initialize the shift depending on the model: ddCOSMO/ddCPCM has an internal shift, 
   ! ddPCM and ddLPB have a symmetric shift
   if (input%ddx_model == ddx_solvation_model%cosmo .or. &
      & input%ddx_model == ddx_solvation_model%cpcm) then
      self%ddx_input%shift = -1.0_wp
   else 
      self%ddx_input%shift = 0.0_wp
   end if

   ! Initialize the rest of the input
   self%ddx_input%conv = input%conv
   self%ddx_input%eta = input%eta
   self%ddx_input%nang = input%nang
   self%ddx_input%lmax = input%lmax
   self%ddx_input%max_iter = input%max_iter
   self%ddx_input%jacobi_ndiis = input%jacobi_ndiis
   self%ddx_input%pm = input%pm
   self%ddx_input%pl = input%pl
   self%ddx_input%incore = input%incore
   self%ddx_input%enable_fmm = input%enable_fmm

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

   ! Create new instance of the solvation model
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

   integer :: model

   call taint(cache, ptr)

   ! Electrostatics at the cavity points
   ! Electric potential
   ptr%ddx_electrostatics%do_phi = .true.
   ! Electric field
   ptr%ddx_electrostatics%do_e = .true.
   ! Electric field gradient
   ptr%ddx_electrostatics%do_g = .true.

   ! Put COSMO and CPCM label back to 1 
   if (self%ddx_input%ddx_model == ddx_solvation_model%cosmo .or. &
      & self%ddx_input%ddx_model == ddx_solvation_model%cpcm) then
      model = 1
   else
      model = self%ddx_input%ddx_model
   end if

   write(*,*) 'gs:', self%ddx_input%nang
   write(*,*) 'fmm:', self%ddx_input%enable_fmm 
   call ddinit(model, mol%nat, mol%xyz, &
      & self%rvdw, self%dielectric_const, ptr%ddx, &
      & ptr%ddx_error, force=1, ngrid=self%ddx_input%nang, &
      & lmax=self%ddx_input%lmax, nproc=self%ddx_input%nproc, &
      & eta=self%ddx_input%eta, kappa=self%ddx_input%kappa, &
      & shift=self%ddx_input%shift, maxiter=self%ddx_input%max_iter, &
      & jacobi_ndiis=self%ddx_input%jacobi_ndiis, pm=self%ddx_input%pm, &
      & pl=self%ddx_input%pl, incore=self%ddx_input%incore, &
      & enable_fmm=self%ddx_input%enable_fmm)
   call check_error(ptr%ddx_error)

   call allocate_state(ptr%ddx%params, ptr%ddx%constants, &
      ptr%ddx_state, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   if (allocated(ptr%multipoles)) then
      deallocate(ptr%multipoles)
   end if
   allocate(ptr%multipoles(1, mol%nat), source=0.0_wp)

   if (allocated(ptr%jmat)) then
      deallocate(ptr%jmat)
   end if
   allocate(ptr%jmat(ptr%ddx%constants%ncav, mol%nat), source=0.0_wp)
   call get_coulomb_matrix(mol%xyz, ptr%ddx%constants%ccav, ptr%jmat)

   if (allocated(ptr%force)) then
      deallocate(ptr%force)
   end if
   allocate(ptr%force(3, mol%nat), source=0.0_wp)

   if (allocated(ptr%ddx_pot)) then
      deallocate(ptr%ddx_pot)
   end if
   allocate(ptr%ddx_pot(mol%nat), source=0.0_wp)

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
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   type(ddx_cache), pointer :: ptr

   call view(cache, ptr)


   ! Recalculate the solution of the ddX system with the new charges after diagonalization
   ! This solution cannot be reused in the potential due to intermediate mixing

   ptr%multipoles(1, :) = wfn%qat(:, 1) / sqrt(4.0_wp*pi)
   call multipole_electrostatics(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%multipoles, 0, ptr%ddx_electrostatics, ptr%ddx_error)
   call multipole_psi(ptr%ddx%params, ptr%multipoles, 0, ptr%ddx_state%psi)
   
   call setup(ptr%ddx%params,ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, ptr%ddx_electrostatics, &
      & ptr%ddx_state%psi, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   call solve(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, self%ddx_input%conv, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Add solvation energy to total energy
   if (self%ddx_input%ddx_model == ddx_solvation_model%lpb) then
      energies(:) = energies + self%feps * 0.5_wp * sum(ptr%ddx_state%x_lpb(:,:,1) * ptr%ddx_state%psi, 1)
   else
      energies(:) = energies + self%feps * 0.5_wp * sum(ptr%ddx_state%xs * ptr%ddx_state%psi, 1) 
   end if

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
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   type(ddx_cache), pointer :: ptr

   call view(cache, ptr)

   ! Solution of the ddX system (direct and adjoint) with the mixed charges
   ! This solution cannot be reused in the energy calculation due to intermediate diagonalization

   ptr%multipoles(1, :) = wfn%qat(:, 1) / sqrt(4.0_wp*pi)
   call multipole_electrostatics(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%multipoles, 0, ptr%ddx_electrostatics, ptr%ddx_error)
   call multipole_psi(ptr%ddx%params, ptr%multipoles, 0, ptr%ddx_state%psi)

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
   ptr%ddx_pot = 0.0_wp
   call gemv(ptr%jmat, ptr%ddx_state%zeta, ptr%ddx_pot(:), alpha=-1.0_wp, beta=1.0_wp, trans='t') 
   ! Scale with 0.5 and feps, and get second contribution to potential
   if (self%ddx_input%ddx_model == ddx_solvation_model%lpb) then
      ptr%ddx_pot(:) = 0.5_wp * self%feps * (ptr%ddx_pot(:) + sqrt(4.0_wp*pi) * ptr%ddx_state%x_lpb(1, :, 1))
   else
      ptr%ddx_pot(:) = 0.5_wp * self%feps * (ptr%ddx_pot(:) + sqrt(4.0_wp*pi) * ptr%ddx_state%xs(1, :))
   end if
 
   ! Add potential to overall potential for new SCF step 
   pot%vat(:,1) = pot%vat(:,1) + ptr%ddx_pot(:)

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
   ! type(ddx_type) :: ptr%ddx

   ! needed for gradient
   real(wp) :: eta
   real(wp), allocatable :: rvdw(:), grid(:,:) 
   real(wp), allocatable :: xyz(:,:) 


   real(wp), allocatable :: v_tensor(:,:,:,:), s_tensor(:,:,:,:), t_tensor(:,:,:), vv_tensor(:,:,:), &
                        & chi_tensor(:,:,:)     
   real(wp), allocatable :: drvdwdr(:,:,:), dgriddr(:,:,:,:,:), dvdr(:,:,:,:,:,:), dtdr(:,:,:,:,:), &
                        & dwdr(:,:,:,:,:), duidr(:,:,:,:), dfidr(:,:,:,:), dphidr(:,:,:), &
                        & dgdr(:,:,:,:), dsdr(:,:,:,:,:,:), dlxdr(:,:), dgdr_con(:,:), dedr(:,:)

   integer :: nsph, isph, ibasis, k, asph, bsph, igrid

   call view(cache, ptr)

   ptr%force = 0.0_wp
   
   ptr%multipoles(1, :) = wfn%qat(:, 1) / sqrt(4.0_wp*pi)
   call multipole_electrostatics(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%multipoles, 0, ptr%ddx_electrostatics, ptr%ddx_error)
   call multipole_psi(ptr%ddx%params, ptr%multipoles, 0, ptr%ddx_state%psi)

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
   
   call solvation_force_terms(ptr%ddx%params, ptr%ddx%constants, &
     & ptr%ddx%workspace, ptr%ddx_state, ptr%ddx_electrostatics, ptr%force, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   call multipole_force_terms(ptr%ddx%params, ptr%ddx%constants, ptr%ddx%workspace, &
     ptr%ddx_state, 0, ptr%multipoles, ptr%force, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   allocate(drvdwdr(3, ptr%ddx%params%nsph, ptr%ddx%params%nsph), source=0.0_wp)
   open(unit=97, file='drdr_0.dat', form='unformatted', status='replace', action='write')
      write(97) drvdwdr
   close(97)
   open(unit=98, file='drdr.dat', access='stream', form='unformatted', status='old')
      read(98) drvdwdr
   close(98)
   ! Convert to Bohr
   drvdwdr = drvdwdr / 0.529177210903

   allocate(rvdw(ptr%ddx%params%nsph), grid(3, ptr%ddx%params%ngrid), xyz(3, ptr%ddx%params%nsph))

   eta = ptr%ddx%params%eta

   nsph =  ptr%ddx%params%nsph
   xyz = mol%xyz 
   rvdw = self%rvdw
   grid = ptr%ddx%constants%cgrid


   !%%%%%%%%%%%%%%%%%%%%%% GET METRICS %%%%%%%%%%%%%%%%%%%%%%
   allocate(v_tensor(nsph, nsph, ptr%ddx%params%ngrid, 3), source=0.0_wp)
   allocate(s_tensor(nsph, nsph, ptr%ddx%params%ngrid, 3), source=0.0_wp)
   allocate(t_tensor(nsph, nsph, ptr%ddx%params%ngrid), source=0.0_wp)
   allocate(vv_tensor(nsph, nsph, ptr%ddx%params%ngrid), source=0.0_wp)
   allocate(chi_tensor(ptr%ddx%params%nsph, ptr%ddx%params%nsph, ptr%ddx%params%ngrid), source=0.0_wp)
   do asph = 1, ptr%ddx%params%nsph
      do bsph = 1, ptr%ddx%params%nsph
         do igrid = 1, ptr%ddx%params%ngrid

            v_tensor(asph, bsph, igrid, :) = xyz(:, asph) &
               & + rvdw(asph) * grid(:,igrid) &
               & - xyz(:, bsph)

            vv_tensor(asph, bsph, igrid) = sqrt(dot_product(v_tensor(asph, bsph, igrid, :), v_tensor(asph, bsph, igrid, :)))
            t_tensor(asph, bsph, igrid) = vv_tensor(asph, bsph, igrid) / rvdw(bsph)
            s_tensor(asph, bsph, igrid,:) = v_tensor(asph, bsph, igrid, :) / vv_tensor(asph, bsph, igrid)

            if (t_tensor(asph, bsph, igrid) <= (1.0_wp-eta)) then
               chi_tensor(asph, bsph, igrid) = 1.0_wp
            
            else if (((1.0_wp-eta) < t_tensor(asph, bsph, igrid)) .and. (t_tensor(asph, bsph, igrid) < 1.0_wp)) then
               chi_tensor(asph, bsph, igrid) = eta**(-5) * (1.0_wp-t_tensor(asph, bsph, igrid))**3 &
                  & * ( 6.0_wp*t_tensor(asph, bsph, igrid)**2 + (15.0_wp*eta-12.0_wp)*t_tensor(asph, bsph, igrid) &
                  & + 10.0_wp*eta**2 - 15.0_wp*eta + 6.0_wp )
            
            else if (t_tensor(asph, bsph, igrid) >= 1.0_wp) then
               chi_tensor(asph, bsph, igrid) = 0.0_wp

            end if

         end do
      end do
   end do  
   

   !%%%%%%%%%%%%%%%%%%%%%% GET GRADIENTS %%%%%%%%%%%%%%%%%%%%%%
   allocate(dgriddr(nsph, 3, nsph, ptr%ddx%params%ngrid, 3), source = 0.0_wp)
   call get_dgriddr(ptr%ddx, xyz, rvdw, grid, drvdwdr, dgriddr)
   
   allocate(dvdr(nsph, 3, nsph, nsph, ptr%ddx%params%ngrid, 3), source=0.0_wp)
   call get_dvdr(ptr%ddx, rvdw, grid, drvdwdr, dgriddr, dvdr)

   allocate(dtdr(nsph, 3, nsph, nsph, ptr%ddx%params%ngrid), source=0.0_wp)
   call get_dtdr(ptr%ddx, rvdw, v_tensor, vv_tensor, drvdwdr, dvdr, dtdr)

   allocate(dsdr(nsph, 3, nsph, nsph, ptr%ddx%params%ngrid, 3), source=0.0_wp)
   call get_dsdr(ptr%ddx, v_tensor, vv_tensor, dvdr, dsdr)

   allocate(dwdr(nsph, 3, nsph, nsph, ptr%ddx%params%ngrid), source=0.0_wp)
   allocate(dfidr(nsph, 3, nsph,  ptr%ddx%params%ngrid), source=0.0_wp)
   allocate(duidr(nsph, 3, nsph, ptr%ddx%params%ngrid), source=0.0_wp)
   call get_duidr(ptr%ddx, t_tensor, chi_tensor, dtdr, dwdr, dfidr, duidr)

   allocate(dphidr(ptr%ddx%params%nsph, 3, ptr%ddx%constants%ncav), source=0.0_wp)
   call get_dphidr(ptr%ddx, wfn%qat(:,1), v_tensor, vv_tensor, dvdr, dphidr)

   allocate(dgdr(nsph, 3, ptr%ddx%constants%nbasis, nsph), source=0.0_wp)
   call get_dgdr(ptr%ddx, ptr%ddx_electrostatics, wfn%qat(:,1), duidr, dphidr, dgdr)
   
   allocate(dlxdr(3, nsph), source=0.0_wp)
   do isph = 1, nsph
       call contract_grad_l_dynrad(ptr%ddx, isph,  ptr%ddx_state%xs, &
           & ptr%ddx_state % sgrid, ptr%ddx%workspace % tmp_vylm(:, 1), &
           & ptr%ddx%workspace % tmp_vdylm(:, :, 1), ptr%ddx%workspace % tmp_vplm(:, 1), &
           & ptr%ddx%workspace % tmp_vcos(:, 1), ptr%ddx%workspace % tmp_vsin(:, 1), &
           & dlxdr(:, isph), dtdr=dtdr, dsdr=dsdr, dfidr=dfidr) 
   end do


   !%%%%%%%%%%%%%%%%%%%%%% PUT TOGETHER %%%%%%%%%%%%%%%%%%%%%%
   allocate(dgdr_con(3, nsph), source=0.0_wp)
   do isph = 1, nsph
      do k = 1, 3
         do ibasis = 1, ptr%ddx%constants%nbasis
            do asph = 1, nsph
               dgdr_con(k, isph) = dgdr_con(k, isph) + &
                  & dgdr(isph, k, ibasis, asph) * ptr%ddx_state%s(ibasis, asph)
            end do
         end do
      end do
   end do

   allocate(dedr(3, nsph), source=0.0_wp)
   dedr = 0.5_wp * self%feps * (dgdr_con - dlxdr)
   
   ! Calculate the gradient of the solvation energy
   ptr%force = self%feps * ptr%force 

   ! Add the gradient of the solvation energy to the total gradient
   gradient = gradient + dedr 

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


end module tblite_solvation_ddx
