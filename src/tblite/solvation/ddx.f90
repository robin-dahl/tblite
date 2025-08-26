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

   use tblite_integral_trafo, only: transform0


   implicit none
   private

   public :: ddx_solvation, new_ddx, ddx_input, ddx_cache, ddx_solvation_model


   !> Possible solvation models to be used within the dd framework
   type :: enum_ddx_solvation_model
      ! COSMO and CPCM temporary defined as 11 and 12 to get correct feps in the first step
      ! Labels are dumped to 1 before passing the model input to the ddX routine,
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
      !> Number of grid points for each atom (=110)
      integer :: nang = grid_size(8)
      !> Maximum angular momentum of basis functions
      integer :: lmax = 2 ! 2
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
      integer :: pm = 8
      !> Maximal degree of local spherical harmonics
      integer :: pl = 8
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
      real(wp), allocatable :: adpmat(:, :, :)
      real(wp), allocatable :: aqpmat(:, :, :)
      !> ddX potential
      real(wp), allocatable :: ddx_pot(:)
      real(wp), allocatable :: ddx_dppot(:, :)
      real(wp), allocatable :: ddx_qppot(:, :)
      !> Solvation energy as returned by ddx
      real(wp) :: esolv
      !> ddx multipole, (1, mol%nat) -> To change: 1 comes from the multipole order, will be increased to 9
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

   integer :: iat, izp
   real(wp) :: feps_param 

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
   allocate(ptr%multipoles(9, mol%nat), source=0.0_wp) !9

   if (allocated(ptr%jmat)) then
      deallocate(ptr%jmat)
   end if
   allocate(ptr%jmat(ptr%ddx%constants%ncav, mol%nat), source=0.0_wp)
   call get_coulomb_matrix(mol%xyz, ptr%ddx%constants%ccav, ptr%jmat)
   if (allocated(ptr%adpmat)) then
      deallocate(ptr%adpmat)
   end if
   allocate(ptr%adpmat(3, ptr%ddx%constants%ncav, mol%nat), source=0.0_wp)
   call get_adp_matrix(mol%xyz, ptr%ddx%constants%ccav, ptr%adpmat)
   if (allocated(ptr%aqpmat)) then
      deallocate(ptr%aqpmat)
   end if
   allocate(ptr%aqpmat(5, ptr%ddx%constants%ncav, mol%nat), source=0.0_wp)
   call get_aqp_matrix(mol%xyz, ptr%ddx%constants%ccav, ptr%aqpmat)

   if (allocated(ptr%force)) then
      deallocate(ptr%force)
   end if
   allocate(ptr%force(3, mol%nat), source=0.0_wp)

   if (allocated(ptr%ddx_pot)) then
      deallocate(ptr%ddx_pot)
   end if
   allocate(ptr%ddx_pot(mol%nat), source=0.0_wp)

   if (allocated(ptr%ddx_dppot)) then
      deallocate(ptr%ddx_dppot)
   end if
   allocate(ptr%ddx_dppot(3, mol%nat), source=0.0_wp)

   if (allocated(ptr%ddx_qppot)) then
      deallocate(ptr%ddx_qppot)
   end if
   allocate(ptr%ddx_qppot(6, mol%nat), source=0.0_wp)

   call multipole_electrostatics(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%multipoles, 2, ptr%ddx_electrostatics, ptr%ddx_error)
   call multipole_psi(ptr%ddx%params, ptr%multipoles, 2, ptr%ddx_state%psi)

   call setup(ptr%ddx%params,ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, ptr%ddx_electrostatics, &
      & ptr%ddx_state%psi, ptr%ddx_error)
   call check_error(ptr%ddx_error)

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

   real(wp) :: fac, trans_fac(5)

   real(wp), parameter :: s3 = sqrt(3.0_wp)
   real(wp), parameter :: s3_4 = s3 * 0.5_wp
   real(wp), parameter :: dtrafo(5,6) =  sqrt(5.0_wp/(4.0_wp*pi)) * reshape([ &
         !   m=-2     m=-1      m=0      m=+1     m=+2
         & 0.0_wp,   0.0_wp,  -0.5_wp,  0.0_wp,   s3_4,   &  ! xx
         &    s3 ,   0.0_wp,   0.0_wp,  0.0_wp,   0.0_wp, &  ! xy
         & 0.0_wp,   0.0_wp,  -0.5_wp,  0.0_wp,  -s3_4,   &  ! yy
         & 0.0_wp,   0.0_wp,   0.0_wp,     s3 ,   0.0_wp, &  ! xz
         & 0.0_wp,      s3 ,   0.0_wp,  0.0_wp,   0.0_wp, &  ! yz
         & 0.0_wp,   0.0_wp,   1.0_wp,  0.0_wp,   0.0_wp  &  ! zz
   ], shape(dtrafo))


   real(wp), parameter :: tol = 1.0d-10
   integer :: i, k

   call view(cache, ptr)


   ! Recalculate the solution of the ddX system with the new charges after diagonalization
   ! This solution cannot be reused in the potential due to intermediate mixing

   ! Monopole
   ptr%multipoles(1, :) = wfn%qat(:, 1) / sqrt(4.0_wp*pi)
   ! Dipoles
   ! This prefactor should be correct, according to https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
   fac = sqrt(3.0_wp/(4.0_wp*pi))
   ptr%multipoles(2, :) = wfn%dpat(1,:,1) * fac
   ptr%multipoles(3, :) = wfn%dpat(2,:,1) * fac
   ptr%multipoles(4, :) = wfn%dpat(3,:,1) * fac
   
   ! Quadrupoles
   do i = 1, mol%nat
      ptr%multipoles(5:9,i) = matmul(dtrafo, wfn%qpat(:,i,1))
   end do

   call multipole_electrostatics(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%multipoles, 2, ptr%ddx_electrostatics, ptr%ddx_error) !2
   call multipole_psi(ptr%ddx%params, ptr%multipoles, 2, ptr%ddx_state%psi) !2

   call setup(ptr%ddx%params,ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, ptr%ddx_electrostatics, &
      & ptr%ddx_state%psi, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   call solve(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%ddx_state, self%ddx_input%conv, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Add solvation energy to total energy
   if (self%ddx_input%ddx_model == ddx_solvation_model%lpb) then
      energies(:) = energies + 0.5_wp * self%feps * sum(ptr%ddx_state%x_lpb(:,:,1) * ptr%ddx_state%psi, 1)
   else
      energies(:) = energies + 0.5_wp * self%feps * sum(ptr%ddx_state%xs * ptr%ddx_state%psi, 1)
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
   integer :: k, i
   real(wp) :: fac, trans_fac(5), ddx_qppot_trans(5,mol%nat)

   real(wp), parameter :: s3 = sqrt(3.0_wp)
   real(wp), parameter :: s3_4 = s3 * 0.5_wp
   real(wp), parameter :: dtrafo(5,6) = sqrt(5.0_wp/(4.0_wp*pi)) * reshape([ &
         !   m=-2     m=-1      m=0      m=+1     m=+2
         & 0.0_wp,   0.0_wp,  -0.5_wp,  0.0_wp,   s3_4,   &  ! xx
         &    s3 ,   0.0_wp,   0.0_wp,  0.0_wp,   0.0_wp, &  ! xy
         & 0.0_wp,   0.0_wp,  -0.5_wp,  0.0_wp,  -s3_4,   &  ! yy
         & 0.0_wp,   0.0_wp,   0.0_wp,     s3 ,   0.0_wp, &  ! xz
         & 0.0_wp,      s3 ,   0.0_wp,  0.0_wp,   0.0_wp, &  ! yz
         & 0.0_wp,   0.0_wp,   1.0_wp,  0.0_wp,   0.0_wp  &  ! zz
         ], shape(dtrafo))

   real(wp), parameter :: pinv(6,5) = transpose( sqrt(4.0_wp*pi/5.0_wp) * reshape([ &
        & 0.000000_wp, 0.000000_wp, -1.0_wp/3.0_wp, 0.000000_wp,  1.0_wp/s3,   & 
        & 1.0_wp/s3,   0.000000_wp,  0.000000_wp,   0.000000_wp,  0.000000_wp, &  
        & 0.000000_wp, 0.000000_wp, -1.0_wp/3.0_wp, 0.000000_wp, -1.0_wp/s3,   & 
        & 0.000000_wp, 0.000000_wp,  0.000000_wp,   1.0_wp/s3,    0.000000_wp, & 
        & 0.000000_wp, 1.0_wp/s3,    0.000000_wp,   0.000000_wp,  0.000000_wp, & 
        & 0.000000_wp, 0.000000_wp,  2.0_wp/3.0_wp, 0.000000_wp,  0.000000_wp  & 
        ], [5,6]) )

   call view(cache, ptr)

   ! Solution of the ddX system (direct and adjoint) with the mixed charges
   ! This solution cannot be reused in the energy calculation due to intermediate diagonalization

   ! Transformation into spherical harmonics, according to https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
   ! Monopole
   ptr%multipoles(1, :) = wfn%qat(:, 1) / sqrt(4.0_wp*pi)
   ! Dipoles
   fac = sqrt(3.0_wp/(4.0_wp*pi))
   ptr%multipoles(2, :) = wfn%dpat(1,:,1) * fac
   ptr%multipoles(3, :) = wfn%dpat(2,:,1) * fac
   ptr%multipoles(4, :) = wfn%dpat(3,:,1) * fac
   ! Quadrupoles
   do i = 1, mol%nat
      ptr%multipoles(5:9,i) = matmul(dtrafo, wfn%qpat(:,i,1)) 
   end do

! write(*,*) wfn%qpat(:,:,1) - matmul(transpose(dtrafo), ptr%multipoles(5:9,:))
! stop
   call multipole_electrostatics(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%multipoles, 2, ptr%ddx_electrostatics, ptr%ddx_error) 
   call multipole_psi(ptr%ddx%params, ptr%multipoles, 2, ptr%ddx_state%psi)

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

   !%%%%%%%%%% MONOPOLE POTENTIAL %%%%%%%%%% 
   ptr%ddx_pot = 0.0_wp
   ! Contract with the Coulomb matrix
   call gemv(ptr%jmat, ptr%ddx_state%zeta, ptr%ddx_pot(:), alpha=-1.0_wp, beta=1.0_wp, trans='t') 
   ! Scale with 0.5*feps, and get second contribution to potential
   if (self%ddx_input%ddx_model == ddx_solvation_model%lpb) then
      ptr%ddx_pot(:) = 0.5_wp * self%feps * (ptr%ddx_pot(:) + sqrt(4.0_wp*pi) * ptr%ddx_state%x_lpb(1, :, 1))
   else
      ptr%ddx_pot(:) = 0.5_wp * self%feps * (ptr%ddx_pot(:) + sqrt(4.0_wp*pi) * 1.0_wp/(ptr%ddx%params%rsph(:)**(0)) * ptr%ddx_state%xs(1, :))
   end if
   ! Add potential to overall potential for new SCF step 
   pot%vat(:,1) = pot%vat(:,1) + ptr%ddx_pot(:)

   !%%%%%%%%%% DIPOLE POTENTIAL %%%%%%%%%%
   ptr%ddx_dppot = 0.0_wp
   do k = 1, 3
      call gemv(ptr%adpmat(k,:,:), ptr%ddx_state%zeta, ptr%ddx_dppot(k,:), alpha=-1.0_wp, beta=1.0_wp, trans='t') 
      ptr%ddx_dppot(k,:) = 0.5_wp * self%feps * (ptr%ddx_dppot(k,:) + fac*4.0_wp*pi/3.0_wp * 1.0_wp/(ptr%ddx%params%rsph(:)**(1)) * ptr%ddx_state%xs(k+1, :))
   end do
   pot%vdp(:,:,1) = pot%vdp(:,:,1) + ptr%ddx_dppot(:,:)

   !%%%%%%%%%%% QUADRUPOLE POTENTIAL %%%%%%%%%%
   ptr%ddx_qppot = 0.0_wp
   ddx_qppot_trans = 0.0_wp
   trans_fac = [ sqrt(15.0_wp/(4.0_wp*pi)),              & ! m = -2,  * Theta_xy
              sqrt(15.0_wp/(4.0_wp*pi)),     & ! m = -1,  * Theta_yz   
              1.0_wp/4.0_wp*sqrt(5.0_wp/(pi)),    & ! m =  0,  * (2 Tzz - Txx - Tyy)
              sqrt(15.0_wp/(4.0_wp*pi)),     & ! m = +1,  * Theta_xz  
              1.0_wp/4.0_wp*sqrt(15.0_wp/(pi)) ]      ! m = +2,  * (Txx - Tyy)
   do k = 1, 5
      call gemv(ptr%aqpmat(k,:,:), ptr%ddx_state%zeta, ddx_qppot_trans(k,:), alpha=-1.0_wp, beta=1.0_wp, trans='t') 
      ddx_qppot_trans(k,:) = 0.5_wp * self%feps * (ddx_qppot_trans(k,:) + trans_fac(k) * 4.0_wp*pi/5.0_wp * 1.0_wp/(ptr%ddx%params%rsph(:)**2) * ptr%ddx_state%xs(k+4, :))
   end do
   ! Transform the quadrupole potential back to Cartesian form with Moore-Penrose pseudo-inverse of the transformation matrix
   do i = 1, mol%nat
      ptr%ddx_qppot(:,i) = matmul(pinv, ddx_qppot_trans(:,i))
   end do
   pot%vqp(:,:,1) = pot%vqp(:,:,1) + ptr%ddx_qppot(:,:)

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

   ptr%force = 0.0_wp

   ptr%multipoles(1, :) = wfn%qat(:, 1) / sqrt(4.0_wp*pi)
   call multipole_electrostatics(ptr%ddx%params, ptr%ddx%constants, &
      & ptr%ddx%workspace, ptr%multipoles, 1, ptr%ddx_electrostatics, ptr%ddx_error)
   call multipole_psi(ptr%ddx%params, ptr%multipoles, 1, ptr%ddx_state%psi)

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
      ptr%ddx_state, 1, ptr%multipoles, ptr%force, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   ! Calculate the gradient of the solvation energy
   ptr%force = self%feps * ptr%force 

   ! Add the gradient of the solvation energy to the total gradient
   gradient =  gradient + ptr%force

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

subroutine get_adp_matrix(xyz, ccav, adpmat)
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: ccav(:, :)
   real(wp), intent(inout) :: adpmat(:, :, :)

   integer :: ic, j
   real(wp) :: vec(3), vec2(3), d2, d

   adpmat(:, :, :) = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(ccav, xyz, adpmat) private(ic, j, vec, vec2, d2, d)
   do ic = 1, size(ccav, 2)
      do j = 1, size(xyz, 2)
         vec(:) = ccav(:, ic) - xyz(:, j)
         d2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         d = sqrt(d2)
         ! yzx vector that maps to l+m=-1,0,1 (z,y,x)
         vec2(1) = vec(2)
         vec2(2) = vec(3)
         vec2(3) = vec(1)
         adpmat(:, ic, j) = vec2(:) / (d**3)
      end do
   end do

end subroutine get_adp_matrix


subroutine get_aqp_matrix(xyz, ccav, aqpmat)
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: ccav(:, :)
   real(wp), intent(inout) :: aqpmat(:, :, :)

   integer :: ic, j
   real(wp) :: vec(3), vec2(3), d2, d, rrT(3,3), rrTcomp(6), rtrans(5)

   real(wp), parameter :: s3 = sqrt(3.0_wp)
   real(wp), parameter :: s3_4 = s3 * 0.5_wp

   real(wp), parameter :: dtrafo(5,6) = sqrt(5.0_wp/(4.0_wp*pi)) * reshape([ &
      !   m=-2     m=-1      m=0      m=+1     m=+2
      & 0.0_wp,   0.0_wp,  -0.5_wp,  0.0_wp,   s3_4,   &  ! xx
      &    s3 ,   0.0_wp,   0.0_wp,  0.0_wp,   0.0_wp, &  ! xy
      & 0.0_wp,   0.0_wp,  -0.5_wp,  0.0_wp,  -s3_4,   &  ! yy
      & 0.0_wp,   0.0_wp,   0.0_wp,     s3 ,   0.0_wp, &  ! xz
      & 0.0_wp,      s3 ,   0.0_wp,  0.0_wp,   0.0_wp, &  ! yz
      & 0.0_wp,   0.0_wp,   1.0_wp,  0.0_wp,   0.0_wp  &  ! zz
      ], shape(dtrafo))

      real(wp) :: Y2(5), err


   aqpmat(:, :, :) = 0.0_wp
   ! $omp parallel do default(none) schedule(runtime) collapse(2) &
   ! $omp shared(ccav, xyz, aqpmat) private(ic, j, vec, vec2, rrT, rrTcomp, rtrans, d2, d, c22s, c21, c20, c22c, c6)
   do ic = 1, size(ccav, 2)
      do j = 1, size(xyz, 2)
         vec(:) = ccav(:, ic) - xyz(:, j)
         d2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         d = sqrt(d2)

         rrT = matmul(reshape(vec, [3,1]), reshape(vec, [1,3]))

         ! Loose trace 
         rrT(1,1) = rrT(1,1) - d2/3.0_wp
         rrT(2,2) = rrT(2,2) - d2/3.0_wp
         rrT(3,3) = rrT(3,3) - d2/3.0_wp

         ! aqpmat(:, ic, j) = rtrans(:) / (d**5)
         ! pack Cartesian, trace-free rr^T into your order
         rrTcomp(1) = rrT(1,1)   ! xx
         rrTcomp(2) = rrT(2,1)   ! xy
         rrTcomp(3) = rrT(2,2)   ! yy
         rrTcomp(4) = rrT(3,1)   ! xz
         rrTcomp(5) = rrT(3,2)   ! yz
         rrTcomp(6) = rrT(3,3)   ! zz
               
         ! --- Voigt scale the off-diagonals before transforming ---
          rrTcomp(2) = sqrt(2.0_wp) * rrTcomp(2)   ! xy
          rrTcomp(4) = sqrt(2.0_wp) * rrTcomp(4)   ! xz
          rrTcomp(5) = sqrt(2.0_wp) * rrTcomp(5)   ! yz
               
         ! spherical kernel (5) via dtrafo
         rtrans = matmul(dtrafo, rrTcomp)
         aqpmat(:, ic, j) = 1.0_wp * rtrans(:) / (d**5)

         ! call Y2m_real_from_vec(vec, Y2)
         ! err = maxval( abs( aqpmat(:,ic,j) * d**3 / 1.0_wp - Y2(:) ) )

         ! print *, 'error: ', err

      end do
   end do

end subroutine get_aqp_matrix

pure subroutine Y2m_real_from_vec(vec, Y2)  ! vec(3) -> Y2(5) in order [-2,-1,0,+1,+2]
  use, intrinsic :: iso_fortran_env, only: wp => real64
  real(wp), intent(in)  :: vec(3)
  real(wp), intent(out) :: Y2(5)
  real(wp), parameter :: pi = acos(-1.0_wp)
  real(wp), parameter :: c15_4pi  = sqrt(15.0_wp/(4.0_wp*pi))
  real(wp), parameter :: c5_16pi  = sqrt( 5.0_wp/(16.0_wp*pi))
  real(wp), parameter :: c15_16pi = sqrt(15.0_wp/(16.0_wp*pi))
  real(wp) :: r, x, y, z, invr

  r = sqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))
  if (r == 0.0_wp) then
     Y2 = 0.0_wp
     return
  end if
  invr = 1.0_wp / r
  x = vec(1) * invr
  y = vec(2) * invr
  z = vec(3) * invr

  Y2(1) = c15_4pi  * (x*y)             ! m = -2
  Y2(2) = c15_4pi  * (y*z)             ! m = -1
  Y2(3) = c5_16pi  * (3.0_wp*z*z - 1.0_wp)  ! m =  0
  Y2(4) = c15_4pi  * (z*x)             ! m = +1
  Y2(5) = c15_16pi * (x*x - y*y)       ! m = +2
end subroutine



  subroutine inv_spd_5x5_chol(A, Ainv)
    implicit none
    real(wp), intent(in)  :: A(5,5)
    real(wp), intent(out) :: Ainv(5,5)
    real(wp) :: L(5,5), z(5), x(5), s
    integer  :: i, j, k

    L = 0.0_wp
    do j = 1, 5
       do i = j, 5
          s = A(i,j); do k = 1, j-1; s = s - L(i,k)*L(j,k); end do
          if (i == j) then
             if (s <= 0.0_wp) stop "inv_spd_5x5_chol: matrix not SPD"
             L(i,j) = sqrt(s)
          else
             L(i,j) = s / L(j,j)
          end if
       end do
    end do

    do j = 1, 5
       do i = 1, 5
          s = merge(1.0_wp, 0.0_wp, i==j)
          do k = 1, i-1; s = s - L(i,k)*z(k); end do
          z(i) = s / L(i,i)
       end do
       do i = 5, 1, -1
          s = z(i); do k = i+1, 5; s = s - L(k,i)*x(k); end do
          x(i) = s / L(i,i)
       end do
       do i = 1, 5
          Ainv(i,j) = x(i)
       end do
    end do
  end subroutine inv_spd_5x5_chol



end module tblite_solvation_ddx
