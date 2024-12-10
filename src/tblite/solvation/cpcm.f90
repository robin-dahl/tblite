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

!> @file tblite/solvation/cpcm.f90
!> Provides a polarizable continuum model

!> Implicit solvation model based on a polarizable dielectric continuum
module tblite_solvation_cpcm
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use tblite_blas, only : dot, gemv
   use tblite_container_cache, only : container_cache
   use tblite_mesh_lebedev, only : grid_size, get_angular_grid, list_bisection
   use tblite_scf_info, only : scf_info, atom_resolved
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   ! use tblite_solvation_cpcm_dd
   use tblite_solvation_data, only : get_vdw_rad_cosmo
   use tblite_solvation_type, only : solvation_type


   use iso_fortran_env, only : output_unit

   use tblite_disp_d4, only: get_eeq_charges


   use ddx, only: ddx_type, ddx_error_type, check_error, ddinit, ddx_state_type, allocate_state
   use ddx, only: fill_guess, fill_guess_adjoint, ddrun, solvation_force_terms, solve_adjoint
   use ddx_core, only: ddx_electrostatics_type, allocate_electrostatics
   use ddx_cosmo, only: cosmo_solve, cosmo_solve_adjoint

   use ddx_multipolar_solutes, only: multipole_electrostatics

   implicit none
   private

   public :: cpcm_solvation, new_cpcm, cpcm_input, cpcm_cache


   !> Input for CPCM solvation
   type :: cpcm_input
      !> Dielectric constant
      real(wp) :: dielectric_const
      !> Number of grid points for each atom
      integer :: nang = grid_size(6)
      !> Scaling of van-der-Waals radii
      real(wp) :: rscale = 1.0_wp
      !> Accuracy for iterative solver
      real(wp) :: conv = 1.0e-8_wp
      !> Regularization parameter
      real(wp) :: eta = 2.0_wp
      !> Maximum angular momentum of basis functions
      integer :: lmax = 6
      !> Van-der-Waals radii for all species
      real(wp), allocatable :: rvdw(:)
   end type cpcm_input


   !> Definition of polarizable continuum model
   type, extends(solvation_type) :: cpcm_solvation
      !> Actual domain decomposition calculator
      ! type(domain_decomposition) :: dd

      !> Dielectric function
      real(wp) :: keps
      !> Dielctric constant
      real(wp) :: dielectric_const
      !> Van-der-Waal radii for all atoms
      real(wp), allocatable :: rvdw(:)


      
      real(wp) :: ddx_tol = 1.0e-11_wp

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
   end type cpcm_solvation

   !> Provide constructor for CPCM solvation
   interface cpcm_solvation
      module procedure :: create_cpcm
   end interface


   !> Restart data for CPCM calculation
   type :: cpcm_cache
      !> Actual domain decomposition calculator
      ! type(domain_decomposition) :: dd
      type(ddx_type) :: xdd 

      type(ddx_state_type) :: ddx_state
      
      type(ddx_error_type) :: ddx_error

      type(ddx_electrostatics_type) :: ddx_electrostatics

      ! !> Electrostatic potential phi(ncav)
      ! ! ddx_state%phi_cav
      ! real(wp), allocatable :: phi(:) !phi_cav ?

      ! !> Psi vector psi(nylm, n)
      ! ! ddx_state%psi
      ! real(wp), allocatable :: psi(:, :)

      ! !> CPCM solution sigma(nylm, n)
      ! ! ddx_state%xs
      ! real(wp), allocatable :: sigma(:, :)

      ! !> CPCM adjoint solution s(nylm, n)
      ! ! ddx_state%s
      ! real(wp), allocatable :: s(:, :)

      !> Interaction matrix with surface charges jmat(ncav, nat)
      real(wp), allocatable :: jmat(:, :)

      real(wp) :: esolv

      real(wp), allocatable :: zeta(:)

      real(wp), allocatable :: force(:, :)

      real(wp), allocatable :: multipoles(:,:)

      integer :: do_force = 0
   

   end type cpcm_cache


   !> Identifier for container
   character(len=*), parameter :: label = "polarizable continuum model"

   real(wp), parameter :: alpha_alpb = 0.571412_wp
   integer, parameter :: ndiis = 25


contains


!> Create new electric field container
subroutine new_cpcm(self, mol, input, error)
   !> Instance of the solvation model
   type(cpcm_solvation), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input for CPCM solvation
   type(cpcm_input), intent(in) :: input
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, izp
   ! real(wp), allocatable :: ang_grid(:, :), ang_weight(:)

   self%label = label

   ! Radii for all atoms
   allocate(self%rvdw(mol%nat))
   if (allocated(input%rvdw)) then
      self%rvdw(:) = input%rscale*input%rvdw(mol%id)
   else
      do iat = 1, mol%nat
         izp = mol%num(mol%id(iat))
         self%rvdw(iat) = input%rscale*get_vdw_rad_cosmo(izp)
      end do
   end if

   write(*,*) 'rvdw = ', self%rvdw

   ! Epssilon
   self%dielectric_const = input%dielectric_const

   write(*,*) 'dielectric_const = ', self%dielectric_const

   ! keps
   self%keps = -0.5_wp * (1.0_wp/self%dielectric_const - 1.0_wp) / (1.0_wp + alpha_alpb)

   write(*,*) 'keps = ', self%keps

end subroutine new_cpcm


!> Type constructor for CPCM splvation
function create_cpcm(mol, input) result(self)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input for CPCM solvation
   type(cpcm_input), intent(in) :: input
   !> Instance of the solvation model
   type(cpcm_solvation) :: self

   type(error_type), allocatable :: error

   call new_cpcm(self, mol, input, error)
end function create_cpcm


!> Update cache from container
subroutine update(self, mol, cache)
   !> Instance of the solvation model
   class(cpcm_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Number of grid points for each atom
   integer :: nang = grid_size(6)
   !> Scaling of van-der-Waals radii
   real(wp) :: rscale = 1.0_wp
   !> Maximum angular momentum of basis functions
   integer :: lmax = 6

   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   type(cpcm_cache), pointer :: ptr
   call taint(cache, ptr)

   ptr%do_force = 1
   ptr%ddx_electrostatics%do_phi = 1

   ! ddinit
   ! 1-cosmo, 2-pcm
   call ddinit(1, mol%nat, mol%xyz, self%rvdw, self%dielectric_const, ptr%xdd, &
      & ptr%ddx_error, ngrid=302, lmax=6, nproc=8, force=ptr%do_force)
   call check_error(ptr%ddx_error)
   write(*,*) '----------CHECKPOINT: ddinit done----------'

   print *, 'ptr%xdd%params%force', ptr%xdd%params%force
   allocate(ptr%force(3, ptr%xdd%params%nsph))

   call allocate_state(ptr%xdd%params, ptr%xdd%constants,  &
      ptr%ddx_state, ptr%ddx_error)
   call check_error(ptr%ddx_error)
   write(*,*) '----------CHECKPOINT: allocate state done----------'

   allocate(ptr%multipoles(1, ptr%xdd%params%nsph))
   write(*,*) '----------CHECKPOINT: allocate multipoles done----------'
 
   allocate(ptr%jmat(ptr%xdd%constants%ncav, mol%nat))

   call get_coulomb_matrix(mol%xyz, ptr%xdd%constants%ccav, ptr%jmat)
   write(*,*) '----------CHECKPOINT: got J matrix----------'

   allocate(ptr%zeta(ptr%xdd%params%ngrid * ptr%xdd%params%nsph))
   write(*,*) '----------CHECKPOINT: got zeta----------'

   ! Activate fast multipole method
   ptr%xdd%params%fmm = .true.


end subroutine update


!!! not needed anymore..?
!> Get electric field energy
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the solvation model
   class(cpcm_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Solvation free energy
   real(wp), intent(inout) :: energies(:)

   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   type(cpcm_cache), pointer :: ptr
   call taint(cache, ptr)

   energies(:) = energies + self%keps * sum(ptr%ddx_state%xs * ptr%ddx_state%psi , 1)
   write(*,*) 'energy = ', sum(energies(:)) / 627.5095_wp

end subroutine get_energy


!> Get electric field potential
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the solvation model
   class(cpcm_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   type(cpcm_cache), pointer :: ptr
   call taint(cache, ptr)

   !> Calculate the multipole moments
   ! With maximum angular momentum of the multipoles = 1:
   ! We have charges, so we need to convert them to monopoles.
   ! For charges the conversion is simply a scaling by 1/sqrt(4*pi).
   ptr%multipoles(1, :) = wfn%qat(:, 1) / ( sqrt( 4 * pi ))

   !> Multipole electrostatics
   ! Compute the electrostatic properties for a multipolar distributions of 
   ! arbitrary order, provided that they are given in real spherical harmonics.
   call multipole_electrostatics(ptr%xdd%params, ptr%xdd%constants, ptr%xdd%workspace, &
      ptr%multipoles, 0, ptr%ddx_electrostatics, ptr%ddx_error)
   call check_error(ptr%ddx_error)

   !> Calculate Representation of the solute density in spherical harmonics
   !! (\f$ \Psi \f$). It is used as RHS for the adjoint linear system.
   !! Dimension (nbasis, nsph).
   call get_psi(wfn%qat(:, 1), ptr%ddx_state%psi)

   !> Calculate all solvation terms
   call ddrun(ptr%xdd, ptr%ddx_state, ptr%ddx_electrostatics, ptr%ddx_state%psi, self%ddx_tol, &
      & ptr%esolv, ptr%ddx_error, ptr%force)
   call check_error(ptr%ddx_error)

   write(*,*) ''
   write(*,*) 'qat = ', wfn%qat(:, 1)
   write(*,*) 'esolv = ', ptr%esolv

   ! We need the weights w, the switching function ui, the sp harm expansion v (v + w -> vw), and the solution to the COSMO eq S
   ! allocate(ptr%zeta(ptr%xdd%constants%ncav))
   call get_zeta(ptr, self%keps)

   ! Contract with the Coulomb matrix
   call gemv(ptr%jmat, ptr%zeta, pot%vat(:, 1), alpha=-1.0_wp, beta=1.0_wp, trans='t')

   !> Caclulate the potential
   pot%vat(:, 1) = pot%vat(:, 1) + (self%keps * sqrt(4*pi)) * ptr%ddx_state%xs(1, :)

end subroutine get_potential


!> Get electric field gradient
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the solvation model
   class(cpcm_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the solvation free energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the solvation free energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: ii, iat, ig
   real(wp), allocatable :: gx(:, :), zeta(:), ef(:, :)
   type(cpcm_cache), pointer :: ptr

   call view(cache, ptr)

   write(*,*) '----------CHECKPOINT: get_gradient----------'

   ! allocate(gx(3, mol%nat), zeta(ptr%dd%ncav), ef(3, max(mol%nat, ptr%dd%ncav)))
   allocate(gx(3, mol%nat), zeta(ptr%xdd%constants%ncav), ef(3, max(mol%nat, ptr%xdd%constants%ncav)))

   ! call solve_cosmo_adjoint(ptr%dd, ptr%ddx_state%psi, ptr%ddx_state%s, .true., &
   !    & accuracy=ptr%dd%conv*1e-3_wp)

   call solve_adjoint(ptr%xdd%params, ptr%xdd%constants, ptr%xdd%workspace, ptr%ddx_state, self%ddx_tol, &
      & ptr%ddx_error)

   ! reset Phi
   ! call get_phi(wfn%qat(:, 1), ptr%jmat, ptr%ddx_state%phi_cav)
   call get_phi(wfn%qat(:, 1), ptr%jmat, ptr%ddx_electrostatics%phi_cav)

   ! now call the routine that computes the ddcosmo specific contributions to the forces.
   !call get_deriv(ptr%dd, self%keps, ptr%ddx_state%phi_cav, ptr%ddx_state%xs, ptr%ddx_state%s, gx)

   ! form the "zeta" intermediate
   !call get_zeta(ptr%dd, self%keps, ptr%ddx_state%s, zeta)
   !call get_zeta(ptr, self%keps)

   ! 1. solute's electric field at the cav points times zeta:
   !    compute the electric field
   ! call efld(mol%nat, wfn%qat(:, 1), ptr%dd%xyz, ptr%dd%ncav, ptr%dd%ccav, ef)
   !call efld(mol%nat, wfn%qat(:, 1), mol%xyz, ptr%xdd%constants%ncav, ptr%xdd%constants%ccav, ef)


   call solvation_force_terms(ptr%xdd%params, ptr%xdd%constants, ptr%xdd%workspace, ptr%ddx_state, &
      & ptr%ddx_electrostatics, ptr%force, ptr%ddx_error)

   write(*,*) 'force = ', ptr%force


   ! contract it with the zeta intermediate
   ! ii = 0
   ! do iat = 1, mol%nat
   !    do ig = 1, ptr%xdd%constants%ngrid
   !       if (self%xdd%constants%ui(ig, iat) > 0.0_wp) then
   !          ii = ii + 1
   !          gx(:, iat) = gx(:, iat) + ptr%zeta(ii)*ef(:, ii)
   !       end if
   !    end do
   ! end do

   ! 2. "zeta's" electric field at the nuclei times the charges.
   !    compute the "electric field"
   ! call efld(ptr%xxd%constants%ncav, ptr%zeta, ptr%xxd%constants%ccav, mol%nat, mol%xyz, ef)

   ! ! contract it with the solute's charges.
   ! do iat = 1, mol%nat
   !    gx(:, iat) = gx(:, iat) + ef(:, iat)*wfn%qat(iat, 1)
   ! end do

   ! gradient(:, :) = gradient(:, :) + gx
end subroutine get_gradient


!> Return dependency on density
pure function variable_info(self) result(info)
   !> Instance of the solvation model
   class(cpcm_solvation), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=atom_resolved)
end function variable_info


subroutine taint(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(cpcm_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(cpcm_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

subroutine view(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(cpcm_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(cpcm_cache)
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
   real(wp), parameter :: fac = sqrt(4*pi)

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


!> Computes the electric field produced by the sources src (nsrc point charges
!  with coordinates csrc) at the ntrg target points ctrg:
subroutine efld(nsrc, src, csrc, ntrg, ctrg, ef)
   integer, intent(in) :: nsrc, ntrg
   real(wp), intent(in) :: src(:)
   real(wp), intent(in) :: csrc(:, :)
   real(wp), intent(in) :: ctrg(:, :)
   real(wp), intent(inout) :: ef(:, :)

   integer :: i, j
   real(wp) :: vec(3), r2, rr, r3, f
   real(wp), parameter :: zero=0.0_wp

   ef(:, :) = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp reduction(+:ef) shared(ntrg, nsrc, ctrg, csrc, src) &
   !$omp private(j, i, f, vec, r2, rr, r3)
   do j = 1, ntrg
      do i = 1, nsrc
         vec(:) = ctrg(:, j) - csrc(:, i)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         rr = sqrt(r2)
         r3 = r2*rr
         f = src(i)/r3
         ef(:, j) = ef(:, j) + f*vec
      end do
   end do

end subroutine efld


subroutine get_zeta(self, keps)
   ! type(domain_decomposition), intent(in) :: self
   type(cpcm_cache), intent(inout) :: self
   real(wp), intent(in) :: keps
   !real(wp), intent(inout) :: zeta(:) ! [self%ncav]

   integer :: its, iat, ii
   real, dimension(74, 2) :: ui_temp

   ! write(*,*) 'ngrid = ', self%xdd%params%ngrid yes
   ! write(*,*) 'keps = ', keps yes
   ! write(*,*) 'ui = ', 
   ! call write_2d_matrix(self%xdd%constants%ui, name='ui', step=5) NO
   ! write(*,*) 'vgrid = ', self%xdd%constants%vgrid yes
   ! write(*,*) 'wgrid = ', self%xdd%constants%wgrid yes 
   !  write(*,*) 's zeta = ', self%ddx_state%s 

   ! => S, XS, AND UI DIFFERENT TO LEGACY CODE


   ii = 0
   do iat = 1, self%xdd%params%nsph
      do its = 1, self%xdd%params%ngrid
         if (self%xdd%constants%ui(its, iat) > 0.0_wp) then
            ii = ii + 1
            self%zeta(ii) = keps * self%xdd%constants%wgrid(its) * self%xdd%constants%ui(its, iat) &
                & * dot_product(self%xdd%constants%vgrid(:, its), self%ddx_state%s(:, iat))
         end if
      end do
   end do



 end subroutine get_zeta



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

   if (present(name)) write(iunit,'(/,"vector printed:",1x,a)') name

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

   if (present(name)) write(iunit,'(/,"matrix printed:",1x,a)') name

   do i = 1, d2, istep
       l = min(i+istep-1,d2)
       write(iunit,'(/,6x)',advance='no')
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



end module tblite_solvation_cpcm
