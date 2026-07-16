! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later

!> Construction of selectable integral backends.
module tblite_integral_handler
   use mctc_env, only : error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_integral_type, only : integral_type
   use tblite_integral_native, only : native_integral_type
#if TBLITE_HAS_LIBCINT
   use tblite_integral_libcint, only : libcint_integral_type
#endif
   implicit none
   private

   public :: new_integral_handler

contains

subroutine new_integral_handler(handler, mol, basis, error, use_libcint)
   class(integral_type), allocatable, intent(out) :: handler
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: basis
   type(error_type), allocatable, intent(out) :: error
   logical, intent(in), optional :: use_libcint
   logical :: libcint

   libcint = .false.
   if (present(use_libcint)) libcint = use_libcint
   if (libcint) then
      if (any(mol%periodic)) then
         call fatal_error(error, "libcint integral backend does not yet support periodic systems")
         return
      end if
#if TBLITE_HAS_LIBCINT
      allocate(libcint_integral_type :: handler)
#else
      call fatal_error(error, "libcint integral backend is not available in this build")
      return
#endif
   else
      allocate(native_integral_type :: handler)
   end if
   call handler%initialize_integral(mol, basis)
end subroutine new_integral_handler

end module tblite_integral_handler
