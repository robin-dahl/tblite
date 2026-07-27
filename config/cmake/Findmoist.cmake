# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later

set(_lib "moist")
set(_pkg "MOIST")
set(_url "https://github.com/lukaswittmann/moist")
set(_rev "v0.6.0-alpha.1")

if(NOT DEFINED "${_pkg}_FIND_METHOD")
  if(DEFINED "${PROJECT_NAME}-dependency-method")
    set("${_pkg}_FIND_METHOD" "${${PROJECT_NAME}-dependency-method}")
  else()
    set("${_pkg}_FIND_METHOD" "cmake" "pkgconf" "subproject" "fetch")
  endif()
  set("_${_pkg}_FIND_METHOD")
endif()

# Keep the dependency-only subproject lean.
set(MOIST_API OFF CACHE BOOL "" FORCE)
set(MOIST_TESTS OFF CACHE BOOL "" FORCE)

include("${CMAKE_CURRENT_LIST_DIR}/tblite-utils.cmake")

# Avoid preparing source-only dependencies when an installed moist CMake
# package is already available.
if("cmake" IN_LIST ${_pkg}_FIND_METHOD AND NOT TARGET moist::moist)
  find_package(moist CONFIG QUIET)
endif()

# moist's CMake build expects Meson wrap dependencies to have populated these
# source directories.  A CMake FetchContent checkout does not process wrap
# files, so provide the missing targets before configuring moist.
if("fetch" IN_LIST ${_pkg}_FIND_METHOD AND NOT TARGET moist::moist)
  include(FetchContent)

  if(NOT TARGET test-drive)
    set(TESTDRIVE_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
      test-drive
      GIT_REPOSITORY "https://github.com/fortran-lang/test-drive"
      GIT_TAG "v0.6.1"
    )
    FetchContent_MakeAvailable(test-drive)
  endif()

  if(NOT TARGET fclap)
    FetchContent_Declare(
      fclap
      GIT_REPOSITORY "https://github.com/chselz/fclap"
      GIT_TAG "8f499d50596aea1726fdd806e93400d1547a6f1e"
    )
    FetchContent_MakeAvailable(fclap)
  endif()

  if(NOT TARGET jonquil)
    FetchContent_Declare(
      jonquil
      GIT_REPOSITORY "https://github.com/lukaswittmann/jonquil"
      GIT_TAG "fdb94458e536266850d63cde04b36bf21872a48a"
    )
    FetchContent_MakeAvailable(jonquil)
  endif()
endif()

tblite_find_package("${_lib}" "${${_pkg}_FIND_METHOD}" "${_url}" "${_rev}")

if(TARGET moist AND NOT TARGET moist::moist)
  add_library(moist::moist ALIAS moist)
endif()

if(DEFINED "_${_pkg}_FIND_METHOD")
  unset("${_pkg}_FIND_METHOD")
  unset("_${_pkg}_FIND_METHOD")
endif()
unset(_lib)
unset(_pkg)
unset(_url)
unset(_rev)
