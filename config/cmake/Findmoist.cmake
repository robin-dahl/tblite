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
