# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later

set(_lib "cint")
set(_pkg "CINT")
set(_url "https://github.com/sunqm/libcint.git")
set(_rev "v6.1.3")

if(NOT DEFINED "${_pkg}_FIND_METHOD")
  if(DEFINED "${PROJECT_NAME}-dependency-method")
    set("${_pkg}_FIND_METHOD" "${${PROJECT_NAME}-dependency-method}")
  else()
    set("${_pkg}_FIND_METHOD" "cmake" "pkgconf" "subproject" "fetch")
  endif()
  set("_${_pkg}_FIND_METHOD")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/tblite-utils.cmake")

tblite_find_package("${_lib}" "${${_pkg}_FIND_METHOD}" "${_url}" "${_rev}")

if(DEFINED "_${_pkg}_FIND_METHOD")
  unset("${_pkg}_FIND_METHOD")
  unset("_${_pkg}_FIND_METHOD")
endif()
unset(_lib)
unset(_pkg)
unset(_url)
unset(_rev)
