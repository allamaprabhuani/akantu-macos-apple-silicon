#===============================================================================
# Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# This file is part of Akantu
# 
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
# 
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.30)
  cmake_policy(SET CMP0167 NEW)
  set(Boost_NO_BOOST_CMAKE OFF CACHE BOOL "" FORCE)
else()
  set(Boost_NO_BOOST_CMAKE ON CACHE BOOL "" FORCE)
endif()

package_declare(Boost EXTERNAL
  NOT_OPTIONAL
  DESCRIPTION "Package handling boost components"
  EXTRA_PACKAGE_OPTIONS TARGET Boost::headers
)


mask_package_options(Boost)
mark_as_advanced(AKANTU_USE_BOOST)

package_on_enabled_script(Boost
  "if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER \"4.8\")
  set(_boost_version \${Boost_MAJOR_VERSION}.\${Boost_MINOR_VERSION})
  if(_boost_version VERSION_LESS 1.58 AND _boost_version VERSION_GREATER 1.53)
    package_set_compile_flags(Boost CXX -DBOOST_SPIRIT_USE_PHOENIX_V3)
  else()
    package_unset_compile_flags(Boost CXX)
  endif()
endif()
")

package_set_package_system_dependency(Boost deb-src libboost-dev)

# only build dependency
# package_set_package_system_dependency(Boost deb libboost)
# package_set_package_system_dependency(Boost rmp boost)
