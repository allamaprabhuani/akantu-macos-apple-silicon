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
# @section DESCRIPTION
#
#===============================================================================


add_mesh(hertz hertz.geo 2 1)
add_mesh(cohesive-contact cohesive-contact.geo 2 1)


register_example(contact_explicit_dynamic
  SOURCES contact_explicit_dynamic.cc
  DEPENDS hertz
  FILES_TO_COPY material.dat
  )

register_example(contact_explicit_static
  SOURCES contact_explicit_static.cc
  DEPENDS hertz
  FILES_TO_COPY material.dat
  )


register_example(cohesive_contact_explicit_dynamic
  SOURCES cohesive_contact_explicit_dynamic.cc
  DEPENDS cohesive-contact
  FILES_TO_COPY material-cohesive.dat
  )
