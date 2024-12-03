#===============================================================================
# Copyright (©) 2021-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


add_mesh(compression_mesh compression.geo DIM 2 ORDER 1)

register_example(compression
  SCRIPT compression.py
  PYTHON
  FILES_TO_COPY compression.dat
  DEPENDS compression_mesh)

register_example(detection_explicit.py
  SCRIPT detection-explicit.py
  PYTHON
  FILES_TO_COPY detection-explicit.dat
    detection-explicit.msh
  )
