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


add_mesh(plate_mesh_mm_py plate.geo DIM 2 ORDER 1)

register_example(multi_material_py
  SCRIPT plate.py
  PYTHON
  FILES_TO_COPY material.dat
  DEPENDS plate_mesh_mm_py)
