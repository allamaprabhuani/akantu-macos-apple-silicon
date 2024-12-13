#!/usr/bin/env python3
__copyright__ = (
    "Copyright (©) 2021-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import akantu as aka


MAX_STEPS = 20000
max_displacement = 1e-2

damping_interval = 10
damping_ratio = 0.9

spatial_dimension = 2


aka.parseInput('compression.dat')

mesh = aka.Mesh(spatial_dimension)
mesh.read('compression.msh')

coupler = aka.CouplerSolidContact(mesh)

solid = coupler.getSolidMechanicsModel()
contact = coupler.getContactMechanicsModel()

material_selector = aka.MeshDataMaterialSelectorString("physical_names", solid)
solid.setMaterialSelector(material_selector)

coupler.initFull(_analysis_method=aka._explicit_lumped_mass)

surface_selector = aka.PhysicalSurfaceSelector(mesh)
detector = contact.getContactDetector()
detector.setSurfaceSelector(surface_selector)

solid.applyBC(aka.FixedValue(0.0, aka._x), "sides")

time_step = solid.getStableTimeStep()
time_step *= 0.1
coupler.setTimeStep(time_step)

coupler.setBaseName("compression")
coupler.addDumpFieldVector("displacement")
coupler.addDumpFieldVector("contact_force")
coupler.addDumpFieldVector("external_force")
coupler.addDumpFieldVector("internal_force")
coupler.addDumpField("gaps")
coupler.addDumpField("areas")
coupler.addDumpField("blocked_dofs")
coupler.addDumpField("grad_u")
coupler.addDumpField("stress")

coupler.dump()

velocity = solid.getVelocity()

increment = max_displacement / MAX_STEPS


for s in range(0, MAX_STEPS):

    print("Step : ", s)

    solid.applyBC(aka.IncrementValue(-increment, aka._y), "loading")
    solid.applyBC(aka.IncrementValue(increment, aka._y), "fixed")

    coupler.solveStep()

    if s % damping_interval == 0:
        velocity *= damping_ratio

    if s % 100 == 0:
        coupler.dump()
