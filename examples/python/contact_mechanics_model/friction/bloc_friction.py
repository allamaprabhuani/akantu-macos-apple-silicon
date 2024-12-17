#!/usr/bin/env python3
__copyright__ = (
    "Copyright (©) 2021-2024 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"

import akantu as aka

aka.parseInput('material.dat')

mesh = aka.Mesh(2)
mesh.read('bloc.msh')

coupler = aka.CouplerSolidContact(mesh)

solid = coupler.getSolidMechanicsModel()
contact = coupler.getContactMechanicsModel()

material_selector = aka.MeshDataMaterialSelectorString("physical_names", solid)
solid.setMaterialSelector(material_selector)

coupler.initFull(_analysis_method=aka._explicit_lumped_mass)

surface_selector = aka.PhysicalSurfaceSelector(mesh)
detector = contact.getContactDetector()
detector.setSurfaceSelector(surface_selector)


time_step = solid.getStableTimeStep()
time_step *= 0.1
coupler.setTimeStep(time_step)
print(f"Time step: {time_step}")

coupler.setBaseName("bloc_dumped")
coupler.addDumpFieldVector("displacement")
coupler.addDumpFieldVector("velocity")
coupler.addDumpFieldVector("normals")
coupler.addDumpFieldVector("tangents")
coupler.addDumpFieldVector("contact_force")
coupler.addDumpFieldVector("external_force")
coupler.addDumpFieldVector("internal_force")
coupler.addDumpField("areas")
coupler.addDumpField("stress")
coupler.addDumpField("blocked_dofs")

solid.applyBC(aka.FixedValue(0.0, aka._x), "XFixed")
solid.applyBC(aka.FixedValue(0.0, aka._y), "YFixed")
solid.applyBC(aka.FixedValue(0.0, aka._x), "loading")
solid.applyBC(aka.FixedValue(0.0, aka._y), "loading")

velocity = solid.getVelocity()
gaps = contact.getGaps()

coupler.dump()

for s in range(0, 10000):
    solid.applyBC(aka.IncrementValue(-1.0 / 10000, aka._y), "loading")

    coupler.solveStep()

    for i, (gap, vel) in enumerate(zip(gaps, velocity)):
        if gap > 0:
            velocity[i] *= 0.99

    if s % 100 == 0:
        coupler.dump()
        print(f"Step {s}\t", end="\r", flush=True)

print("Compression done !")

for s in range(0, 10000):
    solid.applyBC(aka.IncrementValue(2.0 / 10000, aka._x), "loading")

    coupler.solveStep()

    for i, (gap, vel) in enumerate(zip(gaps, velocity)):
        if gap > 0:
            velocity[i] *= 0.99

    if s % 100 == 0:
        coupler.dump()
        print(f"Step {s}\t", end="\r", flush=True)

print("Sliding done !")

for s in range(0, 10000):
    coupler.solveStep()

    for i, (gap, vel) in enumerate(zip(gaps, velocity)):
        if gap > 0:
            velocity[i] *= 0.99

    if s % 100 == 0:
        coupler.dump()
        print(f"Step {s}\t", end="\r", flush=True)
