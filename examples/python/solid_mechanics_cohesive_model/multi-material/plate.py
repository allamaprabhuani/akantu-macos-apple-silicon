#!/usr/bin/env python3
__copyright__ = (
    "Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import akantu as aka
import numpy as np


def initialize_model(model):
    """Initialize the model."""
    cohesive_selector = aka.MaterialCohesiveRulesSelector(model, {
        ("Top", "Bottom"): "tough",
        ("Top", "Top"): "soft",
        ("Bottom", "Bottom"): "soft",
    })

    model.setMaterialSelector(cohesive_selector)

    model.initFull(_analysis_method=aka._static, _is_extrinsic=True)
    model.initNewSolver(aka._explicit_lumped_mass)


def set_dumpers(model):
    """Set the dumpers."""
    model.setBaseName("plate")
    model.addDumpFieldVector("displacement")
    model.addDumpFieldVector("external_force")
    model.addDumpField("strain")
    model.addDumpField("stress")
    model.addDumpField("blocked_dofs")

    model.setBaseNameToDumper("cohesive elements", "cohesive")
    model.addDumpFieldVectorToDumper("cohesive elements", "displacement")
    model.addDumpFieldToDumper("cohesive elements", "damage")
    model.addDumpFieldVectorToDumper("cohesive elements", "tractions")
    model.addDumpFieldVectorToDumper("cohesive elements", "opening")


def do_static_step(model):
    """Do a static step."""
    solver = model.getNonLinearSolver("static")
    solver.set("max_iterations", 100)
    solver.set("threshold", 1e-10)
    solver.set("convergence_type", aka.SolveConvergenceCriteria.residual)

    model.solveStep("static")
    model.dump()
    model.dump("cohesive elements")


def solve(material_file, mesh_file, traction):
    """Define the main problem."""
    aka.parseInput(material_file)
    spatial_dimension = 2

    # -------------------------------------------------------------------------
    # Initialization
    # -------------------------------------------------------------------------
    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    model = aka.SolidMechanicsModelCohesive(mesh)

    initialize_model(model)
    set_dumpers(model)

    # -------------------------------------------------------------------------
    # Boundary conditions
    # -------------------------------------------------------------------------
    model.applyBC(aka.FixedValue(0.0, aka._x), "XBlocked")
    model.applyBC(aka.FixedValue(0.0, aka._y), "YBlocked")

    trac = np.zeros(spatial_dimension)
    trac[int(aka._y)] = traction

    print("Solve for traction ", traction)

    model.getExternalForce()[:] = 0
    model.applyBC(aka.FromTraction(trac), "Traction")

    do_static_step(model)

    model.setTimeStep(model.getStableTimeStep() * 0.1)

    maxsteps = 1000

    for i in range(0, maxsteps):
        print("{0}/{1}".format(i, maxsteps))
        model.checkCohesiveStress()
        model.solveStep("explicit_lumped")
        if i % 10 == 0:
            model.dump()
            model.dump("cohesive elements")


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    mesh_file = "plate.msh"
    material_file = "material.dat"

    traction = 0.095
    solve(material_file, mesh_file, traction)
