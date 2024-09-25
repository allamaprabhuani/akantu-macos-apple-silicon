#!/usr/bin/env python3
# pylint: disable=missing-module-docstring
# pylint: disable=missing-function-docstring
__copyright__ = (
    "Copyright (©) 2022-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import numpy as np
import akantu as aka

spatial_dimension = 2


class LinearCohesive(aka.MaterialCohesive):
    """Material Cohesive Linear."""

    def __init__(self, model, _id):
        """Construct the material and register the parameters to the parser."""
        super().__init__(model, _id)
        super().registerParamReal(
            "G_c", aka._pat_readable | aka._pat_parsable, "Fracture energy"
        )
        super().registerParamReal(
            "beta", aka._pat_readable | aka._pat_parsable, ""
        )
        self.registerInternalReal("delta_max", 1)
        self.beta = 0.
        self.sigma_c = 0.
        self.delta_c = 0.
        self.G_c = 0.

    def initMaterial(self):
        """Initialize the parameters for the material."""
        super().initMaterial()
        self.sigma_c = self.getReal("sigma_c")
        self.G_c = self.getReal("G_c")
        self.beta = self.getReal("beta")
        self.delta_c = 2 * self.G_c / self.sigma_c

    def checkInsertion(self, _check_only):
        """Check if need to insert a cohesive element."""
        model = self.getModel()
        facets = self.getFacetFilter()
        inserter = model.getElementInserter()

        for type_facet in facets.elementTypes(dim=(spatial_dimension - 1)):
            facet_filter = facets(type_facet)
            nb_facet = facet_filter.shape[0]
            if nb_facet == 0:
                continue

            fe_engine = model.getFEEngine("FacetsFEEngine")

            facets_check = inserter.getCheckFacets(type_facet)
            insertion = inserter.getInsertionFacets(type_facet)

            nb_quad_facet = fe_engine.getNbIntegrationPoints(type_facet)
            normals = fe_engine.getNormalsOnIntegrationPoints(type_facet)
            facets_stresses = model.getStressOnFacets(type_facet)

            tangents = model.getTangents(type_facet)

            for facet in facet_filter:
                if not facets_check[facet]:
                    continue

                ref_stress = 0

                for q in range(nb_quad_facet):
                    current_quad = facet * nb_quad_facet + q
                    normal = normals[current_quad, :].ravel()
                    tangent = tangents[current_quad, :].ravel()

                    quad_stresses = facets_stresses[current_quad, :].reshape(
                        2, spatial_dimension, spatial_dimension)

                    stress_1 = quad_stresses[0]
                    stress_2 = quad_stresses[1]

                    avg_stress = stress_1 + stress_2 / 2.0
                    traction = avg_stress.dot(normal)

                    n = traction.dot(normal)
                    n = max(0, n)
                    t = traction.dot(tangent)

                    ref_stress = max(
                        ref_stress, np.sqrt(n * n + t * t / (self.beta**2))
                    )

                if ref_stress > self.sigma_c:
                    insertion[facet] = True

    # constitutive law
    def computeTraction(self, el_type, ghost_type):
        """Compute the traction for a given opening."""
        normals = self.getNormals(el_type, ghost_type)
        openings = self.getOpening(el_type, ghost_type)
        tractions = self.getTraction(el_type, ghost_type)

        delta_max = self.getInternalReal("delta_max")(el_type)

        for el in range(normals.shape[0]):
            normal = normals[el].ravel()
            opening = openings[el].ravel()

            delta_n = opening.dot(normal) * normal
            delta_s = opening - delta_n

            delta = self.beta * np.linalg.norm(delta_s) ** 2 + \
                np.linalg.norm(delta_n) ** 2

            delta_max[el] = max(delta, delta_max[el])

            tractions[el, :] = (
                (delta * delta_s + delta_n)
                * self.sigma_c
                / delta
                * (1 - delta / self.delta_c)
            )


def allocator(_dim, _unused, model, _id):  # NOQA
    """Register the material to the material factory."""
    return LinearCohesive(model, _id)


mat_factory = aka.MaterialFactory.getInstance()
mat_factory.registerAllocator("local_cohesive", allocator)

# -------------------------------------------------------------------------
# Initialization
# -------------------------------------------------------------------------
aka.parseInput("local_material.dat")
mesh = aka.Mesh(spatial_dimension)
mesh.read("plate.msh")

model = aka.SolidMechanicsModelCohesive(mesh)
model.initFull(_analysis_method=aka._static, _is_extrinsic=True)
model.initNewSolver(aka._explicit_lumped_mass)

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

# -------------------------------------------------------------------------
# Boundary conditions
# -------------------------------------------------------------------------
model.applyBC(aka.FixedValue(0.0, aka._x), "XBlocked")
model.applyBC(aka.FixedValue(0.0, aka._y), "YBlocked")

trac = np.zeros(spatial_dimension)
traction = 0.095
trac[int(aka._y)] = traction

model.getExternalForce()[:] = 0
model.applyBC(aka.FromTraction(trac), "Traction")

print("Solve for traction ", traction)
solver = model.getNonLinearSolver("static")
solver.set("max_iterations", 100)
solver.set("threshold", 1e-10)
solver.set("convergence_type", aka.SolveConvergenceCriteria.residual)

model.solveStep("static")
model.dump()
model.dump("cohesive elements")

model.setTimeStep(model.getStableTimeStep() * 0.1)

maxsteps = 100

for i in range(0, maxsteps):
    print("{0}/{1}".format(i, maxsteps))
    model.checkCohesiveStress()
    model.solveStep("explicit_lumped")
    if i % 10 == 0:
        model.dump()
        model.dump("cohesive elements")
