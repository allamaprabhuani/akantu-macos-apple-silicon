multi-material (2D)
'''''''''''''''''''

:Sources:

   .. collapse:: plate.py (click to expand)

      .. literalinclude:: examples/python/solid_mechanics_cohesive_model/cohesive/plate.py
         :language: python
         :lines: 9-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/python/solid_mechanics_cohesive_model/cohesive/material.dat
         :language: text

:Location:

   ``examples/python/solid_mechanics_cohesive_model/`` `multi-material <https://gitlab.com/akantu/akantu/-/blob/master/examples/python/solid_mechanics_cohesive_model/multi-material/>`_


This example is the same as example as the :ref:`ex-python-cohesive` example
The main difference is that in this example there are multiple cohesive laws. The way they are selected is by using a `aka.MaterialCohesiveRulesSelector` ::

  cohesive_selector = aka.MaterialCohesiveRulesSelector(model, {
    ("Top", "Bottom"): "tough",
    ("Top", "Top"): "soft",
    ("Bottom", "Bottom"): "soft",
  })

  model.setMaterialSelector(cohesive_selector)
