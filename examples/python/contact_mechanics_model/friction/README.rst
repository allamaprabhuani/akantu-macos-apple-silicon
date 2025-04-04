Friction (2D)
'''''''''''''
:Sources:

   .. collapse:: bloc_friction.py (click to expand)

      .. literalinclude:: examples/python/contact_mechanics_model/friction/bloc_friction.py
         :language: python
         :lines: 9-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/python/contact_mechanics_model/friction/material.dat
         :language: text

:Location:

   ``examples/python/contact_mechanics_model/`` `friction <https://gitlab.com/akantu/akantu/-/blob/master/examples/python/contact_mechanics_model/friction>`_

In `friction.py` it is shown how to couple the Solid Mechanics Model with the Contact Mechanics Model, applying the friction law.
The example simulate the contact between two blocks with a friction coefficient of 0.3.
It is similar to the 'bloc_friction.cc' example in C++.

.. figure:: examples/c++/contact_mechanics_model/friction/images/bloc_friction.gif
            :align: center
            :width: 100%

            Friction of a bloc against a wall with mu = 0.3.

