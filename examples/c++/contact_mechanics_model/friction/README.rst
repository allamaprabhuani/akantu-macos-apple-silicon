Friction (2D)
'''''''''''''

:Sources:

   .. collapse:: bloc_friction.cc (click to expand)

      .. literalinclude:: examples/c++/contact_mechanics_model/friction/bloc_friction.cc
         :language: c++
         :lines: 20-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/c++/contact_mechanics_model/friction/material.dat
         :language: text

:Location:

   ``examples/c++/contact_mechanics_model/`` `friction <https://gitlab.com/akantu/akantu/-/blob/master/examples/c++/contact_mechanics_model/friction>`_


In ``bloc_friction.cc``, a simple example of 2D friction is presented, using the Contact Mechanics model in a dynamic situation. 
A block is compressed against a wall, before undergoing a constrained lateral displacement from its upper surface. 
A friction coefficient of 0.3 is imposed at the interface, creating a stick-and-slip phenomenon.

The material and contact parameters are set up in the ``material.dat`` file. 
The wall and the block, lower and upper respectively, have the same properties. 
The ``contact_detector`` is used to define the numerical parameters for calculating the contact.

Finally, ``contact_resolution`` is used for the interface properties, and the penalties to be applied for numerical resolution. 
Here, a linear penalty is applied. 
``mu`` defines the coefficient of friction, and can be freely modified. 
``epsilon`` represents numerical penalties to be applied in order to maintain contact at the interface, and must be set by the user during the simulation. 
``epsilon_n`` is used to prevent the two interfaces from interpenetrating, and ``epsilon_t`` to detect friction. 
``is_master_deformable`` is set to false to ignore wall deformation, in order to reduce calculation times.

.. figure:: examples/c++/contact_mechanics_model/friction/images/bloc_friction.gif
            :align: center
            :width: 100%

            Friction of a bloc against a wall with mu = 0.3.

