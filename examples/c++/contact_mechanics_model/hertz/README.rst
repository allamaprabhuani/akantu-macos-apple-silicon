Hertzian contact (2D)
'''''''''''''''''''''

:Sources:

   .. collapse:: contact_explicit_dynamic.cc (click to expand)

      .. literalinclude:: examples/c++/contact_mechanics_model/hertz/contact_explicit_dynamic.cc
         :language: c++
         :lines: 20-

   .. collapse:: contact_explicit_static.cc (click to expand)

      .. literalinclude:: examples/c++/contact_mechanics_model/hertz/contact_explicit_static.cc
         :language: c++
         :lines: 20-

   .. collapse:: material.dat (click to expand)

      .. literalinclude:: examples/c++/contact_mechanics_model/hertz/material.dat
         :language: text

:Location:

   ``examples/c++/contact_mechanics_model/`` `hertz <https://gitlab.com/akantu/akantu/-/blob/master/examples/c++/contact_mechanics_model/hertz>`_


`contact_explicit_static` and `contact_explicit_dynamic` are solving a 2D Hertz contact patch test using the ``CouplerSolidContact``.
The two examples follow what is described extensively in section :ref:`sect-cmm-coupling-with-smm`. The only main difference between `contact_explicit_static` and `contact_explicit_dynamic` is the solver used::
    
    // for contact_explicit_static
    coupler.initFull(_analysis_method = _static);  
    // for contact_explicit_dynamic
    coupler.initFull(_analysis_method = _explicit_lumped_mass);  

The ``material.dat`` file contain a ``contact_detector`` and a ``contact_resolution penalty_linear`` section as explains in section :ref:`sect-cmm-contact-detection`.

.. figure:: examples/c++/contact_mechanics_model/hertz/images/hertz.svg
            :align: center
            :width: 60%

.. figure:: examples/c++/contact_mechanics_model/hertz/images/hertz.png
            :align: center
            :width: 60%
