Examples Layout
===============

The repository examples are grouped by usage mode:

.. code-block:: text

   examples/
   ├── idl/
   │   ├── RenderExampleMW.pro
   │   ├── RenderExampleEUV.pro
   │   ├── InterpolateEBTELexample.pro
   │   └── compile_local_idl
   └── python/
       ├── cli/
       │   ├── RenderExampleMW.py
       │   └── RenderExampleEUV.py
       └── sdk/
           ├── sdk_render_mw.py
           └── sdk_render_euv.py

Python CLI examples
-------------------

Use these for quick rendering runs and to mirror the documented command-line
workflow:

.. code-block:: bash

   PYTHONPATH=src python examples/python/cli/RenderExampleMW.py --help
   PYTHONPATH=src python examples/python/cli/RenderExampleEUV.py --help

Python SDK examples
-------------------

Use these as templates for embedding rendering into your own application code
with typed SDK options/results:

.. code-block:: bash

   PYTHONPATH=src python examples/python/sdk/sdk_render_mw.py --help
   PYTHONPATH=src python examples/python/sdk/sdk_render_euv.py --help

IDL examples
------------

Compile local IDL routines for development/testing:

.. code-block:: idl

   @/path/to/gximagecomputing/examples/idl/compile_local_idl
