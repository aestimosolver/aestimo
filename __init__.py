from . import aestimo
from . import aestimo_eh
from . import config
from . import database
from . import intersubband_optical_transitions

import os

localpath = lambda fname: os.path.abspath(
    os.path.join(os.path.dirname(__file__), fname)
)

# load module docstring
__doc__ = open(localpath("README.md")).read()

__version__ = aestimo.__version__
