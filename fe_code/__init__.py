"""
Module fe_code
==============
"""

import sys

from .structure import Structure
from .node import Node

if sys.version_info < (3, 6):
    raise RuntimeError("The fiber beam-column module requires at least Python 3.6!")

MSG = "-------------------------------------------------------\n"
MSG += "+++++++++ Fiber Beam-Column Program initiated +++++++++\n"
MSG += "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
MSG += "\n"
MSG += "                 Author: Mahmoud Zidan                 \n"
MSG += "-------------------------------------------------------\n\n"
print(MSG)
