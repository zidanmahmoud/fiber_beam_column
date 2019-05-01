
import sys
if sys.version_info < (3, 6):
    raise RuntimeError("The fiber beam-column module requires at least Python 3.6!")
from .node import Node

s =  "-------------------------------------------------------\n"
s += "+++++++++ Fiber Beam-Column Program initiated +++++++++\n"
s += "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
s += "\n"
s += "                 Author: Mahmoud Zidan                 \n"
s += "-------------------------------------------------------\n\n"
print(s)