
from fe_code.node import Node
from fe_code.fiber_beam import FiberBeam
from fe_code.dof import Dof

# print(Dof(1, "x"))

# for dof in FiberBeam(Node(1, 0, 0, 0), Node(2, 0, 0, 0)).dofs:
#     print(dof)
print(FiberBeam(Node(1, 0, 0, 0), Node(2, 0, 0, 0)).dofs)
