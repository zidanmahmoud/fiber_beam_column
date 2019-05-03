from abc import ABC, abstractclassmethod

class Element(ABC):

    @abstractclassmethod
    def calculate_global_stiffness_matrix(self):
        pass

    @abstractclassmethod
    def dofs(self):
        pass
