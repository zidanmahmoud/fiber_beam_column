from abc import ABC, abstractclassmethod


class Material(ABC):
    @abstractclassmethod
    def get_material_tangent_modulus(self):
        pass
