"""
Module contains the abstract class Material
"""

from abc import ABC, abstractclassmethod


class UniaxialIncrementalMaterial(ABC):
    """
    Material abstract class to be used in fiber-beam-column element
    """

    @abstractclassmethod
    def update_strain(self, fiber_strain):
        pass

    @abstractclassmethod
    def _set_trial_state(self, new_strain):
        pass

    @abstractclassmethod
    def finalize_load_step(self):
        pass
