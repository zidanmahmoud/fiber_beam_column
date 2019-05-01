import numpy as np

class Structure:
    def __init__(self):
        self.elements = dict()

    def add_fiber_beam_element(self, element):
        self.elements[id] = element