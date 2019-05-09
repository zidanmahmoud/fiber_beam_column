import sys
sys.path.append('C:/Users/Fei/OngoingProjects/Force_based_fiber_beam_element/small_disp_fiberbeam/test_code_1/fiber_beam_model')
# import necessary modules
import numpy as np
import matplotlib.pyplot as plt
from fiber_beam_model import *
material = KentParkModel(1, 41.37, 1.048, -0.0096, 128.205)

material.initialize_material_strain_incr()
strain = 0.
stress = 0.
tangent_modulus = -material.get_material_tangent_modulus()
for nt in range(4):
    material.update_change_in_material_strain_incr(-0.002)
    if material.check_reversal():
        material.update_model_parameters()
    material.update_material_strain_incr()
    material.update_material_strain()
    material.update_material_stress()
    material.update_material_tangent_modulus()
    strain = np.append(strain, material.get_material_strain())
    stress = np.append(stress, material.get_material_stress())
    tangent_modulus = np.append(tangent_modulus, material.get_material_tangent_modulus())
plt.plot(-strain, -stress, 'bo')
'''
strain = strain[-1]
stress = stress[-1]
tangent_modulus = tangent_modulus[-1]
for nt in range(10):
    material.update_change_in_material_strain_incr(0.0004)
    if material.check_reversal():
        material.update_model_parameters()
    material.update_material_strain_incr()
    material.update_material_strain()
    material.update_material_stress()
    material.update_material_tangent_modulus()
    strain = np.append(strain, material.get_material_strain())
    stress = np.append(stress, material.get_material_stress())
    tangent_modulus = np.append(tangent_modulus, material.get_material_tangent_modulus())
plt.plot(-strain, -stress, 'r-')

strain = strain[-1]
stress = stress[-1]
tangent_modulus = tangent_modulus[-1]
for nt in range(20):
    material.update_change_in_material_strain_incr(-0.00033)
    if material.check_reversal():
        material.update_model_parameters()
    material.update_material_strain_incr()
    material.update_material_strain()
    material.update_material_stress()
    material.update_material_tangent_modulus()
    strain = np.append(strain, material.get_material_strain())
    stress = np.append(stress, material.get_material_stress())
    tangent_modulus = np.append(tangent_modulus, material.get_material_tangent_modulus())
'''
#plt.plot(-strain, -stress, 'g-')
plt.grid(True)
plt.show()
