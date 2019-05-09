import sys
sys.path.append('C:/Users/Fei/OngoingProjects/Force_based_fiber_beam_element/small_disp_fiberbeam/test_code_1/fiber_beam_model')
# import necessary modules
import numpy as np
import matplotlib.pyplot as plt
from fiber_beam_model import *
material = KentParkModel(1, 41.37, 1.048, -0.0096, 128.205)
ep_left = -0.005
ep_right = 0.03
#plt.plot(np.linspace(ep_left,ep_right,10), 580/60*np.linspace(ep_left,ep_right,10)+0.98, 'k--')
#plt.plot(np.linspace(ep_left,ep_right,10), 580/60*np.linspace(ep_left,ep_right,10)-0.98, 'k--')
material.initialize_material_strain_incr()
strain = 0.
stress = 0.
tangent_modulus = material.get_material_tangent_modulus()
for nt in range(10):
    material.update_change_in_material_strain_incr(-0.00154944)
    if material.check_reversal(1):
        material.update_model_parameters()
    material.save_to_last_NR_change_in_material_strain_incr(1)
    material.update_material_strain_incr()
    material.update_material_strain()
    material.update_material_stress()
    material.update_material_tangent_modulus()
    strain = np.append(strain, material.get_material_strain())
    stress = np.append(stress, material.get_material_stress())
    tangent_modulus = np.append(tangent_modulus, material.get_material_tangent_modulus())
plt.plot(strain, stress, 'b-')
#plt.plot(strain, tangent_modulus, 'go')
'''
strain = strain[-1]
stress = stress[-1]
tangent_modulus = tangent_modulus[-1]
for nt in range(1):
    material.update_change_in_material_strain_incr(-0.00256968)
    if material.check_reversal():
        material.update_model_parameters()
    material.update_material_strain_incr()
    material.update_material_strain()
    material.update_material_stress()
    material.update_material_tangent_modulus()
    strain = np.append(strain, material.get_material_strain())
    stress = np.append(stress, material.get_material_stress())
    tangent_modulus = np.append(tangent_modulus, material.get_material_tangent_modulus())
plt.plot(strain, stress, 'bo')
plt.plot(strain, tangent_modulus, 'yo')


strain = strain[-1]
stress = stress[-1]
for nt in range(22):
    material.update_change_in_material_strain_incr(0.0005)
    if material.check_reversal():
        material.update_model_parameters()
    material.update_material_strain_incr()
    material.update_material_strain()
    material.update_material_stress()
    material.update_material_tangent_modulus()
    strain = np.append(strain, material.get_material_strain())
    stress = np.append(stress, material.get_material_stress())
    tangent_modulus = material.get_material_tangent_modulus()
    plt.plot(strain, stress/60, 'b-')

strain = strain[-1]
stress = stress[-1]
for nt in range(22):
    material.update_change_in_material_strain_incr(-0.0005)
    if material.check_reversal():
        material.update_model_parameters()
    material.update_material_strain_incr()
    material.update_material_strain()
    material.update_material_stress()
    material.update_material_tangent_modulus()
    strain = np.append(strain, material.get_material_strain())
    stress = np.append(stress, material.get_material_stress())
    tangent_modulus = material.get_material_tangent_modulus()
    plt.plot(strain, stress/60, 'r-')

strain = strain[-1]
stress = stress[-1]
for nt in range(42):
    material.update_change_in_material_strain_incr(0.0005)
    if material.check_reversal():
        material.update_model_parameters()
    material.update_material_strain_incr()
    material.update_material_strain()
    material.update_material_stress()
    material.update_material_tangent_modulus()
    strain = np.append(strain, material.get_material_strain())
    stress = np.append(stress, material.get_material_stress())
    tangent_modulus = material.get_material_tangent_modulus()
    plt.plot(strain, stress/60, 'b-')

strain = strain[-1]
stress = stress[-1]
for nt in range(42):
    material.update_change_in_material_strain_incr(-0.0005)
    if material.check_reversal():
        material.update_model_parameters()
    material.update_material_strain_incr()
    material.update_material_strain()
    material.update_material_stress()
    material.update_material_tangent_modulus()
    strain = np.append(strain, material.get_material_strain())
    stress = np.append(stress, material.get_material_stress())
    tangent_modulus = material.get_material_tangent_modulus()
    plt.plot(strain, stress/60, 'r-')
'''
plt.grid(True)
plt.show()
