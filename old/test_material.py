import sys
sys.path.append('C:/Users/Fei/OngoingProjects/Force_based_fiber_beam_element/small_disp_fiberbeam/test_code_1/fiber_beam_model')
# import necessary modules
import numpy as np
from fiber_beam_model import *
lin_iso_hardening = IsotropicLinearHardenging('linear_isotropic_hardening', 1., 1., 1.)
lin_iso_hardening_1 = IsotropicLinearHardenging('linear_isotropic_hardening', 1., 1., 1.)
list_of_materials = []
list_of_materials.append(lin_iso_hardening)
list_of_materials.append(lin_iso_hardening_1)
print(len(list_of_materials))
lin_iso_hardening.update_change_in_material_strain_incr(2)
lin_iso_hardening.initialize_material_strain_incr()
lin_iso_hardening.update_material_strain_incr()
lin_iso_hardening.update_material_strain()
lin_iso_hardening.initialize_material_PEEQ_incr()
lin_iso_hardening.update_material_PEEQ_incr()

lin_iso_hardening.update_material_stress()
lin_iso_hardening.update_material_tangent_modulus()
lin_iso_hardening.update_material_PEEQ()

lin_iso_hardening.save_to_material_strain_last_loadstep()
lin_iso_hardening.save_to_material_PEEQ_last_loadstep()
lin_iso_hardening.save_to_material_stress_lasttime_j()

print(lin_iso_hardening.get_material_strain())
print(lin_iso_hardening.get_material_stress())
print(lin_iso_hardening.get_material_tangent_modulus())
print(lin_iso_hardening.get_material_PEEQ())

lin_iso_hardening.update_change_in_material_strain_incr(-2)
lin_iso_hardening.initialize_material_strain_incr()
lin_iso_hardening.update_material_strain_incr()
lin_iso_hardening.update_material_strain()
lin_iso_hardening.initialize_material_PEEQ_incr()
lin_iso_hardening.update_material_PEEQ_incr()

lin_iso_hardening.update_material_stress()
lin_iso_hardening.update_material_tangent_modulus()
lin_iso_hardening.update_material_PEEQ()

lin_iso_hardening.save_to_material_strain_last_loadstep()
lin_iso_hardening.save_to_material_PEEQ_last_loadstep()
lin_iso_hardening.save_to_material_stress_lasttime_j()

print(lin_iso_hardening.get_material_strain())
print(lin_iso_hardening.get_material_stress())
print(lin_iso_hardening.get_material_tangent_modulus())
print(lin_iso_hardening.get_material_PEEQ())
