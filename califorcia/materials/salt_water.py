from ..models import ElectrolyteModel
from .water_11 import epsilon as water_11_model

# Electrolyte model (Salt Water)
# kappa_D: Debye screening length [1/m]. 
kappa_D = 1.0e8 # [1/m] Placeholder value, updated in simulations

epsilon = ElectrolyteModel(water_11_model, kappa_D)
materialclass = epsilon.materialclass
solvent_model = epsilon.solvent_model
