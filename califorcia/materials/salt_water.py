from ..models import LorentzModel, DebyeModel, CombinedModel, ElectrolyteModel

# Parameters for water solvent
# Lorentz part
f1, w1, g1 = 0.77, 2.79545e+16, 2.05101e+16
lorentz = LorentzModel(f1, w1, g1)

# Debye part
epsD, epsInf, tau = 78.3, 1.8430, 7.76690e-12
debye = DebyeModel(epsD, epsInf, tau)

# Combined solvent model (Water)
# Note: CombinedModel sums contributions. Lorentz returns 1 + sum. Debye returns epsInf + ...
# We need to be careful with the constant term.
# Debye returns epsInf + ... where epsInf is already the high-frequency limit.
# Lorentz returns 1 + ...
# If we combine them: lorentz.epsilon(xi) + (debye.epsilon(xi) - 1.0)
water_model = CombinedModel([lorentz, debye])

# Electrolyte model (Salt Water)
# kappa_D: Debye screening length [1/m]. 
# For 0.1M NaCl, kappa_D is approx 10^9 m^-1. 
# We'll use a representative value or define it here.
kappa_D = 1.0e8 # [1/m] Placeholder value

epsilon = ElectrolyteModel(water_model, kappa_D)
materialclass = epsilon.materialclass
solvent_model = epsilon.solvent_model
