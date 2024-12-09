import pyddx
import numpy as np
 
#tobohr = 1 / 0.52917721092
tobohr = 1.88973
charges = np.array([
     -0.50000000000000000,      -0.50000000000000000
])
rvdw = tobohr * np.array([
    3.2503289343471575,        2.4566439620065728
])
centres = tobohr * np.array([
    [ 0.00000,  0.00000, 0.00000],
    [ 0.00000,  0.00000, 1.00000],
]).T
 
print(pyddx.banner())
 
model = pyddx.Model("cpcm", charges, centres, rvdw, solvent_epsilon=78.3553)
 
# Compute solute contributions (here just charges)
solute_multipoles = charges.reshape(1, -1) / np.sqrt(4 * np.pi)
solute_field = model.multipole_electrostatics(solute_multipoles)
solute_psi = model.multipole_psi(solute_multipoles)
 
# Solve the problem
state = pyddx.State(model, solute_field["phi"], solute_psi)
state.fill_guess()
state.solve()
state.fill_guess_adjoint()
state.solve_adjoint()
 
# Show results
energy = 0.5 * np.sum(state.x * solute_psi)
force = state.solvation_force_terms()
print(energy)
print(force)
