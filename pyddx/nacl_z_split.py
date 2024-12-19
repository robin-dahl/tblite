import pyddx
import numpy as np

tobohr = 1 / 0.52917721092
#charges = np.array([-0.79715541, 0.79715541])
#charges = np.array([-0.82282044, 0.82282044])
charges = np.array([-1.0000, 0.0000])

#rvdw = np.array([
#    3.8739385554719026,  3.4015070243167926])
rvdw = tobohr* np.array([
    3.00,  1.00])

centres = tobohr * np.array([
    [ 0.0, 0.0, 0.0],
    [ 0.0, 0.0, 3.0]]).T

model = pyddx.Model("cosmo", centres, rvdw, solvent_epsilon=78.359999999999999, n_lebedev=302,lmax=6)

# Compute solute contributions (here just charges)
solute_multipoles = charges.reshape(1, -1) / np.sqrt(4 * np.pi)
solute_field = model.multipole_electrostatics(solute_multipoles)
solute_psi = model.multipole_psi(solute_multipoles)

#print(solute_field)

# Solve the problem
state = pyddx.State(model, solute_psi, solute_field["phi"])
state.fill_guess()
state.solve()
state.fill_guess_adjoint()
state.solve_adjoint()

# Compute energy and forces
energy = 0.5 * np.sum(state.x * solute_psi)
print(f'xs: {state.x}')
print(f'energy: {energy}')
force_solvation = state.solvation_force_terms(solute_field)
print(f'force_solvation: {force_solvation}')
force_solute = state.multipole_force_terms(solute_multipoles);
print(f'force_solute: {force_solute}')
print(f'tot: {force_solute + force_solvation}')
