import pyddx
import numpy as np

tobohr = 1 / 0.52917721092

charges = np.array([
   #   -0.50000000000000000,      -0.50000000000000000
   -0.79715541, 0.79715541
])
rvdw =  np.array([
    3.8739385554719026,        3.4015070243167926
])
centres = tobohr *  np.array([
    [ 0.00000,  0.00000, 0.00000],
    [ 0.00000,  0.00000, 2.50000],
]).T

model = pyddx.Model("cosmo", centres, rvdw, solvent_epsilon=78.359999999999999,
n_lebedev=302,lmax=6)

# Compute solute contributions (here just charges)
solute_multipoles = charges.reshape(1, -1) / np.sqrt(4 * np.pi)
solute_field = model.multipole_electrostatics(solute_multipoles)
solute_psi = model.multipole_psi(solute_multipoles)

print(f'sol psi : {solute_psi}')
#print(f'sol field : {solute_field}')


# Solve the problem
state = pyddx.State(model, solute_psi, solute_field["phi"])
state.fill_guess()
state.solve()
state.fill_guess_adjoint()
state.solve_adjoint()

# Compute energy and forces
energy = 0.5 * np.sum(state.x * solute_psi)
force = state.solvation_force_terms(solute_field)
force += state.multipole_force_terms(solute_multipoles);

print(f'multipole force term {state.multipole_force_terms(solute_multipoles)}')

# Show results
print('energy',energy)
print('force',force)
#print('x',state.x)
#print('s',state.s)

# dimension of x
#print('dimension of x',state.x.shape)

