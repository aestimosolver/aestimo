"content = '''# Poisson Loop
'''damping:An adjustable parameter (0 < damping < 1) is typically set to 0.5 at low carrier densities. With increasing
carrier densities, a smaller value of it is needed for rapid convergence.'''
damping = 0.2    #averaging factor between iterations to smooth convergence.
Stern_damping=True#the extrapolated-convergence-factor method instead of the fixed-convergence-factor method
max_iterations=80 #maximum number of iterations.
convergence_test=1e-4 #convergence is reached when the ground state energy (meV) is stable to within this number between iterations.

# Valence band
predic_correc=True#predictor corrector method
anti_crossing_length=0.0001 # the lower lenght limit to consider anti-crossing (nm), works with old versions
amort_wave_0=1.5#ratio of half well's width for wavefunction to penetration into the the left adjacent barrier
amort_wave_1=1.5#ratio of half well's width for wavefunction to penetration into the the right adjacent barrier
strain =True # for aestimo_numpy_eh
piezo=False # directly calculationg the induced electric field,for old poisson solver, works with old versions
piezo1=True #indirectly using interface charges.
quantum_effect=True#temporary

# Output Files
parameters = True
electricfield_out = True
potential_out = True
sigma_out = True
probability_out = True
states_out = True
Drift_Diffusion_out=True

# Figures
wavefunction_scalefactor = 400 # scales wavefunctions when plotting QW diagrams
'''

with open('config.py', 'w', encoding='utf-8') as f:
    f.write(content)
print('config.py has been successfully cleaned.')"
