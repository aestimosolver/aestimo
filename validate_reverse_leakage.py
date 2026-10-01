import numpy as np
import aestimo
import json

def run_test_point(voltage, tat_field=1e8, grid_step=0.2):
    print(f"\nTesting V = {voltage}V, tat_field = {tat_field:.1e}, grid_step = {grid_step}nm")
    
    # Load base config from the example
    with open('examples/untitled_project.json', 'r') as f:
        config = json.load(f)
    
    # Apply necessary overrides for the test
    config['vmin'] = voltage
    config['vmax'] = voltage
    config['vstep'] = 0.0
    config['tat_field'] = tat_field
    config['grid_step'] = grid_step
    config['damping'] = 0.2
    config['comp_scheme'] = 9
    config['G_optical'] = 0.0
    config['enable_polarization'] = True
    config['__file__'] = f'test_{voltage}'

    # Convert layers to material matrix as required by aestimo
    mat_list = []
    for l in config.get('layers', []):
        mat_list.append([float(l['thickness']), l['material'], float(l.get('mole', 0.0)), 
                         float(l.get('mole_y', 0.0)), float(l['doping']), l['doping_type'], l['type']])
    config['material'] = mat_list
    config.pop('layers', None)

    try:
        result = aestimo.run_aestimo(config, drawFigures=False, show=False)
        # Extract the current from the result object
        # Based on aestimo.py, av_curr is stored in the result object
        curr = result.av_curr[0] if hasattr(result, 'av_curr') else 'N/A'
        print(f"Resulting Current: {curr} A/m^2")
        return curr
    except Exception as e:
        print(f"Simulation failed: {e}")
        return None

if __name__ == '__main__':
    test_points = [0.0, -0.2, -0.8]
    final_results = {}
    
    for v in test_points:
        res = run_test_point(v)
        final_results[v] = res
    
    print('\n' + '='*30)
    print('FINAL VALIDATION SUMMARY')
    print('='*30)
    for v, c in final_results.items():
        print(f'V = {v:>4}V  =>  Current = {c}')
