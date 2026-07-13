import sys
import json
import time
sys.path.insert(0, '.')
import aestimo

def test_run(polarization, bc_right, solver):
    with open('examples/untitled_project.json', 'r') as f:
        config = json.load(f)
    
    config['enable_polarization'] = polarization
    config['bc_right'] = bc_right
    config['vmin'] = 0.0
    config['vmax'] = 0.1
    config['vstep'] = 0.1
    config['G_optical'] = 0.0
    config['comp_scheme'] = solver
    config['__file__'] = f'test_out_{polarization}_{bc_right}_{solver}'
    
    # Need to convert layer list
    mat_list = []
    for l in config['layers']:
        mat_list.append([float(l['thickness']), l['material'], float(l['mole']), float(l.get('mole_y', 0.0)), float(l['doping']), l['doping_type'], l['type']])
    config['material'] = mat_list
    config['Fapplied'] = float(config.get('field', 0.0))
    config['Each_Step'] = float(config.get('vstep', 0.05))
    config['device_area'] = float(config.get('area', 1e-4))
    if 'tat_field' in config: config['tat_field'] = float(config['tat_field'])
    if 'mat_sys' in config: config['mat_type'] = config['mat_sys']
    
    config['maxgridpoints'] = 200000
    config['max_iterations'] = 200  # fast fail
    
    print(f"\n--- Testing Pol={polarization}, bc_right={bc_right}, solver={solver} ---")
    start = time.time()
    import traceback
    try:
        aestimo.run_aestimo(config, drawFigures=False, show=False)
        print(f"SUCCESS in {time.time()-start:.2f}s")
    except Exception as e:
        print(f"FAILED in {time.time()-start:.2f}s: {e}")
        traceback.print_exc()

if __name__ == "__main__":
    test_run(True, 0.6, 7)
