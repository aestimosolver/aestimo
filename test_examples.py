import os
import sys
import subprocess
import time
from pathlib import Path

def run_example(file_path):
    print(f"Testing {file_path.name}...")
    start_time = time.time()
    try:
        # Run with -d for figures but we don't want to show them in a loop if they block
        # Actually most examples call run_aestimo and might have plt.show()
        # We should use Agg backend to avoid blocking
        env = os.environ.copy()
        env["MPLBACKEND"] = "Agg"
        env["PYTHONPATH"] = os.getcwd()
        
        result = subprocess.run(
            [sys.executable, str(file_path)],
            capture_output=True,
            text=True,
            timeout=120, # 2 minutes timeout per example
            env=env
        )
        
        elapsed = time.time() - start_time
        if result.returncode == 0:
            print(f"PASSED {file_path.name} ({elapsed:.2f}s)")
            return True, None
        else:
            print(f"FAILED {file_path.name} ({elapsed:.2f}s)")
            return False, result.stderr
    except subprocess.TimeoutExpired:
        print(f"TIMEOUT {file_path.name}")
        return False, "Timeout"
    except Exception as e:
        print(f"ERROR {file_path.name}: {e}")
        return False, str(e)

def main():
    examples_dir = Path("examples")
    python_files = sorted([f for f in examples_dir.glob("*.py") if f.name != "samples.py"])
    
    total = len(python_files)
    passed = 0
    failed = []
    
    print(f"Starting tests for {total} examples...\n", flush=True)
    
    for i, file_path in enumerate(python_files):
        print(f"[{i+1}/{total}] ", end="", flush=True)
        success, error = run_example(file_path)
        if success:
            passed += 1
        else:
            failed.append((file_path.name, error))
            
    print(f"\n--- Test Results ---", flush=True)
    print(f"Total:  {total}", flush=True)
    print(f"Passed: {passed}", flush=True)
    print(f"Failed: {len(failed)}", flush=True)
    
    if failed:
        print("\nFailed Examples Details:")
        for name, error in failed:
            print(f"\n--- {name} ---")
            print(error)
            print("-" * (len(name) + 8))
            
if __name__ == "__main__":
    main()
