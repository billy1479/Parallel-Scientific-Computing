#!/usr/bin/env python3
import subprocess
import os
import argparse
from tqdm import tqdm

def generate_and_test_timestep(dt, executable="step-0", n=4, min_mass=1.0, max_mass=100.0, 
                              final_time=5.0, snapshots=0.1):
    """Generate a configuration and test stability with given timestep"""
    # Use Python script to generate configuration
    gen_cmd = [
        "python3", "./paste.txt",
        "--final-time", str(final_time),
        "--snapshots", str(snapshots),
        "--dt", str(dt),
        "--executable-name", executable,
        "--min-mass", str(min_mass),
        "--max-mass", str(max_mass),
        "--dim", "3",
        "--N", str(n),
        "--scenario", "random-grid"
    ]
    
    print(f"Generating configuration with dt={dt}")
    try:
        result = subprocess.run(gen_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error generating configuration: {result.stderr}")
            return False, "Generation failed"
    except Exception as e:
        print(f"Exception during generation: {e}")
        return False, str(e)
    
    # Fix the shell script - add shebang and make it executable
    script_path = f"{executable}.sh"
    try:
        with open(script_path, 'r') as f:
            content = f.read()
        
        with open(script_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write(content)
        
        os.chmod(script_path, 0o755)
    except Exception as e:
        print(f"Error fixing shell script: {e}")
        return False, str(e)
    
    # Run the simulation
    print(f"Running simulation with dt={dt}")
    try:
        result = subprocess.run([f"./{script_path}"], capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            return False, "Simulation crashed"
        return True, "Simulation completed successfully"
    except subprocess.TimeoutExpired:
        return False, "Simulation timed out"
    except Exception as e:
        return False, str(e)

def main():
    parser = argparse.ArgumentParser(description='Test time step stability')
    parser.add_argument('--min-dt', type=float, default=1e-7, help='Minimum time step')
    parser.add_argument('--max-dt', type=float, default=1e-4, help='Maximum time step')
    parser.add_argument('--steps', type=int, default=10, help='Number of steps to test')
    parser.add_argument('--executable', default='step-0', help='Executable name')
    parser.add_argument('--n', type=int, default=4, help='Grid size (N×N×N particles)')
    args = parser.parse_args()
    
    # Create logarithmically spaced time steps
    import numpy as np
    time_steps = np.logspace(np.log10(args.min_dt), np.log10(args.max_dt), args.steps)
    
    results = []
    last_stable = None
    
    print(f"Testing {args.steps} time steps from {args.min_dt} to {args.max_dt}")
    
    for dt in time_steps:
        dt_str = f"{dt:.10f}"
        success, message = generate_and_test_timestep(dt_str, args.executable, args.n)
        results.append((dt, success, message))
        
        print(f"Time step {dt_str}: {'Stable' if success else 'Unstable'} - {message}")
        
        if success:
            last_stable = dt
        elif last_stable is not None:
            # We found the transition point between stable and unstable
            print(f"\nTransition found! Last stable time step: {last_stable:.10f}")
            break
    
    # Report results
    if last_stable is not None:
        print(f"\nLargest stable time step: {last_stable:.10f}")
    else:
        stable_steps = [dt for dt, success, _ in results if success]
        if stable_steps:
            print(f"\nLargest stable time step: {max(stable_steps):.10f}")
        else:
            print("\nNo stable time steps found")
    
    print("\nAll test results:")
    for dt, success, message in results:
        print(f"dt={dt:.10f}: {'Stable' if success else 'Unstable'} - {message}")

if __name__ == "__main__":
    main()