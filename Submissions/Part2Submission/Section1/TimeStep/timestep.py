import subprocess
import re
import numpy as np
import os
import argparse
from tqdm import tqdm

def run_simulation(executable, config, time_step, duration, plot_interval=0.0):
    """Run the n-body simulation with given parameters and return the output."""
    # Command format: ./executable plot-time final-time dt objects
    command = [executable, str(plot_interval), str(duration), str(time_step)] + config.split()
    
    try:
        # Capture both stdout and stderr
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        return result.returncode == 0, result.stdout + result.stderr
    except Exception as e:
        return False, str(e)

def analyze_stability(output, max_velocity_threshold=1000.0, min_distance_threshold=1e-6):
    """Analyze simulation output to determine if it was stable."""
    # Check for errors
    if "error" in output.lower():
        return False, "Error in simulation"
    
    # Extract maximum velocity
    max_v_pattern = r"v_max=\s*([0-9.e+-]+)"
    max_v_matches = re.findall(max_v_pattern, output)
    
    if max_v_matches:
        try:
            max_v = float(max_v_matches[-1])  # Get last reported max velocity
            if max_v > max_velocity_threshold:
                return False, f"Maximum velocity too high: {max_v}"
        except ValueError:
            pass  # Invalid float value, continue checking
    
    # Extract minimum distance
    min_dx_pattern = r"dx_min=\s*([0-9.e+-]+)"
    min_dx_matches = re.findall(min_dx_pattern, output)
    
    if min_dx_matches:
        try:
            min_dx = float(min_dx_matches[-1])  # Get last reported min distance
            if min_dx < min_distance_threshold:
                return False, f"Bodies too close: {min_dx}"
        except ValueError:
            pass  # Invalid float value
    
    # Check if simulation completed
    if "Number of remaining objects" not in output:
        return False, "Simulation did not complete"
        
    return True, "Simulation stable"

def find_largest_stable_time_step(executable, config, test_range, duration, tolerance=1e-6):
    """Find largest stable time step within the given range."""
    time_steps = np.array(test_range)
    stability_results = []
    
    print(f"Testing {len(time_steps)} time steps from {min(time_steps)} to {max(time_steps)}...")
    
    for dt in tqdm(time_steps):
        success, output = run_simulation(executable, config, dt, duration)
        is_stable, reason = analyze_stability(output)
        stability_results.append((dt, is_stable, reason))
        
    # Find largest stable time step
    stable_steps = [dt for dt, is_stable, _ in stability_results if is_stable]
    
    if not stable_steps:
        print("No stable time steps found in the given range.")
        print("Stability results:")
        for dt, is_stable, reason in stability_results:
            print(f"dt={dt}: {'Stable' if is_stable else 'Unstable'} - {reason}")
        return None
    
    largest_stable = max(stable_steps)
    
    # Print detailed results
    print("\nStability results:")
    for dt, is_stable, reason in sorted(stability_results):
        print(f"dt={dt}: {'Stable' if is_stable else 'Unstable'} - {reason}")
    
    return largest_stable

def binary_search_stable_time_step(executable, config, lower_bound, upper_bound, duration, tolerance=1e-6, max_iterations=20):
    """Use binary search to find the largest stable time step."""
    print(f"Finding largest stable time step between {lower_bound} and {upper_bound}...")
    
    iteration = 0
    
    # Ensure lower bound is stable
    _, output = run_simulation(executable, config, lower_bound, duration)
    is_stable, reason = analyze_stability(output)
    if not is_stable:
        print(f"Lower bound {lower_bound} is not stable: {reason}")
        return None
    
    # Ensure upper bound is unstable
    _, output = run_simulation(executable, config, upper_bound, duration)
    is_stable, reason = analyze_stability(output)
    if is_stable:
        print(f"Upper bound {upper_bound} is stable. Increase the upper bound.")
        return upper_bound
    
    while (upper_bound - lower_bound) > tolerance and iteration < max_iterations:
        mid_point = (lower_bound + upper_bound) / 2
        print(f"Testing dt={mid_point}")
        
        _, output = run_simulation(executable, config, mid_point, duration)
        is_stable, reason = analyze_stability(output)
        
        if is_stable:
            print(f"dt={mid_point} is stable: {reason}")
            lower_bound = mid_point
        else:
            print(f"dt={mid_point} is unstable: {reason}")
            upper_bound = mid_point
            
        iteration += 1
    
    return lower_bound

def main():
    parser = argparse.ArgumentParser(description='Find the largest stable time step for N-body simulation.')
    parser.add_argument('--executable', default='./nbody', help='Path to the N-body simulation executable')
    parser.add_argument('--config', required=True, help='Configuration string for the simulation')
    parser.add_argument('--duration', type=float, default=10.0, help='Duration to run each simulation')
    parser.add_argument('--search-type', choices=['range', 'binary'], default='binary', 
                        help='Search method to use')
    parser.add_argument('--min-dt', type=float, default=0.0001, help='Minimum time step to test')
    parser.add_argument('--max-dt', type=float, default=0.1, help='Maximum time step to test')
    parser.add_argument('--steps', type=int, default=20, help='Number of steps to test (for range search)')
    parser.add_argument('--tolerance', type=float, default=1e-6, help='Tolerance for binary search')
    
    args = parser.parse_args()
    
    # Example configurations from the code comments
    example_configs = {
        "one-body": "0.0 0.0 0.0  1.0 0.0 0.0  1.0",
        "two-bodies-orbit": "0.0 0.0 0.0  0.0 0.0 0.0 1.0  1.0 0.0 0.0  0.0 1.0 0.0 0.1",
        "two-bodies-spiral": "0.0 0.0 0.0  0.0 0.5 0.0 1.0  1.0 0.0 0.0  0.0 -0.5 0.0 1.0",
        "three-bodies": "3.0 0.0 0.0  0.0 1.0 0.0  0.4  0.0 0.0 0.0  0.0 0.0 0.0  0.2  2.0 0.0 0.0  0.0 0.0 0.0  1.0"
    }
    
    # If config is a key in example_configs, use that config
    config = example_configs.get(args.config, args.config)
    
    if args.search_type == 'range':
        # Create logarithmic range of time steps
        test_range = np.logspace(np.log10(args.min_dt), np.log10(args.max_dt), args.steps)
        result = find_largest_stable_time_step(args.executable, config, test_range, args.duration)
    else:  # binary search
        result = binary_search_stable_time_step(args.executable, config, args.min_dt, args.max_dt, 
                                             args.duration, args.tolerance)
    
    if result is not None:
        print(f"\nLargest stable time step found: {result}")
    else:
        print("\nFailed to find a stable time step.")

if __name__ == "__main__":
    main()