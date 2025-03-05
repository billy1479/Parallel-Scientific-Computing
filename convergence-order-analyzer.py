#!/usr/bin/env python3
import subprocess
import numpy as np
import os
import re
import argparse
import matplotlib.pyplot as plt
from stable_timestep_calculator import calculate_stable_timestep, read_particle_data
import math

def extract_final_positions(output):
    """Extract final positions of particles from simulation output."""
    # Find the position of the first particle in the final step
    match = re.search(r"Position of first remaining object: ([-+]?[0-9]*\.?[0-9]+), ([-+]?[0-9]*\.?[0-9]+), ([-+]?[0-9]*\.?[0-9]+)", output)
    if match:
        x = float(match.group(1))
        y = float(match.group(2))
        z = float(match.group(3))
        return np.array([x, y, z])
    return None

def generate_config(dt, input_file, output_file):
    """Generate a configuration with the given time step."""
    # Read the original file
    with open(input_file, 'r') as f:
        content = f.read()
    
    # Extract parameters
    params = content.split()
    
    # Replace time step (assuming it's the third parameter)
    if len(params) > 3:
        params[3] = str(dt)
    
    # Write new configuration
    with open(output_file, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write(" ".join(params))
    
    os.chmod(output_file, 0o755)
    return True

def run_simulation(script_path, timeout=300):
    """Run the simulation and return output."""
    try:
        result = subprocess.run([script_path], capture_output=True, text=True, timeout=timeout)
        return result.returncode == 0, result.stdout + result.stderr
    except subprocess.TimeoutExpired:
        return False, "Simulation timed out"
    except Exception as e:
        return False, str(e)

def measure_convergence_order(base_timestep, num_steps=5, refinement_factor=2, 
                            input_file="step-0.sh", timeout=300):
    """Measure convergence order by running simulations with decreasing time steps."""
    # Generate time steps
    time_steps = [base_timestep / (refinement_factor**i) for i in range(num_steps)]
    
    print(f"Testing {num_steps} time steps from {time_steps[0]} to {time_steps[-1]}")
    
    results = []
    positions = []
    
    for i, dt in enumerate(time_steps):
        print(f"\nRunning simulation {i+1}/{num_steps} with dt={dt}")
        
        # Generate configuration
        config_file = f"step-0-dt-{i}.sh"
        generate_config(dt, input_file, config_file)
        
        # Run simulation
        success, output = run_simulation(f"./{config_file}", timeout)
        
        if not success:
            print(f"Error running simulation with dt={dt}: {output[:100]}...")
            continue
        
        # Extract final position
        position = extract_final_positions(output)
        if position is not None:
            positions.append(position)
            results.append(dt)
            print(f"Final position: {position}")
        else:
            print(f"Failed to extract final position for dt={dt}")
    
    # Calculate errors using smallest time step as reference
    if len(positions) < 2:
        print("Not enough successful simulations to calculate convergence order")
        return None
    
    reference_position = positions[-1]  # Smallest time step
    errors = [np.linalg.norm(pos - reference_position) for pos in positions]
    
    # Calculate convergence order
    orders = []
    for i in range(len(errors)-1):
        if errors[i] > 0 and errors[i+1] > 0 and results[i] > 0 and results[i+1] > 0:
            order = math.log(errors[i]/errors[i+1]) / math.log(results[i]/results[i+1])
            orders.append(order)
    
    avg_order = sum(orders) / len(orders) if orders else None
    
    # Plot results
    plt.figure(figsize=(10, 6))
    
    # Error vs time step
    plt.subplot(1, 2, 1)
    plt.loglog(results, errors, 'o-')
    plt.xlabel('Time Step (dt)')
    plt.ylabel('Error')
    plt.title('Error vs Time Step')
    plt.grid(True)
    
    # Error vs time step with fit line
    plt.subplot(1, 2, 2)
    log_dt = np.log10(results)
    log_error = np.log10(errors)
    
    # Linear fit
    if len(log_dt) > 1:
        fit = np.polyfit(log_dt, log_error, 1)
        plt.plot(log_dt, log_error, 'o', label='Data')
        plt.plot(log_dt, np.polyval(fit, log_dt), '-', label=f'Fit: Slope = {fit[0]:.2f}')
        plt.xlabel('Log(Time Step)')
        plt.ylabel('Log(Error)')
        plt.title('Convergence Order Calculation')
        plt.legend()
        plt.grid(True)
    
    plt.tight_layout()
    plt.savefig('convergence_order.png')
    print(f"Plots saved to convergence_order.png")
    
    # Print results
    print("\nConvergence Order Analysis Results:")
    print("-----------------------------------")
    print(f"Time Steps: {results}")
    print(f"Errors: {errors}")
    print(f"Individual Orders: {orders}")
    print(f"Average Convergence Order: {avg_order}")
    
    if avg_order is not None:
        print(f"\nExpected order for Explicit Euler: 1.0")
        print(f"Measured order: {avg_order:.2f}")
        
        if 0.8 <= avg_order <= 1.2:
            print("✓ Results are consistent with a first-order method (Explicit Euler)")
        else:
            print("⚠ Measured order differs from expected first-order behavior")
            if avg_order < 0.8:
                print("  - Order is lower than expected, possibly due to error accumulation or instabilities")
            else:
                print("  - Order is higher than expected, possibly due to special case or error measurement issues")
    
    return avg_order

def main():
    parser = argparse.ArgumentParser(description='Measure convergence order for N-body simulation')
    parser.add_argument('--file', type=str, default='step-0.sh', 
                        help='Path to shell script with particle data')
    parser.add_argument('--base-dt', type=float, default=None, 
                        help='Base time step (if not provided, will be calculated)')
    parser.add_argument('--safety', type=float, default=0.2, 
                        help='Safety factor for calculated time step')
    parser.add_argument('--steps', type=int, default=5, 
                        help='Number of time step refinements to test')
    parser.add_argument('--factor', type=float, default=2.0, 
                        help='Refinement factor between time steps')
    parser.add_argument('--timeout', type=int, default=300,
                        help='Simulation timeout in seconds')
    
    args = parser.parse_args()
    
    # Calculate base time step if not provided
    if args.base_dt is None:
        print("Calculating stable time step...")
        particles = read_particle_data(args.file)
        if particles:
            result = calculate_stable_timestep(particles, args.safety)
            base_dt = result['dt']
            print(f"Calculated stable time step: {base_dt}")
        else:
            print("Error: No particle data found. Please specify --base-dt manually.")
            return
    else:
        base_dt = args.base_dt
        print(f"Using provided base time step: {base_dt}")
    
    # Measure convergence order
    measure_convergence_order(base_dt, args.steps, args.factor, args.file, args.timeout)

if __name__ == "__main__":
    main()