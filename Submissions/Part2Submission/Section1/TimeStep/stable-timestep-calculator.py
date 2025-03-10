#!/usr/bin/env python3
import sys
import os
import math
import numpy as np
import argparse

def read_particle_data(filepath):
    """Read particle data from a shell script file."""
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Extract parameters from the shell script
    params = content.split()
    
    # Skip executable name and first 3 parameters (tPlotDelta, tFinal, timeStepSize)
    particle_params = params[4:]
    
    # Parse particle data - each particle has 7 parameters
    particles = []
    for i in range(0, len(particle_params), 7):
        if i+6 < len(particle_params):
            particle = {
                'position': [float(particle_params[i]), float(particle_params[i+1]), float(particle_params[i+2])],
                'velocity': [float(particle_params[i+3]), float(particle_params[i+4]), float(particle_params[i+5])],
                'mass': float(particle_params[i+6])
            }
            particles.append(particle)
    
    return particles

def distance(p1, p2):
    """Calculate Euclidean distance between two particles."""
    return math.sqrt(
        (p1['position'][0] - p2['position'][0])**2 +
        (p1['position'][1] - p2['position'][1])**2 +
        (p1['position'][2] - p2['position'][2])**2
    )

def velocity_magnitude(p):
    """Calculate velocity magnitude of a particle."""
    return math.sqrt(
        p['velocity'][0]**2 +
        p['velocity'][1]**2 +
        p['velocity'][2]**2
    )

def calculate_acceleration(p1, p2):
    """Calculate gravitational acceleration between two particles."""
    G = 1.0  # Gravitational constant
    dist = distance(p1, p2)
    
    if dist == 0:
        return 0
    
    # For n-body gravity, the acceleration of p1 due to p2 is: a = G * m2 / r^2
    return G * p2['mass'] / (dist * dist)

def calculate_stable_timestep(particles, safety_factor=0.2):
    """Calculate a stable time step for the given particles."""
    if not particles:
        return None
    
    min_distance = float('inf')
    max_velocity = 0
    max_acceleration = 0
    
    # Calculate minimum distance between any two particles
    for i in range(len(particles)):
        for j in range(i+1, len(particles)):
            dist = distance(particles[i], particles[j])
            min_distance = min(min_distance, dist)
    
    # Find maximum velocity
    for p in particles:
        max_velocity = max(max_velocity, velocity_magnitude(p))
    
    # Calculate maximum acceleration
    for i in range(len(particles)):
        for j in range(len(particles)):
            if i != j:
                acc = calculate_acceleration(particles[i], particles[j])
                max_acceleration = max(max_acceleration, acc)
    
    # If all particles are stationary, use a reasonable default
    if max_velocity == 0:
        dt_collision = 0.1
    else:
        dt_collision = min_distance / max_velocity
    
    # If no gravitational acceleration, use a reasonable default
    if max_acceleration == 0:
        dt_gravity = 0.1
    else:
        dt_gravity = math.sqrt(min_distance / max_acceleration)
    
    # Choose the smaller of the two constraints, with safety factor
    dt = safety_factor * min(dt_collision, dt_gravity)
    
    return {
        'dt': dt,
        'min_distance': min_distance,
        'max_velocity': max_velocity,
        'max_acceleration': max_acceleration,
        'dt_collision': dt_collision,
        'dt_gravity': dt_gravity
    }

def main():
    parser = argparse.ArgumentParser(description='Calculate a stable time step for N-body simulation')
    parser.add_argument('--file', type=str, default='step-0.sh', 
                        help='Path to shell script with particle data')
    parser.add_argument('--safety', type=float, default=0.2, 
                        help='Safety factor (eta) between 0.1 and 0.5')
    args = parser.parse_args()
    
    # Check if file exists
    if not os.path.exists(args.file):
        print(f"Error: File {args.file} not found")
        return
    
    # Read particle data
    print(f"Reading particle data from {args.file}...")
    particles = read_particle_data(args.file)
    print(f"Found {len(particles)} particles")
    
    # Calculate stable time step
    result = calculate_stable_timestep(particles, args.safety)
    
    print("\n=============================================")
    print("Stable Time Step Calculation:")
    print("=============================================")
    print(f"Minimum distance between particles: {result['min_distance']:.6e}")
    print(f"Maximum velocity of any particle: {result['max_velocity']:.6e}")
    print(f"Maximum acceleration due to gravity: {result['max_acceleration']:.6e}")
    print("---------------------------------------------")
    print(f"Time step based on collision constraint: {result['dt_collision']:.6e}")
    print(f"Time step based on gravitational constraint: {result['dt_gravity']:.6e}")
    print("---------------------------------------------")
    print(f"Final stable time step (with safety factor {args.safety}): {result['dt']:.6e}")
    print("=============================================")
    
    # Determine which constraint is limiting
    if result['dt_collision'] < result['dt_gravity']:
        print("Stability is limited by the collision constraint")
    else:
        print("Stability is limited by the gravitational constraint")
    
    return result['dt']

if __name__ == "__main__":
    main()