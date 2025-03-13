#!/usr/bin/env python3
import math

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