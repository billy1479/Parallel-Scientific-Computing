# Step 3.1

## Report

- The main calculation is done in the force calculation loop which has a triangular iteration pattern as the outer loop iterates N times, and the inner loop iterates N-1, N-2, ... and so on as i increases in the outer loop
    - This leads to an O(N^2) time complexity as a worst case
    - The force caluation function itself has O(1) as it is a fixed speed
        - This is called 6 times for each iteration -> still linear
    - NEED TO RUN CODE WITH VARIOUS INPUT SIZES AND PLOT GRAPH TO SHOW THIS COMPLEXITY

- Improvement -> Barnes-Hut Algorithm
    - Read about it in lecture slides
    - Basic description I guess will suffice?

- Time step for step-0
    - Seems simple to implement via code
    - Will need to mathematically justify?
    - Fetched from Claude .... 
    Let me explain how the lecture content helps us understand and justify the stability analysis approach for this N-body simulation.

From the lecture slides, we learn that for explicit Euler methods, stability analysis is crucial and can be understood through a few key principles:

1. Physical vs. Numerical Stability
The slides show that a system can be physically stable (like the example λ ≤ 0 in the Dahlquist test equation) while the numerical method might be unstable depending on the timestep. This is exactly what we see in N-body simulations - while the physical system conserves energy, the numerical solution can become unstable.

2. Local Stability Analysis
From the slides, we see that stability can be analyzed by looking at how errors propagate step by step:
```
|1 + λh| ≤ 1
```
This is the key stability criterion shown in the lecture.

For the N-body code, we can apply this same principle by considering the force calculations. The force between two bodies scales as:
```cpp
const double distance3 = distance * distance * distance;
return (x[i][direction]-x[j][direction]) * mass[i]*mass[j] / distance3;
```

When bodies get very close, this term becomes very large, effectively giving us a large negative λ in the stability criterion. Following the lecture's stability analysis:

1. The timestep h must satisfy |1 + λh| ≤ 1 for stability
2. For large negative λ (close bodies), this means we need a very small h
3. This explains why we need adaptive timesteps based on the minimum distance between bodies:
```cpp
minDx = std::min(minDx, distance);
```

The theoretical largest stable timestep would be:
```cpp
h ≤ 2/|λ| ∝ minDx^3/(G * max_mass)
```

This matches what the lecture calls a "stiff" problem - where the timestep is constrained by stability rather than accuracy requirements. As shown in the lecture graphs, when the timestep is too large, the numerical solution "blows up" despite the physical system being stable.

The lecture also confirms that this is a fundamental limitation of explicit Euler - there are other timestepping schemes without this restriction, but they typically require implicit steps which would be computationally expensive for N-body problems.

This theoretical analysis explains and justifies why the code:
1. Tracks minimum distances between bodies
2. Uses these to constrain timestep size
3. Requires smaller timesteps when bodies get close together

Would you like me to explain any part of this analysis in more detail? I can particularly elaborate on how the stability criterion from the lecture maps to specific terms in the N-body force calculations.

- 