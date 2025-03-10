# Step 3.1

## Report

### Section 1

- The main calculation is done in the force calculation loop which has a triangular iteration pattern as the outer loop iterates N times, and the inner loop iterates N-1, N-2, ... and so on as i increases in the outer loop
    - This leads to an O(N^2) time complexity as a worst case
    - The force calculation function itself has O(1) as it is a fixed speed
        - This is called 6 times for each iteration -> still linear
    - NEED TO RUN CODE WITH VARIOUS INPUT SIZES AND PLOT GRAPH TO SHOW THIS COMPLEXITY

- Improvement -> Barnes-Hut Algorithm
    - Read about it in lecture slides
    - Basic description I guess will suffice?

- Time step for step-0
    - Explicit euler method used
    - Code to implement the optimal time step
        - Order of time stepping scheme? 

- Convergence order
    - Code in place to hard measure it
        - On tower = 1.08928
            - For p around 1, that is good
    - Estimated order -> need to estimate this based off time stepping -> not sure of method
    - Estimations based off theory
        - A theoretical order of convergence of 1.0, similar to Symplectic Euler
        - Less favorable error behavior over long time periods compared to Symplectic Euler
        - Energy errors that grow systematically rather than oscillating around a constant value
        - Artifacts in orbital motion that may affect the convergence analysis, especially for longer test times

- Tolerance for error
    - Tolerance is proportional h^p where h is time step used and p is convergence order (1.08 roughly)

- Optimisations for code
    - Not done yet
    - Talk about data structures

### Section 2 - Step-3
- Implemented basic approach from claude in step-3
- Need to run a full test to see if the collisions are working

- Convergence changes
    - 

### Section 3 - Step-4
- Steal adpative time step method from step-0 currently there and use this
- Need to implement it and test it to see if it works
- Then do write up


