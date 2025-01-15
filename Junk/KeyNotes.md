# Step-1 - Vectorisation
- Force update - reduced function calls by storing intermediate value in variables
- Force update - openmp simd for vectorisation with collapse(2) to break down the nested loop
- Force calculation - openmp declare simd vectorisation
- SIMD enabled version of square root, min and max function from std

# Step-2 - Parallelisation
- Brought SIMD enabled maths fucntions from step-1
- enabled parallel for schedule(dynamic) on force update loop
- enabled parallelisation for the subloop from above
- Ensuring variable visiblty via shared() and private() for loop definitions -> see lecture slides for more information in report -> variable i made private variable to do this
- All-to-one synch for the for loop updating forces for force0, force1, force2 