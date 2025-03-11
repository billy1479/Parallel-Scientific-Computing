#pragma once
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

class OctreeNode {
private:
    // Spatial boundaries
    double centerX, centerY, centerZ;
    double halfWidth;  // Half the width of this cube

    // Physical properties
    double totalMass;  // Total mass of bodies in this node
    double centerOfMassX, centerOfMassY, centerOfMassZ;  // Center of mass

    // Children nodes (octants)
    OctreeNode* children[8];

    // Bodies contained directly in this node (only for leaf nodes)
    std::vector<int> bodyIndices;  // Indices into the simulation's body arrays

    // Flag indicating if this is an internal node (has children)
    bool isInternal;
    
    // Depth tracking to prevent infinite recursion
    int depth;
    static const int MAX_DEPTH = 20;

    // Constants
    static const double THETA;  // Threshold for approximation

    // Minimum size for an octree node - prevents infinite subdivision
    static const double MIN_SIZE;

public:
    // Constructor
    OctreeNode(double cx, double cy, double cz, double hw, int d = 0);
    
    // Destructor
    ~OctreeNode();

    // Insert a body into the octree
    void insert(int bodyIdx, double** x, double* mass);

    // Compute force on a body
    void computeForce(int bodyIdx, double** x, double* force, double* mass);

    // Calculate which octant a body belongs to (0-7)
    int getOctant(double posX, double posY, double posZ);

    // Check if a node is far enough for approximation
    bool isFarEnough(double posX, double posY, double posZ);

    // Calculate force between a body and a node
    void calculateForce(int bodyIdx, double** x, double* force, double* mass);

    // Update center of mass for this node
    void updateCenterOfMass(double** x, double* mass);
};