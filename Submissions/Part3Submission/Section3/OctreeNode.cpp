// OctreeNode.cpp
#include "OctreeNode.h"

// Define the static THETA constant
const double OctreeNode::THETA = 0.5;

OctreeNode::OctreeNode(double cx, double cy, double cz, double hw) :
    centerX(cx), centerY(cy), centerZ(cz), halfWidth(hw),
    totalMass(0.0), centerOfMassX(0.0), centerOfMassY(0.0), centerOfMassZ(0.0),
    isInternal(false) {
    
    // Initialize children pointers to null
    for (int i = 0; i < 8; i++) {
        children[i] = nullptr;
    }
}

OctreeNode::~OctreeNode() {
    // Clean up children
    for (int i = 0; i < 8; i++) {
        if (children[i] != nullptr) {
            delete children[i];
            children[i] = nullptr;
        }
    }
}

void OctreeNode::insert(int bodyIdx, double** x, double* mass) {
    // If this node doesn't contain any bodies yet
    if (totalMass == 0.0) {
        bodyIndices.push_back(bodyIdx);
        totalMass = mass[bodyIdx];
        centerOfMassX = x[bodyIdx][0];
        centerOfMassY = x[bodyIdx][1];
        centerOfMassZ = x[bodyIdx][2];
        return;
    }

    // If this is a leaf node with a body already, split it
    if (!isInternal && !bodyIndices.empty()) {
        int existingBodyIdx = bodyIndices[0];
        bodyIndices.clear();

        // Mark as internal node
        isInternal = true;

        // Determine which octant the existing body belongs to
        int octant = getOctant(x[existingBodyIdx][0], x[existingBodyIdx][1], x[existingBodyIdx][2]);
        
        // Create that octant and insert the existing body
        double newHalfWidth = halfWidth / 2.0;
        double newCenterX = centerX + (octant & 1 ? newHalfWidth : -newHalfWidth);
        double newCenterY = centerY + (octant & 2 ? newHalfWidth : -newHalfWidth);
        double newCenterZ = centerZ + (octant & 4 ? newHalfWidth : -newHalfWidth);
        
        children[octant] = new OctreeNode(newCenterX, newCenterY, newCenterZ, newHalfWidth);
        children[octant]->insert(existingBodyIdx, x, mass);
    }

    // Now insert the new body
    if (isInternal) {
        // Determine which octant the new body belongs to
        int octant = getOctant(x[bodyIdx][0], x[bodyIdx][1], x[bodyIdx][2]);
        
        // Create the octant if it doesn't exist
        if (children[octant] == nullptr) {
            double newHalfWidth = halfWidth / 2.0;
            double newCenterX = centerX + (octant & 1 ? newHalfWidth : -newHalfWidth);
            double newCenterY = centerY + (octant & 2 ? newHalfWidth : -newHalfWidth);
            double newCenterZ = centerZ + (octant & 4 ? newHalfWidth : -newHalfWidth);
            
            children[octant] = new OctreeNode(newCenterX, newCenterY, newCenterZ, newHalfWidth);
        }
        
        // Insert the body into the appropriate octant
        children[octant]->insert(bodyIdx, x, mass);
        
        // Update the center of mass of this node
        updateCenterOfMass(x, mass);
    }
}

int OctreeNode::getOctant(double posX, double posY, double posZ) {
    int octant = 0;
    if (posX >= centerX) octant |= 1;
    if (posY >= centerY) octant |= 2;
    if (posZ >= centerZ) octant |= 4;
    return octant;
}

void OctreeNode::updateCenterOfMass(double** x, double* mass) {
    totalMass = 0.0;
    centerOfMassX = centerOfMassY = centerOfMassZ = 0.0;
    
    // For leaf nodes
    if (!isInternal) {
        for (int bodyIdx : bodyIndices) {
            totalMass += mass[bodyIdx];
            centerOfMassX += mass[bodyIdx] * x[bodyIdx][0];
            centerOfMassY += mass[bodyIdx] * x[bodyIdx][1];
            centerOfMassZ += mass[bodyIdx] * x[bodyIdx][2];
        }
    } 
    // For internal nodes
    else {
        for (int i = 0; i < 8; i++) {
            if (children[i] != nullptr) {
                totalMass += children[i]->totalMass;
                centerOfMassX += children[i]->totalMass * children[i]->centerOfMassX;
                centerOfMassY += children[i]->totalMass * children[i]->centerOfMassY;
                centerOfMassZ += children[i]->totalMass * children[i]->centerOfMassZ;
            }
        }
    }
    
    // Compute actual center of mass
    if (totalMass > 0) {
        centerOfMassX /= totalMass;
        centerOfMassY /= totalMass;
        centerOfMassZ /= totalMass;
    }
}

bool OctreeNode::isFarEnough(double posX, double posY, double posZ) {
    // Calculate distance to center of mass
    double dx = centerOfMassX - posX;
    double dy = centerOfMassY - posY;
    double dz = centerOfMassZ - posZ;
    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    // If distance is zero, we need direct calculation
    if (distance == 0.0) return false;
    
    // Check if ratio of size to distance is less than threshold
    return (halfWidth / distance < THETA);
}

void OctreeNode::computeForce(int bodyIdx, double** x, double* force, double* mass) {
    // If this is a leaf node with the same body, skip
    if (!isInternal && bodyIndices.size() == 1 && bodyIndices[0] == bodyIdx) {
        return;
    }
    
    // Check if we can use approximation
    if (!isInternal || isFarEnough(x[bodyIdx][0], x[bodyIdx][1], x[bodyIdx][2])) {
        // Use approximation - compute force between body and center of mass
        calculateForce(bodyIdx, x, force, mass);
    } 
    // Otherwise, recursively compute forces from children
    else {
        for (int i = 0; i < 8; i++) {
            if (children[i] != nullptr) {
                children[i]->computeForce(bodyIdx, x, force, mass);
            }
        }
    }
}

void OctreeNode::calculateForce(int bodyIdx, double** x, double* force, double* mass) {
    // Skip if node is empty
    if (totalMass == 0.0) return;
    
    // Calculate distance
    double dx = centerOfMassX - x[bodyIdx][0];
    double dy = centerOfMassY - x[bodyIdx][1];
    double dz = centerOfMassZ - x[bodyIdx][2];
    
    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    // Avoid division by zero
    if (distance < 1e-10) return;
    
    // F = G * m1 * m2 / r^2 * (r_vec / r)
    // G is absorbed in the mass values in this simulation
    double force_magnitude = (mass[bodyIdx] * totalMass) / (distance * distance * distance);
    
    force[0] += force_magnitude * dx;
    force[1] += force_magnitude * dy;
    force[2] += force_magnitude * dz;
}