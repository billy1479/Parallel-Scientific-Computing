#include "OctreeNode.h"

// Define the static constants
const double OctreeNode::THETA = 0.2;
const double OctreeNode::MIN_SIZE = 1e-6;  // Prevents infinite subdivision

OctreeNode::OctreeNode(double cx, double cy, double cz, double hw, int d) :
    centerX(cx), centerY(cy), centerZ(cz), 
    halfWidth(hw > MIN_SIZE ? hw : MIN_SIZE),  // Ensure minimum size
    totalMass(0.0), centerOfMassX(0.0), centerOfMassY(0.0), centerOfMassZ(0.0),
    isInternal(false), depth(d) {
    
    // Initialize children pointers to null
    for (int i = 0; i < 8; i++) {
        children[i] = nullptr;
    }
}

OctreeNode::~OctreeNode() {
    // Safe deletion of children
    for (int i = 0; i < 8; i++) {
        if (children[i] != nullptr) {
            try {
                OctreeNode* temp = children[i];
                children[i] = nullptr;  // Clear pointer before deletion to avoid cyclic issues
                delete temp;
            } catch (std::exception& e) {
                std::cerr << "Exception while deleting child " << i << ": " << e.what() << std::endl;
            } catch (...) {
                std::cerr << "Unknown exception while deleting child " << i << std::endl;
            }
        }
    }
}
int OctreeNode::getOctant(double posX, double posY, double posZ) {
    // Handle numerical precision issues near boundaries
    double epsilon = 1e-10;
    
    int octant = 0;
    if (posX >= centerX - epsilon) octant |= 1;
    if (posY >= centerY - epsilon) octant |= 2;
    if (posZ >= centerZ - epsilon) octant |= 4;
    return octant;
}

void OctreeNode::insert(int bodyIdx, double** x, double* mass) {
    try {
        // Safety check for null pointers
        if (x == nullptr || mass == nullptr) {
            std::cerr << "Error: Null pointer in insert method" << std::endl;
            return;
        }
        
        // If this node doesn't contain any bodies yet
        if (totalMass <= 0.0) {
            bodyIndices.push_back(bodyIdx);
            totalMass = mass[bodyIdx];
            centerOfMassX = x[bodyIdx][0];
            centerOfMassY = x[bodyIdx][1];
            centerOfMassZ = x[bodyIdx][2];
            return;
        }

        // If this is a leaf node with bodies, and we haven't reached max depth
        if (!isInternal && !bodyIndices.empty() && depth < MAX_DEPTH) {
            // Save existing bodies
            std::vector<int> oldBodies = bodyIndices;
            bodyIndices.clear();

            // Mark as internal node
            isInternal = true;

            // Re-insert all old bodies
            for (size_t i = 0; i < oldBodies.size(); i++) {
                int oldBodyIdx = oldBodies[i];
                
                // Handle each old body
                double posX = x[oldBodyIdx][0];
                double posY = x[oldBodyIdx][1];
                double posZ = x[oldBodyIdx][2];
                
                int octant = getOctant(posX, posY, posZ);
                
                // Create new octant if needed
                if (children[octant] == nullptr) {
                    // Calculate new center coordinates for the child node
                    double newHalfWidth = halfWidth / 2.0;
                    double newCenterX = centerX + ((octant & 1) ? newHalfWidth : -newHalfWidth);
                    double newCenterY = centerY + ((octant & 2) ? newHalfWidth : -newHalfWidth);
                    double newCenterZ = centerZ + ((octant & 4) ? newHalfWidth : -newHalfWidth);
                    
                    // Create child node with regular new operator
                    children[octant] = new OctreeNode(newCenterX, newCenterY, newCenterZ, 
                                                      newHalfWidth, depth + 1);
                }
                
                // Insert old body into the child
                if (children[octant] != nullptr) {
                    children[octant]->insert(oldBodyIdx, x, mass);
                }
            }
        }

        // Now handle the new body
        if (isInternal) {
            // Get position of new body
            double posX = x[bodyIdx][0];
            double posY = x[bodyIdx][1];
            double posZ = x[bodyIdx][2];
            
            // Determine which octant the new body belongs to
            int octant = getOctant(posX, posY, posZ);
            
            // Create new octant if needed
            if (children[octant] == nullptr) {
                double newHalfWidth = halfWidth / 2.0;
                double newCenterX = centerX + ((octant & 1) ? newHalfWidth : -newHalfWidth);
                double newCenterY = centerY + ((octant & 2) ? newHalfWidth : -newHalfWidth);
                double newCenterZ = centerZ + ((octant & 4) ? newHalfWidth : -newHalfWidth);
                
                // Create child node with regular new operator
                children[octant] = new OctreeNode(newCenterX, newCenterY, newCenterZ, 
                                                  newHalfWidth, depth + 1);
            }
            
            // Insert the body into the appropriate octant
            if (children[octant] != nullptr) {
                children[octant]->insert(bodyIdx, x, mass);
            }
            
            // Update the center of mass for this node
            updateCenterOfMass(x, mass);
        } else {
            // This is still a leaf node, just add the body
            bodyIndices.push_back(bodyIdx);
            updateCenterOfMass(x, mass);
        }
    } catch (std::exception& e) {
        std::cerr << "Exception in insert: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception in insert" << std::endl;
    }
}

void OctreeNode::updateCenterOfMass(double** x, double* mass) {
    try {
        // Reset
        totalMass = 0.0;
        centerOfMassX = centerOfMassY = centerOfMassZ = 0.0;
        
        // For leaf nodes with direct bodies
        if (!isInternal) {
            for (size_t i = 0; i < bodyIndices.size(); i++) {
                int bodyIdx = bodyIndices[i];
                double bodyMass = mass[bodyIdx];
                
                totalMass += bodyMass;
                centerOfMassX += bodyMass * x[bodyIdx][0];
                centerOfMassY += bodyMass * x[bodyIdx][1];
                centerOfMassZ += bodyMass * x[bodyIdx][2];
            }
        } 
        // For internal nodes, aggregate from children
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
    } catch (std::exception& e) {
        std::cerr << "Exception in updateCenterOfMass: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception in updateCenterOfMass" << std::endl;
    }
}

bool OctreeNode::isFarEnough(double posX, double posY, double posZ) {
    try {
        // Calculate distance to center of mass
        double dx = centerOfMassX - posX;
        double dy = centerOfMassY - posY;
        double dz = centerOfMassZ - posZ;
        double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
        
        // If distance is essentially zero, we need direct calculation
        if (distance < 1e-10) return false;
        
        // Check if ratio of size to distance is less than threshold
        return (halfWidth / distance < THETA);
    } catch (std::exception& e) {
        std::cerr << "Exception in isFarEnough: " << e.what() << std::endl;
        return false;
    } catch (...) {
        std::cerr << "Unknown exception in isFarEnough" << std::endl;
        return false;
    }
}

void OctreeNode::computeForce(int bodyIdx, double** x, double* force, double* mass) {
    try {
        // Skip empty nodes
        if (totalMass <= 0.0) return;
        
        // If this is a leaf node with the same body, skip
        if (!isInternal && bodyIndices.size() == 1 && bodyIndices[0] == bodyIdx) {
            return;
        }
        
        // Get position of the body
        double posX = x[bodyIdx][0];
        double posY = x[bodyIdx][1];
        double posZ = x[bodyIdx][2];
        
        // Check if we can use approximation
        if (!isInternal || isFarEnough(posX, posY, posZ)) {
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
    } catch (std::exception& e) {
        std::cerr << "Exception in computeForce: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception in computeForce" << std::endl;
    }
}

void OctreeNode::calculateForce(int bodyIdx, double** x, double* force, double* mass) {
    try {
        // Skip if node is empty
        if (totalMass <= 0.0) return;
        
        // Calculate distance
        double dx = centerOfMassX - x[bodyIdx][0];
        double dy = centerOfMassY - x[bodyIdx][1];
        double dz = centerOfMassZ - x[bodyIdx][2];
        
        double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
        
        // Avoid division by zero - add a small softening factor
        if (distance < 1e-10) {
            // Bodies are too close, apply a small random perturbation
            // This helps avoid infinite forces
            return;
        }
        
        // F = G * m1 * m2 / r^2 * (r_vec / r)
        // G is absorbed in the mass values in this simulation
        double force_magnitude = (mass[bodyIdx] * totalMass) / (distance * distance * distance);
        
        force[0] += force_magnitude * dx;
        force[1] += force_magnitude * dy;
        force[2] += force_magnitude * dz;
    } catch (std::exception& e) {
        std::cerr << "Exception in calculateForce: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception in calculateForce" << std::endl;
    }
}