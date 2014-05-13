#include "ReferenceThreeNeighborList.h"
#include <set>
#include <map>
#include <cmath>
#include <iostream>
#include <cassert>

using namespace std;

namespace OpenMM {

typedef std::vector<AtomIndex> AtomList;

static double periodicDifference(double val1, double val2, double period) {
    double diff = val1-val2;
    double base = floor(diff/period+0.5)*period;
    return diff-base;
}

// squared distance between two points
static double compPairDistanceSquared(const RealVec& pos1, const RealVec& pos2, const RealVec& periodicBoxSize, bool usePeriodic) {
    double dx, dy, dz;
    if (!usePeriodic) {
        dx = pos2[0] - pos1[0];
        dy = pos2[1] - pos1[1];
        dz = pos2[2] - pos1[2];
    }
    else {
        dx = periodicDifference(pos2[0], pos1[0], periodicBoxSize[0]);
        dy = periodicDifference(pos2[1], pos1[1], periodicBoxSize[1]);
        dz = periodicDifference(pos2[2], pos1[2], periodicBoxSize[2]);
    }
    return dx*dx + dy*dy + dz*dz;
}

class VoxelIndex 
{
public:
    VoxelIndex(int xx, int yy, int zz) : x(xx), y(yy), z(zz) {}

    // operator<() needed for map
    bool operator<(const VoxelIndex& other) const {
        if      (x < other.x) return true;
        else if (x > other.x) return false;
        else if (y < other.y) return true;
        else if (y > other.y) return false;
        else if (z < other.z) return true;
        else return false;
    }

    int x;
    int y;
    int z;
};


typedef std::pair<const RealVec*, AtomIndex> VoxelItem;
typedef std::vector< VoxelItem > Voxel;

class VoxelHash
{
public:
    VoxelHash(double vsx, double vsy, double vsz, const RealVec& periodicBoxSize, bool usePeriodic) :
            voxelSizeX(vsx), voxelSizeY(vsy), voxelSizeZ(vsz), periodicBoxSize(periodicBoxSize), usePeriodic(usePeriodic) {
        if (usePeriodic) {
            nx = (int) floor(periodicBoxSize[0]/voxelSizeX+0.5);
            ny = (int) floor(periodicBoxSize[1]/voxelSizeY+0.5);
            nz = (int) floor(periodicBoxSize[2]/voxelSizeZ+0.5);
        }
    }

    void insert(const AtomIndex& item, const RealVec& location)
    {
        VoxelIndex voxelIndex = getVoxelIndex(location);
        if ( voxelMap.find(voxelIndex) == voxelMap.end() ) voxelMap[voxelIndex] = Voxel(); 
        Voxel& voxel = voxelMap.find(voxelIndex)->second;
        voxel.push_back( VoxelItem(&location, item) );
    }


    VoxelIndex getVoxelIndex(const RealVec& location) const {
        double xperiodic, yperiodic, zperiodic;
        if (!usePeriodic) {
            xperiodic = location[0];
            yperiodic = location[1];
            zperiodic = location[2];
        }
        else {
            xperiodic = location[0]-periodicBoxSize[0]*floor(location[0]/periodicBoxSize[0]);
            yperiodic = location[1]-periodicBoxSize[1]*floor(location[1]/periodicBoxSize[1]);
            zperiodic = location[2]-periodicBoxSize[2]*floor(location[2]/periodicBoxSize[2]);
        }
        int x = int(floor(xperiodic / voxelSizeX));
        int y = int(floor(yperiodic / voxelSizeY));
        int z = int(floor(zperiodic / voxelSizeZ));
        
        return VoxelIndex(x, y, z);
    }

    void getNeighbors(
            ThreeNeighborList& neighbors,
            const VoxelItem& referencePoint, 
            double maxDistance, 
            double minDistance) const 
    {

        // Loop over neighboring voxels
        // TODO use more clever selection of neighboring voxels
        assert(maxDistance > 0);
        assert(minDistance >= 0);
        assert(voxelSizeX > 0);
        assert(voxelSizeY > 0);
        assert(voxelSizeZ > 0);

        const AtomIndex atomI = referencePoint.second;
        const RealVec& locationI = *referencePoint.first;
        
        double maxDistanceSquared = maxDistance * maxDistance;
        double minDistanceSquared = minDistance * minDistance;

        int dIndexX = int(maxDistance / voxelSizeX) + 1; // How may voxels away do we have to look?
        int dIndexY = int(maxDistance / voxelSizeY) + 1;
        int dIndexZ = int(maxDistance / voxelSizeZ) + 1;
        VoxelIndex centerVoxelIndex = getVoxelIndex(locationI);
        int lastx = centerVoxelIndex.x+dIndexX;
        int lasty = centerVoxelIndex.y+dIndexY;
        int lastz = centerVoxelIndex.z+dIndexZ;
        if (usePeriodic) {
            lastx = min(lastx, centerVoxelIndex.x-dIndexX+nx-1);
            lasty = min(lasty, centerVoxelIndex.y-dIndexY+ny-1);
            lastz = min(lastz, centerVoxelIndex.z-dIndexZ+nz-1);
        }

        std::vector<AtomIndex> nearbyAtoms;
        std::vector<const RealVec*> nearbyAtomsLocation;

        // Loop through nearby Voxels and create vectors of nearby Atoms and their locations
        for (int x = centerVoxelIndex.x - dIndexX; x <= lastx; ++x)
        {
            for (int y = centerVoxelIndex.y - dIndexY; y <= lasty; ++y)
            {
                for (int z = centerVoxelIndex.z - dIndexZ; z <= lastz; ++z)
                {
                    VoxelIndex voxelIndex(x, y, z);
                    if (usePeriodic) {
                        voxelIndex.x = (x+nx)%nx;
                        voxelIndex.y = (y+ny)%ny;
                        voxelIndex.z = (z+nz)%nz;
                    }
                    if (voxelMap.find(voxelIndex) == voxelMap.end()) continue; // no such voxel; skip
                    const Voxel& voxel = voxelMap.find(voxelIndex)->second;
                    for (Voxel::const_iterator itemIter = voxel.begin(); itemIter != voxel.end(); ++itemIter)
                    {
                        const AtomIndex atomJ = itemIter->second;
                        const RealVec& locationJ = *itemIter->first;
                        
                        // Ignore self hits
                        if (atomI == atomJ) continue;
                        
                        double dSquared = compPairDistanceSquared(locationI, locationJ, periodicBoxSize, usePeriodic);
                        if (dSquared > maxDistanceSquared) continue;
                        if (dSquared < minDistanceSquared) continue;

                        nearbyAtoms.push_back( atomJ );
                        nearbyAtomsLocation.push_back( itemIter->first );
                    } // for atom indices
                }  // for z
            }  // for y
        } // for x


        if (nearbyAtoms.size() >= 2)
        {
            // Loop through each combination of the reference atom (atomI) with each pair of nearby atoms atomJ and atomQ
            for (unsigned int atomJindex=0; atomJindex<(nearbyAtoms.size()-1); atomJindex++)
            {
                for (unsigned int atomQindex=atomJindex+1; atomQindex<nearbyAtoms.size(); atomQindex++)
                {
                    // check only distance Q-J, we checked I-J and I-Q in the previous loop
                    double dSquared = compPairDistanceSquared(*nearbyAtomsLocation[atomJindex], *nearbyAtomsLocation[atomQindex], periodicBoxSize, usePeriodic);
                    if (dSquared > maxDistanceSquared) continue;
                    if (dSquared < minDistanceSquared) continue;

                    AtomTriplet triplet = {atomI, nearbyAtoms[atomJindex], nearbyAtoms[atomQindex]};
                    neighbors.push_back(triplet);

                }
            }
        }
    }

private:
    double voxelSizeX, voxelSizeY, voxelSizeZ;
    int nx, ny, nz;
    const RealVec& periodicBoxSize;
    const bool usePeriodic;
    std::map<VoxelIndex, Voxel> voxelMap;
};


// O(n) neighbor list method using voxel hash data structure
void OPENMM_EXPORT computeNeighborListVoxelHash(
                              ThreeNeighborList& neighborList,
                              int nAtoms,
                              const AtomLocationList& atomLocations, 
                              const RealVec& periodicBoxSize,
                              bool usePeriodic,
                              double maxDistance,
                              double minDistance                             )
{
    neighborList.clear();

    double edgeSizeX, edgeSizeY, edgeSizeZ;
    if (!usePeriodic)
        edgeSizeX = edgeSizeY = edgeSizeZ = maxDistance; // TODO - adjust this as needed
    else {
        edgeSizeX = periodicBoxSize[0]/floor(periodicBoxSize[0]/maxDistance);
        edgeSizeY = periodicBoxSize[1]/floor(periodicBoxSize[1]/maxDistance);
        edgeSizeZ = periodicBoxSize[2]/floor(periodicBoxSize[2]/maxDistance);
    }
    VoxelHash voxelHash(edgeSizeX, edgeSizeY, edgeSizeZ, periodicBoxSize, usePeriodic);
    for (AtomIndex atomJ = 0; atomJ < (AtomIndex) nAtoms; ++atomJ) // use "j", because j > i for pairs
    {
        // 1) Find other atoms that are close to this one
        const RealVec& location = atomLocations[atomJ];
        voxelHash.getNeighbors( 
            neighborList, 
            VoxelItem(&location, atomJ),
            maxDistance, 
            minDistance);
            
        // 2) Add this atom to the voxelHash
        voxelHash.insert(atomJ, location);
    }
}

} // namespace OpenMM
