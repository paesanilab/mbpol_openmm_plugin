#ifndef OPENMM_REFERENCE_THREE_NEIGHBORLIST_H_
#define OPENMM_REFERENCE_THREE_NEIGHBORLIST_H_

#include "RealVec.h"
#include "openmm/internal/windowsExport.h"
#include <set>
#include <vector>

namespace OpenMM {

typedef std::vector<RealVec> AtomLocationList;
typedef unsigned int AtomIndex;

struct AtomTriplet
{
    AtomIndex  first, second, third;
};

typedef std::vector<AtomTriplet>  ThreeNeighborList;

// O(n) neighbor list method using voxel hash data structure
// parameter neighborList is automatically clear()ed before 
// neighbors are added
void OPENMM_EXPORT computeNeighborListVoxelHash(
                              ThreeNeighborList& neighborList,
                              int nAtoms,
                              const AtomLocationList& atomLocations, 
                              const RealVec& periodicBoxSize,
                              bool usePeriodic,
                              double maxDistance,
                              double minDistance = 0.0
                             );

} // namespace OpenMM

#endif // OPENMM_REFERENCE_NEIGHBORLIST_H_
