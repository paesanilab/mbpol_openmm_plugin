/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/internal/AssertionUtilities.h"
#include "ReferenceThreeNeighborList.h"
#include "sfmt/SFMT.h"
#include <cassert>
#include <iostream>
#include <vector>

using namespace std;
using namespace OpenMM;

void testNeighborList()
{
    vector<RealVec> particleList(3);
    particleList[0] = RealVec(13.6, 0, 0);
    particleList[1] = RealVec(0, 0, 0);
    particleList[2] = RealVec(2, 0, 0);
    
    ThreeNeighborList neighborList;

    RealVec boxSize;
    
    computeNeighborListVoxelHash(neighborList, particleList.size(), particleList, boxSize, false, 13.7, 0.01);
    assert(neighborList.size() == 1);
    
    computeNeighborListVoxelHash(neighborList, particleList.size(), particleList, boxSize, false, 13.5, 0.01);
    assert(neighborList.size() == 0);
}

void testNeighborListFourAtoms()
{
    vector<RealVec> particleList(4);
    particleList[0] = RealVec(14, 0, 0);
    particleList[1] = RealVec(10, 0, 0);
    particleList[2] = RealVec(10, 0, 0);
    particleList[3] = RealVec(10, 0, 0);

    ThreeNeighborList neighborList;

    RealVec boxSize;

    computeNeighborListVoxelHash(neighborList, particleList.size(), particleList, boxSize, false, 12, 0);
    assert(neighborList.size() == 4);

    computeNeighborListVoxelHash(neighborList, particleList.size(), particleList, boxSize, false, 1, 0);
    assert(neighborList.size() == 1);

}


double periodicDifference(double val1, double val2, double period) {
    double diff = val1-val2;
    double base = floor(diff/period+0.5)*period;
    return diff-base;
}

double distance2(RealVec& pos1, RealVec& pos2, const RealVec& periodicBoxSize) {
    double dx = periodicDifference(pos1[0], pos2[0], periodicBoxSize[0]);
    double dy = periodicDifference(pos1[1], pos2[1], periodicBoxSize[1]);
    double dz = periodicDifference(pos1[2], pos2[2], periodicBoxSize[2]);
    return dx*dx+dy*dy+dz*dz;
}

void verifyNeighborList(ThreeNeighborList& list, int numParticles, vector<RealVec>& positions, const RealVec& periodicBoxSize, double cutoff) {

    int count = 0;
    for (int i = 0; i < numParticles-2; i++)
        for (int j = i+1; j < numParticles-1; j++)
            for (int q = j+1; q < numParticles; q++)
            if ((distance2(positions[i], positions[j], periodicBoxSize) <= cutoff*cutoff) and
                (distance2(positions[i], positions[q], periodicBoxSize) <= cutoff*cutoff) and
                (distance2(positions[j], positions[q], periodicBoxSize) <= cutoff*cutoff))
                count++;

    ASSERT(count == list.size());

    for (int i = 0; i < (int) list.size(); i++) {
        int particle1 = list[i].first;
        int particle2 = list[i].second;
        int particle3 = list[i].third;
        ASSERT(distance2(positions[particle1], positions[particle2], periodicBoxSize) <= cutoff*cutoff);
        ASSERT(distance2(positions[particle1], positions[particle3], periodicBoxSize) <= cutoff*cutoff);
        ASSERT(distance2(positions[particle2], positions[particle3], periodicBoxSize) <= cutoff*cutoff);
    }

}

void testPeriodic(double cutoff) {
    const int numParticles = 40;

    const RealVec periodicBoxSize(20.0, 15.0, 22.0);
    vector<RealVec> particleList(numParticles);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    for (int i = 0; i <numParticles; i++) {
        particleList[i][0] = (RealOpenMM) (genrand_real2(sfmt)*periodicBoxSize[0]*3);
        particleList[i][1] = (RealOpenMM) (genrand_real2(sfmt)*periodicBoxSize[1]*3);
        particleList[i][2] = (RealOpenMM) (genrand_real2(sfmt)*periodicBoxSize[2]*3);
    }
    ThreeNeighborList neighborList;
    computeNeighborListVoxelHash(neighborList, numParticles, particleList, periodicBoxSize, true, cutoff);

    verifyNeighborList(neighborList, numParticles, particleList, periodicBoxSize, cutoff);
}

int main() 
{
try {
    testNeighborList();

    testNeighborListFourAtoms();

    testPeriodic(10.);

    cout << "Test Passed" << endl;
    return 0;
}
catch (...) {
    cerr << "*** ERROR: Test Failed ***" << endl;
    return 1;
}
}

