/* -------------------------------------------------------------------------- *
 *                               OpenMMMBPol                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
 * Authors:                                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "MBPolCpuKernels.h"
#include "openmm/cpu/CpuPlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/MBPolElectrostaticsForce.h"
#include "openmm/internal/MBPolElectrostaticsForceImpl.h"
#include "openmm/internal/MBPolTwoBodyForceImpl.h"
#include "openmm/internal/MBPolThreeBodyForceImpl.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include <iostream>

#include <cmath>
#ifdef _MSC_VER
#include <windows.h>
#endif

using namespace  OpenMM;
using namespace MBPolPlugin;
using namespace std;

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->velocities);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

static RealVec& extractBoxSize(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *(RealVec*) data->periodicBoxSize;
}
#if !(OPENMM_MAJOR_VERSION == 6 && OPENMM_MINOR_VERSION <= 2)
static RealVec* extractBoxVectors(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return (RealVec*) data->periodicBoxVectors;
}
#endif

/* -------------------------------------------------------------------------- *
 *                             MBPolElectrostatics                                *
 * -------------------------------------------------------------------------- */

CpuCalcMBPolElectrostaticsForceKernel::CpuCalcMBPolElectrostaticsForceKernel(std::string name, const Platform& platform, CpuPlatform::PlatformData& data) :
         CalcMBPolElectrostaticsForceKernel(name, platform), data(data), numElectrostatics(0), mutualInducedMaxIterations(200), mutualInducedTargetEpsilon(1.0e-03),
                                                         usePme(false),alphaEwald(0.0), cutoffDistance(1.0) {  

}

CpuCalcMBPolElectrostaticsForceKernel::~CpuCalcMBPolElectrostaticsForceKernel() {
}

void CpuCalcMBPolElectrostaticsForceKernel::initialize(const OpenMM::System& system, const MBPolElectrostaticsForce& force) {

    numElectrostatics   = force.getNumElectrostatics();

    charges.resize(numElectrostatics);
    moleculeIndices.resize(numElectrostatics);
	atomTypes.resize(numElectrostatics);
    dampingFactors.resize(numElectrostatics);
    polarity.resize(numElectrostatics);

    int tholeIndex       = 0;
    for( int ii = 0; ii < numElectrostatics; ii++ ){

        // multipoles

        double charge, dampingFactorD, polarityD;
        std::vector<double> dipolesD;
		int moleculeIndex, atomType;
        force.getElectrostaticsParameters(ii, charge,
                                     moleculeIndex, atomType, dampingFactorD, polarityD );

        charges[ii]                        = static_cast<RealOpenMM>(charge);

        dampingFactors[ii]                 = static_cast<RealOpenMM>(dampingFactorD);
        polarity[ii]                       = static_cast<RealOpenMM>(polarityD);

		moleculeIndices[ii] = moleculeIndex;
		atomTypes[ii] = atomType;

    }

    mutualInducedMaxIterations = force.getMutualInducedMaxIterations();
    mutualInducedTargetEpsilon = force.getMutualInducedTargetEpsilon();

    includeChargeRedistribution = force.getIncludeChargeRedistribution();
    tholeParameters = force.getTholeParameters();

    // PME

    nonbondedMethod  = force.getNonbondedMethod();
    if( nonbondedMethod == MBPolElectrostaticsForce::PME ){
        usePme     = true;
        alphaEwald = force.getAEwald();
        cutoffDistance = force.getCutoffDistance();
        force.getPmeGridDimensions(pmeGridDimension);
        if (pmeGridDimension[0] == 0 || alphaEwald == 0.0) {
            NonbondedForce nb;
            nb.setEwaldErrorTolerance(force.getEwaldErrorTolerance());
            nb.setCutoffDistance(force.getCutoffDistance());
            int gridSizeX, gridSizeY, gridSizeZ;
            NonbondedForceImpl::calcPMEParameters(system, nb, alphaEwald, gridSizeX, gridSizeY, gridSizeZ);
            pmeGridDimension[0] = gridSizeX;
            pmeGridDimension[1] = gridSizeY;
            pmeGridDimension[2] = gridSizeZ;
            std::cout << "Computed PME parameters for MBPolElectrostaticsForce, alphaEwald:" <<
                    alphaEwald << " pmeGrid: " <<  gridSizeX << "," <<  gridSizeY << ","<<  gridSizeZ << std::endl;
        }    
    } else {
        usePme = false;
    }
    return;
}

MBPolCpuElectrostaticsForce* CpuCalcMBPolElectrostaticsForceKernel::setupMBPolCpuElectrostaticsForce(ContextImpl& context )
{

    // mbpolCpuElectrostaticsForce is set to MBPolCpuGeneralizedKirkwoodForce if MBPolGeneralizedKirkwoodForce is present
    // mbpolCpuElectrostaticsForce is set to MBPolCpuPmeElectrostaticsForce if 'usePme' is set
    // mbpolCpuElectrostaticsForce is set to MBPolCpuElectrostaticsForce otherwise

    MBPolCpuElectrostaticsForce* mbpolCpuElectrostaticsForce = NULL;
    if( usePme ) {

         MBPolCpuPmeElectrostaticsForce* mbpolCpuPmeElectrostaticsForce = new MBPolCpuPmeElectrostaticsForce( data.threads, data.threadForce );
         mbpolCpuPmeElectrostaticsForce->setAlphaEwald( alphaEwald );
         mbpolCpuPmeElectrostaticsForce->setCutoffDistance( cutoffDistance );
         mbpolCpuPmeElectrostaticsForce->setPmeGridDimensions( pmeGridDimension );
         RealVec& box = extractBoxSize(context);
         double minAllowedSize = 1.999999*cutoffDistance;
         if (box[0] < minAllowedSize || box[1] < minAllowedSize || box[2] < minAllowedSize){
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
         }
         mbpolCpuPmeElectrostaticsForce->setPeriodicBoxSize(box);
         mbpolCpuElectrostaticsForce = static_cast<MBPolCpuElectrostaticsForce*>(mbpolCpuPmeElectrostaticsForce);

    } else {
         mbpolCpuElectrostaticsForce = new MBPolCpuElectrostaticsForce( MBPolCpuElectrostaticsForce::NoCutoff, data.threads, data.threadForce );
    }

    mbpolCpuElectrostaticsForce->setMutualInducedDipoleTargetEpsilon( mutualInducedTargetEpsilon );
    mbpolCpuElectrostaticsForce->setMaximumMutualInducedDipoleIterations( mutualInducedMaxIterations );

    mbpolCpuElectrostaticsForce->setIncludeChargeRedistribution(includeChargeRedistribution);
    if (tholeParameters.size() > 0)
        mbpolCpuElectrostaticsForce->setTholeParameters(tholeParameters);

    return mbpolCpuElectrostaticsForce;

}

double CpuCalcMBPolElectrostaticsForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    MBPolCpuElectrostaticsForce* mbpolCpuElectrostaticsForce = setupMBPolCpuElectrostaticsForce( context );

    vector<RealVec>& posData   = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy          = mbpolCpuElectrostaticsForce->calculateForceAndEnergy( posData, charges, moleculeIndices, atomTypes, tholes,
                                                                                         dampingFactors, polarity,
                                                                                         forceData);

    delete mbpolCpuElectrostaticsForce;

    return static_cast<double>(energy);
}

void CpuCalcMBPolElectrostaticsForceKernel::getElectrostaticPotential(ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                                                        std::vector< double >& outputElectrostaticPotential ){

    MBPolCpuElectrostaticsForce* mbpolCpuElectrostaticsForce = setupMBPolCpuElectrostaticsForce( context );
    vector<RealVec>& posData                                     = extractPositions(context);
    vector<RealVec> grid( inputGrid.size() );
    vector<RealOpenMM> potential( inputGrid.size() );
    for( unsigned int ii = 0; ii < inputGrid.size(); ii++ ){
        grid[ii] = inputGrid[ii];
    }
    mbpolCpuElectrostaticsForce->calculateElectrostaticPotential( posData, charges,moleculeIndices, atomTypes,  tholes,
                                                                    dampingFactors, polarity,
                                                                    grid, potential );

    outputElectrostaticPotential.resize( inputGrid.size() );
    for( unsigned int ii = 0; ii < inputGrid.size(); ii++ ){
        outputElectrostaticPotential[ii] = potential[ii];
    }

    delete mbpolCpuElectrostaticsForce;

    return;
}

void CpuCalcMBPolElectrostaticsForceKernel::getSystemElectrostaticsMoments(ContextImpl& context, std::vector< double >& outputElectrostaticsMoments){

    // retrieve masses

    const OpenMM::System& system             = context.getSystem();
    vector<RealOpenMM> masses;
    for (int i = 0; i <  system.getNumParticles(); ++i) {
        masses.push_back( static_cast<RealOpenMM>(system.getParticleMass(i)) );
    }    

    MBPolCpuElectrostaticsForce* mbpolCpuElectrostaticsForce = setupMBPolCpuElectrostaticsForce( context );
    vector<RealVec>& posData                                     = extractPositions(context);
    mbpolCpuElectrostaticsForce->calculateMBPolSystemElectrostaticsMoments( masses, posData, charges, moleculeIndices, atomTypes, tholes,
                                                                          dampingFactors, polarity,
                                                                          outputElectrostaticsMoments );

    delete mbpolCpuElectrostaticsForce;

    return;
}

void CpuCalcMBPolElectrostaticsForceKernel::copyParametersToContext(ContextImpl& context, const MBPolElectrostaticsForce& force) {
    if (numElectrostatics != force.getNumElectrostatics())
        throw OpenMMException("updateParametersInContext: The number of multipoles has changed");

    // Record the values.

    int tholeIndex = 0;
    for (int i = 0; i < numElectrostatics; ++i) {
        double charge, dampingFactorD, polarityD;
        std::vector<double> tholeD;
		int moleculeIndex, atomType;
        force.getElectrostaticsParameters(i, charge, moleculeIndex, atomType, dampingFactorD, polarityD);
        moleculeIndices[i] = moleculeIndex;
        atomTypes[i] = atomType;
        charges[i] = (RealOpenMM) charge;
        dampingFactors[i] = (RealOpenMM) dampingFactorD;
        polarity[i] = (RealOpenMM) polarityD;
    }
}
