/* -------------------------------------------------------------------------- *
 *                                OpenMMMBPol                                *
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

#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "openmm/MBPolElectrostaticsForce.h"
#include "openmm/internal/MBPolElectrostaticsForceImpl.h"
#include <stdio.h>

using namespace  OpenMM;
using namespace MBPolPlugin;
using std::string;
using std::vector;

MBPolElectrostaticsForce::MBPolElectrostaticsForce() : nonbondedMethod(NoCutoff), pmeBSplineOrder(5), cutoffDistance(0.9), ewaldErrorTol(1e-4), mutualInducedMaxIterations(200),
                                               mutualInducedTargetEpsilon(1.0e-07), scalingDistanceCutoff(100.0), electricConstant(138.9354558456), aewald(0.0), includeChargeRedistribution(true) {
    pmeGridDimension.resize(3);
    pmeGridDimension[0] = pmeGridDimension[1] = pmeGridDimension[2];
    const double defaultTholeParameters[5] = { 0.4, 0.4, 0.055, 0.626, 0.055 };
    tholeParameters.assign(&defaultTholeParameters[0], &defaultTholeParameters[0]+5);
}

MBPolElectrostaticsForce::NonbondedMethod MBPolElectrostaticsForce::getNonbondedMethod( void ) const {
    return nonbondedMethod;
}

void MBPolElectrostaticsForce::setNonbondedMethod( MBPolElectrostaticsForce::NonbondedMethod method) {
    nonbondedMethod = method;
}

double MBPolElectrostaticsForce::getCutoffDistance( void ) const {
    return cutoffDistance;
}

void MBPolElectrostaticsForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

double MBPolElectrostaticsForce::getAEwald() const { 
    return aewald; 
} 

void MBPolElectrostaticsForce::setIncludeChargeRedistribution( bool chargeRedistribution ) {
    includeChargeRedistribution = chargeRedistribution;
}

bool MBPolElectrostaticsForce::getIncludeChargeRedistribution( void ) const {
    return includeChargeRedistribution;
}
void MBPolElectrostaticsForce::setAEwald(double inputAewald ) { 
    aewald = inputAewald; 
} 
 
int MBPolElectrostaticsForce::getPmeBSplineOrder( void ) const { 
    return pmeBSplineOrder; 
} 
 
void MBPolElectrostaticsForce::getPmeGridDimensions( std::vector<int>& gridDimension ) const { 
    if( gridDimension.size() < 3 ){
        gridDimension.resize(3);
    }
    if( pmeGridDimension.size() > 2 ){
        gridDimension[0] = pmeGridDimension[0];
        gridDimension[1] = pmeGridDimension[1];
        gridDimension[2] = pmeGridDimension[2];
    } else {
        gridDimension[0] = gridDimension[1] = gridDimension[2] = 0;
    }
    return;
} 
 
void MBPolElectrostaticsForce::setPmeGridDimensions( const std::vector<int>& gridDimension ) {
    pmeGridDimension.resize(3);
    pmeGridDimension[0] = gridDimension[0];
    pmeGridDimension[1] = gridDimension[1];
    pmeGridDimension[2] = gridDimension[2];
    return;
} 

int MBPolElectrostaticsForce::getMutualInducedMaxIterations( void ) const {
    return mutualInducedMaxIterations;
}

void MBPolElectrostaticsForce::setMutualInducedMaxIterations( int inputMutualInducedMaxIterations ) {
    mutualInducedMaxIterations = inputMutualInducedMaxIterations;
}

double MBPolElectrostaticsForce::getMutualInducedTargetEpsilon( void ) const {
    return mutualInducedTargetEpsilon;
}

void MBPolElectrostaticsForce::setMutualInducedTargetEpsilon( double inputMutualInducedTargetEpsilon ) {
    mutualInducedTargetEpsilon = inputMutualInducedTargetEpsilon;
}

double MBPolElectrostaticsForce::getEwaldErrorTolerance() const {
    return ewaldErrorTol;
}

void MBPolElectrostaticsForce::setEwaldErrorTolerance(double tol) {
    ewaldErrorTol = tol;
}

int MBPolElectrostaticsForce::addElectrostatics( double charge,
                                       int moleculeIndex, int atomType, double dampingFactor, double polarity) {
    multipoles.push_back(ElectrostaticsInfo( charge, moleculeIndex, atomType, dampingFactor, polarity));
    return multipoles.size()-1;
}

void MBPolElectrostaticsForce::getElectrostaticsParameters(int index, double& charge,
                                                  int& moleculeIndex, int& atomType, double& dampingFactor, double& polarity ) const {
    charge                      = multipoles[index].charge;

    moleculeIndex               = multipoles[index].moleculeIndex;
    atomType                    = multipoles[index].atomType;
    dampingFactor               = multipoles[index].dampingFactor;
    polarity                    = multipoles[index].polarity;
}

void MBPolElectrostaticsForce::setElectrostaticsParameters(int index, double charge,
                                                  int moleculeIndex, int atomType, double dampingFactor, double polarity ) {

    multipoles[index].charge                      = charge;

    multipoles[index].dampingFactor               = dampingFactor;
    multipoles[index].polarity                    = polarity;
    multipoles[index].moleculeIndex = moleculeIndex;
    multipoles[index].atomType = atomType;

}

void MBPolElectrostaticsForce::getElectrostaticPotential( const std::vector< Vec3 >& inputGrid, Context& context, std::vector< double >& outputElectrostaticPotential ){
    dynamic_cast<MBPolElectrostaticsForceImpl&>(getImplInContext(context)).getElectrostaticPotential(getContextImpl(context), inputGrid, outputElectrostaticPotential);
}

ForceImpl* MBPolElectrostaticsForce::createImpl()  const {
    return new MBPolElectrostaticsForceImpl(*this);
}

void MBPolElectrostaticsForce::updateParametersInContext(Context& context) {
    dynamic_cast<MBPolElectrostaticsForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
