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

MBPolElectrostaticsForce::MBPolElectrostaticsForce() : nonbondedMethod(NoCutoff), polarizationType(Mutual), pmeBSplineOrder(5), cutoffDistance(0.9), ewaldErrorTol(1e-4), mutualInducedMaxIterations(200),
                                               mutualInducedTargetEpsilon(1.0e-07), scalingDistanceCutoff(100.0), electricConstant(138.9354558456), aewald(0.0), includeChargeRedistribution(true) {
    pmeGridDimension.resize(3);
    pmeGridDimension[0] = pmeGridDimension[1] = pmeGridDimension[2];
}

MBPolElectrostaticsForce::NonbondedMethod MBPolElectrostaticsForce::getNonbondedMethod( void ) const {
    return nonbondedMethod;
}

void MBPolElectrostaticsForce::setNonbondedMethod( MBPolElectrostaticsForce::NonbondedMethod method) {
    nonbondedMethod = method;
}

MBPolElectrostaticsForce::PolarizationType MBPolElectrostaticsForce::getPolarizationType( void ) const {
    return polarizationType;
}

void MBPolElectrostaticsForce::setPolarizationType( MBPolElectrostaticsForce::PolarizationType type ) {
    polarizationType = type;
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
                                       int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, const std::vector<double>& thole, double dampingFactor, double polarity) {
    std::vector<double> zeroDipole(3);
    std::vector<double> zeroQuadrupole(9);

    std::fill(zeroDipole.begin(), zeroDipole.end(), 0.);
    std::fill(zeroQuadrupole.begin(), zeroQuadrupole.end(), 0.);

    multipoles.push_back(ElectrostaticsInfo( charge, zeroDipole, zeroQuadrupole,  0, multipoleAtomZ,  multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity));
    return multipoles.size()-1;
}

void MBPolElectrostaticsForce::getElectrostaticsParameters(int index, double& charge, std::vector<double>& molecularDipole, std::vector<double>& molecularQuadrupole, 
                                                  int& axisType, int& multipoleAtomZ, int& multipoleAtomX, int& multipoleAtomY, std::vector<double>& thole, double& dampingFactor, double& polarity ) const {
    charge                      = multipoles[index].charge;

    molecularDipole.resize( 3 );
    molecularDipole[0]          = multipoles[index].molecularDipole[0];
    molecularDipole[1]          = multipoles[index].molecularDipole[1];
    molecularDipole[2]          = multipoles[index].molecularDipole[2];

    molecularQuadrupole.resize( 9 );
    molecularQuadrupole[0]      = multipoles[index].molecularQuadrupole[0];
    molecularQuadrupole[1]      = multipoles[index].molecularQuadrupole[1];
    molecularQuadrupole[2]      = multipoles[index].molecularQuadrupole[2];
    molecularQuadrupole[3]      = multipoles[index].molecularQuadrupole[3];
    molecularQuadrupole[4]      = multipoles[index].molecularQuadrupole[4];
    molecularQuadrupole[5]      = multipoles[index].molecularQuadrupole[5];
    molecularQuadrupole[6]      = multipoles[index].molecularQuadrupole[6];
    molecularQuadrupole[7]      = multipoles[index].molecularQuadrupole[7];
    molecularQuadrupole[8]      = multipoles[index].molecularQuadrupole[8];

    axisType                    = multipoles[index].axisType;
    multipoleAtomZ              = multipoles[index].multipoleAtomZ;
    multipoleAtomX              = multipoles[index].multipoleAtomX;
    multipoleAtomY              = multipoles[index].multipoleAtomY;

    thole.resize( 5 );
    for (int i=0; i<5;i++)
    thole[i] = multipoles[index].thole[i];
    dampingFactor               = multipoles[index].dampingFactor;
    polarity                    = multipoles[index].polarity;
}

void MBPolElectrostaticsForce::setElectrostaticsParameters(int index, double charge, const std::vector<double>& molecularDipole, const std::vector<double>& molecularQuadrupole, 
                                                  int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, const std::vector<double>&  thole, double dampingFactor, double polarity ) {

    multipoles[index].charge                      = charge;

    multipoles[index].molecularDipole[0]          = molecularDipole[0];
    multipoles[index].molecularDipole[1]          = molecularDipole[1];
    multipoles[index].molecularDipole[2]          = molecularDipole[2];

    multipoles[index].molecularQuadrupole[0]      = molecularQuadrupole[0];
    multipoles[index].molecularQuadrupole[1]      = molecularQuadrupole[1];
    multipoles[index].molecularQuadrupole[2]      = molecularQuadrupole[2];
    multipoles[index].molecularQuadrupole[3]      = molecularQuadrupole[3];
    multipoles[index].molecularQuadrupole[4]      = molecularQuadrupole[4];
    multipoles[index].molecularQuadrupole[5]      = molecularQuadrupole[5];
    multipoles[index].molecularQuadrupole[6]      = molecularQuadrupole[6];
    multipoles[index].molecularQuadrupole[7]      = molecularQuadrupole[7];
    multipoles[index].molecularQuadrupole[8]      = molecularQuadrupole[8];

    multipoles[index].axisType                    = axisType;
    multipoles[index].multipoleAtomZ              = multipoleAtomZ;
    multipoles[index].multipoleAtomX              = multipoleAtomX;
    multipoles[index].multipoleAtomY              = multipoleAtomY;
    for (int i; i<5;i++)
        multipoles[index].thole [i]                      = thole[i];
    multipoles[index].dampingFactor               = dampingFactor;
    multipoles[index].polarity                    = polarity;

}

void MBPolElectrostaticsForce::setCovalentMap(int index, CovalentType typeId, const std::vector<int>& covalentAtoms ) {

    std::vector<int>& covalentList = multipoles[index].covalentInfo[typeId];
    covalentList.resize( covalentAtoms.size() );
    for( unsigned int ii = 0; ii < covalentAtoms.size(); ii++ ){
       covalentList[ii] = covalentAtoms[ii];
    }
}

void MBPolElectrostaticsForce::getCovalentMap(int index, CovalentType typeId, std::vector<int>& covalentAtoms ) const {

    // load covalent atom index entries for atomId==index and covalentId==typeId into covalentAtoms

    std::vector<int> covalentList = multipoles[index].covalentInfo[typeId];
    covalentAtoms.resize( covalentList.size() );
    for( unsigned int ii = 0; ii < covalentList.size(); ii++ ){
       covalentAtoms[ii] = covalentList[ii];
    }
}

void MBPolElectrostaticsForce::getCovalentMaps(int index, std::vector< std::vector<int> >& covalentLists ) const {

    covalentLists.resize( CovalentEnd );
    for( unsigned int jj = 0; jj < CovalentEnd; jj++ ){
        std::vector<int> covalentList = multipoles[index].covalentInfo[jj];
        std::vector<int> covalentAtoms;
        covalentAtoms.resize( covalentList.size() );
        for( unsigned int ii = 0; ii < covalentList.size(); ii++ ){
           covalentAtoms[ii] = covalentList[ii];
        }
        covalentLists[jj] = covalentAtoms;
    }
}

void MBPolElectrostaticsForce::getElectrostaticPotential( const std::vector< Vec3 >& inputGrid, Context& context, std::vector< double >& outputElectrostaticPotential ){
    dynamic_cast<MBPolElectrostaticsForceImpl&>(getImplInContext(context)).getElectrostaticPotential(getContextImpl(context), inputGrid, outputElectrostaticPotential);
}

void MBPolElectrostaticsForce::getSystemElectrostaticsMoments(Context& context, std::vector< double >& outputElectrostaticsMonents ){
    dynamic_cast<MBPolElectrostaticsForceImpl&>(getImplInContext(context)).getSystemElectrostaticsMoments(getContextImpl(context), outputElectrostaticsMonents);
}

ForceImpl* MBPolElectrostaticsForce::createImpl()  const {
    return new MBPolElectrostaticsForceImpl(*this);
}

void MBPolElectrostaticsForce::updateParametersInContext(Context& context) {
    dynamic_cast<MBPolElectrostaticsForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
