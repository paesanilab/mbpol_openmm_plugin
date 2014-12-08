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
#include "openmm/MBPolDispersionForce.h"
#include "openmm/internal/MBPolDispersionForceImpl.h"

using namespace  OpenMM;
using namespace MBPolPlugin;
using std::string;
using std::vector;

using std::map;
using std::make_pair;
using std::pair;


MBPolDispersionForce::MBPolDispersionForce() : nonbondedMethod(CutoffNonPeriodic), cutoff(1.0e+10) {
}

int MBPolDispersionForce::addParticle(string atomElement ) {
    parameters.push_back(DispersionInfo(atomElement));
    return parameters.size()-1;
}

void MBPolDispersionForce::addDispersionParameters(string firstElement, string secondElement, double c6, double d6){
    c6d6Data[make_pair(firstElement, secondElement)] = make_pair(c6, d6);
}

c6d6Datatype MBPolDispersionForce::getDispersionParameters( void ) const{
    return c6d6Data;
}

int MBPolDispersionForce::getNumMolecules() const {
    return parameters.size();
}

void MBPolDispersionForce::getParticleParameters(int particleIndex, string & atomElement ) const {
    atomElement = parameters[particleIndex].atomElement;
}

void MBPolDispersionForce::setParticleParameters(int particleIndex, string atomElement  ) {
      parameters[particleIndex].atomElement = atomElement;
}

void MBPolDispersionForce::setCutoff( double inputCutoff ){
    cutoff = inputCutoff;
}

double MBPolDispersionForce::getCutoff( void ) const {
    return cutoff;
}

MBPolDispersionForce::NonbondedMethod MBPolDispersionForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void MBPolDispersionForce::setNonbondedMethod(NonbondedMethod method) {
    nonbondedMethod = method;
}

ForceImpl* MBPolDispersionForce::createImpl() const {
    return new MBPolDispersionForceImpl(*this);
}

void MBPolDispersionForce::updateParametersInContext(Context& context) {
    dynamic_cast<MBPolDispersionForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
