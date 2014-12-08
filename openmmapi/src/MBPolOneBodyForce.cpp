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
#include "openmm/MBPolOneBodyForce.h"
#include "openmm/internal/MBPolOneBodyForceImpl.h"

using namespace  OpenMM;
using namespace MBPolPlugin;

MBPolOneBodyForce::MBPolOneBodyForce() {
}

int MBPolOneBodyForce::addOneBody(const std::vector<int> & particleIndices    ) {
    stretchBends.push_back(OneBodyInfo(particleIndices));
    return stretchBends.size()-1;
}

void MBPolOneBodyForce::getOneBodyParameters(int particleIndex, std::vector<int>& particleIndices ) const {
    particleIndices     = stretchBends[particleIndex].particleIndices;
}
void MBPolOneBodyForce::setOneBodyParameters(int index, std::vector<int>& particleIndices  ) {
    stretchBends[index].particleIndices =particleIndices;
}

MBPolOneBodyForce::NonbondedMethod MBPolOneBodyForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void MBPolOneBodyForce::setNonbondedMethod(NonbondedMethod method) {
    nonbondedMethod = method;
}

ForceImpl* MBPolOneBodyForce::createImpl() const {
    return new MBPolOneBodyForceImpl(*this);
}

void MBPolOneBodyForce::updateParametersInContext(Context& context) {
    dynamic_cast<MBPolOneBodyForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
