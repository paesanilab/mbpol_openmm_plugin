/* -------------------------------------------------------------------------- *
 *                               OpenMMMBPol                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/MBPolDispersionForceImpl.h"
#include "openmm/internal/MBPolConstants.h"
#include "openmm/mbpolKernels.h"
#include <map>
#include <cmath>

using namespace  OpenMM;
using namespace MBPolPlugin;
using namespace std;


using std::pair;
using std::make_pair;
using std::vector;
using std::set;


double tang_toennies(const int n, const double& x)
{
    assert(n >= 0);

    int nn = n;

    double sum = 1.0 + x/nn;
    while (--nn != 0)
        sum = 1.0 + sum*x/nn;

    double tt = 1.0 - sum*std::exp(-x);

    if (std::fabs(tt) < 1.0e-8) {

        double term(1);
        for (nn = n; nn != 0; --nn)
            term *= x/nn;

        sum = 0.0;
        for (nn = n + 1; nn < 1000; ++nn) {
            term *= x/nn;
            sum += term;

            if (std::fabs(term/sum) < 1.0e-8)
                break;
        }

        tt = sum*std::exp(-x);
    }

    return tt;
}

double tang_toennies_long_range(const double& x)
{
    const int n = 6;
    const double if6 = 1.0/Factorial<6>::value;
    double tang_toennies_lr = tang_toennies(n, x) / std::pow(x, (n - 3)) + std::exp(-x)*(x*x*(3.0 + x) + 6*(1.0 + x)) * if6 ;

   tang_toennies_lr = tang_toennies_lr/(n - 3);
   return tang_toennies_lr;
}

double energy_long_range_correction(double cutoff, double c6, double d6)
{
    double el12 = -c6/kcal_permol_Aminus6_to_kJ_permol_nmminus6 * pow(d6/10., 3) * tang_toennies_long_range(d6/10.*cutoff);
    return el12;
}

MBPolDispersionForceImpl::MBPolDispersionForceImpl(const MBPolDispersionForce& owner) : owner(owner) {
}

MBPolDispersionForceImpl::~MBPolDispersionForceImpl() {
}

void MBPolDispersionForceImpl::initialize(ContextImpl& context) {
    const OpenMM::System& system = context.getSystem();

    // check that cutoff < 0.5*boxSize

    if (owner.getNonbondedMethod() == MBPolDispersionForce::CutoffPeriodic) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double cutoff = owner.getCutoff();
        if (cutoff > 0.5*boxVectors[0][0] || cutoff > 0.5*boxVectors[1][1] || cutoff > 0.5*boxVectors[2][2])
            throw OpenMMException("MBPolDispersionForce: The cutoff distance cannot be greater than half the periodic box size.");
    }   

    kernel = context.getPlatform().createKernel(CalcMBPolDispersionForceKernel::Name(), context);
    kernel.getAs<CalcMBPolDispersionForceKernel>().initialize(context.getSystem(), owner);
}

double MBPolDispersionForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcMBPolDispersionForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}


double MBPolDispersionForceImpl::calcDispersionCorrection(const System& system, const MBPolDispersionForce& force) {
    if (force.getNonbondedMethod() == MBPolDispersionForce::NoCutoff || force.getNonbondedMethod() == MBPolDispersionForce::CutoffNonPeriodic)
        return 0.0;

    // Identify all particle classes (defined by sigma and epsilon), and count the number of
    // particles in each class.

    map<string, int> classCounts;
    for (int i = 0; i < force.getNumParticles(); i++) {
        string atomElement;
        force.getParticleParameters(i, atomElement);
        if (atomElement != "M")
        {
        map<string, int>::iterator entry = classCounts.find(atomElement);
        if (entry == classCounts.end())
            classCounts[atomElement] = 1;
        else
            entry->second++;
        }
    }

    // Loop over all pairs of classes to compute the coefficient.

    double cutoff = force.getCutoff();
    double energy = 0.;

    c6d6Datatype c6d6Data = force.getDispersionParameters();

    for (map<string, int>::const_iterator class1 = classCounts.begin(); class1 != classCounts.end(); ++class1)
        for (map<string, int>::const_iterator class2 = class1; class2 != classCounts.end(); ++class2) {
            c6d6Datatype ::iterator entry = c6d6Data.find(make_pair(class1->first, class2->first));
            pair<double, double> c6d6;
            if (entry == c6d6Data.end())
                entry = c6d6Data.find(make_pair(class2->first, class1->first));
            c6d6 = entry->second;
            double count = (double) class1->second;
            count *= (double) class2->second;
            double energy_correction = energy_long_range_correction(cutoff, c6d6.first, c6d6.second);
            energy += 2 * count * M_PI * energy_correction;
        }
    return energy;
}


std::vector<std::string> MBPolDispersionForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcMBPolDispersionForceKernel::Name());
    return names;
}

void MBPolDispersionForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcMBPolDispersionForceKernel>().copyParametersToContext(context, owner);
}
