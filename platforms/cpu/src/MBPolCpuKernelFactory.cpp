/* -------------------------------------------------------------------------- *
 *                              OpenMMMBPol                                  *
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

#include "MBPolCpuKernelFactory.h"
#include "MBPolCpuKernels.h"
#include "openmm/cpu/CpuPlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"
#include <iostream>

using namespace  OpenMM;
using namespace MBPolPlugin;

extern "C" OPENMM_EXPORT void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
    std::cout << "Executing registerKernelFactories of CPU MBPol" << std::endl;
    for( int ii = 0; ii < Platform::getNumPlatforms(); ii++ ){
        Platform& platform = Platform::getPlatform(ii);
        if( platform.getName() == "CPU" ){

             MBPolCpuKernelFactory* factory = new MBPolCpuKernelFactory();

             //platform.registerKernelFactory(CalcMBPolOneBodyForceKernel::Name(),           factory);
             //platform.registerKernelFactory(CalcMBPolTwoBodyForceKernel::Name(),                   factory);
             //platform.registerKernelFactory(CalcMBPolThreeBodyForceKernel::Name(),                   factory);
             platform.registerKernelFactory(CalcMBPolElectrostaticsForceKernel::Name(),             factory);
             std::cout << "Registered CPU CalcMBPolElectrostaticsForceKernel::Name" << std::endl;
        }
    }
}

extern "C" OPENMM_EXPORT void registerMBPolCpuKernelFactories() {
    try {
        Platform::getPlatformByName("CPU");
    }
    catch (...) {
        Platform::registerPlatform(new CpuPlatform());
    }
    registerKernelFactories();
}

KernelImpl* MBPolCpuKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    CpuPlatform::PlatformData& data = CpuPlatform::getPlatformData(context);

    // create MBPolCpuData object if contextToMBPolDataMap does not contain
    // key equal to current context
    //if (name == CalcMBPolOneBodyForceKernel::Name())
    //    return new ReferenceCalcMBPolOneBodyForceKernel(name, platform, context.getSystem());

    //if (name == CalcMBPolTwoBodyForceKernel::Name())
    //    return new ReferenceCalcMBPolTwoBodyForceKernel(name, platform, context.getSystem());

    //if (name == CalcMBPolThreeBodyForceKernel::Name())
    //        return new ReferenceCalcMBPolThreeBodyForceKernel(name, platform, context.getSystem());

    if (name == CalcMBPolElectrostaticsForceKernel::Name())
        return new CpuCalcMBPolElectrostaticsForceKernel(name, platform, data);

    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
