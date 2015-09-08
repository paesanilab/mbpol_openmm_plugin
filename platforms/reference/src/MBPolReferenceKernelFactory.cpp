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

#include "MBPolReferenceKernelFactory.h"
#include "MBPolReferenceKernels.h"
#include "openmm/reference/ReferencePlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace  OpenMM;
using namespace MBPolPlugin;

#if defined(WIN32)
    #include <windows.h>
    extern "C" void initMBPolReferenceKernels();
    BOOL WINAPI DllMain(HANDLE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved) {
        if (ul_reason_for_call == DLL_PROCESS_ATTACH)
            initMBPolReferenceKernels();
        return TRUE;
    }
#else
    extern "C" void __attribute__((constructor)) initMBPolReferenceKernels();
#endif

extern "C" void initMBPolReferenceKernels() {
    for( int ii = 0; ii < Platform::getNumPlatforms(); ii++ ){
        Platform& platform = Platform::getPlatform(ii);
        if( platform.getName() == "Reference" ){

             MBPolReferenceKernelFactory* factory = new MBPolReferenceKernelFactory();

             platform.registerKernelFactory(CalcMBPolOneBodyForceKernel::Name(),           factory);
             platform.registerKernelFactory(CalcMBPolTwoBodyForceKernel::Name(),                   factory);
             platform.registerKernelFactory(CalcMBPolThreeBodyForceKernel::Name(),                   factory);
             platform.registerKernelFactory(CalcMBPolDispersionForceKernel::Name(),                   factory);
             platform.registerKernelFactory(CalcMBPolElectrostaticsForceKernel::Name(),             factory);
        }
    }
}

KernelImpl* MBPolReferenceKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    ReferencePlatform::PlatformData& referencePlatformData = *static_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());

    // create MBPolReferenceData object if contextToMBPolDataMap does not contain
    // key equal to current context
    if (name == CalcMBPolOneBodyForceKernel::Name())
        return new ReferenceCalcMBPolOneBodyForceKernel(name, platform, context.getSystem());

    if (name == CalcMBPolTwoBodyForceKernel::Name())
        return new ReferenceCalcMBPolTwoBodyForceKernel(name, platform, context.getSystem());

    if (name == CalcMBPolThreeBodyForceKernel::Name())
            return new ReferenceCalcMBPolThreeBodyForceKernel(name, platform, context.getSystem());

    if (name == CalcMBPolElectrostaticsForceKernel::Name())
        return new ReferenceCalcMBPolElectrostaticsForceKernel(name, platform, context.getSystem());

    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
