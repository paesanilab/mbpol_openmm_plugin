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

/**
 * This tests the Reference implementation of ReferenceMBPolOneBodyForce.
 */

#include "openmm/internal/AssertionUtilities.h"
//#include "MBPolTinkerParameterFile.h"
const double DegreesToRadians = 3.14159265/180.0;
#include "openmm/Context.h"
#include "OpenMMMBPol.h"
#include "openmm/System.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>

using namespace  OpenMM;
using namespace MBPolPlugin;

const double TOL = 1e-4;
#define PI_M               3.141592653589
#define RADIAN            57.29577951308

void testOneOneBody( FILE* log ) {

    System system;
    int numberOfParticles = 3;
    for( int ii = 0; ii < numberOfParticles; ii++ ){
        system.addParticle(1.0);
    }

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    MBPolOneBodyForce* mbpolOneBodyForce = new MBPolOneBodyForce();

    double abLength         = 0.144800000E+01;
    double cbLength         = 0.101500000E+01;
    double angleOneBody = 0.108500000E+03*DegreesToRadians;
    //double kOneBody     = 0.750491578E-01;
    double kOneBody     = 1.0;

    std::vector<int> particleIndices;
    particleIndices.push_back(0);
    particleIndices.push_back(1);
    particleIndices.push_back(2);

    mbpolOneBodyForce->addOneBody(particleIndices);

    system.addForce(mbpolOneBodyForce);
    Context context(system, integrator, Platform::getPlatformByName( "Reference"));

    std::vector<Vec3> positions(numberOfParticles);

    positions[0]             = Vec3( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  );
    positions[1]             = Vec3( -6.218989773e-01, -6.009430735e-01,  1.572437625e+00  );
    positions[2]             = Vec3( -2.017613812e+00, -4.190350349e-01,  2.239642849e+00  );

    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            positions[i][j] *= 1e-1;
        }
    }

    context.setPositions(positions);
    
    State state                      = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces   = state.getForces();
    double cal2joule = 4.184;

    std::vector<Vec3> expectedForces(numberOfParticles);

    expectedForces[0]         = Vec3(-27.48162433,     8.92495995,   2.80995323 );
    expectedForces[1]         = Vec3( 30.78909844,   -11.48714187,  -0.27204770 );
    expectedForces[2]         = Vec3( -3.30747410,     2.56218193,  -2.53790553 );

    // MBPol gives the gradients, we use forces in OpenMM, need to flip sign
    for (int i=0; i<numberOfParticles; i++) {
            for (int j=0; j<3; j++) {
                expectedForces[i][j] *= -1;
            }
        }

    double tolerance=1e-4;
    double expectedEnergy = 0.55975882;

    std::cout  << std::endl << "Forces:" << std::endl;

    for (int i=0; i<numberOfParticles; i++) {
           for (int j=0; j<3; j++) {
            forces[i][j] /= cal2joule*10;
           }
       }

    for (int i=0; i<numberOfParticles; i++) {
         std::cout << forces[i] << " Kcal/mol/A " << std::endl;
    }

    std::cout  << std::endl << "Expected forces:" << std::endl;

    for (int i=0; i<numberOfParticles; i++) {
         std::cout << expectedForces[i] << " Kcal/mol/A " << std::endl;
    }

    double energy = state.getPotentialEnergy();

    std::cout << "Energy: " << energy/cal2joule << " Kcal/mol "<< std::endl;
    std::cout << "Expected energy: " << expectedEnergy << " Kcal/mol "<< std::endl;

    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        ASSERT_EQUAL_VEC( expectedForces[ii], forces[ii], tolerance );
    }
    ASSERT_EQUAL_TOL( expectedEnergy, energy/cal2joule, tolerance );

}

int main( int numberOfArguments, char* argv[] ) {

    try {
        std::cout << "TestReferenceMBPolOneBodyForce running test..." << std::endl;

        FILE* log = NULL;
        //FILE* log = stderr;
        //FILE* log = fopen( "MBPolOneBodyForce1.log", "w" );;
        testOneOneBody( log );
#ifdef MBPOL_DEBUG
        if( log && log != stderr )
            (void) fclose( log );
#endif

    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    //std::cout << "PASS - Test succeeded." << std::endl;
    std::cout << "Done" << std::endl;
    return 0;
}
