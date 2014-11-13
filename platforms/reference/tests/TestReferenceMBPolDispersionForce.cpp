/* -------------------------------------------------------------------------- *
 *                                   OpenMMMBPol                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs                                                   *
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
#include "openmm/Context.h"
#include "OpenMMMBPol.h"
#include "openmm/System.h"

#include "openmm/MBPolDispersionForce.h"

#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#define ASSERT_EQUAL_TOL_MOD(expected, found, tol, testname) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC_MOD(expected, found, tol,testname) {ASSERT_EQUAL_TOL_MOD((expected)[0], (found)[0], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[1], (found)[1], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[2], (found)[2], (tol),(testname));};

using namespace  OpenMM;
using namespace MBPolPlugin;

const double TOL = 1e-4;

void testDispersion( double boxDimension, bool addPositionOffset ) {

    std::string testName      = "testMBPolDispersionForceTest";

    System system;
    // Dispersion Force
    MBPolDispersionForce* dispersionForce = new MBPolDispersionForce();
    dispersionForce->setCutoff( 10 );

    if( boxDimension > 0.0 ){
        Vec3 a( boxDimension, 0.0, 0.0 );
        Vec3 b( 0.0, boxDimension, 0.0 );
        Vec3 c( 0.0, 0.0, boxDimension );
        system.setDefaultPeriodicBoxVectors( a, b, c );
        dispersionForce->setNonbondedMethod(MBPolDispersionForce::CutoffPeriodic);
        dispersionForce->setUseDispersionCorrection(true);
    } else {
        dispersionForce->setNonbondedMethod(MBPolDispersionForce::CutoffNonPeriodic);
    }

    int numberOfWaterMolecules = 3;
    unsigned int particlesPerMolecule = 3;
    int numberOfParticles          = numberOfWaterMolecules * particlesPerMolecule;

    for( unsigned int jj = 0; jj < numberOfParticles; jj += particlesPerMolecule ){
        system.addParticle( 1.5999000e+01 );
        dispersionForce->addParticle( "O");

        system.addParticle( 1.0080000e+00 );
        dispersionForce->addParticle( "H");

        system.addParticle( 1.0080000e+00 );
        dispersionForce->addParticle( "H");

    }


    const double C6_HH = 2.009358600184719e+01; // kcal/mol * A^(-6)
    const double C6_OH = 8.349556669872743e+01; // kcal/mol * A^(-6)
    const double C6_OO = 2.373212214147944e+02; // kcal/mol * A^(-6)

    const double d6_HH =  9.406475169954112e+00; // A^(-1)
    const double d6_OH =  9.775202425217957e+00; // A^(-1)
    const double d6_OO =  9.295485815062264e+00; // A^(-1)

    // <!-- Units: c6 [kJ mol^{-1} nm^{-6}], d6 [nm^{-1}] -->
    dispersionForce->addDispersionParameters("O", "O", 9.92951990e+08, 9.29548582e+01);
    dispersionForce->addDispersionParameters("O", "H", 3.49345451e+08, 9.77520243e+01);
    dispersionForce->addDispersionParameters("H", "H", 8.40715638e+07, 9.40647517e+01);

    system.addForce(dispersionForce);

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    std::vector<Vec3> positions(numberOfParticles);
    std::vector<Vec3> expectedForces(numberOfParticles);
    double expectedEnergy;

    positions[0]             = Vec3( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  );
    positions[1]             = Vec3( -6.218989773e-01, -6.009430735e-01,  1.572437625e+00  );
    positions[2]             = Vec3( -2.017613812e+00, -4.190350349e-01,  2.239642849e+00  );

    positions[3]             = Vec3( -1.763651687e+00, -3.816594649e-01, -1.300353949e+00  );
    positions[4]             = Vec3( -1.903851736e+00, -4.935677617e-01, -3.457810126e-01  );
    positions[5]             = Vec3( -2.527904158e+00, -7.613550077e-01, -1.733803676e+00  );

    positions[6]             = Vec3( -5.588472140e-01,  2.006699172e+00, -1.392786582e-01  );
    positions[7]             = Vec3( -9.411558180e-01,  1.541226676e+00,  6.163293071e-01  );
    positions[8]            = Vec3( -9.858551734e-01,  1.567124294e+00, -8.830970941e-01  );

    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            positions[i][j] *= 1e-1;
        }
    }

    if (addPositionOffset) {
        // move second molecule 1 box dimension in Y direction
        positions[3][1] += boxDimension;
        positions[4][1] += boxDimension;
        positions[5][1] += boxDimension;
    }

    expectedForces[0]     = Vec3( 0.08571609, -2.87858990,  9.55466024 );
    expectedForces[1]     = Vec3( 0.30171745, -0.45439770,  0.69698990 );
    expectedForces[2]     = Vec3(-0.07264503, -0.10500907,  0.33143112 );
    expectedForces[3]     = Vec3(-1.47718406, -3.16347777, -2.38624713 );
    expectedForces[4]     = Vec3(-1.88830235, -1.76186242, -6.63181710 );
    expectedForces[5]     = Vec3(-0.10127396, -0.11424219, -0.15574083 );
    expectedForces[6]     = Vec3( 0.81603792,  1.80760267, -0.23906766 );
    expectedForces[7]     = Vec3( 1.23858368,  3.92938304, -1.37060120 );
    expectedForces[8]     = Vec3( 1.09735028,  2.74059334,  0.20039267 );

    // gradients => forces
    for( unsigned int ii = 0; ii < expectedForces.size(); ii++ ){
        expectedForces[ii] *= -1;
    }
    expectedEnergy        = -6.84471477;

    std::string platformName;
    #define AngstromToNm 0.1    
    #define CalToJoule   4.184

    platformName = "Reference";
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);
    context.applyConstraints(1e-6); // update position of virtual sites

    State state                      = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces         = state.getForces();

    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        forces[ii] /= CalToJoule*10;
    }    

    double tolerance = 1.0e-03;


    double energy = state.getPotentialEnergy() / CalToJoule;

    std::cout << "Energy: " << energy << " Kcal/mol "<< std::endl;
    std::cout << "Expected energy: " << expectedEnergy << " Kcal/mol "<< std::endl;

    std::cout  << std::endl << "Forces:" << std::endl;

    for (int i=0; i<numberOfParticles; i++) {
           std::cout << "Force atom " << i << ": " << expectedForces[i] << " Kcal/mol/A <mbpol>" << std::endl;
           std::cout << "Force atom " << i << ": " << forces[i] << " Kcal/mol/A <openmm-mbpol>" << std::endl << std::endl;
       }

       std::cout << "Comparison of energy and forces with tolerance: " << tolerance << std::endl << std::endl;



   ASSERT_EQUAL_TOL( expectedEnergy, energy, tolerance );

   for( unsigned int ii = 0; ii < forces.size(); ii++ ){
       ASSERT_EQUAL_VEC( expectedForces[ii], forces[ii], tolerance );
   }
   std::cout << "Test Successful: " << testName << std::endl << std::endl;

}

int main( int numberOfArguments, char* argv[] ) {

    try {
        FILE* log = NULL;

        std::cout << "TestReferenceMBPolDispersionForce" << std::endl;
        double boxDimension = 0;
        testDispersion( boxDimension, false );

        std::cout << "TestReferenceMBPolDispersionForce Periodic boundary conditions" << std::endl;
        boxDimension = 50;
        testDispersion( boxDimension, false );

        std::cout << "TestReferenceMBPolDispersionForce Periodic boundary conditions with offset" << std::endl;
        boxDimension = 50;
        testDispersion( boxDimension, true );

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
