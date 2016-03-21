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

/**
 * This tests the Reference implementation of ReferenceMBPolThreeBodyForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMMBPol.h"
#include "openmm/System.h"
#include "openmm/MBPolTwoBodyForce.h"
#include "MBPolReferenceElectrostaticsForce.h"
#include "openmm/VirtualSite.h"

#include "openmm/MBPolThreeBodyForce.h"
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

void runTest( double boxDimension ) {

    std::string testName      = "testMBPolIntegrationTest";

    System system;

    double virtualSiteWeightO = 0.573293118;
    double virtualSiteWeightH = 0.213353441;
    MBPolElectrostaticsForce* mbpolElectrostaticsForce        = new MBPolElectrostaticsForce();;
    mbpolElectrostaticsForce->setForceGroup(0);

    std::vector<double> zeroDipole(3);
    std::vector<double> zeroQuadrupole(9);
    std::vector<double> thole(5);

    std::fill(zeroDipole.begin(), zeroDipole.end(), 0.);
    std::fill(zeroQuadrupole.begin(), zeroQuadrupole.end(), 0.);

    thole[TCC] = 0.4;
    thole[TCD] = 0.4;
    thole[TDD] = 0.055;
    thole[TDDOH]  = 0.626;
    thole[TDDHH] = 0.055;

    // One body interaction
    MBPolOneBodyForce* mbpolOneBodyForce = new MBPolOneBodyForce();
    mbpolOneBodyForce->setForceGroup(1);

    // Two body interaction
    MBPolTwoBodyForce* mbpolTwoBodyForce = new MBPolTwoBodyForce();
    mbpolTwoBodyForce->setForceGroup(2);
    double cutoff = .9;
    mbpolTwoBodyForce->setCutoff( cutoff );

    // Three body interaction
    MBPolThreeBodyForce* mbpolThreeBodyForce = new MBPolThreeBodyForce();
    mbpolThreeBodyForce->setCutoff( cutoff );
    mbpolThreeBodyForce->setForceGroup(3);

    // Dispersion Force
    MBPolDispersionForce* dispersionForce = new MBPolDispersionForce();
    dispersionForce->setCutoff( cutoff );
    dispersionForce->setForceGroup(4);

    if( boxDimension > 0.0 ){
        Vec3 a( boxDimension, 0.0, 0.0 );
        Vec3 b( 0.0, boxDimension, 0.0 );
        Vec3 c( 0.0, 0.0, boxDimension );
        system.setDefaultPeriodicBoxVectors( a, b, c );

        mbpolElectrostaticsForce->setNonbondedMethod( MBPolElectrostaticsForce::PME );
        double ewaldErrorTol = 1e-3;
        mbpolElectrostaticsForce->setEwaldErrorTolerance(ewaldErrorTol);

        mbpolTwoBodyForce->setNonbondedMethod(MBPolTwoBodyForce::CutoffPeriodic);

        mbpolThreeBodyForce->setNonbondedMethod(MBPolThreeBodyForce::CutoffPeriodic);

        dispersionForce->setNonbondedMethod(MBPolDispersionForce::CutoffPeriodic);
        dispersionForce->setUseDispersionCorrection(true);


    } else {
        mbpolElectrostaticsForce->setNonbondedMethod( MBPolElectrostaticsForce::NoCutoff );

        mbpolTwoBodyForce->setNonbondedMethod(MBPolTwoBodyForce::CutoffNonPeriodic);

        mbpolThreeBodyForce->setNonbondedMethod(MBPolThreeBodyForce::CutoffNonPeriodic);

        dispersionForce->setNonbondedMethod(MBPolDispersionForce::CutoffNonPeriodic);

    }

    int numberOfWaterMolecules = 3;
    unsigned int particlesPerMolecule = 4;
    int numberOfParticles          = numberOfWaterMolecules * particlesPerMolecule;

    std::vector<int> particleIndices(particlesPerMolecule);
	int waterMoleculeIndex=0;
    for( unsigned int jj = 0; jj < numberOfParticles; jj += particlesPerMolecule ){
        system.addParticle( 1.5999000e+01 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 0. ); // Virtual Site
        system.setVirtualSite(jj+3, new ThreeParticleAverageSite(jj, jj+1, jj+2,
                                                                   virtualSiteWeightO, virtualSiteWeightH,virtualSiteWeightH));


        particleIndices[0] = jj;
        particleIndices[1] = jj+1;
        particleIndices[2] = jj+2;
        mbpolElectrostaticsForce->addElectrostatics( -5.1966000e-01,
                                            waterMoleculeIndex, 0, 0.001310, 0.001310 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01,
                                            waterMoleculeIndex, 1, 0.000294, 0.000294 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01,
                                            waterMoleculeIndex, 1, 0.000294, 0.000294 );
        mbpolElectrostaticsForce->addElectrostatics(  0.,
                                            waterMoleculeIndex, 2, 0.001310,  0.);
		waterMoleculeIndex++;

        mbpolOneBodyForce->addOneBody(particleIndices);
        mbpolTwoBodyForce->addParticle( particleIndices);
        mbpolThreeBodyForce->addParticle( particleIndices);
        dispersionForce->addParticle( "O");
        dispersionForce->addParticle( "H");
        dispersionForce->addParticle( "H");
        dispersionForce->addParticle( "M");

    }

    // <!-- Units: c6 [kJ mol^{-1} nm^{-6}], d6 [nm^{-1}] -->
    dispersionForce->addDispersionParameters("O", "O", 9.92951990e+08, 9.29548582e+01);
    dispersionForce->addDispersionParameters("O", "H", 3.49345451e+08, 9.77520243e+01);
    dispersionForce->addDispersionParameters("H", "H", 8.40715638e+07, 9.40647517e+01);

    system.addForce(mbpolElectrostaticsForce);
    system.addForce(mbpolOneBodyForce);
    system.addForce(mbpolTwoBodyForce);
    system.addForce(mbpolThreeBodyForce);
    system.addForce(dispersionForce);

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    std::vector<Vec3> positions(numberOfParticles);
    std::vector<Vec3> expectedForces(numberOfParticles);
    double expectedEnergy;

    positions[0]             = Vec3( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  );
    positions[1]             = Vec3( -6.218989773e-01, -6.009430735e-01,  1.572437625e+00  );
    positions[2]             = Vec3( -2.017613812e+00, -4.190350349e-01,  2.239642849e+00  );
    positions[3]             = Vec3( 0,0,0  );

    positions[4]             = Vec3( -1.763651687e+00, -3.816594649e-01, -1.300353949e+00  );
    positions[5]             = Vec3( -1.903851736e+00, -4.935677617e-01, -3.457810126e-01  );
    positions[6]             = Vec3( -2.527904158e+00, -7.613550077e-01, -1.733803676e+00  );
    positions[7]             = Vec3( 0,0,0  );

    positions[8]             = Vec3( -5.588472140e-01,  2.006699172e+00, -1.392786582e-01  );
    positions[9]             = Vec3( -9.411558180e-01,  1.541226676e+00,  6.163293071e-01  );
    positions[10]            = Vec3( -9.858551734e-01,  1.567124294e+00, -8.830970941e-01  );
    positions[11]            = Vec3( 0,0,0  );

    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            positions[i][j] *= 1e-1;
        }
    }

    expectedForces[0]     = Vec3( -25.82906782,  11.96555757,  -3.59633570);
    expectedForces[1]     = Vec3(  26.99237581, -11.11111590,   0.85799351);
    expectedForces[2]     = Vec3(  -1.36390020,   1.80325388,  -2.80944334);
    expectedForces[4]     = Vec3(   1.51466430,   1.05236599,  -4.05062141);
    expectedForces[5]     = Vec3(   0.80024930,  -0.23872567,   9.89802544);
    expectedForces[6]     = Vec3(  -0.11705181,  -0.00181039,  -1.35555449);
    expectedForces[8]     = Vec3(   4.33262397,   1.50639482,  -4.28267317);
    expectedForces[9]     = Vec3(  -3.35396424,  -5.45992971,   6.41676701);
    expectedForces[10]    = Vec3(  -2.97592932,   0.48400940,  -1.07815786);

    // gradients => forces
    for( unsigned int ii = 0; ii < expectedForces.size(); ii++ ){
        expectedForces[ii] *= -1;
    }

    if (boxDimension == 1.8)
    {
        expectedEnergy        = -8.92353;
    } else { // cluster or large box
        expectedEnergy        = -8.78893485;
    }

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
        if ((ii+1) % 4 == 0) { // Set virtual site force to 0
            forces[ii] *= 0;
        }
    }    

    double tolerance = 1.0e-02;

    vector<string> forceLabels;
    forceLabels.push_back("Electrostatics");
    forceLabels.push_back("OneBody");
    forceLabels.push_back("TwoBody");
    forceLabels.push_back("ThreeBody");
    forceLabels.push_back("Dispersion");
    double energy = state.getPotentialEnergy() / CalToJoule;

    std::cout << "Total Energy: " << energy << " Kcal/mol "<< std::endl;
    std::cout << "Expected energy: " << expectedEnergy << " Kcal/mol "<< std::endl;

    for ( unsigned int ii = 0; ii < 5; ii++ ){
        state                      = context.getState(State::Energy, false, pow(2, ii));
        std::cout << "Energy " << forceLabels[ii] << ": " << state.getPotentialEnergy() / CalToJoule << " Kcal/mol "<< std::endl;
    }

//    std::cout  << std::endl << "Forces:" << std::endl;
//
//    for (int i=0; i<numberOfParticles; i++) {
//           std::cout << "Force atom " << i << ": " << expectedForces[i] << " Kcal/mol/A <mbpol>" << std::endl;
//           std::cout << "Force atom " << i << ": " << forces[i] << " Kcal/mol/A <openmm-mbpol>" << std::endl << std::endl;
//       }
//
//       std::cout << "Comparison of energy and forces with tolerance: " << tolerance << std::endl << std::endl;
//
//
//
     ASSERT_EQUAL_TOL( expectedEnergy, energy, tolerance );
//
//   for( unsigned int ii = 0; ii < forces.size(); ii++ ){
//       ASSERT_EQUAL_VEC( expectedForces[ii], forces[ii], tolerance );
//   }
   std::cout << "Test Successful: " << testName << std::endl << std::endl;

}

int main( int numberOfArguments, char* argv[] ) {

    try {
        std::cout << "TestReferenceMBPolIntegrationTest running test..." << std::endl;
        double boxDimension = 0;
        runTest( boxDimension );

        boxDimension = 1.8;
        runTest( boxDimension );

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
