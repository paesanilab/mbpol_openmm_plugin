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
    double cutoff = 0.9;
    mbpolTwoBodyForce->setCutoff( cutoff );

    // Three body interaction
    MBPolThreeBodyForce* mbpolThreeBodyForce = new MBPolThreeBodyForce();
    mbpolThreeBodyForce->setCutoff( cutoff );
    mbpolThreeBodyForce->setForceGroup(3);

    // Dispersion Force
    MBPolDispersionForce* dispersionForce = new MBPolDispersionForce();
    dispersionForce->setCutoff( cutoff );
    dispersionForce->setForceGroup(4);

    vector<string> forceLabels;
    forceLabels.push_back("Electrostatics");
    forceLabels.push_back("OneBody");
    forceLabels.push_back("TwoBody");
    forceLabels.push_back("ThreeBody");
    forceLabels.push_back("Dispersion");

    if( boxDimension > 0.0 ){
        Vec3 a( boxDimension, 0.0, 0.0 );
        Vec3 b( 0.0, boxDimension, 0.0 );
        Vec3 c( 0.0, 0.0, boxDimension );
        system.setDefaultPeriodicBoxVectors( a, b, c );

        mbpolElectrostaticsForce->setNonbondedMethod( MBPolElectrostaticsForce::PME );
        double ewaldErrorTol = 1e-3;
        mbpolElectrostaticsForce->setEwaldErrorTolerance(ewaldErrorTol);

        mbpolOneBodyForce->setNonbondedMethod(MBPolOneBodyForce::Periodic);

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

    int numberOfWaterMolecules = 256;
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

    std::vector<Vec3> positions;

    positions.push_back(Vec3(  2.4848,  5.1657,  0.0968));
    positions.push_back(Vec3(  1.5828,  4.9797,  0.3590));
    positions.push_back(Vec3(  2.6553,  6.0465,  0.4314));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -7.2150,  5.1657,  0.0968));
    positions.push_back(Vec3( -8.1170,  4.9797,  0.3590));
    positions.push_back(Vec3( -7.0445,  6.0465,  0.4314));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  2.4848, -4.5342,  0.0968));
    positions.push_back(Vec3(  1.5828, -4.7201,  0.3590));
    positions.push_back(Vec3(  2.6553, -3.6534,  0.4314));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -7.2150, -4.5342,  0.0968));
    positions.push_back(Vec3( -8.1170, -4.7201,  0.3590));
    positions.push_back(Vec3( -7.0445, -3.6534,  0.4314));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  2.4848,  5.1657, -9.6031));
    positions.push_back(Vec3(  1.5828,  4.9797, -9.3409));
    positions.push_back(Vec3(  2.6553,  6.0465, -9.2685));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -7.2150,  5.1657, -9.6031));
    positions.push_back(Vec3( -8.1170,  4.9797, -9.3409));
    positions.push_back(Vec3( -7.0445,  6.0465, -9.2685));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  2.4848, -4.5342, -9.6031));
    positions.push_back(Vec3(  1.5828, -4.7201, -9.3409));
    positions.push_back(Vec3(  2.6553, -3.6534, -9.2685));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -7.2150, -4.5342, -9.6031));
    positions.push_back(Vec3( -8.1170, -4.7201, -9.3409));
    positions.push_back(Vec3( -7.0445, -3.6534, -9.2685));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.3953,  3.6373,  5.9278));
    positions.push_back(Vec3(  7.1025,  3.1702,  6.7108));
    positions.push_back(Vec3(  6.8111,  3.3316,  5.2335));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.3045,  3.6373,  5.9278));
    positions.push_back(Vec3( -2.5973,  3.1702,  6.7108));
    positions.push_back(Vec3( -2.8887,  3.3316,  5.2335));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.3953, 13.3371,  5.9278));
    positions.push_back(Vec3(  7.1025, 12.8701,  6.7108));
    positions.push_back(Vec3(  6.8111, 13.0315,  5.2335));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.3045, 13.3371,  5.9278));
    positions.push_back(Vec3( -2.5973, 12.8701,  6.7108));
    positions.push_back(Vec3( -2.8887, 13.0315,  5.2335));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.3953,  3.6373, -3.7720));
    positions.push_back(Vec3(  7.1025,  3.1702, -2.9891));
    positions.push_back(Vec3(  6.8111,  3.3316, -4.4664));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.3045,  3.6373, -3.7720));
    positions.push_back(Vec3( -2.5973,  3.1702, -2.9891));
    positions.push_back(Vec3( -2.8887,  3.3316, -4.4664));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.3953, 13.3371, -3.7720));
    positions.push_back(Vec3(  7.1025, 12.8701, -2.9891));
    positions.push_back(Vec3(  6.8111, 13.0315, -4.4664));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.3045, 13.3371, -3.7720));
    positions.push_back(Vec3( -2.5973, 12.8701, -2.9891));
    positions.push_back(Vec3( -2.8887, 13.0315, -4.4664));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.4351,  5.9345,  7.9831));
    positions.push_back(Vec3(  1.1996,  6.3758,  8.3541));
    positions.push_back(Vec3(  0.2585,  6.4049,  7.1680));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.1349,  5.9345,  7.9831));
    positions.push_back(Vec3( 10.8994,  6.3758,  8.3541));
    positions.push_back(Vec3(  9.9583,  6.4049,  7.1680));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.4351, -3.7653,  7.9831));
    positions.push_back(Vec3(  1.1996, -3.3241,  8.3541));
    positions.push_back(Vec3(  0.2585, -3.2949,  7.1680));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.1349, -3.7653,  7.9831));
    positions.push_back(Vec3( 10.8994, -3.3241,  8.3541));
    positions.push_back(Vec3(  9.9583, -3.2949,  7.1680));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.4351,  5.9345, -1.7168));
    positions.push_back(Vec3(  1.1996,  6.3758, -1.3457));
    positions.push_back(Vec3(  0.2585,  6.4049, -2.5319));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.1349,  5.9345, -1.7168));
    positions.push_back(Vec3( 10.8994,  6.3758, -1.3457));
    positions.push_back(Vec3(  9.9583,  6.4049, -2.5319));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.4351, -3.7653, -1.7168));
    positions.push_back(Vec3(  1.1996, -3.3241, -1.3457));
    positions.push_back(Vec3(  0.2585, -3.2949, -2.5319));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.1349, -3.7653, -1.7168));
    positions.push_back(Vec3( 10.8994, -3.3241, -1.3457));
    positions.push_back(Vec3(  9.9583, -3.2949, -2.5319));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.4290,  9.1817,  0.0397));
    positions.push_back(Vec3(  0.8129,  8.7082,  0.7781));
    positions.push_back(Vec3(  1.0378,  9.0322, -0.6840));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.1289,  9.1817,  0.0397));
    positions.push_back(Vec3( 10.5128,  8.7082,  0.7781));
    positions.push_back(Vec3( 10.7376,  9.0322, -0.6840));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.4290, -0.5181,  0.0397));
    positions.push_back(Vec3(  0.8129, -0.9917,  0.7781));
    positions.push_back(Vec3(  1.0378, -0.6677, -0.6840));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.1289, -0.5181,  0.0397));
    positions.push_back(Vec3( 10.5128, -0.9917,  0.7781));
    positions.push_back(Vec3( 10.7376, -0.6677, -0.6840));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.4290,  9.1817,  9.7396));
    positions.push_back(Vec3(  0.8129,  8.7082, 10.4780));
    positions.push_back(Vec3(  1.0378,  9.0322,  9.0158));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.1289,  9.1817,  9.7396));
    positions.push_back(Vec3( 10.5128,  8.7082, 10.4780));
    positions.push_back(Vec3( 10.7376,  9.0322,  9.0158));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.4290, -0.5181,  9.7396));
    positions.push_back(Vec3(  0.8129, -0.9917, 10.4780));
    positions.push_back(Vec3(  1.0378, -0.6677,  9.0158));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.1289, -0.5181,  9.7396));
    positions.push_back(Vec3( 10.5128, -0.9917, 10.4780));
    positions.push_back(Vec3( 10.7376, -0.6677,  9.0158));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  1.1274,  2.0947,  0.3177));
    positions.push_back(Vec3(  1.1563,  2.5594, -0.5190));
    positions.push_back(Vec3(  0.6007,  1.3148,  0.1410));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.8273,  2.0947,  0.3177));
    positions.push_back(Vec3( 10.8561,  2.5594, -0.5190));
    positions.push_back(Vec3( 10.3006,  1.3148,  0.1410));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  1.1274, 11.7945,  0.3177));
    positions.push_back(Vec3(  1.1563, 12.2592, -0.5190));
    positions.push_back(Vec3(  0.6007, 11.0146,  0.1410));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.8273, 11.7945,  0.3177));
    positions.push_back(Vec3( 10.8561, 12.2592, -0.5190));
    positions.push_back(Vec3( 10.3006, 11.0146,  0.1410));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  1.1274,  2.0947, -9.3822));
    positions.push_back(Vec3(  1.1563,  2.5594,-10.2189));
    positions.push_back(Vec3(  0.6007,  1.3148, -9.5588));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.8273,  2.0947, -9.3822));
    positions.push_back(Vec3( 10.8561,  2.5594,-10.2189));
    positions.push_back(Vec3( 10.3006,  1.3148, -9.5588));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  1.1274, 11.7945, -9.3822));
    positions.push_back(Vec3(  1.1563, 12.2592,-10.2189));
    positions.push_back(Vec3(  0.6007, 11.0146, -9.5588));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.8273, 11.7945, -9.3822));
    positions.push_back(Vec3( 10.8561, 12.2592,-10.2189));
    positions.push_back(Vec3( 10.3006, 11.0146, -9.5588));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  3.9522,  5.0413,  7.2224));
    positions.push_back(Vec3(  4.4465,  4.2249,  7.2998));
    positions.push_back(Vec3(  3.4270,  5.0807,  8.0222));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.7477,  5.0413,  7.2224));
    positions.push_back(Vec3( -5.2534,  4.2249,  7.2998));
    positions.push_back(Vec3( -6.2728,  5.0807,  8.0222));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  3.9522, 14.7411,  7.2224));
    positions.push_back(Vec3(  4.4465, 13.9247,  7.2998));
    positions.push_back(Vec3(  3.4270, 14.7805,  8.0222));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.7477, 14.7411,  7.2224));
    positions.push_back(Vec3( -5.2534, 13.9247,  7.2998));
    positions.push_back(Vec3( -6.2728, 14.7805,  8.0222));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  3.9522,  5.0413, -2.4774));
    positions.push_back(Vec3(  4.4465,  4.2249, -2.4000));
    positions.push_back(Vec3(  3.4270,  5.0807, -1.6777));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.7477,  5.0413, -2.4774));
    positions.push_back(Vec3( -5.2534,  4.2249, -2.4000));
    positions.push_back(Vec3( -6.2728,  5.0807, -1.6777));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  3.9522, 14.7411, -2.4774));
    positions.push_back(Vec3(  4.4465, 13.9247, -2.4000));
    positions.push_back(Vec3(  3.4270, 14.7805, -1.6777));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.7477, 14.7411, -2.4774));
    positions.push_back(Vec3( -5.2534, 13.9247, -2.4000));
    positions.push_back(Vec3( -6.2728, 14.7805, -1.6777));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.6495,  1.1950,  3.3320));
    positions.push_back(Vec3( -0.1767,  0.7328,  3.4753));
    positions.push_back(Vec3(  0.7311,  1.2514,  2.3796));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.3494,  1.1950,  3.3320));
    positions.push_back(Vec3(  9.5232,  0.7328,  3.4753));
    positions.push_back(Vec3( 10.4309,  1.2514,  2.3796));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.6495, 10.8949,  3.3320));
    positions.push_back(Vec3( -0.1767, 10.4326,  3.4753));
    positions.push_back(Vec3(  0.7311, 10.9512,  2.3796));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.3494, 10.8949,  3.3320));
    positions.push_back(Vec3(  9.5232, 10.4326,  3.4753));
    positions.push_back(Vec3( 10.4309, 10.9512,  2.3796));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.6495,  1.1950, -6.3679));
    positions.push_back(Vec3( -0.1767,  0.7328, -6.2246));
    positions.push_back(Vec3(  0.7311,  1.2514, -7.3203));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.3494,  1.1950, -6.3679));
    positions.push_back(Vec3(  9.5232,  0.7328, -6.2246));
    positions.push_back(Vec3( 10.4309,  1.2514, -7.3203));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.6495, 10.8949, -6.3679));
    positions.push_back(Vec3( -0.1767, 10.4326, -6.2246));
    positions.push_back(Vec3(  0.7311, 10.9512, -7.3203));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.3494, 10.8949, -6.3679));
    positions.push_back(Vec3(  9.5232, 10.4326, -6.2246));
    positions.push_back(Vec3( 10.4309, 10.9512, -7.3203));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.3579,  9.6839,  9.0892));
    positions.push_back(Vec3(  7.4911, 10.6131,  9.2784));
    positions.push_back(Vec3(  8.0097,  9.2336,  9.6268));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.3420,  9.6839,  9.0892));
    positions.push_back(Vec3( -2.2088, 10.6131,  9.2784));
    positions.push_back(Vec3( -1.6901,  9.2336,  9.6268));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.3579, -0.0159,  9.0892));
    positions.push_back(Vec3(  7.4911,  0.9132,  9.2784));
    positions.push_back(Vec3(  8.0097, -0.4663,  9.6268));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.3420, -0.0159,  9.0892));
    positions.push_back(Vec3( -2.2088,  0.9132,  9.2784));
    positions.push_back(Vec3( -1.6901, -0.4663,  9.6268));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.3579,  9.6839, -0.6107));
    positions.push_back(Vec3(  7.4911, 10.6131, -0.4214));
    positions.push_back(Vec3(  8.0097,  9.2336, -0.0730));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.3420,  9.6839, -0.6107));
    positions.push_back(Vec3( -2.2088, 10.6131, -0.4214));
    positions.push_back(Vec3( -1.6901,  9.2336, -0.0730));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.3579, -0.0159, -0.6107));
    positions.push_back(Vec3(  7.4911,  0.9132, -0.4214));
    positions.push_back(Vec3(  8.0097, -0.4663, -0.0730));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.3420, -0.0159, -0.6107));
    positions.push_back(Vec3( -2.2088,  0.9132, -0.4214));
    positions.push_back(Vec3( -1.6901, -0.4663, -0.0730));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  2.4531,  7.9731,  8.3301));
    positions.push_back(Vec3(  2.5038,  8.5399,  7.5601));
    positions.push_back(Vec3(  3.3435,  7.9661,  8.6822));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -7.2468,  7.9731,  8.3301));
    positions.push_back(Vec3( -7.1960,  8.5399,  7.5601));
    positions.push_back(Vec3( -6.3564,  7.9661,  8.6822));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  2.4531, -1.7268,  8.3301));
    positions.push_back(Vec3(  2.5038, -1.1599,  7.5601));
    positions.push_back(Vec3(  3.3435, -1.7338,  8.6822));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -7.2468, -1.7268,  8.3301));
    positions.push_back(Vec3( -7.1960, -1.1599,  7.5601));
    positions.push_back(Vec3( -6.3564, -1.7338,  8.6822));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  2.4531,  7.9731, -1.3697));
    positions.push_back(Vec3(  2.5038,  8.5399, -2.1398));
    positions.push_back(Vec3(  3.3435,  7.9661, -1.0177));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -7.2468,  7.9731, -1.3697));
    positions.push_back(Vec3( -7.1960,  8.5399, -2.1398));
    positions.push_back(Vec3( -6.3564,  7.9661, -1.0177));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  2.4531, -1.7268, -1.3697));
    positions.push_back(Vec3(  2.5038, -1.1599, -2.1398));
    positions.push_back(Vec3(  3.3435, -1.7338, -1.0177));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -7.2468, -1.7268, -1.3697));
    positions.push_back(Vec3( -7.1960, -1.1599, -2.1398));
    positions.push_back(Vec3( -6.3564, -1.7338, -1.0177));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  2.2120,  7.0752,  3.6984));
    positions.push_back(Vec3(  3.0194,  7.4579,  3.3542));
    positions.push_back(Vec3(  2.4572,  6.1850,  3.9518));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 11.9118,  7.0752,  3.6984));
    positions.push_back(Vec3( 12.7193,  7.4579,  3.3542));
    positions.push_back(Vec3( 12.1571,  6.1850,  3.9518));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  2.2120, -2.6246,  3.6984));
    positions.push_back(Vec3(  3.0194, -2.2420,  3.3542));
    positions.push_back(Vec3(  2.4572, -3.5149,  3.9518));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 11.9118, -2.6246,  3.6984));
    positions.push_back(Vec3( 12.7193, -2.2420,  3.3542));
    positions.push_back(Vec3( 12.1571, -3.5149,  3.9518));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  2.2120,  7.0752, -6.0014));
    positions.push_back(Vec3(  3.0194,  7.4579, -6.3456));
    positions.push_back(Vec3(  2.4572,  6.1850, -5.7480));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 11.9118,  7.0752, -6.0014));
    positions.push_back(Vec3( 12.7193,  7.4579, -6.3456));
    positions.push_back(Vec3( 12.1571,  6.1850, -5.7480));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  2.2120, -2.6246, -6.0014));
    positions.push_back(Vec3(  3.0194, -2.2420, -6.3456));
    positions.push_back(Vec3(  2.4572, -3.5149, -5.7480));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 11.9118, -2.6246, -6.0014));
    positions.push_back(Vec3( 12.7193, -2.2420, -6.3456));
    positions.push_back(Vec3( 12.1571, -3.5149, -5.7480));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.2917,  3.9954,  2.0164));
    positions.push_back(Vec3( -0.6359,  4.0107,  1.7792));
    positions.push_back(Vec3(  0.6688,  3.3081,  1.4667));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  9.9915,  3.9954,  2.0164));
    positions.push_back(Vec3(  9.0640,  4.0107,  1.7792));
    positions.push_back(Vec3( 10.3686,  3.3081,  1.4667));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.2917, 13.6952,  2.0164));
    positions.push_back(Vec3( -0.6359, 13.7105,  1.7792));
    positions.push_back(Vec3(  0.6688, 13.0079,  1.4667));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  9.9915, 13.6952,  2.0164));
    positions.push_back(Vec3(  9.0640, 13.7105,  1.7792));
    positions.push_back(Vec3( 10.3686, 13.0079,  1.4667));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.2917,  3.9954, -7.6834));
    positions.push_back(Vec3( -0.6359,  4.0107, -7.9206));
    positions.push_back(Vec3(  0.6688,  3.3081, -8.2332));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  9.9915,  3.9954, -7.6834));
    positions.push_back(Vec3(  9.0640,  4.0107, -7.9206));
    positions.push_back(Vec3( 10.3686,  3.3081, -8.2332));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.2917, 13.6952, -7.6834));
    positions.push_back(Vec3( -0.6359, 13.7105, -7.9206));
    positions.push_back(Vec3(  0.6688, 13.0079, -8.2332));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  9.9915, 13.6952, -7.6834));
    positions.push_back(Vec3(  9.0640, 13.7105, -7.9206));
    positions.push_back(Vec3( 10.3686, 13.0079, -8.2332));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  4.6919,  8.8870,  9.5723));
    positions.push_back(Vec3(  4.2352,  9.5845, 10.0432));
    positions.push_back(Vec3(  5.3876,  9.3370,  9.0924));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.0079,  8.8870,  9.5723));
    positions.push_back(Vec3( -5.4646,  9.5845, 10.0432));
    positions.push_back(Vec3( -4.3122,  9.3370,  9.0924));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  4.6919, -0.8128,  9.5723));
    positions.push_back(Vec3(  4.2352, -0.1153, 10.0432));
    positions.push_back(Vec3(  5.3876, -0.3629,  9.0924));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.0079, -0.8128,  9.5723));
    positions.push_back(Vec3( -5.4646, -0.1153, 10.0432));
    positions.push_back(Vec3( -4.3122, -0.3629,  9.0924));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  4.6919,  8.8870, -0.1276));
    positions.push_back(Vec3(  4.2352,  9.5845,  0.3434));
    positions.push_back(Vec3(  5.3876,  9.3370, -0.6075));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.0079,  8.8870, -0.1276));
    positions.push_back(Vec3( -5.4646,  9.5845,  0.3434));
    positions.push_back(Vec3( -4.3122,  9.3370, -0.6075));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  4.6919, -0.8128, -0.1276));
    positions.push_back(Vec3(  4.2352, -0.1153,  0.3434));
    positions.push_back(Vec3(  5.3876, -0.3629, -0.6075));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.0079, -0.8128, -0.1276));
    positions.push_back(Vec3( -5.4646, -0.1153,  0.3434));
    positions.push_back(Vec3( -4.3122, -0.3629, -0.6075));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  3.2053,  4.4465,  4.4811));
    positions.push_back(Vec3(  3.6305,  4.4697,  5.3387));
    positions.push_back(Vec3(  2.4987,  3.8081,  4.5805));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -6.4945,  4.4465,  4.4811));
    positions.push_back(Vec3( -6.0694,  4.4697,  5.3387));
    positions.push_back(Vec3( -7.2012,  3.8081,  4.5805));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  3.2053, 14.1464,  4.4811));
    positions.push_back(Vec3(  3.6305, 14.1695,  5.3387));
    positions.push_back(Vec3(  2.4987, 13.5080,  4.5805));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -6.4945, 14.1464,  4.4811));
    positions.push_back(Vec3( -6.0694, 14.1695,  5.3387));
    positions.push_back(Vec3( -7.2012, 13.5080,  4.5805));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  3.2053,  4.4465, -5.2188));
    positions.push_back(Vec3(  3.6305,  4.4697, -4.3611));
    positions.push_back(Vec3(  2.4987,  3.8081, -5.1193));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -6.4945,  4.4465, -5.2188));
    positions.push_back(Vec3( -6.0694,  4.4697, -4.3611));
    positions.push_back(Vec3( -7.2012,  3.8081, -5.1193));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  3.2053, 14.1464, -5.2188));
    positions.push_back(Vec3(  3.6305, 14.1695, -4.3611));
    positions.push_back(Vec3(  2.4987, 13.5080, -5.1193));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -6.4945, 14.1464, -5.2188));
    positions.push_back(Vec3( -6.0694, 14.1695, -4.3611));
    positions.push_back(Vec3( -7.2012, 13.5080, -5.1193));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.9564,  0.0934,  3.0756));
    positions.push_back(Vec3(  7.6836, -0.7516,  2.7172));
    positions.push_back(Vec3(  7.6458,  0.0817,  3.9813));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -1.7434,  0.0934,  3.0756));
    positions.push_back(Vec3( -2.0163, -0.7516,  2.7172));
    positions.push_back(Vec3( -2.0540,  0.0817,  3.9813));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.9564,  9.7932,  3.0756));
    positions.push_back(Vec3(  7.6836,  8.9483,  2.7172));
    positions.push_back(Vec3(  7.6458,  9.7815,  3.9813));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -1.7434,  9.7932,  3.0756));
    positions.push_back(Vec3( -2.0163,  8.9483,  2.7172));
    positions.push_back(Vec3( -2.0540,  9.7815,  3.9813));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.9564,  0.0934, -6.6242));
    positions.push_back(Vec3(  7.6836, -0.7516, -6.9826));
    positions.push_back(Vec3(  7.6458,  0.0817, -5.7186));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -1.7434,  0.0934, -6.6242));
    positions.push_back(Vec3( -2.0163, -0.7516, -6.9826));
    positions.push_back(Vec3( -2.0540,  0.0817, -5.7186));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.9564,  9.7932, -6.6242));
    positions.push_back(Vec3(  7.6836,  8.9483, -6.9826));
    positions.push_back(Vec3(  7.6458,  9.7815, -5.7186));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -1.7434,  9.7932, -6.6242));
    positions.push_back(Vec3( -2.0163,  8.9483, -6.9826));
    positions.push_back(Vec3( -2.0540,  9.7815, -5.7186));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  3.9655,  1.5890,  0.8078));
    positions.push_back(Vec3(  3.0229,  1.7528,  0.7698));
    positions.push_back(Vec3(  4.3082,  2.2867,  1.3668));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.7343,  1.5890,  0.8078));
    positions.push_back(Vec3( -6.6769,  1.7528,  0.7698));
    positions.push_back(Vec3( -5.3917,  2.2867,  1.3668));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  3.9655, 11.2888,  0.8078));
    positions.push_back(Vec3(  3.0229, 11.4526,  0.7698));
    positions.push_back(Vec3(  4.3082, 11.9866,  1.3668));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.7343, 11.2888,  0.8078));
    positions.push_back(Vec3( -6.6769, 11.4526,  0.7698));
    positions.push_back(Vec3( -5.3917, 11.9866,  1.3668));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  3.9655,  1.5890, -8.8921));
    positions.push_back(Vec3(  3.0229,  1.7528, -8.9300));
    positions.push_back(Vec3(  4.3082,  2.2867, -8.3330));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.7343,  1.5890, -8.8921));
    positions.push_back(Vec3( -6.6769,  1.7528, -8.9300));
    positions.push_back(Vec3( -5.3917,  2.2867, -8.3330));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  3.9655, 11.2888, -8.8921));
    positions.push_back(Vec3(  3.0229, 11.4526, -8.9300));
    positions.push_back(Vec3(  4.3082, 11.9866, -8.3330));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.7343, 11.2888, -8.8921));
    positions.push_back(Vec3( -6.6769, 11.4526, -8.9300));
    positions.push_back(Vec3( -5.3917, 11.9866, -8.3330));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  8.9358,  7.5054,  5.8835));
    positions.push_back(Vec3(  8.3038,  8.1984,  6.0764));
    positions.push_back(Vec3(  9.7390,  7.9714,  5.6495));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -0.7640,  7.5054,  5.8835));
    positions.push_back(Vec3( -1.3960,  8.1984,  6.0764));
    positions.push_back(Vec3(  0.0391,  7.9714,  5.6495));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  8.9358, -2.1944,  5.8835));
    positions.push_back(Vec3(  8.3038, -1.5015,  6.0764));
    positions.push_back(Vec3(  9.7390, -1.7285,  5.6495));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -0.7640, -2.1944,  5.8835));
    positions.push_back(Vec3( -1.3960, -1.5015,  6.0764));
    positions.push_back(Vec3(  0.0391, -1.7285,  5.6495));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  8.9358,  7.5054, -3.8163));
    positions.push_back(Vec3(  8.3038,  8.1984, -3.6234));
    positions.push_back(Vec3(  9.7390,  7.9714, -4.0503));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -0.7640,  7.5054, -3.8163));
    positions.push_back(Vec3( -1.3960,  8.1984, -3.6234));
    positions.push_back(Vec3(  0.0391,  7.9714, -4.0503));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  8.9358, -2.1944, -3.8163));
    positions.push_back(Vec3(  8.3038, -1.5015, -3.6234));
    positions.push_back(Vec3(  9.7390, -1.7285, -4.0503));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -0.7640, -2.1944, -3.8163));
    positions.push_back(Vec3( -1.3960, -1.5015, -3.6234));
    positions.push_back(Vec3(  0.0391, -1.7285, -4.0503));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.0834,  6.2049,  9.0922));
    positions.push_back(Vec3(  6.9459,  7.1017,  9.3983));
    positions.push_back(Vec3(  7.6433,  6.2976,  8.3211));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.6165,  6.2049,  9.0922));
    positions.push_back(Vec3( -2.7540,  7.1017,  9.3983));
    positions.push_back(Vec3( -2.0565,  6.2976,  8.3211));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.0834, -3.4950,  9.0922));
    positions.push_back(Vec3(  6.9459, -2.5982,  9.3983));
    positions.push_back(Vec3(  7.6433, -3.4023,  8.3211));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.6165, -3.4950,  9.0922));
    positions.push_back(Vec3( -2.7540, -2.5982,  9.3983));
    positions.push_back(Vec3( -2.0565, -3.4023,  8.3211));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.0834,  6.2049, -0.6076));
    positions.push_back(Vec3(  6.9459,  7.1017, -0.3015));
    positions.push_back(Vec3(  7.6433,  6.2976, -1.3788));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.6165,  6.2049, -0.6076));
    positions.push_back(Vec3( -2.7540,  7.1017, -0.3015));
    positions.push_back(Vec3( -2.0565,  6.2976, -1.3788));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.0834, -3.4950, -0.6076));
    positions.push_back(Vec3(  6.9459, -2.5982, -0.3015));
    positions.push_back(Vec3(  7.6433, -3.4023, -1.3788));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.6165, -3.4950, -0.6076));
    positions.push_back(Vec3( -2.7540, -2.5982, -0.3015));
    positions.push_back(Vec3( -2.0565, -3.4023, -1.3788));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.6777,  2.6707,  8.7697));
    positions.push_back(Vec3(  7.6784,  3.1349,  9.6072));
    positions.push_back(Vec3(  8.5402,  2.8504,  8.3948));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.0222,  2.6707,  8.7697));
    positions.push_back(Vec3( -2.0214,  3.1349,  9.6072));
    positions.push_back(Vec3( -1.1597,  2.8504,  8.3948));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.6777, 12.3705,  8.7697));
    positions.push_back(Vec3(  7.6784, 12.8347,  9.6072));
    positions.push_back(Vec3(  8.5402, 12.5503,  8.3948));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.0222, 12.3705,  8.7697));
    positions.push_back(Vec3( -2.0214, 12.8347,  9.6072));
    positions.push_back(Vec3( -1.1597, 12.5503,  8.3948));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.6777,  2.6707, -0.9301));
    positions.push_back(Vec3(  7.6784,  3.1349, -0.0926));
    positions.push_back(Vec3(  8.5402,  2.8504, -1.3050));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.0222,  2.6707, -0.9301));
    positions.push_back(Vec3( -2.0214,  3.1349, -0.0926));
    positions.push_back(Vec3( -1.1597,  2.8504, -1.3050));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.6777, 12.3705, -0.9301));
    positions.push_back(Vec3(  7.6784, 12.8347, -0.0926));
    positions.push_back(Vec3(  8.5402, 12.5503, -1.3050));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.0222, 12.3705, -0.9301));
    positions.push_back(Vec3( -2.0214, 12.8347, -0.0926));
    positions.push_back(Vec3( -1.1597, 12.5503, -1.3050));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  5.0203,  2.8164,  8.1340));
    positions.push_back(Vec3(  5.8974,  2.5309,  8.3908));
    positions.push_back(Vec3(  4.4888,  2.6972,  8.9215));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -4.6796,  2.8164,  8.1340));
    positions.push_back(Vec3( -3.8024,  2.5309,  8.3908));
    positions.push_back(Vec3( -5.2110,  2.6972,  8.9215));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  5.0203, 12.5162,  8.1340));
    positions.push_back(Vec3(  5.8974, 12.2308,  8.3908));
    positions.push_back(Vec3(  4.4888, 12.3970,  8.9215));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -4.6796, 12.5162,  8.1340));
    positions.push_back(Vec3( -3.8024, 12.2308,  8.3908));
    positions.push_back(Vec3( -5.2110, 12.3970,  8.9215));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  5.0203,  2.8164, -1.5659));
    positions.push_back(Vec3(  5.8974,  2.5309, -1.3091));
    positions.push_back(Vec3(  4.4888,  2.6972, -0.7783));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -4.6796,  2.8164, -1.5659));
    positions.push_back(Vec3( -3.8024,  2.5309, -1.3091));
    positions.push_back(Vec3( -5.2110,  2.6972, -0.7783));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  5.0203, 12.5162, -1.5659));
    positions.push_back(Vec3(  5.8974, 12.2308, -1.3091));
    positions.push_back(Vec3(  4.4888, 12.3970, -0.7783));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -4.6796, 12.5162, -1.5659));
    positions.push_back(Vec3( -3.8024, 12.2308, -1.3091));
    positions.push_back(Vec3( -5.2110, 12.3970, -0.7783));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  5.2823,  3.1509,  3.0410));
    positions.push_back(Vec3(  4.5593,  3.6794,  3.3797));
    positions.push_back(Vec3(  5.2594,  2.3493,  3.5644));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -4.4175,  3.1509,  3.0410));
    positions.push_back(Vec3( -5.1406,  3.6794,  3.3797));
    positions.push_back(Vec3( -4.4405,  2.3493,  3.5644));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  5.2823, 12.8507,  3.0410));
    positions.push_back(Vec3(  4.5593, 13.3793,  3.3797));
    positions.push_back(Vec3(  5.2594, 12.0492,  3.5644));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -4.4175, 12.8507,  3.0410));
    positions.push_back(Vec3( -5.1406, 13.3793,  3.3797));
    positions.push_back(Vec3( -4.4405, 12.0492,  3.5644));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  5.2823,  3.1509, -6.6588));
    positions.push_back(Vec3(  4.5593,  3.6794, -6.3202));
    positions.push_back(Vec3(  5.2594,  2.3493, -6.1355));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -4.4175,  3.1509, -6.6588));
    positions.push_back(Vec3( -5.1406,  3.6794, -6.3202));
    positions.push_back(Vec3( -4.4405,  2.3493, -6.1355));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  5.2823, 12.8507, -6.6588));
    positions.push_back(Vec3(  4.5593, 13.3793, -6.3202));
    positions.push_back(Vec3(  5.2594, 12.0492, -6.1355));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -4.4175, 12.8507, -6.6588));
    positions.push_back(Vec3( -5.1406, 13.3793, -6.3202));
    positions.push_back(Vec3( -4.4405, 12.0492, -6.1355));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  1.6714,  9.5315,  5.6251));
    positions.push_back(Vec3(  2.4847,  9.0569,  5.4516));
    positions.push_back(Vec3(  1.4063,  9.8739,  4.7711));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 11.3712,  9.5315,  5.6251));
    positions.push_back(Vec3( 12.1846,  9.0569,  5.4516));
    positions.push_back(Vec3( 11.1062,  9.8739,  4.7711));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  1.6714, -0.1684,  5.6251));
    positions.push_back(Vec3(  2.4847, -0.6429,  5.4516));
    positions.push_back(Vec3(  1.4063,  0.1740,  4.7711));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 11.3712, -0.1684,  5.6251));
    positions.push_back(Vec3( 12.1846, -0.6429,  5.4516));
    positions.push_back(Vec3( 11.1062,  0.1740,  4.7711));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  1.6714,  9.5315, -4.0747));
    positions.push_back(Vec3(  2.4847,  9.0569, -4.2482));
    positions.push_back(Vec3(  1.4063,  9.8739, -4.9288));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 11.3712,  9.5315, -4.0747));
    positions.push_back(Vec3( 12.1846,  9.0569, -4.2482));
    positions.push_back(Vec3( 11.1062,  9.8739, -4.9288));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  1.6714, -0.1684, -4.0747));
    positions.push_back(Vec3(  2.4847, -0.6429, -4.2482));
    positions.push_back(Vec3(  1.4063,  0.1740, -4.9288));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 11.3712, -0.1684, -4.0747));
    positions.push_back(Vec3( 12.1846, -0.6429, -4.2482));
    positions.push_back(Vec3( 11.1062,  0.1740, -4.9288));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.8786,  3.3263,  4.9125));
    positions.push_back(Vec3(  0.2092,  4.0105,  4.8893));
    positions.push_back(Vec3(  0.6296,  2.7230,  4.2119));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.5785,  3.3263,  4.9125));
    positions.push_back(Vec3(  9.9090,  4.0105,  4.8893));
    positions.push_back(Vec3( 10.3295,  2.7230,  4.2119));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.8786, 13.0261,  4.9125));
    positions.push_back(Vec3(  0.2092, 13.7103,  4.8893));
    positions.push_back(Vec3(  0.6296, 12.4228,  4.2119));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.5785, 13.0261,  4.9125));
    positions.push_back(Vec3(  9.9090, 13.7103,  4.8893));
    positions.push_back(Vec3( 10.3295, 12.4228,  4.2119));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.8786,  3.3263, -4.7874));
    positions.push_back(Vec3(  0.2092,  4.0105, -4.8105));
    positions.push_back(Vec3(  0.6296,  2.7230, -5.4880));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.5785,  3.3263, -4.7874));
    positions.push_back(Vec3(  9.9090,  4.0105, -4.8105));
    positions.push_back(Vec3( 10.3295,  2.7230, -5.4880));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.8786, 13.0261, -4.7874));
    positions.push_back(Vec3(  0.2092, 13.7103, -4.8105));
    positions.push_back(Vec3(  0.6296, 12.4228, -5.4880));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.5785, 13.0261, -4.7874));
    positions.push_back(Vec3(  9.9090, 13.7103, -4.8105));
    positions.push_back(Vec3( 10.3295, 12.4228, -5.4880));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  5.7789,  6.9883,  6.4703));
    positions.push_back(Vec3(  5.2997,  6.1599,  6.4998));
    positions.push_back(Vec3(  6.6819,  6.7365,  6.2752));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -3.9209,  6.9883,  6.4703));
    positions.push_back(Vec3( -4.4002,  6.1599,  6.4998));
    positions.push_back(Vec3( -3.0179,  6.7365,  6.2752));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  5.7789, -2.7115,  6.4703));
    positions.push_back(Vec3(  5.2997, -3.5400,  6.4998));
    positions.push_back(Vec3(  6.6819, -2.9634,  6.2752));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -3.9209, -2.7115,  6.4703));
    positions.push_back(Vec3( -4.4002, -3.5400,  6.4998));
    positions.push_back(Vec3( -3.0179, -2.9634,  6.2752));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  5.7789,  6.9883, -3.2295));
    positions.push_back(Vec3(  5.2997,  6.1599, -3.2001));
    positions.push_back(Vec3(  6.6819,  6.7365, -3.4246));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -3.9209,  6.9883, -3.2295));
    positions.push_back(Vec3( -4.4002,  6.1599, -3.2001));
    positions.push_back(Vec3( -3.0179,  6.7365, -3.4246));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  5.7789, -2.7115, -3.2295));
    positions.push_back(Vec3(  5.2997, -3.5400, -3.2001));
    positions.push_back(Vec3(  6.6819, -2.9634, -3.4246));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -3.9209, -2.7115, -3.2295));
    positions.push_back(Vec3( -4.4002, -3.5400, -3.2001));
    positions.push_back(Vec3( -3.0179, -2.9634, -3.4246));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  5.1477,  1.8605,  5.3043));
    positions.push_back(Vec3(  4.5913,  2.2471,  5.9808));
    positions.push_back(Vec3(  4.7901,  0.9821,  5.1724));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -4.5521,  1.8605,  5.3043));
    positions.push_back(Vec3( -5.1086,  2.2471,  5.9808));
    positions.push_back(Vec3( -4.9097,  0.9821,  5.1724));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  5.1477, 11.5604,  5.3043));
    positions.push_back(Vec3(  4.5913, 11.9470,  5.9808));
    positions.push_back(Vec3(  4.7901, 10.6820,  5.1724));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -4.5521, 11.5604,  5.3043));
    positions.push_back(Vec3( -5.1086, 11.9470,  5.9808));
    positions.push_back(Vec3( -4.9097, 10.6820,  5.1724));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  5.1477,  1.8605, -4.3956));
    positions.push_back(Vec3(  4.5913,  2.2471, -3.7190));
    positions.push_back(Vec3(  4.7901,  0.9821, -4.5275));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -4.5521,  1.8605, -4.3956));
    positions.push_back(Vec3( -5.1086,  2.2471, -3.7190));
    positions.push_back(Vec3( -4.9097,  0.9821, -4.5275));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  5.1477, 11.5604, -4.3956));
    positions.push_back(Vec3(  4.5913, 11.9470, -3.7190));
    positions.push_back(Vec3(  4.7901, 10.6820, -4.5275));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -4.5521, 11.5604, -4.3956));
    positions.push_back(Vec3( -5.1086, 11.9470, -3.7190));
    positions.push_back(Vec3( -4.9097, 10.6820, -4.5275));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.5758,  6.7918,  1.5804));
    positions.push_back(Vec3(  0.5892,  5.8353,  1.6214));
    positions.push_back(Vec3(  0.9161,  7.0733,  2.4300));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.2756,  6.7918,  1.5804));
    positions.push_back(Vec3( 10.2890,  5.8353,  1.6214));
    positions.push_back(Vec3( 10.6159,  7.0733,  2.4300));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.5758, -2.9080,  1.5804));
    positions.push_back(Vec3(  0.5892, -3.8646,  1.6214));
    positions.push_back(Vec3(  0.9161, -2.6266,  2.4300));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.2756, -2.9080,  1.5804));
    positions.push_back(Vec3( 10.2890, -3.8646,  1.6214));
    positions.push_back(Vec3( 10.6159, -2.6266,  2.4300));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.5758,  6.7918, -8.1194));
    positions.push_back(Vec3(  0.5892,  5.8353, -8.0785));
    positions.push_back(Vec3(  0.9161,  7.0733, -7.2698));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.2756,  6.7918, -8.1194));
    positions.push_back(Vec3( 10.2890,  5.8353, -8.0785));
    positions.push_back(Vec3( 10.6159,  7.0733, -7.2698));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.5758, -2.9080, -8.1194));
    positions.push_back(Vec3(  0.5892, -3.8646, -8.0785));
    positions.push_back(Vec3(  0.9161, -2.6266, -7.2698));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.2756, -2.9080, -8.1194));
    positions.push_back(Vec3( 10.2890, -3.8646, -8.0785));
    positions.push_back(Vec3( 10.6159, -2.6266, -7.2698));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.3787,  5.4440,  3.7297));
    positions.push_back(Vec3(  7.8518,  5.2925,  4.5483));
    positions.push_back(Vec3(  6.8476,  4.6563,  3.6100));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.3212,  5.4440,  3.7297));
    positions.push_back(Vec3( -1.8481,  5.2925,  4.5483));
    positions.push_back(Vec3( -2.8522,  4.6563,  3.6100));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.3787, -4.2558,  3.7297));
    positions.push_back(Vec3(  7.8518, -4.4074,  4.5483));
    positions.push_back(Vec3(  6.8476, -5.0436,  3.6100));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.3212, -4.2558,  3.7297));
    positions.push_back(Vec3( -1.8481, -4.4074,  4.5483));
    positions.push_back(Vec3( -2.8522, -5.0436,  3.6100));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.3787,  5.4440, -5.9701));
    positions.push_back(Vec3(  7.8518,  5.2925, -5.1516));
    positions.push_back(Vec3(  6.8476,  4.6563, -6.0898));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.3212,  5.4440, -5.9701));
    positions.push_back(Vec3( -1.8481,  5.2925, -5.1516));
    positions.push_back(Vec3( -2.8522,  4.6563, -6.0898));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.3787, -4.2558, -5.9701));
    positions.push_back(Vec3(  7.8518, -4.4074, -5.1516));
    positions.push_back(Vec3(  6.8476, -5.0436, -6.0898));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.3212, -4.2558, -5.9701));
    positions.push_back(Vec3( -1.8481, -4.4074, -5.1516));
    positions.push_back(Vec3( -2.8522, -5.0436, -6.0898));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.9779,  0.3681,  5.9563));
    positions.push_back(Vec3(  8.4547,  1.1694,  6.1740));
    positions.push_back(Vec3(  7.4025,  0.2196,  6.7071));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -1.7219,  0.3681,  5.9563));
    positions.push_back(Vec3( -1.2451,  1.1694,  6.1740));
    positions.push_back(Vec3( -2.2974,  0.2196,  6.7071));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.9779, 10.0679,  5.9563));
    positions.push_back(Vec3(  8.4547, 10.8692,  6.1740));
    positions.push_back(Vec3(  7.4025,  9.9195,  6.7071));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -1.7219, 10.0679,  5.9563));
    positions.push_back(Vec3( -1.2451, 10.8692,  6.1740));
    positions.push_back(Vec3( -2.2974,  9.9195,  6.7071));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.9779,  0.3681, -3.7435));
    positions.push_back(Vec3(  8.4547,  1.1694, -3.5258));
    positions.push_back(Vec3(  7.4025,  0.2196, -2.9927));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -1.7219,  0.3681, -3.7435));
    positions.push_back(Vec3( -1.2451,  1.1694, -3.5258));
    positions.push_back(Vec3( -2.2974,  0.2196, -2.9927));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.9779, 10.0679, -3.7435));
    positions.push_back(Vec3(  8.4547, 10.8692, -3.5258));
    positions.push_back(Vec3(  7.4025,  9.9195, -2.9927));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -1.7219, 10.0679, -3.7435));
    positions.push_back(Vec3( -1.2451, 10.8692, -3.5258));
    positions.push_back(Vec3( -2.2974,  9.9195, -2.9927));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.5975,  3.0690,  7.5894));
    positions.push_back(Vec3(  0.5145,  4.0047,  7.4034));
    positions.push_back(Vec3(  1.2649,  2.7606,  6.9759));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.2974,  3.0690,  7.5894));
    positions.push_back(Vec3( 10.2144,  4.0047,  7.4034));
    positions.push_back(Vec3( 10.9647,  2.7606,  6.9759));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.5975, 12.7689,  7.5894));
    positions.push_back(Vec3(  0.5145, 13.7045,  7.4034));
    positions.push_back(Vec3(  1.2649, 12.4605,  6.9759));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.2974, 12.7689,  7.5894));
    positions.push_back(Vec3( 10.2144, 13.7045,  7.4034));
    positions.push_back(Vec3( 10.9647, 12.4605,  6.9759));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.5975,  3.0690, -2.1105));
    positions.push_back(Vec3(  0.5145,  4.0047, -2.2964));
    positions.push_back(Vec3(  1.2649,  2.7606, -2.7240));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.2974,  3.0690, -2.1105));
    positions.push_back(Vec3( 10.2144,  4.0047, -2.2964));
    positions.push_back(Vec3( 10.9647,  2.7606, -2.7240));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  0.5975, 12.7689, -2.1105));
    positions.push_back(Vec3(  0.5145, 13.7045, -2.2964));
    positions.push_back(Vec3(  1.2649, 12.4605, -2.7240));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( 10.2974, 12.7689, -2.1105));
    positions.push_back(Vec3( 10.2144, 13.7045, -2.2964));
    positions.push_back(Vec3( 10.9647, 12.4605, -2.7240));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.3346,  4.0084,  1.2247));
    positions.push_back(Vec3(  7.1550,  4.9423,  1.1128));
    positions.push_back(Vec3(  6.5083,  3.6410,  1.5394));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.3652,  4.0084,  1.2247));
    positions.push_back(Vec3( -2.5449,  4.9423,  1.1128));
    positions.push_back(Vec3( -3.1916,  3.6410,  1.5394));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.3346, 13.7083,  1.2247));
    positions.push_back(Vec3(  7.1550, 14.6421,  1.1128));
    positions.push_back(Vec3(  6.5083, 13.3409,  1.5394));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.3652, 13.7083,  1.2247));
    positions.push_back(Vec3( -2.5449, 14.6421,  1.1128));
    positions.push_back(Vec3( -3.1916, 13.3409,  1.5394));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.3346,  4.0084, -8.4751));
    positions.push_back(Vec3(  7.1550,  4.9423, -8.5870));
    positions.push_back(Vec3(  6.5083,  3.6410, -8.1604));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.3652,  4.0084, -8.4751));
    positions.push_back(Vec3( -2.5449,  4.9423, -8.5870));
    positions.push_back(Vec3( -3.1916,  3.6410, -8.1604));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.3346, 13.7083, -8.4751));
    positions.push_back(Vec3(  7.1550, 14.6421, -8.5870));
    positions.push_back(Vec3(  6.5083, 13.3409, -8.1604));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -2.3652, 13.7083, -8.4751));
    positions.push_back(Vec3( -2.5449, 14.6421, -8.5870));
    positions.push_back(Vec3( -3.1916, 13.3409, -8.1604));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  4.4060,  8.9037,  5.0120));
    positions.push_back(Vec3(  4.8385,  8.1184,  5.3483));
    positions.push_back(Vec3(  4.4503,  8.8130,  4.0598));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.2939,  8.9037,  5.0120));
    positions.push_back(Vec3( -4.8614,  8.1184,  5.3483));
    positions.push_back(Vec3( -5.2495,  8.8130,  4.0598));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  4.4060, -0.7962,  5.0120));
    positions.push_back(Vec3(  4.8385, -1.5814,  5.3483));
    positions.push_back(Vec3(  4.4503, -0.8868,  4.0598));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.2939, -0.7962,  5.0120));
    positions.push_back(Vec3( -4.8614, -1.5814,  5.3483));
    positions.push_back(Vec3( -5.2495, -0.8868,  4.0598));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  4.4060,  8.9037, -4.6879));
    positions.push_back(Vec3(  4.8385,  8.1184, -4.3515));
    positions.push_back(Vec3(  4.4503,  8.8130, -5.6401));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.2939,  8.9037, -4.6879));
    positions.push_back(Vec3( -4.8614,  8.1184, -4.3515));
    positions.push_back(Vec3( -5.2495,  8.8130, -5.6401));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  4.4060, -0.7962, -4.6879));
    positions.push_back(Vec3(  4.8385, -1.5814, -4.3515));
    positions.push_back(Vec3(  4.4503, -0.8868, -5.6401));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.2939, -0.7962, -4.6879));
    positions.push_back(Vec3( -4.8614, -1.5814, -4.3515));
    positions.push_back(Vec3( -5.2495, -0.8868, -5.6401));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.7425,  7.4227,  2.0653));
    positions.push_back(Vec3(  7.5631,  6.7355,  2.7075));
    positions.push_back(Vec3(  8.6027,  7.1995,  1.7088));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -1.9574,  7.4227,  2.0653));
    positions.push_back(Vec3( -2.1368,  6.7355,  2.7075));
    positions.push_back(Vec3( -1.0972,  7.1995,  1.7088));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.7425, -2.2771,  2.0653));
    positions.push_back(Vec3(  7.5631, -2.9644,  2.7075));
    positions.push_back(Vec3(  8.6027, -2.5003,  1.7088));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -1.9574, -2.2771,  2.0653));
    positions.push_back(Vec3( -2.1368, -2.9644,  2.7075));
    positions.push_back(Vec3( -1.0972, -2.5003,  1.7088));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.7425,  7.4227, -7.6346));
    positions.push_back(Vec3(  7.5631,  6.7355, -6.9924));
    positions.push_back(Vec3(  8.6027,  7.1995, -7.9911));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -1.9574,  7.4227, -7.6346));
    positions.push_back(Vec3( -2.1368,  6.7355, -6.9924));
    positions.push_back(Vec3( -1.0972,  7.1995, -7.9911));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  7.7425, -2.2771, -7.6346));
    positions.push_back(Vec3(  7.5631, -2.9644, -6.9924));
    positions.push_back(Vec3(  8.6027, -2.5003, -7.9911));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -1.9574, -2.2771, -7.6346));
    positions.push_back(Vec3( -2.1368, -2.9644, -6.9924));
    positions.push_back(Vec3( -1.0972, -2.5003, -7.9911));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  4.6574,  7.8884,  2.5117));
    positions.push_back(Vec3(  4.5382,  8.3654,  1.6901));
    positions.push_back(Vec3(  5.0262,  7.0449,  2.2484));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.0425,  7.8884,  2.5117));
    positions.push_back(Vec3( -5.1616,  8.3654,  1.6901));
    positions.push_back(Vec3( -4.6737,  7.0449,  2.2484));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  4.6574, -1.8114,  2.5117));
    positions.push_back(Vec3(  4.5382, -1.3345,  1.6901));
    positions.push_back(Vec3(  5.0262, -2.6549,  2.2484));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.0425, -1.8114,  2.5117));
    positions.push_back(Vec3( -5.1616, -1.3345,  1.6901));
    positions.push_back(Vec3( -4.6737, -2.6549,  2.2484));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  4.6574,  7.8884, -7.1881));
    positions.push_back(Vec3(  4.5382,  8.3654, -8.0098));
    positions.push_back(Vec3(  5.0262,  7.0449, -7.4515));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.0425,  7.8884, -7.1881));
    positions.push_back(Vec3( -5.1616,  8.3654, -8.0098));
    positions.push_back(Vec3( -4.6737,  7.0449, -7.4515));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3(  4.6574, -1.8114, -7.1881));
    positions.push_back(Vec3(  4.5382, -1.3345, -8.0098));
    positions.push_back(Vec3(  5.0262, -2.6549, -7.4515));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));
    positions.push_back(Vec3( -5.0425, -1.8114, -7.1881));
    positions.push_back(Vec3( -5.1616, -1.3345, -8.0098));
    positions.push_back(Vec3( -4.6737, -2.6549, -7.4515));
    positions.push_back(Vec3(  0.0000,  0.0000,  0.0000));


    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            positions[i][j] *= 1e-1;
        }
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

    double tolerance = 1.0e-03;


    double energy = state.getPotentialEnergy() / CalToJoule;

    double expectedEnergy = -2270.88890;

    std::cout << "Total Energy: " << energy << " Kcal/mol "<< std::endl;
    std::cout << "Expected energy: " << expectedEnergy << " Kcal/mol "<< std::endl;

    for ( unsigned int ii = 0; ii < 5; ii++ ){
        state                      = context.getState(State::Energy, false, pow(2, ii));
        std::cout << "Energy " << forceLabels[ii] << ": " << state.getPotentialEnergy() / CalToJoule << " Kcal/mol "<< std::endl;
    }
   std::cout << "Test Successful: " << testName << std::endl << std::endl;

}

int main( int numberOfArguments, char* argv[] ) {

    try {
        std::cout << "TestReferenceMBPolIntegrationTest running test..." << std::endl;

        double boxDimension = 19.3996888399961804/10.;
        runTest( boxDimension );

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
