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
const double cal2joule = 4.184;

void testWater14( ) {

    std::string testName      = "testMBPolWater14Test";


    // Electrostatics

    double virtualSiteWeightO = 0.573293118;
    double virtualSiteWeightH = 0.213353441;

    MBPolElectrostaticsForce* mbpolElectrostaticsForce        = new MBPolElectrostaticsForce();;
    mbpolElectrostaticsForce->setNonbondedMethod( MBPolElectrostaticsForce::NoCutoff );

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


    // Two body interaction
    MBPolTwoBodyForce* mbpolTwoBodyForce = new MBPolTwoBodyForce();
    double cutoff = 1e10;
    mbpolTwoBodyForce->setCutoff( cutoff );
    mbpolTwoBodyForce->setNonbondedMethod(MBPolTwoBodyForce::CutoffNonPeriodic);

    // Three body interaction
    MBPolThreeBodyForce* mbpolThreeBodyForce = new MBPolThreeBodyForce();
    mbpolThreeBodyForce->setCutoff( cutoff );
    mbpolThreeBodyForce->setNonbondedMethod(MBPolThreeBodyForce::CutoffNonPeriodic);

    // Dispersion Force
    MBPolDispersionForce* dispersionForce = new MBPolDispersionForce();
    dispersionForce->setCutoff( cutoff );
    dispersionForce->setNonbondedMethod(MBPolDispersionForce::CutoffNonPeriodic);

    // Setup system

    int numberOfWaterMolecules = 14;
    unsigned int particlesPerMolecule = 4;
    int numberOfParticles          = numberOfWaterMolecules * particlesPerMolecule;

    std::vector<int> particleIndices(particlesPerMolecule);
    System system;
    int waterMoleculeIndex=0;
    for( unsigned int jj = 0; jj < numberOfParticles; jj += particlesPerMolecule ){
        system.addParticle( 1.5999000e+01 ); // Mass
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 0. ); // Virtual Site
        system.setVirtualSite(jj+3, new ThreeParticleAverageSite(jj, jj+1, jj+2,
                                                                   virtualSiteWeightO, virtualSiteWeightH,virtualSiteWeightH));

        particleIndices[0] = jj;
        particleIndices[1] = jj+1;
        particleIndices[2] = jj+2;

        // Charge, dipoles and quadrupoles (zero in MBPol)
        mbpolElectrostaticsForce->addElectrostatics( -5.1966000e-01,
                                            waterMoleculeIndex, 0, 0.001310, 0.001310 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01,
                                            waterMoleculeIndex, 1, 0.000294, 0.000294 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01,
                                            waterMoleculeIndex, 1, 0.000294, 0.000294 );
        mbpolElectrostaticsForce->addElectrostatics(  0.,
                                            waterMoleculeIndex, 2, 0.001310,  0.);
        waterMoleculeIndex++;
        mbpolOneBodyForce->addOneBody( particleIndices );
        mbpolTwoBodyForce->addParticle( particleIndices );
        mbpolThreeBodyForce->addParticle( particleIndices );
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


    // Atom positions O H H, [A], no virtual sites

    std::vector<Vec3> positions;
    positions.push_back(Vec3( -2.349377641189e-01,  1.798934467398e-01,  1.896881820756e-01 ));
    positions.push_back(Vec3(  1.788456192811e-01, -4.351349633402e-01, -3.765224894244e-01 ));
    positions.push_back(Vec3(  2.195852756811e-01,  8.816778044978e-02,  1.073177788976e+00 ));
    positions.push_back(Vec3( -2.899375600289e+00,  4.533801552398e-01,  4.445704119756e-01 ));
    positions.push_back(Vec3( -1.891275615534e+00,  3.753590916398e-01,  2.490296312756e-01 ));
    positions.push_back(Vec3( -3.372062885819e+00,  2.555551963398e-01, -3.733422905244e-01 ));
    positions.push_back(Vec3(  8.428220990811e-01,  4.870807363398e-01,  3.155605857676e+00 ));
    positions.push_back(Vec3(  2.205240762811e-01,  1.158785963140e+00,  3.577474762076e+00 ));
    positions.push_back(Vec3(  1.087801046381e+00, -1.129208591502e-01,  3.875459300076e+00 ));
    positions.push_back(Vec3(  9.587606425811e-01,  1.739780448140e+00, -2.350038144714e+00 ));
    positions.push_back(Vec3(  1.608290357181e+00,  1.084550067240e+00, -2.527800037751e+00 ));
    positions.push_back(Vec3(  6.371470514811e-01,  1.907637977740e+00, -1.459391248124e+00 ));
    positions.push_back(Vec3( -5.402595366189e-01,  3.436858949240e+00, -3.364294798244e-01 ));
    positions.push_back(Vec3( -9.927075362189e-01,  4.233122552140e+00, -8.308930347244e-01 ));
    positions.push_back(Vec3( -1.277172354239e+00,  2.836232717540e+00, -3.850631209244e-01 ));
    positions.push_back(Vec3(  8.739387783811e-01, -1.626958848900e+00, -2.396798577774e+00 ));
    positions.push_back(Vec3(  1.408213796381e+00, -2.079240526160e+00, -1.701669540424e+00 ));
    positions.push_back(Vec3(  1.515572976181e+00, -1.166885340300e+00, -2.938074529504e+00 ));
    positions.push_back(Vec3(  2.976994446581e+00, -4.267919297502e-01,  1.269092190476e+00 ));
    positions.push_back(Vec3(  2.533682331181e+00,  2.556735979978e-02,  2.039154286576e+00 ));
    positions.push_back(Vec3(  3.128503830381e+00,  2.419881502398e-01,  5.967620551756e-01 ));
    positions.push_back(Vec3( -1.629810395359e+00,  8.289601626398e-01, -3.056166050894e+00 ));
    positions.push_back(Vec3( -6.994519097189e-01,  8.049764115398e-01, -2.943175234164e+00 ));
    positions.push_back(Vec3( -1.874145275955e+00, -8.466235253022e-02, -3.241518244384e+00 ));
    positions.push_back(Vec3( -2.595625033599e+00, -1.646832204690e+00,  2.192150970676e+00 ));
    positions.push_back(Vec3( -2.906585612949e+00, -9.062350233512e-01,  1.670815817176e+00 ));
    positions.push_back(Vec3( -2.669295291189e+00, -1.431610358210e+00,  3.153621682076e+00 ));
    positions.push_back(Vec3(  1.818457619381e+00, -2.794187137660e+00,  1.107580890756e-01 ));
    positions.push_back(Vec3(  2.097391780281e+00, -1.959930744260e+00,  5.631032745756e-01 ));
    positions.push_back(Vec3(  2.529698452281e+00, -3.454934434160e+00,  3.013241934756e-01 ));
    positions.push_back(Vec3( -4.545107775189e-01, -3.318742298960e+00,  1.492762161876e+00 ));
    positions.push_back(Vec3(  1.737151477811e-01, -3.171411334960e+00,  8.066594968756e-01 ));
    positions.push_back(Vec3( -1.041679032169e+00, -2.555367771060e+00,  1.600821086376e+00 ));
    positions.push_back(Vec3( -5.833794722189e-01,  3.052540693140e+00,  2.637054256676e+00 ));
    positions.push_back(Vec3( -1.419205322489e+00,  2.586471523640e+00,  2.720532465476e+00 ));
    positions.push_back(Vec3( -5.746468928189e-01,  3.408128243540e+00,  1.785534806076e+00 ));
    positions.push_back(Vec3(  3.344764010581e+00,  1.359764839140e+00, -6.682810128244e-01 ));
    positions.push_back(Vec3(  2.813836912981e+00,  2.168029298140e+00, -5.595663217244e-01 ));
    positions.push_back(Vec3(  4.237583121481e+00,  1.647152701540e+00, -8.700689423244e-01 ));
    positions.push_back(Vec3( -2.004856348130e+00, -1.781589543180e+00, -2.722539611524e+00 ));
    positions.push_back(Vec3( -1.060302432829e+00, -1.994959751260e+00, -2.496572811216e+00 ));
    positions.push_back(Vec3( -2.596156174389e+00, -2.566367260560e+00, -2.697712987614e+00 ));

    // insert null position for virtual site
    std::vector<Vec3>::iterator beginOfNextWater;
    for (int w=0; w<numberOfWaterMolecules; w++) {
            beginOfNextWater = positions.begin() + w * particlesPerMolecule + (particlesPerMolecule-1);
            positions.insert(beginOfNextWater, Vec3( 0., 0., 0.));
    }

    // A => nm
    for (int i=0; i<numberOfParticles; i++) {
            positions[i] *= 1e-1;
    }

    // Setup context
    std::string platformName = "Reference";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );
    context.setPositions(positions);
    context.applyConstraints(1e-6); // update position of virtual sites

    // Run computation
    State state                      = context.getState(State::Forces | State::Energy);
    double energy = state.getPotentialEnergy() / cal2joule;
    std::vector<Vec3> forces         = state.getForces();

    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        forces[ii] /= cal2joule*10;
        if ((ii+1) % particlesPerMolecule == 0) { // Set virtual site force to 0
            forces[ii] *= 0;
        }
    }    

    // Testing

    double tolerance = 1.0e-03;

    std::vector<Vec3> expectedForces;
    expectedForces.push_back(Vec3(  -9.03018602, -25.22018024, -53.51380581 ));
    expectedForces.push_back(Vec3( -12.52482773,  24.89400935,  25.80226200 ));
    expectedForces.push_back(Vec3(  17.02879706,  -3.55643920,  26.68530047 ));
    expectedForces.push_back(Vec3( -36.55910156,  -0.84885078,   6.80733866 ));
    expectedForces.push_back(Vec3(  49.90118398,   1.80065759,  -2.04507653 ));
    expectedForces.push_back(Vec3( -11.18936272,  -2.00312066,  -2.88289342 ));
    expectedForces.push_back(Vec3(  23.42068352, -19.56587056, -27.10935524 ));
    expectedForces.push_back(Vec3( -33.33921978,  25.15157828,  25.27458819 ));
    expectedForces.push_back(Vec3(   6.20288232,  -4.54221160,   7.17792182 ));
    expectedForces.push_back(Vec3(  24.85926509, -33.25263883,  25.99134252 ));
    expectedForces.push_back(Vec3( -13.17758428,  12.09180034, -25.71547223 ));
    expectedForces.push_back(Vec3( -12.83642431,  17.04564612,  -2.37510580 ));
    expectedForces.push_back(Vec3(  32.21731121, -64.94236003,  33.92641093 ));
    expectedForces.push_back(Vec3( -40.11559851,  49.01185816, -32.26359729 ));
    expectedForces.push_back(Vec3(   8.79597794,  21.26568970,  -2.02614751 ));
    expectedForces.push_back(Vec3( -13.79875732,   7.94084709, -21.43603790 ));
    expectedForces.push_back(Vec3(  11.96232403,  -4.89564516,  11.75815846 ));
    expectedForces.push_back(Vec3(   2.32126747,   2.38141894,   3.09469541 ));
    expectedForces.push_back(Vec3(  21.91680081,   3.37395463, -36.63017621 ));
    expectedForces.push_back(Vec3( -14.91664520,  12.66435924,  31.71106001 ));
    expectedForces.push_back(Vec3(  -0.68988308,  -8.07541704,   1.68094803 ));
    expectedForces.push_back(Vec3(  31.08549012,  -4.09618771,   6.41621800 ));
    expectedForces.push_back(Vec3( -29.01317556,  -6.29083303,  -2.86799406 ));
    expectedForces.push_back(Vec3(   2.50721872,  -2.02990166,  -4.01517086 ));
    expectedForces.push_back(Vec3(  -4.63075645,  10.83217003, -41.17139179 ));
    expectedForces.push_back(Vec3(   4.37795918, -15.42375203,   2.85804473 ));
    expectedForces.push_back(Vec3(   0.47708162,   2.90833544,  37.27058864 ));
    expectedForces.push_back(Vec3( -32.52126814,   1.18150031,  -0.50441484 ));
    expectedForces.push_back(Vec3(  -4.01706902,  16.11448974,   2.85777915 ));
    expectedForces.push_back(Vec3(  24.74097363, -25.17193723,   4.77221810 ));
    expectedForces.push_back(Vec3(  21.55479185,   7.98770439, -24.32863803 ));
    expectedForces.push_back(Vec3( -11.76596457,  -9.89888550,  12.91177869 ));
    expectedForces.push_back(Vec3(  -3.02047642,   6.21275158,   4.97109387 ));
    expectedForces.push_back(Vec3(   3.48368086,  24.91392899, -56.33474548 ));
    expectedForces.push_back(Vec3(   0.48324390,  -4.99402937,   4.84782778 ));
    expectedForces.push_back(Vec3(  -1.36307280, -17.27929672,  51.34101302 ));
    expectedForces.push_back(Vec3(   4.86574918, -13.78940550,   2.82892934 ));
    expectedForces.push_back(Vec3(  -9.64565940,  10.90797883,   4.81889484 ));
    expectedForces.push_back(Vec3(   1.75254400,  -2.20460013,  -1.03773916 ));
    expectedForces.push_back(Vec3(  -4.45254552,  23.26887634,  -5.28665816 ));
    expectedForces.push_back(Vec3(  34.72253327,   4.24026292,   8.12021875 ));
    expectedForces.push_back(Vec3( -30.07018137, -18.10825505,  -2.38021111 ));

    // insert null position for virtual site
    for (int w=0; w<numberOfWaterMolecules; w++) {
            beginOfNextWater = expectedForces.begin() + w * particlesPerMolecule + (particlesPerMolecule-1);
            expectedForces.insert(beginOfNextWater, Vec3( 0., 0., 0.));
    }


    // gradients => forces
    for( unsigned int ii = 0; ii < expectedForces.size(); ii++ ){
        expectedForces[ii] *= -1;
    }
    const double expectedEnergy = -59.10662152;

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
        std::cout << "testMBPolWater14Test running test..." << std::endl;

        testWater14();

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
