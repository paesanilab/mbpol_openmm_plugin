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

    std::vector<double> thole(5);

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


    } else {
        mbpolElectrostaticsForce->setNonbondedMethod( MBPolElectrostaticsForce::NoCutoff );

        mbpolTwoBodyForce->setNonbondedMethod(MBPolTwoBodyForce::CutoffNonPeriodic);

        mbpolThreeBodyForce->setNonbondedMethod(MBPolThreeBodyForce::CutoffNonPeriodic);

    }

    int numberOfWaterMolecules = 50;
    unsigned int particlesPerMolecule = 4;
    int numberOfParticles          = numberOfWaterMolecules * particlesPerMolecule;

    std::vector<int> particleIndices(particlesPerMolecule);
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

        mbpolElectrostaticsForce->addElectrostatics( -5.1966000e-01, jj+1, jj+2, jj+3,
                                            thole, 0.001310, 0.001310 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+2, jj+3,
                                            thole, 0.000294, 0.000294 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+1, jj+3,
                                            thole, 0.000294, 0.000294 );
        mbpolElectrostaticsForce->addElectrostatics(  0., jj, jj+1, jj+2,
                                                    thole,  0.001310,  0.);

        mbpolOneBodyForce->addOneBody(particleIndices);
        mbpolTwoBodyForce->addParticle( particleIndices);
        mbpolThreeBodyForce->addParticle( particleIndices);
    }

    system.addForce(mbpolElectrostaticsForce);
    system.addForce(mbpolOneBodyForce);
    system.addForce(mbpolTwoBodyForce);
    system.addForce(mbpolThreeBodyForce);

    LangevinIntegrator integrator(0.0, 0.1, 0.01);

    std::vector<Vec3> positions;
    std::vector<Vec3> expectedForces;
    double expectedEnergy;

    positions.push_back(Vec3( 2.85736094, -5.55081227, -8.12183429));
    positions.push_back(Vec3( 2.44281114, -6.43473556, -8.07382840));
    positions.push_back(Vec3( 3.82473974, -5.70456735, -8.19899732));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-6.79525511,  4.75911250, -3.88813105));
    positions.push_back(Vec3(-7.35706790,  5.55740118, -3.95831371));
    positions.push_back(Vec3(-5.87912662,  5.16802048, -3.91717400));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-6.82404288,  5.61379358,  1.27419175));
    positions.push_back(Vec3(-7.69644118,  5.77578295,  1.50988331));
    positions.push_back(Vec3(-6.26169320,  6.01020910,  1.94248670));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 8.67093184, -1.30608328, -0.78315366));
    positions.push_back(Vec3(-8.41988721, -0.93466363, -0.94475345));
    positions.push_back(Vec3( 8.07488886, -0.52795751, -0.86983334));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-3.11505935,  1.26377015,  0.77958938));
    positions.push_back(Vec3(-2.39788475,  0.72018722,  0.38979635));
    positions.push_back(Vec3(-2.94280463,  2.10811114,  0.41140302));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-8.51671124, -5.76530765,  7.74108808));
    positions.push_back(Vec3(-7.87908408, -5.49804129,  8.41407895));
    positions.push_back(Vec3(-8.20968797, -6.66268648,  7.64028143));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 5.72789613,  3.42420549,  7.41025239));
    positions.push_back(Vec3( 4.79017835,  3.57078919,  7.05138027));
    positions.push_back(Vec3( 6.11411964,  2.78476544,  6.78143304));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 8.76294937,  5.89517641,  2.05903411));
    positions.push_back(Vec3( 8.09068002,  5.41991739,  1.55581783));
    positions.push_back(Vec3( 8.51030487,  5.58693918,  2.93583106));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-4.01902800,  2.58990689,  6.43969265));
    positions.push_back(Vec3(-3.80830690,  1.68932165,  6.61639957));
    positions.push_back(Vec3(-3.76627940,  2.64611575,  5.55292014));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-1.47564412, -4.55130620,  6.57499487));
    positions.push_back(Vec3(-1.69328116, -4.87917171,  5.69599997));
    positions.push_back(Vec3(-2.20225769, -3.98711832,  6.87550629));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-1.00583365, -4.25714790, -7.72257630));
    positions.push_back(Vec3(-1.74254870, -4.50830189, -8.27768838));
    positions.push_back(Vec3(-0.22387925, -4.28434773, -8.30970008));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 7.63942120,  5.68542007,  7.48493667));
    positions.push_back(Vec3( 8.12635201,  5.62356557,  8.29358686));
    positions.push_back(Vec3( 6.98714092,  4.99478707,  7.60070111));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-2.57445193, -2.27550335, -6.82016356));
    positions.push_back(Vec3(-3.01048097, -2.46346312, -6.01645354));
    positions.push_back(Vec3(-2.06302609, -3.04494539, -7.10580361));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-0.73809516, -3.26261303,  2.29195696));
    positions.push_back(Vec3(-1.41699363, -2.59985513,  2.37577931));
    positions.push_back(Vec3(-0.28243934, -2.99088243,  1.47665394));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 8.04669592, -3.53299536,  7.25871135));
    positions.push_back(Vec3( 8.64762436, -4.28342494,  7.18787287));
    positions.push_back(Vec3( 7.21966794, -3.86200623,  6.83274460));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 4.73231920,  0.89945619,  6.00518809));
    positions.push_back(Vec3( 4.16333424,  0.19746055,  5.74015493));
    positions.push_back(Vec3( 4.85350849,  0.90199381,  6.95332104));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 0.38132916, -6.48615762, -0.51778543));
    positions.push_back(Vec3( 1.18897024, -6.58841661, -0.04769810));
    positions.push_back(Vec3(-0.31470223, -6.42686269,  0.13477770));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 3.19664593,  3.90293994,  6.31179027));
    positions.push_back(Vec3( 3.31465591,  3.31689070,  5.51699030));
    positions.push_back(Vec3( 2.35184523,  3.65650511,  6.66857253));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 3.24242536, -0.29510218,  2.69252789));
    positions.push_back(Vec3( 2.90652749, -1.15200217,  3.05611355));
    positions.push_back(Vec3( 3.86564281, -0.50059484,  2.00820199));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 3.15568037,  2.42073277,  3.70930928));
    positions.push_back(Vec3( 3.46473047,  2.94011188,  2.96806037));
    positions.push_back(Vec3( 3.18478646,  1.48926347,  3.36543611));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 1.29878135, -5.16673665, -5.46192296));
    positions.push_back(Vec3( 0.52205982, -4.78776156, -5.95162466));
    positions.push_back(Vec3( 1.99200780, -5.07662460, -6.12619141));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-3.35585795, -3.72637017,  8.63303497));
    positions.push_back(Vec3(-4.27885756, -4.01574704,  8.55707282));
    positions.push_back(Vec3(-3.29851093, -3.21760072, -8.55333993));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-2.19671563, -5.35826327,  3.93258988));
    positions.push_back(Vec3(-2.11836651, -6.19616924,  3.44496059));
    positions.push_back(Vec3(-1.68674034, -4.72187859,  3.46384511));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-6.65899761,  0.53684591, -1.29708545));
    positions.push_back(Vec3(-5.76638371,  0.31433734, -1.48828204));
    positions.push_back(Vec3(-6.68406703,  1.50756159, -1.07335310));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-6.85220610,  3.30828341, -1.11118607));
    positions.push_back(Vec3(-6.77567044,  4.01628835, -0.43787068));
    positions.push_back(Vec3(-6.75505486,  3.71757689, -1.96491732));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-8.86011346, -4.52978939, -6.61011395));
    positions.push_back(Vec3(-8.43019795, -3.83953486, -6.14510572));
    positions.push_back(Vec3( 8.39682542, -4.10703722, -7.00192033));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 8.67471742,  5.23444868,  4.54951834));
    positions.push_back(Vec3(-8.86362133,  4.41579768,  4.73810684));
    positions.push_back(Vec3( 8.35413658,  5.53522809,  5.38330475));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-8.16826797, -3.03109217,  2.23475195));
    positions.push_back(Vec3(-7.22439594, -2.70017827,  2.18395484));
    positions.push_back(Vec3(-8.25867610, -3.88391387,  1.77975103));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-5.37464373,  5.27828675, -6.60151544));
    positions.push_back(Vec3(-5.03051536,  5.12193108, -7.53165144));
    positions.push_back(Vec3(-6.07214139,  4.66161835, -6.52377251));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-1.15863413, -0.02712881, -0.74645475));
    positions.push_back(Vec3(-0.33640036, -0.09450170, -1.28835651));
    positions.push_back(Vec3(-1.51349409, -0.99287082, -0.71105594));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 7.26135968, -2.93403499, -7.98047603));
    positions.push_back(Vec3( 7.54629570, -2.86994664, -8.97100092));
    positions.push_back(Vec3( 6.42033007, -2.53152624, -8.10296553));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 7.20753665,  8.67553343, -8.83592966));
    positions.push_back(Vec3( 7.83620503,  8.19764654, -8.25444552));
    positions.push_back(Vec3( 6.88465000,  8.11352021,  8.46700018));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-2.15712544, -0.64664741,  2.92772353));
    positions.push_back(Vec3(-1.59212430, -0.11061402,  3.53909752));
    positions.push_back(Vec3(-2.72898808,  0.01961931,  2.44685316));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 2.81698367, -6.91414203,  1.18924772));
    positions.push_back(Vec3( 3.39819981, -6.14801631,  1.37953093));
    positions.push_back(Vec3( 3.09238980, -7.58400287,  1.78922348));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 4.81467328, -4.62899093,  1.45423079));
    positions.push_back(Vec3( 4.59422750, -3.73422731,  1.63452067));
    positions.push_back(Vec3( 5.61507004, -4.42579068,  0.99596812));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 8.75794610, -5.49447558,  1.25168066));
    positions.push_back(Vec3( 8.20341244, -4.78594023,  0.87522915));
    positions.push_back(Vec3( 8.64547410, -6.28986479,  0.68477902));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 3.57214060,  1.66971234, -2.58053403));
    positions.push_back(Vec3( 3.02224731,  2.33710707, -2.07786247));
    positions.push_back(Vec3( 3.55804955,  1.90577853, -3.49525766));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-3.90044231,  4.10603895,  8.82889889));
    positions.push_back(Vec3(-2.96877974,  3.87182172, -8.86405099));
    positions.push_back(Vec3(-3.94908288,  3.60455914,  8.01354917));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 7.39383381, -3.06983109,  0.58970994));
    positions.push_back(Vec3( 7.93634144, -2.67087773,  1.30787959));
    positions.push_back(Vec3( 7.73328228, -2.59664226, -0.16548384));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 6.30897456, -6.01026482, -8.11560112));
    positions.push_back(Vec3( 6.82392069, -5.21489272, -8.31120577));
    positions.push_back(Vec3( 6.90621900, -6.68609686, -8.41239908));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-2.85739773,  0.23550221,  6.62437325));
    positions.push_back(Vec3(-3.37302485, -0.52315639,  6.29505802));
    positions.push_back(Vec3(-1.93616178,  0.00125324,  6.64684388));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-6.88532817, -4.61025442, -8.44015108));
    positions.push_back(Vec3(-6.78261563, -3.65960508, -8.41085663));
    positions.push_back(Vec3(-7.64357546, -4.76709544, -7.78712038));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 1.58560591, -0.08666679, -1.90052646));
    positions.push_back(Vec3( 1.90678915, -0.95213307, -1.67008498));
    positions.push_back(Vec3( 2.29581362,  0.27373393, -2.39548101));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-1.68382787, -2.70280182, -0.28470772));
    positions.push_back(Vec3(-1.80337488, -3.28352879, -1.06707946));
    positions.push_back(Vec3(-2.58550155, -2.70484793,  0.09078284));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-0.88609249, -6.39591875, -3.23511204));
    positions.push_back(Vec3(-0.31246175, -6.52063484, -2.42324532));
    positions.push_back(Vec3(-0.51189779, -6.91109922, -3.92189582));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 1.54361367, -2.30440442,  3.77268499));
    positions.push_back(Vec3( 2.00934692, -2.59544066,  4.56358258));
    positions.push_back(Vec3( 0.69560000, -2.73854999,  3.76637164));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 0.82205203, -3.84463700,  8.31466504));
    positions.push_back(Vec3( 0.32816106, -4.13197608,  7.59224083));
    positions.push_back(Vec3( 1.53033815, -4.48116252,  8.43551439));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-8.46624219,  5.76863416, -8.38053447));
    positions.push_back(Vec3(-7.90807629,  6.28385825, -7.77168073));
    positions.push_back(Vec3(-7.83748479,  5.32961922, -8.97240869));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3(-1.71399162, -3.92050611, -2.82231781));
    positions.push_back(Vec3(-1.36022862, -4.72786151, -3.29342441));
    positions.push_back(Vec3(-2.23738279, -3.40325836, -3.41361371));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));
    positions.push_back(Vec3( 1.46313225,  2.55586533, -0.87404907));
    positions.push_back(Vec3( 0.81502005,  3.14063862, -0.50488322));
    positions.push_back(Vec3( 1.10224998,  1.68909320, -0.78617541));
    positions.push_back(Vec3( 0.00000000,  0.00000000,  0.00000000));



    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            positions[i][j] *= 1e-1;
        }
    }

    // gradients => forces
    for( unsigned int ii = 0; ii < expectedForces.size(); ii++ ){
        expectedForces[ii] *= -1;
    }
    expectedEnergy        = -244.37507;

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

    vector<std::string> forceLabels;
    forceLabels.push_back("Electrostatics");
    forceLabels.push_back("OneBody");
    forceLabels.push_back("TwoBody");
    forceLabels.push_back("ThreeBody");
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

        double boxDimension = 1.8;
        runTest( boxDimension );

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
