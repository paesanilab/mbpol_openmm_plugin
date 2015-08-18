/* -------------------------------------------------------------------------- *
 *                                   OpenMMMBPol                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs                                                     *
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
 * This tests the Reference implementation of ReferenceMBPolTwoBodyForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMMBPol.h"
#include "openmm/System.h"
#include "openmm/MBPolTwoBodyForce.h"
#include "MBPolReferenceTwoBodyForce.h"
#include "openmm/LangevinIntegrator.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#include "openmm/reference/RealVec.h"

#define ASSERT_EQUAL_TOL_MOD(expected, found, tol, testname) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC_MOD(expected, found, tol,testname) {ASSERT_EQUAL_TOL_MOD((expected)[0], (found)[0], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[1], (found)[1], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[2], (found)[2], (tol),(testname));};

using namespace OpenMM;
using namespace MBPolPlugin;

const double TOL = 1e-4;\

void testTwoWaterOneChloride(double boxDimension, bool addPositionOffset) {
	std::string testName = "testMBPol2BodyInteraction";

	System system;
	int numberOfParticles = 7;
	MBPolTwoBodyForce* mbpolTwoBodyForce = new MBPolTwoBodyForce();
	double cutoff = 10;
	mbpolTwoBodyForce->setCutoff(cutoff);

	unsigned int particlesPerMolecule = 3;

	if (boxDimension > 0.0) {
		Vec3 a(boxDimension, 0.0, 0.0);
		Vec3 b(0.0, boxDimension, 0.0);
		Vec3 c(0.0, 0.0, boxDimension);
		system.setDefaultPeriodicBoxVectors(a, b, c);
		mbpolTwoBodyForce->setNonbondedMethod(
				MBPolTwoBodyForce::CutoffPeriodic);
	} else {
		mbpolTwoBodyForce->setNonbondedMethod(
				MBPolTwoBodyForce::CutoffNonPeriodic);
	}
	//add all water
	std::vector<int> particleIndices(particlesPerMolecule);
	for (unsigned int jj = 0; jj < numberOfParticles - 1; jj +=
			particlesPerMolecule) {
		system.addParticle(1.5999000e+01);
		system.addParticle(1.0080000e+00);
		system.addParticle(1.0080000e+00);
		particleIndices[0] = jj;
		particleIndices[1] = jj + 1;
		particleIndices[2] = jj + 2;
		mbpolTwoBodyForce->addParticle(particleIndices);
	}
	//add chloride
	system.addParticle(35);
	particleIndices.resize(1);
	particleIndices[0] = numberOfParticles - 1;

	mbpolTwoBodyForce->addParticle(particleIndices);

	LangevinIntegrator integrator(0.0, 0.1, 0.01);

	std::vector<Vec3> positions(numberOfParticles);
	std::vector<Vec3> expectedForces(numberOfParticles);
	double expectedEnergy;

	positions[Oa]  = Vec3( 0.2071123419,        0.1431350648,        0.2575948809);
	positions[Ha1] = Vec3( 0.3968154394,       -0.4466985822,       -0.5164780569);
	positions[Ha2] = Vec3(-0.0790971024,        0.9461358941,       -0.1868734592);

	positions[Ob]  = Vec3( 3.1469464639,        0.2657051814,       -0.5347590767);
	positions[Hb1] = Vec3( 2.3878654247,        0.3762497422,        0.0532565299);
	positions[Hb2] = Vec3( 2.7252172819,       -0.2612759028,       -1.2421153352);

	positions[Cl2] = Vec3( 1.0072476456,       -1.3374671605,       -2.3154322027);

	for (int i = 0; i < numberOfParticles; i++) {
		for (int j = 0; j < 3; j++) {
			positions[i][j] *= 1e-1;
		}
	}

	expectedEnergy = 19.6295;

	system.addForce(mbpolTwoBodyForce);

	std::string platformName;
#define AngstromToNm 0.1
#define CalToJoule   4.184

	platformName = "Reference";
	Context context(system, integrator,
			Platform::getPlatformByName(platformName));

	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	std::vector<Vec3> forces = state.getForces();

	for (unsigned int ii = 0; ii < forces.size(); ii++) {
		forces[ii] /= CalToJoule * 10;
		expectedForces[ii] *= -1; // gradient -> forces
	}

	double tolerance = 1.0e-03;

	double energy = state.getPotentialEnergy() / CalToJoule;

	std::cout << "Energy: " << energy << " Kcal/mol " << std::endl;
	std::cout << "Expected energy: " << expectedEnergy << " Kcal/mol "
			<< std::endl;

	std::cout << "Comparison of energy with tolerance: " << tolerance
			<< std::endl << std::endl;

	ASSERT_EQUAL_TOL(expectedEnergy, energy, tolerance);

	std::cout << "Test Successful: " << testName << std::endl << std::endl;

}

void testTwoBodyChloride(double boxDimension, bool addPositionOffset) {
	std::string testName = "testMBPol2BodyInteraction";

	System system;
	int numberOfParticles = 4;
	MBPolTwoBodyForce* mbpolTwoBodyForce = new MBPolTwoBodyForce();
	double cutoff = 10;
	mbpolTwoBodyForce->setCutoff(cutoff);

	unsigned int particlesPerMolecule = 3;

	if (boxDimension > 0.0) {
		Vec3 a(boxDimension, 0.0, 0.0);
		Vec3 b(0.0, boxDimension, 0.0);
		Vec3 c(0.0, 0.0, boxDimension);
		system.setDefaultPeriodicBoxVectors(a, b, c);
		mbpolTwoBodyForce->setNonbondedMethod(
				MBPolTwoBodyForce::CutoffPeriodic);
	} else {
		mbpolTwoBodyForce->setNonbondedMethod(
				MBPolTwoBodyForce::CutoffNonPeriodic);
	}
	//add all water
	std::vector<int> particleIndices(particlesPerMolecule);
	for (unsigned int jj = 0; jj < numberOfParticles - 1; jj +=
			particlesPerMolecule) {
		system.addParticle(1.5999000e+01);
		system.addParticle(1.0080000e+00);
		system.addParticle(1.0080000e+00);
		particleIndices[0] = jj;
		particleIndices[1] = jj + 1;
		particleIndices[2] = jj + 2;
		mbpolTwoBodyForce->addParticle(particleIndices);
	}
	//add chloride
	system.addParticle(35);
	particleIndices.resize(1);
	particleIndices[0] = numberOfParticles - 1;
//    particleIndices[1] = -1;
//    particleIndices[2] = -1;
	mbpolTwoBodyForce->addParticle(particleIndices);

	LangevinIntegrator integrator(0.0, 0.1, 0.01);

	std::vector<Vec3> positions(numberOfParticles);
	std::vector<Vec3> expectedForces(numberOfParticles);
	double expectedEnergy;

	positions[0] = Vec3(6.028800000e+00, -4.021520000e+00, 1.384500000e+00);
	positions[1] = Vec3(5.895700000e+00, -3.282780000e+00, 7.820000000e-01);
	positions[2] = Vec3(6.859700000e+00, -4.452620000e+00, 1.204600000e+00);

	positions[3] = Vec3(5.004300000e+00, -1.784700000e+00, 3.409400000e+00);

	for (int i = 0; i < numberOfParticles; i++) {
		for (int j = 0; j < 3; j++) {
			positions[i][j] *= 1e-1;
		}
	}

//    if (addPositionOffset) {
//        // move second molecule 1 box dimension in Y direction
//        positions[3][1] += boxDimension;
//        positions[4][1] += boxDimension;
//        positions[5][1] += boxDimension;
//    }

	expectedForces[0] = Vec3(-3.278090251e+00, 4.781931804e+00,
			2.443283214e+00);
	expectedForces[1] = Vec3(5.700852530e-01, -1.298619977e+00,
			1.232417788e+00);
	expectedForces[2] = Vec3(2.181168186e-01, 1.161499952e+00, 5.698707152e-01);
	expectedForces[3] = Vec3(2.489888179e+00, -4.644811779e+00,
			-4.245571717e+00);

	expectedEnergy = 2.027705104;

	system.addForce(mbpolTwoBodyForce);

	std::string platformName;
#define AngstromToNm 0.1
#define CalToJoule   4.184

	platformName = "Reference";
	Context context(system, integrator,
			Platform::getPlatformByName(platformName));

	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	std::vector<Vec3> forces = state.getForces();

	for (unsigned int ii = 0; ii < forces.size(); ii++) {
		forces[ii] /= CalToJoule * 10;
		expectedForces[ii] *= -1; // gradient -> forces
	}

	double tolerance = 1.0e-03;

	double energy = state.getPotentialEnergy() / CalToJoule;

	std::cout << "Energy: " << energy << " Kcal/mol " << std::endl;
	std::cout << "Expected energy: " << expectedEnergy << " Kcal/mol "
			<< std::endl;

	std::cout << std::endl << "Forces:" << std::endl;

	for (int i = 0; i < numberOfParticles; i++) {
		std::cout << "Force atom " << i << ": " << expectedForces[i]
				<< " Kcal/mol/A <mbpol>" << std::endl;
		std::cout << "Force atom " << i << ": " << forces[i]
				<< " Kcal/mol/A <openmm-mbpol>" << std::endl << std::endl;
	}

	std::cout << "Comparison of energy and forces with tolerance: " << tolerance
			<< std::endl << std::endl;

	ASSERT_EQUAL_TOL(expectedEnergy, energy, tolerance);

	for (unsigned int ii = 0; ii < forces.size(); ii++) {
		ASSERT_EQUAL_VEC(expectedForces[ii], forces[ii], tolerance);
	}
	std::cout << "Test Successful: " << testName << std::endl << std::endl;

}

void testTwoBody(double boxDimension, bool addPositionOffset) {

	std::string testName = "testMBPol2BodyInteraction";

	System system;
	int numberOfParticles = 6;
	MBPolTwoBodyForce* mbpolTwoBodyForce = new MBPolTwoBodyForce();
	double cutoff = 10;
	mbpolTwoBodyForce->setCutoff(cutoff);

	unsigned int particlesPerMolecule = 3;

	if (boxDimension > 0.0) {
		Vec3 a(boxDimension, 0.0, 0.0);
		Vec3 b(0.0, boxDimension, 0.0);
		Vec3 c(0.0, 0.0, boxDimension);
		system.setDefaultPeriodicBoxVectors(a, b, c);
		mbpolTwoBodyForce->setNonbondedMethod(
				MBPolTwoBodyForce::CutoffPeriodic);
	} else {
		mbpolTwoBodyForce->setNonbondedMethod(
				MBPolTwoBodyForce::CutoffNonPeriodic);
	}

	std::vector<int> particleIndices(particlesPerMolecule);
	for (unsigned int jj = 0; jj < numberOfParticles; jj +=
			particlesPerMolecule) {
		system.addParticle(1.5999000e+01);
		system.addParticle(1.0080000e+00);
		system.addParticle(1.0080000e+00);
		particleIndices[0] = jj;
		particleIndices[1] = jj + 1;
		particleIndices[2] = jj + 2;
		mbpolTwoBodyForce->addParticle(particleIndices);
	}

	LangevinIntegrator integrator(0.0, 0.1, 0.01);

	std::vector<Vec3> positions(numberOfParticles);
	std::vector<Vec3> expectedForces(numberOfParticles);
	double expectedEnergy;
	positions[0] = Vec3(-1.516074336e+00, -2.023167650e-01, 1.454672917e+00);
	positions[1] = Vec3(-6.218989773e-01, -6.009430735e-01, 1.572437625e+00);
	positions[2] = Vec3(-2.017613812e+00, -4.190350349e-01, 2.239642849e+00);

	positions[3] = Vec3(-1.763651687e+00, -3.816594649e-01, -1.300353949e+00);
	positions[4] = Vec3(-1.903851736e+00, -4.935677617e-01, -3.457810126e-01);
	positions[5] = Vec3(-2.527904158e+00, -7.613550077e-01, -1.733803676e+00);

	for (int i = 0; i < numberOfParticles; i++) {
		for (int j = 0; j < 3; j++) {
			positions[i][j] *= 1e-1;
		}
	}

	if (addPositionOffset) {
		// move second molecule 1 box dimension in Y direction
		positions[3][1] += boxDimension;
		positions[4][1] += boxDimension;
		positions[5][1] += boxDimension;
	}

	expectedForces[0] = Vec3(-4.85337479, -4.47836379, -20.08989563);
	expectedForces[1] = Vec3(-0.31239868, 0.52518586, -1.88893830);
	expectedForces[2] = Vec3(0.00886712, 0.73323536, -1.81715325);
	expectedForces[3] = Vec3(-0.65181727, -0.72947395, 5.88973293);
	expectedForces[4] = Vec3(4.82340981, 3.20090213, 16.49522051);
	expectedForces[5] = Vec3(0.98531382, 0.74851439, 1.41103374);

	expectedEnergy = 6.14207815;

	system.addForce(mbpolTwoBodyForce);
	std::string platformName;
#define AngstromToNm 0.1
#define CalToJoule   4.184

	platformName = "Reference";
	Context context(system, integrator,
			Platform::getPlatformByName(platformName));

	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	std::vector<Vec3> forces = state.getForces();

	for (unsigned int ii = 0; ii < forces.size(); ii++) {
		forces[ii] /= CalToJoule * 10;
		expectedForces[ii] *= -1; // gradient -> forces
	}

	double tolerance = 1.0e-03;

	double energy = state.getPotentialEnergy() / CalToJoule;

	std::cout << "Energy: " << energy << " Kcal/mol " << std::endl;
	std::cout << "Expected energy: " << expectedEnergy << " Kcal/mol "
			<< std::endl;

	std::cout << std::endl << "Forces:" << std::endl;

	for (int i = 0; i < numberOfParticles; i++) {
		std::cout << "Force atom " << i << ": " << expectedForces[i]
				<< " Kcal/mol/A <mbpol>" << std::endl;
		std::cout << "Force atom " << i << ": " << forces[i]
				<< " Kcal/mol/A <openmm-mbpol>" << std::endl << std::endl;
	}

	std::cout << "Comparison of energy and forces with tolerance: " << tolerance
			<< std::endl << std::endl;

	ASSERT_EQUAL_TOL(expectedEnergy, energy, tolerance);

	for (unsigned int ii = 0; ii < forces.size(); ii++) {
		ASSERT_EQUAL_VEC(expectedForces[ii], forces[ii], tolerance);
	}
	std::cout << "Test Successful: " << testName << std::endl << std::endl;

}

void testImageMolecules(bool runTestWithAtomImaging, bool addPositionOffset) {

	double boxDimension = 10.;
	RealVec box(boxDimension, boxDimension, boxDimension);

	unsigned int numberOfParticles = 6;
	std::vector<RealVec> particlePositions(numberOfParticles);

	particlePositions[0] = RealVec(0., 0., 0.);
	particlePositions[1] = RealVec(0., 1., 0.);
	particlePositions[2] = RealVec(0., 0., 0.);

	particlePositions[3] = RealVec(4.5, 0., 0.);
	particlePositions[4] = RealVec(-4.5, 0., 0.);
	particlePositions[5] = RealVec(4.5, 0., 0.);

	std::vector<RealVec> expectedParticlePositions(numberOfParticles);

	expectedParticlePositions[0] = RealVec(0., 0., 0.);
	expectedParticlePositions[1] = RealVec(0., 1., 0.);
	expectedParticlePositions[2] = RealVec(0., 0., 0.);

	expectedParticlePositions[3] = RealVec(4.5, 0., 0.);
	expectedParticlePositions[4] = RealVec(5.5, 0., 0.);
	expectedParticlePositions[5] = RealVec(4.5, 0., 0.);

	std::vector<RealVec> originalParticlePositions(numberOfParticles);
	originalParticlePositions = particlePositions;

	if (addPositionOffset) {
		// move second molecule 1 box dimension in Y direction
		particlePositions[3][0] += boxDimension;
		particlePositions[4][0] += boxDimension;
		particlePositions[5][0] += boxDimension;
	}

	std::vector<std::vector<int> > allParticleIndices(6);

	allParticleIndices[0].push_back(0);
	allParticleIndices[0].push_back(1);
	allParticleIndices[0].push_back(2);

	allParticleIndices[3].push_back(3);
	allParticleIndices[3].push_back(4);
	allParticleIndices[3].push_back(5);

	double imagedPositions[numberOfParticles * 2];
	imageMolecules(box, particlePositions);

	if (runTestWithAtomImaging) {
		// Check that making periodic images of everything with respect to the first oxygen fails
		imageParticles(box, particlePositions[0], particlePositions[4]);
		imageParticles(box, particlePositions[0], particlePositions[5]);
	}

	RealVec tempPosition;
	for (int i = 0; i < numberOfParticles; i++) {
		std::cout << "Position atom " << i << ": "
				<< expectedParticlePositions[i] << " A" << std::endl;
		std::cout << "Position atom " << i << ": " << particlePositions[i]
				<< " A" << std::endl;
	}

	for (int i = 0; i < numberOfParticles; i++) {
		ASSERT_EQUAL_VEC(expectedParticlePositions[i], particlePositions[i],
				1e-6);
	}
}

void testImageMoleculesOneWaterChloride(bool runTestWithAtomImaging,
		bool addPositionOffset) {

	double boxDimension = 10.;
	RealVec box(boxDimension, boxDimension, boxDimension);

	unsigned int numberOfParticles = 4;
	std::vector<RealVec> particlePositions(numberOfParticles);

	particlePositions[0] = RealVec(0., 1., 0.);
	particlePositions[1] = RealVec(0., 0., 0.);
	particlePositions[2] = RealVec(0., 2., 0.);

	particlePositions[3] = RealVec(3, -4.5, 0.);

	std::vector<RealVec> expectedParticlePositions(numberOfParticles);

	expectedParticlePositions[0] = RealVec(0., 1., 0.);
	expectedParticlePositions[1] = RealVec(0., 0., 0.);
	expectedParticlePositions[2] = RealVec(0., 2., 0.);

	expectedParticlePositions[3] = RealVec(3, 5.5, 0.);

	if (addPositionOffset) {
		// move second molecule 1 box dimension in Y direction
		particlePositions[3][0] += boxDimension;
	}

	double imagedPositions[numberOfParticles * 2];
	imageMolecules(box, particlePositions);

	if (runTestWithAtomImaging) {
		// Check that making periodic images of everything with respect to the first oxygen fails
		imageParticles(box, particlePositions[0], particlePositions[4]);
		imageParticles(box, particlePositions[0], particlePositions[5]);
	}

	RealVec tempPosition;
	for (int i = 0; i < numberOfParticles; i++) {
		std::cout << "Position atom " << i << ": "
				<< expectedParticlePositions[i] << " A" << std::endl;
		std::cout << "Position atom " << i << ": " << particlePositions[i]
				<< " A" << std::endl;
	}

	for (int i = 0; i < numberOfParticles; i++) {
		ASSERT_EQUAL_VEC(expectedParticlePositions[i], particlePositions[i],
				1e-6);
	}
}

void testImageMoleculesTwoWaterChloride(bool runTestWithAtomImaging,
		bool addPositionOffset) {

	double boxDimension = 10.;
	RealVec box(boxDimension, boxDimension, boxDimension);

	unsigned int numberOfParticles = 7;
	std::vector<RealVec> particlePositions(numberOfParticles);

	particlePositions[0] = RealVec(0., 1., 0.);
	particlePositions[1] = RealVec(0., 0., 0.);
	particlePositions[2] = RealVec(0., 2., 0.);

	particlePositions[3] = RealVec(0, -4.5, 0.);
	particlePositions[4] = RealVec(-1, -4.5, 0.);
	particlePositions[5] = RealVec(1, -4.5, 0.);

	particlePositions[6] = RealVec(3, -4.5, 0.);

	std::vector<RealVec> expectedParticlePositions(numberOfParticles);

	expectedParticlePositions[0] = RealVec(0., 1., 0.);
	expectedParticlePositions[1] = RealVec(0., 0., 0.);
	expectedParticlePositions[2] = RealVec(0., 2., 0.);

	expectedParticlePositions[3] = RealVec(0, 5.5, 0.);
	expectedParticlePositions[4] = RealVec(-1, 5.5, 0.);
	expectedParticlePositions[5] = RealVec(1, 5.5, 0.);

	expectedParticlePositions[6] = RealVec(3, 5.5, 0.);

	if (addPositionOffset) {
		// move second molecule 1 box dimension in Y direction
		particlePositions[3][0] += boxDimension;
		particlePositions[4][0] += boxDimension;
		particlePositions[5][0] += boxDimension;
	}

	double imagedPositions[numberOfParticles * 2];
	imageMolecules(box, particlePositions);

	if (runTestWithAtomImaging) {
		// Check that making periodic images of everything with respect to the first oxygen fails
		imageParticles(box, particlePositions[0], particlePositions[4]);
		imageParticles(box, particlePositions[0], particlePositions[5]);
	}

	RealVec tempPosition;
	for (int i = 0; i < numberOfParticles; i++) {
		std::cout << "Position atom " << i << ": "
				<< expectedParticlePositions[i] << " A" << std::endl;
		std::cout << "Position atom " << i << ": " << particlePositions[i]
				<< " A" << std::endl;
	}

	for (int i = 0; i < numberOfParticles; i++) {
		ASSERT_EQUAL_VEC(expectedParticlePositions[i], particlePositions[i],
				1e-6);
	}
}

void testImageMoleculesThreeWaterChloride(bool runTestWithAtomImaging,
		bool addPositionOffset) {

	double boxDimension = 10.;
	RealVec box(boxDimension, boxDimension, boxDimension);

	unsigned int numberOfParticles = 10;
	std::vector<RealVec> particlePositions(numberOfParticles);

	particlePositions[0] = RealVec(0., 1., 0.);
	particlePositions[1] = RealVec(0., 0., 0.);
	particlePositions[2] = RealVec(0., 2., 0.);

	particlePositions[3] = RealVec(0, -4.5, 0.);
	particlePositions[4] = RealVec(-1, -4.5, 0.);
	particlePositions[5] = RealVec(1, -4.5, 0.);

	particlePositions[6] = RealVec(-2, 4., 0.);
	particlePositions[7] = RealVec(-2, 3., 0.);
	particlePositions[8] = RealVec(-2, -5., 0.);

	particlePositions[9] = RealVec(3, -4.5, 0.);

	std::vector<RealVec> expectedParticlePositions(numberOfParticles);

	expectedParticlePositions[0] = RealVec(0., 1., 0.);
	expectedParticlePositions[1] = RealVec(0., 0., 0.);
	expectedParticlePositions[2] = RealVec(0., 2., 0.);

	expectedParticlePositions[3] = RealVec(0, 5.5, 0.);
	expectedParticlePositions[4] = RealVec(-1, 5.5, 0.);
	expectedParticlePositions[5] = RealVec(1, 5.5, 0.);

	expectedParticlePositions[6] = RealVec(-2, 4., 0.);
	expectedParticlePositions[7] = RealVec(-2, 3., 0.);
	expectedParticlePositions[8] = RealVec(-2, 5., 0.);

	expectedParticlePositions[9] = RealVec(3, 5.5, 0.);

	if (addPositionOffset) {
		// move second molecule 1 box dimension in Y direction
		particlePositions[3][0] += boxDimension;
		particlePositions[4][0] += boxDimension;
		particlePositions[5][0] += boxDimension;
	}

	double imagedPositions[numberOfParticles * 2];
	imageMolecules(box, particlePositions);

	if (runTestWithAtomImaging) {
		// Check that making periodic images of everything with respect to the first oxygen fails
		imageParticles(box, particlePositions[0], particlePositions[4]);
		imageParticles(box, particlePositions[0], particlePositions[5]);
	}

	RealVec tempPosition;
	for (int i = 0; i < numberOfParticles; i++) {
		std::cout << "Position atom " << i << ": "
				<< expectedParticlePositions[i] << " A" << std::endl;
		std::cout << "Position atom " << i << ": " << particlePositions[i]
				<< " A" << std::endl;
	}

	for (int i = 0; i < numberOfParticles; i++) {
		ASSERT_EQUAL_VEC(expectedParticlePositions[i], particlePositions[i],
				1e-6);
	}
}

int main(int numberOfArguments, char* argv[]) {

	try {
		std::cout << "TestReferenceMBPolTwoBodyForce running test..."
				<< std::endl;

		double boxDimension = 0;

		std::cout << "TestReferenceMBPolTwoBodyForce Cluster" << std::endl;
		testTwoBody(boxDimension, false);

		bool runTestWithAtomImaging = false;
		std::cout
				<< "TestReferenceMBPolTwoBodyForce  Image Molecules with two water"
				<< std::endl;
		testImageMolecules(runTestWithAtomImaging, false);
		// shift molecule of 1 boxDimension
		testImageMolecules(runTestWithAtomImaging, true);

		std::cout
				<< "TestReferenceMBPolTwoBodyForce  Periodic boundary conditions"
				<< std::endl;
		boxDimension = 50;
		testTwoBody(boxDimension, false);

		std::cout
				<< "TestReferenceMBPolTwoBodyForce  Periodic boundary conditions with boxDimension offset on second water molecule"
				<< std::endl;
		boxDimension = 50;
		testTwoBody(boxDimension, true);

		std::cout << "TestReferenceMBPolTwoBodyForce  With chloride"
				<< std::endl;
		boxDimension = 0;
		testTwoBodyChloride(boxDimension, false);

		std::cout << "TestReferenceMBPolTwoBodyForce  With two water and a chloride"
				<< std::endl;
		boxDimension = 0;
		testTwoWaterOneChloride(boxDimension, false);

		std::cout
				<< "TestReferenceMBPolTwoBodyForce  Image Molecules with one water and chloride"
				<< std::endl;
		testImageMoleculesOneWaterChloride(runTestWithAtomImaging, false);
		// shift molecule of 1 boxDimension
		testImageMoleculesOneWaterChloride(runTestWithAtomImaging, true);

		std::cout
				<< "TestReferenceMBPolTwoBodyForce  Image Molecules with two water and chloride"
				<< std::endl;
		testImageMoleculesTwoWaterChloride(runTestWithAtomImaging, false);
		// shift molecule of 1 boxDimension
		testImageMoleculesTwoWaterChloride(runTestWithAtomImaging, true);

		std::cout
				<< "TestReferenceMBPolTwoBodyForce  Image Molecules with three water and chloride"
				<< std::endl;
		testImageMoleculesThreeWaterChloride(runTestWithAtomImaging, false);
		// shift molecule of 1 boxDimension
		testImageMoleculesThreeWaterChloride(runTestWithAtomImaging, true);

	} catch (const std::exception& e) {
		std::cout << "exception: " << e.what() << std::endl;
		std::cout << "FAIL - ERROR.  Test failed." << std::endl;
		return 1;
	}

	std::cout << "Done" << std::endl;
	return 0;
}
