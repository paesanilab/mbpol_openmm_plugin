/*
 * TestCudaMBPolElectrostaticsForce.cpp
 *
 *  Created on: Sep 30, 2015
 *      Author: sebastian
 *
 * This tests the CUDA implementation of ReferenceMBPolElectrostaticsForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMMBPol.h"
#include "openmm/System.h"
#include "openmm/MBPolElectrostaticsForce.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/VirtualSite.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#define ASSERT_EQUAL_TOL_MOD(expected, found, tol, testname) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC_MOD(expected, found, tol,testname) {ASSERT_EQUAL_TOL_MOD((expected)[0], (found)[0], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[1], (found)[1], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[2], (found)[2], (tol),(testname));};

// #define COMPUTE_FINITE_DIFFERENCES_FORCES

using namespace OpenMM;
using namespace MBPolPlugin;
using namespace std;
const double TOL = 1e-4;
const double cal2joule = 4.184;

extern "C" OPENMM_EXPORT void registerMBPolCudaKernelFactories();

static void testWater3VirtualSite() {

	std::string testName = "testWater3VirtualSite";
	std::cout << "Test START: " << testName << std::endl;

	int numberOfParticles = 12;
	double cutoff = 0.70;

	std::vector<double> outputElectrostaticsMoments;
	std::vector<Vec3> inputGrid;
	std::vector<double> outputGridPotential;

	// beginning of Electrostatics setup
	MBPolElectrostaticsForce::NonbondedMethod nonbondedMethod =
			MBPolElectrostaticsForce::NoCutoff;

	System system;
	// box dimensions

	// double boxDimension                               = 1.8643;
	// Vec3 a( boxDimension, 0.0, 0.0 );
	// Vec3 b( 0.0, boxDimension, 0.0 );
	// Vec3 c( 0.0, 0.0, boxDimension );
	// system.setDefaultPeriodicBoxVectors( a, b, c );

	MBPolElectrostaticsForce* mbpolElectrostaticsForce =
			new MBPolElectrostaticsForce();
	;
	mbpolElectrostaticsForce->setNonbondedMethod(nonbondedMethod);
	//mbpolElectrostaticsForce->setPolarizationType( polarizationType );
	//mbpolElectrostaticsForce->setCutoffDistance( cutoff );
	//mbpolElectrostaticsForce->setMutualInducedTargetEpsilon( 1.0e-06 );
	//mbpolElectrostaticsForce->setMutualInducedMaxIterations( 500 );
	//mbpolElectrostaticsForce->setAEwald( 5.4459052e+00 );
	//mbpolElectrostaticsForce->setEwaldErrorTolerance( 1.0e-04 );

	double virtualSiteWeightO = 0.573293118;
	double virtualSiteWeightH = 0.213353441;
	for (unsigned int jj = 0; jj < numberOfParticles; jj += 4) {
		system.addParticle(1.5999000e+01);
		system.addParticle(1.0080000e+00);
		system.addParticle(1.0080000e+00);
		system.addParticle(0.); // Virtual Site
		system.setVirtualSite(jj + 3,
				new ThreeParticleAverageSite(jj, jj + 1, jj + 2,
						virtualSiteWeightO, virtualSiteWeightH,
						virtualSiteWeightH));

	}

	std::vector<double> zeroDipole(3);
	std::vector<double> zeroQuadrupole(9);
	std::vector<double> thole(5);

	std::fill(zeroDipole.begin(), zeroDipole.end(), 0.);
	std::fill(zeroQuadrupole.begin(), zeroQuadrupole.end(), 0.);

	thole[TCC] = 0.4;
	thole[TCD] = 0.4;
	thole[TDD] = 0.055;
	thole[TDDOH] = 0.626;
	thole[TDDHH] = 0.055;

	for (unsigned int jj = 0; jj < numberOfParticles; jj += 4) {
		mbpolElectrostaticsForce->addElectrostatics(-5.1966000e-01, jj + 1,
				jj + 2, jj + 3, thole, 0.001310, 0.001310);
		mbpolElectrostaticsForce->addElectrostatics(2.5983000e-01, jj, jj + 2,
				jj + 3, thole, 0.000294, 0.000294);
		mbpolElectrostaticsForce->addElectrostatics(2.5983000e-01, jj, jj + 1,
				jj + 3, thole, 0.000294, 0.000294);
		mbpolElectrostaticsForce->addElectrostatics(0., jj, jj + 1, jj + 2,
				thole, 0.001310, 0.);
	}

	system.addForce(mbpolElectrostaticsForce);

	static std::vector<Vec3> positions; // Static to work around bug in Visual Studio that makes compilation very very slow.
	positions.resize(numberOfParticles);

	positions[0] = Vec3(-1.516074336e+00, -2.023167650e-01, 1.454672917e+00);
	positions[1] = Vec3(-6.218989773e-01, -6.009430735e-01, 1.572437625e+00);
	positions[2] = Vec3(-2.017613812e+00, -4.190350349e-01, 2.239642849e+00);
	positions[3] = Vec3(-1.43230412, -0.33360265, 1.64727446);

	positions[4] = Vec3(-1.763651687e+00, -3.816594649e-01, -1.300353949e+00);
	positions[5] = Vec3(-1.903851736e+00, -4.935677617e-01, -3.457810126e-01);
	positions[6] = Vec3(-2.527904158e+00, -7.613550077e-01, -1.733803676e+00);
	positions[7] = Vec3(-1.95661974, -0.48654484, -1.18917052);

	positions[8] = Vec3(-5.588472140e-01, 2.006699172e+00, -1.392786582e-01);
	positions[9] = Vec3(-9.411558180e-01, 1.541226676e+00, 6.163293071e-01);
	positions[10] = Vec3(-9.858551734e-01, 1.567124294e+00, -8.830970941e-01);
	positions[11] = Vec3(-0.73151769, 1.8136042, -0.13676332);

	for (int i = 0; i < numberOfParticles; i++) {
		for (int j = 0; j < 3; j++) {
			positions[i][j] *= 1e-1;
		}
	}

	std::string platformName;
	platformName = "CUDA";
	LangevinIntegrator integrator(0.0, 0.1, 0.01);
	// constructing the context is causing an exception
	Context context(system, integrator,
			Platform::getPlatformByName(platformName));

	context.setPositions(positions);
	context.applyConstraints(1e-4); // update position of virtual site
	double tolerance = 1.0e-04;

//    // test energy and forces

	State state = context.getState(State::Forces | State::Energy);
	std::vector<Vec3> forces = state.getForces();
	double energy = state.getPotentialEnergy();
	double cal2joule = 4.184;

	double expectedEnergy = -15.818784 * cal2joule;
	std::cout << "Energy: " << energy / cal2joule << " Kcal/mol " << std::endl;
	std::cout << "Expected energy: " << expectedEnergy / cal2joule
			<< " Kcal/mol " << std::endl;

	std::vector<Vec3> expectedForces(numberOfParticles);
	expectedForces[0] = Vec3(2.38799956, 0.126835228, 8.86189407);
	expectedForces[1] = Vec3(-4.21263312, -0.72316292, 3.37076777);
	expectedForces[2] = Vec3(2.19240288, -2.24806806, 1.96210789);
	expectedForces[4] = Vec3(3.59486021, -2.16710895, 3.57138432);
	expectedForces[5] = Vec3(-4.54547068, -4.58639226, -17.4258666);
	expectedForces[6] = Vec3(-3.27239433, -1.96722979, 1.1170853);
	expectedForces[8] = Vec3(-1.44387205, -3.22471108, -2.61329967);
	expectedForces[9] = Vec3(3.35011312, 6.07136704, -0.197008793);
	expectedForces[10] = Vec3(1.94899441, 8.7184708, 1.35293571);

	// gradient -> forces
	for (int i = 0; i < numberOfParticles; i++) {
		for (int j = 0; j < 3; j++) {
			forces[i][j] /= cal2joule * 10;
			if ((i + 1) % 4 == 0) { // Set virtual site force to 0
				forces[i][j] = 0;
			}
		}
	}

	for (int i = 0; i < numberOfParticles; i++) {
		for (int j = 0; j < 3; j++) {
			expectedForces[i][j] *= -1;
		}
	}
	std::cout << std::endl << "Forces:" << std::endl;

	const double eps = 1.0e-4;

	double x_orig;

	std::vector<Vec3> finiteDifferenceForces(numberOfParticles);
//	for (int i = 0; i < numberOfParticles; i++) {
//		finiteDifferenceForces.push_back(Vec3(0., 0., 0.));
//	}
//	for (int i = 0; i < numberOfParticles; i++) {
//		for (int xyz = 0; xyz < 3; xyz++) {
//			x_orig = positions[i][xyz];
//
//			positions[i][xyz] = x_orig + eps;
//			context.setPositions(positions);
//			context.applyConstraints(1e-4); // update position of virtual site
//			state = context.getState(State::Energy);
//			const double Ep = state.getPotentialEnergy();
//
//			positions[i][xyz] = x_orig + 2 * eps;
//			context.setPositions(positions);
//			context.applyConstraints(1e-4); // update position of virtual site
//			state = context.getState(State::Energy);
//			const double E2p = state.getPotentialEnergy();
//
//			positions[i][xyz] = x_orig - eps;
//			context.setPositions(positions);
//			context.applyConstraints(1e-4); // update position of virtual site
//			state = context.getState(State::Energy);
//			const double Em = state.getPotentialEnergy();
//
//			positions[i][xyz] = x_orig - 2 * eps;
//			context.setPositions(positions);
//			context.applyConstraints(1e-4); // update position of virtual site
//			state = context.getState(State::Energy);
//			const double E2m = state.getPotentialEnergy();
//
//			finiteDifferenceForces[i][xyz] = (8 * (Ep - Em) - (E2p - E2m))
//					/ (12 * eps);
//			positions[i][xyz] = x_orig;
//		}
//
//	}

	for (int i = 0; i < numberOfParticles; i++) {
		for (int j = 0; j < 3; j++) {
			finiteDifferenceForces[i][j] /= -1 * cal2joule * 10;
		}

	}

	for (int i = 0; i < numberOfParticles; i++) {
		std::cout << "Force atom " << i << ": " << expectedForces[i]
				<< " Kcal/mol/A <mbpol>" << std::endl;
		std::cout << "Force atom " << i << ": " << forces[i]
				<< " Kcal/mol/A <openmm-mbpol>" << std::endl;
		std::cout << "Force atom " << i << ": " << finiteDifferenceForces[i]
				<< " Kcal/mol/A <openmm-mbpol finite differences>" << std::endl
				<< std::endl;
	}

	std::cout << "Comparison of energy and forces with tolerance: " << tolerance
			<< std::endl << std::endl;

	ASSERT_EQUAL_TOL_MOD(expectedEnergy, energy, tolerance, testName);

	for (unsigned int ii = 0; ii < forces.size(); ii++) {
		ASSERT_EQUAL_VEC_MOD(expectedForces[ii], forces[ii], tolerance,
				testName);
	}

	std::cout << "Test Successful: " << testName << std::endl << std::endl;

	return;
}
int main(int numberOfArguments, char* argv[]) {
	registerMBPolCudaKernelFactories();
	try {
		std::cout << "TestReferenceMBPolElectrostaticsForce running test..."
				<< std::endl;
		testWater3VirtualSite();
	} catch (const std::exception& e) {
		std::cout << "exception: " << e.what() << std::endl;
		std::cout << "FAIL - ERROR.  Test failed." << std::endl;
		return 1;
	}
	std::cout << "Done" << std::endl;
	return 0;
}
