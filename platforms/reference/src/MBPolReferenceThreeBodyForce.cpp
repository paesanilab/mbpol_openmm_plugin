/* Portions copyright (c) 2006 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
#include "MBPolReferenceThreeBodyForce.h"
#include "MBPolReferenceTwoBodyForce.h"
#include <algorithm>
#include <cctype>
#include "mbpol_3body_constants.h"
#include "poly-3b-v2x.h"
#include "poly-3b-h2o-cl-v2x.h"
#include <list>
#include <iostream>

using std::vector;
using OpenMM::RealVec;

MBPolReferenceThreeBodyForce::MBPolReferenceThreeBodyForce() :
		_nonbondedMethod(NoCutoff), _cutoff(1.0e+10) {

	_periodicBoxDimensions = RealVec(0.0, 0.0, 0.0);
}

MBPolReferenceThreeBodyForce::NonbondedMethod MBPolReferenceThreeBodyForce::getNonbondedMethod(
		void) const {
	return _nonbondedMethod;
}

void MBPolReferenceThreeBodyForce::setNonbondedMethod(
		MBPolReferenceThreeBodyForce::NonbondedMethod nonbondedMethod) {
	_nonbondedMethod = nonbondedMethod;
}

void MBPolReferenceThreeBodyForce::setCutoff(double cutoff) {
	_cutoff = cutoff;
}

double MBPolReferenceThreeBodyForce::getCutoff(void) const {
	return _cutoff;
}

void MBPolReferenceThreeBodyForce::setPeriodicBox(const RealVec& box) {
	_periodicBoxDimensions = box;
}

RealVec MBPolReferenceThreeBodyForce::getPeriodicBox(void) const {
	return _periodicBoxDimensions;
}

double var(const double& k, const double& r0, const OpenMM::RealVec& a1,
		const OpenMM::RealVec& a2) {
	const double dx[3] = { (a1[0] - a2[0]) * nm_to_A, (a1[1] - a2[1]) * nm_to_A,
			(a1[2] - a2[2]) * nm_to_A };

	const double dsq = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
	const double d = std::sqrt(dsq);

	return std::exp(-k * (d - r0));
}

void g_var(const double& g, const double& k, const double& r0,
		const OpenMM::RealVec& a1, const OpenMM::RealVec& a2,
		OpenMM::RealVec& g1, OpenMM::RealVec& g2) {
	const double dx[3] = { (a1[0] - a2[0]) * nm_to_A, (a1[1] - a2[1]) * nm_to_A,
			(a1[2] - a2[2]) * nm_to_A };

	const double dsq = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
	const double d = std::sqrt(dsq);

	double gg = -k * g * std::exp(-k * (d - r0)) / d;

	const double cal2joule = 4.184;

	gg *= cal2joule * 10. * -1;

	for (int i = 0; i < 3; ++i) {
		g1[i] += gg * dx[i];
		g2[i] -= gg * dx[i];
	}
}

double threebody_f_switch(const double& r, double& g) {
	if (r > r3f) {
		g = 0.0;
		return 0.0;
	} else if (r > r3i) {
		const double t1 = M_PI / (r3f - r3i);
		const double x = (r - r3i) * t1;
		g = -std::sin(x) * t1 / 2.0;
		return (1.0 + std::cos(x)) / 2.0;
	} else {
		g = 0.0;
		return 1.0;
	}
}
double threebody_f_switch_chloride(const double& r, double& g) {
	if (r > r3f_chloride) {
		g = 0.0;
		return 0.0;
	} else if (r > r3i_chloride) {
		const double t1 = M_PI / (r3f_chloride - r3i_chloride);
		const double x = (r - r3i_chloride) * t1;
		g = -std::sin(x) * t1 / 2.0;
		return (1.0 + std::cos(x)) / 2.0;
	} else {
		g = 0.0;
		return 1.0;
	}
}

RealOpenMM MBPolReferenceThreeBodyForce::calculateTripletIxn(int siteI,
		int siteJ, int siteQ, const std::vector<RealVec>& particlePositions,
		const std::vector<std::vector<int> >& allParticleIndices,
		vector<RealVec>& forces) const {

	// siteI and siteJ are indices in a oxygen-only array, in order to get the position of an oxygen, we need:
	// allParticleIndices[siteI][0]
	// first hydrogen: allParticleIndices[siteI][1]
	// second hydrogen: allParticleIndices[siteI][2]
	// same for the second water molecule
	const double cal2joule = 4.184;

	//if the sites are all waters
	if ((allParticleIndices[siteJ][0] == allParticleIndices[siteJ][1] - 1
			&& allParticleIndices[siteJ][0] == allParticleIndices[siteJ][2] - 2
			&& allParticleIndices[siteI][0] == allParticleIndices[siteI][1] - 1
			&& allParticleIndices[siteI][0] == allParticleIndices[siteI][2] - 2
			&& allParticleIndices[siteQ][0] == allParticleIndices[siteQ][1] - 1
			&& allParticleIndices[siteQ][0] == allParticleIndices[siteQ][2] - 2)) {
		std::vector<RealVec> allPositions;

		// the iterator constructor can also be used to construct from arrays:
		int sites_array[] = { siteI, siteJ, siteQ };
		std::list<int> sites(sites_array,
				sites_array + sizeof(sites_array) / sizeof(int));

		for (std::list<int>::iterator it = sites.begin(); it != sites.end();
				it++) {
			for (unsigned int i = 0; i < 3; i++)
				allPositions.push_back(
						particlePositions[allParticleIndices[*it][i]]);
		}

		if (_nonbondedMethod == CutoffPeriodic)
			imageMolecules(_periodicBoxDimensions, allPositions);

		double x[36];

		x[0] = var(kHH_intra, dHH_intra, allPositions[Ha1], allPositions[Ha2]);
		x[1] = var(kHH_intra, dHH_intra, allPositions[Hb1], allPositions[Hb2]);
		x[2] = var(kHH_intra, dHH_intra, allPositions[Hc1], allPositions[Hc2]);
		x[3] = var(kOH_intra, dOH_intra, allPositions[Oa], allPositions[Ha1]);
		x[4] = var(kOH_intra, dOH_intra, allPositions[Oa], allPositions[Ha2]);
		x[5] = var(kOH_intra, dOH_intra, allPositions[Ob], allPositions[Hb1]);
		x[6] = var(kOH_intra, dOH_intra, allPositions[Ob], allPositions[Hb2]);
		x[7] = var(kOH_intra, dOH_intra, allPositions[Oc], allPositions[Hc1]);
		x[8] = var(kOH_intra, dOH_intra, allPositions[Oc], allPositions[Hc2]);

		x[9] = var(kHH, dHH, allPositions[Ha1], allPositions[Hb1]);
		x[10] = var(kHH, dHH, allPositions[Ha1], allPositions[Hb2]);
		x[11] = var(kHH, dHH, allPositions[Ha1], allPositions[Hc1]);
		x[12] = var(kHH, dHH, allPositions[Ha1], allPositions[Hc2]);
		x[13] = var(kHH, dHH, allPositions[Ha2], allPositions[Hb1]);
		x[14] = var(kHH, dHH, allPositions[Ha2], allPositions[Hb2]);
		x[15] = var(kHH, dHH, allPositions[Ha2], allPositions[Hc1]);
		x[16] = var(kHH, dHH, allPositions[Ha2], allPositions[Hc2]);
		x[17] = var(kHH, dHH, allPositions[Hb1], allPositions[Hc1]);
		x[18] = var(kHH, dHH, allPositions[Hb1], allPositions[Hc2]);
		x[19] = var(kHH, dHH, allPositions[Hb2], allPositions[Hc1]);
		x[20] = var(kHH, dHH, allPositions[Hb2], allPositions[Hc2]);
		x[21] = var(kOH, dOH, allPositions[Oa], allPositions[Hb1]);
		x[22] = var(kOH, dOH, allPositions[Oa], allPositions[Hb2]);
		x[23] = var(kOH, dOH, allPositions[Oa], allPositions[Hc1]);
		x[24] = var(kOH, dOH, allPositions[Oa], allPositions[Hc2]);
		x[25] = var(kOH, dOH, allPositions[Ob], allPositions[Ha1]);
		x[26] = var(kOH, dOH, allPositions[Ob], allPositions[Ha2]);
		x[27] = var(kOH, dOH, allPositions[Ob], allPositions[Hc1]);
		x[28] = var(kOH, dOH, allPositions[Ob], allPositions[Hc2]);
		x[29] = var(kOH, dOH, allPositions[Oc], allPositions[Ha1]);
		x[30] = var(kOH, dOH, allPositions[Oc], allPositions[Ha2]);
		x[31] = var(kOH, dOH, allPositions[Oc], allPositions[Hb1]);
		x[32] = var(kOH, dOH, allPositions[Oc], allPositions[Hb2]);
		x[33] = var(kOO, dOO, allPositions[Oa], allPositions[Ob]);
		x[34] = var(kOO, dOO, allPositions[Oa], allPositions[Oc]);
		x[35] = var(kOO, dOO, allPositions[Ob], allPositions[Oc]);

		double g[36];
		double retval = poly_3b_v2x::eval(thefit, x, g);

		RealVec rab, rac, rbc;
		double drab(0), drac(0), drbc(0);

		rab = (allPositions[Oa] - allPositions[Ob]) * nm_to_A;
		drab += rab.dot(rab);

		rac = (allPositions[Oa] - allPositions[Oc]) * nm_to_A;
		drac += rac.dot(rac);

		rbc = (allPositions[Ob] - allPositions[Oc]) * nm_to_A;
		drbc += rbc.dot(rbc);

		drab = std::sqrt(drab);
		drac = std::sqrt(drac);
		drbc = std::sqrt(drbc);

		if ((drab < 2) or (drac < 2) or (drbc < 2))
			return 0.;
		double gab, gac, gbc;

		const double sab = threebody_f_switch(drab, gab);
		const double sac = threebody_f_switch(drac, gac);
		const double sbc = threebody_f_switch(drbc, gbc);

		const double s = sab * sac + sab * sbc + sac * sbc;

		for (int n = 0; n < 36; ++n)
			g[n] *= s;

		std::vector<RealVec> allForces;
		allForces.resize(allPositions.size());

		g_var(g[0], kHH_intra, dHH_intra, allPositions[Ha1], allPositions[Ha2],
				allForces[Ha1], allForces[Ha2]);
		g_var(g[1], kHH_intra, dHH_intra, allPositions[Hb1], allPositions[Hb2],
				allForces[Hb1], allForces[Hb2]);
		g_var(g[2], kHH_intra, dHH_intra, allPositions[Hc1], allPositions[Hc2],
				allForces[Hc1], allForces[Hc2]);
		g_var(g[3], kOH_intra, dOH_intra, allPositions[Oa], allPositions[Ha1],
				allForces[Oa], allForces[Ha1]);
		g_var(g[4], kOH_intra, dOH_intra, allPositions[Oa], allPositions[Ha2],
				allForces[Oa], allForces[Ha2]);
		g_var(g[5], kOH_intra, dOH_intra, allPositions[Ob], allPositions[Hb1],
				allForces[Ob], allForces[Hb1]);
		g_var(g[6], kOH_intra, dOH_intra, allPositions[Ob], allPositions[Hb2],
				allForces[Ob], allForces[Hb2]);
		g_var(g[7], kOH_intra, dOH_intra, allPositions[Oc], allPositions[Hc1],
				allForces[Oc], allForces[Hc1]);
		g_var(g[8], kOH_intra, dOH_intra, allPositions[Oc], allPositions[Hc2],
				allForces[Oc], allForces[Hc2]);
		g_var(g[9], kHH, dHH, allPositions[Ha1], allPositions[Hb1],
				allForces[Ha1], allForces[Hb1]);
		g_var(g[10], kHH, dHH, allPositions[Ha1], allPositions[Hb2],
				allForces[Ha1], allForces[Hb2]);
		g_var(g[11], kHH, dHH, allPositions[Ha1], allPositions[Hc1],
				allForces[Ha1], allForces[Hc1]);
		g_var(g[12], kHH, dHH, allPositions[Ha1], allPositions[Hc2],
				allForces[Ha1], allForces[Hc2]);
		g_var(g[13], kHH, dHH, allPositions[Ha2], allPositions[Hb1],
				allForces[Ha2], allForces[Hb1]);
		g_var(g[14], kHH, dHH, allPositions[Ha2], allPositions[Hb2],
				allForces[Ha2], allForces[Hb2]);
		g_var(g[15], kHH, dHH, allPositions[Ha2], allPositions[Hc1],
				allForces[Ha2], allForces[Hc1]);
		g_var(g[16], kHH, dHH, allPositions[Ha2], allPositions[Hc2],
				allForces[Ha2], allForces[Hc2]);
		g_var(g[17], kHH, dHH, allPositions[Hb1], allPositions[Hc1],
				allForces[Hb1], allForces[Hc1]);
		g_var(g[18], kHH, dHH, allPositions[Hb1], allPositions[Hc2],
				allForces[Hb1], allForces[Hc2]);
		g_var(g[19], kHH, dHH, allPositions[Hb2], allPositions[Hc1],
				allForces[Hb2], allForces[Hc1]);
		g_var(g[20], kHH, dHH, allPositions[Hb2], allPositions[Hc2],
				allForces[Hb2], allForces[Hc2]);
		g_var(g[21], kOH, dOH, allPositions[Oa], allPositions[Hb1],
				allForces[Oa], allForces[Hb1]);
		g_var(g[22], kOH, dOH, allPositions[Oa], allPositions[Hb2],
				allForces[Oa], allForces[Hb2]);
		g_var(g[23], kOH, dOH, allPositions[Oa], allPositions[Hc1],
				allForces[Oa], allForces[Hc1]);
		g_var(g[24], kOH, dOH, allPositions[Oa], allPositions[Hc2],
				allForces[Oa], allForces[Hc2]);
		g_var(g[25], kOH, dOH, allPositions[Ob], allPositions[Ha1],
				allForces[Ob], allForces[Ha1]);
		g_var(g[26], kOH, dOH, allPositions[Ob], allPositions[Ha2],
				allForces[Ob], allForces[Ha2]);
		g_var(g[27], kOH, dOH, allPositions[Ob], allPositions[Hc1],
				allForces[Ob], allForces[Hc1]);
		g_var(g[28], kOH, dOH, allPositions[Ob], allPositions[Hc2],
				allForces[Ob], allForces[Hc2]);
		g_var(g[29], kOH, dOH, allPositions[Oc], allPositions[Ha1],
				allForces[Oc], allForces[Ha1]);
		g_var(g[30], kOH, dOH, allPositions[Oc], allPositions[Ha2],
				allForces[Oc], allForces[Ha2]);
		g_var(g[31], kOH, dOH, allPositions[Oc], allPositions[Hb1],
				allForces[Oc], allForces[Hb1]);
		g_var(g[32], kOH, dOH, allPositions[Oc], allPositions[Hb2],
				allForces[Oc], allForces[Hb2]);
		g_var(g[33], kOO, dOO, allPositions[Oa], allPositions[Ob],
				allForces[Oa], allForces[Ob]);
		g_var(g[34], kOO, dOO, allPositions[Oa], allPositions[Oc],
				allForces[Oa], allForces[Oc]);
		g_var(g[35], kOO, dOO, allPositions[Ob], allPositions[Oc],
				allForces[Ob], allForces[Oc]);

		// gradients of the switching function

		gab *= (sac + sbc) * retval / drab;
		gac *= (sab + sbc) * retval / drac;
		gbc *= (sab + sac) * retval / drbc;

		retval *= s;

		for (int n = 0; n < 3; ++n) {
			allForces[Oa][n] += (gab * rab[n] + gac * rac[n]) * cal2joule
					* -nm_to_A;
			allForces[Ob][n] += (gbc * rbc[n] - gab * rab[n]) * cal2joule
					* -nm_to_A;
			allForces[Oc][n] -= (gac * rac[n] + gbc * rbc[n]) * cal2joule
					* -nm_to_A;
		}

		unsigned int j = 0;
		for (std::list<int>::iterator it = sites.begin(); it != sites.end();
				it++) {
			for (unsigned int i = 0; i < 3; i++) {
				forces[allParticleIndices[*it][i]] += allForces[j];
				j++;
			}
		}

		RealOpenMM energy = retval * cal2joule;

		return energy;
	} else {
		//THERE IS AN ION
		//determine which site is the ion
//		// Determine which site is the ion
//		if (!(allParticleIndices[siteJ][0] == allParticleIndices[siteJ][1] - 1
//				&& allParticleIndices[siteJ][0]
//						== allParticleIndices[siteJ][2] - 2)) {
//			// if siteJ is the ion, swap siteJ and siteI
//			int temp = siteJ;
//			siteJ = siteI;
//			siteI = temp;
//
//		} else if (!(allParticleIndices[siteQ][0]
//				== allParticleIndices[siteQ][1] - 1
//				&& allParticleIndices[siteQ][0]
//						== allParticleIndices[siteQ][2] - 2)) {
//			// if siteJ is the ion, swap siteQ and siteI
//			int temp = siteQ;
//			siteQ = siteI;
//			siteI = temp;
//
//		}
//		// from here on the ion will be in site I
		std::vector<RealVec> allPositions;

		// pushes back all positions of particles in waters
		for (unsigned int i = 0; i < 3; i++) {
			allPositions.push_back(
					particlePositions[allParticleIndices[siteQ][i]]);
		}
		for (unsigned int i = 0; i < 3; i++) {
			allPositions.push_back(
					particlePositions[allParticleIndices[siteJ][i]]);
		}
		// only push back the Cl position
		allPositions.push_back(particlePositions[allParticleIndices[siteI][0]]);

		if (_nonbondedMethod == CutoffPeriodic)
			imageMolecules(_periodicBoxDimensions, allPositions);

		double x[21];

		x[0] = var(kHH_intra_chloride, dHH_intra_chloride, allPositions[Ha1],
				allPositions[Ha2]);
		x[1] = var(kHH_intra_chloride, dHH_intra_chloride, allPositions[Hb1],
				allPositions[Hb2]);
		x[2] = var(kOH_intra_chloride, dOH_intra_chloride, allPositions[Oa],
				allPositions[Ha1]);
		x[3] = var(kOH_intra_chloride, dOH_intra_chloride, allPositions[Oa],
				allPositions[Ha2]);
		x[4] = var(kOH_intra_chloride, dOH_intra_chloride, allPositions[Ob],
				allPositions[Hb1]);
		x[5] = var(kOH_intra_chloride, dOH_intra_chloride, allPositions[Ob],
				allPositions[Hb2]);

		x[6] = var(kHH_chloride, dHH_chloride, allPositions[Ha1],
				allPositions[Hb1]);
		x[7] = var(kHH_chloride, dHH_chloride, allPositions[Ha1],
				allPositions[Hb2]);
		x[8] = var(kHH_chloride, dHH_chloride, allPositions[Ha2],
				allPositions[Hb1]);
		x[9] = var(kHH_chloride, dHH_chloride, allPositions[Ha2],
				allPositions[Hb2]);
		x[10] = var(kOH_chloride, dOH_chloride, allPositions[Oa],
				allPositions[Hb1]);
		x[11] = var(kOH_chloride, dOH_chloride, allPositions[Oa],
				allPositions[Hb2]);
		x[12] = var(kOH_chloride, dOH_chloride, allPositions[Ob],
				allPositions[Ha1]);
		x[13] = var(kOH_chloride, dOH_chloride, allPositions[Ob],
				allPositions[Ha2]);
		x[14] = var(kOO_chloride, dOO_chloride, allPositions[Oa],
				allPositions[Ob]);
		x[15] = var(kClH, dClH, allPositions[Cl2], allPositions[Ha1]);
		x[16] = var(kClH, dClH, allPositions[Cl2], allPositions[Ha2]);
		x[17] = var(kClH, dClH, allPositions[Cl2], allPositions[Hb1]);
		x[18] = var(kClH, dClH, allPositions[Cl2], allPositions[Hb2]);
		x[19] = var(kClO, dClO, allPositions[Cl2], allPositions[Oa]);
		x[20] = var(kClO, dClO, allPositions[Cl2], allPositions[Ob]);

		double g[21];
		double retval = h2o_cl::poly_3b_h2o_cl_v2x::eval(thefit_chloride, x, g);

		RealVec rab, rac, rbc;
		double drab(0), drac(0), drbc(0);

		rab = (allPositions[Oa] - allPositions[Ob]) * nm_to_A;
		drab += rab.dot(rab);

		rac = (allPositions[Oa] - allPositions[Cl2]) * nm_to_A;
		drac += rac.dot(rac);

		rbc = (allPositions[Ob] - allPositions[Cl2]) * nm_to_A;
		drbc += rbc.dot(rbc);

		drab = std::sqrt(drab);
		drac = std::sqrt(drac);
		drbc = std::sqrt(drbc);

		if ((drab < 2) or (drac < 2) or (drbc < 2))
			return 0.;
		double gab, gac, gbc;

		const double sab = threebody_f_switch_chloride(drab, gab);
		const double sac = threebody_f_switch_chloride(drac, gac);
		const double sbc = threebody_f_switch_chloride(drbc, gbc);

		const double s = sab * sac + sab * sbc + sac * sbc;
		for (int n = 0; n < 21; ++n)
			g[n] *= s;
		std::vector<RealVec> allForces;
		allForces.resize(allPositions.size());

		g_var(g[0], kHH_intra_chloride, dHH_intra_chloride, allPositions[Ha1],
				allPositions[Ha2], allForces[Ha1], allForces[Ha2]);
		g_var(g[1], kHH_intra_chloride, dHH_intra_chloride, allPositions[Hb1],
				allPositions[Hb2], allForces[Hb1], allForces[Hb2]);
		g_var(g[2], kOH_intra_chloride, dOH_intra_chloride, allPositions[Oa],
				allPositions[Ha1], allForces[Oa], allForces[Ha1]);
		g_var(g[3], kOH_intra_chloride, dOH_intra_chloride, allPositions[Oa],
				allPositions[Ha2], allForces[Oa], allForces[Ha2]);
		g_var(g[4], kOH_intra_chloride, dOH_intra_chloride, allPositions[Ob],
				allPositions[Hb1], allForces[Ob], allForces[Hb1]);
		g_var(g[5], kOH_intra_chloride, dOH_intra_chloride, allPositions[Ob],
				allPositions[Hb2], allForces[Ob], allForces[Hb2]);

		g_var(g[6], kHH_chloride, dHH_chloride, allPositions[Ha1],
				allPositions[Hb1], allForces[Ha1], allForces[Hb1]);
		g_var(g[7], kHH_chloride, dHH_chloride, allPositions[Ha1],
				allPositions[Hb2], allForces[Ha1], allForces[Hb2]);
		g_var(g[8], kHH_chloride, dHH_chloride, allPositions[Ha2],
				allPositions[Hb1], allForces[Ha2], allForces[Hb1]);
		g_var(g[9], kHH_chloride, dHH_chloride, allPositions[Ha2],
				allPositions[Hb2], allForces[Ha2], allForces[Hb2]);
		g_var(g[10], kOH_chloride, dOH_chloride, allPositions[Oa],
				allPositions[Hb1], allForces[Oa], allForces[Hb1]);
		g_var(g[11], kOH_chloride, dOH_chloride, allPositions[Oa],
				allPositions[Hb2], allForces[Oa], allForces[Hb2]);
		g_var(g[12], kOH_chloride, dOH_chloride, allPositions[Ob],
				allPositions[Ha1], allForces[Ob], allForces[Ha1]);
		g_var(g[13], kOH_chloride, dOH_chloride, allPositions[Ob],
				allPositions[Ha2], allForces[Ob], allForces[Ha2]);
		g_var(g[14], kOO_chloride, dOO_chloride, allPositions[Oa],
				allPositions[Ob], allForces[Oa], allForces[Ob]);
		g_var(g[15], kClH, dClH, allPositions[Cl2], allPositions[Ha1],
				allForces[Cl2], allForces[Ha1]);
		g_var(g[16], kClH, dClH, allPositions[Cl2], allPositions[Ha2],
				allForces[Cl2], allForces[Ha2]);
		g_var(g[17], kClH, dClH, allPositions[Cl2], allPositions[Hb1],
				allForces[Cl2], allForces[Hb1]);
		g_var(g[18], kClH, dClH, allPositions[Cl2], allPositions[Hb2],
				allForces[Cl2], allForces[Hb2]);
		g_var(g[19], kClO, dClO, allPositions[Cl2], allPositions[Oa],
				allForces[Cl2], allForces[Oa]);
		g_var(g[20], kClO, dClO, allPositions[Cl2], allPositions[Ob],
				allForces[Cl2], allForces[Ob]);

		// gradients of the switching function

		gab *= (sac + sbc) * retval / drab;
		gac *= (sab + sbc) * retval / drac;
		gbc *= (sab + sac) * retval / drbc;

		retval *= s;

		for (int n = 0; n < 3; ++n) {
			allForces[Oa][n] += (gab * rab[n] + gac * rac[n]) * cal2joule
					* -nm_to_A;
			allForces[Ob][n] += (gbc * rbc[n] - gab * rab[n]) * cal2joule
					* -nm_to_A;
			allForces[Cl2][n] -= (gac * rac[n] + gbc * rbc[n]) * cal2joule
					* -nm_to_A;
		}

//		unsigned int j = 0;
//		for (std::list<int>::iterator it = sites.begin(); it != sites.end();
//				it++) {
//			for (unsigned int i = 0; i < 3; i++) {
//				forces[allParticleIndices[*it][i]] += allForces[j];
//				j++;
//			}
//		}
		// water molecules
		forces[allParticleIndices[siteQ][0]] += allForces[Oa];
		forces[allParticleIndices[siteQ][1]] += allForces[Ha1];
		forces[allParticleIndices[siteQ][2]] += allForces[Ha2];
		forces[allParticleIndices[siteJ][0]] += allForces[Ob];
		forces[allParticleIndices[siteJ][1]] += allForces[Hb1];
		forces[allParticleIndices[siteJ][2]] += allForces[Hb2];
		// Cl
		forces[allParticleIndices[siteI][0]] += allForces[Cl2];

		RealOpenMM energy = retval * cal2joule;

		return energy;

	}

}

RealOpenMM MBPolReferenceThreeBodyForce::calculateForceAndEnergy(
		int numParticles, const vector<RealVec>& particlePositions,
		const std::vector<std::vector<int> >& allParticleIndices,
		const ThreeNeighborList& neighborList, vector<RealVec>& forces) const {

	// loop over neighbor list
	//    (1) calculate pair vdw ixn
	//    (2) accumulate forces: if particle is a site where interaction position != particle position,
	//        then call addReducedForce() to apportion force to particle and its covalent partner
	//        based on reduction factor

	RealOpenMM energy = 0.;
	for (unsigned int ii = 0; ii < neighborList.size(); ii++) {

		MBPolPlugin::AtomTriplet triplet = neighborList[ii];
		int siteI = triplet.first;
		int siteJ = triplet.second;
		int siteQ = triplet.third;

		energy += calculateTripletIxn(siteI, siteJ, siteQ, particlePositions,
				allParticleIndices, forces);

	}

	return energy;
}
