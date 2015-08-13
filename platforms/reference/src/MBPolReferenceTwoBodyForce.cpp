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

#include <iostream>

#include "MBPolReferenceTwoBodyForce.h"
#include <algorithm>
#include <cctype>
#include "mbpol_2body_constants.h"
#include "poly-2b-v6x.h"
#include "poly-2b-h2o-cl-v2x.h"
#include "openmm/internal/MBPolConstants.h"

using std::vector;
using OpenMM::RealVec;
using namespace MBPolPlugin;

MBPolReferenceTwoBodyForce::MBPolReferenceTwoBodyForce() :
		_nonbondedMethod(NoCutoff), _cutoff(1.0e+10) {

	_periodicBoxDimensions = RealVec(0.0, 0.0, 0.0);
}

MBPolReferenceTwoBodyForce::NonbondedMethod MBPolReferenceTwoBodyForce::getNonbondedMethod(
		void) const {
	return _nonbondedMethod;
}

void MBPolReferenceTwoBodyForce::setNonbondedMethod(
		MBPolReferenceTwoBodyForce::NonbondedMethod nonbondedMethod) {
	_nonbondedMethod = nonbondedMethod;
}

void MBPolReferenceTwoBodyForce::setCutoff(double cutoff) {
	_cutoff = cutoff;
}

double MBPolReferenceTwoBodyForce::getCutoff(void) const {
	return _cutoff;
}

void MBPolReferenceTwoBodyForce::setPeriodicBox(const RealVec& box) {
	_periodicBoxDimensions = box;
}

RealVec MBPolReferenceTwoBodyForce::getPeriodicBox(void) const {
	return _periodicBoxDimensions;
}

void imageParticles(const RealVec& box, const RealVec & referenceParticle,
		RealVec& particleToImage) {
	/* std::cout << "Got to imageParticles" << std::endl; */
	// Periodic boundary conditions imaging of particleToImage with respect to referenceParticle
	double distance, factor;
	for (unsigned int i = 0; i < 3; i++) {
		distance = referenceParticle[i] - particleToImage[i];
		factor = std::floor(distance / box[i] + 0.5);
		particleToImage[i] += box[i] * factor;
	}
}

void imageMolecules(const RealVec& box, std::vector<RealVec>& allPositions) {
	/* std::cout << "Got to imageMolecules" << std::endl; */

	// Take first oxygen as central atom
	// image its two hydrogens with respect of the first oxygen
	imageParticles(box, allPositions[Oa], allPositions[Ha1]);
	imageParticles(box, allPositions[Oa], allPositions[Ha2]);

	int lastIndex = allPositions.size() - 1;
	if (lastIndex >= Ob) // at least Two molecules
			{
		// Now image the first particle of the second molecule
		// (the Cl or the oxygen of the first second molecule)

		imageParticles(box, allPositions[Oa], allPositions[Ob]);

		// if the second particle is not an ion
		if (lastIndex > Cl) {
			// Image the hydrogen of the second molecule with respect to the oxygen of the second molecule
			imageParticles(box, allPositions[Ob], allPositions[Hb1]);
			imageParticles(box, allPositions[Ob], allPositions[Hb2]);
		}
		if (lastIndex >= Oc) // at least Three molecules
				{
			// Now image the first particle of the third molecule
			// (the Cl or the oxygen of the first third molecule)

			imageParticles(box, allPositions[Oa], allPositions[Oc]);

			// if the third particle is not an ion
			if (lastIndex > Cl2) {
				// Image the hydrogen of the third molecule with respect to the oxygen of the third molecule
				imageParticles(box, allPositions[Oc], allPositions[Hc1]);
				imageParticles(box, allPositions[Oc], allPositions[Hc2]);
			}
			// if there is 3 particles and an ion
			if (lastIndex == Cl3) {
				// Now image the ion
				imageParticles(box, allPositions[Oa], allPositions[Cl3]);
			}

		}
	}
}

RealOpenMM MBPolReferenceTwoBodyForce::calculatePairIxn(int siteI, int siteJ,
		const std::vector<RealVec>& particlePositions,
		const std::vector<std::vector<int> >& allParticleIndices,
		vector<RealVec>& forces) const {
	/* std::cout << "Got to calculatePairIxn" << std::endl; */

	/* std::cout << siteJ << ": ";
	 for (unsigned int i = 0; i < 3; i++)
	 /* std::cout << allParticleIndices[siteJ][i] << " ";
	 /* std::cout << std::endl; */
	/* std::cout << siteI << ": ";
	 for (unsigned int i = 0; i < 3; i++)
	 /* std::cout << allParticleIndices[siteI][i] << " ";
	 /* std::cout << std::endl; */

	/* std::cout << "is check true? :"
	 << (allParticleIndices[siteJ][0] == allParticleIndices[siteJ][1] - 1
	 && allParticleIndices[siteJ][0]
	 == allParticleIndices[siteJ][2] - 2
	 && allParticleIndices[siteI][0]
	 == allParticleIndices[siteI][1] - 1
	 && allParticleIndices[siteI][0]
	 == allParticleIndices[siteI][2] - 2) << std::endl; */

	// if statement to determine if W-W interaction, based on the fact that waters will have
	// multiple particles that are increasing in index
	if ((allParticleIndices[siteJ][0] == allParticleIndices[siteJ][1] - 1
			&& allParticleIndices[siteJ][0] == allParticleIndices[siteJ][2] - 2
			&& allParticleIndices[siteI][0] == allParticleIndices[siteI][1] - 1
			&& allParticleIndices[siteI][0] == allParticleIndices[siteI][2] - 2)) {
		// siteI and siteJ are indices in a oxygen-only array, in order to get the position of an oxygen, we need:
		// allParticleIndices[siteI][0]
		// first hydrogen: allParticleIndices[siteI][1]
		// second hydrogen: allParticleIndices[siteI][2]
		// same for the second water molecule
		// offsets
		std::vector<RealVec> allPositions;

		for (unsigned int i = 0; i < 3; i++)
			allPositions.push_back(
					particlePositions[allParticleIndices[siteI][i]] * nm_to_A);

		for (unsigned int i = 0; i < 3; i++)
			allPositions.push_back(
					particlePositions[allParticleIndices[siteJ][i]] * nm_to_A);

		std::vector<RealVec> extraPoints;
		extraPoints.resize(4);

		if (_nonbondedMethod == CutoffPeriodic)
			imageMolecules(_periodicBoxDimensions * nm_to_A, allPositions);

		RealVec dOO = allPositions[Oa] - allPositions[Ob];

		const double rOOsq = dOO[0] * dOO[0] + dOO[1] * dOO[1]
				+ dOO[2] * dOO[2];
		const double rOO = std::sqrt(rOOsq);

		if (rOO > r2f)
			return 0.0;

		if (rOO < 2.)
			return 0.0;

		// the extra-points

		monomer ma, mb;

		ma.setup(allPositions[Oa], allPositions[Ha1], allPositions[Ha2],
				in_plane_gamma, out_of_plane_gamma, extraPoints[Xa1],
				extraPoints[Xa2]);

		mb.setup(allPositions[Ob], allPositions[Hb1], allPositions[Hb2],
				in_plane_gamma, out_of_plane_gamma, extraPoints[Xb1],
				extraPoints[Xb2]);

		// variables

		double v[31]; // stored separately (gets passed to poly::eval)

		variable ctxt[31];

		v[0] = ctxt[0].v_exp(d0_intra, k_HH_intra, allPositions[Ha1],
				allPositions[Ha2]);
		v[1] = ctxt[1].v_exp(d0_intra, k_HH_intra, allPositions[Hb1],
				allPositions[Hb2]);

		v[2] = ctxt[2].v_exp(d0_intra, k_OH_intra, allPositions[Oa],
				allPositions[Ha1]);
		v[3] = ctxt[3].v_exp(d0_intra, k_OH_intra, allPositions[Oa],
				allPositions[Ha2]);
		v[4] = ctxt[4].v_exp(d0_intra, k_OH_intra, allPositions[Ob],
				allPositions[Hb1]);
		v[5] = ctxt[5].v_exp(d0_intra, k_OH_intra, allPositions[Ob],
				allPositions[Hb2]);

		v[6] = ctxt[6].v_coul(d0_inter, k_HH_coul, allPositions[Ha1],
				allPositions[Hb1]);
		v[7] = ctxt[7].v_coul(d0_inter, k_HH_coul, allPositions[Ha1],
				allPositions[Hb2]);
		v[8] = ctxt[8].v_coul(d0_inter, k_HH_coul, allPositions[Ha2],
				allPositions[Hb1]);
		v[9] = ctxt[9].v_coul(d0_inter, k_HH_coul, allPositions[Ha2],
				allPositions[Hb2]);

		v[10] = ctxt[10].v_coul(d0_inter, k_OH_coul, allPositions[Oa],
				allPositions[Hb1]);
		v[11] = ctxt[11].v_coul(d0_inter, k_OH_coul, allPositions[Oa],
				allPositions[Hb2]);
		v[12] = ctxt[12].v_coul(d0_inter, k_OH_coul, allPositions[Ob],
				allPositions[Ha1]);
		v[13] = ctxt[13].v_coul(d0_inter, k_OH_coul, allPositions[Ob],
				allPositions[Ha2]);

		v[14] = ctxt[14].v_coul(d0_inter, k_OO_coul, allPositions[Oa],
				allPositions[Ob]);

		v[15] = ctxt[15].v_exp(d0_inter, k_XH_main, extraPoints[Xa1],
				allPositions[Hb1]);
		v[16] = ctxt[16].v_exp(d0_inter, k_XH_main, extraPoints[Xa1],
				allPositions[Hb2]);
		v[17] = ctxt[17].v_exp(d0_inter, k_XH_main, extraPoints[Xa2],
				allPositions[Hb1]);
		v[18] = ctxt[18].v_exp(d0_inter, k_XH_main, extraPoints[Xa2],
				allPositions[Hb2]);
		v[19] = ctxt[19].v_exp(d0_inter, k_XH_main, extraPoints[Xb1],
				allPositions[Ha1]);
		v[20] = ctxt[20].v_exp(d0_inter, k_XH_main, extraPoints[Xb1],
				allPositions[Ha2]);
		v[21] = ctxt[21].v_exp(d0_inter, k_XH_main, extraPoints[Xb2],
				allPositions[Ha1]);
		v[22] = ctxt[22].v_exp(d0_inter, k_XH_main, extraPoints[Xb2],
				allPositions[Ha2]);

		v[23] = ctxt[23].v_exp(d0_inter, k_XO_main, allPositions[Oa],
				extraPoints[Xb1]);
		v[24] = ctxt[24].v_exp(d0_inter, k_XO_main, allPositions[Oa],
				extraPoints[Xb2]);
		v[25] = ctxt[25].v_exp(d0_inter, k_XO_main, allPositions[Ob],
				extraPoints[Xa1]);
		v[26] = ctxt[26].v_exp(d0_inter, k_XO_main, allPositions[Ob],
				extraPoints[Xa2]);

		v[27] = ctxt[27].v_exp(d0_inter, k_XX_main, extraPoints[Xa1],
				extraPoints[Xb1]);
		v[28] = ctxt[28].v_exp(d0_inter, k_XX_main, extraPoints[Xa1],
				extraPoints[Xb2]);
		v[29] = ctxt[29].v_exp(d0_inter, k_XX_main, extraPoints[Xa2],
				extraPoints[Xb1]);
		v[30] = ctxt[30].v_exp(d0_inter, k_XX_main, extraPoints[Xa2],
				extraPoints[Xb2]);

		double g[31];
		const double E_poly = poly_2b_v6x_eval(thefit, v, g);

		std::vector<RealVec> allForces;
		allForces.resize(allPositions.size());
		/* std::cout << "there are " << allForces.size() << " forces" << std::endl; */

		std::vector<RealVec> extraForces;
		extraForces.resize(extraPoints.size());
		/* std::cout << "there are " << extraForces.size() << " extra forces" << std::endl; */

		ctxt[0].grads(g[0], allForces[Ha1], allForces[Ha2]);
		ctxt[1].grads(g[1], allForces[Hb1], allForces[Hb2]);

		ctxt[2].grads(g[2], allForces[Oa], allForces[Ha1]);
		ctxt[3].grads(g[3], allForces[Oa], allForces[Ha2]);
		ctxt[4].grads(g[4], allForces[Ob], allForces[Hb1]);
		ctxt[5].grads(g[5], allForces[Ob], allForces[Hb2]);

		ctxt[6].grads(g[6], allForces[Ha1], allForces[Hb1]);
		ctxt[7].grads(g[7], allForces[Ha1], allForces[Hb2]);
		ctxt[8].grads(g[8], allForces[Ha2], allForces[Hb1]);
		ctxt[9].grads(g[9], allForces[Ha2], allForces[Hb2]);

		ctxt[10].grads(g[10], allForces[Oa], allForces[Hb1]);
		ctxt[11].grads(g[11], allForces[Oa], allForces[Hb2]);
		ctxt[12].grads(g[12], allForces[Ob], allForces[Ha1]);
		ctxt[13].grads(g[13], allForces[Ob], allForces[Ha2]);

		ctxt[14].grads(g[14], allForces[Oa], allForces[Ob]);

		ctxt[15].grads(g[15], extraForces[Xa1], allForces[Hb1]);
		ctxt[16].grads(g[16], extraForces[Xa1], allForces[Hb2]);
		ctxt[17].grads(g[17], extraForces[Xa2], allForces[Hb1]);
		ctxt[18].grads(g[18], extraForces[Xa2], allForces[Hb2]);
		ctxt[19].grads(g[19], extraForces[Xb1], allForces[Ha1]);
		ctxt[20].grads(g[20], extraForces[Xb1], allForces[Ha2]);
		ctxt[21].grads(g[21], extraForces[Xb2], allForces[Ha1]);
		ctxt[22].grads(g[22], extraForces[Xb2], allForces[Ha2]);

		ctxt[23].grads(g[23], allForces[Oa], extraForces[Xb1]);
		ctxt[24].grads(g[24], allForces[Oa], extraForces[Xb2]);
		ctxt[25].grads(g[25], allForces[Ob], extraForces[Xa1]);
		ctxt[26].grads(g[26], allForces[Ob], extraForces[Xa2]);

		ctxt[27].grads(g[27], extraForces[Xa1], extraForces[Xb1]);
		ctxt[28].grads(g[28], extraForces[Xa1], extraForces[Xb2]);
		ctxt[29].grads(g[29], extraForces[Xa2], extraForces[Xb1]);
		ctxt[30].grads(g[30], extraForces[Xa2], extraForces[Xb2]);

		// distribute gradients w.r.t. the X-points

		ma.grads(extraForces[Xa1], extraForces[Xa2], in_plane_gamma,
				out_of_plane_gamma, allForces[Oa], allForces[Ha1],
				allForces[Ha2]);

		mb.grads(extraForces[Xb1], extraForces[Xb2], in_plane_gamma,
				out_of_plane_gamma, allForces[Ob], allForces[Hb1],
				allForces[Hb2]);

		// the switch

		double gsw;
		double sw = f_switch(rOO, gsw);

		double cal2joule = 4.184;

		// first water molecule
		forces[allParticleIndices[siteI][0]] += allForces[Oa] * sw * cal2joule
				* -10.;
		forces[allParticleIndices[siteI][1]] += allForces[Ha1] * sw * cal2joule
				* -10.;
		forces[allParticleIndices[siteI][2]] += allForces[Ha2] * sw * cal2joule
				* -10.;
		// second water molecule
		forces[allParticleIndices[siteJ][0]] += allForces[Ob] * sw * cal2joule
				* -10.;
		forces[allParticleIndices[siteJ][1]] += allForces[Hb1] * sw * cal2joule
				* -10.;
		forces[allParticleIndices[siteJ][2]] += allForces[Hb2] * sw * cal2joule
				* -10.;

		// gradient of the switch
		gsw *= E_poly / rOO;
		for (int i = 0; i < 3; ++i) {
			const double d = gsw * dOO[i];
			forces[allParticleIndices[siteI][0]][i] += d * cal2joule * -10.;
			forces[allParticleIndices[siteJ][0]][i] -= d * cal2joule * -10.;
		}

		RealOpenMM energy = sw * E_poly * cal2joule;
		return energy;

	} else { // if a w-Cl interaction
		// for the following code, an ion is in siteI so if the ion comes in as siteJ switch the sites
		if (allParticleIndices[siteI][0] == allParticleIndices[siteI][1] - 1
				&& allParticleIndices[siteI][0]
						== allParticleIndices[siteI][2] - 2) {
			// if siteI is the water and it is known that there is an ion, switch the sites (not tested)
			int temp = siteI;
			siteI = siteJ;
			siteJ = temp;
		}

		std::vector<RealVec> allPositions;

		// pushes back all positions of particles in water
		for (unsigned int i = 0; i < 3; i++) {
			allPositions.push_back(
					particlePositions[allParticleIndices[siteJ][i]] * nm_to_A);
			/* std::cout << "pushed back particle " << allParticleIndices[siteJ][i] << std::endl; */
		}
		// only push back the Cl position
		allPositions.push_back(
				particlePositions[allParticleIndices[siteI][0]] * nm_to_A);
		/* std::cout << "pushed back particle " << allParticleIndices[siteI][0] << std::endl; */
		std::vector<RealVec> extraPoints;

		extraPoints.resize(2);
		// FIXME support ions in image molecules
		//		if (_nonbondedMethod == CutoffPeriodic)
		//			imageMolecules(_periodicBoxDimensions * nm_to_A, allPositions);

		RealVec dOCl = allPositions[Oa] - allPositions[Cl];
		/* std::cout << "Is there a way to check this?" << dOCl << " visual check: looks good"<<  std::endl; */

		const double rOClsq = dOCl[0] * dOCl[0] + dOCl[1] * dOCl[1]
				+ dOCl[2] * dOCl[2];
		const double rOCl = std::sqrt(rOClsq);

		if (rOCl > r2f_chloride)
			return 0.0;

		// the extra-points
		monomer m;
		m.setup(allPositions[Oa], allPositions[Ha1], allPositions[Ha2],
				in_plane_gamma_chloride, out_of_plane_gamma_chloride,
				extraPoints[Xa1], extraPoints[Xa2]);

		//variables
		double v[8]; // stored separately (gets passed to poly::eval)
		variable ctxt[8];

		v[0] = ctxt[0].v_exp(d0_intra, k_HH_intra_chloride, allPositions[Ha1],
				allPositions[Ha2]);

		v[1] = ctxt[1].v_exp(d0_intra, k_OH_intra_chloride, allPositions[Oa],
				allPositions[Ha1]);
		v[2] = ctxt[2].v_exp(d0_intra, k_OH_intra_chloride, allPositions[Oa],
				allPositions[Ha2]);

		v[3] = ctxt[3].v_coul(d0_inter, k_ClH_coul, allPositions[Cl],
				allPositions[Ha1]);
		v[4] = ctxt[4].v_coul(d0_inter, k_ClH_coul, allPositions[Cl],
				allPositions[Ha2]);

		v[5] = ctxt[5].v_coul(d0_inter, k_ClO_coul, allPositions[Cl],
				allPositions[Oa]);

		v[6] = ctxt[6].v_exp(d0_inter, k_cllp_main, allPositions[Cl],
				extraPoints[Xa1]);
		v[7] = ctxt[7].v_exp(d0_inter, k_cllp_main, allPositions[Cl],
				extraPoints[Xa2]);

		double g[31];
		const double E_poly = h2o_cl::poly_2b_h2o_cl_v2x::eval(thefit_chloride,
				v, g);
		/* std::cout << "E_poly = " << E_poly << std::endl; */

		std::vector<RealVec> allForces;
		allForces.resize(allPositions.size());
		/* std::cout << "there are " << allForces.size() << " forces" << std::endl; */

		std::vector<RealVec> extraForces;
		extraForces.resize(extraPoints.size());
		/* std::cout << "there are " << extraForces.size() << " extra forces" << std::endl; */

		ctxt[0].grads(g[0], allForces[Ha1], allForces[Ha2]);

		ctxt[1].grads(g[1], allForces[Oa], allForces[Ha1]);
		ctxt[2].grads(g[2], allForces[Oa], allForces[Ha2]);

		ctxt[3].grads(g[3], allForces[Cl], allForces[Ha1]);
		ctxt[4].grads(g[4], allForces[Cl], allForces[Ha2]);

		ctxt[5].grads(g[5], allForces[Cl], allForces[Oa]);

		ctxt[6].grads(g[6], allForces[Cl], extraForces[Xa1]);
		ctxt[7].grads(g[7], allForces[Cl], extraForces[Xa2]);

		// distribute gradients w.r.t. the X-points
		m.grads(extraForces[Xa1], extraForces[Xa2], in_plane_gamma_chloride,
				out_of_plane_gamma_chloride, allForces[Oa], allForces[Ha1],
				allForces[Ha2]);
		// the switch
		double gsw;
		const double sw = f_switch(rOCl, gsw);
		/* std::cout << "sw = " << sw << std::endl; */
		double cal2joule = 4.184;

		// water molecule
		forces[allParticleIndices[siteJ][0]] += allForces[Oa] * sw * cal2joule
				* -10.;
		forces[allParticleIndices[siteJ][1]] += allForces[Ha1] * sw * cal2joule
				* -10.;
		forces[allParticleIndices[siteJ][2]] += allForces[Ha2] * sw * cal2joule
				* -10.;
		// Cl
		forces[allParticleIndices[siteI][0]] += allForces[Cl] * sw * cal2joule
				* -10.;
		/* std::cout << "forces has a size of " << forces.size() << std::endl; */
		// gradient of the switch
		gsw *= E_poly / rOCl;
		for (int i = 0; i < 3; ++i) {
			const double d = gsw * dOCl[i];
			forces[allParticleIndices[siteI][0]][i] += d * cal2joule * -10.;
			forces[allParticleIndices[siteJ][0]][i] -= d * cal2joule * -10.;
		}

		RealOpenMM energy = sw * E_poly * cal2joule;
		return energy;
	}
}

RealOpenMM MBPolReferenceTwoBodyForce::calculateForceAndEnergy(int numParticles,
		const vector<RealVec>& particlePositions,
		const std::vector<std::vector<int> >& allParticleIndices,
		const NeighborList& neighborList, vector<RealVec>& forces) const {
	/* std::cout << "Got to calculateForceAndEnergy" << std::endl; */

	// loop over neighbor list
	//    (1) calculate pair TwoBody ixn
	//    (2) accumulate forces: if particle is a site where interaction position != particle position,
	//        then call addReducedForce() to apportion force to particle and its covalent partner
	//        based on reduction factor
	RealOpenMM energy = 0.;
	for (unsigned int ii = 0; ii < neighborList.size(); ii++) {

		OpenMM::AtomPair pair = neighborList[ii];
		int siteI = pair.first;
		int siteJ = pair.second;

		energy += calculatePairIxn(siteI, siteJ, particlePositions,
				allParticleIndices, forces);

	}

	return energy;
}
