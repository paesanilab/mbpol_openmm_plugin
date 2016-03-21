extern "C" __global__ void recordInducedDipoles(const long long* __restrict__ fieldBuffers, const long long* __restrict__ fieldPolarBuffers,
		real* __restrict__ inducedDipole, real* __restrict__ inducedDipolePolar, const float* __restrict__ polarizability) {
	for (int atom = blockIdx.x*blockDim.x+threadIdx.x; atom < NUM_ATOMS; atom += gridDim.x*blockDim.x) {
		real scale = polarizability[atom]/(real) 0x100000000;
		inducedDipole[3*atom] = scale*fieldBuffers[atom];
		inducedDipole[3*atom+1] = scale*fieldBuffers[atom+PADDED_NUM_ATOMS];
		inducedDipole[3*atom+2] = scale*fieldBuffers[atom+PADDED_NUM_ATOMS*2];
		inducedDipolePolar[3*atom] = scale*fieldPolarBuffers[atom];
		inducedDipolePolar[3*atom+1] = scale*fieldPolarBuffers[atom+PADDED_NUM_ATOMS];
		inducedDipolePolar[3*atom+2] = scale*fieldPolarBuffers[atom+PADDED_NUM_ATOMS*2];
	}
}

/**
 * Normalize a vector and return what its magnitude was.
 */
inline __device__ real normVector(real3& v) {
	real n = SQRT(dot(v, v));
	v *= (n > 0 ? RECIP(n) : 0);
	return n;
}

/**
 * Compute the electrostatic potential at each of a set of points.
 */
extern "C" __global__ void computePotentialAtPoints(const real4* __restrict__ posq,
		const real* __restrict__ inducedDipole, const real4* __restrict__ points,
		real* __restrict__ potential, int numPoints, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
	extern __shared__ real4 localPosq[];
	real3* localInducedDipole = (real3*)&localPosq[blockDim.x];
	for (int basePoint = blockIdx.x*blockDim.x; basePoint < numPoints; basePoint += gridDim.x*blockDim.x) {
		int point = basePoint+threadIdx.x;
		real4 pointPos = points[point];
		real p = 0;
		for (int baseAtom = 0; baseAtom < NUM_ATOMS; baseAtom += blockDim.x) {
			int atom = baseAtom+threadIdx.x;

			// Load data into shared memory.

			if (atom < NUM_ATOMS) {
				localPosq[threadIdx.x] = posq[atom];
				localInducedDipole[threadIdx.x] = make_real3(inducedDipole[3*atom], inducedDipole[3*atom+1], inducedDipole[3*atom+2]);
			}
			__syncthreads();

			// Loop over atoms and compute the potential at this point.

			if (point < numPoints) {
				int end = min(blockDim.x, NUM_ATOMS-baseAtom);
				for (int i = 0; i < end; i++) {
					real3 delta = trimTo3(localPosq[i]-pointPos);
#ifdef USE_PERIODIC
					APPLY_PERIODIC_TO_DELTA(delta)
#endif
					real r2 = dot(delta, delta);
					real rInv = RSQRT(r2);
					p += localPosq[i].w*rInv;
					real rr2 = rInv*rInv;
					real rr3 = rInv*rr2;
					real scu = dot(localInducedDipole[i], delta);
					p -= (scu)*rr3;
				}
			}
			__syncthreads();
		}
		potential[point] = p;
	}
}

extern "C" __global__ void computeWaterCharge(
        real4* __restrict__ posq,
        real3* __restrict__ chargeDerivatives, unsigned int numMultipoles,
        const int* __restrict__ moleculeIndex, const int* __restrict__ atomType) {

        const unsigned int moleculeId = blockIdx.x*blockDim.x+threadIdx.x;

        // TODO make this flexible based on moleculeIndex and atomType
        const unsigned int O  = moleculeId * 4;
        const unsigned int H1 = O + 1;
        const unsigned int H2 = O + 2;
        const unsigned int M  = O + 3;

        if (M <= NUM_ATOMS) {
            const real Bohr_A = 0.52917721092; // CODATA 2010
            // M-site positioning (TTM2.1-F)
            const real gammaM = 0.426706882;

            const real gamma1 = 1.0 - gammaM;
            const real gamma2 = gammaM/2;
            const real ath0 = 1.82400520401572996557;
            const real costhe = -0.24780227221366464506;
            const real reoh = 0.958649;
            const real b1D = 1.0;
            const real a = 0.2999e0;
            const real b = -0.6932e0;
            const real c0 = 1.0099e0;
            const real c1 = -0.1801e0;
            const real c2 = 0.0892e0;


            const real e =  1.602176565e-19; // C CODATA 2010

            // interaction energy of 2 unit charges 1A apart
            const real E_cc = 1.0e-7*(c0*e*c0*e)/1.0e-10; // in J
            const real Na = 6.02214129e+23; // CODATA 2010
            const real kcal_J = 4184.0;

            const real CHARGECON = SQRT(E_cc*Na/kcal_J);

            const unsigned int idxD0[84] = {
                   1, 1, 1, 2, 1, 1, 1, 2, 2, 3, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4,
                   1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 1, 1, 1, 1, 1,
                   1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6, 1, 1, 1, 1,
                   1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5,
                   5, 6, 6, 7
            };

            const unsigned int idxD1[84] = {
                   1, 1, 2, 1, 1, 2, 3, 1, 2, 1, 1, 2, 3, 4, 1, 2, 3, 1, 2, 1,
                   1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2, 3, 1, 2, 1, 1, 2, 3, 4, 5,
                   6, 1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2, 3, 1, 2, 1, 1, 2, 3, 4,
                   5, 6, 7, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2,
                   3, 1, 2, 1
            };

            const unsigned int idxD2[84] = {
                   1, 2, 1, 1, 3, 2, 1, 2, 1, 1, 4, 3, 2, 1, 3, 2, 1, 2, 1, 1,
                   5, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2, 1, 2, 1, 1, 6, 5, 4, 3, 2,
                   1, 5, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2, 1, 2, 1, 1, 7, 6, 5, 4,
                   3, 2, 1, 6, 5, 4, 3, 2, 1, 5, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2,
                   1, 2, 1, 1
            };


            const real coefD[84] = {
                  -2.1689686086730e-03, 1.4910379754728e-02, 5.3546078430060e-02,
                  -7.4055995388666e-02,-3.7764333017616e-03, 1.4089887256484e-01,
                  -6.2584207687264e-02,-1.1260393113022e-01,-5.7824159269319e-02,
                   1.4360743650655e-02,-1.5469680141070e-02,-1.3036350092795e-02,
                   2.7515837781556e-02, 1.4098478875076e-01,-2.7663168397781e-02,
                  -5.2378176254797e-03,-1.0237198381792e-02, 8.9571999265473e-02,
                   7.2920263098603e-03,-2.6873260551686e-01, 2.0220870325864e-02,
                  -7.0764766270927e-02, 1.2140640273760e-01, 2.0978491966341e-02,
                  -1.9443840512668e-01, 4.0826835370618e-02,-4.5365190474650e-02,
                   6.2779900072132e-02,-1.3194351021000e-01,-1.4673032718563e-01,
                   1.1894031277247e-01,-6.4952851564679e-03, 8.8503610374493e-02,
                   1.4899437409291e-01, 1.3962841511565e-01,-2.6459446720450e-02,
                  -5.0128914532773e-02, 1.8329676428116e-01,-1.5559089125095e-01,
                  -4.0176879767592e-02, 3.6192059996636e-01, 1.0202887240343e-01,
                   1.9318668580051e-01,-4.3435977107932e-01,-4.2080828803311e-02,
                   1.9144626027273e-01,-1.7851138969948e-01, 1.0524533875070e-01,
                  -1.7954071602185e-02, 5.2022455612120e-02,-2.8891891146828e-01,
                  -4.7452036576319e-02,-1.0939400546289e-01, 3.5916564473568e-01,
                  -2.0162789820172e-01,-3.5838629543696e-01, 5.6706523551202e-03,
                   1.3849337488211e-01,-4.1733982195604e-01, 4.1641570764241e-01,
                  -1.2243429796296e-01, 4.7141730971228e-02,-1.8224510249551e-01,
                  -1.8880981556620e-01,-3.1992359561800e-01,-1.8567550546587e-01,
                   6.1850530431280e-01,-6.1142756235141e-02,-1.6996135584933e-01,
                   5.4252879499871e-01, 6.6128603899427e-01, 1.2107016404639e-02,
                  -1.9633639729189e-01, 2.7652059420824e-03,-2.2684111109778e-01,
                  -4.7924491598635e-01, 2.4287790137314e-01,-1.4296023329441e-01,
                   8.9664665907006e-02,-1.4003228575602e-01,-1.3321543452254e-01,
                  -1.8340983193745e-01, 2.3426707273520e-01, 1.5141050914514e-01
            };

            real3 ROH1 = (trimTo3(posq[H1]) - trimTo3(posq[ O])) * 10.;
            real3 ROH2 = (trimTo3(posq[H2]) - trimTo3(posq[ O])) * 10.;
            real3 RHH  = (trimTo3(posq[H1]) - trimTo3(posq[H2])) * 10.;
            real dROH1 = SQRT(dot(ROH1, ROH1));
            real dROH2 = SQRT(dot(ROH2, ROH2));
            real dRHH  = SQRT(dot( RHH,  RHH));

            const real costh =
                dot(ROH1, ROH2)/(dROH1*dROH2);

            const real efac = EXP(-b1D*(POW((dROH1 - reoh), 2)
                                             + POW((dROH2 - reoh), 2)));

            const real x1 = (dROH1 - reoh)/reoh;
            const real x2 = (dROH2 - reoh)/reoh;
            const real x3 = costh - costhe;

            real fmat[3][16];

            for (unsigned int i = 0; i < 3; ++i) {
                fmat[i][0] = 0.0;
                fmat[i][1] = 1.0;
            }

            for (unsigned int j = 2; j < 16; ++j) {
                fmat[0][j] = fmat[0][j - 1]*x1;
                fmat[1][j] = fmat[1][j - 1]*x2;
                fmat[2][j] = fmat[2][j - 1]*x3;
            }

            // Calculate the dipole moment

            real p1(0), p2(0);
            real pl1 = costh;
            real pl2 = 0.5*(3*pl1*pl1 - 1.0);

            real dp1dr1(0);
            real dp1dr2(0);
            real dp1dcabc(0);
            real dp2dr1(0);
            real dp2dr2(0);
            real dp2dcabc(0);

            for (unsigned int j = 1; j < 84; ++j) {
                const unsigned int inI = idxD0[j];
                const unsigned int inJ = idxD1[j];
                const unsigned int inK = idxD2[j];

                p1 += coefD[j]*fmat[0][inI]*fmat[1][inJ]*fmat[2][inK];
                p2 += coefD[j]*fmat[0][inJ]*fmat[1][inI]*fmat[2][inK];

                dp1dr1 +=
                    coefD[j]*(inI - 1)*fmat[0][inI - 1]*fmat[1][inJ]*fmat[2][inK];
                dp1dr2 +=
                    coefD[j]*(inJ - 1)*fmat[0][inI]*fmat[1][inJ - 1]*fmat[2][inK];
                dp1dcabc +=
                    coefD[j]*(inK - 1)*fmat[0][inI]*fmat[1][inJ]*fmat[2][inK - 1];
                dp2dr1 +=
                    coefD[j]*(inJ - 1)*fmat[0][inJ - 1]*fmat[1][inI]*fmat[2][inK];
                dp2dr2 +=
                    coefD[j]*(inI - 1)*fmat[0][inJ]*fmat[1][inI - 1]*fmat[2][inK];
                dp2dcabc +=
                    coefD[j]*(inK - 1)*fmat[0][inJ]*fmat[1][inI]*fmat[2][inK - 1];
            }

            const real xx = Bohr_A;
            const real xx2 = xx*xx;
            dp1dr1 /= reoh/xx;
            dp1dr2 /= reoh/xx;
            dp2dr1 /= reoh/xx;
            dp2dr2 /= reoh/xx;

            const real pc0 =
                a*(POW(dROH1, b) + POW(dROH2, b))*(c0 + pl1*c1 + pl2*c2);

            const real dpc0dr1 =
                a*b*POW(dROH1, b - 1)*(c0 + pl1*c1 + pl2*c2)*xx2;
            const real dpc0dr2 =
                a*b*POW(dROH2, b - 1)*(c0 + pl1*c1 + pl2*c2)*xx2;
            real dpc0dcabc =
                a*(POW(dROH1, b) + POW(dROH2, b))*(c1 + 0.5*(6.0*pl1)*c2)*xx;

            const real defacdr1 = -2.0*b1D*(dROH1 - reoh)*efac*xx;
            const real defacdr2 = -2.0*b1D*(dROH2 - reoh)*efac*xx;

            dp1dr1 = dp1dr1*efac + p1*defacdr1 + dpc0dr1;
            dp1dr2 = dp1dr2*efac + p1*defacdr2 + dpc0dr2;
            dp1dcabc = dp1dcabc*efac + dpc0dcabc;
            dp2dr1 = dp2dr1*efac + p2*defacdr1 + dpc0dr1;
            dp2dr2 = dp2dr2*efac + p2*defacdr2 + dpc0dr2;
            dp2dcabc = dp2dcabc*efac + dpc0dcabc;

            p1 = coefD[0] + p1*efac + pc0*xx; // q^H1 in TTM2-F
            p2 = coefD[0] + p2*efac + pc0*xx; // q^H2 paper

            real chargeO= -(p1 + p2);  // Oxygen
            real chargeH1 = p1; // Hydrogen-1
            real chargeH2 = p2;  // Hydrogen-2
            real gamma2div1 = gamma2/gamma1;

            posq[O].w = 0.;
            posq[H1].w = chargeH1 + gamma2div1*(chargeH1 + chargeH2);
            posq[H2].w = chargeH2 + gamma2div1*(chargeH1 + chargeH2);
            posq[M].w = chargeO/gamma1;

            dp1dr1 /= xx;
            dp1dr2 /= xx;
            dp2dr1 /= xx;
            dp2dr2 /= xx;

            const real f1q1r13 = (dp1dr1 - (dp1dcabc*costh/dROH1))/dROH1;
            const real f1q1r23 = dp1dcabc/(dROH1*dROH2);
            const real f2q1r23 = (dp1dr2 - (dp1dcabc*costh/dROH2))/dROH2;
            const real f2q1r13 = dp1dcabc/(dROH2*dROH1);
            const real f1q2r13 = (dp2dr1 - (dp2dcabc*costh/dROH1))/dROH1;
            const real f1q2r23 = dp2dcabc/(dROH1*dROH2);
            const real f2q2r23 = (dp2dr2 - (dp2dcabc*costh/dROH2))/dROH2;
            const real f2q2r13 = dp2dcabc/(dROH2*dROH1);

            // first index is atom w.r.t. to which the derivative is
            // second index is the charge being differentiated

            real3 chargeDerivativesH1vsH1;
            real3 chargeDerivativesH1vsH2;
            real3 chargeDerivativesH1vsO;
            //gradient of charge h1(second index) wrt displacement of h1(first index)
            chargeDerivativesH1vsH1 = f1q1r13*ROH1 + f1q1r23*ROH2;
            chargeDerivativesH1vsH2 = f2q1r13*ROH1 + f2q1r23*ROH2;
            chargeDerivativesH1vsO = -(chargeDerivativesH1vsH1+chargeDerivativesH1vsH2);

            real3 chargeDerivativesH2vsH1;
            real3 chargeDerivativesH2vsH2;
            real3 chargeDerivativesH2vsO;

            //    //gradient of charge h1(second index) wrt displacement of h1(first index)
            chargeDerivativesH2vsH1 = f1q2r13*ROH1 + f1q2r23*ROH2;
            chargeDerivativesH2vsH2 = f2q2r13*ROH1 + f2q2r23*ROH2;
            chargeDerivativesH2vsO = -(chargeDerivativesH2vsH1+chargeDerivativesH2vsH2);

            real3 chargeDerivativesOvsH1;
            real3 chargeDerivativesOvsH2;
            real3 chargeDerivativesOvsO;


            chargeDerivativesOvsH1 = -(chargeDerivativesH1vsH1+ chargeDerivativesH2vsH1);
            chargeDerivativesOvsH2 =  -(chargeDerivativesH1vsH2+ chargeDerivativesH2vsH2);
            chargeDerivativesOvsO =  -(chargeDerivativesH1vsO+ chargeDerivativesH2vsO);

            chargeDerivatives[3*M]   = make_real3(0, 0, 0);
            chargeDerivatives[3*M+1] = make_real3(0, 0, 0);
            chargeDerivatives[3*M+2] = make_real3(0, 0, 0);

            real3 sumH1 = gamma2div1*(chargeDerivativesH1vsH1+chargeDerivativesH2vsH1);
            real3 sumH2 = gamma2div1*(chargeDerivativesH1vsH2+chargeDerivativesH2vsH2);
            real3 sumO = gamma2div1*(chargeDerivativesH1vsO+chargeDerivativesH2vsO);

            // vs H1
            chargeDerivatives[3*H1] = chargeDerivativesH1vsH1 + sumH1;
            chargeDerivatives[3*H2] = chargeDerivativesH1vsH2 + sumH2;
            chargeDerivatives[3*O]  = chargeDerivativesH1vsO + sumO;

            // vs H2
            chargeDerivatives[3*H1+1] = chargeDerivativesH2vsH1 + sumH1;
            chargeDerivatives[3*H2+1] = chargeDerivativesH2vsH2 + sumH2;
            chargeDerivatives[3*O+1]  = chargeDerivativesH2vsO +  sumO;

            // vs M
            chargeDerivatives[3*H1+2] = chargeDerivativesOvsH1 - 2*sumH1;
            chargeDerivatives[3*H2+2] = chargeDerivativesOvsH2 - 2*sumH2;
            chargeDerivatives[3*O+2 ] = chargeDerivativesOvsO  - 2*sumO;

            //// convert from q/A to q/nm
            for (unsigned int i = O; i <= H2; ++i) {
                for (unsigned int s = 0; s < 3; ++s) {
                    chargeDerivatives[3*i+s] *= 10;
                }

            }

        }

}

extern "C" __global__ void computeChargeDerivativesForces(
        const real3* __restrict__ chargeDerivatives, unsigned int numMultipoles,
        unsigned long long* __restrict__ forceBuffers, const long long* __restrict__ potentialBuffers) {

        const unsigned int moleculeId = blockIdx.x*blockDim.x+threadIdx.x;

        // TODO make this flexible based on moleculeIndex and atomType
        const unsigned int O  = moleculeId * 4;
        const unsigned int H1 = O + 1;
        const unsigned int H2 = O + 2;
        const unsigned int M  = O + 3;

        const real scale = RECIP((real) 0x100000000);

        if (M <= numMultipoles) {

            real3 f;
            for( unsigned int a = O; a <= H2; a++ ){
                f = chargeDerivatives[3*a]   * ((real) potentialBuffers[H1])*scale +
                     chargeDerivatives[3*a+1] * ((real) potentialBuffers[H2])*scale +
                     chargeDerivatives[3*a+2] * ((real) potentialBuffers[M])*scale;
                f *= - ENERGY_SCALE_FACTOR;
                forceBuffers[a]                    += static_cast<unsigned long long>((long long) (f.x*0x100000000));
                forceBuffers[a+PADDED_NUM_ATOMS]   += static_cast<unsigned long long>((long long) (f.y*0x100000000));
                forceBuffers[a+2*PADDED_NUM_ATOMS] += static_cast<unsigned long long>((long long) (f.z*0x100000000));
            }

        }
}
