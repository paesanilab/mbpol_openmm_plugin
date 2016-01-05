__device__ void computeOneInteractionF1(AtomData& atom1, volatile AtomData& atom2, float dScale, float pScale, float mScale, real& energy, real3& outputForce) {

    // FIXME thole copy in unique location
    const enum TholeIndices { TCC, TCD, TDD, TDDOH, TDDHH };
    const float thole[5] =  { 0.4, 0.4, 0.4,   0.4,   0.4 };
	// thole[TDD] = 0.055;
	// thole[TDDOH] = 0.626;
	// thole[TDDHH] = 0.055;

	// deltas
    real3 delta = make_real3(atom2.posq.x - atom1.posq.x, atom2.posq.y - atom1.posq.y, atom2.posq.z - atom1.posq.z);

    real r2 = dot(delta, delta);
    real r = SQRT(r2);
    real rr1 = RECIP(r);
    real rr2 = rr1*rr1;
    real rr3 = rr1*rr2;
    real rr5 = 3*rr3*rr2;
    real rr7 = 5*rr5*rr2;

    real sci3 = dot(atom1.inducedDipole, delta);
    real sci4 = dot(atom2.inducedDipole, delta);

    real scip2 = dot(atom1.inducedDipole,atom2.inducedDipolePolar) +
                 dot(atom2.inducedDipole,atom1.inducedDipolePolar);

    real scip3 = dot(atom1.inducedDipolePolar, delta);

    real scip4 = dot(atom2.inducedDipolePolar, delta);

    real gli1 = atom2.posq.w*sci3 - atom1.posq.w*sci4;

    real glip1 = atom2.posq.w*scip3 - atom1.posq.w*scip4;

    real gl0 = atom1.posq.w*atom2.posq.w;

    // if isSameWater set gl0 to zero
    bool isSameWater = atom1.moleculeIndex == atom2.moleculeIndex;
    gl0  *= !isSameWater;
    gli1 *= !isSameWater;
    glip1 *= !isSameWater;

    // real scale1CC = getAndScaleInverseRs( particleI, particleK, r, true, 1, TCC);
    // real scale3CD = getAndScaleInverseRs( particleI, particleK, r, true, 3, TCD);

    float damp      = POW(atom1.damp*atom2.damp, 1.0f/6.0f); // AA in MBPol

    real do_scaling = (damp != 0.0) & ( damp > -50.0 ); // damp or not

    real ratio       = POW(r/damp, 4.); // rA4 in MBPol

    real dampForExpCC = -1 * thole[TCC] * ratio;
    // EXP(ttm::gammln(3.0/4.0)) = 1.2254167024651776
    real scale3CC = ( 1.0 - do_scaling*EXP(dampForExpCC) ); // needed for force
    real scale1CC = scale3CC + do_scaling*(POW(thole[TCC], 1.0f/4.0f)*(r/damp)*1.2254167024651776*gammq(3.0/4.0, -dampForExpCC));

    real dampForExpCD = -1 * thole[TCD] * ratio;
    real scale3CD = ( 1.0 - do_scaling*EXP(dampForExpCD) );

    real em = rr1 * gl0 * scale1CC;
    real ei = 0.5f * gli1 * rr3 * scale3CD;
    energy = em+ei;

    energy *= ENERGY_SCALE_FACTOR;

    // RealOpenMM scale3CC = getAndScaleInverseRs( particleI, particleK, r, true, 3, TCC);
    // RealOpenMM scale5CD = getAndScaleInverseRs( particleI, particleK, r, true, 5, TCD);
    // RealOpenMM scale5DD = getAndScaleInverseRs( particleI, particleK, r, true, 5, TDD);
    // RealOpenMM scale7DD = getAndScaleInverseRs( particleI, particleK, r, true, 7, TDD);

    real scale5CD = scale3CD - do_scaling * (4./3.) * thole[TCD] * EXP(dampForExpCD) * ratio;

    int tdd = TDD;
    if (isSameWater) {
        if ((atom1.atomType == 0) | (atom2.atomType == 0)) { // one is oxygen
            tdd = TDDOH;
        } else { // both hydrogens
            tdd = TDDHH;
        }
    }
    real dampForExpDD = -1 * thole[tdd] * ratio;
    real scale5DD =  1.0 - do_scaling * EXP(dampForExpDD) *  (1. + (4./3.) * thole[tdd] * ratio);
    real scale7DD = scale5DD - do_scaling * ((4./15.) * thole[tdd] * (4. * thole[tdd] * ratio - 1.) * EXP(dampForExpDD) / POW(damp, 4.0f) * POW(r, 4));

    real gf1 = rr3*gl0*scale3CC;

    real3 ftm2 = gf1*delta;

    // intermediate variables for the induced components

    real gfi1 = 0.5 * rr5 * gli1 * scale5CD +  // charge - induced dipole
                0.5 * rr5 * glip1* scale5CD  +  // charge - induced dipole
                0.5 * rr5 * scip2* scale5DD  +  //induced dipole - induced dipole
              - 0.5 * rr7 * (sci3*scip4+scip3*sci4) *scale7DD; // induced dipole - induced dipole

    real3 ftm2i = gfi1*delta;

    ftm2i += ( atom1.inducedDipolePolar *  sci4 + // iPdipole_i * idipole_k
               atom1.inducedDipole      * scip4 +
               atom2.inducedDipolePolar *  sci3 + // iPdipole_k * idipole_i
               atom2.inducedDipole * scip3  ) * 0.5 * rr5 * scale5DD;

    // Same water atoms have no induced-dipole/charge interaction
    if (not( isSameWater )) {
    ftm2i += (
                   ( atom1.inducedDipole +
             atom1.inducedDipolePolar )*-atom2.posq.w +
                   ( atom2.inducedDipole +
             atom2.inducedDipolePolar )* atom1.posq.w
                 ) * 0.5 * rr3 * scale3CD;
    }

    outputForce = - ENERGY_SCALE_FACTOR * (ftm2+ftm2i);
}
