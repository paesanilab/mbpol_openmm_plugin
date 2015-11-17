
__device__ real getAndScaleInverseRs(AtomData& atom1, volatile AtomData& atom2, real r, bool justScale, int interactionOrder, int interactionType)
{	
	// a way to see if they are the same molecule
	float thole[] = {0.4, 0.4, 0.055, 0.626, 0.055};
	enum TholeIndices { TCC, TCD, TDD, TDDOH, TDDHH };
	
	real rrI;
    if (justScale) {
        rrI         = 1.;
    } else {
        real rI             =  1.0/r;
        real r2I            =  rI*rI;
        rrI                    = rI;
        for( unsigned int ii  = 3; ii <= interactionOrder; ii=ii+2 ){
           rrI *= (ii-2)*r2I;
        }
    }
	
	float thole = 0.4;
	real pGamma = thole;
    
    if (interactionType == TDD) {
    	 bool oneIsOxygen = ;
    	 if (oneIsOxygen) {
                pgamma = thole[TDDOH];
            } 
            else {
                pgamma = thole[TDDHH];
            }
            else {
            pgamma = thole[TDD];
            }
        }
    }
        
	// get and scale move to diffrent kernel func
    real pdamp = atom1.damp*atom2.damp;
    if (pdamp != 0) {
   
        real ratio = r/pdamp;
        //FIXME: 
        float thole = 0.4;
        float pGamma = thole;

        real damp = ratio*ratio*ratio*pGamma;
        real dampExp = EXP(-damp);
        real damp1 = damp + 1;
        real damp2 = damp*damp;

        scale3 = 1 - dampExp;
        scale5 = 1 - damp1*dampExp;
        scale7 = 1 - (damp1 + 0.6f*damp2)*dampExp;

        real factor = 3*damp*dampExp*rr2;
        real factor7 = -0.2f + 0.6f*damp;
        
        ddsc3_0 = factor*xr;
        ddsc5_0 = ddsc3_0*damp;
        ddsc7_0 = ddsc5_0*factor7;

        ddsc3_1 = factor*yr;
        ddsc5_1 = ddsc3_1*damp;
        ddsc7_1 = ddsc5_1*factor7;

        ddsc3_2 = factor*zr;
        ddsc5_2 = ddsc3_2*damp;
        ddsc7_2 = ddsc5_2*factor7;


    }
	return 0;
}

__device__ void computeOneInteractionF1(AtomData& atom1, volatile AtomData& atom2, float dScale, float pScale, float mScale, real& energy, real3& outputForce) {    

    const float uScale = 1;
    real ddsc3_0 = 0;
    real ddsc3_1 = 0;
    real ddsc3_2 = 0;

    real ddsc5_0 = 0;
    real ddsc5_1 = 0;
    real ddsc5_2 = 0;

    real ddsc7_0 = 0;
    real ddsc7_1 = 0;
    real ddsc7_2 = 0;
    
	// deltas 
    real xr = atom2.posq.x - atom1.posq.x;
    real yr = atom2.posq.y - atom1.posq.y;
    real zr = atom2.posq.z - atom1.posq.z;
    
    real r2 = xr*xr + yr*yr + zr*zr;
    real r = SQRT(r2);
    real rr1 = RECIP(r);
    real rr2 = rr1*rr1;
    real rr3 = rr1*rr2;
    real rr5 = 3*rr3*rr2;
    real rr7 = 5*rr5*rr2;

    real scale3 = 1;
    real scale5 = 1;
    real scale7 = 1;
    
	// get and scale move to diffrent kernel func
    real pdamp = atom1.damp*atom2.damp;
    if (pdamp != 0) {
   
        real ratio = r/pdamp;
        //FIXME: 
        float thole = 0.4;
        float pGamma = thole;

        real damp = ratio*ratio*ratio*pGamma;
        real dampExp = EXP(-damp);
        real damp1 = damp + 1;
        real damp2 = damp*damp;

        scale3 = 1 - dampExp;
        scale5 = 1 - damp1*dampExp;
        scale7 = 1 - (damp1 + 0.6f*damp2)*dampExp;

        real factor = 3*damp*dampExp*rr2;
        real factor7 = -0.2f + 0.6f*damp;
        
        ddsc3_0 = factor*xr;
        ddsc5_0 = ddsc3_0*damp;
        ddsc7_0 = ddsc5_0*factor7;

        ddsc3_1 = factor*yr;
        ddsc5_1 = ddsc3_1*damp;
        ddsc7_1 = ddsc5_1*factor7;

        ddsc3_2 = factor*zr;
        ddsc5_2 = ddsc3_2*damp;
        ddsc7_2 = ddsc5_2*factor7;


    }
      
    real scale3i = rr3*scale3*uScale;
    real scale5i = rr5*scale5*uScale;

    real dsc3 = rr3*scale3*dScale;
    real psc3 = rr3*scale3*pScale;

    real dsc5 = rr5*scale5*dScale;
    real psc5 = rr5*scale5*pScale;

    real dsc7 = rr7*scale7*dScale;
    real psc7 = rr7*scale7*pScale;

    real sc2 = atom1.dipole.x*atom2.dipole.x + atom1.dipole.y*atom2.dipole.y + atom1.dipole.z*atom2.dipole.z;

    real sc4 = atom2.dipole.x*xr + atom2.dipole.y*yr + atom2.dipole.z*zr;

    real sc3 = atom1.dipole.x*xr + atom1.dipole.y*yr + atom1.dipole.z*zr;
    
    real sci1 = atom1.inducedDipole.x*atom2.dipole.x + atom1.inducedDipole.y*atom2.dipole.y + atom1.inducedDipole.z*atom2.dipole.z +
                atom2.inducedDipole.x*atom1.dipole.x + atom2.inducedDipole.y*atom1.dipole.y + atom2.inducedDipole.z*atom1.dipole.z;
        
    real sci3 = atom1.inducedDipole.x*xr + atom1.inducedDipole.y*yr + atom1.inducedDipole.z*zr;
    real sci4 = atom2.inducedDipole.x*xr + atom2.inducedDipole.y*yr + atom2.inducedDipole.z*zr;
    
    real scip1 = atom1.inducedDipolePolar.x*atom2.dipole.x + atom1.inducedDipolePolar.y*atom2.dipole.y + atom1.inducedDipolePolar.z*atom2.dipole.z +
                 atom2.inducedDipolePolar.x*atom1.dipole.x + atom2.inducedDipolePolar.y*atom1.dipole.y + atom2.inducedDipolePolar.z*atom1.dipole.z;

    real scip2 = atom1.inducedDipole.x*atom2.inducedDipolePolar.x + atom1.inducedDipole.y*atom2.inducedDipolePolar.y + atom1.inducedDipole.z*atom2.inducedDipolePolar.z +
                 atom2.inducedDipole.x*atom1.inducedDipolePolar.x + atom2.inducedDipole.y*atom1.inducedDipolePolar.y + atom2.inducedDipole.z*atom1.inducedDipolePolar.z;
    
    real scip3 = ((atom1.inducedDipolePolar.x)*(xr) + (atom1.inducedDipolePolar.y)*(yr) + (atom1.inducedDipolePolar.z)*(zr));

    real scip4 = ((atom2.inducedDipolePolar.x)*(xr) + (atom2.inducedDipolePolar.y)*(yr) + (atom2.inducedDipolePolar.z)*(zr));

    real gli1 = atom2.posq.w*sci3 - atom1.posq.w*sci4;
    
    real gli6 = sci1;
    real glip1 = atom2.posq.w*scip3 - atom1.posq.w*scip4;
    real glip6 = scip1;
    real gli2 = -sc3*sci4 - sci3*sc4;
    
    real glip2 = -sc3*scip4 - scip3*sc4;
    real factor3 = rr3*((gli1  +  gli6)*pScale + (glip1  + glip6)*dScale);
    real factor5 = rr5*(gli2*pScale + glip2*dScale);
    
    real ftm2i_0 = -0.5f*(factor3*ddsc3_0 + factor5*ddsc5_0);
    real ftm2i_1 = -0.5f*(factor3*ddsc3_1 + factor5*ddsc5_1);
    real ftm2i_2 = -0.5f*(factor3*ddsc3_2 + factor5*ddsc5_2);
      
    real gl0 = atom1.posq.w*atom2.posq.w;
    real gl1 = atom2.posq.w*sc3 - atom1.posq.w*sc4;
    real gl2 = -sc3*sc4;
    real gl6 = sc2;
    
    bool isSameWater = (atom1.posq.x == atom2.posq.x && atom1.posq.y == atom2.posq.y && atom1.posq.z == atom2.posq.z && atom1.posq.w == atom2.posq.w);
    if (isSameWater) {
    	gl0 = 0;
    	//gli0 = 0;
    	//glip0 = 0;
    }
    real gf1 = rr3*gl0 + rr5*(gl1+gl6) + rr7*gl2;

    real gf2 = -atom2.posq.w*rr3 + sc4*rr5;
    real gf5 = 2*(-atom2.posq.w*rr5+sc4*rr7);

    real gf3 =  atom1.posq.w*rr3 + sc3*rr5;
    real gf6 = 2*(-atom1.posq.w*rr5-sc3*rr7);

    real em = mScale*(rr1*gl0 + rr3*(gl1+gl6) + rr5*gl2);
    real ei = 0.5f*((gli1+gli6)*psc3 + gli2*psc5);
    energy = em+ei;
    
    real ftm2_0 = mScale*(gf1*xr + gf2*atom1.dipole.x + gf3*atom2.dipole.x);
    real ftm2_1 = mScale*(gf1*yr + gf2*atom1.dipole.y + gf3*atom2.dipole.y);
    real ftm2_2 = mScale*(gf1*zr + gf2*atom1.dipole.z + gf3*atom2.dipole.z);

    real gfi1 = rr2*(1.5f*((gli1+gli6)*psc3 + (glip1+glip6)*dsc3 + scip2*scale3i) + 2.5f*(gli2*psc5 + glip2*dsc5 - (sci3*scip4+scip3*sci4)*scale5i));
    ftm2i_0 += gfi1*xr;
    ftm2i_1 += gfi1*yr;
    ftm2i_2 += gfi1*zr;

    real gfi5 = (sci4*psc7 + scip4*dsc7);

    real gfi6 = -(sci3*psc7 + scip3*dsc7);

    ftm2i_0 += 0.5f*(-atom2.posq.w*(atom1.inducedDipole.x*psc3 + atom1.inducedDipolePolar.x*dsc3) +
               sc4*(atom1.inducedDipole.x*psc5 + atom1.inducedDipolePolar.x*dsc5)) +
      
               0.5f*(atom1.posq.w*(atom2.inducedDipole.x*psc3+atom2.inducedDipolePolar.x*dsc3) +
               sc3*(atom2.inducedDipole.x*psc5 +atom2.inducedDipolePolar.x*dsc5)) +

               scale5i*(sci4*atom1.inducedDipolePolar.x+scip4*atom1.inducedDipole.x +
                        sci3*atom2.inducedDipolePolar.x+scip3*atom2.inducedDipole.x)*0.5f +
      
               0.5f*(sci4*psc5+scip4*dsc5)*atom1.dipole.x +
               0.5f*(sci3*psc5+scip3*dsc5)*atom2.dipole.x;
      
    ftm2i_1 += 0.5f*(-atom2.posq.w*(atom1.inducedDipole.y*psc3 + atom1.inducedDipolePolar.y*dsc3) +
               sc4*(atom1.inducedDipole.y*psc5 + atom1.inducedDipolePolar.y*dsc5)) +

               (atom1.posq.w*(atom2.inducedDipole.y*psc3+atom2.inducedDipolePolar.y*dsc3) +
                    sc3*(atom2.inducedDipole.y*psc5+atom2.inducedDipolePolar.y*dsc5))*0.5f +
                    scale5i*(sci4*atom1.inducedDipolePolar.y+scip4*atom1.inducedDipole.y + sci3*atom2.inducedDipolePolar.y+scip3*atom2.inducedDipole.y)*0.5f +

               0.5f*(sci4*psc5+scip4*dsc5)*atom1.dipole.y +
               0.5f*(sci3*psc5+scip3*dsc5)*atom2.dipole.y;
      
    ftm2i_2 += 0.5f*(-atom2.posq.w*(atom1.inducedDipole.z*psc3 + atom1.inducedDipolePolar.z*dsc3) +
               sc4*(atom1.inducedDipole.z*psc5 + atom1.inducedDipolePolar.z*dsc5)) +

               (atom1.posq.w*(atom2.inducedDipole.z*psc3+atom2.inducedDipolePolar.z*dsc3) +
                    sc3*(atom2.inducedDipole.z*psc5+atom2.inducedDipolePolar.z*dsc5))*0.5f +
                    scale5i*(sci4*atom1.inducedDipolePolar.z+scip4*atom1.inducedDipole.z +
                    sci3*atom2.inducedDipolePolar.z+scip3*atom2.inducedDipole.z)*0.5f +

               0.5f*(sci4*psc5+scip4*dsc5)*atom1.dipole.z +
               0.5f*(sci3*psc5+scip3*dsc5)*atom2.dipole.z;

#ifdef DIRECT_POLARIZATION
    real gfd = 0.5*(3*rr2*scip2*scale3i - 5*rr2*(scip3*sci4+sci3*scip4)*scale5i);
    real temp5 = 0.5*scale5i;
    real fdir_0 = gfd*xr + temp5*(sci4*atom1.inducedDipolePolar.x + scip4*atom1.inducedDipole.x + sci3*atom2.inducedDipolePolar.x + scip3*atom2.inducedDipole.x);
    real fdir_1 = gfd*yr + temp5*(sci4*atom1.inducedDipolePolar.y + scip4*atom1.inducedDipole.y + sci3*atom2.inducedDipolePolar.y + scip3*atom2.inducedDipole.y);
    real fdir_2 = gfd*zr + temp5*(sci4*atom1.inducedDipolePolar.z + scip4*atom1.inducedDipole.z + sci3*atom2.inducedDipolePolar.z + scip3*atom2.inducedDipole.z);
    ftm2i_0 -= fdir_0;
    ftm2i_1 -= fdir_1;
    ftm2i_2 -= fdir_2;
#else
    real scaleF = 0.5f*uScale;
    real inducedFactor3 = scip2*rr3*scaleF;
    real inducedFactor5 = (sci3*scip4+scip3*sci4)*rr5*scaleF;
    real findmp_0 = inducedFactor3*ddsc3_0 - inducedFactor5*ddsc5_0;
    real findmp_1 = inducedFactor3*ddsc3_1 - inducedFactor5*ddsc5_1;
    real findmp_2 = inducedFactor3*ddsc3_2 - inducedFactor5*ddsc5_2;
    ftm2i_0 -= findmp_0;
    ftm2i_1 -= findmp_1;
    ftm2i_2 -= findmp_2;
#endif

    outputForce.x = -(ftm2_0+ftm2i_0);
    outputForce.y = -(ftm2_1+ftm2i_1);
    outputForce.z = -(ftm2_2+ftm2i_2);
}
