extern "C" __global__ void computeLabFrameMoments(real4* __restrict__ posq, int4* __restrict__ multipoleParticles, float* __restrict__ molecularDipoles,
		 real* __restrict__ labFrameDipoles) {
	// get coordinates of this atom and the z & x axis atoms
	// compute the vector between the atoms and 1/sqrt(d2), d2 is distance between
	// this atom and the axis atom

	// this atom is referred to as the k-atom in notes below

	// code common to ZThenX and Bisector

	for (int atom = blockIdx.x*blockDim.x+threadIdx.x; atom < NUM_ATOMS; atom += gridDim.x*blockDim.x) {
		int4 particles = multipoleParticles[atom];
		if (particles.x >= 0 && particles.z >= 0) {
			real4 thisParticlePos = posq[atom];
			real4 posZ = posq[particles.z];
			real3 vectorZ = make_real3(posZ.x-thisParticlePos.x, posZ.y-thisParticlePos.y, posZ.z-thisParticlePos.z);
			real4 posX = posq[particles.x];
			real3 vectorX = make_real3(posX.x-thisParticlePos.x, posX.y-thisParticlePos.y, posX.z-thisParticlePos.z);
			int axisType = particles.w;

			/*
			 z-only
			 (1) norm z
			 (2) select random x
			 (3) x = x - (x.z)z
			 (4) norm x
			 
			 z-then-x
			 (1) norm z
			 (2) norm x (not needed)
			 (3) x = x - (x.z)z
			 (4) norm x
			 
			 bisector
			 (1) norm z
			 (2) norm x 
			 (3) z = x + z
			 (4) norm z
			 (5) x = x - (x.z)z 
			 (6) norm x 
			 
			 z-bisect
			 (1) norm z
			 (2) norm x 
			 (3) norm y 
			 (3) x = x + y
			 (4) norm x
			 (5) x = x - (x.z)z 
			 (6) norm x 
			 
			 3-fold
			 (1) norm z
			 (2) norm x 
			 (3) norm y 
			 (4) z = x + y + z
			 (5) norm z
			 (6) x = x - (x.z)z 
			 (7) norm x 
			 
			 */

			// branch based on axis type
			vectorZ = normalize(vectorZ);

			if (axisType == 1) {

				// bisector

				vectorX = normalize(vectorX);
				vectorZ += vectorX;
				vectorZ = normalize(vectorZ);
			}
			else if (axisType == 2 || axisType == 3) {

				// z-bisect

				if (particles.y >= 0 && particles.y < NUM_ATOMS) {
					real4 posY = posq[particles.y];
					real3 vectorY = make_real3(posY.x-thisParticlePos.x, posY.y-thisParticlePos.y, posY.z-thisParticlePos.z);
					vectorY = normalize(vectorY);
					vectorX = normalize(vectorX);
					if (axisType == 2) {
						vectorX += vectorY;
						vectorX = normalize(vectorX);
					}
					else {

						// 3-fold

						vectorZ += vectorX + vectorY;
						vectorZ = normalize(vectorZ);
					}
				}

			}
			else if (axisType >= 4)
			vectorX = make_real3((real) 0.1f);

			// x = x - (x.z)z

			vectorX -= dot(vectorZ, vectorX)*vectorZ;
			vectorX = normalize(vectorX);
			real3 vectorY = cross(vectorZ, vectorX);

			// use identity rotation matrix for unrecognized axis types

			if (axisType < 0 || axisType > 4) {

				vectorX.x = 1;
				vectorX.y = 0;
				vectorX.z = 0;

				vectorY.x = 0;
				vectorY.y = 1;
				vectorY.z = 0;

				vectorZ.x = 0;
				vectorZ.y = 0;
				vectorZ.z = 1;
			}

			// Check the chirality and see whether it needs to be reversed

			bool reverse = false;
			if (axisType != 0 && particles.x >= 0 && particles.y >=0 && particles.z >= 0) {
				real4 posY = posq[particles.y];
				real delta[4][3];

				delta[0][0] = thisParticlePos.x - posY.x;
				delta[0][1] = thisParticlePos.y - posY.y;
				delta[0][2] = thisParticlePos.z - posY.z;

				delta[1][0] = posZ.x - posY.x;
				delta[1][1] = posZ.y - posY.y;
				delta[1][2] = posZ.z - posY.z;

				delta[2][0] = posX.x - posY.x;
				delta[2][1] = posX.y - posY.y;
				delta[2][2] = posX.z - posY.z;

				delta[3][0] = delta[1][1]*delta[2][2] - delta[1][2]*delta[2][1];
				delta[3][1] = delta[2][1]*delta[0][2] - delta[2][2]*delta[0][1];
				delta[3][2] = delta[0][1]*delta[1][2] - delta[0][2]*delta[1][1];

				real volume = delta[3][0]*delta[0][0] + delta[3][1]*delta[1][0] + delta[3][2]*delta[2][0];
				reverse = (volume < 0);
			}

			// Transform the dipole

			unsigned int offset = 3*atom;
			real molDipole[3];
			molDipole[0] = molecularDipoles[offset];
			molDipole[1] = molecularDipoles[offset+1];
			molDipole[2] = molecularDipoles[offset+2];
			if (reverse)
			molDipole[1] *= -1;
			labFrameDipoles[offset] = molDipole[0]*vectorX.x + molDipole[1]*vectorY.x + molDipole[2]*vectorZ.x;
			labFrameDipoles[offset+1] = molDipole[0]*vectorX.y + molDipole[1]*vectorY.y + molDipole[2]*vectorZ.y;
			labFrameDipoles[offset+2] = molDipole[0]*vectorX.z + molDipole[1]*vectorY.z + molDipole[2]*vectorZ.z;
		}
		else {
			labFrameDipoles[3*atom] = molecularDipoles[3*atom];
			labFrameDipoles[3*atom+1] = molecularDipoles[3*atom+1];
			labFrameDipoles[3*atom+2] = molecularDipoles[3*atom+2];
		}
	}
}

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
extern "C" __global__ void computePotentialAtPoints(const real4* __restrict__ posq, const real* __restrict__ labFrameDipole,
		const real* __restrict__ labFrameQuadrupole, const real* __restrict__ inducedDipole, const real4* __restrict__ points,
		real* __restrict__ potential, int numPoints, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
	extern __shared__ real4 localPosq[];
	real3* localDipole = (real3*) &localPosq[blockDim.x];
	real3* localInducedDipole = (real3*) &localDipole[blockDim.x];
	real* localQuadrupole = (real*) &localInducedDipole[blockDim.x];
	for (int basePoint = blockIdx.x*blockDim.x; basePoint < numPoints; basePoint += gridDim.x*blockDim.x) {
		int point = basePoint+threadIdx.x;
		real4 pointPos = points[point];
		real p = 0;
		for (int baseAtom = 0; baseAtom < NUM_ATOMS; baseAtom += blockDim.x) {
			int atom = baseAtom+threadIdx.x;

			// Load data into shared memory.

			if (atom < NUM_ATOMS) {
				localPosq[threadIdx.x] = posq[atom];
				localDipole[threadIdx.x] = make_real3(labFrameDipole[3*atom], labFrameDipole[3*atom+1], labFrameDipole[3*atom+2]);
				localInducedDipole[threadIdx.x] = make_real3(inducedDipole[3*atom], inducedDipole[3*atom+1], inducedDipole[3*atom+2]);
				localQuadrupole[5*threadIdx.x] = labFrameQuadrupole[5*atom];
				localQuadrupole[5*threadIdx.x+1] = labFrameQuadrupole[5*atom+1];
				localQuadrupole[5*threadIdx.x+2] = labFrameQuadrupole[5*atom+2];
				localQuadrupole[5*threadIdx.x+3] = labFrameQuadrupole[5*atom+3];
				localQuadrupole[5*threadIdx.x+4] = labFrameQuadrupole[5*atom+4];
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
					real scd = dot(localDipole[i], delta);
					real scu = dot(localInducedDipole[i], delta);
					p -= (scd+scu)*rr3;
					real rr5 = 3*rr3*rr2;
					real scq = delta.x*dot(delta, make_real3(localQuadrupole[5*i+0], localQuadrupole[5*i+1], localQuadrupole[5*i+2])) +
					delta.y*dot(delta, make_real3(localQuadrupole[5*i+1], localQuadrupole[5*i+3], localQuadrupole[5*i+4])) +
					delta.z*dot(delta, make_real3(localQuadrupole[5*i+2], localQuadrupole[5*i+4], -localQuadrupole[5*i]-localQuadrupole[5*i+3]));
					p += scq*rr5;
				}
			}
			__syncthreads();
		}
		potential[point] = p;
	}
}
