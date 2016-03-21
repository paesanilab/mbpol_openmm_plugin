#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

typedef struct {
    real3 pos, force, dipole, inducedDipole, inducedDipolePolar;
    real q;
    float damp;
    real potential;
    int moleculeIndex;
    int atomType;
} AtomData;

__device__ void computeOneInteractionF1(AtomData& atom1, volatile AtomData& atom2, real4 delta, real4 bn, real bn5, float forceFactor, real3& force, real& energy, real2& potential);
__device__ void computeOneInteractionF2(AtomData& atom1, volatile AtomData& atom2, real4 delta, real4 bn, float forceFactor, real3& force, real& energy, real2& potential);

inline __device__ void loadAtomData(AtomData& data, int atom, const real4* __restrict__ posq,
        const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar,
        const float* __restrict__ damping, const int* __restrict__ moleculeIndex, const int* __restrict__ atomType) {
    real4 atomPosq = posq[atom];
    data.pos = make_real3(atomPosq.x, atomPosq.y, atomPosq.z);
    data.q = atomPosq.w;
    data.inducedDipole.x = inducedDipole[atom*3];
    data.inducedDipole.y = inducedDipole[atom*3+1];
    data.inducedDipole.z = inducedDipole[atom*3+2];
    data.inducedDipolePolar.x = inducedDipolePolar[atom*3];
    data.inducedDipolePolar.y = inducedDipolePolar[atom*3+1];
    data.inducedDipolePolar.z = inducedDipolePolar[atom*3+2];
    data.damp = damping[atom];
    data.moleculeIndex = moleculeIndex[atom];
    data.atomType = atomType[atom];
}

__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, bool hasExclusions, float forceFactor,
                                      real& energy, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    real4 delta;
    delta.x = atom2.pos.x - atom1.pos.x;
    delta.y = atom2.pos.y - atom1.pos.y;
    delta.z = atom2.pos.z - atom1.pos.z;
    real energy_before = energy;

    // periodic box

    APPLY_PERIODIC_TO_DELTA(delta)
    delta.w = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
    if (delta.w > CUTOFF_SQUARED)
        return;

    real r = SQRT(delta.w);
    real ralpha = EWALD_ALPHA*r;

    real alsq2 = 2*EWALD_ALPHA*EWALD_ALPHA;
    real alsq2n = 0;
    if (EWALD_ALPHA > 0)
        alsq2n = RECIP(SQRT_PI*EWALD_ALPHA);
    real exp2a = EXP(-(ralpha*ralpha));

    real rr1 = RECIP(r);
    delta.w = rr1;
#ifdef USE_DOUBLE_PRECISION
    const real erfcAlphaR = erfc(ralpha);
#else
    // This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.  They cite the following as
    // the original source: C. Hastings, Jr., Approximations for Digital Computers (1955).  It has a maximum
    // error of 1.5e-7.

    const real t = RECIP(1.0f+0.3275911f*ralpha);
    const real erfcAlphaR = (0.254829592f+(-0.284496736f+(1.421413741f+(-1.453152027f+1.061405429f*t)*t)*t)*t)*t*exp2a;
#endif
    real bn0 = erfcAlphaR*rr1;
    energy += forceFactor*atom1.q*atom2.q*bn0;
    real rr2 = rr1*rr1;
    alsq2n *= alsq2;

    real4 bn;
    bn.x = (bn0+alsq2n*exp2a)*rr2;

    alsq2n *= alsq2;
    bn.y = (3*bn.x+alsq2n*exp2a)*rr2;

    alsq2n *= alsq2;
    bn.z = (5*bn.y+alsq2n*exp2a)*rr2;

    bn.w = bn0; // bn4 is not used, so we store bn0

    alsq2n *= alsq2;
    real bn5 = (9*bn.w+alsq2n*exp2a)*rr2;

    real3 force;
    real2 potential = make_real2(0., 0.);

    computeOneInteractionF1(atom1, atom2, delta, bn, bn5, forceFactor, force, energy, potential);
    computeOneInteractionF2(atom1, atom2, delta, bn, forceFactor, force, energy, potential);

    atom1.force += force;
    atom1.potential += potential.x;
    if (forceFactor == 1) {
        atom2.force -= force;
        atom2.potential += potential.y;
    }
}

/**
 * Compute the self energy and self torque.
 */
__device__ void computeSelfEnergyAndTorque(AtomData& atom1, real& energy) {
    real fterm = -EWALD_ALPHA/SQRT_PI;
    real cii = atom1.q*atom1.q;
    real selfEnergy = cii;
    selfEnergy *= fterm;

    energy += selfEnergy;
    #ifdef INCLUDE_CHARGE_REDISTRIBUTION
        atom1.potential += atom1.q * fterm * 2;
    #endif

}

/**
 * Compute electrostatic interactions.
 */
extern "C" __global__ void computeElectrostatics(
        unsigned long long* __restrict__ forceBuffers, unsigned long long* __restrict__ potentialBuffers,real* __restrict__ energyBuffer,
        const real4* __restrict__ posq,
        const ushort2* __restrict__ exclusionTiles, unsigned int startTileIndex, unsigned int numTileIndices,
#ifdef USE_CUTOFF
        const int* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, const real4* __restrict__ blockCenter,
        const unsigned int* __restrict__ interactingAtoms,
#endif
        const real* __restrict__ inducedDipole,
        const real* __restrict__ inducedDipolePolar,
        const float* __restrict__ damping, const int* __restrict__ moleculeIndex, const int* __restrict__ atomType) {

    const unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    const unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE;
    const unsigned int tgx = threadIdx.x & (TILE_SIZE-1);
    const unsigned int tbx = threadIdx.x - tgx;
    real energy = 0;
    __shared__ AtomData localData[THREAD_BLOCK_SIZE];

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        AtomData data;
        unsigned int atom1 = x*TILE_SIZE + tgx;
        loadAtomData(data, atom1, posq, inducedDipole, inducedDipolePolar, damping, moleculeIndex, atomType);
        data.force = make_real3(0);
        data.potential = 0.;
        if (x == y) {
            // This tile is on the diagonal.

            localData[threadIdx.x].pos = data.pos;
            localData[threadIdx.x].q = data.q;
            localData[threadIdx.x].inducedDipole = data.inducedDipole;
            localData[threadIdx.x].inducedDipolePolar = data.inducedDipolePolar;
            localData[threadIdx.x].damp = data.damp;
            localData[threadIdx.x].moleculeIndex = data.moleculeIndex;
            localData[threadIdx.x].atomType = data.atomType;

            // Compute forces.

            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+j;
                if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    float d = 1.;
                    float p = 1.;
                    float m = 1.;
                    computeOneInteraction(data, localData[tbx+j], true, 0.5f, energy, periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
                }
            }
            if (atom1 < NUM_ATOMS)
                computeSelfEnergyAndTorque(data, energy);
            data.force *= -ENERGY_SCALE_FACTOR;
            atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
            atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
            atomicAdd(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
            atomicAdd(&potentialBuffers[atom1], static_cast<unsigned long long>((long long) (data.potential*0x100000000)));
        }
        else {
            // This is an off-diagonal tile.

            unsigned int j = y*TILE_SIZE + tgx;
            loadAtomData(localData[threadIdx.x], j, posq, inducedDipole, inducedDipolePolar, damping, moleculeIndex, atomType);
            localData[threadIdx.x].force = make_real3(0);
            localData[threadIdx.x].potential = 0.;
            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = y*TILE_SIZE+tj;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    float d = 1.;
                    float p = 1.;
                    float m = 1.;
                    computeOneInteraction(data, localData[tbx+tj], true, 1, energy, periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
            data.force *= -ENERGY_SCALE_FACTOR;
            localData[threadIdx.x].force *= -ENERGY_SCALE_FACTOR;
            unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
            atomicAdd(&potentialBuffers[offset], static_cast<unsigned long long>((long long) (data.potential*0x100000000)));
            offset = y*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));
            atomicAdd(&potentialBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].potential*0x100000000)));
        }
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_CUTOFF
    const unsigned int numTiles = interactionCount[0];
    int pos = (int) (numTiles > maxTiles ? startTileIndex+warp*(long long)numTileIndices/totalWarps : warp*(long long)numTiles/totalWarps);
    int end = (int) (numTiles > maxTiles ? startTileIndex+(warp+1)*(long long)numTileIndices/totalWarps : (warp+1)*(long long)numTiles/totalWarps);
#else
    const unsigned int numTiles = numTileIndices;
    int pos = (int) (startTileIndex+warp*(long long)numTiles/totalWarps);
    int end = (int) (startTileIndex+(warp+1)*(long long)numTiles/totalWarps);
#endif
    int skipBase = 0;
    int currentSkipIndex = tbx;
    __shared__ int atomIndices[THREAD_BLOCK_SIZE];
    __shared__ volatile int skipTiles[THREAD_BLOCK_SIZE];
    skipTiles[threadIdx.x] = -1;
    
    while (pos < end) {
        bool includeTile = true;

        // Extract the coordinates of this tile.
        
        int x, y;
#ifdef USE_CUTOFF
        if (numTiles <= maxTiles)
            x = tiles[pos];
        else
#endif
        {
            y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
                y += (x < y ? -1 : 1);
                x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
            }

            // Skip over tiles that have exclusions, since they were already processed.

            while (skipTiles[tbx+TILE_SIZE-1] < pos) {
                if (skipBase+tgx < NUM_TILES_WITH_EXCLUSIONS) {
                    ushort2 tile = exclusionTiles[skipBase+tgx];
                    skipTiles[threadIdx.x] = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
                }
                else
                    skipTiles[threadIdx.x] = end;
                skipBase += TILE_SIZE;            
                currentSkipIndex = tbx;
            }
            while (skipTiles[currentSkipIndex] < pos)
                currentSkipIndex++;
            includeTile = (skipTiles[currentSkipIndex] != pos);
        }
        if (includeTile) {
            unsigned int atom1 = x*TILE_SIZE + tgx;

            // Load atom data for this tile.

            AtomData data;
            loadAtomData(data, atom1, posq, inducedDipole, inducedDipolePolar, damping, moleculeIndex, atomType);
            data.force = make_real3(0);
            data.potential = 0.;
#ifdef USE_CUTOFF
            unsigned int j = (numTiles <= maxTiles ? interactingAtoms[pos*TILE_SIZE+tgx] : y*TILE_SIZE + tgx);
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[threadIdx.x] = j;
            loadAtomData(localData[threadIdx.x], j, posq, inducedDipole, inducedDipolePolar, damping, moleculeIndex, atomType);
            localData[threadIdx.x].force = make_real3(0);
            localData[threadIdx.x].potential = 0.;

            // Compute forces.

            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = atomIndices[tbx+tj];
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    computeOneInteraction(data, localData[tbx+tj], false, 1, energy, periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
            data.force *= -ENERGY_SCALE_FACTOR;
            localData[threadIdx.x].force *= -ENERGY_SCALE_FACTOR;

            // Write results.

            unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (data.force.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.force.z*0x100000000)));
            atomicAdd(&potentialBuffers[offset], static_cast<unsigned long long>((long long) (data.potential*0x100000000)));
#ifdef USE_CUTOFF
            offset = atomIndices[threadIdx.x];
#else
            offset = y*TILE_SIZE + tgx;
#endif
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.x*0x100000000)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.y*0x100000000)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].force.z*0x100000000)));
            atomicAdd(&potentialBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].potential*0x100000000)));
        }
        pos++;
    }
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy*ENERGY_SCALE_FACTOR;
}
