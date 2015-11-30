#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

typedef struct {
    real4 posq;
    real3 field, fieldPolar, dipole;
    float damp;
} AtomData;

inline __device__ void loadAtomData(AtomData& data, int atom, const real4* __restrict__ posq, const real* __restrict__ labFrameDipole, const float* __restrict__ damping) {
    data.posq = posq[atom];
    data.dipole.x = labFrameDipole[atom*3];
    data.dipole.y = labFrameDipole[atom*3+1];
    data.dipole.z = labFrameDipole[atom*3+2];
    data.damp = damping[atom];
}

#ifdef USE_EWALD
__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real3 deltaR, float dScale, float pScale, real3* fields) {
    real r2 = dot(deltaR, deltaR);
    if (r2 <= CUTOFF_SQUARED) {
        // calculate the error function damping terms

        real r = SQRT(r2);
        real ralpha = EWALD_ALPHA*r;
        real exp2a = EXP(-(ralpha*ralpha));
#ifdef USE_DOUBLE_PRECISION
        const real erfcAlphaR = erfc(ralpha);
#else
        // This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.  They cite the following as
        // the original source: C. Hastings, Jr., Approximations for Digital Computers (1955).  It has a maximum
        // error of 1.5e-7.

        const real t = RECIP(1.0f+0.3275911f*ralpha);
        const real erfcAlphaR = (0.254829592f+(-0.284496736f+(1.421413741f+(-1.453152027f+1.061405429f*t)*t)*t)*t)*t*exp2a;
#endif
        real bn0 = erfcAlphaR/r;
        real alsq2 = 2*EWALD_ALPHA*EWALD_ALPHA;
        real alsq2n = RECIP(SQRT_PI*EWALD_ALPHA);
        alsq2n *= alsq2;
        real bn1 = (bn0+alsq2n*exp2a)/r2;
        alsq2n *= alsq2;
        real bn2 = (3*bn1+alsq2n*exp2a)/r2;
        alsq2n *= alsq2;
        real bn3 = (5*bn2+alsq2n*exp2a)/r2;

        // compute the error function scaled and unscaled terms

        real scale3 = 1;
        real damp = atom1.damp*atom2.damp;
        if (damp != 0) {
            real ratio = (r/damp);
            ratio = ratio*ratio*ratio;
            //FIXME: 
            // use minimum TCC (charge-charge) thole factor
            float thole = 0.4;
        	real pgamma = thole;
            damp = -pgamma*ratio;
            if (damp > -50) {
                real expdamp = EXP(damp);
                scale3 = 1 - expdamp;
            }
        }
        real dsc3 = dScale*scale3;

        real psc3 = pScale*scale3;

        real r3 = r*r2;
        real drr3 = (1-dsc3)/r3;

        real prr3 = (1-psc3)/r3;

        real3 fim = -deltaR*(bn1*atom2.posq.w);
        real3 fkm = deltaR*(bn1*atom1.posq.w);
        real3 fid = -deltaR*(drr3*atom2.posq.w);
        real3 fkd = deltaR*(drr3*atom1.posq.w);
        real3 fip = -deltaR*(prr3*atom2.posq.w);
        real3 fkp = deltaR*(prr3*atom1.posq.w);

        // increment the field at each site due to this interaction

        fields[0] = fim-fid;
        fields[1] = fim-fip;
        fields[2] = fkm-fkd;
        fields[3] = fkm-fkp;
    }
    else {
        fields[0] = make_real3(0);
        fields[1] = make_real3(0);
        fields[2] = make_real3(0);
        fields[3] = make_real3(0);
    }
}
#else
__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real3 deltaR, float dScale, float pScale, real3* fields) {
    real rI = RSQRT(dot(deltaR, deltaR));
    real r = RECIP(rI);
    real r2I = rI*rI;

    real rr3 = rI*r2I;
    real rr5 = 3*rr3*r2I;
    real rr7 = 5*rr5*r2I;
 
    // get scaling factors, if needed
    
    float damp = atom1.damp*atom2.damp;
    real dampExp;
    if (damp != 0) {

        // get scaling factors
      
        real ratio = r/damp;
        //FIXME: 
        float thole = 0.4;
        float pGamma = thole; 
        damp = ratio*ratio*ratio*pGamma;
        dampExp = EXP(-damp);
    }
    else
        dampExp = 0;
      
    rr3 *= 1 - dampExp;
    rr5 *= 1 - (1+damp)*dampExp;
    rr7 *= 1 - (1+damp+(0.6f*damp*damp))*dampExp;
       
    real dir = dot(atom1.dipole, deltaR);
    real dkr = dot(atom2.dipole, deltaR);

    real factor = -rr3*atom2.posq.w + rr5*dkr;
    real3 field1 = deltaR*factor - rr3*atom2.dipole;
    factor = rr3*atom1.posq.w + rr5*dir;
    real3 field2 = deltaR*factor - rr3*atom1.dipole;

    fields[0] = dScale*field1;
    fields[1] = pScale*field1;
    fields[2] = dScale*field2;
    fields[3] = pScale*field2;
}
#endif

__device__ real computeDScaleFactor(unsigned int polarizationGroup, int index) {
    return (polarizationGroup & 1<<index ? 0 : 1);
}

__device__ float computePScaleFactor(uint2 covalent, unsigned int polarizationGroup, int index) {
    int mask = 1<<index;
    bool x = (covalent.x & mask);
    bool y = (covalent.y & mask);
    bool p = (polarizationGroup & mask);
    return (x && y ? 0.0f : (x && p ? 0.5f : 1.0f));
}

/**
 * Compute nonbonded interactions.
 */
extern "C" __global__ void computeFixedField(
        unsigned long long* __restrict__ fieldBuffers, unsigned long long* __restrict__ fieldPolarBuffers, const real4* __restrict__ posq,
        const uint2* __restrict__ covalentFlags, const unsigned int* __restrict__ polarizationGroupFlags, const ushort2* __restrict__ exclusionTiles,
        unsigned int startTileIndex, unsigned int numTileIndices,
#ifdef USE_CUTOFF
        const int* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, const real4* __restrict__ blockCenter,
        const unsigned int* __restrict__ interactingAtoms,
#endif
        const real* __restrict__ labFrameDipole, const float* __restrict__ damping) {
    const unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    const unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE;
    const unsigned int tgx = threadIdx.x & (TILE_SIZE-1);
    const unsigned int tbx = threadIdx.x - tgx;
    __shared__ AtomData localData[THREAD_BLOCK_SIZE];

    // First loop: process tiles that contain exclusions.
    
    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const ushort2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        AtomData data;
        data.field = make_real3(0);
        data.fieldPolar = make_real3(0);
        unsigned int atom1 = x*TILE_SIZE + tgx;
        loadAtomData(data, atom1, posq, labFrameDipole, damping);
        uint2 covalent = covalentFlags[pos*TILE_SIZE+tgx];
        unsigned int polarizationGroup = polarizationGroupFlags[pos*TILE_SIZE+tgx];
        if (x == y) {
            // This tile is on the diagonal.

            const unsigned int localAtomIndex = threadIdx.x;
            localData[localAtomIndex].posq = data.posq;
            localData[localAtomIndex].dipole = data.dipole;
            localData[localAtomIndex].damp = data.damp;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                real3 delta = trimTo3(localData[tbx+j].posq-data.posq);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                int atom2 = y*TILE_SIZE+j;
                if (atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 fields[4];
                    float d = computeDScaleFactor(polarizationGroup, j);
                    float p = computePScaleFactor(covalent, polarizationGroup, j);
                    computeOneInteraction(data, localData[tbx+j], delta, d, p, fields);
                    data.field += fields[0];
                    data.fieldPolar += fields[1];
                }
            }
        }
        else {
            // This is an off-diagonal tile.

            const unsigned int localAtomIndex = threadIdx.x;
            unsigned int j = y*TILE_SIZE + tgx;
            loadAtomData(localData[localAtomIndex], j, posq, labFrameDipole, damping);
            localData[localAtomIndex].field = make_real3(0);
            localData[localAtomIndex].fieldPolar = make_real3(0);
            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                real3 delta = trimTo3(localData[tbx+tj].posq-data.posq);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                int atom2 = y*TILE_SIZE+tj;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 fields[4];
                    float d = computeDScaleFactor(polarizationGroup, tj);
                    float p = computePScaleFactor(covalent, polarizationGroup, tj);
                    computeOneInteraction(data, localData[tbx+tj], delta, d, p, fields);
                    data.field += fields[0];
                    data.fieldPolar += fields[1];
                    localData[tbx+tj].field += fields[2];
                    localData[tbx+tj].fieldPolar += fields[3];
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
        }
        
        // Write results.
        
        unsigned int offset = x*TILE_SIZE + tgx;
        atomicAdd(&fieldBuffers[offset], static_cast<unsigned long long>((long long) (data.field.x*0x100000000)));
        atomicAdd(&fieldBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.y*0x100000000)));
        atomicAdd(&fieldBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.z*0x100000000)));
        atomicAdd(&fieldPolarBuffers[offset], static_cast<unsigned long long>((long long) (data.fieldPolar.x*0x100000000)));
        atomicAdd(&fieldPolarBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.y*0x100000000)));
        atomicAdd(&fieldPolarBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.z*0x100000000)));
        if (x != y) {
            offset = y*TILE_SIZE + tgx;
            atomicAdd(&fieldBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.x*0x100000000)));
            atomicAdd(&fieldBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.y*0x100000000)));
            atomicAdd(&fieldBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.z*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.x*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.y*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.z*0x100000000)));
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
            data.field = make_real3(0);
            data.fieldPolar = make_real3(0);
            loadAtomData(data, atom1, posq, labFrameDipole, damping);
#ifdef USE_CUTOFF
            unsigned int j = (numTiles <= maxTiles ? interactingAtoms[pos*TILE_SIZE+tgx] : y*TILE_SIZE + tgx);
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[threadIdx.x] = j;
            const unsigned int localAtomIndex = threadIdx.x;
            loadAtomData(localData[localAtomIndex], j, posq, labFrameDipole, damping);
            localData[localAtomIndex].field = make_real3(0);
            localData[localAtomIndex].fieldPolar = make_real3(0);

            // Compute the full set of interactions in this tile.

            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                real3 delta = trimTo3(localData[tbx+tj].posq-data.posq);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                int atom2 = atomIndices[tbx+tj];
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS) {
                    real3 fields[4];
                    computeOneInteraction(data, localData[tbx+tj], delta, 1, 1, fields);
                    data.field += fields[0];
                    data.fieldPolar += fields[1];
                    localData[tbx+tj].field += fields[2];
                    localData[tbx+tj].fieldPolar += fields[3];
                }
                tj = (tj + 1) & (TILE_SIZE - 1);
            }

            // Write results.

            unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&fieldBuffers[offset], static_cast<unsigned long long>((long long) (data.field.x*0x100000000)));
            atomicAdd(&fieldBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.y*0x100000000)));
            atomicAdd(&fieldBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.z*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset], static_cast<unsigned long long>((long long) (data.fieldPolar.x*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.y*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.z*0x100000000)));
#ifdef USE_CUTOFF
            offset = atomIndices[threadIdx.x];
#else
            offset = y*TILE_SIZE + tgx;
#endif
            atomicAdd(&fieldBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.x*0x100000000)));
            atomicAdd(&fieldBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.y*0x100000000)));
            atomicAdd(&fieldBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.z*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.x*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.y*0x100000000)));
            atomicAdd(&fieldPolarBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.z*0x100000000)));
        }
        pos++;
    }
}
