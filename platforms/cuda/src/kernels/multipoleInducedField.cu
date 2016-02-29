#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

typedef struct {
    real3 pos;
    real3 field, fieldPolar, inducedDipole, inducedDipolePolar;
    float damp;
    int moleculeIndex;
    int atomType;
} AtomData;

inline __device__ void loadAtomData(AtomData& data, int atom, const real4* __restrict__ posq, const real* __restrict__ inducedDipole,
        const real* __restrict__ inducedDipolePolar, const float* __restrict__ damping, const int* __restrict__ moleculeIndex, const int* __restrict__ atomType) {
    real4 atomPosq = posq[atom];
    data.pos = make_real3(atomPosq.x, atomPosq.y, atomPosq.z);
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

inline __device__ void zeroAtomData(AtomData& data) {
    data.field = make_real3(0);
    data.fieldPolar = make_real3(0);
}

#ifdef USE_EWALD
__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real3 deltaR, bool isSelfInteraction) {
    if (isSelfInteraction)
        return;
    real scale1, scale2;
    real r2 = dot(deltaR, deltaR);
    bool isSameWater = atom1.moleculeIndex == atom2.moleculeIndex;
    if (r2 < CUTOFF_SQUARED) {
        real rI = RSQRT(r2);
        real r = RECIP(rI);

        // calculate the error function damping terms

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
        // FIXME thole copy in unique location
        const enum TholeIndices { TCC, TCD, TDD, TDDOH, TDDHH };
        const float thole[5] =  { 0.4, 0.4, 0.4,   0.4,   0.4 };
        real bn0 = erfcAlphaR*rI;
        real alsq2 = 2*EWALD_ALPHA*EWALD_ALPHA;
        real alsq2n = RECIP(SQRT_PI*EWALD_ALPHA);
        alsq2n *= alsq2;
        real bn1 = (bn0+alsq2n*exp2a)*rI*rI;

        alsq2n *= alsq2;
        real bn2 = (3*bn1+alsq2n*exp2a)*rI*rI;

        // compute the error function scaled and unscaled terms

        //RealOpenMM scale3 = getAndScaleInverseRs(particleI, particleJ, r, true, 3, TDD);
        //RealOpenMM scale5 = getAndScaleInverseRs(particleI, particleJ, r, true, 5, TDD);
        real damp      = POW(atom1.damp*atom2.damp, 1.0f/6.0f); // AA in MBPol

        bool do_scaling = (damp != 0.0) & ( damp > -50.0 ); // damp or not

        real ratio       = POW(r/damp, 4); // rA4 in MBPol

        // FIXME identify if we need to use TDDOH and so on
        int tdd = TDD;
        if ((isSameWater) && (atom1.atomType != 2) && (atom2.atomType != 2)) {
            if ((atom1.atomType == 0) | (atom2.atomType == 0)) { // one is oxygen
                tdd = TDDOH;
            } else { // both hydrogens
                tdd = TDDHH;
            }
        }
        real pgamma = thole[tdd];
        real dampForExp = -1 * pgamma * ratio;

        real scale3 = 1.0;
        if (do_scaling)
            scale3 -= EXP(dampForExp);
        real scale5 = scale3;
        if (do_scaling)
            scale5 -= (4./3.) * pgamma * EXP(dampForExp) * ratio;

        real r3 = (r*r2);
        real r5 = (r3*r2);
        real rr3 = (1-scale3)/r3;
        real rr5 = 3*(1-scale5)/r5;

        scale1 = rr3 - bn1;
        scale2 = bn2 - rr5;
    }
    else {
        scale1 = 0;
        scale2 = 0;
    }
    real dDotDelta = scale2*dot(deltaR, atom2.inducedDipole);
    atom1.field += scale1*atom2.inducedDipole + dDotDelta*deltaR;
    dDotDelta = scale2*dot(deltaR, atom2.inducedDipolePolar);
    atom1.fieldPolar += scale1*atom2.inducedDipolePolar + dDotDelta*deltaR;
    dDotDelta = scale2*dot(deltaR, atom1.inducedDipole);
    atom2.field += scale1*atom1.inducedDipole + dDotDelta*deltaR;
    dDotDelta = scale2*dot(deltaR, atom1.inducedDipolePolar);
    atom2.fieldPolar += scale1*atom1.inducedDipolePolar + dDotDelta*deltaR;
}
#else
__device__ void computeOneInteraction(AtomData& atom1, AtomData& atom2, real3 deltaR, bool isSelfInteraction) {
    if (isSelfInteraction)
        return;
    // FIXME thole copy in unique location
    const enum TholeIndices { TCC, TCD, TDD, TDDOH, TDDHH };
    const float thole[5] =  { 0.4, 0.4, 0.4,   0.4,   0.4 };

    // RealOpenMM scale3 = getAndScaleInverseRs( particleI, particleJ, r, false, 3, TDD);
    // RealOpenMM scale5 = getAndScaleInverseRs( particleI, particleJ, r, false, 5, TDD);

    real rI = RSQRT(dot(deltaR, deltaR));
    real r = RECIP(rI);
    real r2I = rI*rI;
    real rr3 = -rI*r2I;
    real rr5 = -3*rr3*r2I;

    real damp      = pow(atom1.damp*atom2.damp, 1.0f/6.0f); // AA in MBPol

    bool do_scaling = (damp != 0.0) & ( damp > -50.0 ); // damp or not

    real ratio       = pow(r/damp, 4); // rA4 in MBPol

    int tdd = TDD;
    bool isSameWater = atom1.moleculeIndex == atom2.moleculeIndex;
    if ((isSameWater) && (atom1.atomType != 2) && (atom2.atomType != 2)) {
        if ((atom1.atomType == 0) | (atom2.atomType == 0)) { // one is oxygen
            tdd = TDDOH;
        } else { // both hydrogens
            tdd = TDDHH;
        }
    }
    real pgamma = thole[TDD];
    real dampForExp = -1 * pgamma * ratio;

    real rr3_factor = 1.0;
    if (do_scaling)
        rr3_factor -= EXP(dampForExp);
    rr3 *= rr3_factor;

    real rr5_factor = rr3_factor;
    if (do_scaling)
        rr5_factor -= (4./3.) * pgamma * EXP(dampForExp) * ratio;
    rr5 *= rr5_factor;

    real dDotDelta = rr5*dot(deltaR, atom2.inducedDipole);
    atom1.field += rr3*atom2.inducedDipole + dDotDelta*deltaR;

    dDotDelta = rr5*dot(deltaR, atom2.inducedDipolePolar);
    atom1.fieldPolar += rr3*atom2.inducedDipolePolar + dDotDelta*deltaR;
    dDotDelta = rr5*dot(deltaR, atom1.inducedDipole);
    atom2.field += rr3*atom1.inducedDipole + dDotDelta*deltaR;
    dDotDelta = rr5*dot(deltaR, atom1.inducedDipolePolar);
    atom2.fieldPolar += rr3*atom1.inducedDipolePolar + dDotDelta*deltaR;
}
#endif

/**
 * Compute the mutual induced field.
 */
extern "C" __global__ void computeInducedField(
        unsigned long long* __restrict__ field, unsigned long long* __restrict__ fieldPolar, const real4* __restrict__ posq, const ushort2* __restrict__ exclusionTiles, 
        const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar, unsigned int startTileIndex, unsigned int numTileIndices,
#ifdef USE_CUTOFF
        const int* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize,
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, const real4* __restrict__ blockCenter, const unsigned int* __restrict__ interactingAtoms,
#endif
        const float* __restrict__ damping, const int* __restrict__ moleculeIndex, const int* __restrict__ atomType) {
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
        zeroAtomData(data);
        unsigned int atom1 = x*TILE_SIZE + tgx;
        loadAtomData(data, atom1, posq, inducedDipole, inducedDipolePolar, damping, moleculeIndex, atomType);
        if (x == y) {
            // This tile is on the diagonal.

            localData[threadIdx.x].pos = data.pos;
            localData[threadIdx.x].inducedDipole = data.inducedDipole;
            localData[threadIdx.x].inducedDipolePolar = data.inducedDipolePolar;
            localData[threadIdx.x].damp = data.damp;
            localData[threadIdx.x].moleculeIndex = data.moleculeIndex;
            localData[threadIdx.x].atomType = data.atomType;

            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                real3 delta = localData[tbx+j].pos-data.pos;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                int atom2 = y*TILE_SIZE+j;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
                    computeOneInteraction(data, localData[tbx+j], delta, atom1 == atom2);
            }
        }
        else {
            // This is an off-diagonal tile.

            loadAtomData(localData[threadIdx.x], y*TILE_SIZE+tgx, posq, inducedDipole, inducedDipolePolar, damping, moleculeIndex, atomType);

            zeroAtomData(localData[threadIdx.x]);
            unsigned int tj = tgx;
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                real3 delta = localData[tbx+tj].pos-data.pos;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                int atom2 = y*TILE_SIZE+j;
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
                    computeOneInteraction(data, localData[tbx+tj], delta, false);
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
        }

        // Write results.

        unsigned int offset = x*TILE_SIZE + tgx;
        atomicAdd(&field[offset], static_cast<unsigned long long>((long long) (data.field.x*0x100000000)));
        atomicAdd(&field[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.y*0x100000000)));
        atomicAdd(&field[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.z*0x100000000)));
        atomicAdd(&fieldPolar[offset], static_cast<unsigned long long>((long long) (data.fieldPolar.x*0x100000000)));
        atomicAdd(&fieldPolar[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.y*0x100000000)));
        atomicAdd(&fieldPolar[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.z*0x100000000)));

        if (x != y) {
            offset = y*TILE_SIZE + tgx;
            atomicAdd(&field[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.x*0x100000000)));
            atomicAdd(&field[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.y*0x100000000)));
            atomicAdd(&field[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.z*0x100000000)));
            atomicAdd(&fieldPolar[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.x*0x100000000)));
            atomicAdd(&fieldPolar[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.y*0x100000000)));
            atomicAdd(&fieldPolar[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.z*0x100000000)));

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
            zeroAtomData(data);

            loadAtomData(data, atom1, posq, inducedDipole, inducedDipolePolar, damping, moleculeIndex, atomType);

#ifdef USE_CUTOFF
            unsigned int j = (numTiles <= maxTiles ? interactingAtoms[pos*TILE_SIZE+tgx] : y*TILE_SIZE + tgx);
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[threadIdx.x] = j;

            loadAtomData(localData[threadIdx.x], j, posq, inducedDipole, inducedDipolePolar, damping, moleculeIndex, atomType);

            zeroAtomData(localData[threadIdx.x]);

            // Compute the full set of interactions in this tile.

            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                real3 delta = localData[tbx+tj].pos-data.pos;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                int atom2 = atomIndices[tbx+tj];
                if (atom1 < NUM_ATOMS && atom2 < NUM_ATOMS)
                    computeOneInteraction(data, localData[tbx+tj], delta, false);
                tj = (tj + 1) & (TILE_SIZE - 1);
            }

            // Write results.

            unsigned int offset = x*TILE_SIZE + tgx;
            atomicAdd(&field[offset], static_cast<unsigned long long>((long long) (data.field.x*0x100000000)));
            atomicAdd(&field[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.y*0x100000000)));
            atomicAdd(&field[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.field.z*0x100000000)));
            atomicAdd(&fieldPolar[offset], static_cast<unsigned long long>((long long) (data.fieldPolar.x*0x100000000)));
            atomicAdd(&fieldPolar[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.y*0x100000000)));
            atomicAdd(&fieldPolar[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (data.fieldPolar.z*0x100000000)));

#ifdef USE_CUTOFF
            offset = atomIndices[threadIdx.x];
#else
            offset = y*TILE_SIZE + tgx;
#endif
            atomicAdd(&field[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.x*0x100000000)));
            atomicAdd(&field[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.y*0x100000000)));
            atomicAdd(&field[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].field.z*0x100000000)));
            atomicAdd(&fieldPolar[offset], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.x*0x100000000)));
            atomicAdd(&fieldPolar[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.y*0x100000000)));
            atomicAdd(&fieldPolar[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (localData[threadIdx.x].fieldPolar.z*0x100000000)));

        }
        pos++;
    }
}

extern "C" __global__ void recordInducedDipolesForDIIS(const long long* __restrict__ fixedField, const long long* __restrict__ fixedFieldPolar,
        const long long* __restrict__ fixedFieldS, const long long* __restrict__ inducedField, const long long* __restrict__ inducedFieldPolar,
        const real* __restrict__ inducedDipole, const real* __restrict__ inducedDipolePolar, const float* __restrict__ polarizability, float2* __restrict__ errors,
        real* __restrict__ prevDipoles, real* __restrict__ prevDipolesPolar, real* __restrict__ prevErrors, int iteration, bool recordPrevErrors, real* __restrict__ matrix) {
    extern __shared__ real2 buffer[];
#ifdef USE_EWALD
    const real ewaldScale = (4/(real) 3)*(EWALD_ALPHA*EWALD_ALPHA*EWALD_ALPHA)/SQRT_PI;
#else
    const real ewaldScale = 0;
#endif
    const real fieldScale = 1/(real) 0x100000000;
    real sumErrors = 0;
    real sumPolarErrors = 0;
    for (int atom = blockIdx.x*blockDim.x + threadIdx.x; atom < NUM_ATOMS; atom += blockDim.x*gridDim.x) {
        real scale = polarizability[atom];
        for (int component = 0; component < 3; component++) {
            int dipoleIndex = 3*atom+component;
            int fieldIndex = atom+component*PADDED_NUM_ATOMS;
            if (iteration >= MAX_PREV_DIIS_DIPOLES) {
                // We have filled up the buffer for previous dipoles, so shift them all over by one.
                
                for (int i = 1; i < MAX_PREV_DIIS_DIPOLES; i++) {
                    int index1 = dipoleIndex+(i-1)*NUM_ATOMS*3;
                    int index2 = dipoleIndex+i*NUM_ATOMS*3;
                    prevDipoles[index1] = prevDipoles[index2];
                    prevDipolesPolar[index1] = prevDipolesPolar[index2];
                    if (recordPrevErrors)
                        prevErrors[index1] = prevErrors[index2];
                }
            }
            
            // Compute the new dipole, and record it along with the error.
            
            real oldDipole = inducedDipole[dipoleIndex];
            real oldDipolePolar = inducedDipolePolar[dipoleIndex];
            long long fixedS = (fixedFieldS == NULL ? (long long) 0 : fixedFieldS[fieldIndex]);
            real newDipole = scale*((fixedField[fieldIndex]+fixedS+inducedField[fieldIndex])*fieldScale+ewaldScale*oldDipole);
            real newDipolePolar = scale*((fixedFieldPolar[fieldIndex]+fixedS+inducedFieldPolar[fieldIndex])*fieldScale+ewaldScale*oldDipolePolar);
            int storePrevIndex = dipoleIndex+min(iteration, MAX_PREV_DIIS_DIPOLES-1)*NUM_ATOMS*3;
            prevDipoles[storePrevIndex] = newDipole;
            prevDipolesPolar[storePrevIndex] = newDipolePolar;
            if (recordPrevErrors)
                prevErrors[storePrevIndex] = newDipole-oldDipole;
            sumErrors += (newDipole-oldDipole)*(newDipole-oldDipole);
            sumPolarErrors += (newDipolePolar-oldDipolePolar)*(newDipolePolar-oldDipolePolar);
        }
    }
    
    // Sum the errors over threads and store the total for this block.
    
    buffer[threadIdx.x] = make_real2(sumErrors, sumPolarErrors);
    __syncthreads();
    for (int offset = 1; offset < blockDim.x; offset *= 2) {
        if (threadIdx.x+offset < blockDim.x && (threadIdx.x&(2*offset-1)) == 0) {
            buffer[threadIdx.x].x += buffer[threadIdx.x+offset].x;
            buffer[threadIdx.x].y += buffer[threadIdx.x+offset].y;
        }
        __syncthreads();
    }
    if (threadIdx.x == 0)
        errors[blockIdx.x] = make_float2((float) buffer[0].x, (float) buffer[0].y);
    
    if (iteration >= MAX_PREV_DIIS_DIPOLES && recordPrevErrors && blockIdx.x == 0) {
        // Shift over the existing matrix elements.
        
        for (int i = 0; i < MAX_PREV_DIIS_DIPOLES-1; i++) {
            if (threadIdx.x < MAX_PREV_DIIS_DIPOLES-1)
                matrix[threadIdx.x+i*MAX_PREV_DIIS_DIPOLES] = matrix[(threadIdx.x+1)+(i+1)*MAX_PREV_DIIS_DIPOLES];
            __syncthreads();
        }
    }
}

extern "C" __global__ void computeDIISMatrix(real* __restrict__ prevErrors, int iteration, real* __restrict__ matrix) {
    extern __shared__ real sumBuffer[];
    int j = min(iteration, MAX_PREV_DIIS_DIPOLES-1);
    for (int i = blockIdx.x; i <= j; i += gridDim.x) {
        // All the threads in this thread block work together to compute a single matrix element.

        real sum = 0;
        for (int index = threadIdx.x; index < NUM_ATOMS*3; index += blockDim.x)
            sum += prevErrors[index+i*NUM_ATOMS*3]*prevErrors[index+j*NUM_ATOMS*3];
        sumBuffer[threadIdx.x] = sum;
        __syncthreads();
        for (int offset = 1; offset < blockDim.x; offset *= 2) { 
            if (threadIdx.x+offset < blockDim.x && (threadIdx.x&(2*offset-1)) == 0)
                sumBuffer[threadIdx.x] += sumBuffer[threadIdx.x+offset];
            __syncthreads();
        }
        if (threadIdx.x == 0) {
            matrix[i+MAX_PREV_DIIS_DIPOLES*j] = sumBuffer[0];
            if (i != j)
                matrix[j+MAX_PREV_DIIS_DIPOLES*i] = sumBuffer[0];
        }
    }
}

extern "C" __global__ void updateInducedFieldByDIIS(real* __restrict__ inducedDipole, real* __restrict__ inducedDipolePolar, 
        const real* __restrict__ prevDipoles, const real* __restrict__ prevDipolesPolar, const float* __restrict__ coefficients, int numPrev) {
    for (int index = blockIdx.x*blockDim.x + threadIdx.x; index < 3*NUM_ATOMS; index += blockDim.x*gridDim.x) {
        real sum = 0;
        real sumPolar = 0;
        for (int i = 0; i < numPrev; i++) {
            sum += coefficients[i]*prevDipoles[i*3*NUM_ATOMS+index];
            sumPolar += coefficients[i]*prevDipolesPolar[i*3*NUM_ATOMS+index];
        }
        inducedDipole[index] = sum;
        inducedDipolePolar[index] = sumPolar;
    }
}
