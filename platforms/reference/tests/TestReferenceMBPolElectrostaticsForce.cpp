/* -------------------------------------------------------------------------- *
 *                                   OpenMMMBPol                             *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,  *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/**
 * This tests the Reference implementation of ReferenceMBPolElectrostaticsForce.
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMMBPol.h"
#include "openmm/System.h"
#include "openmm/MBPolElectrostaticsForce.h"
#include "MBPolReferenceElectrostaticsForce.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/VirtualSite.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#define ASSERT_EQUAL_TOL_MOD(expected, found, tol, testname) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC_MOD(expected, found, tol,testname) {ASSERT_EQUAL_TOL_MOD((expected)[0], (found)[0], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[1], (found)[1], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[2], (found)[2], (tol),(testname));};

// #define COMPUTE_FINITE_DIFFERENCES_FORCES

using namespace  OpenMM;
using namespace MBPolPlugin;
using namespace std;
const double TOL = 1e-4;
const double cal2joule = 4.184;

// getAndScaleInverseRs is protected, so we need to create a wrapping class for testing it.
class WrappedMBPolReferenceElectrostaticsForce : public MBPolReferenceElectrostaticsForce {
    public:
    void wrapGetAndScaleInverseRs(
            RealOpenMM dampI, RealOpenMM dampJ,
            RealOpenMM tholeI, RealOpenMM tholeJ,
            RealOpenMM r, bool justScale, RealOpenMM & damp, MapIntRealOpenMM& rrI
            )   {

                    std::vector<ElectrostaticsParticleData> particleData;
                    particleData.resize(2);
                    particleData[0].dampingFactor = dampI;
                    particleData[1].dampingFactor = dampJ;
                    for (int i=0; i<5; i++) {
                        particleData[0].thole[i] = tholeI;
                        particleData[1].thole[i] = tholeJ;
                    }

                    for (int order=1; order <=7; order+=2) {
                        rrI[order] = getAndScaleInverseRs(particleData[0], particleData[1], r, justScale, order, TCC);
                    }
                }
};

static void testGetAndScaleInverseRs( ) {

    std::string testName      = "testGetAndScaleInverseRs";

    RealOpenMM damp=10.;
    //RealOpenMM dampO=0.306988;
    //RealOpenMM dampH=0.28135;
    RealOpenMM dampO=0.001310;
    RealOpenMM dampH=0.000294;
    MapIntRealOpenMM rrI;
    RealOpenMM r=9.860634018e-02; // from Water3 test
    RealOpenMM thole=0.400;


    WrappedMBPolReferenceElectrostaticsForce* mbpolReferenceElectrostaticsForce = new WrappedMBPolReferenceElectrostaticsForce();;
    mbpolReferenceElectrostaticsForce->wrapGetAndScaleInverseRs( dampO, dampH,
                          thole, thole, r, false, damp, rrI);

    ASSERT_EQUAL_TOL_MOD(9.33047, rrI[1], 1e-5, testName); // from this plugin after integration testing with mbpol on water3
    ASSERT_EQUAL_TOL_MOD(5.324612470e+02, rrI[3], 1e-5, testName); // from mbpol
    // mbpol multiplies by constant factor (3) later, MBPOL in this function
    ASSERT_EQUAL_TOL_MOD(4.747626558e+03*3., rrI[5], 1e-5, testName); // from mbpol
    ASSERT_EQUAL_TOL_MOD(             -2.13404e+07, rrI[7], 1e-5, testName); // from this plugin after integration testing with mbpol on water3

}

static void testGetAndScaleInverseRsInterMulecolar() {

    std::string testName      = "testGetAndScaleInverseRsInterMulecolar";

    RealOpenMM damp=0.;
    RealOpenMM dampO=0.001310;
    MapIntRealOpenMM rrI;
    RealOpenMM r=2.771936396e+00*1e-1; // from Water3 test
    RealOpenMM thole=0.400;

    WrappedMBPolReferenceElectrostaticsForce* mbpolReferenceElectrostaticsForce = new WrappedMBPolReferenceElectrostaticsForce();;
    mbpolReferenceElectrostaticsForce->wrapGetAndScaleInverseRs( dampO, dampO,
                          thole, thole, r, false, damp, rrI);

    ASSERT_EQUAL_TOL_MOD(3.607586381e-01*1e1, rrI[1], 1e-5, testName); // from mbpol
    ASSERT_EQUAL_TOL_MOD(4.695157736e-02*1e3, rrI[3], 1e-5, testName); // from mbpol
    ASSERT_EQUAL_TOL_MOD(6.110587933e-03*1e5*3., rrI[5], 1e-5, testName); // from mbpol
    ASSERT_EQUAL_TOL_MOD(119289., rrI[7], 1e-5, testName); // from this plugin after integration testing with mbpol on water3
}

class WrappedMBPolReferenceElectrostaticsForceForIndDipole : public MBPolReferenceElectrostaticsForce {
    public:
    void wrapCalculateInducedDipolePairIxns()   {
        string testName = "computeInducedDipoles";
        std::cout << "wrapCalculateInducedDipolePairIxns" << std::endl;

        int numberOfParticles = 2;
        std::vector<RealVec> positions(numberOfParticles);
        positions[0]             = RealVec( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  );
        positions[1]             = RealVec( -1.763651687e+00, -3.816594649e-01, -1.300353949e+00  );

        for (int i=0; i<numberOfParticles; i++) {
             for (int j=0; j<3; j++) {
                 positions[i][j] *= 1e-1;
             }
         }

        std::vector<RealOpenMM> charges, dipoles, tholes, dampingFactors, polarity;
        std::vector<RealOpenMM> quadrupoles;
        std::vector<int> intZeros;

        for (int i=0; i<numberOfParticles; i++){
            charges.push_back(-5.1966000e-01);
            for (int j=0; j<5; j++){
                tholes.push_back(0.4);
            }
            dampingFactors.push_back(0.001310);
            polarity.push_back(0.001310);
            for (int j=0; j<3; j++){
                dipoles.push_back(0.);
            }
            for (int j=0; j<6; j++){
                quadrupoles.push_back(0.);
            }
            intZeros.push_back(0);
        }

        std::vector<ElectrostaticsParticleData> particleData;
        _numParticles = numberOfParticles;
        loadParticleData(positions, charges,
                              tholes, dampingFactors, polarity, intZeros, intZeros, intZeros, particleData );

        _fixedElectrostaticsField.resize( numberOfParticles );
        _fixedElectrostaticsFieldPolar.resize( numberOfParticles );
        _fixedElectrostaticsField[0]      = RealVec(-6.040604308e-03*1e2, -4.375756834e-03*1e2, -6.721950569e-02*1e2);
        _fixedElectrostaticsFieldPolar[0] = RealVec(0., 0., 0);
        _fixedElectrostaticsField[1]      = RealVec(6.040604308e-03*1e2, 4.375756834e-03*1e2, 6.721950569e-02*1e2);
        _fixedElectrostaticsFieldPolar[1] = RealVec(0., 0., 0);

        for( unsigned int ii = 0; ii < _numParticles; ii++ ){
            _fixedElectrostaticsField[ii]      *= particleData[ii].polarity;
            _fixedElectrostaticsFieldPolar[ii] *= particleData[ii].polarity;
        }
        _inducedDipole.resize( numberOfParticles );
        _inducedDipolePolar.resize( numberOfParticles );
        std::vector<UpdateInducedDipoleFieldStruct> updateInducedDipoleField;
        updateInducedDipoleField.push_back( UpdateInducedDipoleFieldStruct( &_fixedElectrostaticsField,       &_inducedDipole ) );
        updateInducedDipoleField.push_back( UpdateInducedDipoleFieldStruct( &_fixedElectrostaticsFieldPolar,  &_inducedDipolePolar ) );

//        std::cout << "initializeInducedDipoles" << std::endl;

        initializeInducedDipoles( updateInducedDipoleField );

//        for( unsigned int ii = 0; ii < numberOfParticles; ii++ ){
//            std::cout << updateInducedDipoleField[0].inducedDipoles[0][ii] << std::endl;
//        }

//        std::cout << "calculateInducedDipolePairIxns" << std::endl;

        convergeInduceDipoles( particleData, updateInducedDipoleField );

//        for( unsigned int ii = 0; ii < numberOfParticles; ii++ ){
//            std::cout << "******** Particle " << ii << std::endl;
//            //"[ protoncharge nm ]"
//            std::cout << "inducedDipoles:     " << updateInducedDipoleField[0].inducedDipoles[0][ii]*18.2226*10 << "[ protoncharge A sqrt(A kcal/mol) ]"<< std::endl;
//            std::cout << "fixedElectrostaticsField:" << updateInducedDipoleField[0].fixedElectrostaticsField[0][ii] << "[ Kj/mol/nm ]" << std::endl;
//            std::cout << "inducedDipoleField: " << updateInducedDipoleField[0].inducedDipoleField[ii] << "[ Kj/mol/nm ]" << std::endl;
//        }

//        std::cout << "END of wrapCalculateInducedDipolePairIxns" << std::endl;

        std::vector<Vec3> expectedInducedDipoles(numberOfParticles);
        expectedInducedDipoles[0] = Vec3(-7.046394571e-03, -5.104341822e-03, -7.841188329e-02);
        expectedInducedDipoles[1] = Vec3( 7.046394571e-03,  5.104341822e-03,  7.841188329e-02);
        for (int i=0; i<numberOfParticles; i++) {
            for (int j=0; j<3; j++) {
                expectedInducedDipoles[i][j] *= 1e-1;
            }
        }

        double tolerance = 1e-7;
        for (int i=0; i<numberOfParticles; i++) {
            ASSERT_EQUAL_VEC_MOD(expectedInducedDipoles[i], updateInducedDipoleField[0].inducedDipoles[0][i], tolerance, testName);
        }
    }
};

class WrappedMBPolReferenceElectrostaticsForceForPmeDipole: public MBPolReferencePmeElectrostaticsForce {
    public:
    void wrapCalculateInducedDipolePairIxns()   {
        string testName = "computeInducedDipolesPme";
        std::cout << testName << std::endl;

        int numberOfParticles = 2;
        std::vector<RealVec> positions(numberOfParticles);
        positions[0]             = RealVec( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  );
        positions[1]             = RealVec( -1.763651687e+00, -3.816594649e-01, -1.300353949e+00  );

        for (int i=0; i<numberOfParticles; i++) {
             for (int j=0; j<3; j++) {
                positions[i][j] *= 1e-1;
             }
         }

        std::vector<RealOpenMM> charges, dipoles, tholes, dampingFactors, polarity;
        std::vector<RealOpenMM> quadrupoles;
        std::vector<int> intZeros;

        for (int i=0; i<numberOfParticles; i++){
            charges.push_back(-5.1966000e-01);
            for (int j=0; j<5; j++){
                tholes.push_back(0.4);
            }
            dampingFactors.push_back(0.001310);
            polarity.push_back(0.001310);
            intZeros.push_back(0);
        }

        std::vector<ElectrostaticsParticleData> particleData;
        _numParticles = numberOfParticles;
        loadParticleData(positions, charges,
                              tholes, dampingFactors, polarity, intZeros, intZeros, intZeros, particleData );

        _fixedElectrostaticsField.resize( numberOfParticles );
        _fixedElectrostaticsFieldPolar.resize( numberOfParticles );
        _fixedElectrostaticsField[0]      = RealVec(-6.040604308e-03*1e2, -4.375756834e-03*1e2, -6.721950569e-02*1e2);
        _fixedElectrostaticsFieldPolar[0] = RealVec(0., 0., 0);
        _fixedElectrostaticsField[1]      = RealVec(6.040604308e-03*1e2, 4.375756834e-03*1e2, 6.721950569e-02*1e2);
        _fixedElectrostaticsFieldPolar[1] = RealVec(0., 0., 0);

        for( unsigned int ii = 0; ii < _numParticles; ii++ ){
            _fixedElectrostaticsField[ii]      *= particleData[ii].polarity;
            _fixedElectrostaticsFieldPolar[ii] *= particleData[ii].polarity;
        }
        _inducedDipole.resize( numberOfParticles );
        _inducedDipolePolar.resize( numberOfParticles );
        std::vector<UpdateInducedDipoleFieldStruct> updateInducedDipoleField;
        updateInducedDipoleField.push_back( UpdateInducedDipoleFieldStruct( &_fixedElectrostaticsField,       &_inducedDipole ) );
        updateInducedDipoleField.push_back( UpdateInducedDipoleFieldStruct( &_fixedElectrostaticsFieldPolar,  &_inducedDipolePolar ) );

//        std::cout << "initializeInducedDipoles" << std::endl;

        //initializeInducedDipoles( updateInducedDipoleField );
        _inducedDipole.resize( _numParticles );
        _inducedDipolePolar.resize( _numParticles );

        for( unsigned int ii = 0; ii < _numParticles; ii++ ){
            _inducedDipole[ii]       = _fixedElectrostaticsField[ii];
            _inducedDipolePolar[ii]  = _fixedElectrostaticsFieldPolar[ii];
        }

        resizePmeArrays();

//      for( unsigned int ii = 0; ii < numberOfParticles; ii++ ){
//          std::cout << updateInducedDipoleField[0].inducedDipoles[0][ii] << std::endl;
//      }

//        std::cout << "calculateInducedDipolePairIxns" << std::endl;

        convergeInduceDipoles( particleData, updateInducedDipoleField );

      for( unsigned int ii = 0; ii < numberOfParticles; ii++ ){
          std::cout << "******** Particle " << ii << std::endl;
          //"[ protoncharge nm ]"
          std::cout << "inducedDipoles:     " << updateInducedDipoleField[0].inducedDipoles[0][ii]*18.2226*10 << "[ protoncharge A sqrt(A kcal/mol) ]"<< std::endl;
          std::cout << "fixedElectrostaticsField:" << updateInducedDipoleField[0].fixedElectrostaticsField[0][ii] << "[ Kj/mol/nm ]" << std::endl;
          std::cout << "inducedDipoleField: " << updateInducedDipoleField[0].inducedDipoleField[ii] << "[ Kj/mol/nm ]" << std::endl;
      }

//        std::cout << "END of wrapCalculateInducedDipolePairIxns" << std::endl;

        std::vector<Vec3> expectedInducedDipoles(numberOfParticles);
        expectedInducedDipoles[0] = Vec3(-7.046394571e-03, -5.104341822e-03, -7.841188329e-02);
        expectedInducedDipoles[1] = Vec3( 7.046394571e-03,  5.104341822e-03,  7.841188329e-02);
        for (int i=0; i<numberOfParticles; i++) {
            for (int j=0; j<3; j++) {
                expectedInducedDipoles[i][j] *= 1e-1;
            }
        }

        double tolerance = 1e-7;
        for (int i=0; i<numberOfParticles; i++) {
            ASSERT_EQUAL_VEC_MOD(expectedInducedDipoles[i], updateInducedDipoleField[0].inducedDipoles[0][i], tolerance, testName);
        }
    }
};

class WrappedMBPolReferenceElectrostaticsForceForComputeWaterCharge : public MBPolReferenceElectrostaticsForce {
    public:
    void testComputeWaterCharge()   {
        string testName = "testComputeWaterCharge";

        int numberOfParticles = 4;
        std::vector<RealVec> positions(numberOfParticles);
        positions[0]             = Vec3( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  );
        positions[1]             = Vec3( -6.218989773e-01, -6.009430735e-01,  1.572437625e+00  );
        positions[2]             = Vec3( -2.017613812e+00, -4.190350349e-01,  2.239642849e+00  );
        positions[3]             = Vec3( -1.43230412, -0.33360265,  1.64727446 );

        for (int i=0; i<numberOfParticles; i++) {
             for (int j=0; j<3; j++) {
                positions[i][j] *= 1e-1;
             }
         }

        ElectrostaticsParticleData particleO;
        ElectrostaticsParticleData particleH1;
        ElectrostaticsParticleData particleH2;
        ElectrostaticsParticleData particleM;
        for (int j=0; j<3; j++) {
            particleO.position[j]  = positions[0][j];
            particleH1.position[j] = positions[1][j];
            particleH2.position[j] = positions[2][j];
            particleM.position[j]  = positions[3][j];
        }

        computeWaterCharge(particleO, particleH1, particleH2, particleM);

        ASSERT_EQUAL_TOL_MOD(0., particleO.charge, 1e-5, testName);
        ASSERT_EQUAL_TOL_MOD(0.573599422, particleH1.charge, 1e-5, testName);
        ASSERT_EQUAL_TOL_MOD(0.577197137, particleH2.charge, 1e-5, testName);
        ASSERT_EQUAL_TOL_MOD(-1.15079656, particleM.charge, 1e-5, testName);

//        std::cout << "Charges" << std::endl;
//
//        std::cout << "O: " << particleO.charge << std::endl;
//        std::cout << "H1: " << particleH1.charge << std::endl;
//        std::cout << "H2: " << particleH2.charge << std::endl;
//        std::cout << "M: " << particleM.charge << std::endl;
//
//        std::cout << "Derivatives" << std::endl;
//
//        std::cout << "H1 vs H1: " << particleH1.chargeDerivatives[0] << std::endl;
//        std::cout << "H1 vs H2: " << particleH1.chargeDerivatives[1] << std::endl;
//        std::cout << "H1 vs M : " << particleH1.chargeDerivatives[2] << std::endl;
//
//        std::cout << "H2 vs H1: " << particleH2.chargeDerivatives[0] << std::endl;
//        std::cout << "H2 vs H2: " << particleH2.chargeDerivatives[1] << std::endl;
//        std::cout << "H2 vs M : " << particleH2.chargeDerivatives[2] << std::endl;
//
//        std::cout << "M vs H1: " << particleM.chargeDerivatives[0] << std::endl;
//        std::cout << "M vs H2: " << particleM.chargeDerivatives[1] << std::endl;
//        std::cout << "M vs M : " << particleM.chargeDerivatives[2] << std::endl;
//        std::cout << "O vs H1: " << particleO.chargeDerivatives[0] << std::endl;
//        std::cout << "O vs H2: " << particleO.chargeDerivatives[1] << std::endl;
//        std::cout << "O vs M : " << particleO.chargeDerivatives[2] << std::endl;

        std::vector<Vec3> expectedChargeDerivatives(9);
        expectedChargeDerivatives[0]         = Vec3( -0.224842979, 0.157051233, -0.139425246  );
        expectedChargeDerivatives[1]         = Vec3( -0.118671613, 0.106113269, -0.118471774 );
        expectedChargeDerivatives[2]         = Vec3(  0.343514592, -0.263164503, 0.25789702 );
        expectedChargeDerivatives[3]         = Vec3( -0.00533173093, 0.0989902789, -0.187436499  );
        expectedChargeDerivatives[4]         = Vec3( 0.065462366, 0.123151092, -0.285810407 );
        expectedChargeDerivatives[5]         = Vec3( -0.060130635, -0.222141371, 0.473246906  );
        expectedChargeDerivatives[6]         = Vec3(  0.23017471, -0.256041512, 0.326861745  );
        expectedChargeDerivatives[7]         = Vec3( 0.0532092469, -0.229264361, 0.404282181  );
        expectedChargeDerivatives[8]         = Vec3( -0.283383957, 0.485305874, -0.731143926  );

        for (int i=0; i<9; i++) {
            for (int j=0; j<3; j++) {
                expectedChargeDerivatives[i][j] *= 10;
            }
        }

        double tolerance = 1e-6;
        int start_i = 0;
        for (int i=0; i<3; i++) {
            ASSERT_EQUAL_VEC_MOD(particleH1.chargeDerivatives[i], expectedChargeDerivatives[start_i+i], tolerance, testName);
        }
        start_i = 3;
        for (int i=0; i<3; i++) {
            ASSERT_EQUAL_VEC_MOD(particleH2.chargeDerivatives[i], expectedChargeDerivatives[start_i+i], tolerance, testName);
        }
        start_i = 6;
        for (int i=0; i<3; i++) {
            ASSERT_EQUAL_VEC_MOD(particleO.chargeDerivatives[i], expectedChargeDerivatives[start_i+i], tolerance, testName);
        }
    }
};

static void testWater3VirtualSite() {

    std::string testName      = "testWater3VirtualSite";
    std::cout << "Test START: " << testName << std::endl;

    int numberOfParticles     = 12;
    double cutoff             = 0.70;

    std::vector<double> outputElectrostaticsMoments;
    std::vector< Vec3 > inputGrid;
    std::vector< double > outputGridPotential;


    // beginning of Electrostatics setup
    MBPolElectrostaticsForce::NonbondedMethod nonbondedMethod = MBPolElectrostaticsForce::NoCutoff;

    System system;
    // box dimensions

    // double boxDimension                               = 1.8643;
    // Vec3 a( boxDimension, 0.0, 0.0 );
    // Vec3 b( 0.0, boxDimension, 0.0 );
    // Vec3 c( 0.0, 0.0, boxDimension );
    // system.setDefaultPeriodicBoxVectors( a, b, c );

    MBPolElectrostaticsForce* mbpolElectrostaticsForce        = new MBPolElectrostaticsForce();;
    mbpolElectrostaticsForce->setNonbondedMethod( nonbondedMethod );
    //mbpolElectrostaticsForce->setPolarizationType( polarizationType );
    //mbpolElectrostaticsForce->setCutoffDistance( cutoff );
    //mbpolElectrostaticsForce->setMutualInducedTargetEpsilon( 1.0e-06 );
    //mbpolElectrostaticsForce->setMutualInducedMaxIterations( 500 );
    //mbpolElectrostaticsForce->setAEwald( 5.4459052e+00 );
    //mbpolElectrostaticsForce->setEwaldErrorTolerance( 1.0e-04 );

    double virtualSiteWeightO = 0.573293118;
    double virtualSiteWeightH = 0.213353441;
    for( unsigned int jj = 0; jj < numberOfParticles; jj += 4 ){
        system.addParticle( 1.5999000e+01 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 0. ); // Virtual Site
        system.setVirtualSite(jj+3, new ThreeParticleAverageSite(jj, jj+1, jj+2,
                                                           virtualSiteWeightO, virtualSiteWeightH,virtualSiteWeightH));

    }

    std::vector<double> zeroDipole(3);
    std::vector<double> zeroQuadrupole(9);
    std::vector<double> thole(5);

    std::fill(zeroDipole.begin(), zeroDipole.end(), 0.);
    std::fill(zeroQuadrupole.begin(), zeroQuadrupole.end(), 0.);

    thole[TCC] = 0.4;
    thole[TCD] = 0.4;
    thole[TDD] = 0.055;
    thole[TDDOH]  = 0.626;
    thole[TDDHH] = 0.055;

    for( unsigned int jj = 0; jj < numberOfParticles; jj += 4 ){
        mbpolElectrostaticsForce->addElectrostatics( -5.1966000e-01, jj+1, jj+2, jj+3,
                                            thole, 0.001310, 0.001310 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+2, jj+3,
                                            thole, 0.000294, 0.000294 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+1, jj+3,
                                            thole, 0.000294, 0.000294 );
        mbpolElectrostaticsForce->addElectrostatics(  0., jj, jj+1, jj+2,
                                                    thole,  0.001310,  0.);
    }

    system.addForce(mbpolElectrostaticsForce);

    static std::vector<Vec3> positions; // Static to work around bug in Visual Studio that makes compilation very very slow.
    positions.resize(numberOfParticles);

    positions[0]             = Vec3( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  );
    positions[1]             = Vec3( -6.218989773e-01, -6.009430735e-01,  1.572437625e+00  );
    positions[2]             = Vec3( -2.017613812e+00, -4.190350349e-01,  2.239642849e+00  );
    positions[3]             = Vec3( -1.43230412, -0.33360265,  1.64727446 );


    positions[4]             = Vec3( -1.763651687e+00, -3.816594649e-01, -1.300353949e+00  );
    positions[5]             = Vec3( -1.903851736e+00, -4.935677617e-01, -3.457810126e-01  );
    positions[6]             = Vec3( -2.527904158e+00, -7.613550077e-01, -1.733803676e+00  );
    positions[7]             = Vec3( -1.95661974, -0.48654484, -1.18917052 );


    positions[8]             = Vec3( -5.588472140e-01,  2.006699172e+00, -1.392786582e-01  );
    positions[9]             = Vec3( -9.411558180e-01,  1.541226676e+00,  6.163293071e-01  );
    positions[10]            = Vec3( -9.858551734e-01,  1.567124294e+00, -8.830970941e-01  );
    positions[11]            = Vec3( -0.73151769,  1.8136042 , -0.13676332 );

    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            positions[i][j] *= 1e-1;
        }
    }

    std::string platformName;
    platformName = "Reference";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);
    context.applyConstraints(1e-4); // update position of virtual site

    double tolerance          = 1.0e-04;

//    // test energy and forces
//
    State state                = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces   = state.getForces();
    double energy              = state.getPotentialEnergy();
    double cal2joule = 4.184;

    double expectedEnergy = -15.818784*cal2joule;
    std::cout << "Energy: " << energy/cal2joule << " Kcal/mol "<< std::endl;
    std::cout << "Expected energy: " << expectedEnergy/cal2joule << " Kcal/mol "<< std::endl;

    std::vector<Vec3> expectedForces(numberOfParticles);
    expectedForces[0]         = Vec3(  2.38799956, 0.126835228,   8.86189407  );
    expectedForces[1]         = Vec3( -4.21263312, -0.72316292,   3.37076777  );
    expectedForces[2]         = Vec3(  2.19240288, -2.24806806,   1.96210789  );
    expectedForces[4]         = Vec3(  3.59486021, -2.16710895,   3.57138432  );
    expectedForces[5]         = Vec3( -4.54547068, -4.58639226,  -17.4258666  );
    expectedForces[6]         = Vec3( -3.27239433, -1.96722979,    1.1170853  );
    expectedForces[8]         = Vec3( -1.44387205, -3.22471108,  -2.61329967  );
    expectedForces[9]         = Vec3(  3.35011312,  6.07136704, -0.197008793  );
    expectedForces[10]         = Vec3(  1.94899441,   8.7184708,   1.35293571  );

    // gradient -> forces
    for (int i=0; i<numberOfParticles; i++) {
           for (int j=0; j<3; j++) {
            forces[i][j] /= cal2joule*10;
            if ((i+1) % 4 == 0) { // Set virtual site force to 0
                forces[i][j] = 0;
            }
           }
       }

    for (int i=0; i<numberOfParticles; i++) {
           for (int j=0; j<3; j++) {
               expectedForces[i][j] *= -1;
           }
       }
    std::cout  << std::endl << "Forces:" << std::endl;


    const double eps = 1.0e-4;

    double x_orig;

    std::vector<Vec3> finiteDifferenceForces(numberOfParticles);
    for (int i=0; i<numberOfParticles; i++) {
        finiteDifferenceForces.push_back(Vec3( 0.,  0., 0.  ));
    }
    for (int i=0; i<numberOfParticles; i++) {
        for (int xyz=0; xyz<3; xyz++) {
            x_orig = positions[i][xyz];

            positions[i][xyz] = x_orig + eps;
            context.setPositions(positions);
            context.applyConstraints(1e-4); // update position of virtual site
            state                = context.getState(State::Energy);
            const double Ep  = state.getPotentialEnergy();

            positions[i][xyz] = x_orig + 2*eps;
            context.setPositions(positions);
            context.applyConstraints(1e-4); // update position of virtual site
            state                = context.getState(State::Energy);
            const double E2p  = state.getPotentialEnergy();

            positions[i][xyz] = x_orig - eps;
            context.setPositions(positions);
            context.applyConstraints(1e-4); // update position of virtual site
            state                = context.getState(State::Energy);
            const double Em   = state.getPotentialEnergy();

            positions[i][xyz] = x_orig - 2*eps;
            context.setPositions(positions);
            context.applyConstraints(1e-4); // update position of virtual site
            state                = context.getState(State::Energy);
            const double E2m   = state.getPotentialEnergy();

            finiteDifferenceForces[i][xyz] = (8*(Ep - Em) - (E2p - E2m))/(12*eps);
            positions[i][xyz] = x_orig;
        }

    }

    for (int i=0; i<numberOfParticles; i++) {
           for (int j=0; j<3; j++) {
            finiteDifferenceForces[i][j] /= -1*cal2joule*10;
           }

       }

    for (int i=0; i<numberOfParticles; i++) {
        std::cout << "Force atom " << i << ": " << expectedForces[i] << " Kcal/mol/A <mbpol>" << std::endl;
        std::cout << "Force atom " << i << ": " << forces[i] << " Kcal/mol/A <openmm-mbpol>" << std::endl;
        std::cout << "Force atom " << i << ": " << finiteDifferenceForces[i] << " Kcal/mol/A <openmm-mbpol finite differences>" << std::endl << std::endl;
    }

    std::cout << "Comparison of energy and forces with tolerance: " << tolerance << std::endl << std::endl;

    ASSERT_EQUAL_TOL_MOD( expectedEnergy, energy, tolerance, testName );

    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        ASSERT_EQUAL_VEC_MOD( expectedForces[ii], forces[ii], tolerance, testName );
    }

    std::cout << "Test Successful: " << testName << std::endl << std::endl;


    return;
}

class WrappedMBPolReferenceElectrostaticsForceForCalculateElectrostaticPairIxn : public MBPolReferenceElectrostaticsForce {
    public:

    void testCalculateElectrostaticPairIxn () {

        setIncludeChargeRedistribution(false);

        string testName = "testCalculateElectrostaticPairIxn";
        std::cout << "Test START: " << testName << std::endl;

        int numberOfParticles = 2;
        std::vector<ElectrostaticsParticleData> particleData;
        particleData.resize(numberOfParticles);
        particleData[0].dampingFactor = 0.001310;
        particleData[1].dampingFactor = 0.001310;
        particleData[0].polarity = 0.001310;
        particleData[1].polarity = 0.001310;
        RealOpenMM thole = 0.4;
        for (int i=0; i<5; i++) {
        particleData[0].thole[i] = thole;
        particleData[1].thole[i] = thole;
        }

        particleData[0].charge = -5.1966000e-01;
        particleData[1].charge = -5.1966000e-01;
        std::vector<RealVec> positions(numberOfParticles);
        positions[0]             = RealVec( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  );
        positions[1]             = RealVec( -1.763651687e+00, -3.816594649e-01, -1.300353949e+00  );

        for (int i=0; i<numberOfParticles; i++) {
            particleData[1].particleIndex = i;
            particleData[i].multipoleAtomZs = -1;
            particleData[i].multipoleAtomXs = -1;
            particleData[i].multipoleAtomYs = -1;
            for (int j=0; j<3; j++) {
                particleData[i].position[j] = positions[i][j] * 1e-1;
            }
         }


        _inducedDipole.push_back(Vec3(-0.128403739, -0.0930144581, -1.4288696));
        _inducedDipole.push_back(Vec3(0.128403739, 0.0930144581, 1.4288696));
        for (int i=0; i<numberOfParticles; i++) {
                _inducedDipole[i] *= 1e-1;
                _inducedDipolePolar.push_back(_inducedDipole[i]);
        }


        std::vector<RealVec> forces;
        for (int i=0; i<numberOfParticles; i++) {
            forces.push_back(Vec3( 0.,  0., 0.  ));
        }
        RealOpenMM energy = calculateElectrostaticPairIxn( particleData, 0, 1,forces);

        std::vector<RealVec> expectedForces;
        expectedForces.push_back(Vec3(3.11046, 2.25319, 34.6131));
        expectedForces.push_back(Vec3(-3.11046, -2.25319, -34.6131));

        double tolerance = 1e-5;

        for( unsigned int ii = 0; ii < forces.size(); ii++ ){
            forces[ii] /= (cal2joule * 10);
//            ASSERT_EQUAL_VEC_MOD( expectedForces[ii], forces[ii], tolerance, testName );
            std::cout << expectedForces[ii] << " Kcal/mol/A <expected>" << std::endl;
            std::cout << forces[ii] << " Kcal/mol/A" << std::endl;
        }
        energy /= cal2joule;
//        ASSERT_EQUAL_TOL_MOD( 0.0634441, energy, tolerance, testName );
        std::cout << "Energy: " << energy << " Kcal/mol "<< std::endl;
        std::cout << "Energy: " << 0.0634441 << " Kcal/mol <expected>"<< std::endl;
        std::cout << "Test END: " << testName << std::endl << std::endl;

    }
};
static void testWater3() {

    std::string testName      = "testWater3";

    int numberOfParticles     = 9;
    double cutoff             = 0.70;

    std::vector<double> outputElectrostaticsMoments;
    std::vector< Vec3 > inputGrid;
    std::vector< double > outputGridPotential;


    // beginning of Electrostatics setup
    MBPolElectrostaticsForce::NonbondedMethod nonbondedMethod = MBPolElectrostaticsForce::NoCutoff;

    System system;
    // box dimensions
    MBPolElectrostaticsForce* mbpolElectrostaticsForce        = new MBPolElectrostaticsForce();;
    mbpolElectrostaticsForce->setNonbondedMethod( nonbondedMethod );
    mbpolElectrostaticsForce->setIncludeChargeRedistribution(false);

    unsigned int particlesPerMolecule = 3;

    for( unsigned int jj = 0; jj < numberOfParticles; jj += particlesPerMolecule ){
        system.addParticle( 1.5999000e+01 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 1.0080000e+00 );
    }

    std::vector<double> zeroDipole(3);
    std::vector<double> zeroQuadrupole(9);
    std::vector<double> thole(5);
    std::fill(zeroDipole.begin(), zeroDipole.end(), 0.);
    std::fill(zeroQuadrupole.begin(), zeroQuadrupole.end(), 0.);

    thole[TCC] = 0.4;
    thole[TCD] = 0.4;
    thole[TDD] = 0.055;
    thole[TDDOH]  = 0.626;
    thole[TDDHH] = 0.055;

    for( unsigned int jj = 0; jj < numberOfParticles; jj += particlesPerMolecule ){
        mbpolElectrostaticsForce->addElectrostatics( -5.1966000e-01, jj+1, jj+2, -1,
                                            thole, 0.001310, 0.001310 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+2, -1,
                                            thole, 0.000294, 0.000294 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+1, -1,
                                            thole, 0.000294, 0.000294 );
    }

    system.addForce(mbpolElectrostaticsForce);

    static std::vector<Vec3> positions; // Static to work around bug in Visual Studio that makes compilation very very slow.
    positions.resize(numberOfParticles);

    positions[0]             = Vec3( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  );
    positions[1]             = Vec3( -6.218989773e-01, -6.009430735e-01,  1.572437625e+00  );
    positions[2]             = Vec3( -2.017613812e+00, -4.190350349e-01,  2.239642849e+00  );

    positions[3]             = Vec3( -1.763651687e+00, -3.816594649e-01, -1.300353949e+00  );
    positions[4]             = Vec3( -1.903851736e+00, -4.935677617e-01, -3.457810126e-01  );
    positions[5]             = Vec3( -2.527904158e+00, -7.613550077e-01, -1.733803676e+00  );

    positions[6]             = Vec3( -5.588472140e-01,  2.006699172e+00, -1.392786582e-01  );
    positions[7]             = Vec3( -9.411558180e-01,  1.541226676e+00,  6.163293071e-01  );
    positions[8]             = Vec3( -9.858551734e-01,  1.567124294e+00, -8.830970941e-01  );

    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            positions[i][j] *= 1e-1;
        }
    }

    std::string platformName;
    platformName = "Reference";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);
    context.applyConstraints(1e-4); // update position of virtual site

    double tolerance          = 1.0e-04;

//    // test energy and forces
//
    State state                = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces   = state.getForces();
    double energy              = state.getPotentialEnergy();
//    std::cout << "Forces" << std::endl;


//    std::vector<Vec3> expectedForces(numberOfParticles);
//    expectedForces[0]         = Vec3( -1.029233628e-01,  1.752006876e-01, -2.394228296e-01  );
//    expectedForces[1]         = Vec3(  1.238286503e-01, -9.713944883e-02,  9.278441270e-02  );
//    expectedForces[2]         = Vec3( -1.992936921e-02, -8.084103617e-02,  1.660930712e-01  );
//    expectedForces[3]         = Vec3(  2.181116801e-01,  1.127169979e-01, -1.998507867e-01  );
//    expectedForces[4]         = Vec3( -1.021411513e-01, -6.244910893e-02,  1.595471969e-01  );
//    expectedForces[5]         = Vec3( -1.214347018e-01, -6.329887574e-02,  2.105405984e-02  );
//    expectedForces[6]         = Vec3(  1.708442625e-01,  1.860776100e-01,  2.249030303e-02  );
//    expectedForces[7]         = Vec3( -7.205290616e-02, -7.830256131e-02,  4.942309713e-02  );
//    expectedForces[8]         = Vec3( -9.430310162e-02, -9.196426456e-02, -7.211852443e-02  );
//    for (int i=0; i<numberOfParticles; i++) {
//        for (int j=0; j<3; j++) {
//            expectedForces[i][j] *= cal2joule*10;
//        }
//    }


    //for( unsigned int ii = 0; ii < forces.size(); ii++ ){
    //    ASSERT_EQUAL_VEC_MOD( expectedForces[ii], forces[ii], tolerance, testName );
    //}

    for (int i=0; i<numberOfParticles; i++) {
           for (int j=0; j<3; j++) {
            forces[i][j] /= cal2joule*10;
           }
       }

    std::cout << "Test start: " << testName << std::endl;

    std::cout  << std::endl << "Forces:" << std::endl;

    for (int i=0; i<numberOfParticles; i++) {
         std::cout << forces[i] << " Kcal/mol/A " << std::endl;
    }
    // Energy elec+ind(kcal/mol): -2.134083549e-02
    double expectedEnergy = -19.6545*cal2joule;
    // ASSERT_EQUAL_TOL_MOD( expectedEnergy, energy, tolerance, testName );
    std::cout << "Energy: " << energy/cal2joule << " Kcal/mol "<< std::endl;
    std::cout << "Expected energy: " << expectedEnergy/cal2joule << " Kcal/mol "<< std::endl;
    const double eps = 1.0e-4;

    double x_orig;

    std::vector<Vec3> finiteDifferenceForces(numberOfParticles);
    for (int i=0; i<numberOfParticles; i++) {
        finiteDifferenceForces.push_back(Vec3( 0.,  0., 0.  ));
    }
    for (int i=0; i<numberOfParticles; i++) {
        for (int xyz=0; xyz<3; xyz++) {
            x_orig = positions[i][xyz];

            positions[i][xyz] = x_orig + eps;
            context.setPositions(positions);
            state                = context.getState(State::Energy);
            const double Ep  = state.getPotentialEnergy();

            positions[i][xyz] = x_orig + 2*eps;
            context.setPositions(positions);
            state                = context.getState(State::Energy);
            const double E2p  = state.getPotentialEnergy();

            positions[i][xyz] = x_orig - eps;
            context.setPositions(positions);
            state                = context.getState(State::Energy);
            const double Em   = state.getPotentialEnergy();

            positions[i][xyz] = x_orig - 2*eps;
            context.setPositions(positions);
            state                = context.getState(State::Energy);
            const double E2m   = state.getPotentialEnergy();

            finiteDifferenceForces[i][xyz] = (8*(Ep - Em) - (E2p - E2m))/(12*eps);
            positions[i][xyz] = x_orig;
        }

    }

    // Flip sign to convert gradient -> forces
    for (int i=0; i<numberOfParticles; i++) {
           for (int j=0; j<3; j++) {
            finiteDifferenceForces[i][j] /= -1*cal2joule*10;
           }

       }
    std::cout  << std::endl << "Finite difference forces:" << std::endl;


    for (int i=0; i<numberOfParticles; i++) {
         std::cout << finiteDifferenceForces[i] << " Kcal/mol/A " << std::endl;
    }
    std::cout << "Test END: " << testName << std::endl << std::endl;

    return;
}


static void testWater3VirtualSitePMEHugeBox() {

    std::string testName      = "testWater3VirtualSitePMEHugeBox";
    std::cout << "Test START: " << testName << std::endl;

    int numberOfParticles     = 4*3;
    double cutoff             = 10.;

    std::vector<double> outputElectrostaticsMoments;
    std::vector< Vec3 > inputGrid;
    std::vector< double > outputGridPotential;
    double mutualInducedDipoleTargetEpsilon = 1e-4;

    // beginning of Electrostatics setup
    MBPolElectrostaticsForce::NonbondedMethod nonbondedMethod = MBPolElectrostaticsForce::PME;

    System system;

    double boxDimension                               = 50;
    Vec3 a( boxDimension, 0.0, 0.0 );
    Vec3 b( 0.0, boxDimension, 0.0 );
    Vec3 c( 0.0, 0.0, boxDimension );
    system.setDefaultPeriodicBoxVectors( a, b, c );

    MBPolElectrostaticsForce* mbpolElectrostaticsForce        = new MBPolElectrostaticsForce();;
    mbpolElectrostaticsForce->setNonbondedMethod( nonbondedMethod );
    mbpolElectrostaticsForce->setCutoffDistance( cutoff );
    mbpolElectrostaticsForce->setIncludeChargeRedistribution(true);
    mbpolElectrostaticsForce->setMutualInducedTargetEpsilon(mutualInducedDipoleTargetEpsilon);
    //mbpolElectrostaticsForce->setIncludeChargeRedistribution(false);

    // disable Ewald by setting alpha to very low value
    mbpolElectrostaticsForce->setAEwald( 1e-15 );

    std::vector<int> pmeGridDimension( 3 );
    int inputPmeGridDimension = 20;
    pmeGridDimension[0] = pmeGridDimension[1] = pmeGridDimension[2] = inputPmeGridDimension;
    mbpolElectrostaticsForce->setPmeGridDimensions( pmeGridDimension );

    double virtualSiteWeightO = 0.573293118;
    double virtualSiteWeightH = 0.213353441;
    for( unsigned int jj = 0; jj < numberOfParticles; jj += 4 ){
        system.addParticle( 1.5999000e+01 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 0. ); // Virtual Site
        system.setVirtualSite(jj+3, new ThreeParticleAverageSite(jj, jj+1, jj+2,
                                                           virtualSiteWeightO, virtualSiteWeightH,virtualSiteWeightH));

    }

    std::vector<double> zeroDipole(3);
    std::vector<double> zeroQuadrupole(9);
    std::vector<double> thole(5);

    std::fill(zeroDipole.begin(), zeroDipole.end(), 0.);
    std::fill(zeroQuadrupole.begin(), zeroQuadrupole.end(), 0.);

    thole[TCC] = 0.4;
    thole[TCD] = 0.4;
    thole[TDD] = 0.055;
    thole[TDDOH]  = 0.626;
    thole[TDDHH] = 0.055;

    // testing with no polarizability

    for( unsigned int jj = 0; jj < numberOfParticles; jj += 4 ){
        mbpolElectrostaticsForce->addElectrostatics( -5.1966000e-01, jj+1, jj+2, jj+3,
                                            thole, 0.001310, 0.000000 );
                                            //thole, 0.001310, 0.001310 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+2, jj+3,
                                            thole, 0.000294, 0.000000 );
                                            //thole, 0.000294, 0.000294 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+1, jj+3,
                                            thole, 0.000294, 0.000000 );
                                            //thole, 0.000294, 0.000294 );
        mbpolElectrostaticsForce->addElectrostatics(  0., jj, jj+1, jj+2,
                                                    thole,  0.001310,  0.);
    }

    system.addForce(mbpolElectrostaticsForce);

    static std::vector<Vec3> positions; // Static to work around bug in Visual Studio that makes compilation very very slow.
    positions.resize(numberOfParticles);

    positions[0]             = Vec3( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  );
    positions[1]             = Vec3( -6.218989773e-01, -6.009430735e-01,  1.572437625e+00  );
    positions[2]             = Vec3( -2.017613812e+00, -4.190350349e-01,  2.239642849e+00  );
    positions[3]             = Vec3( -1.43230412, -0.33360265,  1.64727446 );

    if (numberOfParticles > 4) {
    positions[4]             = Vec3( -1.763651687e+00, -3.816594649e-01, -1.300353949e+00  );
    positions[5]             = Vec3( -1.903851736e+00, -4.935677617e-01, -3.457810126e-01  );
    positions[6]             = Vec3( -2.527904158e+00, -7.613550077e-01, -1.733803676e+00  );
    positions[7]             = Vec3( -1.95661974, -0.48654484, -1.18917052 );
    }

    if (numberOfParticles > 8) {
    positions[8]             = Vec3( -5.588472140e-01,  2.006699172e+00, -1.392786582e-01  );
    positions[9]             = Vec3( -9.411558180e-01,  1.541226676e+00,  6.163293071e-01  );
    positions[10]            = Vec3( -9.858551734e-01,  1.567124294e+00, -8.830970941e-01  );
    positions[11]            = Vec3( -0.73151769,  1.8136042 , -0.13676332 );
    }

    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            positions[i][j] *= 1e-1;
        }
    }

    std::string platformName;
    platformName = "Reference";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);
    context.applyConstraints(1e-4); // update position of virtual site

    double tolerance          = 1.0e-04;

//    // test energy and forces
//
    State state                = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces   = state.getForces();
    double energy              = state.getPotentialEnergy();
    double cal2joule = 4.184;

    double expectedEnergy = -13.0493*cal2joule;
    std::cout << "Energy: " << energy/cal2joule << " Kcal/mol "<< std::endl;
    std::cout << "Expected energy: " << expectedEnergy/cal2joule << " Kcal/mol "<< std::endl;

    std::cout  << std::endl << "AEwald:" << mbpolElectrostaticsForce->getAEwald() << std::endl;
    mbpolElectrostaticsForce->getPmeGridDimensions(pmeGridDimension);
    std::cout  << std::endl << "PmeGridDimensions:" << pmeGridDimension[0] << std::endl;


    std::vector<Vec3> expectedForces(4*3);

    expectedForces[0]         = Vec3( -1.4269, 0.448043, -5.04341   );
    expectedForces[1]         = Vec3( 2.9682, 0.701196, -3.20899    );
    expectedForces[2]         = Vec3( -1.63706, 1.94763, -1.8898    );
    expectedForces[4]         = Vec3( -2.92071, 1.4094, -2.05605    );
    expectedForces[5]         = Vec3( 3.38821, 2.85987, 11.0781     );
    expectedForces[6]         = Vec3( 2.3639, 1.41592, -0.900715    );
    expectedForces[8]         = Vec3( 0.537103, 2.18093, 1.70919    );
    expectedForces[9]         = Vec3( -2.14789, -4.38903, 0.972416 );
    expectedForces[10]         = Vec3(-1.12484, -6.57397, -0.660752  );

    // gradient -> forces
    for (int i=0; i<numberOfParticles; i++) {
           for (int j=0; j<3; j++) {
            forces[i][j] /= cal2joule*10;
            if ((i+1) % 4 == 0) { // Set virtual site force to 0
                forces[i][j] = 0;
            }
           }
       }

    std::cout  << std::endl << "Forces:" << std::endl;

    const double eps = 1.0e-4;

    double x_orig;

    std::vector<Vec3> finiteDifferenceForces(numberOfParticles);
    for (int i=0; i<numberOfParticles; i++) {
        finiteDifferenceForces.push_back(Vec3( 0.,  0., 0.  ));
    }

    for (int i=0; i<numberOfParticles; i++) {
    #ifdef COMPUTE_FINITE_DIFFERENCES_FORCES
        for (int xyz=0; xyz<3; xyz++) {
            x_orig = positions[i][xyz];

            positions[i][xyz] = x_orig + eps;
            context.setPositions(positions);
            context.applyConstraints(1e-4); // update position of virtual site
            state                = context.getState(State::Energy);
            const double Ep  = state.getPotentialEnergy();

            positions[i][xyz] = x_orig + 2*eps;
            context.setPositions(positions);
            context.applyConstraints(1e-4); // update position of virtual site
            state                = context.getState(State::Energy);
            const double E2p  = state.getPotentialEnergy();

            positions[i][xyz] = x_orig - eps;
            context.setPositions(positions);
            context.applyConstraints(1e-4); // update position of virtual site
            state                = context.getState(State::Energy);
            const double Em   = state.getPotentialEnergy();

            positions[i][xyz] = x_orig - 2*eps;
            context.setPositions(positions);
            context.applyConstraints(1e-4); // update position of virtual site
            state                = context.getState(State::Energy);
            const double E2m   = state.getPotentialEnergy();

            finiteDifferenceForces[i][xyz] = (8*(Ep - Em) - (E2p - E2m))/(12*eps);
        finiteDifferenceForces[i][xyz] /= -1*cal2joule*10;
            positions[i][xyz] = x_orig;
        }
    #endif
        std::cout << "Force atom " << i << ": " << forces[i] << " Kcal/mol/A <openmm-mbpol>" << std::endl;
        std::cout << "Force atom " << i << ": " << expectedForces[i] << " Kcal/mol/A <precomputerd finite differences>" << std::endl;
#ifdef COMPUTE_FINITE_DIFFERENCES_FORCES
        std::cout << "Force atom " << i << ": " << finiteDifferenceForces[i] << " Kcal/mol/A <openmm-mbpol finite differences>" << std::endl;
#endif
        std::cout << std::endl;
    }

    std::cout << "Comparison of energy and forces with tolerance: " << tolerance << std::endl << std::endl;

    ASSERT_EQUAL_TOL_MOD( expectedEnergy, energy, tolerance, testName );

    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        ASSERT_EQUAL_VEC_MOD( expectedForces[ii], forces[ii], tolerance, testName );
    }

    std::cout << "Test Successful: " << testName << std::endl << std::endl;


    return;
}

static void testWater3VirtualSitePMESmallBox() {

    double cutoff             = 10.;
    cutoff = .9;
    int numberOfParticles     = 4*3;
    double boxDimension                               = 1.8;
    bool zeroPolarizability = false;
    bool includeChargeRedistribution = true;


    std::string testName      = "testWater3VirtualSitePMESmallBox";
    std::cout << "Test START: " << testName << std::endl;


    std::vector<double> outputElectrostaticsMoments;
    std::vector< Vec3 > inputGrid;
    std::vector< double > outputGridPotential;


    // beginning of Electrostatics setup
    MBPolElectrostaticsForce::NonbondedMethod nonbondedMethod = MBPolElectrostaticsForce::PME;

    System system;

    Vec3 a( boxDimension, 0.0, 0.0 );
    Vec3 b( 0.0, boxDimension, 0.0 );
    Vec3 c( 0.0, 0.0, boxDimension );
    system.setDefaultPeriodicBoxVectors( a, b, c );

    MBPolElectrostaticsForce* mbpolElectrostaticsForce        = new MBPolElectrostaticsForce();;
    mbpolElectrostaticsForce->setNonbondedMethod( nonbondedMethod );
    mbpolElectrostaticsForce->setCutoffDistance( cutoff );
    mbpolElectrostaticsForce->setIncludeChargeRedistribution(includeChargeRedistribution);


    // setting alpha of Ewald to zero triggers automatic estimation of alpha and grid sized based on error tolerance

    if (boxDimension > 3) {
        mbpolElectrostaticsForce->setAEwald( 1e-15 );
        std::vector<int> pmeGridDimension( 3 );
        int inputPmeGridDimension = 20;
        pmeGridDimension[0] = pmeGridDimension[1] = pmeGridDimension[2] = inputPmeGridDimension;
        mbpolElectrostaticsForce->setPmeGridDimensions( pmeGridDimension );
    } else {
        mbpolElectrostaticsForce->setAEwald( 0. );
        mbpolElectrostaticsForce->setEwaldErrorTolerance( 1.0e-03 );
    }

    double virtualSiteWeightO = 0.573293118;
    double virtualSiteWeightH = 0.213353441;
    for( unsigned int jj = 0; jj < numberOfParticles; jj += 4 ){
        system.addParticle( 1.5999000e+01 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 0. ); // Virtual Site
        system.setVirtualSite(jj+3, new ThreeParticleAverageSite(jj, jj+1, jj+2,
                                                           virtualSiteWeightO, virtualSiteWeightH,virtualSiteWeightH));

    }

    std::vector<double> zeroDipole(3);
    std::vector<double> zeroQuadrupole(9);
    std::vector<double> thole(5);

    std::fill(zeroDipole.begin(), zeroDipole.end(), 0.);
    std::fill(zeroQuadrupole.begin(), zeroQuadrupole.end(), 0.);

    thole[TCC] = 0.4;
    thole[TCD] = 0.4;
    thole[TDD] = 0.055;
    thole[TDDOH]  = 0.626;
    thole[TDDHH] = 0.055;

    for( unsigned int jj = 0; jj < numberOfParticles; jj += 4 ){
        if (zeroPolarizability) {
            mbpolElectrostaticsForce->addElectrostatics( -5.1966000e-01, jj+1, jj+2, jj+3,
                                                thole, 0.001310, 0. );
            mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+2, jj+3,
                                                thole, 0.000294, 0. );
            mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+1, jj+3,
                                                thole, 0.000294, 0. );
            std::cout << "POLARIZABILITY SET TO ZERO"<< std::endl;
        } else {
        mbpolElectrostaticsForce->addElectrostatics( -5.1966000e-01, jj+1, jj+2, jj+3,
                                            thole, 0.001310, 0.001310 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+2, jj+3,
                                            thole, 0.000294, 0.000294 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+1, jj+3,
                                            thole, 0.000294, 0.000294 );
        }
        mbpolElectrostaticsForce->addElectrostatics(  0., jj, jj+1, jj+2,
                                                    thole,  0.001310,  0.);
    }

    system.addForce(mbpolElectrostaticsForce);

    static std::vector<Vec3> positions; // Static to work around bug in Visual Studio that makes compilation very very slow.
    positions.resize(numberOfParticles);

    positions[0]             = Vec3( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  );
    positions[1]             = Vec3( -6.218989773e-01, -6.009430735e-01,  1.572437625e+00  );
    positions[2]             = Vec3( -2.017613812e+00, -4.190350349e-01,  2.239642849e+00  );
    positions[3]             = Vec3( -1.43230412, -0.33360265,  1.64727446 );

    if (numberOfParticles > 4) {
        positions[4]             = Vec3( -1.763651687e+00, -3.816594649e-01, -1.300353949e+00  );
        positions[5]             = Vec3( -1.903851736e+00, -4.935677617e-01, -3.457810126e-01  );
        positions[6]             = Vec3( -2.527904158e+00, -7.613550077e-01, -1.733803676e+00  );
        positions[7]             = Vec3( -1.95661974, -0.48654484, -1.18917052 );


        positions[8]             = Vec3( -5.588472140e-01,  2.006699172e+00, -1.392786582e-01  );
        positions[9]             = Vec3( -9.411558180e-01,  1.541226676e+00,  6.163293071e-01  );
        positions[10]            = Vec3( -9.858551734e-01,  1.567124294e+00, -8.830970941e-01  );
        positions[11]            = Vec3( -0.73151769,  1.8136042 , -0.13676332 );
    }

    for (int i=0; i<numberOfParticles; i++) {
        for (int j=0; j<3; j++) {
            positions[i][j] *= 1e-1;
        }
    }

    std::string platformName;
    platformName = "Reference";
    LangevinIntegrator integrator(0.0, 0.1, 0.01);
    Context context(system, integrator, Platform::getPlatformByName( platformName ) );

    context.setPositions(positions);
    context.applyConstraints(1e-4); // update position of virtual site

    double tolerance          = 1.0e-02;

//    // test energy and forces
//

    State state                = context.getState(State::Forces | State::Energy);
    std::vector<Vec3> forces   = state.getForces();
    double energy              = state.getPotentialEnergy();
    double cal2joule = 4.184;

    double expectedEnergy = -66.7426;
    std::cout << "Energy: " << energy/cal2joule << " Kcal/mol "<< std::endl;
    std::cout << "Expected energy: " << expectedEnergy/cal2joule << " Kcal/mol "<< std::endl;

    std::vector<Vec3> expectedForces;

//// finite difference forces no polarizability
//    expectedForces.push_back(Vec3( -1.37135, 0.53378, -5.07505 ));
//    expectedForces.push_back(Vec3( 2.92977, 0.652437, -3.18227 ));
//    expectedForces.push_back(Vec3( -1.65258, 1.91069, -1.88755 ));
//    expectedForces.push_back(Vec3( -0, -0, -0 ));
//    expectedForces.push_back(Vec3( -2.89755, 1.49146, -2.11498 ));
//    expectedForces.push_back(Vec3( 3.37435, 2.81564, 11.0826 ));
//    expectedForces.push_back(Vec3( 2.35644, 1.37722, -0.849565 ));
//    expectedForces.push_back(Vec3( -0, -0, -0 ));
//    expectedForces.push_back(Vec3( 0.563548, 2.23457, 1.6571 ));
//    expectedForces.push_back(Vec3( -2.15644, -4.41274, 0.987988 ));
//    expectedForces.push_back(Vec3( -1.14349, -6.60323, -0.617352 ));
//    expectedForces.push_back(Vec3( -0, -0, -0 ));
//// finite differences for normal case
    expectedForces.push_back(Vec3( -2.32662, -0.0215643, -8.9463 ));
    expectedForces.push_back(Vec3( 4.17346, 0.664051, -3.33561 ));
    expectedForces.push_back(Vec3( -2.21919, 2.20524, -1.9436 ));
    expectedForces.push_back(Vec3( -0, -0, -0 ));
    expectedForces.push_back(Vec3( -3.57208, 2.2761, -3.66722 ));
    expectedForces.push_back(Vec3( 4.54289, 4.54649, 17.4865 ));
    expectedForces.push_back(Vec3( 3.2676, 1.92351, -1.05576 ));
    expectedForces.push_back(Vec3( -0, -0, -0 ));
    expectedForces.push_back(Vec3( 1.47085, 3.31187, 2.5387 ));
    expectedForces.push_back(Vec3( -3.36353, -6.13306, 0.230528 ));
    expectedForces.push_back(Vec3( -1.97036, -8.77306, -1.30534 ));
    expectedForces.push_back(Vec3( -0, -0, -0 ));
//
//    // finite difference forces with no polarizability and no charge redistribution
//    expectedForces.push_back(Vec3( -2.29703, 2.44192, -7.15701 ));
//    expectedForces.push_back(Vec3( 2.22037, -0.960757, 1.1415 ));
//    expectedForces.push_back(Vec3( 0.0101264, -0.494861, 2.03038 ));
//    expectedForces.push_back(Vec3( -0, -0, -0 ));
//    expectedForces.push_back(Vec3( 1.2845, 3.00543, -2.30249 ));
//    expectedForces.push_back(Vec3( -0.0309592, 0.126436, 5.54208 ));
//    expectedForces.push_back(Vec3( -0.248866, -0.380292, 0.50001 ));
//    expectedForces.push_back(Vec3( -0, -0, -0 ));
//    expectedForces.push_back(Vec3( 2.28936, 3.71712, 1.099 ));
//    expectedForces.push_back(Vec3( -1.9654, -3.44716, 0.0428621 ));
//    expectedForces.push_back(Vec3( -1.2611, -4.00799, -0.895414 ));
//    expectedForces.push_back(Vec3( -0, -0, -0 ));

//    // finite difference forces for 1 molecule with no polizability
//    expectedForces.push_back(Vec3( -0.0102713, 0.0168557, -0.0252335 ));
//    expectedForces.push_back(Vec3( -0.000491702, -0.00838144, 0.0143935 ));
//    expectedForces.push_back(Vec3( 0.0126705, -0.00890535, 0.0102875 ));
//    expectedForces.push_back(Vec3( -0, -0, -0 ));

//    // forces with polarizability and charge redistribution HUGE BOX
//    expectedForces.push_back(Vec3(-2.79154, 2.38403, -9.48745 ));
//    expectedForces.push_back(Vec3(2.66988, -0.970286, 1.37622 ));
//    expectedForces.push_back(Vec3(-0.0273603, -0.454075, 2.39462 ));
//    expectedForces.push_back(Vec3(0, 0, 0 ));
//    expectedForces.push_back(Vec3(1.72996, 3.83253, -2.83729 ));
//    expectedForces.push_back(Vec3(0.106845, 0.561475, 7.84569 ));
//    expectedForces.push_back(Vec3(-0.272013, -0.384251, 0.645139 ));
//    expectedForces.push_back(Vec3(0, 0, 0 ));
//    expectedForces.push_back(Vec3(2.8275, 4.2067, 1.50836 ));
//    expectedForces.push_back(Vec3(-2.54637, -4.21511, -0.253532 ));
//    expectedForces.push_back(Vec3(-1.69689, -4.96102, -1.19175 ));
//    expectedForces.push_back(Vec3(0, 0, 0 ));

    // finite difference forces with polarizability and charge redistribution HUGE BOX
//
//    expectedForces.push_back(Vec3(-2.97971, 2.53886, -9.69661 ));
//    expectedForces.push_back(Vec3(2.9417, -1.0957, 1.46631 ));
//    expectedForces.push_back(Vec3(-0.112144, -0.484293, 2.50821 ));
//    expectedForces.push_back(Vec3(-0, -0, -0 ));
//    expectedForces.push_back(Vec3(1.75998, 3.90645, -3.30798 ));
//    expectedForces.push_back(Vec3(0.0199501, 0.442349, 8.29842 ));
//    expectedForces.push_back(Vec3(-0.214005, -0.338206, 0.668481 ));
//    expectedForces.push_back(Vec3(-0, -0, -0 ));
//    expectedForces.push_back(Vec3(3.03821, 4.46561, 1.70258 ));
//    expectedForces.push_back(Vec3(-2.59869, -4.27455, -0.206005 ));
//    expectedForces.push_back(Vec3(-1.8553, -5.16052, -1.43341 ));
//    expectedForces.push_back(Vec3(-0, -0, -0 ));


    // gradient -> forces
    for (int i=0; i<numberOfParticles; i++) {
           for (int j=0; j<3; j++) {
            forces[i][j] /= cal2joule*10;
            if ((i+1) % 4 == 0) { // Set virtual site force to 0
                forces[i][j] = 0;
            }
           }
       }

    std::cout  << std::endl << "Forces:" << std::endl;


    const double eps = 1.0e-4;
#ifdef COMPUTE_FINITE_DIFFERENCES_FORCES

    double x_orig;

    std::vector<Vec3> finiteDifferenceForces(numberOfParticles);
    for (int i=0; i<numberOfParticles; i++) {
        finiteDifferenceForces.push_back(Vec3( 0.,  0., 0.  ));
    }
    for (int i=0; i<numberOfParticles; i++) {

        std::cout  << "Finite differences: particle " << i << "/" << numberOfParticles << std::endl;

        for (int xyz=0; xyz<3; xyz++) {

            x_orig = positions[i][xyz];

            positions[i][xyz] = x_orig + eps;
            context.setPositions(positions);
            context.applyConstraints(1e-4); // update position of virtual site
            state                = context.getState(State::Energy);
            const double Ep  = state.getPotentialEnergy();

            positions[i][xyz] = x_orig + 2*eps;
            context.setPositions(positions);
            context.applyConstraints(1e-4); // update position of virtual site
            state                = context.getState(State::Energy);
            const double E2p  = state.getPotentialEnergy();

            positions[i][xyz] = x_orig - eps;
            context.setPositions(positions);
            context.applyConstraints(1e-4); // update position of virtual site
            state                = context.getState(State::Energy);
            const double Em   = state.getPotentialEnergy();

            positions[i][xyz] = x_orig - 2*eps;
            context.setPositions(positions);
            context.applyConstraints(1e-4); // update position of virtual site
            state                = context.getState(State::Energy);
            const double E2m   = state.getPotentialEnergy();

            finiteDifferenceForces[i][xyz] = (8*(Ep - Em) - (E2p - E2m))/(12*eps);
            positions[i][xyz] = x_orig;
        }

    }

    for (int i=0; i<numberOfParticles; i++) {
           for (int j=0; j<3; j++) {
            finiteDifferenceForces[i][j] /= -1*cal2joule*10;
           }

       }
#endif

    for (int i=0; i<numberOfParticles; i++) {
        std::cout << "Force atom " << i << ": " << forces[i] << " Kcal/mol/A <openmm-mbpol>" << std::endl;
#ifndef COMPUTE_FINITE_DIFFERENCES_FORCES
        std::cout << "Force atom " << i << ": " << expectedForces[i] << " Kcal/mol/A <finite difference forces with polarizability and charge redistribution SMALL BOX>" << std::endl;
#endif
#ifdef COMPUTE_FINITE_DIFFERENCES_FORCES
        std::cout << "Force atom " << i << ": " << finiteDifferenceForces[i] << " Kcal/mol/A <openmm-mbpol finite differences>" << std::endl ;
#endif
        std::cout << std::endl;
    }

    std::cout << "Comparison of energy and forces with tolerance: " << tolerance << std::endl << std::endl;

    ASSERT_EQUAL_TOL_MOD( expectedEnergy, energy, tolerance, testName );

    for( unsigned int ii = 0; ii < forces.size(); ii++ ){
        ASSERT_EQUAL_VEC_MOD( expectedForces[ii], forces[ii], tolerance, testName );
    }

    std::cout << "Test Successful: " << testName << std::endl << std::endl;


    return;
}

class WrappedMBPolReferencePmeElectrostaticsForceForcalculatePmeDirectElectrostaticPairIxn : public MBPolReferencePmeElectrostaticsForce {
    public:

    void testCalculateElectrostaticPairIxn () {

        setIncludeChargeRedistribution(false);

        string testName = "testCalculatePmeDirectElectrostaticPairIxn";
        std::cout << "Test START: " << testName << std::endl;

        int numberOfParticles = 2;
        std::vector<ElectrostaticsParticleData> particleData;
        particleData.resize(numberOfParticles);
        particleData[0].dampingFactor = 0.001310;
        particleData[1].dampingFactor = 0.001310;
        particleData[0].polarity = 0.001310;
        particleData[1].polarity = 0.001310;
        RealOpenMM thole = 0.4;
        for (int i=0; i<5; i++) {
        particleData[0].thole[i] = thole;
        particleData[1].thole[i] = thole;
        }

        particleData[0].charge = -5.1966000e-01;
        particleData[1].charge = -5.1966000e-01;
        std::vector<RealVec> positions(numberOfParticles);
        positions[0]             = RealVec( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  );
        positions[1]             = RealVec( -1.763651687e+00, -3.816594649e-01, -1.300353949e+00  );

        for (int i=0; i<numberOfParticles; i++) {
            particleData[1].particleIndex = i;
            particleData[i].multipoleAtomZs = -1;
            particleData[i].multipoleAtomXs = -1;
            particleData[i].multipoleAtomYs = -1;
            for (int j=0; j<3; j++) {
                particleData[i].position[j] = positions[i][j] * 1e-1;
            }
        }


        _inducedDipole.push_back(Vec3(-0.128403739, -0.0930144581, -1.4288696));
        _inducedDipole.push_back(Vec3(0.128403739, 0.0930144581, 1.4288696));
        for (int i=0; i<numberOfParticles; i++) {
                _inducedDipole[i] *= 1e-1;
                _inducedDipolePolar.push_back(_inducedDipole[i]);
        }

        std::vector<RealVec> forces;
        std::vector<RealOpenMM> electrostaticPotential;

        for (int i=0; i<numberOfParticles; i++) {
            forces.push_back(Vec3( 0.,  0., 0.  ));
            electrostaticPotential.push_back(0.);
        }
        RealOpenMM energy = calculatePmeDirectElectrostaticPairIxn( particleData, 0, 1,forces, electrostaticPotential);

        std::vector<RealVec> expectedForces;
        expectedForces.push_back(Vec3(3.11046, 2.25319, 34.6131));
        expectedForces.push_back(Vec3(-3.11046, -2.25319, -34.6131));

        double tolerance = 1e-5;

        for( unsigned int ii = 0; ii < forces.size(); ii++ ){
            forces[ii] /= (cal2joule * 10);
//            ASSERT_EQUAL_VEC_MOD( expectedForces[ii], forces[ii], tolerance, testName );
            std::cout << expectedForces[ii] << " Kcal/mol/A <expected>" << std::endl;
            std::cout << forces[ii] << " Kcal/mol/A" << std::endl;
        }
        energy /= cal2joule;
//        ASSERT_EQUAL_TOL_MOD( 0.0634441, energy, tolerance, testName );
        std::cout << "Energy: " << energy << " Kcal/mol "<< std::endl;
        std::cout << "Energy: " << 0.0634441 << " Kcal/mol <expected>"<< std::endl;
        std::cout << "Test END: " << testName << std::endl << std::endl;

    }
};

int main( int numberOfArguments, char* argv[] ) {

    try {
        std::cout << "TestReferenceMBPolElectrostaticsForce running test..." << std::endl;

        testGetAndScaleInverseRs();
        testGetAndScaleInverseRsInterMulecolar();

        WrappedMBPolReferenceElectrostaticsForceForIndDipole* mbpolReferenceElectrostaticsForce = new WrappedMBPolReferenceElectrostaticsForceForIndDipole();
        mbpolReferenceElectrostaticsForce->setMutualInducedDipoleTargetEpsilon(1e-7);
        mbpolReferenceElectrostaticsForce->wrapCalculateInducedDipolePairIxns();

        WrappedMBPolReferenceElectrostaticsForceForPmeDipole* mbpolReferenceElectrostaticsForcePme = new WrappedMBPolReferenceElectrostaticsForceForPmeDipole();
        mbpolReferenceElectrostaticsForcePme->setMutualInducedDipoleTargetEpsilon(1e-7);
        mbpolReferenceElectrostaticsForcePme->setCutoffDistance( 10. );
        mbpolReferenceElectrostaticsForcePme->setAlphaEwald( 1e-15 );
        std::vector<int> pmeGrid(3);
        std::fill(pmeGrid.begin(), pmeGrid.end(), 20.);
        mbpolReferenceElectrostaticsForcePme->setPmeGridDimensions(pmeGrid);
        RealVec boxSize;
        boxSize[0] = boxSize[1] = boxSize[2] = 50;
        mbpolReferenceElectrostaticsForcePme->setPeriodicBoxSize(boxSize);

        mbpolReferenceElectrostaticsForcePme->wrapCalculateInducedDipolePairIxns();
//
        WrappedMBPolReferenceElectrostaticsForceForCalculateElectrostaticPairIxn* wrapperForComputeElectrostaticPairIxn = new WrappedMBPolReferenceElectrostaticsForceForCalculateElectrostaticPairIxn();
        wrapperForComputeElectrostaticPairIxn->testCalculateElectrostaticPairIxn();

        WrappedMBPolReferenceElectrostaticsForceForComputeWaterCharge* wrapperForComputeWaterCharge = new WrappedMBPolReferenceElectrostaticsForceForComputeWaterCharge();
        wrapperForComputeWaterCharge->testComputeWaterCharge();

        testWater3();

        testWater3VirtualSite();

        WrappedMBPolReferencePmeElectrostaticsForceForcalculatePmeDirectElectrostaticPairIxn* mbpolReferenceElectrostaticsForcePmePair = new WrappedMBPolReferencePmeElectrostaticsForceForcalculatePmeDirectElectrostaticPairIxn();
        mbpolReferenceElectrostaticsForcePmePair->setMutualInducedDipoleTargetEpsilon(1e-7);
        mbpolReferenceElectrostaticsForcePmePair->setCutoffDistance( 10. );
        mbpolReferenceElectrostaticsForcePmePair->setAlphaEwald( 1e-15 );
        // std::vector<int> pmeGrid(3);
        std::fill(pmeGrid.begin(), pmeGrid.end(), 20.);
        mbpolReferenceElectrostaticsForcePmePair->setPmeGridDimensions(pmeGrid);
        // RealVec boxSize;
        boxSize[0] = boxSize[1] = boxSize[2] = 50;
        mbpolReferenceElectrostaticsForcePmePair->setPeriodicBoxSize(boxSize);

        mbpolReferenceElectrostaticsForcePmePair->testCalculateElectrostaticPairIxn();

        testWater3VirtualSitePMEHugeBox();

        testWater3VirtualSitePMESmallBox();

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
