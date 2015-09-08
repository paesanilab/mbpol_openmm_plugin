#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "OpenMMMBPol.h"
#include "openmm/System.h"
#include "openmm/VirtualSite.h"
#include "openmm/LocalEnergyMinimizer.h"

#include "openmm/LangevinIntegrator.h"
#include "openmm/VerletIntegrator.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

#define ASSERT_EQUAL_TOL_MOD(expected, found, tol, testname) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << testname << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC_MOD(expected, found, tol,testname) {ASSERT_EQUAL_TOL_MOD((expected)[0], (found)[0], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[1], (found)[1], (tol),(testname)); ASSERT_EQUAL_TOL_MOD((expected)[2], (found)[2], (tol),(testname));};

using namespace  OpenMM;
using namespace MBPolPlugin;

const double TOL = 1e-4;
const double cal2joule = 4.184;

// Handy homebrew PDB writer for quick-and-dirty trajectory output.
void writePdbFrame(std::string filenametag, int frameNum, double time_ps, const OpenMM::State& state, unsigned int particlesPerMolecule)
{
    // Open file
    FILE * pFile;
    char buf[4];
    sprintf( buf, "%04d", frameNum );
    std::string zeroPaddedFrameNum = buf;
    std::string filename = "simulate14WaterCluster_" + filenametag + "_" + zeroPaddedFrameNum + ".pdb";
    pFile = fopen (filename.c_str(),"w");

    // std::cout << "Write Pdb: " << filename << std::endl;

    // Reference atomic positions in the OpenMM State.
    const std::vector<OpenMM::Vec3>& pos_nm = state.getPositions();

    const char* atomNames[] = {" O  ", " H1 ", " H2 ", " VS "}; // cycle through these
    fprintf(pFile, "MODEL     %d\n", frameNum);
    fprintf(pFile, "REMARK 250 time=%.3f picoseconds\n", time_ps);
    int atomCount = 1;
    for (int atom=0; atom < (int)pos_nm.size(); ++atom)
    {
//        if (atom%4 != 3)
//        {
//            fprintf(pFile, "HETATM%5d %4s HOH  %4d    ",        // start of pdb HETATM line
//                atomCount, atomNames[atom%4], 1 + (atomCount-1)/3); // atom number, name, residue #
//            fprintf(pFile, "%8.3f%8.3f%8.3f",                   // middle of pdb HETATM line
//                pos_nm[atom][0]*10, pos_nm[atom][1]*10, pos_nm[atom][2]*10);
//            fprintf(pFile, "  1.00  0.00            \n");       // end of pdb HETATM line
//            atomCount++;
//        }
        fprintf(pFile, "HETATM%5d %4s HOH  %4d    ",        // start of pdb HETATM line
           atomCount, atomNames[atom%particlesPerMolecule], 0); // atom number, name, residue #
            //atomCount, atomNames[atom%particlesPerMolecule], 1 + (atomCount-1)/3); // atom number, name, residue #
        fprintf(pFile, "%8.3f%8.3f%8.3f",                   // middle of pdb HETATM line
            pos_nm[atom][0]*10, pos_nm[atom][1]*10, pos_nm[atom][2]*10);
        fprintf(pFile, "  1.00  0.00            \n");       // end of pdb HETATM line
        atomCount++;
    }

    fprintf(pFile, "ENDMDL\n"); // end of frame
    fclose (pFile);
}

void simulate14WaterCluster() {

    std::string testName = "simulate14WaterCluster";

    // Electrostatics

    double virtualSiteWeightO = 0.573293118;
    double virtualSiteWeightH = 0.213353441;

    MBPolElectrostaticsForce* mbpolElectrostaticsForce        = new MBPolElectrostaticsForce();;
    mbpolElectrostaticsForce->setNonbondedMethod( MBPolElectrostaticsForce::NoCutoff );

    std::vector<double> thole(5);

    thole[TCC] = 0.4;
    thole[TCD] = 0.4;
    thole[TDD] = 0.055;
    thole[TDDOH]  = 0.626;
    thole[TDDHH] = 0.055;

    // One body interaction
    MBPolOneBodyForce* mbpolOneBodyForce = new MBPolOneBodyForce();


    // Two body interaction
    MBPolTwoBodyForce* mbpolTwoBodyForce = new MBPolTwoBodyForce();
    double cutoff = 6.5 / 10.; // 6.5 A => nm
    cutoff = 1e10;
    mbpolTwoBodyForce->setCutoff( cutoff );
    mbpolTwoBodyForce->setNonbondedMethod(MBPolTwoBodyForce::CutoffNonPeriodic);

    // Three body interaction
    MBPolThreeBodyForce* mbpolThreeBodyForce = new MBPolThreeBodyForce();
    //cutoff = 5. / 10.; // 4.5 A => nm
    cutoff = 1e10; // 4.5 A => nm
    mbpolThreeBodyForce->setCutoff( cutoff );
    mbpolThreeBodyForce->setNonbondedMethod(MBPolThreeBodyForce::CutoffNonPeriodic);

    // Atom positions O H H, [A], no virtual sites

    std::vector<Vec3> positions;

    // CONFIGURATION OPTIONS
    bool useWater3 = true;
    unsigned int particlesPerMolecule = 4;

    // restore initial positions for water14, so that we can run constant energy
    // simulation even if energy minimization fails
    bool restoreInitialPositionsAfterMinimization = not useWater3;
    
    if (useWater3) {
        positions.push_back(Vec3( -1.516074336e+00, -2.023167650e-01,  1.454672917e+00  ));
        positions.push_back(Vec3( -6.218989773e-01, -6.009430735e-01,  1.572437625e+00  ));
        positions.push_back(Vec3( -2.017613812e+00, -4.190350349e-01,  2.239642849e+00  ));
        positions.push_back(Vec3( -1.763651687e+00, -3.816594649e-01, -1.300353949e+00  ));
        positions.push_back(Vec3( -1.903851736e+00, -4.935677617e-01, -3.457810126e-01  ));
        positions.push_back(Vec3( -2.527904158e+00, -7.613550077e-01, -1.733803676e+00  ));
        positions.push_back(Vec3( -5.588472140e-01,  2.006699172e+00, -1.392786582e-01  ));
        positions.push_back(Vec3( -9.411558180e-01,  1.541226676e+00,  6.163293071e-01  ));
        positions.push_back(Vec3( -9.858551734e-01,  1.567124294e+00, -8.830970941e-01  ));
    } else {
        // water14
        // 1
        positions.push_back(Vec3( -2.349377641189e-01,  1.798934467398e-01,  1.896881820756e-01 ));
        positions.push_back(Vec3(  1.788456192811e-01, -4.351349633402e-01, -3.765224894244e-01 ));
        positions.push_back(Vec3(  2.195852756811e-01,  8.816778044978e-02,  1.073177788976e+00 ));
        // 2
        positions.push_back(Vec3( -2.899375600289e+00,  4.533801552398e-01,  4.445704119756e-01 ));
        positions.push_back(Vec3( -1.891275615534e+00,  3.753590916398e-01,  2.490296312756e-01 ));
        positions.push_back(Vec3( -3.372062885819e+00,  2.555551963398e-01, -3.733422905244e-01 ));
        // 3
        positions.push_back(Vec3(  8.428220990811e-01,  4.870807363398e-01,  3.155605857676e+00 ));
        positions.push_back(Vec3(  2.205240762811e-01,  1.158785963140e+00,  3.577474762076e+00 ));
        positions.push_back(Vec3(  1.087801046381e+00, -1.129208591502e-01,  3.875459300076e+00 ));
        // 4
        positions.push_back(Vec3(  9.587606425811e-01,  1.739780448140e+00, -2.350038144714e+00 ));
        positions.push_back(Vec3(  1.608290357181e+00,  1.084550067240e+00, -2.527800037751e+00 ));
        positions.push_back(Vec3(  6.371470514811e-01,  1.907637977740e+00, -1.459391248124e+00 ));
        // 5
        positions.push_back(Vec3( -5.402595366189e-01,  3.436858949240e+00, -3.364294798244e-01 ));
        positions.push_back(Vec3( -9.927075362189e-01,  4.233122552140e+00, -8.308930347244e-01 ));
        positions.push_back(Vec3( -1.277172354239e+00,  2.836232717540e+00, -3.850631209244e-01 ));
        // 6
        positions.push_back(Vec3(  8.739387783811e-01, -1.626958848900e+00, -2.396798577774e+00 ));
        positions.push_back(Vec3(  1.408213796381e+00, -2.079240526160e+00, -1.701669540424e+00 ));
        positions.push_back(Vec3(  1.515572976181e+00, -1.166885340300e+00, -2.938074529504e+00 ));
        // 7
        positions.push_back(Vec3(  2.976994446581e+00, -4.267919297502e-01,  1.269092190476e+00 ));
        positions.push_back(Vec3(  2.533682331181e+00,  2.556735979978e-02,  2.039154286576e+00 ));
        positions.push_back(Vec3(  3.128503830381e+00,  2.419881502398e-01,  5.967620551756e-01 ));
        // 8
        positions.push_back(Vec3( -1.629810395359e+00,  8.289601626398e-01, -3.056166050894e+00 ));
        positions.push_back(Vec3( -6.994519097189e-01,  8.049764115398e-01, -2.943175234164e+00 ));
        positions.push_back(Vec3( -1.874145275955e+00, -8.466235253022e-02, -3.241518244384e+00 ));
        // 9
        positions.push_back(Vec3( -2.595625033599e+00, -1.646832204690e+00,  2.192150970676e+00 ));
        positions.push_back(Vec3( -2.906585612949e+00, -9.062350233512e-01,  1.670815817176e+00 ));
        positions.push_back(Vec3( -2.669295291189e+00, -1.431610358210e+00,  3.153621682076e+00 ));
        // 10
        positions.push_back(Vec3(  1.818457619381e+00, -2.794187137660e+00,  1.107580890756e-01 ));
        positions.push_back(Vec3(  2.097391780281e+00, -1.959930744260e+00,  5.631032745756e-01 ));
        positions.push_back(Vec3(  2.529698452281e+00, -3.454934434160e+00,  3.013241934756e-01 ));
        // 11
        positions.push_back(Vec3( -4.545107775189e-01, -3.318742298960e+00,  1.492762161876e+00 ));
        positions.push_back(Vec3(  1.737151477811e-01, -3.171411334960e+00,  8.066594968756e-01 ));
        positions.push_back(Vec3( -1.041679032169e+00, -2.555367771060e+00,  1.600821086376e+00 ));
        // 12
        positions.push_back(Vec3( -5.833794722189e-01,  3.052540693140e+00,  2.637054256676e+00 ));
        positions.push_back(Vec3( -1.419205322489e+00,  2.586471523640e+00,  2.720532465476e+00 ));
        positions.push_back(Vec3( -5.746468928189e-01,  3.408128243540e+00,  1.785534806076e+00 ));
        // 13
        positions.push_back(Vec3(  3.344764010581e+00,  1.359764839140e+00, -6.682810128244e-01 ));
        positions.push_back(Vec3(  2.813836912981e+00,  2.168029298140e+00, -5.595663217244e-01 ));
        positions.push_back(Vec3(  4.237583121481e+00,  1.647152701540e+00, -8.700689423244e-01 ));
        // 14
        positions.push_back(Vec3( -2.004856348130e+00, -1.781589543180e+00, -2.722539611524e+00 ));
        positions.push_back(Vec3( -1.060302432829e+00, -1.994959751260e+00, -2.496572811216e+00 ));
        positions.push_back(Vec3( -2.596156174389e+00, -2.566367260560e+00, -2.697712987614e+00 ));
    }

    // Setup system

    int numberOfWaterMolecules = positions.size()/3; // virtual sites are added later, here
    // positions has only 3 atoms per water molecule
    
    int numberOfParticles          = numberOfWaterMolecules * particlesPerMolecule;

    if (particlesPerMolecule < 4) {
        mbpolElectrostaticsForce->setIncludeChargeRedistribution(false);
    }
    std::vector<int> particleIndices(particlesPerMolecule);
    System system;
    for( unsigned int jj = 0; jj < numberOfParticles; jj += particlesPerMolecule ){
        system.addParticle( 1.5999000e+01 ); // Mass
        system.addParticle( 1.0080000e+00 );
        system.addParticle( 1.0080000e+00 );

        particleIndices[0] = jj;
        particleIndices[1] = jj+1;
        particleIndices[2] = jj+2;

        // Charge, dipoles and quadrupoles (zero in MBPol)
        mbpolElectrostaticsForce->addElectrostatics( -5.1966000e-01, jj+1, jj+2, jj+3,
                                            thole, 0.001310, 0.001310 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+2, jj+3,
                                            thole, 0.000294, 0.000294 );
        mbpolElectrostaticsForce->addElectrostatics(  2.5983000e-01, jj, jj+1, jj+3,
                                            thole, 0.000294, 0.000294 );
        if (particlesPerMolecule == 4) {
            system.addParticle( 0. ); // Virtual Site
            system.setVirtualSite(jj+3, new ThreeParticleAverageSite(jj, jj+1, jj+2,
                                                                   virtualSiteWeightO, virtualSiteWeightH,virtualSiteWeightH));
            mbpolElectrostaticsForce->addElectrostatics(  0., jj, jj+1, jj+2,
                                                        thole,  0.001310,  0.);
        }


        mbpolOneBodyForce->addOneBody(particleIndices);
        mbpolTwoBodyForce->addParticle( particleIndices);
        mbpolThreeBodyForce->addParticle( particleIndices);
        dispersionForce->addParticle( "O");
        dispersionForce->addParticle( "H");
        dispersionForce->addParticle( "H");
        dispersionForce->addParticle( "M");

    }

    // <!-- Units: c6 [kJ mol^{-1} nm^{-6}], d6 [nm^{-1}] -->
    dispersionForce->addDispersionParameters("O", "O", 9.92951990e+08, 9.29548582e+01);
    dispersionForce->addDispersionParameters("O", "H", 3.49345451e+08, 9.77520243e+01);
    dispersionForce->addDispersionParameters("H", "H", 8.40715638e+07, 9.40647517e+01);

    system.addForce(mbpolElectrostaticsForce);
    system.addForce(mbpolOneBodyForce);
    system.addForce(mbpolTwoBodyForce);
    system.addForce(mbpolThreeBodyForce);
    system.addForce(dispersionForce);

    if (particlesPerMolecule == 4) {
        // insert null position for virtual site
        std::vector<Vec3>::iterator beginOfNextWater;
        for (int w=0; w<numberOfWaterMolecules; w++) {
                beginOfNextWater = positions.begin() + w * particlesPerMolecule + (particlesPerMolecule-1);
                positions.insert(beginOfNextWater, Vec3( 0., 0., 0.));
        }
    }
    // A => nm
    for (int i=0; i<numberOfParticles; i++) {
            positions[i] *= 1e-1;
    }

    State state;
    std::vector<Vec3> forces;
    double energy;

    // Setup context
    std::string platformName = "Reference";
    double temperature_K = 100;
    double frictionCoeff_invps = 0.1;
    double stepSize_ps = 0.2 * 1e-3; // 0.2 femtoseconds

    // Simulate at constant energy
    
    VerletIntegrator constantEnergyIntegrator(stepSize_ps);
    Context constantEnergyContext(system, constantEnergyIntegrator, Platform::getPlatformByName( platformName ) );
    
    // Local minimization of energy
    double minimizationTolerance = 5;

    constantEnergyContext.setPositions(positions);
    constantEnergyContext.applyConstraints(1e-6); // update position of virtual sites

    state    = constantEnergyContext.getState(State::Positions | State::Energy);
    energy = state.getPotentialEnergy();
    writePdbFrame("BeforeMin", 0, 0, state, particlesPerMolecule); // output coordinates

    std::cout << "Potential energy before energy minimization: " << energy/cal2joule << " Kcal/mol" << std::endl;

    // Local energy minimization
    LocalEnergyMinimizer::minimize(constantEnergyContext, minimizationTolerance);

    state    = constantEnergyContext.getState(State::Positions | State::Energy);
    energy = state.getPotentialEnergy(); 
    writePdbFrame("AfterMin", 0, 0, state, particlesPerMolecule); // output coordinates
    std::cout << "Potential energy after energy minimization: " << energy/cal2joule << " Kcal/mol" << std::endl << std::endl;

    std::cout << "Constant Energy simulation:" <<  std::endl;

    if (restoreInitialPositionsAfterMinimization) {
        std::cout << "Reset position to initial state" << std::endl;

        constantEnergyContext.setPositions(positions);
        constantEnergyContext.applyConstraints(1e-6); // update position of virtual sites
    }

    // set velocity randomly given temperature
    std::cout << "setVelocitiesToTemperature: " << temperature_K << "K" << std::endl;
    constantEnergyContext.setVelocitiesToTemperature(temperature_K);

    for (int frameNum=1; ;++frameNum) {
        
        OpenMM::State state    = constantEnergyContext.getState(State::Positions | State::Energy);
        energy = state.getPotentialEnergy() + state.getKineticEnergy();
        double potentialEnergy = state.getPotentialEnergy();
        const double  time_ps = state.getTime();
        std::cout << "Time (ps): " << std::setw (6) << time_ps << " | Energy: " << std::setw(9) << energy/cal2joule << " Kcal/mol | Potential Energy: " << std::setw(9) << potentialEnergy/cal2joule << " Kcal/mol" << std::endl;
        writePdbFrame("Econst", frameNum, time_ps, state, particlesPerMolecule); // output coordinates

        if (frameNum >= 10)
        // if (time_ps >= 100.)
            break;

        // Advance state many steps at a time, for efficient use of OpenMM.
        // 500 steps = .1 picoseconds
        constantEnergyIntegrator.step(10);
    }
    
}

int main( int numberOfArguments, char* argv[] ) {

    try {
        std::cout << "simulate14WaterCluster" << std::endl;

        simulate14WaterCluster();

    } catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
