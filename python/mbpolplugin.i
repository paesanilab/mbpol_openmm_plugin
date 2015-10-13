%module mbpolplugin

%include "std_string.i"
%include "typemaps.i"
%include "std_vector.i"
namespace std {
  %template(vectord) vector<double>;
  %template(vectorddd) vector< vector< vector<double> > >;
  %template(vectori) vector<int>;
  %template(vectorii) vector < vector<int> >;
  %template(vectorstring) vector<string>;
};

%import(module="simtk.openmm") "OpenMMSwigHeaders.i"

%{
#include "openmm/MBPolElectrostaticsForce.h"
#include "openmm/MBPolOneBodyForce.h"
#include "openmm/MBPolTwoBodyForce.h"
#include "openmm/MBPolThreeBodyForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%feature("autodoc", "1");
%nodefaultctor;


using namespace OpenMM;

namespace MBPolPlugin {

class MBPolElectrostaticsForce : public OpenMM::Force {
public:
    MBPolElectrostaticsForce();

    int getNumElectrostatics() const;

    enum NonbondedMethod { NoCutoff, PME };

    void setNonbondedMethod(NonbondedMethod method);

    double getCutoffDistance(void) const;

    void setCutoffDistance(double distance);

    void setIncludeChargeRedistribution( bool chargeRedistribution );

    bool getIncludeChargeRedistribution( void ) const;

    // double getAEwald() const;

    // void setAEwald(double aewald);

    // int getPmeBSplineOrder() const;

    // void getPmeGridDimensions(std::vector<int>& gridDimension) const;

    // void setPmeGridDimensions(const std::vector<int>& gridDimension);

    int addElectrostatics(double charge, 
                      int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, const std::vector<double>& thole, double dampingFactor, double polarity);

    void getElectrostaticsParameters(int index, double& charge,
                                 int& axisType, int& multipoleAtomZ, int& multipoleAtomX, int& multipoleAtomY, std::vector<double>& thole, double& dampingFactor, double& polarity) const;

    void setElectrostaticsParameters(int index, double charge,
                                 int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, const std::vector<double>& thole, double dampingFactor, double polarity);

    // void setCovalentMap(int index, CovalentType typeId, const std::vector<int>& covalentAtoms);

    // void getCovalentMap(int index, CovalentType typeId, std::vector<int>& covalentAtoms) const;

    // void getCovalentMaps(int index, std::vector < std::vector<int> >& covalentLists) const;

    int getMutualInducedMaxIterations(void) const;

    void setMutualInducedMaxIterations(int inputMutualInducedMaxIterations);

    double getMutualInducedTargetEpsilon(void) const;

    void setMutualInducedTargetEpsilon(double inputMutualInducedTargetEpsilon);

    // double getEwaldErrorTolerance() const;

    // void setEwaldErrorTolerance(double tol);

    void getElectrostaticPotential(const std::vector< Vec3 >& inputGrid,
                                     Context& context, std::vector< double >& outputElectrostaticPotential);

    void updateParametersInContext(Context& context);
};

class MBPolOneBodyForce : public OpenMM::Force {
public:

    MBPolOneBodyForce();

    int getNumOneBodys() const;

    enum NonbondedMethod { Periodic, NonPeriodic };

    void setNonbondedMethod(NonbondedMethod method);
    NonbondedMethod getNonbondedMethod() const;

    int addOneBody(const std::vector<int> & particleIndices);

    void getOneBodyParameters(int particleIndex, std::vector<int>& particleIndices ) const;
    void setOneBodyParameters(int index, std::vector<int>& particleIndices  );

    void updateParametersInContext(Context& context);

};

class MBPolTwoBodyForce : public Force {
public:
    MBPolTwoBodyForce();

    int getNumParticles() const;

    void setParticleParameters(int particleIndex, std::vector<int>& particleIndices);

    void getParticleParameters(int particleIndex, std::vector<int>& particleIndices) const;

    int addParticle(const std::vector<int> & particleIndices);

    int getNumMolecules(void) const;
    void setCutoff(double cutoff);

    double getCutoff(void) const;


    enum NonbondedMethod { NoCutoff, CutoffPeriodic, CutoffNonPeriodic };

    NonbondedMethod getNonbondedMethod() const;
    void setNonbondedMethod(NonbondedMethod method);

    void updateParametersInContext(Context& context);

};

class MBPolThreeBodyForce : public Force {
public:
    MBPolThreeBodyForce();

    int getNumParticles() const;

    void setParticleParameters(int particleIndex, std::vector<int>& particleIndices);

    void getParticleParameters(int particleIndex, std::vector<int>& particleIndices) const;

    int addParticle(const std::vector<int> & particleIndices);

    int getNumMolecules(void) const;
    void setCutoff(double cutoff);

    double getCutoff(void) const;

    enum NonbondedMethod { NoCutoff, CutoffPeriodic, CutoffNonPeriodic };


    NonbondedMethod getNonbondedMethod() const;
    void setNonbondedMethod(NonbondedMethod method);

    void updateParametersInContext(Context& context);
};


} // namespace

%pythoncode %{
  # when we import * from the python module, we only want to import the
  # actual classes, and not the swigregistration methods, which have already
  # been called, and are now unneeded by the user code, and only pollute the
  # namespace
  __all__ = [k for k in locals().keys() if not (k.endswith('_swigregister') or k.startswith('_'))]
%}
