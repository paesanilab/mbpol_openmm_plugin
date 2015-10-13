import mbpolplugin
from simtk.openmm import app
from simtk import unit

try:
    from exceptions import ValueError, NotImplementedError
except:
    pass

## @private
class MBPolOneBodyForceGenerator:

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []

    @staticmethod
    def parseElement(element, forceField):
        generator = MBPolOneBodyForceGenerator()
        forceField.registerGenerator(generator)

        # <MBPolOneBodyForce>
        #     <Residue class1="OW" class2="HW" class3="HW" />
        # </MBPolOneBodyForce>

        for MBPolOneBodyForce_template in element.findall('MBPolOneBodyForce'):
            types = forceField._findAtomTypes(MBPolOneBodyForce_template.attrib, 3)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

            else:
                outputString = self.__class__ + ": error getting types: %s %s %s" % (
                                    MBPolOneBodyForce_template.attrib['class1'],
                                    MBPolOneBodyForce_template.attrib['class2'],
                                    MBPolOneBodyForce_template.attrib['class3'])
                raise ValueError(outputString)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        methodMap = {app.NoCutoff:mbpolplugin.MBPolOneBodyForce.NonPeriodic,
                     app.PME:mbpolplugin.MBPolOneBodyForce.Periodic,
                     app.CutoffPeriodic:mbpolplugin.MBPolOneBodyForce.Periodic,
                     app.CutoffNonPeriodic:mbpolplugin.MBPolOneBodyForce.NonPeriodic}

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mbpolplugin.MBPolOneBodyForce]

        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for MBPolOneBodyForce')

        if len(existing) == 0:
            force = mbpolplugin.MBPolOneBodyForce()
            sys.addForce(force)
        else:
            force = existing[0]

        force.setNonbondedMethod(methodMap[nonbondedMethod])

        for residue in data.atoms[0].residue.chain.residues():
            global_atoms_indices = [a.index for a in residue.atoms()]

            if len(global_atoms_indices) > 1:
                global_atoms_indices += [-1] * (3 - len(global_atoms_indices))
                v = mbpolplugin.vectori()
                for i in range(3):
                    v.push_back(global_atoms_indices[i])

                force.addOneBody(v)

app.forcefield.parsers["MBPolOneBodyForce"] = MBPolOneBodyForceGenerator.parseElement

## @private
class MBPolTwoBodyForceGenerator:

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []

    @staticmethod
    def parseElement(element, forceField):
        generator = MBPolTwoBodyForceGenerator()
        forceField.registerGenerator(generator)

        # <MBPolTwoBodyForce>
        #     <Residue class1="OW" class2="HW" class3="HW" />
        # </MBPolTwoBodyForce>

        for MBPolTwoBodyForce_template in element.findall('MBPolTwoBodyForce'):
            types = forceField._findAtomTypes(MBPolTwoBodyForce_template, 3)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

            else:
                outputString = self.__class__ + ": error getting types: %s %s %s" % (
                                    MBPolTwoBodyForce_template.attrib['class1'],
                                    MBPolTwoBodyForce_template.attrib['class2'],
                                    MBPolTwoBodyForce_template.attrib['class3'])
                raise ValueError(outputString)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        methodMap = {app.NoCutoff:mbpolplugin.MBPolTwoBodyForce.NoCutoff,
                     app.PME:mbpolplugin.MBPolTwoBodyForce.CutoffPeriodic,
                     app.CutoffPeriodic:mbpolplugin.MBPolTwoBodyForce.CutoffPeriodic,
                     app.CutoffNonPeriodic:mbpolplugin.MBPolTwoBodyForce.CutoffNonPeriodic}

        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for MBPolTwoBodyForce')

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mbpolplugin.MBPolTwoBodyForce]

        if len(existing) == 0:
            force = mbpolplugin.MBPolTwoBodyForce()
            force.setCutoff(float(nonbondedCutoff.value_in_unit(unit.nanometer)))
            sys.addForce(force)
        else:
            force = existing[0]

        force.setNonbondedMethod(methodMap[nonbondedMethod])

        for residue in data.atoms[0].residue.chain.residues():
            global_atoms_indices = [a.index for a in residue.atoms()]

            global_atoms_indices += [-1] * (3 - len(global_atoms_indices))
            v = mbpolplugin.vectori()
            for i in range(3):
                v.push_back(global_atoms_indices[i])

            force.addParticle(v)

app.forcefield.parsers["MBPolTwoBodyForce"] = MBPolTwoBodyForceGenerator.parseElement

## @private
class MBPolThreeBodyForceGenerator:

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []

    @staticmethod
    def parseElement(element, forceField):
        generator = MBPolThreeBodyForceGenerator()
        forceField.registerGenerator(generator)

        # <MBPolThreeBodyForce>
        #     <Residue class1="OW" class2="HW" class3="HW" />
        # </MBPolThreeBodyForce>

        for MBPolThreeBodyForce_template in element.findall('MBPolThreeBodyForce'):
            types = forceField._findAtomTypes(MBPolThreeBodyForce_template.attrib, 3)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

            else:
                outputString = self.__class__ + ": error getting types: %s %s %s" % (
                                    MBPolThreeBodyForce_template.attrib['class1'],
                                    MBPolThreeBodyForce_template.attrib['class2'],
                                    MBPolThreeBodyForce_template.attrib['class3'])
                raise ValueError(outputString)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        methodMap = {app.NoCutoff:mbpolplugin.MBPolThreeBodyForce.NoCutoff,
                     app.PME:mbpolplugin.MBPolThreeBodyForce.CutoffPeriodic,
                     app.CutoffPeriodic:mbpolplugin.MBPolThreeBodyForce.CutoffPeriodic,
                     app.CutoffNonPeriodic:mbpolplugin.MBPolThreeBodyForce.CutoffNonPeriodic}

        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for MBPolThreeBodyForce')

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mbpolplugin.MBPolThreeBodyForce]

        if len(existing) == 0:
            force = mbpolplugin.MBPolThreeBodyForce()
            force.setCutoff(float(nonbondedCutoff.value_in_unit(unit.nanometer)))
            sys.addForce(force)
        else:
            force = existing[0]

        force.setNonbondedMethod(methodMap[nonbondedMethod])

        for residue in data.atoms[0].residue.chain.residues():
            global_atoms_indices = [a.index for a in residue.atoms()]

            global_atoms_indices += [-1] * (3 - len(global_atoms_indices))
            v = mbpolplugin.vectori()
            for i in range(3):
                v.push_back(global_atoms_indices[i])

            force.addParticle(v)


app.forcefield.parsers["MBPolThreeBodyForce"] = MBPolThreeBodyForceGenerator.parseElement

## @private
class MBPolDispersionForceGenerator:

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.parameters = dict()

    @staticmethod
    def parseElement(element, forceField):
        generator = MBPolDispersionForceGenerator()
        forceField.registerGenerator(generator)

        # <MBPolDispersionForce>
        #     <Residue name="HOH" class1="O" class2="H" class3="H" />
        #     <Parameters class1="O" class2="O" c6="9.92951990e+08" d6="9.29548582e+01" />
        #     <Parameters class1="O" class2="H" c6="3.49345451e+08" d6="9.77520243e+01" />
        #     <Parameters class1="H" class2="H" c6="8.40715638e+07" d6="9.40647517e+01" />
        # </MBPolDispersionForce>

        for residue in element.findall('Residue'):
            types = forceField._findAtomTypes(residue.attrib, 3)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

            else:
                outputString = self.__class__ + ": error getting types: %s %s %s" % (
                                    residue.attrib['class1'],
                                    residue.attrib['class2'],
                                    residue.attrib['class3'])
                raise ValueError(outputString)

        for parameter in element.findall('Parameters'):
            generator.parameters[(parameter.attrib['class1'], parameter.attrib['class2'])] = \
                (float(parameter.attrib['c6']), float(parameter.attrib['d6']))

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        methodMap = {app.NoCutoff:mbpolplugin.MBPolDispersionForce.NoCutoff,
                     app.PME:mbpolplugin.MBPolDispersionForce.CutoffPeriodic,
                     app.CutoffPeriodic:mbpolplugin.MBPolDispersionForce.CutoffPeriodic,
                     app.CutoffNonPeriodic:mbpolplugin.MBPolDispersionForce.CutoffNonPeriodic}

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mbpolplugin.MBPolDispersionForce]

        if len(existing) == 0:
            force = mbpolplugin.MBPolDispersionForce()
            force.setCutoff(float(nonbondedCutoff.value_in_unit(unit.nanometer)))
            sys.addForce(force)
        else:
            force = existing[0]

        force.setNonbondedMethod(methodMap[nonbondedMethod])

        for atom in data.atoms:
            if atom.element:
                force.addParticle(atom.element.symbol)
            else:
                # virtual site
                force.addParticle('M')

        for elements, c6d6 in self.parameters.items():
            force.addDispersionParameters(elements[0], elements[1], c6d6[0], c6d6[1])

app.forcefield.parsers["MBPolDispersionForce"] = MBPolDispersionForceGenerator.parseElement

## @private
class MBPolElectrostaticsForceGenerator:

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.typeMap = {}
        self.setIncludeChargeRedistribution = True;

    @staticmethod
    def parseElement(element, forceField):

        # <MBPolElectrostaticsForce>
        #     <Residue class1="OW" class2="HW" class3="HW" thole-charge-charge="0.4" thole-charge-dipole="0.4" thole-dipole-dipole-intermolecules="0.055" thole-dipole-dipole-1-2="0.055" thole-dipole-dipole-1-3="0.626" thole-dipole-dipole-2-3="0.626" /> 
        #     <Atom type="MBPol-O" charge="-5.1966000e-01" damping-factor="0.00131" polarizability="0.00131" />
        #     <Atom type="MBPol-H" charge="2.5983000e-01" damping-factor="0.000294" polarizability="0.000294" />
        #     <Atom type="MBPol-M" charge="0" damping-factor="0.00131" polarizability="0" />
        # <ChargeRedistribution set="False"/>
        # </MBPolElectrostaticsForce>
        generator = MBPolElectrostaticsForceGenerator()
        forceField.registerGenerator(generator)
        # mbpolElectrostaticsForce->setIncludeChargeRedistribution(false);
        for MBPolElectrostaticsForce_template in element.findall('MBPolElectrostaticsForce'):
            types = forceField._findAtomTypes(MBPolElectrostaticsForce_template.attrib, 3)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

            else:
                outputString = self.__class__ + ": error getting types: %s %s %s" % (
                                    MBPolElectrostaticsForce_template.attrib['class1'],
                                    MBPolElectrostaticsForce_template.attrib['class2'],
                                    MBPolElectrostaticsForce_template.attrib['class3'])
                raise ValueError(outputString)

        for residue in element.findall('Residue'):
            
        #     <Residue name="HOH" class1="OW" class2="HW" class3="HW" thole-charge-charge="0.4" thole-charge-dipole="0.4" thole-dipole-dipole-intermolecules="0.055" thole-dipole-dipole-1-2="0.055" thole-dipole-dipole-1-3="0.626" thole-dipole-dipole-2-3="0.626" /> 

            name = residue.attrib["name"]

            v = mbpolplugin.vectord()
            thole_components = ["thole-charge-charge", "thole-charge-dipole", "thole-dipole-dipole-intermolecules", "thole-dipole-dipole-OH", "thole-dipole-dipole-HH", ]
            for each in thole_components:
                v.push_back(float(residue.attrib[each]))

            generator.typeMap[name] = dict(thole = v)

        for atom in element.findall('Atom'):
        #     <Atom type="MBPol-H" charge="2.5983000e-01" damping-factor="0.000294" polarizability="0.000294" />
            types = forceField._findAtomTypes(atom.attrib, 1)
            if None not in types:
    
                for t in types[0]:
                    generator.typeMap[t] = dict(  polarizability = float(atom.attrib['polarizability']),
                                             charge         = float(atom.attrib['charge']),
                                             damping_factor = float(atom.attrib['damping-factor']) )

            else:
                outputString = "MBPolElectrostaticsForceGenerator: error getting type for atom: %s" % (atom.attrib['type'])
                raise ValueError(outputString)
                
        for value in element.findall('setIncludeChargeRedistribution'):
            if value.attrib['chargeRedistribution'] == 'True':
                generator.setIncludeChargeRedistribution = True
            elif value.attrib['chargeRedistribution'] == 'False':
                generator.setIncludeChargeRedistribution = False
         
    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        # CutoffNonPeriodic defaults to NoCutoff
        methodMap = {app.NoCutoff:mbpolplugin.MBPolElectrostaticsForce.NoCutoff,
                     app.CutoffNonPeriodic:mbpolplugin.MBPolElectrostaticsForce.NoCutoff,
                     app.PME:mbpolplugin.MBPolElectrostaticsForce.PME}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for MBPolElectrostaticsForce')

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mbpolplugin.MBPolElectrostaticsForce]

        if len(existing) == 0:
            force = mbpolplugin.MBPolElectrostaticsForce()
            force.setCutoffDistance(float(nonbondedCutoff.value_in_unit(unit.nanometer)))
            #force.setIncludeChargeRedistribution(False)
            sys.addForce(force)
        else:
            force = existing[0]
        
        force.setIncludeChargeRedistribution(self.setIncludeChargeRedistribution)
        force.setNonbondedMethod(methodMap[nonbondedMethod])

        for atom in data.atoms:
            t = data.atomType[atom]
            global_atoms_indices = [a.index for a in atom.residue.atoms() if a.index != atom.index]

            global_atoms_indices += [-1] * (3 - len(global_atoms_indices))

            force.addElectrostatics(self.typeMap[t]['charge'], global_atoms_indices[0], global_atoms_indices[1], global_atoms_indices[2], self.typeMap[atom.residue.name]['thole'], self.typeMap[t]['damping_factor'], self.typeMap[t]['polarizability'])

app.forcefield.parsers["MBPolElectrostaticsForce"] = MBPolElectrostaticsForceGenerator.parseElement
