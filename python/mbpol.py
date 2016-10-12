import mbpolplugin
from simtk.openmm import app
from simtk import unit

ATOM_TYPES = {
    "MBPol-O" : 0,
    "MBPol-H" : 1,
    "MBPol-M" : 2,
}

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

        for i in range(len(data.angles)):
            angle = data.angles[i]
            atom1 = data.atoms[angle[0]]
            atom2 = data.atoms[angle[1]]
            atom3 = data.atoms[angle[2]]
            if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH' and atom3.residue.name == 'HOH':
                # FIXME loop through all residues of MBPolOneBodyForce and match their name
                v = mbpolplugin.vectori()
                v.push_back(atom2.index)
                v.push_back(atom1.index)
                v.push_back(atom3.index)
                force.addOneBody(v);

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

        for i in range(len(data.angles)):
            angle = data.angles[i]
            atom1 = data.atoms[angle[0]]
            atom2 = data.atoms[angle[1]]
            atom3 = data.atoms[angle[2]]
            if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH' and atom3.residue.name == 'HOH':
                # FIXME loop through all residues of MBPolTwoBodyForce and match their name
                v = mbpolplugin.vectori()
                v.push_back(atom2.index)
                v.push_back(atom1.index)
                v.push_back(atom3.index)

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

        for i in range(len(data.angles)):
            angle = data.angles[i]
            atom1 = data.atoms[angle[0]]
            atom2 = data.atoms[angle[1]]
            atom3 = data.atoms[angle[2]]
            if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH' and atom3.residue.name == 'HOH':
                # FIXME loop through all residues of MBPolThreeBodyForce and match their name
                v = mbpolplugin.vectori()
                v.push_back(atom2.index)
                v.push_back(atom1.index)
                v.push_back(atom3.index)

                force.addParticle(v)

app.forcefield.parsers["MBPolThreeBodyForce"] = MBPolThreeBodyForceGenerator.parseElement

# @private
class MBPolElectrostaticsForceGenerator:

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.typeMap = {}
        self.thole = []
    @staticmethod
    def parseElement(element, forceField):
        # <MBPolElectrostaticsForce>
        #     <Residue class1="OW" class2="HW" class3="HW" thole-charge-charge="0.4" thole-charge-dipole="0.4" thole-dipole-dipole-intermolecules="0.055" thole-dipole-dipole-1-2="0.055" thole-dipole-dipole-1-3="0.626" thole-dipole-dipole-2-3="0.626" /> 
        #     <Atom type="MBPol-O" charge="-5.1966000e-01" damping-factor="0.00131" polarizability="0.00131" />
        #     <Atom type="MBPol-H" charge="2.5983000e-01" damping-factor="0.000294" polarizability="0.000294" />
        #     <Atom type="MBPol-M" charge="0" damping-factor="0.00131" polarizability="0" />
        # </MBPolElectrostaticsForce>

        generator = MBPolElectrostaticsForceGenerator()
        forceField.registerGenerator(generator)

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

        v = mbpolplugin.vectord()
        thole_components = ['thole-charge-charge', 'thole-charge-dipole', 'thole-dipole-dipole', 'thole-dipole-dipole-singlebond', 'thole-dipole-dipole']
        for each in thole_components:
            v.push_back(float(element.attrib[each]))
        generator.thole = v
        for residue in element.findall('Residue'):
        #<Residue name="HOH" class1="O" class2="H" class3="H" />
            name = residue.attrib["name"]

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
            sys.addForce(force)
        else:
            force = existing[0]
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setTholeParameters(self.thole)

        for i in range(len(data.angles)):
            angle = data.angles[i]
            # FIXME loop through all residues of MBPolElectrostaticsForce and match their name

            # FIXME cheating! get virtual site index by max(otheratom indices)  + 1
            global_atoms_indices = set([angle[local_atom_index] for local_atom_index in [0,1,2]])
            virtual_site_index = max(global_atoms_indices) + 1
            global_atoms_indices.add(virtual_site_index)

            for atomIndex in global_atoms_indices:
                atom = data.atoms[atomIndex]
                t = data.atomType[atom]
                if t in self.typeMap:
                    force.addElectrostatics(self.typeMap[t]['charge'], i, ATOM_TYPES[t], self.typeMap[t]['damping_factor'], self.typeMap[t]['polarizability'])

                else:
                    raise ValueError('No type for atom %s %s %d' % (atom.name, atom.residue.name, atom.residue.index))
        # add non water
        for atom in data.atoms:
            t = data.atomType[atom]
            if t in self.typeMap and t not in ['MBPol-O', 'MBPol-H', 'MBPol-M']:
                force.addElectrostatics(self.typeMap[t]['charge'], -1, -1, -1, self.typeMap['HOH']['thole'], self.typeMap[t]['damping_factor'], self.typeMap[t]['polarizability'])

app.forcefield.parsers["MBPolElectrostaticsForce"] = MBPolElectrostaticsForceGenerator.parseElement
