from __future__ import print_function

import unittest
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
import mbpol

class TestReferenceMBPolIntegrationWithCl(unittest.TestCase):

    def test_water3andCl(self):
        expected_energy = -47.3516 + 1.81488
        nonbondedMethod = app.CutoffNonPeriodic
        pdb = app.PDBFile("pdb_files/water3AndCl.pdb")
        ###forcefield = app.ForceField("../mbpol_no_custom_dispersion.xml")
        forcefield = app.ForceField("../mbpol.xml")
        nonbondedCutoff = 10*unit.nanometers
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff)
        integrator = mm.LangevinIntegrator(0.0, 0.1, 0.01)

        forces = system.getForces()[:-1]
        for num, force in enumerate(forces):
            force.setForceGroup(num)

        platform = mm.Platform.getPlatformByName('Reference')
        simulation = app.Simulation(pdb.topology, system, integrator, platform)
        simulation.context.setPositions(pdb.positions)
        simulation.context.computeVirtualSites()

        state = simulation.context.getState(getForces=True, getEnergy=True, getPositions=True)
        potential_energy = state.getPotentialEnergy()
        potential_energy.in_units_of(unit.kilocalorie_per_mole)

        total = potential_energy.in_units_of(unit.kilocalorie_per_mole)

        ### Needs to use no dispersion force in the XML to run this
        ### force_labels = ["electrostatics", "onebody", "twobody", "threebody"]
        ### for num, force in enumerate(forces):
        ###     state = simulation.context.getState(getForces=True, getEnergy=True, groups=2**num)
        ###     potential_energy = state.getPotentialEnergy()
        ###     print(force_labels[num] + ":")
        ###     print(potential_energy.in_units_of(unit.kilocalorie_per_mole))
        print("Total:")
        print(total)

        self.assertTrue(abs(potential_energy.in_units_of(unit.kilocalorie_per_mole)._value - expected_energy) < .1)

if __name__ == '__main__':
    unittest.main()
