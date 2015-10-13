from __future__ import print_function

import unittest
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
import mbpol

class TestReferenceMBPolElectrostatics(unittest.TestCase):
    """This tests the Reference implementation of ReferenceMBPolOneBodyForce."""

    def testWater3VirtualSite(self, nonbondedMethod=app.NoCutoff):
         expected_energy = -15.818784
         pdb = app.PDBFile("pdb_files/water3_electrostatics.pdb")
         forcefield = app.ForceField("../mbpol_no_custom_dispersion.xml") 
         # this xml file is the same as mbpol.xml but the Custom dispersion force is never added to the system.  This is needed because it cannot be removed from the system once it is added and thereroe prevents the testing of individual forces within mbpol.xml
         forcefield._forces = forcefield._forces[0:1]
         nonbondedCutoff = .7*unit.nanometer
         if (nonbondedMethod == app.PME):
             # expected_energy = -13.0493 # with no polarizability
             expected_energy = -15.818784 # with polarizability
             nonbondedCutoff = 10*unit.nanometer
             boxDimension = 50
             boxsize = [boxDimension, boxDimension, boxDimension]
             pdb.topology.setUnitCellDimensions( boxsize )

         system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff)
         integrator = mm.LangevinIntegrator(0.0, 0.1, 0.01)
         platform = mm.Platform.getPlatformByName('Reference')
         simulation = app.Simulation(pdb.topology, system, integrator, platform)
         simulation.context.setPositions(pdb.positions)
         simulation.context.computeVirtualSites()
         state = simulation.context.getState(getForces=True, getEnergy=True, getPositions=True)
         potential_energy = state.getPotentialEnergy()
         potential_energy.in_units_of(unit.kilocalorie_per_mole)

         print(potential_energy.in_units_of(unit.kilocalorie_per_mole)._value)

         self.assertTrue(abs(potential_energy.in_units_of(unit.kilocalorie_per_mole)._value - expected_energy) < .1)

    def testWater3VirtualSitePMEHugeBox(self):
         self.testWater3VirtualSite(nonbondedMethod=app.PME)

    def testWater3(self):
        nonbondedCutoff = .7*unit.nanometer
        nonbondedMethod=app.NoCutoff
        #charge dirstribution set in the xml
        pdb = app.PDBFile("../water3_novirtualsite.pdb")
        forcefield = app.ForceField("../mbpol_no_custom_dispersion_no_charge_redistribution.xml")
        # this xml file is the same as mbpol.xml but the Custom dispersion force is never added to the system.  This is needed because it cannot be removed from the system once it is added and therefore prevents the testing of individual forces within mbpol.xml
        forcefield._forces = forcefield._forces[0:1]
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff)
        integrator = mm.LangevinIntegrator(0.0, 0.1, 0.01)
        platform = mm.Platform.getPlatformByName('Reference')
        simulation = app.Simulation(pdb.topology, system, integrator, platform)
        simulation.context.setPositions(pdb.positions)
        state = simulation.context.getState(getForces=True, getEnergy=True, getPositions=True)
        potential_energy = state.getPotentialEnergy()
        print("Energy:")
        print(potential_energy.in_units_of(unit.kilocalorie_per_mole))
        forces = state.getForces()
        for force in forces:
            print(force*0.0239005736)
        expected_energy = -6.98374
        self.assertTrue(abs(potential_energy.in_units_of(unit.kilocalorie_per_mole)._value - expected_energy) < .1)
if __name__ == '__main__':
    unittest.main()
