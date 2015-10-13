from __future__ import print_function

import unittest
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
import mbpol

class TestCustomForce(unittest.TestCase):
    """This tests the Reference implementation of ReferenceMBPolOneBodyForce."""

    def testOneOneBody(self):
        nonbondedMethod=app.PME
        expected_energy = 0.55975882
        pdb = app.PDBFile("pdb_files/water1.pdb")
        forcefield = app.ForceField("../mbpol.xml")
        forcefield._forces = forcefield._forces[1:2]
        nonbondedCutoff = 0.9*unit.nanometer
        
        if (nonbondedMethod == app.PME):
            boxDimension = 1.8
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
        
        #print(potential_energy.in_units_of(unit.kilocalorie_per_mole)._value)
        
        
        self.assertTrue(abs(potential_energy.in_units_of(unit.kilocalorie_per_mole)._value - expected_energy) < .001)
        
if __name__ == '__main__':
    unittest.main()