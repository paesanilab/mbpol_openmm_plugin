from __future__ import print_function

import unittest
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
import mbpol
from prepare_xml import prepare_xml

class TestCustomForce(unittest.TestCase):
    """Test the functionality of V(i-TTM) = V(electrostatics) + V(disp) + V(rep) xml file."""

    def test_one(self, pdb_file="./pdb_files/na_f.pdb", expected_energy=-163.00791517):
        pdb = app.PDBFile(pdb_file)
        nonbondedMethod=app.CutoffNonPeriodic  
        #forcefield = app.ForceField("../i-TTM_50.xml")
        prepare_xml(percent_pol=.5)
        forcefield = app.ForceField("../i-TTM_template.xml")
        nonbondedCutoff = 1e3*unit.nanometer
        if (nonbondedMethod == app.CutoffPeriodic):
            boxsize = [50, 50, 50]
            pdb.topology.setUnitCellDimensions( boxsize )
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbondedMethod, nonBondedCutoff=nonbondedCutoff)
      
        integrator = mm.VerletIntegrator(0.02*unit.femtoseconds)
        platform = mm.Platform.getPlatformByName('Reference')
        simulation = app.Simulation(pdb.topology, system, integrator, platform)
        simulation.context.setPositions(pdb.positions)
        simulation.context.computeVirtualSites()
        state = simulation.context.getState(getForces=True, getEnergy=True, getPositions=True)
        potential_energy = state.getPotentialEnergy()
        potential_energy.in_units_of(unit.kilocalorie_per_mole)
        print("@ {} calculated energy = {}  expected energy = {}".format(pdb_file, potential_energy.in_units_of(unit.kilocalorie_per_mole)._value, expected_energy))
        
        
        self.assertTrue(abs(potential_energy.in_units_of(unit.kilocalorie_per_mole)._value - expected_energy) < .1)
    def test_two(self):
        self.test_one(pdb_file="./pdb_files/water_iodide.pdb", expected_energy=110.32196257)
       
if __name__ == '__main__':
    unittest.main()
