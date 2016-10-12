from __future__ import print_function

import unittest
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
import mbpol

class TestCustomForce(unittest.TestCase):
    """Test the functionality of Custom Dispersion Force xml file."""

    def test_water_br(self, pdb_file="./pdb_files/water_br.pdb", expected_energy=22.20519668):
        pdb = app.PDBFile(pdb_file)
        nonbondedMethod=app.CutoffNonPeriodic  
        forcefield = app.ForceField("../i-TTM_100.xml")
        nonbondedCutoff = 1e3*unit.nanometer
        
        if (nonbondedMethod == app.CutoffPeriodic):
            boxsize = [50, 50, 50]
            pdb.topology.setUnitCellDimensions( boxsize )
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbondedMethod, nonBondedCutoff=nonbondedCutoff)
        system.removeForce(2) # remove dispersion
        system.removeForce(0) # remove electrostatics
        integrator = mm.VerletIntegrator(0.02*unit.femtoseconds)
        platform = mm.Platform.getPlatformByName('Reference')
        simulation = app.Simulation(pdb.topology, system, integrator, platform)
        simulation.context.setPositions(pdb.positions)
        simulation.context.computeVirtualSites()
        state = simulation.context.getState(getForces=True, getEnergy=True, getPositions=True)
        potential_energy = state.getPotentialEnergy()
        print("@ {} calculated energy = {}  expected energy = {}".format(pdb_file, potential_energy.in_units_of(unit.kilocalorie_per_mole)._value, expected_energy))
        self.assertTrue(abs(potential_energy.in_units_of(unit.kilocalorie_per_mole)._value - expected_energy) < .01)
    def test_water_br21(self):
        self.test_water_br(pdb_file="./pdb_files/water_br21.pdb", expected_energy=72.96587861)
    def test_Cl_Na(self):
        self.test_water_br(pdb_file="./pdb_files/cl_na.pdb", expected_energy=122.81085303)
    def test_I_Li(self):
        self.test_water_br(pdb_file="./pdb_files/i_li.pdb", expected_energy=162.97995526)
       
if __name__ == '__main__':
    unittest.main()
