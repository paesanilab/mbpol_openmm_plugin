from __future__ import print_function

import unittest
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
import mbpol

class TestCustomForce(unittest.TestCase):
    """Test the functionality of Custom Dispersion Force xml file."""

    def test_water_br(self, pdb_file="./pdb_files/water_br.pdb", expected_energy=2.20583105):
        pdb = app.PDBFile(pdb_file)
        nonbondedMethod=app.CutoffNonPeriodic  
        #mbpol with i-ttm is just mbpol with no dispersion and ion info right now
        forcefield = app.ForceField("../i-TTM_elec.xml")
        nonbondedCutoff = 1e3*unit.nanometer
        #electrostatics only
        forcefield._forces = [forcefield._forces[0]] 

        if (nonbondedMethod == app.CutoffPeriodic):
            boxsize = [50, 50, 50]
            pdb.topology.setUnitCellDimensions( boxsize )
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbondedMethod, nonBondedCutoff=nonbondedCutoff)
#       for i in range(0, system.getNumForces()):
#           print(system.getForce(i))
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
    def test_water_br1(self):
        self.test_water_br(pdb_file="./pdb_files/water_br21.pdb", expected_energy=2.75913036)
    def test_Cl_Na(self):
        self.test_water_br(pdb_file="./pdb_files/cl_na.pdb", expected_energy=-222.76239105)
    def test_I_Li(self):
        self.test_water_br(pdb_file="./pdb_files/i_li.pdb", expected_energy=-272.02766332)
       
if __name__ == '__main__':
    unittest.main()