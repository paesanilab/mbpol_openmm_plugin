from __future__ import print_function

import unittest
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
import mbpol

class TestReferenceMBPolIntegration(unittest.TestCase):
    def test_water3(self, nonbondedMethod=app.CutoffNonPeriodic):
        expected_energy = -8.78893485
        pdb = app.PDBFile("../water3.pdb")
        forcefield = app.ForceField("../mbpol.xml")
        nonbondedCutoff = 0.9*unit.nanometers
        if (nonbondedMethod == app.PME):
            expected_energy = -8.92353
            boxDimension = 1.9
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
        self.assertTrue(abs(potential_energy.in_units_of(unit.kilocalorie_per_mole)._value - expected_energy) < .1)

    def test_water3_periodic(self):
        self.test_water3(nonbondedMethod=app.PME)

    ###def test_water50_periodic(self):
    ###    nonbondedMethod=app.PME
    ###    expected_energy = -244.37507
    ###    pdb = app.PDBFile("pdb_files/water50.pdb")
    ###    forcefield = app.ForceField("../mbpol.xml")
    ###    nonbondedCutoff = 0.9*unit.nanometer
    ###    
    ###    if (nonbondedMethod == app.PME):
    ###        boxDimension = 1.8
    ###        boxsize = [boxDimension, boxDimension, boxDimension]
    ###        pdb.topology.setUnitCellDimensions( boxsize )
    ###        
    ###    system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff)
    ###        
    ###    integrator = mm.LangevinIntegrator(0.0, 0.1, 0.01)
    ###    platform = mm.Platform.getPlatformByName('Reference')
    ###    simulation = app.Simulation(pdb.topology, system, integrator, platform)
    ###    simulation.context.setPositions(pdb.positions)
    ###    simulation.context.computeVirtualSites()
    ###    state = simulation.context.getState(getForces=True, getEnergy=True, getPositions=True)
    ###    potential_energy = state.getPotentialEnergy()
    ###    potential_energy.in_units_of(unit.kilocalorie_per_mole)
    ###    
    ###    #print(potential_energy.in_units_of(unit.kilocalorie_per_mole)._value)
    ###    self.assertTrue(abs(potential_energy.in_units_of(unit.kilocalorie_per_mole)._value - expected_energy) < 1)
    ###    
    ###def test_water256_periodic(self):
    ###    nonbondedMethod=app.PME
    ###    expected_energy = -2270.88890
    ###    pdb = app.PDBFile("pdb_files/water256_integration_test.pdb")
    ###    forcefield = app.ForceField("../mbpol.xml")

    ###    nonbondedCutoff = 0.9*unit.nanometer
    ###    
    ###    if (nonbondedMethod == app.PME):
    ###        boxDimension = 19.3996888399961804/10.
    ###        boxsize = [boxDimension, boxDimension, boxDimension]
    ###        pdb.topology.setUnitCellDimensions( boxsize )
    ###        
    ###    system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff)
    ###        
    ###    integrator = mm.LangevinIntegrator(0.0, 0.1, 0.01)
    ###    platform = mm.Platform.getPlatformByName('Reference')
    ###    simulation = app.Simulation(pdb.topology, system, integrator, platform)
    ###    simulation.context.setPositions(pdb.positions)
    ###    simulation.context.computeVirtualSites()
    ###    state = simulation.context.getState(getForces=True, getEnergy=True, getPositions=True)
    ###    potential_energy = state.getPotentialEnergy()
    ###    potential_energy.in_units_of(unit.kilocalorie_per_mole)
    ###    
    ###    #print(potential_energy.in_units_of(unit.kilocalorie_per_mole)._value)
    ###    
    ###    
    ###    self.assertTrue(abs(potential_energy.in_units_of(unit.kilocalorie_per_mole)._value - expected_energy) < 20)

if __name__ == '__main__':
    unittest.main()
