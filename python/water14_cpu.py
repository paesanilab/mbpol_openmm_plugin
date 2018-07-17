
# In[ ]:

from __future__ import print_function


# In[ ]:

from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
import mbpol


##### Input system in pdb format

# In[ ]:

pdb = app.PDBFile("water14_cluster.pdb")


##### Define the type of potential, first file defines all elements, only the water model is in the second xml file

# In[ ]:

forcefield = app.ForceField(mbpol.__file__.replace('mbpol.py', 'mbpol.xml'))
# use tip4p
#forcefield = app.ForceField("tip4pfb.xml")


##### Create the System, define an integrator, define the Simulation

# In[ ]:

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic, nonBondedCutoff=1e3*unit.nanometer)
integrator = mm.VerletIntegrator(0.00001*unit.femtoseconds)


# In[ ]:

platform = mm.Platform.getPlatformByName('CPU')
simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
simulation.context.computeVirtualSites()


##### Compute initial energy and forces with getState

# In[ ]:

state = simulation.context.getState(getForces=True, getEnergy=True)
potential_energy = state.getPotentialEnergy()
potential_energy.in_units_of(unit.kilocalorie_per_mole)


# In[ ]:

kilocalorie_per_mole_per_angstrom = unit.kilocalorie_per_mole/unit.angstrom
for f in state.getForces():
    print(f.in_units_of(kilocalorie_per_mole_per_angstrom))


##### Local energy minimization

# In[ ]:

from simtk.openmm import LocalEnergyMinimizer


# In[ ]:

#LocalEnergyMinimizer.minimize(simulation.context, 1e-1)


##### Energy, forces and positions after minimization

# In[ ]:

state = simulation.context.getState(getForces=True, getEnergy=True, getPositions=True)
potential_energy = state.getPotentialEnergy()
potential_energy.in_units_of(unit.kilocalorie_per_mole)


# In[ ]:

kilocalorie_per_mole_per_angstrom = unit.kilocalorie_per_mole/unit.angstrom
for f in state.getForces():
    print(f.in_units_of(kilocalorie_per_mole_per_angstrom))


# In[ ]:

state.getPositions()


##### Run a constant energy simulation (Verlet integrator)

# In[ ]:

simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
# Equilibrate
simulation.step(10)


# Add a `reporter` that prints out the simulation status every 10 steps

# In[ ]:

simulation.reporters.append(app.StateDataReporter(sys.stdout, 10, step=True,
    potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
    speed=True, totalSteps=110, separator='\t'))


# Add a `PDBReporter` that writes molecules positions every 20 steps in a pdb file.

# In[ ]:

simulation.reporters.append(app.PDBReporter('trajectory.pdb', 20))


# Run 100 steps

# In[ ]:

simulation.step(100)


# In[ ]:

#get_ipython().system(u'head trajectory.pdb')


# In[ ]:

#get_ipython().system(u'echo Number of lines: `wc -l trajectory.pdb`')
#

# In[ ]:



