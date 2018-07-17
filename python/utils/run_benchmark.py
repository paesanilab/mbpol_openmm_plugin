from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
import sys
import os.path
import pandas as pd
import datetime
import mbpol

configurations = [("mbpol", "CPU")]
boundaries = ["cluster"]
num_molecules_list = [512]

pdb_filenames = { 512: "pdb/H2O_1",
                  256: "pdb/TIP5P_PIMC_0C_LIQ.0"
                }

steps = 1
timestep = 0.02*unit.femtoseconds
temperature = 300*unit.kelvin

for num_molecules in num_molecules_list:
    for boundary in boundaries:
        for forcefield_configuration, platform_name  in configurations:

            pdb_filename = pdb_filenames[num_molecules]
            if forcefield_configuration.startswith("amoeba"):
                pdb_filename += "_AMOEBA"
            pdb_filename += ".pdb"
            pdb = app.PDBFile(pdb_filename)

            tag = os.path.basename(pdb_filename).replace(".pdb", "")

            if boundary == "PME":
                # PME box size
                if tag.startswith("H2O"):
                    boxSize = 24.8343697144598003
                elif tag.startswith("FM_25C") or tag.startswith("TIP5P_25C"):
                    boxSize = 19.7316565863235596
                else:
                    boxSize = 19.3996888399961804
                boxSize /= 10.
                pdb.topology.setUnitCellDimensions((boxSize,boxSize,boxSize))
                nonbondedMethod = app.PME

            elif boundary == "cluster":

                nonbondedMethod = app.CutoffNonPeriodic if forcefield_configuration == "mbpol" else app.NoCutoff

            forcefield = app.ForceField(forcefield_configuration + ".xml")

            system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbondedMethod, nonbondedCutoff=.9*unit.nanometer)
            integrator = mm.VerletIntegrator(timestep)

            platform = mm.Platform.getPlatformByName(platform_name)
            simulation = app.Simulation(pdb.topology, system, integrator, platform)
            simulation.context.setPositions(pdb.positions)
            simulation.context.computeVirtualSites()

            #state = simulation.context.getState(getForces=True, getEnergy=True)
            simulation.context.setVelocitiesToTemperature(temperature)

            #simulation.reporters.append(app.StateDataReporter("simulation_{tag}_{timestep}_{temperature}.log".format(**locals()).replace(" ",""), 100, step=True, totalEnergy=True,
            #potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
            #speed=True, totalSteps=steps, separator=','))

            start = datetime.datetime.now()
            simulation.step(steps)
            end = datetime.datetime.now()
            print("%s, %s, %s, %d, %.0f" % (forcefield_configuration, platform_name, boundary, num_molecules, (end-start).total_seconds()))
