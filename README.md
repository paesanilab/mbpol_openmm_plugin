MB-pol plugin for OpenMM
=======================

[![Build Status](https://travis-ci.org/paesanilab/mbpol_openmm_plugin.svg?branch=master)](https://travis-ci.org/paesanilab/mbpol_openmm_plugin)
[![Anaconda Badge](https://anaconda.org/paesanilab/mbpol/badges/version.svg)](https://anaconda.org/paesanilab/mbpol)

`mbpol` is a plugin for the `OpenMM` toolkit for molecular simulations using the `MB-pol` many-body water potential. Currently, the `mbpol` plugin is available within the “reference implementation” of the `OpenMM` toolkit, which is appropriate for energy calculations for both clusters and systems in periodic boundary conditions. Multicore and GPU implementations for actual molecular dynamics production runs within `OpenMM` will become available in the near future.

The development of `MB-pol` started with a detailed analysis of the two- and three-body water interactions evaluated at the CCSD(T) level to quantitatively assess the accuracy of current force fields, DFT models, and ab initio based interaction potentials that are commonly used in molecular simulations [1]. On the basis of this analysis, the full-dimensional MB-pol  potential was developed entirely from "first principles" building upon the many-body expansion of the interaction energy of water molecules [2,3]. 

`MB-pol` explicitly treats the one-body (intramolecular distortion energy) term and the short-ranged two- and three-body terms. `MB-pol` can thus be viewed as a classical polarizable potential supplemented by short-range two- and three-body terms that effectively represent quantum-mechanical interactions arising from the overlap of the monomer electron densities. Specifically, at all separations, the total MB-pol two-body term includes (damped) dispersion forces derived from ab initio computed asymptotic expansions of the dispersion energy along with electrostatic contributions due to the interactions between the molecular permanent and induced moments. At short-range, this two-body term is supplemented by a 4th-degree permutationally invariant polynomial that smoothly switches to zero as the oxygen-oxygen separation in the dimer approaches 6.5 Angstrom. Similarly, the `MB-pol` three-body term includes a three-body polarization term at all separations, which is supplemented by a short-range 4th-degree permutationally invariant polynomial that effectively corrects for the deficiencies of a purely classical representation of the three-body interactions in regions where the electron densities of the three monomers overlap. This short-range three-body contribution is smoothly switched off once the oxygen-oxygen separation between any water molecule and the other two water molecules of a trimer reaches a value of 4.5 Angstrom. In `MB-pol`, all induced interactions are described through many-body polarization. `MB-pol` thus contains many-body effects at all monomer separations as well as at all orders, in an explicit way up to the third order and in a mean-field fashion at all higher orders. 

Without containing any empirical parameters, `MB-pol` accurately describes the properties of gas-phase clusters, including the dimer vibration-rotation tunneling spectrum  [2], the second and third virial coefficients [2,3], and cluster structures and energies [3]. Simulations carried out with path-integral molecular dynamics (PIMD) and centroid molecular dynamics (CMD) demonstrate that MB-pol provides a highly accurate description of the liquid phase of water at ambient conditions in comparison with experiment for several structural, thermodynamic, and dynamical properties [4]. Many-body molecular dynamics (MB-MD) simulations carried out with `MB-pol` in combination with many-body representations of the dipole moment (`MB-mu`) and polarizability (`MB-alpha  `), built upon the convergence of the many-body expansion of the electrostatic properties of water [5], predict infrared (IR) and Raman spectra of liquid water in agreement with the experimental results [6].


References:

1. [G.R. Medders, V. Babin, and F. Paesani, J. Chem. Theory Comput. 9, 1103–1114 (2013)](http://pubs.acs.org/doi/abs/10.1021/ct300913g).
2. [V. Babin, C. Leforestier, and F. Paesani, J. Chem. Theory Comput. 9, 5395–5403 (2013)](http://pubs.acs.org/doi/abs/10.1021/ct400863t).
3. [V. Babin, G.R. Medders, and F. Paesani, J. Chem. Theory Comput. 10, 1599–1607 (2014)](http://pubs.acs.org/doi/abs/10.1021/ct500079y).
4. [G.R. Medders, V. Babin, and F. Paesani, J. Chem. Theory Comput. 10, 2906–2910 (2014)](http://pubs.acs.org/doi/abs/10.1021/ct5004115).
5. [G.R. Medders, and F. Paesani, J. Chem. Theory Comput. 9, 4844–4852 (2013)](http://pubs.acs.org/doi/abs/10.1021/ct400696d).
6. [G.R. Medders, and F. Paesani, J. Chem. Theory Comput. 11, 1145–1154 (2015)](http://pubs.acs.org/doi/abs/10.1021/ct501131j).


available components are:

* `MBPolReferenceElectrostaticsForce`
* `MBPolReferenceOneBodyForce`
* `MBPolReferenceTwoBodyForce`
* `MBPolReferenceThreeBodyForce`
* `MBPolReferenceDispersionForce`

The parameters of each component are defined in [`python/mbpol.xml`](https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/python/mbpol.xml).
As of version `1.0`, only the `Reference` platform, i.e. single threaded C++ on CPU, is supported. It currently can simulate clusters of water molecules and water systems in periodic boundary conditions using with particle mesh Ewald (PME).

Included with the mbpol plugin, there are new potential energy functions (i-TTM) describing the interactions between halide ions / alkali metal ions and water molecules, and also ion / ion potential energy surfaces. The i-TTM potentials are derived from fits to electronic structure data and include an explicit treatment of two-body repulsion, electrostatics, and dispersion energy. Many-body effects are represented through classical polarization within an extended Thole-type model. By construction, the i-TTM potentials are compatible with the flexible and fully ab initio MB-pol potential, which has recently been shown to accurately predict the properties of water from the gas to the condensed phase.

The i-TTM potential is expressed as V(i-TTM) = V(electrostatics) + V(disp) + V(rep). The first term [V(electrostatics)] is the electrostatic energy (charge and dipole interactions), which is calculated in the same way as in MBPol. The second term [V(disp)] is the dispersion energy, which has the form V(disp) = f(delta,R) (C6/R^6), where f(delta,R) is the tang-toennies damping function, C6 the dispersion coefficient for that atom pair, and delta is the damping parameter that depends on the pair of atoms, same as in MBPol. The last term [V(rep)] is the repulsive energy, also defined as a sum of pairwise interactions between the different atoms, and can be expressed as V(rep) = A exp(-bR), where A and b are fitting parameters that depend on the pair of atoms. As specified in the literature, it has been shown that b and delta should be set to be the same for the same pair of atoms.

These new potential energy surfaces can be used to run bulk simulations and study the diffusion of the different ions at different concentration of the salt CX in water, being C the cation and X the halide, One can also study spectroscopy of the bulk (Infrared spechtroscopy, IR) or the air/water interface (Sum Frequency Generation, SFG) as a function of the salt concentration.


References:

1. [D. Arismendi-Arrieta, M. Riera, P. Bajaj, R. Prosmiti, and F. Paesani, J. Phys. Chem. B. 120, 1822–1832 (2016)](http://pubs.acs.org/doi/abs/10.1021/acs.jpcb.5b09562)



## How to install

See the [`INSTALL.md`](https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/INSTALL.md) file in the package.

## Automated script builder

`OpenMM` uses a framework approach, users need to create a Python script for setting up and running simulations, the easiest way is to create a script with the online configuration tool <http://builder.openmm.org>, however the builder does not support `MBPol`.
In order to simplify running `MBPol` for new `OpenMM` users we created a tool, `mbpol_builder` that loads a configuration file and outputs a Python script ready to be run.

The script covers these use-cases:

* Cluster simulations: energy computation, geometry optimization and constant energy simulations
* Bulk simulations: cubic box Periodic Mesh Ewald energy computation, constant energy, constant temperature and constant pressure simulations

How to use `mbpol_builder`:

* Create a new folder for the simulation
* Copy a `pdb` format initial positions file in the folder, this plugin provides a sample `water14_cluster.pdb` and `water256_bulk.pdb` in the `python/` folder
* Copy a sample configuration file in the folder, the plugin provides `mbpol_config.ini`
* Customize the configuration file parameters, see the comments
* run `mbpol_builder mbpol_config.ini generated_script_filename.py`
* run `python generated_script_filename.py` to run the simulation

## Example simulation

Simulation of a cluster of 14 water molecules:

* [IPython Notebook example on `nbviewer`](http://nbviewer.ipython.org/gist/zonca/54c7040c1cf3f583930f)
* IPython Notebook example source: <https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/python/water14.ipynb>
* Python Script: <https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/python/water14.py>
* C++ example:
  <https://github.com/paesanilab/mbpol_openmm_plugin/blob/master/examples/Simulate14WaterCluster/Simulate14WaterCluster.cpp>, see the documentation in the `examples/` folder.

## Support

* For installation issues or software bugs open an issue on Github: <https://github.com/paesanilab/mbpol_openmm_plugin/issues>
* For other questions email the main developer Andrea Zonca (username: zonca, domain: sdsc.edu)
