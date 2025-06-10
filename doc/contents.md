# Code Contents {#doc-contents}

[TOC]

# Introduction for new users {#intro_new_users}

* [Installation](installation.md) How to install MPI-AMRVAC.
* [Getting Started](getting_started.md) How to run your first test problem.
* [Features](features.md) An overview of the main features of MPI-AMRVAC.
* [Acknowledgments](acknowledgments.md) Information on collaboration and
financial support.
* [FAQ](faq.md) Frequently asked questions.
* [**Tests**](test.html) An overview of documented test cases.
* @ref demo-movies.md
* @ref contributing.md
* [Documentation](documentation.md) How the documentation works.

# General {#general}

* [Equations](equations.md) The equations and parameters in physics modules.
* [User files](amrvacusr.md) How to create a new problem, specify initial
  conditions and customize functionalities.
* [Parameters](par.md) Description of all parameters in "amrvac.par" parameter file.
* [Auxiliary variables](mpiamrvac_nw.md) Description of the intended use
  for _nw, nwflux, nwaux, nwextra, nwauxio_ parameters.
* [Command line](commandline.md) Help on command-line parameters.
* [Examples](examples.md) Description of various example simulations for which
  parameter files and user modules have been provided.

# Discretization methods and AMR strategy {#discretization}

* [Discretization](discretization.md) The equation and its discretization, the
basic variables, the structure of the grid, boundaries, ghost cells.
* [Time Discretization](time_discretization.md) The way the time integration works.
* [Methods](methods.md) Properties of the discretization methods like TVDLF,
TVDMU, TVD, HLL...
* [Slope limiters](limiter.md) Slope limiters for reconstruction to suppress spurious numerical oscillations
* [AMR aspects](amrstructure.md) Some essential info on global parameters and
the data structures for the block-tree AMR.
* [Using polar/cylindrical/spherical coordinates](axial.md) Some information on simulations using non-Cartesian grids.
* [Small Values](smallvalues.md) Info on positivity fixes for (M)HD runs.

# Additional Physics {#special_sources}

* [Thermal conduction](thermal_conduction.md) Description of solving thermal conduction.
* [Radiative cooling](radiative_cooling.md) Description of adding radiative cooling.
* [Pressure less fluids (dust)](dust.md) Getting started with the dust module.
* [Test particles in (M)HD](particle.md) Description of test particle tracing routines.
* [Reaction-diffusion equations](reaction_diffusion.md) Getting started with the reaction-diffusion module.
* [Advection-Reaction-diffusion equations](advection_reaction_diffusion.md) Using the advection-reaction-diffusion module.
* [Two-fluid equations](twofluid.md) Getting started with the two-fluid module.
* [CAK radiation force](cakforce.md) Getting started with the CAK force for line-driven wind outflows.
* [Adding a new physics module](addmodule.md) Description of how to add your own physics module.

# Source Code {#source_code}

* [Source](source.md) Description of the dimension independent source language,
which is translated to F90 by the VACPP preprocessor.
* [Variable Names](varnames.md) How variable names are formed in the source
files.
* [VACPP](vacpp.md) Making and running the VACPP preprocessor itself.
* [Compilation](compilation.md) Info on compilation, debugging and adding libraries.

# IO and data analysis {#io}

* [File format](fileformat.md) Description of the format of MPI-AMRVAC data file (.dat).
* [Visualisation and analysis in Python with yt](yt_usage.md) without data conversion.
* [Converting data files for visualization](convert.md) Brief notes on how to
convert to Tecplot (.plt), and VTK (.vtu) data files.
* [Slices](slices.md) How to output hypersurfaces (slices) for restart or
visualization.
* [Line of sight views](collapsed.md) How to output collapsed views for
visualisation and analysis (e.g. column densities).
* [Analysis routine](analysis.md) Using the run-time analysis routine.
* [3D Printing](print3D.md) A note on how to generate 3D printed results.
* [Yt visualization](yt_usage.md) The recommended yt usage for visualization.
