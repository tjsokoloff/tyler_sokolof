#!/usr/bin/env python3
#activate TPMS-env
from math import log10

import numpy as np
import openmc
import sys
import glob

# Dakota files
params_file = sys.argv[1]
results_file = sys.argv[2]

## debug
params_file = sys.argv[1]
results_file = sys.argv[2]

print("PARAMETER FILE:", params_file)

with open(params_file, 'r') as f:
    contents = f.read()

print(contents)

# Read fuel radius from Dakota
with open(params_file, 'r') as f:
    f.readline()
    r = float(f.readline().split()[0])

openmc.config['cross_sections'] = "/home/tjsok/TPMS/openmc-workshop/endfb-viii.1-hdf5/cross_sections.xml"


###############################################################################
# Create materials for the problem

uo2 = openmc.Material(name='UO2 fuel at 2.4% wt enrichment')
uo2.set_density('g/cm3', 10.29769)
uo2.add_element('U', 1., enrichment=2.4)
uo2.add_element('O', 2.)

helium = openmc.Material(name='Helium for gap')
helium.set_density('g/cm3', 0.001598)
helium.add_element('He', 2.4044e-4)

zircaloy = openmc.Material(name='Zircaloy 4')
zircaloy.set_density('g/cm3', 6.55)
zircaloy.add_element('Sn', 0.014  , 'wo')
zircaloy.add_element('Fe', 0.00165, 'wo')
zircaloy.add_element('Cr', 0.001  , 'wo')
zircaloy.add_element('Zr', 0.98335, 'wo')

borated_water = openmc.Material(name='Borated water')
borated_water.set_density('g/cm3', 0.740582)
borated_water.add_element('B', 4.0e-5)
borated_water.add_element('H', 5.0e-2)
borated_water.add_element('O', 2.4e-2)
borated_water.add_s_alpha_beta('c_H_in_H2O')

# Collect the materials together and export to XML
materials = openmc.Materials([uo2, helium, zircaloy, borated_water])
materials.export_to_xml()

###############################################################################
# Define problem geometry

# Create cylindrical surfaces
fuel_or_radius = 0.2

print("Fuel radius =", fuel_or_radius)

fuel_or = openmc.ZCylinder(r=fuel_or_radius, name='Fuel OR')
clad_ir = openmc.ZCylinder(r=fuel_or_radius, name='Clad IR')
clad_or = openmc.ZCylinder(r=fuel_or_radius + 0.25, name='Clad OR')

# Create a region represented as the inside of a rectangular prism
pitch = 2.5
box = openmc.model.RectangularPrism(pitch, pitch, boundary_type='reflective')

# Create cells, mapping materials to regions
fuel = openmc.Cell(fill=uo2, region=-fuel_or)
gap = openmc.Cell(fill=helium, region=+fuel_or & -clad_ir)
clad = openmc.Cell(fill=zircaloy, region=+clad_ir & -clad_or)
water = openmc.Cell(fill=borated_water, region=+clad_or & -box)

# Create a geometry and export to XML
geometry = openmc.Geometry([fuel, gap, clad, water])
geometry.export_to_xml()

###############################################################################
# Define problem settings

# Indicate how many particles to run
settings = openmc.Settings()
settings.batches = 200
settings.inactive = 20
settings.particles = 400000

# Create an initial uniform spatial source distribution over fissionable zones
lower_left = (-pitch/2, -pitch/2, -1)
upper_right = (pitch/2, pitch/2, 1)
uniform_dist = openmc.stats.Box(lower_left, upper_right)
settings.source = openmc.IndependentSource(
    space=uniform_dist, constraints={'fissionable': True})

# For source convergence checks, add a mesh that can be used to calculate the
# Shannon entropy
entropy_mesh = openmc.RegularMesh()
entropy_mesh.lower_left = (-fuel_or.r, -fuel_or.r)
entropy_mesh.upper_right = (fuel_or.r, fuel_or.r)
entropy_mesh.dimension = (10, 10)
settings.entropy_mesh = entropy_mesh
settings.export_to_xml()

###############################################################################
# Define tallies

# Create a mesh that will be used for tallying
mesh = openmc.RegularMesh()
mesh.dimension = (100, 100)
mesh.lower_left = (-pitch/2, -pitch/2)
mesh.upper_right = (pitch/2, pitch/2)

# Create a mesh filter that can be used in a tally
mesh_filter = openmc.MeshFilter(mesh)

# Now use the mesh filter in a tally and indicate what scores are desired
mesh_tally = openmc.Tally(name="Mesh tally")
mesh_tally.filters = [mesh_filter]
mesh_tally.scores = ['flux', 'fission', 'nu-fission']

# Let's also create a tally to get the flux energy spectrum. We start by
# creating an energy filter
e_min, e_max = 1e-5, 20.0e6
groups = 500
energies = np.logspace(log10(e_min), log10(e_max), groups + 1)
energy_filter = openmc.EnergyFilter(energies)

spectrum_tally = openmc.Tally(name="Flux spectrum")
spectrum_tally.filters = [energy_filter]
spectrum_tally.scores = ['flux']

# Instantiate a Tallies collection and export to XML
tallies = openmc.Tallies([mesh_tally, spectrum_tally])
tallies.export_to_xml()

###############################################################################
# Create geometry plot

plot = openmc.Plot()
plot.filename = 'pin_cell'
plot.origin = (0.0, 0.0, 0.0)
plot.width = (pitch, pitch)
plot.pixels = (1000, 1000)
plot.color_by = 'material'

plot.colors = {
    uo2: 'orange',
    helium: 'white',
    zircaloy: 'gray',
    borated_water: 'blue'
}

plots = openmc.Plots([plot])
plots.export_to_xml()

#print(f"Fuel radius = {fuel_or.r:.5f} cm")
# Generate plot
#openmc.plot_geometry()

#openmc.run()

#print(f"Fuel radius = {fuel_or.r:.5f} cm")


############################################################################################
#sp_file = openmc.run()
#with openmc.StatePoint(sp_file) as sp:
#    keff = sp.keff.nominal_value
openmc.run()

sp_file = sorted(glob.glob("statepoint.*.h5"))[-1]

with openmc.StatePoint(sp_file) as sp:
	keff = sp.keff.nominal_value

objective = -keff

with open(results_file, 'w') as f:
    f.write(f"{objective}\n")

#print(f"Radius = {fuel_radius:.5f} cm")
#print(f"keff = {keff:.6f}")
