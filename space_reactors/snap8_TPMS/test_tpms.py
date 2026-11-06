import openmc
tpms_period = 2.0
fuel = openmc.Material(name='fuel'); fuel.add_nuclide('U235', 1.0); fuel.set_density('g/cm3', 10)
clad = openmc.Material(name='clad'); clad.add_nuclide('Fe56', 1.0); clad.set_density('g/cm3', 8)
gap  = openmc.Material(name='gap');  gap.add_nuclide('He4', 1.0);   gap.set_density('g/cm3', 1e-4)
cool = openmc.Material(name='cool'); cool.add_nuclide('Na23', 1.0); cool.set_density('g/cm3', 0.8)

cool_frac, clad_frac, fuel_frac = 0.6326, 0.1182, 0.2342
a = cool_frac/2; b = a + clad_frac/2; c = 0.5 - fuel_frac/2; d = 0.5 + fuel_frac/2
tc1, tg1, tf1, tf2, tg2, tc2 = openmc.TPMS.from_relative_densities(
    "Gyroid", [a, b, c, d, 1-b, 1-a], tpms_period)

box = openmc.model.RectangularParallelepiped(-3, 3, -3, 3, -3, 3, boundary_type='reflective')
cells = [
    openmc.Cell(fill=fuel, region=+tf1 & -tf2 & -box),
    openmc.Cell(fill=clad, region=+tc1 & -tg1 & -box),
    openmc.Cell(fill=clad, region=+tg2 & -tc2 & -box),
    openmc.Cell(fill=gap,  region=+tg1 & -tf1 & -box),
    openmc.Cell(fill=gap,  region=+tf2 & -tg2 & -box),
    openmc.Cell(fill=cool, region=(-tc1 | +tc2) & -box),
]
model = openmc.Model(openmc.Geometry(cells), openmc.Materials([fuel, clad, gap, cool]))
model.settings.particles = 100; model.settings.batches = 2
model.export_to_xml()

p = openmc.Plot(); p.basis = 'xy'; p.origin = (0, 0, 0); p.width = (6, 6); p.pixels = (600, 600)
p.color_by = 'material'; p.colors = {fuel: 'red', clad: 'grey', gap: 'yellow', cool: 'blue'}
p.filename = 'tpms_test'
p.to_ipython_image()
