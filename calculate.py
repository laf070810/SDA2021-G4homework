from ROOT import *
import numpy as np


scintillator_layer_num = 67
cell_num_on_x_per_layer = 7
cell_num_on_y_per_layer = 7
cell_width_on_x = 40.4
cell_width_on_y = 40.4
scintillator_thickness = 4.
gap_thickness = 2.
Fe_thickness = 1.
plastic_thickness = 7.

calorimeter_size_x = cell_width_on_x * cell_num_on_x_per_layer
calorimeter_size_y = cell_width_on_y * cell_num_on_y_per_layer
calorimeter_size_z = scintillator_layer_num * scintillator_thickness + \
    (scintillator_layer_num - 1) * gap_thickness + 2 * (Fe_thickness + plastic_thickness)
world_size_x = 1.2 * calorimeter_size_x
world_size_y = 1.2 * calorimeter_size_y
world_size_z = 1.2 * calorimeter_size_z
first_position = -calorimeter_size_z / 2


for i in range(1):
    photon_energy = 0.5 * (i + 1)

    f = TFile(f'data/B4_{i}.root')
    tr = f.Get('B4')

    nf = TFile(f'data/E_and_position_{photon_energy}GeV.root', 'recreate')
    ntr = TTree('E_and_position', '')
    buffer = [np.zeros(1) for _ in range(4)]
    ntr.Branch('Esum', buffer[0], 'Esum/D')
    ntr.Branch('shower_x', buffer[1], 'shower_x/D')
    ntr.Branch('shower_y', buffer[2], 'shower_y/D')
    ntr.Branch('shower_z', buffer[3], 'shower_z/D')

    seed_x = 0
    seed_y = 0
    seed_z = 0
    radius = 100

    num = 0
    for event in tr:
        if num % 100 == 0:
            print(f"{num}/{tr.GetEntries()}")
        num += 1

        Esum = 0
        shower_x = 0
        shower_y = 0
        shower_z = 0
        index = -1
        for iz in range(scintillator_layer_num):
            for iy in range(cell_num_on_y_per_layer):
                for ix in range(cell_num_on_x_per_layer):
                    x_position = -calorimeter_size_x / 2 + ix * cell_width_on_x + cell_width_on_x / 2
                    y_position = -calorimeter_size_y / 2 + iy * cell_width_on_y + cell_width_on_y / 2
                    z_position = first_position + Fe_thickness + plastic_thickness + iz * \
                        (scintillator_thickness + gap_thickness) + scintillator_thickness / 2
                    index += 1

                    if ((x_position - seed_x) ** 2 + (y_position - seed_y) ** 2 > radius ** 2):
                        continue

                    Esum += event.Eabs[index]
                    shower_x += x_position * event.Eabs[index]
                    shower_y += y_position * event.Eabs[index]
                    shower_z += z_position * event.Eabs[index]
        shower_x /= Esum
        shower_y /= Esum
        shower_z /= Esum

        buffer[0][0] = Esum
        buffer[1][0] = shower_x
        buffer[2][0] = shower_y
        buffer[3][0] = shower_z
        ntr.Fill()

    ntr.Write()
    nf.Close()
    f.Close()
