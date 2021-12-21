from ROOT import *
import numpy as np


scintillator_layer_num = 67
cell_num_on_x_per_layer = 10
cell_num_on_y_per_layer = 10
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
world_size_z = 1.2 * (calorimeter_size_z + 2 * 5000)
first_position = -calorimeter_size_z / 2


for i in range(6):
    photon_energy = [1, 2, 5, 10, 50, 100][i]

    f = TFile(f'data/B4_{i}.root')
    tr = f.Get('B4')

    nf = TFile(f'data/E_and_position_{photon_energy}GeV.root', 'recreate')
    ntr = TTree('E_and_position', '')
    buffer = [np.zeros(1) for _ in range(3)]
    ntr.Branch('Esum', buffer[0], 'Esum/D')
    ntr.Branch('shower_x', buffer[1], 'shower_x/D')
    ntr.Branch('shower_y', buffer[2], 'shower_y/D')

    radius = 100

    num = 0
    for event in tr:
        if num % 100 == 0:
            print(f"{num}/{tr.GetEntries()}")
        num += 1

        index = 0
        Esum_xy = np.zeros((cell_num_on_x_per_layer, cell_num_on_y_per_layer))
        for iz in range(scintillator_layer_num):
            for iy in range(cell_num_on_y_per_layer):
                for ix in range(cell_num_on_x_per_layer):
                    Esum_xy[ix, iy] += event.Eabs[index]
                    index += 1

        seed_x_index, seed_y_index = np.unravel_index(np.argmax(Esum_xy), Esum_xy.shape)
        seed_x = -calorimeter_size_x / 2 + seed_x_index * cell_width_on_x + cell_width_on_x / 2
        seed_y = -calorimeter_size_y / 2 + seed_y_index * cell_width_on_y + cell_width_on_y / 2

        Esum = 0
        shower_x = 0
        shower_y = 0
        for iy in range(cell_num_on_y_per_layer):
            for ix in range(cell_num_on_x_per_layer):
                x_position = -calorimeter_size_x / 2 + ix * cell_width_on_x + cell_width_on_x / 2
                y_position = -calorimeter_size_y / 2 + iy * cell_width_on_y + cell_width_on_y / 2

                if ((x_position - seed_x) ** 2 + (y_position - seed_y) ** 2 > radius ** 2):
                    continue

                Esum += Esum_xy[ix, iy]
                shower_x += x_position * Esum_xy[ix, iy]
                shower_y += y_position * Esum_xy[ix, iy]
        shower_x /= Esum
        shower_y /= Esum

        buffer[0][0] = Esum
        buffer[1][0] = shower_x
        buffer[2][0] = shower_y
        ntr.Fill()

    ntr.Write()
    nf.Close()
    f.Close()
