import json
import numpy as np
import numpy.ma as ma
from ROOT import *


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


with open('data/calibration_result.json') as f:
    calibration_result = json.load(f)

f = TFile(f'data/B4_pi0.root')
tr = f.Get('B4')

nf = TFile(f'data/E_and_position_pi0.root', 'recreate')
ntr = TTree('E_and_position', '')
buffer = [np.zeros(1) for _ in range(7)]
ntr.Branch('Esum1', buffer[0], 'Esum1/D')
ntr.Branch('Esum2', buffer[1], 'Esum2/D')
ntr.Branch('incidence1_x', buffer[2], 'incidence1_x/D')
ntr.Branch('incidence1_y', buffer[3], 'incidence1_y/D')
ntr.Branch('incidence2_x', buffer[4], 'incidence2_x/D')
ntr.Branch('incidence2_y', buffer[5], 'incidence2_y/D')
ntr.Branch('M_pi0', buffer[6], 'M_pi0/D')

radius = 100

num = 0
for event in tr:
    if num % 100 == 0:
        print(f"{num}/{tr.GetEntries()}")
    num += 1

    index = 0
    Esum_xy = ma.zeros((cell_num_on_x_per_layer, cell_num_on_y_per_layer))
    E_xy0 = ma.zeros((cell_num_on_x_per_layer, cell_num_on_y_per_layer))
    for iz in range(scintillator_layer_num):
        for iy in range(cell_num_on_y_per_layer):
            for ix in range(cell_num_on_x_per_layer):
                Esum_xy[ix, iy] += event.Eabs[index]
                if iz < 10:
                    E_xy0[ix, iy] += event.Eabs[index]
                index += 1

    seed_x_index, seed_y_index = np.unravel_index(ma.argmax(Esum_xy), Esum_xy.shape)
    seed_x = -calorimeter_size_x / 2 + seed_x_index * cell_width_on_x + cell_width_on_x / 2
    seed_y = -calorimeter_size_y / 2 + seed_y_index * cell_width_on_y + cell_width_on_y / 2

    Esum1 = 0
    for iy in range(cell_num_on_y_per_layer):
        for ix in range(cell_num_on_x_per_layer):
            x_position = -calorimeter_size_x / 2 + ix * cell_width_on_x + cell_width_on_x / 2
            y_position = -calorimeter_size_y / 2 + iy * cell_width_on_y + cell_width_on_y / 2

            if ((x_position - seed_x) ** 2 + (y_position - seed_y) ** 2 > radius ** 2):
                continue

            Esum1 += Esum_xy[ix, iy]
            Esum_xy[ix, iy] = ma.masked
    Esum1 = calibration_result['beta1'] * Esum1 + calibration_result['beta0']

    seed_x_index, seed_y_index = np.unravel_index(ma.argmax(Esum_xy), Esum_xy.shape)
    seed_x = -calorimeter_size_x / 2 + seed_x_index * cell_width_on_x + cell_width_on_x / 2
    seed_y = -calorimeter_size_y / 2 + seed_y_index * cell_width_on_y + cell_width_on_y / 2

    Esum2 = 0
    for iy in range(cell_num_on_y_per_layer):
        for ix in range(cell_num_on_x_per_layer):
            x_position = -calorimeter_size_x / 2 + ix * cell_width_on_x + cell_width_on_x / 2
            y_position = -calorimeter_size_y / 2 + iy * cell_width_on_y + cell_width_on_y / 2

            if Esum_xy.mask[ix, iy]:
                continue
            if ((x_position - seed_x) ** 2 + (y_position - seed_y) ** 2 > radius ** 2):
                continue

            Esum2 += Esum_xy[ix, iy]
            Esum_xy[ix, iy] = ma.masked
    Esum2 = calibration_result['beta1'] * Esum2 + calibration_result['beta0']

    incidence1_x_index, incidence1_y_index = np.unravel_index(ma.argmax(E_xy0), E_xy0.shape)
    incidence1_x = -calorimeter_size_x / 2 + incidence1_x_index * cell_width_on_x + cell_width_on_x / 2
    incidence1_y = -calorimeter_size_y / 2 + incidence1_y_index * cell_width_on_y + cell_width_on_y / 2

    radius2 = 60
    for iy in range(cell_num_on_y_per_layer):
        for ix in range(cell_num_on_x_per_layer):
            x_position = -calorimeter_size_x / 2 + ix * cell_width_on_x + cell_width_on_x / 2
            y_position = -calorimeter_size_y / 2 + iy * cell_width_on_y + cell_width_on_y / 2
            if ((x_position - incidence1_x) ** 2 + (y_position - incidence1_y) ** 2 > radius2 ** 2):
                continue
            E_xy0[ix, iy] = ma.masked

    incidence2_x_index, incidence2_y_index = np.unravel_index(ma.argmax(E_xy0), E_xy0.shape)
    incidence2_x = -calorimeter_size_x / 2 + incidence2_x_index * cell_width_on_x + cell_width_on_x / 2
    incidence2_y = -calorimeter_size_y / 2 + incidence2_y_index * cell_width_on_y + cell_width_on_y / 2

    a1 = np.sqrt(incidence1_x ** 2 + incidence1_y ** 2 + 5000 ** 2)
    a2 = np.sqrt(incidence2_x ** 2 + incidence2_y ** 2 + 5000 ** 2)
    b = np.sqrt((incidence2_x - incidence1_x) ** 2 + (incidence2_y - incidence1_y) ** 2)
    costheta = (a1 ** 2 + a2 ** 2 - b ** 2) / (2 * a1 * a2)

    buffer[0][0] = Esum1
    buffer[1][0] = Esum2
    buffer[2][0] = incidence1_x
    buffer[3][0] = incidence1_y
    buffer[4][0] = incidence2_x
    buffer[5][0] = incidence2_y
    buffer[6][0] = np.sqrt(2 * Esum1 * Esum2 * (1 - costheta))
    ntr.Fill()

ntr.Write()
nf.Close()
f.Close()
