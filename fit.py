import json
import numpy as np
from ROOT import *

Esum_means = [0 for _ in range(6)]
Esum_mean_errors = [0 for _ in range(6)]
Esum_sigmas = [0 for _ in range(6)]
Esum_sigmas_errors = [0 for _ in range(6)]

shower_x_means = [0 for _ in range(6)]
shower_x_mean_errors = [0 for _ in range(6)]
shower_x_sigmas = [0 for _ in range(6)]
shower_x_sigmas_errors = [0 for _ in range(6)]

shower_y_means = [0 for _ in range(6)]
shower_y_mean_errors = [0 for _ in range(6)]
shower_y_sigmas = [0 for _ in range(6)]
shower_y_sigmas_errors = [0 for _ in range(6)]

for i in range(6):
    photon_energy = [1, 2, 5, 10, 50, 100][i]

    f = TFile(f'data/E_and_position_{photon_energy}GeV.root')
    tr = f.Get('E_and_position')
    nf = TFile(f'data/fit_result_{photon_energy}GeV.root', 'recreate')

    w = RooWorkspace('Esum_fit')
    w.factory(
        f'RooCBShape::cbs(Esum[{-100 + 130 * photon_energy}, {100 + 170 * photon_energy}], mean[{-100 + 130 * photon_energy}, {100 + 170 * photon_energy}], sigma[0.1, {100 * np.sqrt(photon_energy)}], alpha[0, 5], n[0, 2])')
    Esum = w.var('Esum')
    mean = w.var('mean')
    sigma = w.var('sigma')
    model = w.pdf('cbs')
    data = RooDataSet('EsumData', 'EsumData', tr, Esum)
    result = model.fitTo(data, RooFit.Save())

    Esum_means[i] = mean.getValV()
    Esum_mean_errors[i] = mean.getError()
    Esum_sigmas[i] = sigma.getValV()
    Esum_sigmas_errors[i] = sigma.getError()
    data.Write()
    w.Write()
    result.Write()

    w = RooWorkspace('shower_x_fit')
    w.factory(
        f'RooGaussian::gauss(shower_x[{-15 + 0.5 * i}, {15 - 0.5 * i}], mean[0, {-15 + 0.5 * i}, {15 - 0.5 * i}], sigma[0, {15 - 0.5 * i}])')
    shower_x = w.var('shower_x')
    mean = w.var('mean')
    sigma = w.var('sigma')
    model = w.pdf('gauss')
    data = RooDataSet('ShowerXData', 'ShowerXData', tr, shower_x)
    result = model.fitTo(data, RooFit.Save())

    shower_x_means[i] = mean.getValV()
    shower_x_mean_errors[i] = mean.getError()
    shower_x_sigmas[i] = sigma.getValV()
    shower_x_sigmas_errors[i] = sigma.getError()
    data.Write()
    w.Write()
    result.Write()

    w = RooWorkspace('shower_y_fit')
    w.factory(
        f'RooGaussian::gauss(shower_y[{-15 + 0.5 * i}, {15 - 0.5 * i}], mean[0, {-15 + 0.5 * i}, {15 - 0.5 * i}], sigma[0, {15 - 0.5 * i}])')
    shower_y = w.var('shower_y')
    mean = w.var('mean')
    sigma = w.var('sigma')
    model = w.pdf('gauss')
    data = RooDataSet('ShowerYData', 'ShowerYData', tr, shower_y)
    result = model.fitTo(data, RooFit.Save())

    shower_y_means[i] = mean.getValV()
    shower_y_mean_errors[i] = mean.getError()
    shower_y_sigmas[i] = sigma.getValV()
    shower_y_sigmas_errors[i] = sigma.getError()
    data.Write()
    w.Write()
    result.Write()

    nf.Close()
    f.Close()

result = {}
result['Esum_means'] = Esum_means
result['Esum_mean_errors'] = Esum_mean_errors
result['Esum_sigmas'] = Esum_sigmas
result['Esum_sigma_errors'] = Esum_sigmas_errors
result['shower_x_means'] = shower_x_means
result['shower_x_mean_errors'] = shower_x_mean_errors
result['shower_x_sigmas'] = shower_x_sigmas
result['shower_x_sigma_errors'] = shower_x_sigmas_errors
result['shower_y_means'] = shower_y_means
result['shower_y_mean_errors'] = shower_y_mean_errors
result['shower_y_sigmas'] = shower_y_sigmas
result['shower_y_sigma_errors'] = shower_y_sigmas_errors
with open('data/calibration_data.json', 'w') as f:
    json.dump(result, f)
