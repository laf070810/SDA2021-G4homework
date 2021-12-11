from ROOT import *
import json

Esum_means = [0 for _ in range(20)]
Esum_mean_errors = [0 for _ in range(20)]
Esum_sigmas = [0 for _ in range(20)]
Esum_sigmas_errors = [0 for _ in range(20)]

shower_x_means = [0 for _ in range(20)]
shower_x_mean_errors = [0 for _ in range(20)]
shower_x_sigmas = [0 for _ in range(20)]
shower_x_sigmas_errors = [0 for _ in range(20)]

shower_y_means = [0 for _ in range(20)]
shower_y_mean_errors = [0 for _ in range(20)]
shower_y_sigmas = [0 for _ in range(20)]
shower_y_sigmas_errors = [0 for _ in range(20)]

shower_z_means = [0 for _ in range(20)]
shower_z_mean_errors = [0 for _ in range(20)]
shower_z_sigmas = [0 for _ in range(20)]
shower_z_sigmas_errors = [0 for _ in range(20)]

for i in range(20):
    photon_energy = 0.5 * (i + 1)

    f = TFile(f'data/E_and_position_{photon_energy}GeV.root')
    tr = f.Get('E_and_position')
    nf = TFile(f'data/fit_result_{photon_energy}GeV.root', 'recreate')

    w = RooWorkspace('Esum_fit')
    w.factory(
        f'RooGaussian::gauss(Esum[{0 + 70 * i}, {400 + 70 * i}], mean[{0 + 70 * i}, {400 + 70 * i}], sigma[30, 0.1, 80])')
    Esum = w.var('Esum')
    mean = w.var('mean')
    sigma = w.var('sigma')
    model = w.pdf('gauss')
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

    w = RooWorkspace('shower_z_fit')
    w.factory(
        f'RooGaussian::gauss(shower_z[{-200 + 3 * i}, {50 + 1.5 * i}], mean[{-200 + 3 * i}, {50 + 1.5 * i}], sigma[50, 0, 100])')
    shower_z = w.var('shower_z')
    mean = w.var('mean')
    sigma = w.var('sigma')
    model = w.pdf('gauss')
    data = RooDataSet('ShowerZData', 'ShowerZData', tr, shower_z)
    result = model.fitTo(data, RooFit.Save())

    shower_z_means[i] = mean.getValV()
    shower_z_mean_errors[i] = mean.getError()
    shower_z_sigmas[i] = sigma.getValV()
    shower_z_sigmas_errors[i] = sigma.getError()
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
result['shower_z_means'] = shower_z_means
result['shower_z_mean_errors'] = shower_z_mean_errors
result['shower_z_sigmas'] = shower_z_sigmas
result['shower_z_sigma_errors'] = shower_z_sigmas_errors
with open('data/calibration.json', 'w') as f:
    json.dump(result, f)
