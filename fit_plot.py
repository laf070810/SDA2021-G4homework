import json
import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from ROOT import *

import matplotlib
matplotlib.use('Agg')

gROOT.ProcessLine(".x ~/lhcbStyle.C")
gROOT.SetBatch(True)


with open('data/calibration_data.json') as f:
    calibration_result = json.load(f)

Esum_means = np.array(calibration_result['Esum_means'])
Esum_mean_errors = np.array(calibration_result['Esum_mean_errors'])
Esum_sigmas = np.array(calibration_result['Esum_sigmas'])
Esum_sigmas_errors = np.array(calibration_result['Esum_sigma_errors'])

shower_x_means = np.array(calibration_result['shower_x_means'])
shower_x_mean_errors = np.array(calibration_result['shower_x_mean_errors'])
shower_x_sigmas = np.array(calibration_result['shower_x_sigmas'])
shower_x_sigmas_errors = np.array(calibration_result['shower_x_sigma_errors'])

shower_y_means = np.array(calibration_result['shower_y_means'])
shower_y_mean_errors = np.array(calibration_result['shower_y_mean_errors'])
shower_y_sigmas = np.array(calibration_result['shower_y_sigmas'])
shower_y_sigmas_errors = np.array(calibration_result['shower_y_sigma_errors'])

photon_energy = np.array([1, 2, 5, 10, 50, 100])

plt.style.use('classic')
plt.figure(figsize=(20, 15))
plt.subplot(221)
plt.errorbar(photon_energy, Esum_means, yerr=Esum_sigmas, fmt='ko')
plt.xlabel(r'$E_{ph} \  (\mathrm{GeV})$', fontweight='bold', fontsize=24, loc='right')
plt.ylabel(r'$\mu(E_{sum}) \  (\mathrm{MeV})$', fontweight='bold', fontsize=24, loc='top')
plt.title(r'$\mu(E_{sum})$', y=1.05, fontweight='bold', fontsize=24)
plt.subplot(222)
plt.plot(photon_energy, Esum_sigmas / Esum_means, 'ko')
plt.xlabel(r'$E_{ph} \  (\mathrm{GeV})$', fontweight='bold', fontsize=24, loc='right')
plt.ylabel(r'$\frac{\sigma(E_{sum})}{\mu(E_{sum})}$', fontweight='bold', fontsize=24, loc='top')
plt.title(r'$\frac{\sigma(E_{sum})}{\mu(E_{sum})}$', y=1.05, fontweight='bold', fontsize=24)
plt.subplot(223)
plt.plot(photon_energy, shower_x_sigmas, 'ko')
plt.xlabel(r'$E_{ph} \  (\mathrm{GeV})$', fontweight='bold', fontsize=24, loc='right')
plt.ylabel(r'$\sigma(x) \  (\mathrm{mm})$', fontweight='bold', fontsize=24, loc='top')
plt.title(r'$\sigma(x)$', y=1.05, fontweight='bold', fontsize=24)
plt.subplot(224)
plt.plot(photon_energy, shower_y_sigmas, 'ko')
plt.xlabel(r'$E_{ph} \  (\mathrm{GeV})$', fontweight='bold', fontsize=24, loc='right')
plt.ylabel(r'$\sigma(y) \  (\mathrm{mm})$', fontweight='bold', fontsize=24, loc='top')
plt.title(r'$\sigma(y)$', y=1.05, fontweight='bold', fontsize=24)
plt.subplots_adjust(hspace=0.4, wspace=0.4)
plt.savefig('hists/calibration.png')

A = np.vstack([Esum_means, np.ones(len(Esum_means))]).T
beta1, beta0 = np.linalg.lstsq(A, photon_energy * 1000)[0]
print('calibration equation: E_gamma = beta1 * E_sum + beta0')
print(f'beta1: {beta1} (1)')
print(f'beta0: {beta0} (MeV)')
print()
calibration_result = {}
calibration_result['beta1'] = beta1
calibration_result['beta0'] = beta0


def func(x, u, y):
    alpha = x[0]
    beta = x[1]
    return np.sqrt((alpha / np.sqrt(u)) ** 2 + beta ** 2) - y


x0 = [1, 0]
res = least_squares(func, x0, args=(photon_energy, Esum_sigmas / Esum_means))
print(res)
print()
calibration_result['alpha'] = res.x[0]
calibration_result['beta'] = res.x[1]
with open('data/calibration_result.json', 'w') as f:
    json.dump(calibration_result, f)

plt.style.use('classic')
plt.figure(figsize=(20, 15))
plt.plot(photon_energy, Esum_sigmas / Esum_means, 'ko')
x = np.linspace(1, 100, 100)
plt.plot(x, np.sqrt((res.x[0] / np.sqrt(x)) ** 2 + res.x[1] ** 2))
plt.xlabel(r'$E_{ph} \  (\mathrm{GeV})$', fontweight='bold', fontsize=24, loc='right')
plt.ylabel(r'$\frac{\sigma(E_{sum})}{\mu(E_{sum})}$', fontweight='bold', fontsize=24, loc='top')
plt.title(r'$\frac{\sigma(E_{sum})}{\mu(E_{sum})}$', y=1.05, fontweight='bold', fontsize=24)
plt.savefig('hists/fit_calibration.png')

A = np.vstack([photon_energy, np.ones(len(photon_energy))]).T
beta1, beta0 = np.linalg.lstsq(A, Esum_sigmas ** 2 / Esum_means)[0]
print(f'Fano factor: {beta1 * 1000} (MeV)')
print(f'beta0: {beta0} (MeV)')

plt.style.use('classic')
plt.figure(figsize=(30, 10))
plt.subplot(121)
plt.plot(photon_energy, Esum_sigmas ** 2, 'ko')
plt.xlabel(r'$E_{ph} \  (\mathrm{GeV})$', fontweight='bold', fontsize=24, loc='right')
plt.ylabel(r'$Var(E_{sum}) \  (\mathrm{MeV}^2)$', fontweight='bold', fontsize=24, loc='top')
plt.title(r'$Var(E_{sum})$', y=1.05, fontweight='bold', fontsize=24)
plt.subplot(122)
plt.plot(photon_energy, Esum_sigmas ** 2 / Esum_means, 'ko', label='Original data')
plt.xlabel(r'$E_{ph} \  (\mathrm{GeV})$', fontweight='bold', fontsize=24, loc='right')
plt.ylabel(r'$\frac{Var(E_{sum})}{\mu(E_{sum})} \  (\mathrm{MeV})$', fontweight='bold', fontsize=24, loc='top')
plt.title(r'$\frac{Var(E_{sum})}{\mu(E_{sum})}$', y=1.05, fontweight='bold', fontsize=24)
plt.plot(photon_energy, beta1 * photon_energy + beta0, label='Fitted line')
plt.subplots_adjust(hspace=0.4, wspace=0.4)
plt.savefig('hists/variance.png')

for i in range(6):
    photon_energy = [1, 2, 5, 10, 50, 100][i]

    f = TFile(f'data/fit_result_{photon_energy}GeV.root')

    w = f.Get('Esum_fit')
    Esum = w.var('Esum')
    model = w.pdf('gauss')
    data = f.Get('EsumData')

    frame = Esum.frame()
    data.plotOn(frame)
    model.plotOn(frame)
    c = TCanvas()
    frame.GetXaxis().SetTitle("E_{sum} (MeV)")
    # frame.GetYaxis().SetTitle("Events/6 MeV")
    frame.Draw()
    c.SaveAs(f"hists/fit_E_{photon_energy}GeV.pdf")

    w = f.Get('shower_x_fit')
    shower_x = w.var('shower_x')
    model = w.pdf('gauss')
    data = f.Get('ShowerXData')

    frame = shower_x.frame()
    data.plotOn(frame)
    model.plotOn(frame)
    c = TCanvas()
    frame.GetXaxis().SetTitle("x (mm)")
    # frame.GetYaxis().SetTitle("Events/6 MeV")
    frame.Draw()
    c.SaveAs(f"hists/fit_ShowerX_{photon_energy}GeV.pdf")

    w = f.Get('shower_y_fit')
    shower_y = w.var('shower_y')
    model = w.pdf('gauss')
    data = f.Get('ShowerYData')

    frame = shower_y.frame()
    data.plotOn(frame)
    model.plotOn(frame)
    c = TCanvas()
    frame.GetXaxis().SetTitle("y (mm)")
    # frame.GetYaxis().SetTitle("Events/6 MeV")
    frame.Draw()
    c.SaveAs(f"hists/fit_ShowerY_{photon_energy}GeV.pdf")

    f.Close()
