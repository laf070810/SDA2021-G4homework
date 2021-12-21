import json
import numpy as np
from ROOT import *

f = TFile(f'data/E_and_position_pi0.root')
tr = f.Get('E_and_position')
nf = TFile(f'data/fit_result_pi0.root', 'recreate')

w = RooWorkspace('M_pi0_fit')
w.factory(
    f'RooGaussian::gauss(M_pi0[80, 200], mean[80, 200], sigma[1, 100])')
M_pi0 = w.var('M_pi0')
mean = w.var('mean')
sigma = w.var('sigma')
model = w.pdf('gauss')
data = RooDataSet('Mpi0Data', 'Mpi0Data', tr, M_pi0)
result = model.fitTo(data, RooFit.Save())

data.Write()
w.Write()
result.Write()

nf.Close()
f.Close()
