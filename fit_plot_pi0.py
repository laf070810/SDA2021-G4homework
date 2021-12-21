from ROOT import *
gROOT.ProcessLine(".x ~/lhcbStyle.C")
gROOT.SetBatch(True)


f = TFile(f'data/fit_result_pi0.root')

w = f.Get('M_pi0_fit')
M_pi0 = w.var('M_pi0')
model = w.pdf('gauss')
data = f.Get('Mpi0Data')

frame = M_pi0.frame()
data.plotOn(frame)
model.plotOn(frame)
c = TCanvas()
frame.GetXaxis().SetTitle("M_{pi^0} (MeV)")
# frame.GetYaxis().SetTitle("Events/6 MeV")
frame.Draw()
c.SaveAs(f"hists/fit_M_pi0.pdf")

f.Close()
