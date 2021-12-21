from ROOT import *
gROOT.ProcessLine(".x ~/lhcbStyle.C")
gROOT.SetBatch(True)


def draw_1Dhists(hists):
    for hist in hists:
        canvas = TCanvas()
        hist.SetStats(False)
        hist.SetMarkerStyle(24)
        hist.SetMarkerColor(kBlue)
        hist.SetLineColor(kBlue)
        hist.SetMinimum(0.)
        hist.SetXTitle(hist.GetTitle())
        hist.Draw("e")
        canvas.SaveAs('hists/' + hist.GetName()[1:] + '.pdf')


f = TFile(f'data/E_and_position_pi0.root')
tr = f.Get('E_and_position')

nBins = 100
hM_pi0 = TH1F("hM_pi0", "M_{pi^{0}} (MeV)", nBins, 20, 200)

num = 0
for event in tr:
    if num % 100 == 0:
        print(f"{num}/{tr.GetEntries()}")
    num += 1

    hM_pi0.Fill(event.M_pi0)

draw_1Dhists([hM_pi0])
