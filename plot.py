from ROOT import *
gROOT.ProcessLine(".x ~/lhcbStyle.C")


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


for i in range(20):
    photon_energy = 0.5 * (i + 1)

    f = TFile(f'data/E_and_position_{photon_energy}GeV.root')
    tr = f.Get('E_and_position')

    nBins = 50
    hE = TH1F(f"hE_{photon_energy}GeV", "E_{sum} (MeV)", nBins, 0 + 70 * i, 400 + 70 * i)
    hShowerX = TH1F(f"hShowerX_{photon_energy}GeV", "x (mm)", nBins, -15 + 0.5 * i, 15 - 0.5 * i)
    hShowerY = TH1F(f"hShowerY_{photon_energy}GeV", "y (mm)", nBins, -15 + 0.5 * i, 15 - 0.5 * i)
    hShowerZ = TH1F(f"hShowerZ_{photon_energy}GeV", "z (mm)", nBins, -200 + 3 * i, 50 + 1.5 * i)

    num = 0
    for event in tr:
        if num % 1000 == 0:
            print(f"{num}/{tr.GetEntries()}")
        num += 1

        hE.Fill(event.Esum)
        hShowerX.Fill(event.shower_x)
        hShowerY.Fill(event.shower_y)
        hShowerZ.Fill(event.shower_z)

    draw_1Dhists([hE, hShowerX, hShowerY, hShowerZ])
