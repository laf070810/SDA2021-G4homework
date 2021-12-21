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


for i in range(6):
    photon_energy = [1, 2, 5, 10, 50, 100][i]

    f = TFile(f'data/E_and_position_{photon_energy}GeV.root')
    tr = f.Get('E_and_position')

    nBins = 50
    hE = TH1F(f"hE_{photon_energy}GeV", "E_{sum} (MeV)", nBins, -100 + 130 * photon_energy, 100 + 170 * photon_energy)
    hShowerX = TH1F(f"hShowerX_{photon_energy}GeV", "x (mm)", nBins, -100, 100)
    hShowerY = TH1F(f"hShowerY_{photon_energy}GeV", "y (mm)", nBins, -100, 100)

    num = 0
    for event in tr:
        if num % 100 == 0:
            print(f"{num}/{tr.GetEntries()}")
        num += 1

        hE.Fill(event.Esum)
        hShowerX.Fill(event.shower_x)
        hShowerY.Fill(event.shower_y)

    draw_1Dhists([hE, hShowerX, hShowerY])
