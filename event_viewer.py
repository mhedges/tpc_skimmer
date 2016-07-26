from ROOT import TFile, TTree, gROOT, TCanvas, TH1F, TH2F

ifile = TFile('test.root')
t = gROOT.FindObject('tr')

h = TH2F('h', 'Track', 80, 0.5, 80.5, 336, 0.5, 336.5)
h2=TH1F('h2','bcids',256,0.5,256.5)

for event in t:
    np = event.npoints
    if max(event.distances) > 3E6:
        for i in xrange(np):
            h.Fill(event.col[i],event.row[i],event.tot[i])
            h2.Fill(event.bcid[i])
        c = TCanvas()
        c.Divide(2,1)
        c.cd(1)
        h.Draw("COLZ")
        c.cd(2)
        h2.Draw()
        print 'Number of distances = ', len(event.distances), len(event.row)
        print 'Number of points = ', np
        raw_input('Well?')
    h.Reset()
    h2.Reset()
