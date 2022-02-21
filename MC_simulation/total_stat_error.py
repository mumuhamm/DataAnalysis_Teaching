from ROOT import TCanvas, TH1F, TSlider, THStack, TPad, TFile, TPaveLabel, TPaveText, TLegend
from ROOT import gROOT, gBenchmark, gRandom
import sys, os, ROOT


ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetPadTopMargin(0.05)
ROOT.gStyle.SetPadRightMargin(0.05)
ROOT.gStyle.SetPadBottomMargin(0.16)
ROOT.gStyle.SetPadLeftMargin(0.20)
ROOT.gStyle.SetPadBorderMode(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

# Create a new canvas, and customize it.
c1 = ROOT.TCanvas( 'c1', 'c1', 200, 10, 600, 600 )
c1.SetGrid();

hs = ROOT.THStack("hs","Stacked 1D histograms")
h1st = ROOT.TH1F("h1st","test hstack",40,-4,4)
h1st.FillRandom("gaus",20000)
h1st.SetFillColor(2)
h1st.SetMarkerStyle(21)
h1st.SetMarkerColor(2) #Red
hs.Add(h1st)
h2st = ROOT.TH1F("h2st","test hstack",40,-4,4)
h2st.FillRandom("gaus",15000)
h2st.SetFillColor(4)
h2st.SetMarkerStyle(21)
h2st.SetMarkerColor(4)#Blue
hs.Add(h2st)
h3st = ROOT.TH1F("h3st","test hstack",40,-4,4)
h3st.FillRandom("gaus",10000)
h3st.SetFillColor(8)
h3st.SetMarkerStyle(21)
h3st.SetMarkerColor(8)#Green
hs.Add(h3st)
hs.SetTitle('Stack histogram')
htot = ROOT.TH1F("htot","error hstack",40,-4,4)
htot.Add(h1st)
htot.Add(h2st)
htot.Add(h3st)
htot.SetMarkerStyle(8)
htot.SetMarkerColor(1)
hs.Draw()
htot.Draw('ELP SAME')
hs.GetHistogram().GetXaxis().SetTitle('G(x)') #ths.GetHistogram()
hs.GetHistogram().GetYaxis().SetTitle('Stacked frequency')
leg = ROOT.TLegend(0.68,0.72,0.98,0.92)
leg.AddEntry(h1st, 'red')
leg.AddEntry(h2st, 'blue')
leg.AddEntry(h3st, 'green')
leg.AddEntry(htot, 'total+Error')
leg.Draw()
c1.Modified()
c1.Update()
c1.SaveAs("staterror.png")

'''

#!/usr/bin/python
import sys, os, ROOT

file = ROOT.TFile( sys.argv[1] )
tree = file.Get("fTree")
minRun = tree.GetMinimum("fRunNumber") - 10
maxRun = tree.GetMaximum("fRunNumber") + 10

# create a canvas and divide it
c1 = ROOT.TCanvas("paddy","paddy",700,700)
c1.Divide(1,2,0,0)

# Histogram for the reset rate in Hz
h1 = ROOT.TH2F("h1","h1",10000,minRun,maxRun,10000,7200,7450)
h1.SetMarkerStyle(24)
h1.SetMarkerSize(0.4)
h1.GetYaxis().SetTitle("Baseline Value (ADC)")

# Histogram for the Nitrogen Level
h2 = ROOT.TH2F("h2","h2",10000,minRun,maxRun,10000,72,98)
h2.SetMarkerStyle(24)
h2.SetMarkerSize(0.4)
h2.GetXaxis().SetTitle("Run Number")
h2.GetYaxis().SetTitle("LN_{2} Level (%)")

#top canvas
c1.cd(1)
ROOT.gPad.SetBottomMargin(0.001)
ROOT.gPad.SetTopMargin(0.01)
ROOT.gPad.SetRightMargin(0.01)
tree.Draw("baselineValueMean:fRunNumber>>h1","","goff")
h1.Draw()

# bottom canvas
c1.cd(2)
ROOT.gPad.SetTopMargin(0.001)
ROOT.gPad.SetRightMargin(0.01)
tree.Draw("fDetLnLevelAtBeginRun:fRunNumber>>h2","","goff")
h2.Draw()

# update the canvas
c1.Update()
raw_input('press enter')

=============

from PlotterToolsDataMC import *
from ROOT import *
from ROOT import TPad
#gStyle.SetOptStat(0)
#gROOT.SetBatch(kTRUE)
#first loop over all variables
for v in Vars:
   hists = []
   hists2 = []
   leg = TLegend(0.6, 0.8, 0.89, 0.89)
   
   leg.SetBorderSize(0)
   Max = -0.
   Max2 = -0.
   for fi,f in enumerate(MC):  ## get all MC plots, because they have to stacked!
      ch = TChain('H4GSel')
      ch.Add(f[0])
      hname = v[1]+'_'+str(fi)
      h = TH1F(hname, v[2], v[3], v[4], v[5])
      ch.Draw(v[0]+'>>'+hname,TCut(genCut)) ## add cut based on what you want to plot, blind or unblind or anything else
      #total1 = h.Integral
      #print total1
      h.Scale(float(f[4]),"nosw2")
      #h.Scale(1/float(h.Integral()))
      h.SetLineColor(f[2])
      h.SetLineWidth(2)
      h.SetFillColor(f[3])
      h.Sumw2()
      mc_copy = h.Clone("copy")
      hists.append([h,ch,f[1]])
      if h.GetMaximum() > Max:
         Max = h.GetMaximum()
   c0 = TCanvas('a','a',800,2000)   ##now starts the drawing part, start by stacking all the MC up + add a ratio plot
   SetOwnership(c0,False)
 s = THStack("s","")
   for fi,hh in enumerate(hists):
      #leg.AddEntry(hh[0], hh[2], 'lf')
      s.Add(hh[0])
      #s.Add(h3)
      hh[0].SetMaximum(Max*1.5)
      hh[0].SetMinimum(0.0001)
      if fi == 0:
         hh[0].Draw('')
      if fi > 0:
         hh[0].Draw('same')
      if fi == 0:
         leg.AddEntry(hh[0],hh[2],'lf')
      if fi == 2:
         leg.AddEntry(hh[0],hh[2],'lf')
      if fi == 5:
         leg.AddEntry(hh[0],hh[2],'lf')
      #h_err = TH1F(h_err, h_err, 35, 100,180)
      #h_err.Add(hh[0])
      #h_err.Add(hh[1])
      #h_err.Add(hh[2])
   #s.Add(h3)
   s.Draw("hist E2")
   s.GetXaxis().SetTitle(v[6])
   s.GetYaxis().SetTitle('Normalized Yields')
   s.GetYaxis().SetTitleOffset(1.6)
   s.SetMaximum(100000)
   
   ===============


gBenchmark.Start( 'hsum' )

# Create some histograms.
total  = TH1F( 'total', 'This is the total distribution', 100, -4, 4 )
main   = TH1F( 'main', 'Main contributor', 100, -4, 4 )
s1     = TH1F( 's1', 'This is the first signal', 100, -4, 4 )
s2     = TH1F( 's2', 'This is the second signal', 100, -4, 4 )
total.Sumw2()   # this makes sure that the sum of squares of weights will be stored

# Set canvas/frame attributes.
total.SetMarkerStyle( 21 )
total.SetMarkerSize( 0.7 )
main.SetFillColor( 16 )
s1.SetFillColor( 42 )
s2.SetFillColor( 46 )

# Initialize random number generator.
gRandom.SetSeed()
gauss, landau = gRandom.Gaus, gRandom.Landau

# for speed, bind and cache the Fill member functions
histos = [ 'total', 'main', 's1', 's2' ]
for name in histos:
   exec('%sFill = %s.Fill' % (name,name))

# Fill histograms randomly
kUPDATE = 500
for i in range( 10000 ):
 # Generate random values.
   xmain = gauss( -1, 1.5 )
   xs1   = gauss( -0.5, 0.5 )
   xs2   = landau( 1, 0.15 )
   mainFill( xmain )

 # Fill histograms.
   s1Fill( xs1, 0.3 )
   s2Fill( xs2, 0.2 )
   totalFill( xmain )
   totalFill( xs1, 0.3 )
   totalFill( xs2, 0.2 )

 # Update display every kUPDATE events.
   if i and (i%kUPDATE) == 0 :
      if i == kUPDATE :
         total.Draw( 'e1p' )
         main.Draw( 'same' )
         s1.Draw( 'same' )
         s2.Draw( 'same' )
         c1.Update()
         slider = TSlider( 'slider', 'test', 4.2, 0, 4.6, total.GetMaximum(), 38 )
         slider.SetFillColor( 46 )

      if slider:
         slider.SetRange( 0, float(i) / 10000. )

      c1.Modified()
      c1.Update()

# Destroy member functions cache.
for name in histos:
   exec('del %sFill' % name)
del histos

# Done, finalized and trigger an update.
slider.SetRange( 0, 1 )
total.Draw( 'sameaxis' ) # to redraw axis hidden by the fill area

gBenchmark.Show( 'hsum' )
'''
