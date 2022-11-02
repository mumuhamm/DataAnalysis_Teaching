import ROOT as r
import os, sys
from ROOT import *

input = TFile.Open( 'input.root', 'READ' )
dummy_data = input.FindObjectAny('dummy_data')
mc_1 = input.FindObjectAny('mc_1')
mc_1.SetFillColor(2)
mc_1 = input.FindObjectAny('mc_1')
mc_2.SetFillColor(4)

##### to illustrate initial state of affairs
mc_stack = THStack()
mc_stack.Add( mc_1 )
mc_stack.Add( mc_2 )

##### mc input to TFractionFitter
mc = TObjArray(2)
mc.Add( mc_1 )
mc.Add( mc_2 )

##### fit
fit = TFractionFitter(dummy_data,mc)
status = fit.Fit()

##### get the new normalisations
value_0 = r.Double(0)
error_0 = r.Double(0)
fit.GetResult(0, value_0, error_0)
# ------------------------------ #
value_1 = r.Double(0)
error_1 = r.Double(0)
fit.GetResult(1, value_1, error_1)
##### I assume data = value_0*mc_1 + value_1*mc_2
corr_mc_1 = mc_1.Clone( 'corr_mc_1' )
corr_mc_2 = mc_2.Clone( 'corr_mc_2' )
corr_mc_1.Scale(value_0)
corr_mc_2.Scale(value_1)

#### get new stack
corr_mc_stack = THStack()
corr_mc_stack.Add( corr_mc_1 )
corr_mc_stack.Add( corr_mc_2 )


############## PLOTTING ####################
##### how things looked at the begining
mc_stack.SetMaximum(60)
mc_stack.Draw('hist')
dummy_data.Draw('Ep same')
gPad.SaveAs('step1.pdf')

##### how things look after the fitting
result = fit.GetPlot()
result.SetFillColor(3)
result.SetMaximum(60)
result.Draw("")
dummy_data.Draw("Ep same")
gPad.SaveAs('step2.pdf')

##### cross check I've extracted the right information
# ---- this should be equal to step2.pdf
corr_mc_stack.SetMaximum(60)
corr_mc_stack.Draw('hist')
dummy_data.Draw('Ep same')
gPad.SaveAs('step3.pdf')
