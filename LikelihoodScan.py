#!/usr/bin/env python
# coding: utf-8

# ### Import necessary modules 

# In[ ]:


import ROOT
from ROOT import RooFit
from ROOT import RooRealVar, RooAddPdf, RooArgSet, RooArgList, TCanvas, TH1F, TH1D, TH2F
import ctypes

from IPython import get_ipython
# ### Setup a model and generate a toy data

# In[ ]:


nbins = 100
RegL    = RooRealVar("RegL","RegL",-1.36, 0.83)
ImgL    = RooRealVar("ImgL","ImgL",-0.14, 0.49)

mean1   = RooRealVar("mean1","Re(g_{R}^{V})", -1.3, .8)
sigma1  = RooRealVar("sigma1","sigma1", 0.001, 1.5)
gaus1   = ROOT.RooGaussian("gaus1","gaus1", RegL, mean1, sigma1)

mean2   = RooRealVar("mean2","Im(g_{R}^{V})", -0.1, 0.4)
sigma2  = RooRealVar("sigma2","sigma2", 0.001, 0.3)
gaus2   = ROOT.RooGaussian("gaus2","gaus2", ImgL, mean2, sigma2)

frac    = RooRealVar("frac", "frac", 0,1)

genmodel= RooAddPdf("genmodel", "genmodel", RooArgList(gaus1, gaus2), RooArgList(frac), True)
data    = genmodel.generate(RooArgSet(RegL, ImgL),100000)


# ### Plot projections of variables from data

# In[ ]:


glr_h = ROOT.RooAbsData.createHistogram(data, "glr_h", RegL)
gli_h = ROOT.RooAbsData.createHistogram(data, "gli_h", ImgL)

glr = RegL.frame(RooFit.Title("#italic{Re}(g^{V}_{L})"), ROOT.RooFit.Bins(nbins))
data.plotOn(glr, RooFit.DataError(ROOT.RooAbsData.SumW2))

gli = ImgL.frame(RooFit.Title("#italic{Im}(g^{V}_{L})"), ROOT.RooFit.Bins(nbins))
data.plotOn(gli, RooFit.DataError(ROOT.RooAbsData.SumW2))

c1 = TCanvas()
glr.Draw()
c2 = TCanvas()
gli.Draw()

'''
# ### Create and fit a model

# In[ ]:


mean11   = RooRealVar("mean11","Re(g_{R}^{V})", -1.3, .8)
sigma11  = RooRealVar("sigma11","sigma1", 0.001, 1.5)
gaus11   = ROOT.RooGaussian("gaus11","gaus1", RegL, mean11, sigma11)
mean22   = RooRealVar("mean22","Im(g_{R}^{V})", -0.1, 0.4)
sigma22  = RooRealVar("sigma22","sigma2", 0.001, 0.3)
gaus22   = ROOT.RooGaussian("gaus22","gaus2", ImgL, mean22, sigma22)

fitmodel = RooAddPdf("fitmodel", "fitmodel", RooArgList(gaus11, gaus22), RooArgList(frac), True);
fitRes = fitmodel.fitTo(data, ROOT.RooFit.Save())
fitRes.Print("v")


# ### Plot after fitting

# In[ ]:


glr1 = RegL.frame(RooFit.Title("#italic{Re}(g^{V}_{L})"), RooFit.Bins(nbins))
data.plotOn(glr1, RooFit.DataError(ROOT.RooAbsData.SumW2))
fitmodel.plotOn(glr1)
fitmodel.paramOn(glr1)
fitmodel.plotOn(glr1, RooFit.LineColor(ROOT.kBlue), RooFit.LineWidth(1))


# ### Create a pull

# In[ ]:


pullglr = RegL.frame(RooFit.Title("glr pull"))
hpull1 = glr1.pullHist()

pullglr.addPlotable(hpull1,"P0")
pullglr.SetMinimum(-3)
pullglr.SetMaximum(+3)
pullglr.SetYTitle("pull")
pullglr.SetMarkerStyle(20)
pullglr.SetNdivisions(10)
chisquare_glr1 = glr1.chiSquare()
print("Chi square of glr fit is :", chisquare_glr1)

c3 = ROOT.TCanvas("c", "c",0,0,600,600)
pad1 = ROOT.TPad("pad1","pad1",0,0.36,1,1)
pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.25)
pad1.SetBottomMargin(0.00001)
pad1.SetBorderMode(0)
pad2.SetTopMargin(0.00001)
pad2.SetBottomMargin(0.1)
pad2.SetBorderMode(0)
pad1.Draw()
pad2.Draw()
pad1.cd()
ROOT.gStyle.SetOptTitle(0)
c3.SetFillColor(0)
c3.SetBorderSize(2)
c3.SetLeftMargin(0.1422222)
c3.SetRightMargin(0.04444445)
glr1.SetStats(0)
glr1.Draw()
pad2.cd()
pullglr.SetStats(0)
pullglr.Draw()
c3.cd()


# In[ ]:


# imaginary part of the gl

gli1 = ImgL.frame(RooFit.Title("#italic{Im}(g^{V}_{L})"), RooFit.Bins(nbins))
data.plotOn(gli1, RooFit.DataError(ROOT.RooAbsData.SumW2))
fitmodel.plotOn(gli1)
fitmodel.paramOn(gli1)
fitmodel.plotOn(gli1, RooFit.LineColor(4), RooFit.LineWidth(1))

pullgli = ImgL.frame(RooFit.Title("gli pull"))
hpull2 = gli1.pullHist()

pullgli.addPlotable(hpull2,"P0")
pullgli.SetMinimum(-3)
pullgli.SetMaximum(+3) 
pullgli.SetYTitle("pull")
pullgli.SetMarkerStyle(20)
pullgli.SetNdivisions(10)
chisquare_gli1 = gli1.chiSquare()
print("Chi square of gli fit is :", chisquare_gli1)

cc = ROOT.TCanvas("cc", "cc",0,0,600,600)
pad11 = ROOT.TPad("pad11","pad11",0,0.36,1,1)
pad22 = ROOT.TPad("pad22","pad22",0,0,1,0.25)
pad11.SetBottomMargin(0.00001)
pad11.SetBorderMode(0)
pad22.SetTopMargin(0.00001)
pad22.SetBottomMargin(0.1)
pad22.SetBorderMode(0)
pad11.Draw()
pad22.Draw()
pad11.cd()
ROOT.gStyle.SetOptTitle(0)
cc.SetFillColor(0)
cc.SetBorderSize(2)
cc.SetLeftMargin(0.1422222)
cc.SetRightMargin(0.04444445)
gli1.SetStats(0)
gli1.Draw()
pad22.cd()
pullgli.SetStats(0)
pullgli.Draw()
cc.cd()


# In[ ]:


glr_var = fitRes.floatParsFinal().find("mean11")
gli_var = fitRes.floatParsFinal().find("mean22")

glr_mean = glr_var.getVal()
glr_err = glr_var.getAsymErrorHi()

gli_mean = gli_var.getVal()
gli_err = gli_var.getAsymErrorHi()

glr_err /= 4.
gli_err /= 4.
nll = fitmodel.createNLL(data, RooFit.NumCPU(8))
nll.SetName("nll")
minuit = ROOT.RooMinuit(nll)
minuit.migrad()


# In[ ]:


amin,edm,errdef,nvpar,nparx = ctypes.c_double(0.), ctypes.c_double(0.), ctypes.c_double(0.), ctypes.c_int(0), ctypes.c_int(0)
nllMin = nll.getVal()
tmp = str(nllMin)
tmp += " )";
likl = ROOT.RooFormulaVar("likl", "L", "exp(-nll  + "+tmp, RooArgList(nll))
manL = ROOT.TH2D("manL", "likelihood contour", 40, glr_mean-(20*glr_err), glr_mean+(20*glr_err),40,gli_mean-(20*gli_err), gli_mean+(20*gli_err))
manL.SetXTitle(mean11.GetTitle())
manL.SetYTitle(mean22.GetTitle())
glr_m = glr_mean-(20*glr_err);
for i in range(1, manL.GetNbinsX()):
  gli_m = gli_mean-(20*gli_err);
  for j in range(1, manL.GetNbinsX()):
    #;++j,gli_m+=gli_err)
    x = manL.GetXaxis().GetBinCenter(i)
    mean11.setVal(x)
    y = manL.GetYaxis().GetBinCenter(j)
    mean22.setVal(y)
    fitter = ROOT.TVirtualFitter.Fitter(glr_h)
    fitter.GetStats(amin,edm,errdef,nvpar,nparx)
    manL.SetBinContent(i, j, amin.value - likl.getVal())
    gli_m += gli_err
  glr_m += glr_err

con = ROOT.TCanvas("con", "con", 0, 0, 800,600)
manL.Draw("cont4z")

m = ROOT.RooMinuit(nll)
m.migrad()
m.hesse()
con1 = ROOT.TCanvas("con1", "con1", 0, 0, 800,600)
p1   = m.contour(mean11, mean22, 1, 2, 3) 
p1.Draw("CONT4Z")


# #### Do something to make all output plots visible

# In[ ]:


get_ipython().run_cell_magic('javascript', '', 'IPython.OutputArea.auto_scroll_threshold = 9999;')


# ### Plot all set of plot at once

# In[ ]:


get_ipython().run_line_magic('jsroot', 'off')
ROOT.gROOT.GetListOfCanvases().Draw()


# In[ ]:


'''

