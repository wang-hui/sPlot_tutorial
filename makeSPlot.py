#!/usr/bin/env python3

import ROOT
from array import array

ROOT.gROOT.SetBatch(True)

file_path = "/storage4/store/user/zhipeng/Data_For_ReWeight_v2/2022E_No1.root"
tree_name = "JpsiTree"

data_file = ROOT.TFile.Open(file_path)
input_tree = data_file.Get(tree_name)

mMin, mMax = 2.5, 3.5
Jmass = ROOT.RooRealVar("mass", "M(J/psi)", mMin, mMax)

ctMin, ctMax = -0.03, 0.12
Jctau = ROOT.RooRealVar("ctau", "c#tau(J/#psi) [cm]", ctMin, ctMax)

# Variable-bin edges need to be passed as C-compatible double arrays.
yBins = array("d", [0.0, 0.3, 0.6, 0.9, 1.2, 1.8, 2.4])
ptBins = array("d", [ 
        10.0, 10.5, 11.0, 11.5, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0,
        18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0,
        28.0, 29.0, 30.0, 32.0, 34.0, 36.0, 38.0, 42.0, 46.0, 50.0,
        60.0, 75.0, 95.0, 120.0, 200.0
        ])

Jpt = ROOT.RooRealVar("pt", "p_{T}(J/#psi) [GeV]", ptBins[0], ptBins[-1])
Jrapidity = ROOT.RooRealVar("rapidity", "|y(J/#psi)|", yBins[0], yBins[-1])

data = ROOT.RooDataSet(
    "data", 
    "data", 
    ROOT.RooArgSet(Jmass, Jctau, Jpt, Jrapidity),
    ROOT.RooFit.Import(input_tree)
)

data = data.reduce(ROOT.RooFit.EventRange(0, 100000))

print(f">>> data.numEntries(): {data.numEntries()}")
# ================================================================
#  Mass PDF: Signal = f_G * DoubleGauss + (1-f_G) * DoubleCB
# ================================================================
mean = ROOT.RooRealVar("mean", "m_{J/\psi}", 3.0964, 3.08, 3.11)

gSigma1 = ROOT.RooRealVar("gSigma1", "gSigma1", 0.015, 0.010, 0.040)
gSigma2 = ROOT.RooRealVar("gSigma2", "gSigma2", 0.023, 0.013, 0.050)
pdf_gauss1 = ROOT.RooGaussian("pdf_gauss1", "pdf_gauss1", Jmass, mean, gSigma1)
pdf_gauss2 = ROOT.RooGaussian("pdf_gauss2", "pdf_gauss2", Jmass, mean, gSigma2)

cbSigmaL = ROOT.RooRealVar("cbSigmaL", "sigma Left", 0.016, 0.001, 0.035)
cbAlphaL = ROOT.RooRealVar("cbAlphaL", "alpha Left", 1.5, 0., 10.0) 
cbNL     = ROOT.RooRealVar("cbNL",     "n Left",     3.5, 1.1, 50.0)
pdf_cb_left = ROOT.RooCBShape("pdf_cb_left", "Left Tail CB", Jmass, mean, cbSigmaL, cbAlphaL, cbNL)

cbSigmaR = ROOT.RooRealVar("cbSigmaR", "sigma Right", 0.016, 0.01, 0.05)
cbAlphaR = ROOT.RooRealVar("cbAlphaR", "alpha Right", 1.5, 0., 5.0)  
cbNR     = ROOT.RooRealVar("cbNR",     "n Right",     3.5, 1.1, 50.0)
pdf_cb_right = ROOT.RooCBShape("pdf_cb_right", "Right Tail CB", Jmass, mean, cbSigmaR, cbAlphaR, cbNR)

fracG1 = ROOT.RooRealVar("fracG1", "f_{G,1}", 0.6, 0.0, 1.0)
pdf_doubleGauss = ROOT.RooAddPdf("pdf_doubleGauss", "pdf_doubleGauss", 
                                 ROOT.RooArgList(pdf_gauss1, pdf_gauss2), ROOT.RooArgList(fracG1))
fracCB1 = ROOT.RooRealVar("fracCB1", "fracCB1", 0.5, 0.0, 1.0)
pdf_doubleCB = ROOT.RooAddPdf("pdf_doubleCB", "pdf_doubleCB", 
                              ROOT.RooArgList(pdf_cb_left, pdf_cb_right), ROOT.RooArgList(fracCB1))

fracGaussInSig = ROOT.RooRealVar("fracGaussInSig", "fracGaussInSig", 0.8, 0.0, 1.0)
pdf_sigMass = ROOT.RooAddPdf("pdf_sigMass", "pdf_sigMass", 
                             ROOT.RooArgList(pdf_doubleGauss, pdf_doubleCB), ROOT.RooArgList(fracGaussInSig))
                             
# Background mass PDF
chebP1 = ROOT.RooRealVar("chebP1", "chebP1", -0.66076225, -2.0, 2.0)
chebP2 = ROOT.RooRealVar("chebP2", "chebP2", 0.030170337 , -2.0, 2.0)
chebP3 = ROOT.RooRealVar("chebP3", "chebP3", 0.053274777, -2.0, 2.0)
pdf_bkgMass = ROOT.RooChebychev("pdf_bkgMass", "pdf_bkgMass", Jmass, ROOT.RooArgList(chebP1, chebP2, chebP3))

# ================================================================
#  Ctau PDF
# ================================================================
# --- Prompt J/psi: Double Gaussian
ctRes_mean = ROOT.RooRealVar("ctRes_mean", "ctRes_mean", 0.0, -0.0001, 0.0001)
ctRes_sig1 = ROOT.RooRealVar("ctRes_sig1", "ctRes_sig1", 2.20688e-03, 0.0005, 0.050)
ctRes_sig2 = ROOT.RooRealVar("ctRes_sig2", "ctRes_sig2", 1.53310e-03, 0.0005, 0.010)
ctRes_sig3 = ROOT.RooRealVar("ctRes_sig3", "ctRes_sig3", 1.53310e-03, 0.0005, 0.010)
ctRes_frac1 = ROOT.RooRealVar("ctRes_frac1", "ctRes_frac1", 8.98616e-01, 0.0, 1.0)
ctRes_frac2 = ROOT.RooRealVar("ctRes_frac2", "ctRes_frac2", 8.98616e-01, 0.0, 1.0)

pdf_ctPrompt_g1 = ROOT.RooGaussian("pdf_ctPrompt_g1", "pdf_ctPrompt_g1", Jctau, ctRes_mean, ctRes_sig1)
pdf_ctPrompt_g2 = ROOT.RooGaussian("pdf_ctPrompt_g2", "pdf_ctPrompt_g2", Jctau, ctRes_mean, ctRes_sig2)
pdf_ctPrompt_g3 = ROOT.RooGaussian("pdf_ctPrompt_g3", "pdf_ctPrompt_g3", Jctau, ctRes_mean, ctRes_sig3)
pdf_ctPrompt_g12 = ROOT.RooAddPdf("pdf_ctPrompt_g12", "pdf_ctPrompt_g12", ROOT.RooArgList(pdf_ctPrompt_g1, pdf_ctPrompt_g2), ROOT.RooArgList(ctRes_frac1))
pdf_ctPrompt = ROOT.RooAddPdf("pdf_ctPrompt", "pdf_ctPrompt", ROOT.RooArgList(pdf_ctPrompt_g12, pdf_ctPrompt_g3), ROOT.RooArgList(ctRes_frac2))

# --- Non-prompt J/psi: 
ctRes_frac3 = ROOT.RooRealVar("ctRes_frac3", "ctRes_frac3", 8.98616e-01, 0.0, 1.0)
ctNonPrompt_tau = ROOT.RooRealVar("ctNonPrompt_tau", "ctNonPrompt_tau", 3.53161e-02, 0.005, 0.15)
resModel_g1 = ROOT.RooGaussModel("resModel_g1", "resModel_g1", Jctau, ctRes_mean, ctRes_sig1)
resModel_g2 = ROOT.RooGaussModel("resModel_g2", "resModel_g2", Jctau, ctRes_mean, ctRes_sig2)
resModel_sum = ROOT.RooAddModel("resModel_sum", "resModel_sum", ROOT.RooArgList(resModel_g1, resModel_g2), ROOT.RooArgList(ctRes_frac3))
pdf_ctNonPrompt = ROOT.RooDecay("pdf_ctNonPrompt", "pdf_ctNonPrompt", Jctau, ctNonPrompt_tau, resModel_sum, ROOT.RooDecay.SingleSided)

# ================================================================
#  Background ctau
# ================================================================
ctBkg_mean = ROOT.RooRealVar("ctBkg_mean", "ctBkg_mean", 0.0, -0.002, 0.002)

ctBkg_sig1 = ROOT.RooRealVar("ctBkg_sig1", "ctBkg_sig1", 0.0015, 0.0005, 0.0035)
ctBkg_sig2 = ROOT.RooRealVar("ctBkg_sig2", "ctBkg_sig2", 0.0050, 0.0036, 0.0200)
ctBkg_frac1 = ROOT.RooRealVar("ctBkg_frac1", "ctBkg_frac1", 0.75, 0.0, 1.0)

resModelBkg_g1 = ROOT.RooGaussModel("resModelBkg_g1", "Bkg Gauss Model 1", Jctau, ctBkg_mean, ctBkg_sig1)
resModelBkg_g2 = ROOT.RooGaussModel("resModelBkg_g2", "Bkg Gauss Model 2", Jctau, ctBkg_mean, ctBkg_sig2)

resModelBkg_sum = ROOT.RooAddModel(
    "resModelBkg_sum", "Bkg Combined Res Model", 
    ROOT.RooArgList(resModelBkg_g1, resModelBkg_g2), 
    ROOT.RooArgList(ctBkg_frac1)
)

pdf_ctBkg_prompt = resModelBkg_sum

ctBkg_tau_right = ROOT.RooRealVar("ctBkg_tau_right", "Bkg Tau Right", 0.03, 0.005, 0.2)
pdf_ctBkg_tailRight = ROOT.RooDecay(
    "pdf_ctBkg_tailRight", "Bkg Right Tail", 
    Jctau, ctBkg_tau_right, resModelBkg_sum, 
    ROOT.RooDecay.SingleSided
)

decay_Inverted = 1 
ctBkg_tau_left = ROOT.RooRealVar("ctBkg_tau_left", "Bkg Tau Left", 0.005, 0.0001, 0.02)
pdf_ctBkg_tailLeft = ROOT.RooDecay(
    "pdf_ctBkg_tailLeft", "Bkg Left Tail", 
    Jctau, ctBkg_tau_left, resModelBkg_sum, 
    decay_Inverted
)

ctBkg_f_prompt = ROOT.RooRealVar("ctBkg_f_prompt", "Bkg fraction of Prompt", 0.4, 0.0, 1.0)
ctBkg_f_right  = ROOT.RooRealVar("ctBkg_f_right",  "Bkg fraction of Right Tail", 0.5, 0.0, 1.0)

pdf_ctBkg = ROOT.RooAddPdf(
    "pdf_ctBkg", "Final Background ctau PDF", 
    ROOT.RooArgList(pdf_ctBkg_prompt, pdf_ctBkg_tailRight, pdf_ctBkg_tailLeft), 
    ROOT.RooArgList(ctBkg_f_prompt, ctBkg_f_right)
)

# ================================================================
#  2D total PDF: 3 components
# ================================================================
pdf_promptJpsi = ROOT.RooProdPdf("pdf_promptJpsi", "Prompt J/#psi (mass x ctau)", ROOT.RooArgList(pdf_sigMass, pdf_ctPrompt))
pdf_nonPromptJpsi = ROOT.RooProdPdf("pdf_nonPromptJpsi", "Non-prompt J/#psi (mass x ctau)", ROOT.RooArgList(pdf_sigMass, pdf_ctNonPrompt))
pdf_bkg = ROOT.RooProdPdf("pdf_bkg", "Combinatorial bkg (mass x ctau)", ROOT.RooArgList(pdf_bkgMass, pdf_ctBkg))

totalEntries = data.numEntries()
nPrompt = ROOT.RooRealVar("nPrompt", "nPrompt", totalEntries * 0.5, 0, totalEntries)
nNonPrompt = ROOT.RooRealVar("nNonPrompt", "nNonPrompt", totalEntries * 0.3, 0,  totalEntries)
nBkg = ROOT.RooRealVar("nBkg", "nBkg", totalEntries * 0.2, 0,  totalEntries)

pdf_total = ROOT.RooAddPdf("pdf_total", "Total 2D PDF", 
                            ROOT.RooArgList(pdf_promptJpsi, pdf_nonPromptJpsi, pdf_bkg), 
                            ROOT.RooArgList(nPrompt, nNonPrompt, nBkg))

cbAlphaL.setVal(3.5727402 )
cbAlphaR.setVal(0.22160347 )
cbNL.setVal(2.4127218 )
cbNR.setVal(33.993631 )
cbSigmaL.setVal(0.024216187 )
cbSigmaR.setVal(0.014339921 )
chebP1.setVal(-0.62088706 )
chebP2.setVal(-0.017403616 )
chebP3.setVal(0.032699596 )
ctBkg_f_prompt.setVal(0.16179959 )
ctBkg_f_right.setVal(0.75756908 )
ctBkg_frac1.setVal(0.46803348 )
ctBkg_mean.setVal(0.00019331703 )
ctBkg_sig1.setVal(0.0014130729 )
ctBkg_sig2.setVal(0.0036017046 )
ctBkg_tau_left.setVal(0.008660065 )
ctBkg_tau_right.setVal(0.038080983 )
ctNonPrompt_tau.setVal(0.037525516 )
ctRes_frac1.setVal(0.064788014 )
ctRes_frac2.setVal(0.6294201 )
ctRes_frac3.setVal(1.7585673e-11 )
ctRes_mean.setVal(7.2035145e-05 )
ctRes_sig1.setVal(0.0058924935 )
ctRes_sig2.setVal(0.0014441716 )
ctRes_sig3.setVal(0.0026387308 )
fracCB1.setVal(0.14238574 )
fracG1.setVal(0.52198676 )
fracGaussInSig.setVal(0.87308713 )
gSigma1.setVal(0.039999983 )
gSigma2.setVal(0.020750125 )
mean.setVal(3.0948056 )
nBkg.setVal(21298.862 )
nNonPrompt.setVal(31862.065 )
nPrompt.setVal(46838.984 )

fitResult = pdf_total.fitTo(data,
                            ROOT.RooFit.Extended(ROOT.kTRUE),
                            ROOT.RooFit.Save(ROOT.kTRUE),
                            ROOT.RooFit.Minimizer("Minuit2", "Migrad"),
                            ROOT.RooFit.Strategy(0),
                            ROOT.RooFit.NumCPU(4)
                            )

fitEDM = fitResult.edm()
fitStatus = fitResult.status()
fitCovQ = fitResult.covQual()

print("==========================================")
print(f"EDM: {fitEDM}")
print(f"Fit status: {fitStatus}")
print(f"Covariance quality: {fitCovQ}")
fitResult.Print()

for par in fitResult.floatParsFinal():
    print(f"{par.GetName()}.setVal({par.getVal():.8g} )")


# ================================================================
#  sPlot: prompt J/psi pT-y distribution
# ================================================================
sData = ROOT.RooStats.SPlot(
    "sData",
    "sData",
    data,
    pdf_total,
    ROOT.RooArgList(nPrompt, nNonPrompt, nBkg)
)

hPromptPtY = ROOT.TH2D(
    "hPromptPtY",
    "Prompt J/#psi; p_{T}(J/#psi) [GeV]; |y(J/#psi)|",
    len(ptBins) - 1,
    ptBins,
    len(yBins) - 1,
    yBins
)

for i in range(data.numEntries()):
    row = data.get(i)
    pt = row.getRealValue("pt")
    rapidity = row.getRealValue("rapidity")
    promptWeight = row.getRealValue("nPrompt_sw")

    hPromptPtY.Fill(pt, rapidity, promptWeight)

hPromptPt = hPromptPtY.ProjectionX("hPromptPt")
hPromptY = hPromptPtY.ProjectionY("hPromptY")

# ================================================================
#  Fit cross-check plots
# ================================================================

plotMass = Jmass.frame(
    ROOT.RooFit.Name("plotMass"),
    ROOT.RooFit.Title("J/#psi mass fit")
)

# Observed mass distribution.
data.plotOn(
    plotMass,
    ROOT.RooFit.Binning(100)
)

# Total 2D fit projected onto mass.
pdf_total.plotOn(
    plotMass,
    ROOT.RooFit.LineColor(ROOT.kRed),
    ROOT.RooFit.LineWidth(2)
)

# Prompt component.
pdf_total.plotOn(
    plotMass,
    ROOT.RooFit.Components("pdf_promptJpsi"),
    ROOT.RooFit.LineColor(ROOT.kGreen + 2),
    ROOT.RooFit.LineStyle(ROOT.kDashed),
    ROOT.RooFit.LineWidth(2)
)

# Non-prompt component.
pdf_total.plotOn(
    plotMass,
    ROOT.RooFit.Components("pdf_nonPromptJpsi"),
    ROOT.RooFit.LineColor(ROOT.kBlue),
    ROOT.RooFit.LineStyle(ROOT.kDashed),
    ROOT.RooFit.LineWidth(2)
)

# Background component.
pdf_total.plotOn(
    plotMass,
    ROOT.RooFit.Components("pdf_bkg"),
    ROOT.RooFit.LineColor(ROOT.kMagenta),
    ROOT.RooFit.LineStyle(ROOT.kDashed),
    ROOT.RooFit.LineWidth(2)
)

plotCtau = Jctau.frame(
    ROOT.RooFit.Name("plotCtau"),
    ROOT.RooFit.Title("J/#psi c#tau fit")
)

# Observed ctau distribution.
data.plotOn(
    plotCtau,
    ROOT.RooFit.Binning(100)
)

# Total 2D fit projected onto ctau.
pdf_total.plotOn(
    plotCtau,
    ROOT.RooFit.LineColor(ROOT.kRed),
    ROOT.RooFit.LineWidth(2)
)

# Prompt component.
pdf_total.plotOn(
    plotCtau,
    ROOT.RooFit.Components("pdf_promptJpsi"),
    ROOT.RooFit.LineColor(ROOT.kGreen + 2),
    ROOT.RooFit.LineStyle(ROOT.kDashed),
    ROOT.RooFit.LineWidth(2)
)

# Non-prompt component.
pdf_total.plotOn(
    plotCtau,
    ROOT.RooFit.Components("pdf_nonPromptJpsi"),
    ROOT.RooFit.LineColor(ROOT.kBlue),
    ROOT.RooFit.LineStyle(ROOT.kDashed),
    ROOT.RooFit.LineWidth(2)
)

# Background component.
pdf_total.plotOn(
    plotCtau,
    ROOT.RooFit.Components("pdf_bkg"),
    ROOT.RooFit.LineColor(ROOT.kMagenta),
    ROOT.RooFit.LineStyle(ROOT.kDashed),
    ROOT.RooFit.LineWidth(2)
)
# ================================================================
#  Save output histograms and fit cross-check plots
# ================================================================
output_file = ROOT.TFile("fitOutput.root", "RECREATE")

hPromptPtY.Write()
hPromptPt.Write()
hPromptY.Write()
plotMass.Write()
plotCtau.Write()

output_file.Close()
