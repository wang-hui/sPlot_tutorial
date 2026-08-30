import ROOT
from DataFormats.FWLite import Events, Handle

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

inputFile = "step0_9.root"
refFile = "fitOutput.root"

# ----------------------------------------------------------------------
# Reference 2D distribution: prompt J/psi pT vs |y|
# ----------------------------------------------------------------------
fRef = ROOT.TFile.Open(refFile)
hRef = fRef.Get("hPromptPtY")
if not hRef:
    raise RuntimeError("Cannot find hPromptPtY")

genPartsMod = "genParticles"
genPartsType = Handle("vector<reco::GenParticle>")

# ----------------------------------------------------------------------
#input gen MC has to be run twice
#first time to get pT, |y| distribution to calculate weights
#second time to apply weights
#First pass:
# ----------------------------------------------------------------------
events = Events(inputFile)

nEvents = 0
nJpsi = 0

hGenPtY = hRef.Clone("hGenPtY")
hGenPtY.Reset()
hGenPtY.SetTitle("Original gen-level J/#psi;p_{T};|y|")

xMin = hRef.GetXaxis().GetXmin()
xMax = hRef.GetXaxis().GetXmax()
yMin = hRef.GetYaxis().GetXmin()
yMax = hRef.GetYaxis().GetXmax()

for event in events:
    nEvents += 1

    event.getByLabel(genPartsMod, genPartsType)
    genParts = genPartsType.product()

    jpsi = None
    for p in genParts:
        if abs(p.pdgId()) == 443:
            jpsi = p
            break

    if jpsi is None:
        continue

    pt = jpsi.pt()
    y = abs(jpsi.y())

    # Only consider the phase space defined by hPromptPtY.
    if not (xMin <= pt < xMax and yMin <= y < yMax):
        continue

    hGenPtY.Fill(pt, y)
    nJpsi += 1

refIntegral = hRef.Integral()
genIntegral = hGenPtY.Integral()

# ----------------------------------------------------------------------
# Construct the 2D weight map.
# ----------------------------------------------------------------------
hWeight = hRef.Clone("hPtYWeight")
hWeight.Reset()
hWeight.SetTitle("p_{T}-|y| reweighting factor;p_{T};|y|")

for ix in range(1, hRef.GetNbinsX() + 1):
    for iy in range(1, hRef.GetNbinsY() + 1):
        refFrac = hRef.GetBinContent(ix, iy) / refIntegral
        genFrac = hGenPtY.GetBinContent(ix, iy) / genIntegral

        # If the source sample has no events in a bin, the event weight
        # is undefined. Such a bin is assigned weight 0.
        if genFrac > 0:
            hWeight.SetBinContent(ix, iy, refFrac / genFrac)
        else:
            hWeight.SetBinContent(ix, iy, 0.0)

# ----------------------------------------------------------------------
# Event-level weight function.
# It takes pt and y and returns the pT-|y| reweighting factor.
# ----------------------------------------------------------------------
def getPtYWeight(pt, y):
    y = abs(y)

    if not (xMin <= pt < xMax and yMin <= y < yMax):
        return 0.0

    ix = hWeight.GetXaxis().FindFixBin(pt)
    iy = hWeight.GetYaxis().FindFixBin(y)

    return hWeight.GetBinContent(ix, iy)

# ----------------------------------------------------------------------
# Second pass:
# test the reweighting by drawing original and weighted pT and |y|.
# ----------------------------------------------------------------------

hPtOriginal = ROOT.TH1D(
    "hPtOriginal",
    "Gen-level J/#psi p_{T};p_{T};Events",
    20, xMin, xMax
)

hPtWeighted = ROOT.TH1D(
    "hPtWeighted",
    "Weighted gen-level J/#psi p_{T};p_{T};Weighted events",
    20, xMin, xMax
)

hYOriginal = ROOT.TH1D(
    "hYOriginal",
    "Gen-level J/#psi |y|;|y|;Events",
    20, yMin, yMax
)
hYWeighted = ROOT.TH1D(
    "hYWeighted",
    "Weighted gen-level J/#psi |y|;|y|;Weighted events",
    20, yMin, yMax
)

events = Events(inputFile)

nSelected = 0
sumWeights = 0.0

for event in events:
    event.getByLabel(genPartsMod, genPartsType)
    genParts = genPartsType.product()

    jpsi = None
    for p in genParts:
        if abs(p.pdgId()) == 443:
            jpsi = p
            break

    if jpsi is None:
        continue

    pt = jpsi.pt()
    y = abs(jpsi.y())

    # Same cut in the first pass.
    if not (xMin <= pt < xMax and yMin <= y < yMax):
        continue

    weight = getPtYWeight(pt, y)

    hPtOriginal.Fill(pt)
    hPtWeighted.Fill(pt, weight)

    hYOriginal.Fill(y)
    hYWeighted.Fill(y, weight)

    nSelected += 1
    sumWeights += weight

print(f"Total input events:       {nEvents}")
print(f"J/psi in reference range: {nJpsi}")
print(f"Selected events: {nSelected}")
print(f"Sum of weights: {sumWeights:.2f}")

# ----------------------------------------------------------------------
# Draw pT comparison.
# ----------------------------------------------------------------------
cPt = ROOT.TCanvas("cPt", "pT comparison", 800, 600)

hPtOriginal.SetLineColor(ROOT.kBlack)
hPtOriginal.SetLineWidth(2)

hPtWeighted.SetLineColor(ROOT.kRed)
hPtWeighted.SetLineWidth(2)

hPtWeighted.Draw("HIST")
hPtOriginal.Draw("HIST SAME")

legPt = ROOT.TLegend(0.62, 0.72, 0.88, 0.88)
legPt.AddEntry(hPtOriginal, "Original", "l")
legPt.AddEntry(hPtWeighted, "Weighted", "l")
legPt.Draw()

cPt.SaveAs("pt_weight_test.png")

# ----------------------------------------------------------------------
# Draw |y| comparison.
# ----------------------------------------------------------------------
cY = ROOT.TCanvas("cY", "|y| comparison", 800, 600)

hYOriginal.SetLineColor(ROOT.kBlack)
hYOriginal.SetLineWidth(2)

hYWeighted.SetLineColor(ROOT.kRed)
hYWeighted.SetLineWidth(2)

hYWeighted.Draw("HIST")
hYOriginal.Draw("HIST SAME")

legY = ROOT.TLegend(0.62, 0.72, 0.88, 0.88)
legY.AddEntry(hYOriginal, "Original", "l")
legY.AddEntry(hYWeighted, "Weighted", "l")
legY.Draw()

cY.SaveAs("y_weight_test.png")

