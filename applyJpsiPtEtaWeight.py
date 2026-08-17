import ROOT
from DataFormats.FWLite import Events, Handle

#ROOT.gROOT.SetBatch(True)

# ============================================================
# Input/output
# ============================================================

inputFile = "step0_9.root"
weightFile = "JpsiPtEtaWeight.root"

genPartsMod = "genParticles"

# ============================================================
# Load the 2D weight histogram
# ============================================================

fWeight = ROOT.TFile.Open(weightFile)
hWeight = fWeight.Get("hPtEtaWeight")

if not hWeight:
    raise RuntimeError("Cannot find hPtEtaWeight")

# ============================================================
# GenParticle handle
# ============================================================

genPartsType = Handle("vector<reco::GenParticle>")

events = Events(inputFile)

# ============================================================
# Output histograms
# ============================================================

hPtBefore = ROOT.TH1D(
    "hPtBefore",
    "J/#psi p_{T};p_{T} [GeV];events",
    10, 0., 200.
)

hPtAfter = ROOT.TH1D(
    "hPtAfter",
    "J/#psi p_{T};p_{T} [GeV];events",
    10, 0., 200.
)

hEtaBefore = ROOT.TH1D(
    "hEtaBefore",
    "J/#psi #eta;#eta;events",
    10, -3., 3.
)

hEtaAfter = ROOT.TH1D(
    "hEtaAfter",
    "J/#psi #eta;#eta;events",
    10, -3., 3.
)

hPtBefore.Sumw2()
hPtAfter.Sumw2()
hEtaBefore.Sumw2()
hEtaAfter.Sumw2()

# ============================================================
# Event loop
# ============================================================

nEvents = 0
nJpsi = 0

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

    nJpsi += 1

    pt = jpsi.pt()
    eta = jpsi.eta()
    absEta = abs(eta)

    # --------------------------------------------------------
    # Before weighting
    # --------------------------------------------------------

    hPtBefore.Fill(pt)
    hEtaBefore.Fill(eta)

    # --------------------------------------------------------
    # Look up pT-|eta| weight
    # --------------------------------------------------------

    ix = hWeight.GetXaxis().FindBin(pt)
    iy = hWeight.GetYaxis().FindBin(absEta)

    weight = 1.0

    if (1 <= ix <= hWeight.GetNbinsX() and
        1 <= iy <= hWeight.GetNbinsY()):

        weight = hWeight.GetBinContent(ix, iy)

    # --------------------------------------------------------
    # After weighting
    # --------------------------------------------------------

    hPtAfter.Fill(pt, weight)
    hEtaAfter.Fill(eta, weight)

print("Events processed =", nEvents)
print("J/psi found      =", nJpsi)

# ============================================================
# Plot pT
# ============================================================

hPtBefore.SetLineColor(ROOT.kBlack)
hPtBefore.SetLineWidth(2)

hPtAfter.SetLineColor(ROOT.kRed)
hPtAfter.SetLineWidth(2)

cPt = ROOT.TCanvas("cPt", "J/psi pT", 800, 600)

hPtBefore.Draw("HIST")
hPtAfter.Draw("HIST SAME")

legPt = ROOT.TLegend(0.55, 0.72, 0.88, 0.88)
legPt.AddEntry(hPtBefore, "Particle gun", "l")
legPt.AddEntry(hPtAfter, "After p_{T}-|#eta| weighting", "l")
legPt.Draw()

cPt.SaveAs("JpsiPt_before_after.png")

# ============================================================
# Plot eta
# ============================================================

hEtaBefore.SetLineColor(ROOT.kBlack)
hEtaBefore.SetLineWidth(2)

hEtaAfter.SetLineColor(ROOT.kRed)
hEtaAfter.SetLineWidth(2)

cEta = ROOT.TCanvas("cEta", "J/psi eta", 800, 600)

hEtaBefore.Draw("HIST")
hEtaAfter.Draw("HIST SAME")

legEta = ROOT.TLegend(0.55, 0.72, 0.88, 0.88)
legEta.AddEntry(hEtaBefore, "Particle gun", "l")
legEta.AddEntry(hEtaAfter, "After p_{T}-|#eta| weighting", "l")
legEta.Draw()

cEta.SaveAs("JpsiEta_before_after.png")

