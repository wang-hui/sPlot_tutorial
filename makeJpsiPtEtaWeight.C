using namespace RooFit;

void makeJpsiPtEtaWeight(
    const char* inputFile = "/nnucms/home/yide/storage2/2022F/Mass7/ParkingDoubleMuonLowMass7/RUN3MINI_2022FMass7/260424_110418/0000/mymultilep_99.root",
    const char* treeName = "mkcands/X_data",
    const char* outputFile = "JpsiPtEtaWeight.root"
)
{
    const double JPSI_MASS = 3.0969;
    const double MUON_MASS = 0.1056583745;

    const double MASS_MIN = 2.5;
    const double MASS_MAX = 3.5;
    const double CTAU_MIN = -0.03;
    const double CTAU_MAX = 0.12;

    const double PT_MIN = 10.;
    const double PT_MAX = 120.;
    const double ETA_MIN = 0.;
    const double ETA_MAX = 2.4;

    const int NPT = 55;
    const int NETA = 24;

    // =========================================================
    // RooFit observables
    // =========================================================

    RooRealVar mass("mass", "m_{#mu#mu} [GeV]", MASS_MIN, MASS_MAX);
    RooRealVar ctau("ctau", "c#tau", CTAU_MIN, CTAU_MAX);
    RooRealVar pt("pt", "p_{T}(J/#psi) [GeV]", PT_MIN, PT_MAX);
    RooRealVar eta("eta", "|#eta(J/#psi)|", ETA_MIN, ETA_MAX);

    RooArgSet vars(mass, ctau, pt, eta);

    // =========================================================
    // Read data with RDataFrame
    // =========================================================

    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df(treeName, inputFile);
    const unsigned int nSlots = df.GetNSlots();

    std::vector<std::unique_ptr<RooDataSet>> slotData;
    slotData.reserve(nSlots);

    for (unsigned int slot = 0; slot < nSlots; ++slot)
        slotData.emplace_back(
            std::make_unique<RooDataSet>(
                Form("data_slot%u", slot),
                Form("data_slot%u", slot),
                vars
            )
        );

    df.ForeachSlot(
        [&](unsigned int slot,
            const ROOT::RVec<float>& muPx,
            const ROOT::RVec<float>& muPy,
            const ROOT::RVec<float>& muPz,
            const ROOT::RVec<float>& mu1Idx,
            const ROOT::RVec<float>& mu2Idx,
            const ROOT::RVec<float>& mumuCtau) {

            for (size_t i = 0; i < mumuCtau.size(); ++i) {

                const int idx1 = static_cast<int>(mu1Idx[i]);
                const int idx2 = static_cast<int>(mu2Idx[i]);

                ROOT::Math::PxPyPzMVector mu1(
                    muPx[idx1], muPy[idx1], muPz[idx1], MUON_MASS
                );

                ROOT::Math::PxPyPzMVector mu2(
                    muPx[idx2], muPy[idx2], muPz[idx2], MUON_MASS
                );

                const auto jpsi = mu1 + mu2;

                const double thisMass = jpsi.M();
                const double thisPt = jpsi.Pt();
                const double thisEta = std::abs(jpsi.Eta());
                const double thisCtau = mumuCtau[i];

                if (thisMass < MASS_MIN || thisMass > MASS_MAX) continue;
                if (thisCtau < CTAU_MIN || thisCtau > CTAU_MAX) continue;
                if (thisPt < PT_MIN || thisPt > PT_MAX) continue;
                if (thisEta < ETA_MIN || thisEta > ETA_MAX) continue;

                mass.setVal(thisMass);
                ctau.setVal(thisCtau);
                pt.setVal(thisPt);
                eta.setVal(thisEta);

                slotData[slot]->add(vars);
            }
        },
        {
            "muPx",
            "muPy",
            "muPz",
            "mumuonlymu1Idx",
            "mumuonlymu2Idx",
            "mumuonlyctau"
        }
    );

    // Merge thread-local datasets.
    RooDataSet data("data", "J/#psi data", vars);

    for (const auto& d : slotData)
        data.append(*d);

    std::cout << "Candidates entering fit = "
              << data.numEntries() << std::endl;

    // =========================================================
    // Mass model
    // =========================================================

    RooRealVar massMean("massMean", "mass mean", 3.0969, 3.085, 3.11);
    RooRealVar massSigma("massSigma", "mass sigma", 0.025, 0.005, 0.10);
    RooRealVar massAlpha("massAlpha", "CB alpha", 1.5, 0.2, 5.);
    RooRealVar massN("massN", "CB n", 5., 0.5, 20.);

    RooCBShape signalMass(
        "signalMass", "J/#psi mass",
        mass, massMean, massSigma, massAlpha, massN
    );

    RooRealVar bkgMassSlope(
        "bkgMassSlope", "background mass slope",
        -1., -10., 0.
    );

    RooExponential bkgMass(
        "bkgMass", "background mass",
        mass, bkgMassSlope
    );

    // =========================================================
    // c*tau model
    // =========================================================

    RooRealVar promptCtauMean(
        "promptCtauMean", "prompt c#tau mean",
        0., -0.005, 0.005
    );

    RooRealVar promptCtauSigma(
        "promptCtauSigma", "prompt c#tau resolution",
        0.01, 0.001, 0.05
    );

    RooGaussian promptCtau(
        "promptCtau", "prompt c#tau",
        ctau, promptCtauMean, promptCtauSigma
    );

    RooRealVar ctauBias("ctauBias", "c#tau bias", 0.);
    ctauBias.setConstant(true);

    RooRealVar ctauRes(
        "ctauRes", "c#tau resolution",
        0.01, 0.001, 0.05
    );

    RooGaussModel resolutionModel(
        "resolutionModel", "resolution",
        ctau, ctauBias, ctauRes
    );

    RooRealVar nonPromptTau(
        "nonPromptTau", "non-prompt lifetime",
        0.045, 0.005, 0.20
    );

    RooDecay nonPromptCtau(
        "nonPromptCtau", "non-prompt c#tau",
        ctau, nonPromptTau, resolutionModel,
        RooDecay::SingleSided
    );

    RooRealVar bkgCtauMean(
        "bkgCtauMean", "background c#tau mean",
        0., -0.02, 0.02
    );

    RooRealVar bkgCtauSigma(
        "bkgCtauSigma", "background c#tau sigma",
        0.015, 0.002, 0.10
    );

    RooGaussian bkgCtauCore(
        "bkgCtauCore", "background c#tau core",
        ctau, bkgCtauMean, bkgCtauSigma
    );

    RooRealVar bkgCtauTau(
        "bkgCtauTau", "background c#tau tail",
        0.03, 0.005, 0.20
    );

    RooDecay bkgCtauTail(
        "bkgCtauTail", "background c#tau tail",
        ctau, bkgCtauTau, resolutionModel,
        RooDecay::DoubleSided
    );

    RooRealVar bkgCtauFraction(
        "bkgCtauFraction", "background core fraction",
        0.7, 0., 1.
    );

    RooAddPdf bkgCtau(
        "bkgCtau", "background c#tau",
        RooArgList(bkgCtauCore, bkgCtauTail),
        RooArgList(bkgCtauFraction)
    );

    // =========================================================
    // 2D fit components
    // =========================================================

    RooProdPdf promptPdf(
        "promptPdf", "prompt J/#psi",
        RooArgSet(signalMass, promptCtau)
    );

    RooProdPdf nonPromptPdf(
        "nonPromptPdf", "non-prompt J/#psi",
        RooArgSet(signalMass, nonPromptCtau)
    );

    RooProdPdf bkgPdf(
        "bkgPdf", "background",
        RooArgSet(bkgMass, bkgCtau)
    );

    // =========================================================
    // Extended model
    // =========================================================

    const double nData = data.numEntries();

    RooRealVar Nprompt(
        "Nprompt", "prompt yield",
        0.5 * nData, 0., 2. * nData
    );

    RooRealVar NnonPrompt(
        "NnonPrompt", "non-prompt yield",
        0.2 * nData, 0., 2. * nData
    );

    RooRealVar Nbkg(
        "Nbkg", "background yield",
        0.3 * nData, 0., 2. * nData
    );

    RooAddPdf model(
        "model",
        "prompt + non-prompt + background",
        RooArgList(promptPdf, nonPromptPdf, bkgPdf),
        RooArgList(Nprompt, NnonPrompt, Nbkg)
    );

    // =========================================================
    // Initial fit
    // =========================================================

    RooFitResult* fitResult = model.fitTo(
        data,
        Extended(true),
        Save(true),
        PrintLevel(-1)
    );

    fitResult->Print("v");

    // =========================================================
    // SPlot
    // =========================================================

    RooArgList yields(Nprompt, NnonPrompt, Nbkg);

    RooStats::SPlot splot(
        "splot", "splot", data, &model, yields
    );

    RooDataSet* sData = splot.GetSDataSet();

    // =========================================================
    // Output objects
    // =========================================================

    TH2D hPtEtaWeight(
        "hPtEtaWeight",
        "J/#psi p_{T}-|#eta| weight; p_{T} [GeV]; |#eta|",
        NPT, PT_MIN, PT_MAX,
        NETA, ETA_MIN, ETA_MAX
    );

    TH1D hPtProfile(
        "hPtProfile",
        "Prompt J/#psi p_{T}; p_{T} [GeV]; normalized yield",
        NPT, PT_MIN, PT_MAX
    );

    TH1D hEtaProfile(
        "hEtaProfile",
        "Prompt J/#psi |#eta|; |#eta|; normalized yield",
        NETA, ETA_MIN, ETA_MAX
    );

    // =========================================================
    // Prompt data distribution in pT-|eta|
    // =========================================================

    TH2D hDataPtEta(
        "hDataPtEta_tmp", "",
        NPT, PT_MIN, PT_MAX,
        NETA, ETA_MIN, ETA_MAX
    );

    for (int i = 0; i < sData->numEntries(); ++i) {

        const RooArgSet* row = sData->get(i);

        const double thisPt = row->getRealValue("pt");
        const double thisEta = row->getRealValue("eta");
        const double promptWeight = row->getRealValue("Nprompt_sw");

        hDataPtEta.Fill(thisPt, thisEta, promptWeight);
    }

    // =========================================================
    // Normalize data distribution
    // =========================================================

    const double dataIntegral = hDataPtEta.Integral();

    if (dataIntegral <= 0.) {
        std::cerr << "ERROR: prompt sWeighted yield <= 0" << std::endl;
        return;
    }

    hDataPtEta.Scale(1. / dataIntegral);

    const double gunBinProbability = 1. / (NPT * NETA);

    for (int ix = 1; ix <= NPT; ++ix) {
        for (int iy = 1; iy <= NETA; ++iy) {

            const double dataProbability =
                hDataPtEta.GetBinContent(ix, iy);

            const double weight =
                dataProbability / gunBinProbability;

            hPtEtaWeight.SetBinContent(ix, iy, weight);
        }
    }

    // =========================================================
    // 1D profiles
    // =========================================================

    for (int ix = 1; ix <= NPT; ++ix) {

        double sum = 0.;

        for (int iy = 1; iy <= NETA; ++iy)
            sum += hDataPtEta.GetBinContent(ix, iy);

        hPtProfile.SetBinContent(ix, sum);
    }

    for (int iy = 1; iy <= NETA; ++iy) {

        double sum = 0.;

        for (int ix = 1; ix <= NPT; ++ix)
            sum += hDataPtEta.GetBinContent(ix, iy);

        hEtaProfile.SetBinContent(iy, sum);
    }

    // =========================================================
    // Mass diagnostic
    // =========================================================

    RooPlot* plotMass = mass.frame(Title("J/#psi mass fit"));

    data.plotOn(plotMass, Binning(100));

    model.plotOn(plotMass, LineColor(kRed));

    model.plotOn(
        plotMass,
        Components(promptPdf),
        LineColor(kGreen + 2),
        LineStyle(kDashed)
    );

    model.plotOn(
        plotMass,
        Components(nonPromptPdf),
        LineColor(kBlue),
        LineStyle(kDashed)
    );

    model.plotOn(
        plotMass,
        Components(bkgPdf),
        LineColor(kMagenta),
        LineStyle(kDashed)
    );

    // =========================================================
    // c*tau diagnostic
    // =========================================================

    RooPlot* plotCtau = ctau.frame(Title("J/#psi c#tau fit"));

    data.plotOn(plotCtau, Binning(100));

    model.plotOn(plotCtau, LineColor(kRed));

    model.plotOn(
        plotCtau,
        Components(promptPdf),
        LineColor(kGreen + 2),
        LineStyle(kDashed)
    );

    model.plotOn(
        plotCtau,
        Components(nonPromptPdf),
        LineColor(kBlue),
        LineStyle(kDashed)
    );

    model.plotOn(
        plotCtau,
        Components(bkgPdf),
        LineColor(kMagenta),
        LineStyle(kDashed)
    );

    // =========================================================
    // Save histograms to root file
    // =========================================================

    TFile fout(outputFile, "RECREATE");

    hPtEtaWeight.Write();
    hPtProfile.Write();
    hEtaProfile.Write();
    plotMass->Write("plotMass");
    plotCtau->Write("plotCtau");

    fout.Close();

    std::cout << "Saved histograms to " << outputFile << std::endl;
}
