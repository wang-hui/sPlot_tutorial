using namespace RooFit;

const double MUON_MASS = 0.1056583745;

TLorentzVector makeP4(float px,float py,float pz){
    TLorentzVector p4;
    p4.SetXYZM(px,py,pz,MUON_MASS);
    return p4;
}

//==================================================
// trigger
//==================================================
bool passEventTrigger(const std::vector<std::string>& trigNames,
                      const std::vector<unsigned int>& trigRes){
    for(size_t i=0;i<trigNames.size();++i){
        if(trigNames[i].find("HLT_DoubleMu4_3_LowMass_v")!=std::string::npos
           && trigRes[i]==1) return true;
    }
    return false;
}

//==================================================
ROOT::RVec<bool> passTrigMatch(const ROOT::RVec<float>& mu1Idx,
                               const ROOT::RVec<float>& mu2Idx,
                               const ROOT::RVec<int>& muTrigMatch){
    ROOT::RVec<bool> pass(mu1Idx.size(),false);
    for(size_t i=0;i<mu1Idx.size();++i){
        int i1=(int)mu1Idx[i];
        int i2=(int)mu2Idx[i];
        pass[i]=(muTrigMatch[i1]==1 && muTrigMatch[i2]==1);
    }
    return pass;
}

ROOT::RVec<bool> passSoftMuon(const ROOT::RVec<float>& mu1Idx,
                              const ROOT::RVec<float>& mu2Idx,
                              const ROOT::RVec<int>& muSoftID){
    ROOT::RVec<bool> pass(mu1Idx.size(),false);
    for(size_t i=0;i<mu1Idx.size();++i){
        int i1=(int)mu1Idx[i];
        int i2=(int)mu2Idx[i];
        pass[i]=(muSoftID[i1]==1 && muSoftID[i2]==1);
    }
    return pass;
}

ROOT::RVec<bool> passVertexCL(const ROOT::RVec<float>& vtxCL){
    return vtxCL>=0.01;
}

ROOT::RVec<bool> passOppositeSign(const ROOT::RVec<float>& mu1Idx,
                                  const ROOT::RVec<float>& mu2Idx,
                                  const ROOT::RVec<float>& muCharge){
    ROOT::RVec<bool> pass(mu1Idx.size(),false);
    for(size_t i=0;i<mu1Idx.size();++i){
        int i1=(int)mu1Idx[i];
        int i2=(int)mu2Idx[i];
        pass[i]=(muCharge[i1]+muCharge[i2])==0;
    }
    return pass;
}

ROOT::RVec<bool> passKinematics(const ROOT::RVec<float>& mu1Idx,
                                const ROOT::RVec<float>& mu2Idx,
                                const ROOT::RVec<float>& px,
                                const ROOT::RVec<float>& py,
                                const ROOT::RVec<float>& pz){
    ROOT::RVec<bool> pass(mu1Idx.size(),false);
    for(size_t i=0;i<mu1Idx.size();++i){
        int i1=(int)mu1Idx[i];
        int i2=(int)mu2Idx[i];

        auto p1=makeP4(px[i1],py[i1],pz[i1]);
        auto p2=makeP4(px[i2],py[i2],pz[i2]);

         bool acc1 = (fabs(p1.Eta()) < 2.4 && p1.Pt() > 5.0);
         bool acc2 = (fabs(p2.Eta()) < 2.4 && p2.Pt() > 5.0);
        
        pass[i]=acc1&&acc2;
    }
    return pass;
}

ROOT::RVec<bool> passSeagull(const ROOT::RVec<float>& mu1Idx,
                             const ROOT::RVec<float>& mu2Idx,
                             const ROOT::RVec<float>& px,
                             const ROOT::RVec<float>& py,
                             const ROOT::RVec<float>& pz){
    ROOT::RVec<bool> pass(mu1Idx.size(),false);

    for(size_t i=0;i<mu1Idx.size();++i){
        int i1=(int)mu1Idx[i];
        int i2=(int)mu2Idx[i];

        auto p1=makeP4(px[i1],py[i1],pz[i1]);
        auto p2=makeP4(px[i2],py[i2],pz[i2]);

        float dphi=p1.Phi()-p2.Phi();
        if(dphi>M_PI)dphi-=2*M_PI;
        if(dphi<-M_PI)dphi+=2*M_PI;

        pass[i]=dphi>0;
    }
    return pass;
}

ROOT::RVec<float> getJPsiRapidity(const ROOT::RVec<float>& px,
                                 const ROOT::RVec<float>& py,
                                 const ROOT::RVec<float>& pz,
                                 const ROOT::RVec<float>& mass)
{
    ROOT::RVec<float> y(px.size());

    for(size_t i=0;i<px.size();++i){

        TLorentzVector p4;

        p4.SetXYZM(px[i],py[i],pz[i],mass[i]);

        y[i]=p4.Rapidity();
    }

    return y;
}

ROOT::RVec<float> getAbsY(const ROOT::RVec<float>& y)
{
    ROOT::RVec<float> out(y.size());

    for(size_t i=0;i<y.size();++i)
        out[i]=std::abs(y[i]);

    return out;
}

//==================================================
// MAIN
//==================================================
void mumumass_efficiency(){
ROOT::EnableImplicitMT();
ROOT::RDataFrame df(
    "mkcands/X_data",
    "/eos/user/y/yingyinw/ntuple_practice/CMSSW_13_0_13/src/UserCode/MultiLepPAT/test/Run_Crab_Script/crab_Official_MC_JpsiPT5_MultiLepPAT_v4/results/*.root");
std::cout<<"Total events = "<<*df.Count()<<std::endl;
//ROOT::RDataFrame df("mkcands/X_data","../ntuple_practice/CMSSW_13_0_13/src/UserCode/MultiLepPAT/test/mymultilep.root");

auto df_gen =
     //remove the filter name so RdataFrame report ignores this filter
     //df.Filter("Offical_Pythia_J_MC_xPx.size()>0", "pass gen")
     df.Filter("Offical_Pythia_J_MC_xPx.size()>0")
     .Define("genJpsiPt",
         "sqrt(Offical_Pythia_J_MC_xPx*Offical_Pythia_J_MC_xPx + Offical_Pythia_J_MC_xPy*Offical_Pythia_J_MC_xPy)")
     .Define("genJpsiY",
         getJPsiRapidity,
         {"Offical_Pythia_J_MC_xPx",
         "Offical_Pythia_J_MC_xPy",
         "Offical_Pythia_J_MC_xPz",
         "Offical_Pythia_J_MC_xM"})
     .Define("genJpsiAbsY",getAbsY,{"genJpsiY"});
double yBins[] = {0.0,0.3,0.6,0.9,1.2,1.8,2.4};
double ptBins[]={10.0,10.5,11.0,11.5,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.0,32.0,34.0,36.0,38.0,42.0,46.0,50.0,60.0,75.0,95.0,120.0,200};
ROOT::RDF::TH2DModel hModel(
    "h2",
    "Dimuon Efficiency;|y|;p_{T} (GeV)",
    sizeof(yBins) / sizeof(double) - 1,
    yBins,
    sizeof(ptBins) / sizeof(double) - 1,
    ptBins
);

auto df_base=df_gen.Define("all_masses","mumuonlyMass")
.Define("jpsiPt","sqrt(mumuonlyPx*mumuonlyPx + mumuonlyPy*mumuonlyPy)")
.Define("jpsiY",getJPsiRapidity,
{"mumuonlyPx","mumuonlyPy","mumuonlyPz","mumuonlyMass"})
.Define("jpsiAbsY",getAbsY,{"jpsiY"})
.Define("hDen_njpsi","jpsiPt.size()");

auto hDen = df_base.Histo2D(hModel, "genJpsiAbsY", "genJpsiPt");
auto hDen_njpsi = df_base.Histo1D(ROOT::RDF::TH1DModel("hDen_njpsi","hDen_njpsi", 10, 0, 10), "hDen_njpsi");

auto df_trig=df_base.Filter(passEventTrigger,
{"TrigNames","TrigRes"}, "trg pass");

auto df1=df_trig.Define("mask_1",passSoftMuon,
{"mumuonlymu1Idx","mumuonlymu2Idx","muIsPatSoftMuon"})
.Filter("Any(mask_1)", "ID");

auto df2=df1.Define("cut_vtx",passVertexCL,{"mumuonlyVtxCL"})
.Define("mask_2","mask_1 && cut_vtx")
.Filter("Any(mask_2)", "vtx");

auto df3=df2.Define("cut_os",passOppositeSign,
{"mumuonlymu1Idx","mumuonlymu2Idx","muCharge"})
.Define("mask_3","mask_2 && cut_os")
.Filter("Any(mask_3)", "OS");

auto df3_trig=df3.Define("cut_trigMatch",passTrigMatch,
{"mumuonlymu1Idx","mumuonlymu2Idx","muIsJpsiTrigMatch"})
.Define("mask_3_trig","mask_3 && cut_trigMatch")
.Filter("Any(mask_3_trig)", "trg match");

auto df4=df3_trig.Define("cut_kin",passKinematics,
{"mumuonlymu1Idx","mumuonlymu2Idx","muPx","muPy","muPz"})
.Define("mask_4","mask_3_trig && cut_kin")
.Filter("Any(mask_4)", "pT&eta");

auto df5=df4.Define("cut_seagull",passSeagull,
{"mumuonlymu1Idx","mumuonlymu2Idx","muPx","muPy","muPz"})
.Define("mask_5","mask_4 && cut_seagull")
.Filter("Any(mask_5)", "seagull");

auto df_result = df5.Define("mask_final",
        "mask_5 && "
        "(all_masses >= 2.5 && all_masses <= 3.5) && "
        "(mumuonlyctau >= -0.03 && mumuonlyctau <= 0.12)")
    .Filter("Any(mask_final)", "m&ctau")
    .Define("good_pt","jpsiPt[mask_final]")
    .Define("good_y","jpsiAbsY[mask_final]")
    .Define("hNum_njpsi","good_pt.size()");

df_result.Report()->Print();

auto hNum = df_result.Histo2D(hModel, "genJpsiAbsY", "genJpsiPt");
auto hNum_njpsi = df_result.Histo1D(ROOT::RDF::TH1DModel("hNum_njpsi","hNum_njpsi", 10, 0, 10), "hNum_njpsi");

TH2D *hEff=(TH2D*)hNum.GetPtr()->Clone("hEff");

//================================
// Binomial error calculation
//================================

hEff->Reset();

for(int ix=1; ix<=hNum->GetNbinsX(); ix++)
{

    for(int iy=1; iy<=hNum->GetNbinsY(); iy++)
    {
        double Nnum_bin =
        hNum->GetBinContent(ix,iy);

        double Nden_bin =
        hDen->GetBinContent(ix,iy);

        if(Nden_bin > 0)
        {
            double eff =
            Nnum_bin / Nden_bin;
            double err =
            sqrt(eff*(1-eff)/Nden_bin);
            hEff->SetBinContent(ix,iy,eff);
            hEff->SetBinError(ix,iy,err);
        }
        else
        {
            hEff->SetBinContent(ix,iy,0);

            hEff->SetBinError(ix,iy,0);
        }
    }
}

//================================
// Save efficiency + error to txt
//================================

std::ofstream outfile("DimuonEfficiency.txt");

outfile
<< "# y_low y_high pt_low pt_high efficiency error"
<< std::endl;

for(int ix=1; ix<=hEff->GetNbinsX(); ix++)
{
    double yLow =hEff->GetXaxis()->GetBinLowEdge(ix);

    double yHigh =hEff->GetXaxis()->GetBinUpEdge(ix);

    for(int iy=1; iy<=hEff->GetNbinsY(); iy++)
    {
        double ptLow =
        hEff->GetYaxis()->GetBinLowEdge(iy);

        double ptHigh =hEff->GetYaxis()->GetBinUpEdge(iy);

        double eff =hEff->GetBinContent(ix,iy);

        double err =hEff->GetBinError(ix,iy);

        outfile
        << std::fixed
        << std::setprecision(3)
        << yLow << " "
        << yHigh << " "
        << ptLow << " "
        << ptHigh << " "
        << eff << " "
        << err
        << std::endl;
    }
}

outfile.close();
std::cout
<<"Efficiency table saved: DimuonEfficiency.txt"
<<std::endl;
//===============================
// calculate integrated efficiency
//===============================
double Nden = hDen->Integral();
double Nnum = hNum->Integral();
std::cout << "==========================" << std::endl;
std::cout << "Denominator integral = " << Nden
          << ", entries = " << hDen->GetEntries()
          << std::endl;

std::cout << "Numerator integral = " << Nnum
          << ", entries = " << hNum->GetEntries()
          << std::endl;

std::cout << "Candidate Efficiency = " << Nnum / Nden
          << ", entries efficiency = " << hNum->GetEntries() / hDen->GetEntries()
          << std::endl;

std::cout << "==========================" << std::endl;

//===============================
// Draw efficiency map
//===============================

TCanvas *c =new TCanvas("c","Dimuon Efficiency",1200,900);

c->SetLogy();
hDen_njpsi->Draw();
c->SaveAs("hDen_njpsi.png");
hNum_njpsi->Draw();
c->SaveAs("hNum_njpsi.png");

gStyle->SetOptStat(0);
gStyle->SetPaintTextFormat("4.3f");
c->SetLeftMargin(0.12);
c->SetRightMargin(0.16);
c->SetBottomMargin(0.12);
c->SetTopMargin(0.08);
hEff->SetStats(0);
//hEff->SetMinimum(0.05);
//hEff->SetMaximum(0.5);
hEff->GetXaxis()->SetTitle("|y|(J/#psi)");
hEff->GetYaxis()->SetTitle("p_{T}(J/#psi) [GeV]");
hEff->GetXaxis()->SetTitleSize(0.05);
hEff->GetYaxis()->SetTitleSize(0.05);
hEff->GetXaxis()->SetLabelSize(0.045);
hEff->GetYaxis()->SetLabelSize(0.045);
hEff->SetMarkerSize(1.1);
hEff->Draw("COLZ TEXT");
//c->SaveAs("DimuonEfficiency.pdf");
c->SaveAs("DimuonEfficiency.png");

TFile fout("DimuonEfficiency.root","RECREATE");
hNum->Write("hNum");
hEff->Write("hEff");
fout.Close();
std::cout<<"Efficiency map saved!"<<std::endl;
}
