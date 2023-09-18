// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file k1analysis.cxx
/// \brief Reconstruction of track-track decay resonance candidates
///
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>

#include "TLorentzVector.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "DCAFitter/DCAFitterN.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/collisionCuts.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCalibration/MeanVertexObject.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct k1analysis {
  SliceCache cache;
  int mRunNumber;
  float d_bz;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::dataformats::MeanVertexObject* mVtx = nullptr;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  using ResoEvents = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::McCollisionLabels>;
  using ResoTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::pidTPCPi, aod::pidTPCKa, aod::pidTOFPi, aod::pidTOFKa, aod::TracksDCA, aod::McTrackLabels>;

  // CCDB options for magnetic field
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};

  Configurable<bool> cfgFatalWhenNull{"cfgFatalWhenNull", true, "Fatal when null"}; // Fatal when null, default true

  /// DCA fitter
  Configurable<double> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};
  Configurable<bool> d_UseWeightedPCA{"d_UseWeightedPCA", false, "Vertices use cov matrices"};
  Configurable<int> useMatCorrType{"useMatCorrType", 2, "0: none, 1: TGeo, 2: LUT"};
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};

  // Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.99, "V0 CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> v0minradius{"v0minradius", 0.0, "v0minradius"};
  Configurable<float> v0maxradius{"v0maxradius", 5.0, "v0maxradius"};
  Configurable<float> maxV0DCAtoPV{"maxV0DCAtoPV", 0.5, "maximum V0 DCA to PV"};

  /// Event cuts
  o2::analysis::CollisonCuts colCuts;
  Configurable<bool> ConfIsRun3{"ConfIsRun3", true, "Running on the Run3 data"}; // Choose if running on converted data or pilot beam
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};

  /// Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  /// Pre-selection cuts
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  /// DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};
  /// PID Selections
  Configurable<double> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 2.0, "TPC nSigma cut for Pion"};                    // TPC
  Configurable<double> cMaxTPCnSigmaPion_bach{"cMaxTPCnSigmaPion_bach", 2.0, "TPC nSigma cut for bachelor Pion"}; // TPC
  Configurable<double> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 2.0, "TOF nSigma cut for Pion"};                    // TOF
  Configurable<double> cMaxTOFnSigmaPion_bach{"cMaxTOFnSigmaPion_bach", 2.0, "TOF nSigma cut for Bachelor Pion"}; // TOF
  // Kaon
  Configurable<std::vector<double>> kaonTPCPIDpTintv{"kaonTPCPIDpTintv", {999.}, "pT intervals for Kaon TPC PID cuts"};
  Configurable<std::vector<double>> kaonTPCPIDcuts{"kaonTPCPIDcuts", {2}, "nSigma list for Kaon TPC PID cuts"};
  Configurable<std::vector<double>> kaonTOFPIDpTintv{"kaonTOFPIDpTintv", {999.}, "pT intervals for Kaon TOF PID cuts"};
  Configurable<std::vector<double>> kaonTOFPIDcuts{"kaonTOFPIDcuts", {2}, "nSigma list for Kaon TOF PID cuts"};

  // bachelor pion TOF PID?
  Configurable<int> cDoTOFPID{"cDoTOFPID", 1, "Do TOF PID"};

  // K(892)0 selection
  Configurable<double> cK892masswindow{"cK892masswindow", 0.1, "K(892)0 inv mass selection window"};
  Configurable<double> cPiPiMin{"cPiPiMin", 0, "Pion pair inv mass selection minimum"};
  Configurable<double> cPiPiMax{"cPiPiMax", 999, "Pion pair inv mass selection maximum"};

  // K1 selection
  Configurable<double> cK1MaxRap{"cK1MaxRap", 0.5, "K1 maximum rapidity"};
  Configurable<double> cK1MinRap{"cK1MinRap", -0.5, "K1 minimum rapidity"};

  // Define o2 fitter, 2-prong, active memory (no need to redefine per event)
  o2::vertexing::DCAFitterN<2> fitter;
  o2::track::TrackParCov lTrk1;
  o2::track::TrackParCov lTrk2;

  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter collisionFilter = nabs(aod::collision::posZ) < ConfEvtZvtx;
  Filter trackCutFilter = requireGlobalTrackInFilter();
  Filter trackFilter = nabs(aod::track::eta) < cfgCutEta; // Eta cut

  void init(o2::framework::InitContext&)
  {
    mRunNumber = 0;
    d_bz = 0;
    colCuts.setCuts(ConfEvtZvtx, ConfEvtTriggerCheck, ConfEvtTriggerSel, ConfEvtOfflineCheck, ConfIsRun3);
    colCuts.init(&histos);

    ccdb->setURL(ccdburl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(cfgFatalWhenNull);
    uint64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    std::vector<double> centBinning = {0., 1., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 80., 90., 100., 200.};
    AxisSpec centAxis = {centBinning, "T0M (%)"};
    AxisSpec ptAxis = {200, 0, 20, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec dcaxyAxis = {50, 0, 0.5, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {50, 0, 0.5, "DCA_{#it{z}} (cm)"};
    AxisSpec dcadauAxis = {200, 0, 2.0, "DCA_{daugthers} (cm)"};
    AxisSpec radiusAxis = {150, 0, 15.0, "radius (cm)"};
    AxisSpec cpaAxis = {500, 0.95, 1.0, "CPA"};
    AxisSpec invMassAxis = {900, 0.6, 1.5, "Invariant Mass (GeV/#it{c}^2)"};        // K(892)0
    AxisSpec invMassAxisReso = {1600, 0.9f, 2.5f, "Invariant Mass (GeV/#it{c}^2)"}; // K1
    AxisSpec invMassAxisScan = {250, 0, 2.5, "Invariant Mass (GeV/#it{c}^2)"};      // For selection
    AxisSpec pidQAAxis = {130, -6.5, 6.5};
    AxisSpec dataTypeAxis = {9, 0, 9, "Histogram types"};
    AxisSpec mcTypeAxis = {4, 0, 4, "Histogram types"};

    // Mass QA (quick check)
    histos.add("k892invmass", "Invariant mass of K(892)0", HistType::kTH1F, {invMassAxis});
    histos.add("k1invmass", "Invariant mass of K1(1270)pm", HistType::kTH1F, {invMassAxisReso});
    histos.add("k1invmass_LS", "Invariant mass of K1(1270)pm", HistType::kTH1F, {invMassAxisReso});
    histos.add("k1invmass_Mix", "Invariant mass of K1(1270)pm", HistType::kTH1F, {invMassAxisReso});
    // DCA QA
    histos.add("QAbefore/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAbefore/trkDCAxy_ka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAbefore/trkDCAxy_pi_bach", "DCAxy distribution of bachelor pion track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAbefore/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAbefore/trkDCAz_ka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAbefore/trkDCAz_pi_bach", "DCAz distribution of bachelor pion track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAafter/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAafter/trkDCAxy_ka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAafter/trkDCAxy_pi_bach", "DCAxy distribution of bachelor pion track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAafter/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAafter/trkDCAz_ka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAafter/trkDCAz_pi_bach", "DCAz distribution of bachelor pion track candidates", HistType::kTH1F, {dcazAxis});
    // pT QA
    histos.add("QAbefore/trkpT_pi", "pT distribution of pion track candidates", HistType::kTH1F, {ptAxis});
    histos.add("QAbefore/trkpT_ka", "pT distribution of kaon track candidates", HistType::kTH1F, {ptAxis});
    histos.add("QAbefore/trkpT_pi_bach", "pT distribution of bachelor pion track candidates", HistType::kTH1F, {ptAxis});
    histos.add("QAafter/trkpT_pi", "pT distribution of pion track candidates", HistType::kTH1F, {ptAxis});
    histos.add("QAafter/trkpT_ka", "pT distribution of kaon track candidates", HistType::kTH1F, {ptAxis});
    histos.add("QAafter/trkpT_pi_bach", "pT distribution of bachelor pion track candidates", HistType::kTH1F, {ptAxis});
    // PID QA before cuts
    histos.add("QAbefore/TOF_TPC_Map_pi", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QAbefore/TOF_Nsigma_pi", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAbefore/TPC_Nsigma_pi", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAbefore/TOF_TPC_Map_ka", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QAbefore/TOF_Nsigma_ka", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAbefore/TPC_Nsigmaka", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAbefore/TOF_TPC_Map_pi_bach", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QAbefore/TOF_Nsigma_pi_bach", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAbefore/TPC_Nsigma_pi_bach", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    // PID QA after cuts
    histos.add("QAafter/TOF_TPC_Map_pi", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QAafter/TOF_Nsigma_pi", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAafter/TPC_Nsigma_pi", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAafter/TOF_TPC_Map_ka", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QAafter/TOF_Nsigma_ka", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAafter/TPC_Nsigmaka", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAafter/TOF_TPC_Map_pi_bach", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QAafter/TOF_Nsigma_pi_bach", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAafter/TPC_Nsigma_pi_bach", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});

    histos.add("QAMCbefore/InvMass_piK_pipi", "Invariant mass of pion + kaon and pion+pion;Invariant Mass (GeV/#it{c}^{2});Invariant Mass (GeV/#it{c}^{2});", {HistType::kTH2F, {invMassAxisScan, invMassAxisScan}});

    // DCA V0
    histos.add("QAbefore/K892Radius", "K892 vtx position", HistType::kTH1F, {radiusAxis});
    histos.add("QAbefore/K892DCAdau", "DCA between K892 daughters", HistType::kTH1F, {dcadauAxis});
    histos.add("QAbefore/K892CPA", "K892 Cosine Pointing Angle", HistType::kTH1F, {cpaAxis});
    histos.add("QAafter/K892Radius", "K892 vtx position", HistType::kTH1F, {radiusAxis});
    histos.add("QAafter/K892DCAdau", "DCA between K892 daughters", HistType::kTH1F, {dcadauAxis});
    histos.add("QAafter/K892CPA", "K892 Cosine Pointing Angle", HistType::kTH1F, {cpaAxis});

    // Invariant mass histograms
    histos.add("hK892invmass_MM", "Invariant mass of K(892)0 (Matter + Matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxis});
    histos.add("hK892invmass_MA", "Invariant mass of K(892)0 (Matter + Anti-matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxis});
    histos.add("hK892invmass_AM", "Invariant mass of K(892)0 (Anti-matter + Matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxis});
    histos.add("hK892invmass_AA", "Invariant mass of K(892)0 (Anti-matter + Anti-matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxis});
    histos.add("hK1invmass_MM", "Invariant mass of K(892)0 + pion (Matter + Matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
    histos.add("hK1invmass_MA", "Invariant mass of K(892)0 + pion (Matter + Anti-matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
    histos.add("hK1invmass_AM", "Invariant mass of K(892)0 + pion (Anti-matter + Matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
    histos.add("hK1invmass_AA", "Invariant mass of K(892)0 + pion (Anti-matter + Anti-matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});

    histos.add("hK1invmass_MM_Mix", "Invariant mass of K(892)0 + pion (Matter + Matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
    histos.add("hK1invmass_MA_Mix", "Invariant mass of K(892)0 + pion (Matter + Anti-matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
    histos.add("hK1invmass_AM_Mix", "Invariant mass of K(892)0 + pion (Anti-matter + Matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
    histos.add("hK1invmass_AA_Mix", "Invariant mass of K(892)0 + pion (Anti-matter + Anti-matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});

    histos.add("QAMC/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAMC/trkDCAxy_ka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAMC/trkDCAxy_pi_bach", "DCAxy distribution of bachelor pion track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAMC/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAMC/trkDCAz_ka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAMC/trkDCAz_pi_bach", "DCAz distribution of bachelor pion track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAMC/TOF_TPC_Map_pi", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QAMC/TOF_Nsigma_pi", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAMC/TPC_Nsigma_pi", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAMC/TOF_TPC_Map_ka", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QAMC/TOF_Nsigma_ka", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAMC/TPC_Nsigmaka", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAMC/TOF_TPC_Map_pi_bach", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QAMC/TOF_Nsigma_pi_bach", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAMC/TPC_Nsigma_pi_bach", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAMC/K892Radius", "K892 vtx position", HistType::kTH1F, {radiusAxis});
    histos.add("QAMC/K892DCAdau", "DCA between K892 daughters", HistType::kTH1F, {dcadauAxis});
    histos.add("QAMC/K892CPA", "K892 Cosine Pointing Angle", HistType::kTH1F, {cpaAxis});
    histos.add("k892invmassMC", "Invariant mass of K(892)0", HistType::kTH1F, {invMassAxis});
    histos.add("k1invmassMC", "Invariant mass of K1(1270)pm", HistType::kTH1F, {invMassAxisReso});

    // initialize O2 2-prong fitter (only once)
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(d_UseAbsDCA);
    fitter.setWeightedFinalPCA(d_UseWeightedPCA);

    // Material correction in the DCA fitter
    o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
    if (useMatCorrType == 1)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    if (useMatCorrType == 2)
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    fitter.setMatCorrType(matCorr);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc) // Simple copy from LambdaKzeroFinder.cxx
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      fitter.setBz(d_bz);
      o2::parameters::GRPMagField grpmag;
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mVtx = ccdb->getForTimeStamp<o2::dataformats::MeanVertexObject>(mVtxPath, bc.timestamp());
    mRunNumber = bc.runNumber();
    // Set magnetic field value once known
    fitter.setBz(d_bz);
    LOGF(info, "Bz set to %f for run: ", d_bz, mRunNumber);
  }

  double massKa = TDatabasePDG::Instance()->GetParticle(kKPlus)->Mass();
  double massPi = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();
  double massK892 = TDatabasePDG::Instance()->GetParticle(313)->Mass();

  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
    // basic track cuts
    if (track.pt() < cMinPtcut)
      return false;
    if (track.dcaXY() > cMaxDCArToPVcut)
      return false;
    if (track.dcaZ() < cMinDCAzToPVcut || track.dcaZ() > cMaxDCAzToPVcut)
      return false;
    return true;
  }

  void processData(soa::Filtered<ResoEvents>::iterator const& collision,
                   soa::Filtered<ResoTracks> const& tracks,
                   aod::BCsWithTimestamps const&,
                   aod::McCollisions const& mcCols, aod::McParticles const& mcParticles)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    initCCDB(bc);
    // Default event selection
    if (!colCuts.isSelected(collision)) {
      LOG(info) << "Event rejected by cuts";
      return;
    }
    colCuts.fillQA(collision);

    o2::dataformats::VertexBase primaryVertex;
    primaryVertex.setPos({collision.posX(), collision.posY(), collision.posZ()});
    primaryVertex.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());

    auto getMothersIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lMothersIndeces{};
      for (auto& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother index lMother: %d", lMother.globalIndex());
        lMothersIndeces.push_back(lMother.globalIndex());
      }
      return lMothersIndeces;
    };
    auto getMothersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lMothersPDGs{};
      for (auto& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother pdgcode lMother: %d", lMother.pdgCode());
        lMothersPDGs.push_back(lMother.pdgCode());
      }
      return lMothersPDGs;
    };

    std::vector<int> mothers1;
    std::vector<int> mothers2;
    std::vector<int> mothers3;
    std::vector<int> motherPDGs1;
    std::vector<int> motherPDGs2;
    std::vector<int> motherPDGs3;

    auto vKaonTPCPIDpTintv = static_cast<std::vector<double>>(kaonTPCPIDpTintv);
    auto vKaonTPCPIDcuts = static_cast<std::vector<double>>(kaonTPCPIDcuts);
    auto vKaonTOFPIDpTintv = static_cast<std::vector<double>>(kaonTOFPIDpTintv);
    auto vKaonTOFPIDcuts = static_cast<std::vector<double>>(kaonTOFPIDcuts);
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonanceK892, lDecayDaughter_bach, lResonanceK1;
    auto lengthOfkaonTPCPIDpTintv = static_cast<int>(vKaonTPCPIDpTintv.size());
    // Loop over tracks
    for (auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {
      // Full index policy is needed to consider all possible combinations
      if (trk1.index() == trk2.index())
        continue; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.
      //// Initialize variables
      // Trk1: Pion, Trk2: Kaon
      auto isTrk1Selected = true;
      auto isTrk2Selected = true;
      auto isTrk1hasTOF = false;
      auto isTrk2hasTOF = false;
      auto trk1ptPi = trk1.pt();
      auto trk1NSigmaPiTPC = trk1.tpcNSigmaPi();
      auto trk1NSigmaPiTOF = -999.;
      auto trk2ptKa = trk2.pt();
      auto trk2NSigmaKaTPC = trk2.tpcNSigmaKa();
      auto trk2NSigmaKaTOF = -999.;

      // apply the track cut
      if (!trackCut(trk1))
        isTrk1Selected = false;
      if (!trackCut(trk2))
        isTrk2Selected = false;

      // TPC PID Selection
      if (std::abs(trk1NSigmaPiTPC) > cMaxTPCnSigmaPion)
        isTrk1Selected = false;
      // For Kaon candidate, we need to apply pT-dependent PID cuts
      if (lengthOfkaonTPCPIDpTintv > 0) {
        for (int i = 0; i < lengthOfkaonTPCPIDpTintv; i++) {
          if (trk2ptKa < vKaonTPCPIDpTintv[i]) {
            if (std::abs(trk2NSigmaKaTPC) > vKaonTPCPIDcuts[i])
              isTrk2Selected = false;
          }
        }
      }
      if (trk1.hasTOF()) {
        isTrk1hasTOF = true;
        trk1NSigmaPiTOF = trk1.tofNSigmaPi();
      }
      if (trk2.hasTOF()) {
        isTrk2hasTOF = true;
        trk2NSigmaKaTOF = trk2.tofNSigmaKa();
        if (lengthOfkaonTPCPIDpTintv > 0) {
          for (int i = 0; i < lengthOfkaonTPCPIDpTintv; i++) {
            if (trk2ptKa < vKaonTOFPIDpTintv[i]) {
              if (std::abs(trk2NSigmaKaTOF) > vKaonTOFPIDcuts[i])
                isTrk2Selected = false;
            }
          }
        }
      }
      //// QA plots before the selection
      //  --- PID QA Pion
      histos.fill(HIST("QAbefore/TPC_Nsigma_pi"), trk1ptPi, trk1NSigmaPiTPC);
      if (isTrk1hasTOF) {
        histos.fill(HIST("QAbefore/TOF_Nsigma_pi"), trk1ptPi, trk1NSigmaPiTOF);
        histos.fill(HIST("QAbefore/TOF_TPC_Map_pi"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
      }
      //  --- PID QA Kaon
      histos.fill(HIST("QAbefore/TPC_Nsigmaka"), trk2ptKa, trk2NSigmaKaTPC);
      if (isTrk2hasTOF) {
        histos.fill(HIST("QAbefore/TOF_Nsigma_ka"), trk2ptKa, trk2NSigmaKaTOF);
        histos.fill(HIST("QAbefore/TOF_TPC_Map_ka"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
      }
      histos.fill(HIST("QAbefore/trkpT_pi"), trk1ptPi);
      histos.fill(HIST("QAbefore/trkpT_ka"), trk2ptKa);

      histos.fill(HIST("QAbefore/trkDCAxy_pi"), trk1.dcaXY());
      histos.fill(HIST("QAbefore/trkDCAxy_ka"), trk2.dcaXY());
      histos.fill(HIST("QAbefore/trkDCAz_pi"), trk1.dcaZ());
      histos.fill(HIST("QAbefore/trkDCAz_ka"), trk2.dcaZ());

      //// Apply the selection
      if (!isTrk1Selected || !isTrk2Selected)
        continue;

      // Try to progate to dca
      lTrk1 = getTrackParCov(trk1);
      lTrk2 = getTrackParCov(trk2);

      int nCand = fitter.process(lTrk1, lTrk2);
      if (nCand == 0) {
        continue;
      }
      const auto& vtx = fitter.getPCACandidate();

      // Fiducial: min radius
      auto sVtxRadius = RecoDecay::sqrtSumOfSquares(vtx[0], vtx[1]);
      auto thisdcav0dau = fitter.getChi2AtPCACandidate();

      std::array<float, 3> pos = {0.};
      std::array<float, 3> pvec0;
      std::array<float, 3> pvec1;
      for (int i = 0; i < 3; i++) {
        pos[i] = vtx[i];
      }
      fitter.getTrack(0).getPxPyPzGlo(pvec0);
      fitter.getTrack(1).getPxPyPzGlo(pvec1);

      auto cpa = RecoDecay::cpa(array{primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, array{pos[0], pos[1], pos[2]}, array{pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]});

      histos.fill(HIST("QAbefore/K892Radius"), sVtxRadius);
      histos.fill(HIST("QAbefore/K892DCAdau"), thisdcav0dau);
      histos.fill(HIST("QAbefore/K892CPA"), cpa);

      // Fiducial: max radius
      if (sVtxRadius > v0maxradius || sVtxRadius < v0minradius) {
        continue;
      }
      // DCA daugthers
      if (thisdcav0dau > dcav0dau) {
        continue;
      }
      if (cpa < v0cospa) {
        continue;
      }

      histos.fill(HIST("QAafter/TPC_Nsigma_pi"), trk1ptPi, trk1NSigmaPiTPC);
      if (isTrk1hasTOF) {
        histos.fill(HIST("QAafter/TOF_Nsigma_pi"), trk1ptPi, trk1NSigmaPiTOF);
        histos.fill(HIST("QAafter/TOF_TPC_Map_pi"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
      }
      //  --- PID QA Kaon
      histos.fill(HIST("QAafter/TPC_Nsigmaka"), trk2ptKa, trk2NSigmaKaTPC);
      if (isTrk1hasTOF) {
        histos.fill(HIST("QAafter/TOF_Nsigma_ka"), trk2ptKa, trk2NSigmaKaTOF);
        histos.fill(HIST("QAafter/TOF_TPC_Map_ka"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
      }
      histos.fill(HIST("QAafter/trkpT_pi"), trk1ptPi);
      histos.fill(HIST("QAafter/trkpT_ka"), trk2ptKa);

      histos.fill(HIST("QAafter/trkDCAxy_pi"), trk1.dcaXY());
      histos.fill(HIST("QAafter/trkDCAxy_ka"), trk2.dcaXY());
      histos.fill(HIST("QAafter/trkDCAz_pi"), trk1.dcaZ());
      histos.fill(HIST("QAafter/trkDCAz_ka"), trk2.dcaZ());

      histos.fill(HIST("QAafter/K892Radius"), sVtxRadius);
      histos.fill(HIST("QAafter/K892DCAdau"), thisdcav0dau);
      histos.fill(HIST("QAafter/K892CPA"), cpa);

      //// Resonance reconstruction
      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massPi);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
      lResonanceK892 = lDecayDaughter1 + lDecayDaughter2;

      histos.fill(HIST("k892invmass"), lResonanceK892.M()); // quick check
      if (trk1.sign() > 0) {                                // Positive pion
        if (trk2.sign() > 0)                                // Positive kaon
          histos.fill(HIST("hK892invmass_MM"), collision.centFT0M(), lResonanceK892.Pt(), lResonanceK892.M());
        else                                                                                                   // Negative kaon
          histos.fill(HIST("hK892invmass_AM"), collision.centFT0M(), lResonanceK892.Pt(), lResonanceK892.M()); // Anti-K(892)0
      } else {                                                                                                 // Negative pion
        if (trk2.sign() > 0)                                                                                   // Positive kaon
          histos.fill(HIST("hK892invmass_MA"), collision.centFT0M(), lResonanceK892.Pt(), lResonanceK892.M()); // K(892)0
        else                                                                                                   // Negative kaon
          histos.fill(HIST("hK892invmass_AA"), collision.centFT0M(), lResonanceK892.Pt(), lResonanceK892.M());
      }

      // Like-sign rejection
      if (trk1.sign() * trk2.sign() > 0)
        continue;

      // Mass window cut
      if (std::abs(lResonanceK892.M() - massK892) > cK892masswindow)
        continue;

      for (auto bTrack : tracks) {
        // ID cut
        if (bTrack.index() == trk1.index() || bTrack.index() == trk2.index())
          continue;
        auto isbTrkSelected = true;
        // Track cut
        if (!trackCut(bTrack))
          isbTrkSelected = false;

        auto bTrkPt = bTrack.pt();
        auto bTrkTPCnSigmaPi = bTrack.tpcNSigmaPi();
        auto bTrack_TOFnSigma = -99.0;
        auto isbTrkhasTOF = false;
        if (bTrack.hasTOF()) {
          bTrack_TOFnSigma = bTrack.tofNSigmaPi();
          isbTrkhasTOF = true;
        }
        // PID selection
        if (std::abs(bTrkTPCnSigmaPi) > cMaxTPCnSigmaPion_bach)
          isbTrkSelected = false;
        if (cDoTOFPID && isbTrkhasTOF) {
          if (std::abs(bTrack_TOFnSigma) > cMaxTOFnSigmaPion_bach)
            isbTrkSelected = false;
        }
        histos.fill(HIST("QAbefore/trkpT_pi_bach"), bTrkPt);
        //  --- PID QA Pion
        histos.fill(HIST("QAbefore/TPC_Nsigma_pi_bach"), bTrkPt, bTrkTPCnSigmaPi);
        if (isbTrkhasTOF) {
          histos.fill(HIST("QAbefore/TOF_Nsigma_pi_bach"), bTrkPt, bTrack_TOFnSigma);
          histos.fill(HIST("QAbefore/TOF_TPC_Map_pi_bach"), bTrack_TOFnSigma, bTrkTPCnSigmaPi);
        }
        histos.fill(HIST("QAbefore/trkDCAxy_pi_bach"), bTrack.dcaXY());
        histos.fill(HIST("QAbefore/trkDCAz_pi_bach"), bTrack.dcaZ());

        if (!isbTrkSelected) // bachelor track selection
          continue;

        histos.fill(HIST("QAafter/trkpT_pi_bach"), bTrkPt);
        //  --- PID QA Pion
        histos.fill(HIST("QAafter/TPC_Nsigma_pi_bach"), bTrkPt, bTrkTPCnSigmaPi);
        if (isbTrkhasTOF) {
          histos.fill(HIST("QAafter/TOF_Nsigma_pi_bach"), bTrkPt, bTrack_TOFnSigma);
          histos.fill(HIST("QAafter/TOF_TPC_Map_pi_bach"), bTrack_TOFnSigma, bTrkTPCnSigmaPi);
        }
        histos.fill(HIST("QAafter/trkDCAxy_pi_bach"), bTrack.dcaXY());
        histos.fill(HIST("QAafter/trkDCAz_pi_bach"), bTrack.dcaZ());

        // K1 reconstruction
        lDecayDaughter_bach.SetXYZM(bTrack.px(), bTrack.py(), bTrack.pz(), massPi);
        lResonanceK1 = lResonanceK892 + lDecayDaughter_bach;

        // Rapidity cut
        if (lResonanceK1.Rapidity() > cK1MaxRap || lResonanceK1.Rapidity() < cK1MinRap)
          continue;

        if (bTrack.sign() > 0) {                              // bachelor pi+
          if (trk2.sign() > 0) {                              // kaon + means K(892)0 is matter.
            histos.fill(HIST("k1invmass"), lResonanceK1.M()); // quick check
            histos.fill(HIST("hK1invmass_MM"), collision.centFT0M(), lResonanceK1.Pt(), lResonanceK1.M());
          } else {
            histos.fill(HIST("k1invmass_LS"), lResonanceK1.M()); // quick check
            histos.fill(HIST("hK1invmass_AM"), collision.centFT0M(), lResonanceK1.Pt(), lResonanceK1.M());
          }
        } else {                                                 // bachelor pi-
          if (trk2.sign() > 0) {                                 // kaon + means K(892)0 is matter.
            histos.fill(HIST("k1invmass_LS"), lResonanceK1.M()); // quick check
            histos.fill(HIST("hK1invmass_MA"), collision.centFT0M(), lResonanceK1.Pt(), lResonanceK1.M());
          } else {
            histos.fill(HIST("k1invmass"), lResonanceK1.M()); // quick check
            histos.fill(HIST("hK1invmass_AA"), collision.centFT0M(), lResonanceK1.Pt(), lResonanceK1.M());
          }
        }
        // MC CASE
        mothers1.clear();
        motherPDGs1.clear();
        mothers2.clear();
        motherPDGs2.clear();
        // only consider tracks with MC particle
        if (!trk1.has_mcParticle())
          continue;
        if (!trk2.has_mcParticle())
          continue;
        auto mcPart1 = mcParticles.rawIteratorAt(trk1.mcParticleId() - mcParticles.offset());
        auto mcPart2 = mcParticles.rawIteratorAt(trk2.mcParticleId() - mcParticles.offset());
        auto pdg1 = abs(mcPart1.pdgCode());
        auto pdg2 = abs(mcPart2.pdgCode());
        // First track is pion from K(892)0
        if (pdg1 != 211)
          continue;
        // Second track is Kaon from K(892)
        if (pdg2 != 321)
          continue;
        if (!mcPart1.has_mothers())
          continue;
        if (!mcPart2.has_mothers())
          continue;
        // Check mother PDG code
        mothers1 = getMothersIndeces(mcPart1);
        motherPDGs1 = getMothersPDGCodes(mcPart1);
        mothers2 = getMothersIndeces(mcPart2);
        motherPDGs2 = getMothersPDGCodes(mcPart2);
        if (motherPDGs1[0] != 313)
          continue;
        if (motherPDGs2[0] != 313)
          continue;
        // Check if the mother is the same
        if (mothers1[0] != mothers2[0])
          continue;
        histos.fill(HIST("k892invmassMC"), lResonanceK892.M()); // quick check

        // K1 reconstruction
        if (!bTrack.has_mcParticle())
          continue;
        mothers3.clear();
        motherPDGs3.clear();
        auto mcPart3 = mcParticles.rawIteratorAt(bTrack.mcParticleId() - mcParticles.offset());
        auto pdg3 = abs(mcPart3.pdgCode());
        if (pdg3 != 211)
          continue;
        if (!mcPart3.has_mothers())
          continue;
        mothers3 = getMothersIndeces(mcPart3);
        motherPDGs3 = getMothersPDGCodes(mcPart3);
        if (motherPDGs3[0] != 10323)
          continue;

        // Check the daughter of mother
        auto mcK1 = mcParticles.rawIteratorAt(mothers3[0] - mcParticles.offset());
        auto mcK1Daughters = mcK1.daughters_as<aod::McParticles>();
        if (mcK1Daughters.size() != 2)
          continue;
        auto mcK1Daughter1 = mcK1Daughters.iteratorAt(0);
        auto mcK1Daughter2 = mcK1Daughters.iteratorAt(1);
        auto mcK1Daughter1index = mcK1Daughter1.globalIndex();
        auto mcK1Daughter2index = mcK1Daughter2.globalIndex();
        if (mothers1[0] != mcK1Daughter1index && mothers1[0] != mcK1Daughter2index)
          continue;

        //  --- PID QA Pion
        histos.fill(HIST("QAMC/TPC_Nsigma_pi"), trk1ptPi, trk1NSigmaPiTPC);
        if (isTrk1hasTOF) {
          histos.fill(HIST("QAMC/TOF_Nsigma_pi"), trk1ptPi, trk1NSigmaPiTOF);
          histos.fill(HIST("QAMC/TOF_TPC_Map_pi"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
        }
        //  --- PID QA Kaon
        histos.fill(HIST("QAMC/TPC_Nsigmaka"), trk2ptKa, trk2NSigmaKaTPC);
        if (isTrk2hasTOF) {
          histos.fill(HIST("QAMC/TOF_Nsigma_ka"), trk2ptKa, trk2NSigmaKaTOF);
          histos.fill(HIST("QAMC/TOF_TPC_Map_ka"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
        }
        //  --- PID QA Pion
        histos.fill(HIST("QAMC/TPC_Nsigma_pi_bach"), bTrkPt, bTrkTPCnSigmaPi);
        if (isbTrkhasTOF) {
          histos.fill(HIST("QAMC/TOF_Nsigma_pi_bach"), bTrkPt, bTrack_TOFnSigma);
          histos.fill(HIST("QAMC/TOF_TPC_Map_pi_bach"), bTrack_TOFnSigma, bTrkTPCnSigmaPi);
        }
        histos.fill(HIST("QAMC/trkDCAxy_pi"), trk1.dcaXY());
        histos.fill(HIST("QAMC/trkDCAxy_ka"), trk2.dcaXY());
        histos.fill(HIST("QAMC/trkDCAz_pi"), trk1.dcaZ());
        histos.fill(HIST("QAMC/trkDCAz_ka"), trk2.dcaZ());
        histos.fill(HIST("QAMC/trkDCAxy_pi_bach"), bTrack.dcaXY());
        histos.fill(HIST("QAMC/trkDCAz_pi_bach"), bTrack.dcaZ());

        histos.fill(HIST("QAMC/K892Radius"), sVtxRadius);
        histos.fill(HIST("QAMC/K892DCAdau"), thisdcav0dau);
        histos.fill(HIST("QAMC/K892CPA"), cpa);

        histos.fill(HIST("k1invmassMC"), lResonanceK1.M()); // quick check
      }
    }
  }
  PROCESS_SWITCH(k1analysis, processData, "Process Event for data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<k1analysis>(cfgc, TaskName{"lf-k1analysis"})};
}
