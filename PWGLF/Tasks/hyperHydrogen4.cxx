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
//
// 
//

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"

#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/TrackSelectorPID.h"
#include "PID/PIDResponse.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/StrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct hyperHydrogen4{
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    AxisSpec vtxZAxis = {100, -20, 20};

    std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};
    AxisSpec centAxis = {centBinning, "V0M (%)"};
    std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 10., 20.};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    histos.add("EventQA/Centrality", "Centrality distribution (V0M)", kTH1F, {centAxis});
    histos.add("EventQA/VtxZBeforeSel", "Vertex distribution in Z;Z (cm)", kTH1F, {vtxZAxis});
    histos.add("EventQA/VtxZAfterSel", "Vertex distribution in Z;Z (cm)", kTH1F, {vtxZAxis});

    // Mass QA (quick check) PDG value: 3.931
    histos.add("hyperHydrogen4invmass", "Invariant mass of Hyper Hydrogen 4", kTH1F, {{1400, 3.4, 4.6, "Invariant Mass (GeV/#it{c}^2)"}});
    histos.add("antihyperHydrogen4invmass", "Invariant mass of Anti-Hyper Hydrogen 4", kTH1F, {{1400, 3.4, 4.6, "Invariant Mass (GeV/#it{c}^2)"}});

    // 3d histogram
    histos.add("h3hyperHydrogen4invmass", "Invariant mass of Hyper Hydrogen 4", kTH3F, {{100, 0.0f, 100.0f}, {150, 0.0f, 15.0f}, {1400, 3.4, 4.6}});
    histos.add("h3antihyperHydrogen4invmass", "Invariant mass of Anti-Hyper Hydrogen 4", kTH3F, {{100, 0.0f, 100.0f}, {150, 0.0f, 15.0f}, {1400, 3.4, 4.6}});
  }

  // Configurables
  Configurable<double> d_bz{"d_bz", -5.0, "bz field"};

  /// DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut1{"cMaxDCArToPVcut1", 0.05, "Track1 DCAr cut to PV Maximum"};
  Configurable<double> cMaxDCArToPVcut2{"cMaxDCArToPVcut2", 0.05, "Track2 DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut1{"cMaxDCAzToPVcut1", 2.0, "Track1 DCAz cut to PV Maximum"};
  Configurable<double> cMaxDCAzToPVcut2{"cMaxDCAzToPVcut2", 2.0, "Track2 DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut1{"cMinDCAzToPVcut1", 0.0, "Track1 DCAz cut to PV Minimum"};
  Configurable<double> cMinDCAzToPVcut2{"cMinDCAzToPVcut2", 0.0, "Track2 DCAz cut to PV Minimum"};

  // Secondary vertex selections
  Configurable<float> cHyperCosPA{"cHyperCosPA", 0.995, "Hyper CosPA to PV"};
  Configurable<float> cMaxDCAHyper{"cMaxDCAHyper", 1.5, "Maximum DCA between hyper daughters"};
  Configurable<float> cMinHyperRadius{"cMinHyperRadius", 5.0, "Minimum hyper radius from PV"};

  /// PID selections
  // TPC
  Configurable<float> pidTPCMinpT{"pidTPCMinpT", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<float> pidTPCMaxpT{"pidTPCMaxpT", 10., "Upper bound of track pT for TPC PID"};
  Configurable<float> nSigmaTPC{"nSigmaTPC", 5., "Nsigma cut on TPC only"};
  Configurable<float> nSigmaTPCCombined{"nSigmaTPCCombined", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF
  Configurable<float> pidTOFMinpT{"pidTOFMinpT", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<float> pidTOFMaxpT{"pidTOFMaxpT", 10., "Upper bound of track pT for TOF PID"};
  Configurable<float> nSigmaTOF{"nSigmaTOF", 5., "Nsigma cut on TOF only"};
  Configurable<float> nSigmaTOFCombined{"nSigmaTOFCombined", 5., "Nsigma cut on TOF combined with TPC"};

  double massPi = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();
  double massAlpha = TDatabasePDG::Instance()->GetParticle(o2::track::PID::Alpha)->Mass();

  void process(const soa::Join<o2::aod::Collisions, o2::aod::EvSels, aod::CentRun2V0Ms>::iterator& inputCollision,
                soa::Join<aod::Track, aod::TrackCov, aod::TrackExtra, aod::pidTPCPi, aod::pidTPCAl, aod::pidTOFAl> const& track,
                aod::V0Datas const& fullV0s)
  {
    // Performing the event selection
    if (!inputCollision.alias()[kINT7]) {
      return;
    }
    if (!inputCollision.sel7()) {
      return;
    } // <- why we need two int7 cuts here? (B.Lim)

    histos.fill(HIST("EventQA/VtxZBeforeSel"), inputCollision.posZ());
    if (abs(inputCollision.posZ()) > 10.f) {
      return;
    } // z- vertext cut
    histos.fill(HIST("EventQA/VtxZAfterSel"), inputCollision.posZ());

    // fill centrality histos
    histos.fill(HIST("EventQA/Centrality"), inputCollision.centRun2V0M());

    /// PID selectors
    // Pion
    TrackSelectorPID selectorPion(kPiPlus);
    selectorPion.setRangePtTPC(pidTPCMinpT, pidTPCMaxpT);
    selectorPion.setRangeNSigmaTPC(-nSigmaTPC, nSigmaTPC);
    // Alpha
    TrackSelectorPID selectorAlpha(selectorPion); // copy the selection criteria
    selectorAlpha.setPDG(o2::track::PID::Alpha);                 // overide only pdg code
    for (auto v0: fullV0s) {
        auto posTrack = v0.posTrack_as<aod::Track>();
        auto negTrack = v0.negTrack_as<aod::Track>();

        // Daughter DCA to PV cuts
        double trk1DCAr = posTrack.dcaXY();
        double trk2DCAr = negTrack.dcaXY();
        double trk1DCAz = posTrack.dcaZ();
        double trk2DCAz = negTrack.dcaZ();
        if ((trk1DCAz < cMinDCAzToPVcut1) || (trk2DCAz < cMinDCAzToPVcut2))
          continue;
        if ((trk1DCAz > cMaxDCAzToPVcut1) || (trk2DCAz > cMaxDCAzToPVcut2))
          continue;
        if ((trk1DCAr > cMaxDCArToPVcut1) || (trk2DCAr > cMaxDCArToPVcut2))
          continue;

        // Un-like sign pair only
        if (posTrack.sign() * negTrack.sign() > 0)
          continue;

        // Track PID cuts
        // TPC only PID
        int pidTrack1PionTPC = selectorPion.getStatusTrackPIDTPC(posTrack);
        int pidTrack1AlphaTPC = selectorAlpha.getStatusTrackPIDTPC(posTrack);
        int pidTrack2PionTPC = selectorPion.getStatusTrackPIDTPC(negTrack);
        int pidTrack2AlphaTPC = selectorAlpha.getStatusTrackPIDTPC(negTrack);

        // TPC+TOF PID
        int pidTrack1Pion = selectorPion.getStatusTrackPIDAll(posTrack);
        int pidTrack1Alpha = selectorAlpha.getStatusTrackPIDAll(posTrack);
        int pidTrack2Pion = selectorPion.getStatusTrackPIDAll(negTrack);
        int pidTrack2Alpha = selectorAlpha.getStatusTrackPIDAll(negTrack);

        auto isAnti = (pidTrack1AlphaTPC < 5) ? false : true;

        // Secondary Vertex selections
        if (v0.dcaV0daughters() > cMaxDCAHyper)
          continue;

        if (v0.v0cosPA() < cHyperCosPA)
          continue;

        if (v0.v0radius() < cMinHyperRadius)
          continue;
        
        auto arrMom = array{
            array{v0.pxpos(), v0.pypos(), v0.pzpos()},
            array{v0.pxneg(), v0.pyneg(), v0.pzneg()}
        };
        auto arrMass = array{massPi, massAlpha};
        auto hyperMass = RecoDecay::M2(arrMom, arrMass);
        
        if (!isAnti) {
            histos.fill(HIST("hyperHydrogen4invmass"), hyperMass);
            histos.fill(HIST("h3hyperHydrogen4invmass"), inputCollision.centRun2V0M(), v0.pt(), hyperMass);
        }
        if (isAnti) {
            histos.fill(HIST("antihyperHydrogen4invmass"), hyperMass);
            histos.fill(HIST("h3antihyperHydrogen4invmass"), inputCollision.centRun2V0M(), v0.pt(), hyperMass);
        }   
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<hyperHydrogen4>(cfgc)};
  return workflow;
}