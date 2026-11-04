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

/// \file hadronPhotonProducer.cxx
/// \author Julius Kinner
/// \brief hadron-photon angular correlation table producer
///
/// Table producer for hadron-photon angular correlation analysis

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/HadronPhotonCorrelation.h"
//
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TableHelper.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TPDGCode.h>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;

// correlation derived data ===================================================================================================================================================================

template <typename T>
concept IsCollisionMCD = requires(T& collision) { collision.mcCollision(); };
template <typename T>
concept IsMcCollision = requires(T& collision) { collision.getGeneratorId(); };

struct HadronPhotonProducer {
  struct : ProducesGroup {
    // reco
    Produces<aod::CorrCollisions> corrCollision;
    Produces<aod::CorrCollisionIds> corrCollisionId;
    Produces<aod::CorrCollisionMcIds> corrCollisionMcId;
    Produces<aod::Triggers> trigger;
    Produces<aod::TriggerTrackIds> triggerTrackId;
    Produces<aod::Hadrons> hadron;
    Produces<aod::HadronTrackIds> hadronTrackId;
    Produces<aod::Pipms> pipm;
    Produces<aod::PipmTrackIds> pipmTrackId;
    Produces<aod::PhotonPcms> photonPcm;
    Produces<aod::PhotonPcmTrackIds> photonPcmTrackId;
    Produces<aod::PhotonPcmV0PhotonIds> photonPcmV0PhotonId;
    Produces<aod::PhotonPcmPairs> photonPcmPair;
    Produces<aod::PhotonPcmPairTrackIds> photonPcmPairTrackId;
    Produces<aod::PhotonPcmPairV0PhotonIds> photonPcmPairV0PhotonId;
    // true
    Produces<aod::McCollisionsExtraCorr> mcCollisionExtraCorr;
    Produces<aod::TrueTriggers> trueTrigger;
  } tableProduction;

  Configurable<double> zPvMax{"zPvMax", 7, "maximum absZ primary-vertex cut"};
  Configurable<int> occupancyMin{"occupancyMin", 0, "minimum occupancy cut"};
  Configurable<int> occupancyMax{"occupancyMax", 2000, "maximum occupancy cut"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "JE framework - event selection"};
  Configurable<bool> skipMBGapEvents{"skipMBGapEvents", true, "skip MB gap events"};
  Configurable<bool> applyRctSelection{"applyRctSelection", true, "apply RCT selection"};
  Configurable<std::string> rctLabel{"rctLabel", "CBT_hadronPid", "RCT selection label"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "JE framework - track selections"};
  Configurable<std::string> triggerMasks{"triggerMasks", "", "JE framework - skimmed data trigger masks (relevent for correlation: fTrackLowPt,fTrackHighPt)"};

  Configurable<double> ptMinHadron{"ptMinHadron", 0.1, "min pt value for hadrons"};

  Configurable<double> piPidLowPt{"piPidLowPt", 0.5, "max pt value for pipm pid without tof"};
  Configurable<double> piPidHighPt{"piPidHighPt", 2.5, "min pt value for pipm pid without tof in relativistic rise of Bethe-Bloch"};
  Configurable<std::vector<double>> nSigmaPiTpcLowPt{"nSigmaPiTpcLowPt", {-2, 2}, "minimum-maximum nSigma for pipm in tpc at low pt"};
  Configurable<std::vector<double>> nSigmaPiTpcMidPt{"nSigmaPiTpcMidPt", {-1, 1}, "minimum-maximum nSigma for pipm in tpc at mid pt"};
  Configurable<std::vector<double>> nSigmaPiTof{"nSigmaPiTof", {-1, 2}, "minimum-maximum nSigma for pipm in tof"};
  Configurable<std::vector<double>> nSigmaPiRelRise{"nSigmaPiRelRise", {0, 2}, "minimum-maximum nSigma pipm tpc at high pt"};

  Configurable<std::vector<double>> pcmAbsDcaZLim{"pcmAbsDcaZLim", {0, 5}, "absolute values of dcaZ limits for pcm"};

  Configurable<bool> oneTriggerPerEvent{"oneTriggerPerEvent", false, "use only one trigger (of highest pt) in an event"};

  // configurables from other tasks

  double absEtaMax = -1;
  double ptTrigMin = -1;

  // derivatives of configurables

  std::vector<int> eventSelectionBits;
  int trackSelection = -1;
  std::vector<int> triggerMaskBits;

  struct TruthLevelInfo {
    TruthLevel truthLevel = TruthLevel::Unknown;
    double ptTruePure = -1;
  };
  struct PhotonImpurityInfo {
    PhotonImpurity impurity = PhotonImpurity::Unknown;
    double ptTrueOrigin = -1;
  };

  // map jetTracks to produced triggers
  std::map<int32_t, int32_t> trackToTriggerMap;

  HistogramRegistry histos{"histogramRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  // for mc
  Service<framework::O2DatabasePDG> pdg{};

  // partitions++
  SliceCache cache;

  Preslice<aod::JetTracks> perColTracks = aod::jtrack::collisionId;
  Preslice<aod::JetParticles> perColMcParticles = aod::jmcparticle::mcCollisionId;

  Preslice<aod::V0PhotonsKF> perColV0Photons = aod::v0photonkf::collisionId;

  // functions ================================================================================================================================================================================

  // get v0 leg tracks {posTrack, negTrack} as jet tracks via JTrackPI. throws in case of problems
  auto getV0LegJetTracks(auto const& v0Photon, auto const& tracks) {
    auto const posTrack = tracks.rawIteratorAt(v0Photon.posTrack().trackId() - tracks.offset());
    auto const negTrack = tracks.rawIteratorAt(v0Photon.negTrack().trackId() - tracks.offset());
    // check index match
    if (posTrack.trackId() != v0Photon.posTrack().trackId() || negTrack.trackId() != v0Photon.negTrack().trackId()) {
      throw std::runtime_error("v0Photon legs not correctly matched to jetTracks");
    }
    return std::tuple(posTrack, negTrack);
  }

  // selections ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // reco only --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // event selection
  bool checkEventSelection(auto const& collision)
  {
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return false;
    }
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, skipMBGapEvents, applyRctSelection, rctLabel)) {
      return false;
    }
    if (std::abs(collision.posZ()) > zPvMax) {
      return false;
    }
    if (collision.trackOccupancyInTimeRange() < occupancyMin || collision.trackOccupancyInTimeRange() > occupancyMax) {
      return false;
    }
    return true;
  }

  // checks track cuts (wow, such explanation)
  bool checkTrackCuts(auto const& track)
  {
    if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
      return false;
    }
    if (!jetderiveddatautilities::applyTrackKinematics(track, ptMinHadron, 1000, -1 * absEtaMax, absEtaMax)) {
      return false;
    }
    return true;
  }

  // checks pipm selection (just pid (no additional track cuts)
  bool checkPipmTPCTOF(auto const& track)
  {
    // too low for tof
    if (track.pt() < piPidLowPt) {
      return (track.tpcNSigmaPi() > nSigmaPiTpcLowPt.value[0] && track.tpcNSigmaPi() < nSigmaPiTpcLowPt.value[1]);
    }
    // Bethe-Bloch overlap (-> tpc + tof)
    if (track.pt() < piPidHighPt) {
      if (track.hasTOF()) { // has to stay inside pt-if due to return-layout of function
        if (track.tpcNSigmaPi() > nSigmaPiTpcMidPt.value[0] && track.tpcNSigmaPi() < nSigmaPiTpcMidPt.value[1] &&
            track.tofNSigmaPi() > nSigmaPiTof.value[0] && track.tofNSigmaPi() < nSigmaPiTof.value[1]) {
          return true;
        }
      }
      return false;
    }
    // Bethe-Bloch rel rise (too high for tof)
    return (track.tpcNSigmaPi() > nSigmaPiRelRise.value[0] && track.tpcNSigmaPi() < nSigmaPiRelRise.value[1]);
  }

  // checks pipm selection (just pid (no additional track cuts))
  bool checkPipmTPC(auto const& track)
  {
    // Bethe-Bloch rel rise
    if (track.pt() > piPidHighPt) {
      if (track.tpcNSigmaPi() > nSigmaPiRelRise.value[0] && track.tpcNSigmaPi() < nSigmaPiRelRise.value[1]) {
        return true;
      }
    }
    return false;
  }

  // reco MC ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // checks if mcParticle is charged
  bool isCharged(auto const& mcParticle)
  {
    auto const pdgParticle = pdg->GetParticle(mcParticle.pdgCode());
    return pdgParticle && pdgParticle->Charge() != 0;
  }
  // checks if fast decaying mcParticle is 'primary'
  bool isDecayPrimary(auto const& mcParticle)
  {
    // identify decaying primary
    return mcParticle.producedByGenerator();
  }

  // check if particle has mother in parent tree
  // returns optional pt of mother
  std::optional<double> checkForMother(auto const& mcParticle, int const pdgCode, bool const checkAntiParticle, bool const includePreHadronisation = true)
  {
    if (!mcParticle.has_mothers()) {
      return std::nullopt;
    }
    auto const mothers = mcParticle.template mothers_as<aod::JetParticles>();
    auto mother = mothers.begin();
    auto motherEnd = mothers.end();
    while (mother != motherEnd) {
      // optional exclude pre-hadonisation
      if (!includePreHadronisation && (*mother).producedByGenerator() && std::abs((*mother).getGenStatusCode()) < 80) {
        return std::nullopt;
      }
      if ((*mother).pdgCode() == pdgCode || (checkAntiParticle && (*mother).pdgCode() == -pdgCode)) {
        return (*mother).pt();
      }
      if (auto const ptMother = checkForMother(*mother, pdgCode, checkAntiParticle, includePreHadronisation)) {
        return ptMother;
      }
      mother++;
    }
    return std::nullopt;
  }

  // checks if tracks come from photon conversion
  template <typename T_track>
  bool isConversionPhoton(T_track const& posTrack, T_track const& negTrack)
  {
    auto const& posParticle = posTrack.mcParticle();
    auto const& negParticle = negTrack.mcParticle();
    if (!posParticle.has_mothers() || !negParticle.has_mothers()) {
      return false;
    }
    // check same mother
    auto const& posMothers = posParticle.template mothers_as<aod::JetParticles>();
    auto const& negMothers = negParticle.template mothers_as<aod::JetParticles>();
    if (posMothers.begin()->globalIndex() != negMothers.begin()->globalIndex()) {
      return false;
    }
    // check photon
    if (posMothers.begin()->pdgCode() != PDG_t::kGamma) {
      return false;
    }

    return true;
  };
  // checks if tracks come from double conversion
  template <typename T_track>
  bool isGGFromDoubleConversion(T_track const& posTrack1, T_track const& negTrack1, T_track const& posTrack2, T_track const& negTrack2, int const pdgCode)
  {
    // check conversion
    if (!isConversionPhoton(posTrack1, negTrack1) || !isConversionPhoton(posTrack2, negTrack2)) {
      return false;
    }
    auto const& photons1 = posTrack1.mcParticle().template mothers_as<aod::JetParticles>();
    auto const& photons2 = posTrack2.mcParticle().template mothers_as<aod::JetParticles>();
    if (!photons1.begin()->has_mothers() || !photons2.begin()->has_mothers()) {
      return false;
    }
    auto const& mothers1 = photons1.begin()->template mothers_as<aod::JetParticles>();
    auto const& mothers2 = photons2.begin()->template mothers_as<aod::JetParticles>();
    // check same mother of photons
    if (mothers1.begin()->globalIndex() != mothers2.begin()->globalIndex()) {
      return false;
    }
    // check particle type
    if (mothers1.begin()->pdgCode() != pdgCode) {
      return false;
    }

    return true;
  };

  // get highest reco truthLevel particle canidate

  template <AssocType AssocTypeV>
  TruthLevelInfo getMcRecoInfoHadronPipm(auto const& collision, auto const& track) {
    if (!track.has_mcParticle()) {
      return {TruthLevel::Reco, -1};
    }
    auto const& mcParticle = track.mcParticle();
    if (mcParticle.mcCollisionId() != collision.mcCollisionId()) {
      return {TruthLevel::RecoMcMatched, -1};
    }
    if constexpr (AssocTypeV == AssocType::Hadron) {
      if (!isCharged(mcParticle) || !mcParticle.isPhysicalPrimary()) {
        return {TruthLevel::RecoAssocEv, -1};
      }
    } else if constexpr (AssocTypeV == AssocType::Pipm) {
      if (std::abs(mcParticle.pdgCode()) != PDG_t::kPiPlus || !mcParticle.isPhysicalPrimary()) {
        return {TruthLevel::RecoAssocEv, -1};
      }
    } else {
      static_assert(AlwaysFalseValue<AssocTypeV>, "assocType not found");
    }
    double const ptTruePure = mcParticle.pt();
    return {TruthLevel::RecoPure, ptTruePure};
  }
  std::tuple<TruthLevelInfo, PhotonImpurityInfo> getMcRecoInfoPhoton(auto const& collision, auto const& v0Photon, auto const& tracks) {
    TruthLevelInfo truthLevelInfo;
    PhotonImpurityInfo impurityInfo;

    auto const [posTrack, negTrack] = getV0LegJetTracks(v0Photon, tracks);
    if (!posTrack.has_mcParticle() || !negTrack.has_mcParticle()) {
      truthLevelInfo.truthLevel = TruthLevel::Reco;
      return {truthLevelInfo, impurityInfo};
    }
    if (posTrack.mcParticle().mcCollisionId() != collision.mcCollisionId() || negTrack.mcParticle().mcCollisionId() != collision.mcCollisionId()) {
      truthLevelInfo.truthLevel = TruthLevel::RecoMcMatched;
      impurityInfo.impurity = PhotonImpurity::WrongEv;
      return {truthLevelInfo, impurityInfo};
    }
    auto const& photons = posTrack.mcParticle().template mothers_as<aod::JetParticles>();
    bool const isPhoton = isConversionPhoton(posTrack, negTrack);
    if (!isPhoton || !photons.begin()->isPhysicalPrimary()) {
      truthLevelInfo.truthLevel = TruthLevel::RecoAssocEv;
      if (!isPhoton) {
        impurityInfo.impurity = PhotonImpurity::MisId;
      } else {
        if (auto const ptMother = checkForMother(*(photons.begin()), PDG_t::kK0Short, true, false)) {
          impurityInfo.impurity = PhotonImpurity::K0Short;
          impurityInfo.ptTrueOrigin = *ptMother;
        } else if (auto const ptMother = checkForMother(*(photons.begin()), PDG_t::kK0Long, true, false)) {
          impurityInfo.impurity = PhotonImpurity::K0Long;
          impurityInfo.ptTrueOrigin = *ptMother;
        } else if (auto const ptMother = checkForMother(*(photons.begin()), PDG_t::kLambda0, true, false)) {
          impurityInfo.impurity = PhotonImpurity::Lambda0;
          impurityInfo.ptTrueOrigin = *ptMother;
        } else {
          impurityInfo.impurity = PhotonImpurity::Other;
        }
      }
      return {truthLevelInfo, impurityInfo};
    }
    truthLevelInfo.truthLevel = TruthLevel::RecoPure;
    truthLevelInfo.ptTruePure = photons.begin()->pt();
    impurityInfo.impurity = PhotonImpurity::Pure;
    return {truthLevelInfo, impurityInfo};
  }
  template <AssocType AssocTypeV>
  std::tuple<TruthLevelInfo, PhotonImpurityInfo> getMcRecoInfoDiPhoton(auto const& collision, auto const& v0Photon1, auto const& v0Photon2, auto const& tracks) {
    TruthLevelInfo truthLevelInfo;
    PhotonImpurityInfo impurityInfo;

    int constexpr PdgCode = []() {
      if constexpr (AssocTypeV == AssocType::Pi0Pcm) {
        return PDG_t::kPi0;
      } else if constexpr (AssocTypeV == AssocType::EtaPcm) {
        return constants::physics::Pdg::kEta;
      } else {
        static_assert(AlwaysFalseValue<AssocTypeV>, "assocType not found");
      }
    }();
    auto const [posTrack1, negTrack1] = getV0LegJetTracks(v0Photon1, tracks);
    auto const [posTrack2, negTrack2] = getV0LegJetTracks(v0Photon2, tracks);
    if (!posTrack1.has_mcParticle() || !negTrack1.has_mcParticle() || !posTrack2.has_mcParticle() || !negTrack2.has_mcParticle()) {
      truthLevelInfo.truthLevel = TruthLevel::Reco;
      return {truthLevelInfo, impurityInfo};
    }
    if (!isGGFromDoubleConversion(posTrack1, negTrack1, posTrack2, negTrack2, PdgCode)) {
      truthLevelInfo.truthLevel = TruthLevel::RecoMcMatched;
      return {truthLevelInfo, impurityInfo};
    }
    if (posTrack1.mcParticle().mcCollisionId() != collision.mcCollisionId()) {
      truthLevelInfo.truthLevel = TruthLevel::PseudoReco;
      if constexpr (AssocTypeV == AssocType::Pi0Pcm) {
        impurityInfo.impurity = PhotonImpurity::WrongEv;
      }
      return {truthLevelInfo, impurityInfo};
    }
    // note: do not dereference 'mothers' and store element in 'const&' (leads to seg fault)
    auto const& photons1 = posTrack1.mcParticle().template mothers_as<aod::JetParticles>();
    auto const& mothersOfPhoton = photons1.begin()->template mothers_as<aod::JetParticles>();
    if (!isDecayPrimary(*(mothersOfPhoton.begin()))) {
      truthLevelInfo.truthLevel = TruthLevel::PseudoRecoAssocEv;
      if constexpr (AssocTypeV == AssocType::Pi0Pcm) {
        if (auto const ptMother = checkForMother(*(photons1.begin()), PDG_t::kK0Short, true, false)) {
          impurityInfo.impurity = PhotonImpurity::K0Short;
          impurityInfo.ptTrueOrigin = *ptMother;
        } else if (auto const ptMother = checkForMother(*(photons1.begin()), PDG_t::kK0Long, true, false)) {
          impurityInfo.impurity = PhotonImpurity::K0Long;
          impurityInfo.ptTrueOrigin = *ptMother;
        } else if (auto const ptMother = checkForMother(*(photons1.begin()), PDG_t::kLambda0, true, false)) {
          impurityInfo.impurity = PhotonImpurity::Lambda0;
          impurityInfo.ptTrueOrigin = *ptMother;
        } else {
          impurityInfo.impurity = PhotonImpurity::Other;
        }
      }
      return {truthLevelInfo, impurityInfo};
    }
    truthLevelInfo.truthLevel = TruthLevel::PseudoRecoPure;
    truthLevelInfo.ptTruePure = mothersOfPhoton.begin()->pt();
    if constexpr (AssocTypeV == AssocType::Pi0Pcm) {
      impurityInfo.impurity = PhotonImpurity::Pure;
    }
    return {truthLevelInfo, impurityInfo};
  }

  // other ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // mc split event selection based in the thrid decimal of mc-true posZ
  int getThirdDecimalTruePosZ(auto const& collision)
  {
    double absPosZ;

    // get third decimal
    if constexpr (IsCollisionMCD<decltype(collision)>) {
      if (!collision.has_mcCollision()) {
        return -1;
      }
      absPosZ = std::abs(collision.template mcCollision_as<aod::JetMcCollisions>().posZ());
    } else if constexpr (IsMcCollision<decltype(collision)>) {
      absPosZ = std::abs(collision.posZ());
    } else {
      static_assert(AlwaysFalseType<decltype(collision)>, "collision type not allowed");
    }
    int const thirdDecimalPosZ = static_cast<int>(1000 * absPosZ) % 10;

    return thirdDecimalPosZ;
  }

  // analysis /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void init(InitContext& initContext)
  {
    // selection JE utilities
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(triggerMasks);

    // get variables from other tasks
    if (!common::core::getTaskOptionValue(initContext, "hadron-photon-correlation", "absEtaMax", absEtaMax, false)) {
      throw std::runtime_error("absEtaMax not set");
    };
    if (!common::core::getTaskOptionValue(initContext, "hadron-photon-correlation", "ptTrigMin", ptTrigMin, false)) {
      throw std::runtime_error("ptTrigMin not set");
    }

    // histograms
    AxisSpec const axisZPv{40, -10, 10, "#it{z}_{pv}"};
    AxisSpec const axisMult{101, -0.5, 100.5, "multiplicity"};
    AxisSpec const axisHadronicRate{1001, -0.5, 1000.5, "hadronicRate"};
    AxisSpec const axisPt{200, 0, 100, "#it{p}_{T}"};
    AxisSpec const axisPhi{72, 0, constants::math::TwoPI, "#varphi"};
    AxisSpec const axisEta{64, -absEtaMax, absEtaMax, "#eta"};
    AxisSpec const axisMgg{240, 0, 0.8, "#it{m}_{#gamma#gamma}"};
    // reco
    histos.add("reco/h2_zPvMult", "h2_zPvMult", kTH2D, {axisZPv, axisMult}, true);
    histos.add("reco/h1_hadronicRate", "h1_hadronicRate", kTH1D, {axisHadronicRate}, true);
    histos.add("reco/h3_ptPhiEta_trigger", "h3_ptPhiEta_trigger", kTH3D, {axisPt, axisPhi, axisEta}, true);
    histos.add("reco/h3_ptPhiEta_hadron", "h3_ptPhiEta_hadron", kTH3D, {axisPt, axisPhi, axisEta}, true);
    histos.add("reco/h3_ptPhiEta_pipm", "h3_ptPhiEta_pipm", kTH3D, {axisPt, axisPhi, axisEta}, true);
    histos.add("reco/h3_ptPhiEta_photonPcm", "h3_ptPhiEta_photonPcm", kTH3D, {axisPt, axisPhi, axisEta}, true);
    histos.add("reco/h3_ptPhiEta_photonPcmPair", "h3_ptPhiEta_photonPcmPair", kTH3D, {axisPt, axisPhi, axisEta}, true);
    // true
    histos.add("true/h2_zPvMult", "h2_zPvMult", kTH2D, {axisZPv, axisMult}, true);
    histos.add("true/h3_ptPhiEta_trigger", "h3_ptPhiEta_trigger", kTH3D, {axisPt, axisPhi, axisEta}, true);
  }

  void fillProcessRecoCollisionTrigger(auto const& collisions, auto const& tracks)
  {
    // clear mapping at start of dataframe (in first process)
    trackToTriggerMap.clear();

    for (auto const& collision : collisions) {
      // event selection
      bool isSelectedEvent = checkEventSelection(collision);

      // group collision
      auto const tracksThisEvent = tracks.sliceBy(perColTracks, collision.globalIndex());

      // mc split event selection
      // default -1 for data reconstruction
      int thirdDecimalTruePosZ = -1;
      if constexpr (IsCollisionMCD<decltype(collision)>) {
        thirdDecimalTruePosZ = getThirdDecimalTruePosZ(collision);
      }
      // maximum track pt
      float maxTrackPt = -1;
      int maxTriggerCollisionId = -1;
      int maxTriggerGlobalIndex = -1;
      float maxTriggerPhi = 0;
      float maxTriggerEta = 0;
      // number global tracks
      int nGlobalTracks = 0;

      // count global tracks (for independence of multiplicity task (uses only JE derieved data))
      for (auto const& track : tracksThisEvent) {
        // track selection
        if (!checkTrackCuts(track)) {
          continue;
        }

        nGlobalTracks++;
        if (track.pt() > maxTrackPt) {
          maxTrackPt = track.pt();
          maxTriggerCollisionId = track.collisionId();
          maxTriggerGlobalIndex = track.globalIndex();
          maxTriggerPhi = track.phi();
          maxTriggerEta = track.eta();
        }

        if (!isSelectedEvent) {
          continue;
        }
        if (track.pt() < ptTrigMin) {
          continue;
        }
        if (oneTriggerPerEvent) {
          continue;
        }

        // trigger info
        histos.fill(HIST("reco/h3_ptPhiEta_trigger"), track.pt(), track.phi(), track.eta());
        tableProduction.trigger(track.collisionId(),
                                track.pt(), track.phi(), track.eta());
        tableProduction.triggerTrackId(track.globalIndex());
        trackToTriggerMap[track.globalIndex()] = tableProduction.trigger.lastIndex();
      }
      bool const isTriggerEvent = maxTrackPt > ptTrigMin;
      if (oneTriggerPerEvent && isSelectedEvent && isTriggerEvent) {
        histos.fill(HIST("reco/h3_ptPhiEta_trigger"), maxTrackPt, maxTriggerPhi, maxTriggerEta);
        tableProduction.trigger(maxTriggerCollisionId,
                                maxTrackPt, maxTriggerPhi, maxTriggerEta);
        tableProduction.triggerTrackId(maxTriggerGlobalIndex);
        trackToTriggerMap[maxTriggerGlobalIndex] = tableProduction.trigger.lastIndex();
      }

      // collision info
      histos.fill(HIST("reco/h2_zPvMult"), collision.posZ(), nGlobalTracks);
      histos.fill(HIST("reco/h1_hadronicRate"), collision.hadronicRate());
      tableProduction.corrCollision(isSelectedEvent, thirdDecimalTruePosZ, isTriggerEvent, maxTrackPt, collision.posZ(), nGlobalTracks);
      tableProduction.corrCollisionId(collision.globalIndex());
      // mc true collsion label
      if constexpr (IsCollisionMCD<decltype(collision)>) {
        tableProduction.corrCollisionMcId(collision.mcCollisionId());
      }
    }
  }

  void processRecoCollisionTrigger(soa::Join<aod::JetCollisions, aod::JCollisionPIs> const& collisions, aod::JetTracks const& tracks)
  {
    fillProcessRecoCollisionTrigger(collisions, tracks);
  }
  PROCESS_SWITCH(HadronPhotonProducer, processRecoCollisionTrigger, "process correlation collision_extra and trigger table (reconstructed)", false);

  void processMcRecoCollisionTrigger(soa::Join<aod::JetCollisionsMCD, aod::JCollisionPIs> const& collisions, aod::JetTracks const& tracks, aod::JetMcCollisions const&)
  {
    fillProcessRecoCollisionTrigger(collisions, tracks);
  }
  PROCESS_SWITCH(HadronPhotonProducer, processMcRecoCollisionTrigger, "process correlation collision_extra and trigger table (reconstructed MC)", false);

  template <typename TrackPidTable>
  void fillProcessRecoHadronPipm(auto const& collisions, auto const& tracks)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!checkEventSelection(collision)) {
        continue;
      }

      // group collision
      auto const tracksThisEvent = tracks.sliceBy(perColTracks, collision.globalIndex());

      // hadron/pipm
      for (auto const& track : tracksThisEvent) {
        std::vector<int> intersectTriggerIds;

        // track selection
        if (!checkTrackCuts(track)) {
          continue;
        }

        // get possible overlap with triggers
        auto const triggerTrackIt = trackToTriggerMap.find(track.globalIndex());
        if (triggerTrackIt != trackToTriggerMap.end()) {
          intersectTriggerIds.push_back(triggerTrackIt->second);
        }

        // mc reco info
        TruthLevelInfo const truthLevelInfoHadron = [&] {
          if constexpr (IsCollisionMCD<decltype(collision)>) {
            return getMcRecoInfoHadronPipm<AssocType::Hadron>(collision, track);
          } else {
            return TruthLevelInfo{};
          }
        }();

        // hadron info
        histos.fill(HIST("reco/h3_ptPhiEta_hadron"), track.pt(), track.phi(), track.eta());
        tableProduction.hadron(track.collisionId(), intersectTriggerIds,
                              track.pt(), track.phi(), track.eta(),
                              static_cast<int>(truthLevelInfoHadron.truthLevel), truthLevelInfoHadron.ptTruePure);
        tableProduction.hadronTrackId(track.globalIndex());

        // pipm selection
        auto const& trackPid = track.template track_as<TrackPidTable>();
        if constexpr (requires { trackPid.tofNSigmaPi(); }) {
          if (!checkPipmTPCTOF(trackPid)) {
            continue;
          }
        } else {
          if (!checkPipmTPC(trackPid)) {
            continue;
          }
        }

        // mc reco info
        TruthLevelInfo const truthLevelInfoPipm = [&] {
          if constexpr (IsCollisionMCD<decltype(collision)>) {
            return getMcRecoInfoHadronPipm<AssocType::Pipm>(collision, track);
          } else {
            return TruthLevelInfo{};
          }
        }();

        // pipm info
        histos.fill(HIST("reco/h3_ptPhiEta_pipm"), track.pt(), track.phi(), track.eta());
        tableProduction.pipm(track.collisionId(), intersectTriggerIds,
                            track.pt(), track.phi(), track.eta(),
                            static_cast<int>(truthLevelInfoPipm.truthLevel), truthLevelInfoPipm.ptTruePure);
        tableProduction.pipmTrackId(track.globalIndex());
      }
    }
  }

  void processRecoHadronPipmTPCTOF(aod::JetCollisions const& collisions,
                                   soa::Join<aod::JetTracks, aod::JTrackPIs> const& tracks, soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTOFPi> const&)
  {
    fillProcessRecoHadronPipm<soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTOFPi>>(collisions, tracks);
  }
  PROCESS_SWITCH(HadronPhotonProducer, processRecoHadronPipmTPCTOF, "process hadron and pipm (TPC-TOF) table (reconstructed)", false);

  void processMcRecoHadronPipmTPCTOF(aod::JetCollisionsMCD const& collisions,
                                     soa::Join<aod::JetTracksMCD, aod::JTrackPIs> const& tracks, soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTOFPi> const&,
                                     aod::JetParticles const&)
  {
    fillProcessRecoHadronPipm<soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTOFPi>>(collisions, tracks);
  }
  PROCESS_SWITCH(HadronPhotonProducer, processMcRecoHadronPipmTPCTOF, "process hadron and pipm (TPC-TOF) table (reconstructed MC)", false);

  void processRecoHadronPipmTPC(aod::JetCollisions const& collisions,
                                soa::Join<aod::JetTracks, aod::JTrackPIs> const& tracks, soa::Join<aod::Tracks, aod::pidTPCPi> const&)
  {
    fillProcessRecoHadronPipm<soa::Join<aod::Tracks, aod::pidTPCPi>>(collisions, tracks);
  }
  PROCESS_SWITCH(HadronPhotonProducer, processRecoHadronPipmTPC, "process hadron and pipm (TPC) table (reconstructed)", false);

  void processMcRecoHadronPipmTPC(aod::JetCollisionsMCD const& collisions,
                                  soa::Join<aod::JetTracksMCD, aod::JTrackPIs> const& tracks, soa::Join<aod::Tracks, aod::pidTPCPi> const&,
                                  aod::JetParticles const&)
  {
    fillProcessRecoHadronPipm<soa::Join<aod::Tracks, aod::pidTPCPi>>(collisions, tracks);
  }
  PROCESS_SWITCH(HadronPhotonProducer, processMcRecoHadronPipmTPC, "process hadron and pipm (TPC) table (reconstructed MC)", false);

  void fillProcessRecoPhotonPcm(auto const& collisions, auto const& v0Photons, auto const& tracks)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!checkEventSelection(collision)) {
        continue;
      }

      // photonsPcm (for some reason collisionId not an index column (?))
      auto const v0PhotonsThisEvent = v0Photons.sliceBy(perColV0Photons, collision.collisionId());

      // photonPcm
      for (auto const& v0Photon : v0PhotonsThisEvent) {
        std::vector<int> intersectTriggerIds;

        // photon selection
        if (std::abs(v0Photon.eta()) > absEtaMax) {
          continue;
        }
        if (std::abs(v0Photon.dcaZtopv()) < pcmAbsDcaZLim.value[0] || std::abs(v0Photon.dcaZtopv()) > pcmAbsDcaZLim.value[1]) {
          continue;
        }

        // get legs and check v0Leg-jetTrack match
        auto const [posTrack, negTrack] = getV0LegJetTracks(v0Photon, tracks);

        // get possible overlap with triggers
        auto const triggerPosTrackIt = trackToTriggerMap.find(posTrack.globalIndex());
        if (triggerPosTrackIt != trackToTriggerMap.end()) {
          intersectTriggerIds.push_back(triggerPosTrackIt->second);
        }
        auto const triggerNegTrackIt = trackToTriggerMap.find(negTrack.globalIndex());
        if (triggerNegTrackIt != trackToTriggerMap.end()) {
          intersectTriggerIds.push_back(triggerNegTrackIt->second);
        }

        // mc reco info
        auto const [truthLevelInfoPhotonPcm, impurityInfoPhotonPcm] = [&] {
          if constexpr (IsCollisionMCD<decltype(collision)>) {
            return getMcRecoInfoPhoton(collision, v0Photon, tracks);
          } else {
            return std::tuple{TruthLevelInfo{}, PhotonImpurityInfo{}};
          }
        }();

        // photonPcm info
        histos.fill(HIST("reco/h3_ptPhiEta_photonPcm"), v0Photon.pt(), v0Photon.phi(), v0Photon.eta());
        tableProduction.photonPcm(v0Photon.collisionId(), intersectTriggerIds,
                                  v0Photon.pt(), v0Photon.phi(), v0Photon.eta(),
                                  static_cast<int>(truthLevelInfoPhotonPcm.truthLevel), truthLevelInfoPhotonPcm.ptTruePure,
                                  static_cast<int>(impurityInfoPhotonPcm.impurity), impurityInfoPhotonPcm.ptTrueOrigin);
        tableProduction.photonPcmTrackId(v0Photon.posTrack().trackId(), v0Photon.negTrack().trackId());
        tableProduction.photonPcmV0PhotonId(v0Photon.globalIndex());
      }

      // photonPcmPair
      for (auto const& [v0Photon1, v0Photon2] : soa::combinations(soa::CombinationsStrictlyUpperIndexPolicy(v0PhotonsThisEvent, v0PhotonsThisEvent))) {
        std::vector<int> intersectTriggerIds;

        // get kinematics
        ROOT::Math::PtEtaPhiMVector const p4V0Pcm1(v0Photon1.pt(), v0Photon1.eta(), v0Photon1.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector const p4V0Pcm2(v0Photon2.pt(), v0Photon2.eta(), v0Photon2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector const p4V0PcmPair = p4V0Pcm1 + p4V0Pcm2;

        // pair selection
        if (std::abs(p4V0PcmPair.Eta()) > absEtaMax) {
          continue;
        }

        // get legs and check v0Leg-jetTrack match
        auto const [posTrack1, negTrack1] = getV0LegJetTracks(v0Photon1, tracks);
        auto const [posTrack2, negTrack2] = getV0LegJetTracks(v0Photon2, tracks);

        // get possible overlap with triggers
        auto const triggerPosTrack1It = trackToTriggerMap.find(posTrack1.globalIndex());
        if (triggerPosTrack1It != trackToTriggerMap.end()) {
          intersectTriggerIds.push_back(triggerPosTrack1It->second);
        }
        auto const triggerNegTrack1It = trackToTriggerMap.find(negTrack1.globalIndex());
        if (triggerNegTrack1It != trackToTriggerMap.end()) {
          intersectTriggerIds.push_back(triggerNegTrack1It->second);
        }
        auto const triggerPosTrack2It = trackToTriggerMap.find(posTrack2.globalIndex());
        if (triggerPosTrack2It != trackToTriggerMap.end()) {
          intersectTriggerIds.push_back(triggerPosTrack2It->second);
        }
        auto const triggerNegTrack2It = trackToTriggerMap.find(negTrack2.globalIndex());
        if (triggerNegTrack2It != trackToTriggerMap.end()) {
          intersectTriggerIds.push_back(triggerNegTrack2It->second);
        }

        // mc reco info
        auto const [truthLevelInfoPi0Pcm, photonImpurityInfoPi0Pcm] = [&] {
          if constexpr (IsCollisionMCD<decltype(collision)>) {
            return getMcRecoInfoDiPhoton<AssocType::Pi0Pcm>(collision, v0Photon1, v0Photon2, tracks);
          } else {
            return std::tuple{TruthLevelInfo{}, PhotonImpurityInfo{}};
          }
        }();
        TruthLevelInfo const truthLevelInfoEtaPcm = [&] {
          if constexpr (IsCollisionMCD<decltype(collision)>) {
            return std::get<0>(getMcRecoInfoDiPhoton<AssocType::EtaPcm>(collision, v0Photon1, v0Photon2, tracks));
          } else {
            return TruthLevelInfo{};
          }
        }();
        double const ptTruePurePhotonPcmPair = std::max(truthLevelInfoPi0Pcm.ptTruePure, truthLevelInfoEtaPcm.ptTruePure);

        // photonPcmPair info
        histos.fill(HIST("reco/h3_ptPhiEta_photonPcmPair"), p4V0PcmPair.Pt(), RecoDecay::constrainAngle(p4V0PcmPair.Phi(), 0), p4V0PcmPair.Eta());
        tableProduction.photonPcmPair(v0Photon1.collisionId(), intersectTriggerIds,
                                      p4V0PcmPair.Pt(), RecoDecay::constrainAngle(p4V0PcmPair.Phi(), 0), p4V0PcmPair.Eta(), p4V0PcmPair.M(),
                                      static_cast<int>(truthLevelInfoPi0Pcm.truthLevel), static_cast<int>(truthLevelInfoEtaPcm.truthLevel), ptTruePurePhotonPcmPair,
                                      static_cast<int>(photonImpurityInfoPi0Pcm.impurity), photonImpurityInfoPi0Pcm.ptTrueOrigin);
        tableProduction.photonPcmPairTrackId(v0Photon1.posTrack().trackId(), v0Photon1.negTrack().trackId(), v0Photon2.posTrack().trackId(), v0Photon2.negTrack().trackId());
        tableProduction.photonPcmPairV0PhotonId(v0Photon1.globalIndex(), v0Photon2.globalIndex());
      }
    }
  }

  void processRecoPhotonPcm(soa::Join<aod::JetCollisions, aod::JCollisionPIs> const& collisions,
                            aod::V0PhotonsKF const& v0Photons, aod::V0Legs const&, soa::Join<aod::JetTracks, aod::JTrackPIs> const& tracks)
  {
    fillProcessRecoPhotonPcm(collisions, v0Photons, tracks);
  }
  PROCESS_SWITCH(HadronPhotonProducer, processRecoPhotonPcm, "process photonPcm table (reconstructed)", false);

  void processMcRecoPhotonPcm(soa::Join<aod::JetCollisionsMCD, aod::JCollisionPIs> const& collisions,
                              aod::V0PhotonsKF const& v0Photons, aod::V0Legs const&, soa::Join<aod::JetTracksMCD, aod::JTrackPIs> const& tracks,
                              aod::JetParticles const&)
  {
    fillProcessRecoPhotonPcm(collisions, v0Photons, tracks);
  }
  PROCESS_SWITCH(HadronPhotonProducer, processMcRecoPhotonPcm, "process photonPcm table (reconstructed MC)", false);

  void processMcCorrTables(aod::JetMcCollisions const& mcCollisions, aod::JetParticles const& mcParticles)
  {
    for (auto const& mcCollision : mcCollisions) {
      // group collision
      auto const mcParticlesThisEvent = mcParticles.sliceBy(perColMcParticles, mcCollision.globalIndex());

      // mc split event selection
      int const thirdDecimalTruePosZ = getThirdDecimalTruePosZ(mcCollision);
      // maximum charged particle pt
      float maxChargedPt = -1;
      int maxTriggerCollisionId = -1;
      int maxTriggerGlobalIndex = -1;
      float maxTriggerPhi = 0;
      float maxTriggerEta = 0;
      // number charged particles in eta range
      int nCharged = 0;

      // particle loop
      for (auto const& mcParticle : mcParticlesThisEvent) {
        // track selection
        auto const pdgParticle = pdg->GetParticle(mcParticle.pdgCode());
        if (!pdgParticle || pdgParticle->Charge() == 0) {
          continue;
        }
        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }
        if (std::abs(mcParticle.eta()) > absEtaMax || mcParticle.pt() < ptMinHadron) {
          continue;
        }

        nCharged++;
        if (mcParticle.pt() > maxChargedPt) {
          maxChargedPt = mcParticle.pt();
          maxTriggerCollisionId = mcParticle.mcCollisionId();
          maxTriggerGlobalIndex = mcParticle.globalIndex();
          maxTriggerPhi = mcParticle.phi();
          maxTriggerEta = mcParticle.eta();
        }

        // trigger selection
        if (mcParticle.pt() < ptTrigMin) {
          continue;
        }
        if (oneTriggerPerEvent) {
          continue;
        }

        // trigger info
        histos.fill(HIST("true/h3_ptPhiEta_trigger"), mcParticle.pt(), mcParticle.phi(), mcParticle.eta());
        tableProduction.trueTrigger(mcParticle.mcCollisionId(), mcParticle.globalIndex(),
                                    mcParticle.pt(), mcParticle.phi(), mcParticle.eta());
      }
      bool const isTriggerEvent = maxChargedPt > ptTrigMin;
      if (oneTriggerPerEvent && isTriggerEvent) {
        histos.fill(HIST("true/h3_ptPhiEta_trigger"), maxChargedPt, maxTriggerPhi, maxTriggerEta);
        tableProduction.trueTrigger(maxTriggerCollisionId, maxTriggerGlobalIndex,
                                    maxChargedPt, maxTriggerPhi, maxTriggerEta);
      }

      // collision info
      histos.fill(HIST("true/h2_zPvMult"), mcCollision.posZ(), nCharged);
      tableProduction.mcCollisionExtraCorr(thirdDecimalTruePosZ, isTriggerEvent, maxChargedPt, nCharged);
    }
  }
  PROCESS_SWITCH(HadronPhotonProducer, processMcCorrTables, "process table production (mc)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& configContext)
{
  return WorkflowSpec{adaptAnalysisTask<HadronPhotonProducer>(configContext)};
}
