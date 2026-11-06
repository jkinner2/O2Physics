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

/// \file hadronPhotonCorrelation.cxx
/// \author Julius Kinner
/// \brief hadron-photon angular correlation analysis
///
/// Analysis for angular correlations between jets and photons via two-particle correlations with charged high-pt triggers
/// Associated hadrons (tracks), pipm, photons (pcm), pi0 (pcm)
/// Also contains checks and monte-carlo (efficiency, purity, mc-true correlation,...)
/// End goal of studying correlations between direct photons and jets

#include "PWGJE/DataModel/HadronPhotonCorrelation.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
//
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TableHelper.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TMCProcess.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <format>
#include <functional>
#include <initializer_list>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

using CorrCollisions = aod::CorrCollisions;
using CorrCollisionsMcD = soa::Join<aod::CorrCollisions, aod::CorrCollisionMcIds>;
using CorrMcCollisions = soa::Join<aod::JetMcCollisions, aod::McCollisionsExtraCorr>;

using PhotonPcmsWithV0 = soa::Join<aod::PhotonPcms, aod::PhotonPcmV0PhotonIds>;
using PhotonPcmPairsWithV0 = soa::Join<aod::PhotonPcmPairs, aod::PhotonPcmPairV0PhotonIds>;

using BinningZPvMult = ColumnBinningPolicy<aod::corr_collision::PosZ, aod::corr_collision::NGlobalTracks>;

// correlation analysis ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// =/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/

const double absEtaMaxDefault = 0.8;

struct HadronPhotonCorrelation {
  // configurables

  // main config
  struct : o2::framework::ConfigurableGroup {
    // general (kenobi)
    Configurable<std::string> pathCcdbAnalysis{"pathCcdbAnalysis", "Users/j/jkinner/analysis/set_in_config", "base path to the ccdb for analysis objects"};
    Configurable<std::string> urlCcdb{"urlCcdb", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<int64_t> noLaterThanCcdbConfig{"noLaterThanCcdbConfig", -1, "latest acceptable timestamp of creation for the object (-1 for task start time)"};

    // shared across analysis tasks
    Configurable<double> absEtaMax{"absEtaMax", 0.8, "maximum absEta cut"};
    Configurable<double> ptTrigMin{"ptTrigMin", 5, "minimum pT of triggers"};

    // analysis settings
    Configurable<int> splitMcEvN{"splitMcEvN", 1, "number of event fractions for selection in MC (max 10) | (select 1 for data)"};
    Configurable<int> splitMcEvSelect{"splitMcEvSelect", 1, "selected event fraction number in MC (starts at 1) | (select 1 for data)"};

    // mc scaling to data
    // caution: applied to true-level very broadly, full mc-true possibly inconsistent (e.g. gammaDir scaled as well), only use for correcting
    Configurable<bool> applyScalingPhotonPt{"applyScalingPhotonPt", false, "apply photon mc scaling factor"};
    Configurable<std::vector<double>> ptRangeScalingPhoton{"ptRangeScalingPhoton", {0.5, 15}, "pt for photon mc scaling factor using tcm, switching to constant outside"};
    Configurable<std::vector<double>> parScalingPhotonPtData{"parScalingPhotonPtData", {}, "tcm params, corrected data, photonPcm"};
    Configurable<std::vector<double>> parScalingPhotonPtMc{"parScalingPhotonPtMc", {}, "tcm params, mc gen, photonPcm"};
    Configurable<bool> applyScalingPi0Pt{"applyScalingPi0Pt", false, "apply pi0 mc scaling factor"};
    Configurable<std::vector<double>> ptRangeScalingPi0{"ptRangeScalingPi0", {0.5, 10}, "pt for pi0 mc scaling factor using tcm, switching to constant outside"};
    Configurable<std::vector<double>> parScalingPi0PtData{"parScalingPi0PtData", {}, "tcm params, corrected data, pi0Pcm"};
    Configurable<std::vector<double>> parScalingPi0PtMc{"parScalingPi0PtMc", {}, "tcm params, mc gen, pi0Pcm"};
    Configurable<bool> applyScalingEtaPt{"applyScalingEtaPt", false, "apply eta mc scaling factor"};
    Configurable<std::vector<double>> ptRangeScalingEta{"ptRangeScalingEta", {1.0, 10}, "pt for eta mc scaling factor using tcm, switching to constant outside"};
    Configurable<std::vector<double>> parScalingEtaPtData{"parScalingEtaPtData", {}, "tcm params, corrected data, etaPcm"};
    Configurable<std::vector<double>> parScalingEtaPtMc{"parScalingEtaPtMc", {}, "tcm params, mc gen, etaPcm"};
    // scalings and corrections applied to reco
    Configurable<bool> applyScalingK0ShortPt{"applyScalingK0ShortPt", false, "apply K0Short mc scaling factor"};
    Configurable<bool> applyTriggerPtCorrection{"applyTriggerPtCorrection", false, "whether to apply on-the-fly total correction factor for triggers"};

    Configurable<bool> doTrigEvMixing{"doTrigEvMixing", false, "whether to use trigger events for trigger mixing"};

    Configurable<int> nTriggerSavedForMixing{"nTriggerSavedForMixing", 8192, "number of triggers that are saved for mixing with other events"};
    Configurable<int> nTriggerBinMinThreshold{"nTriggerBinMinThreshold", 4, "threshold minimum number of triggers in mixing bin required"};
    Configurable<double> iTriggerStartFraction{"iTriggerStartFraction", 0.5, "ratio of bin size for randomised iTrigger start"};
    Configurable<double> nTriggerScaleCoefficient{"nTriggerScaleCoefficient", 1, "coeffcient a in mixing-number pt power scaling"};
    Configurable<double> nTriggerScaleExponent{"nTriggerScaleExponent", 2, "exponent a in mixing-number pt power scaling"};
    Configurable<int> nMixingAt0McTrue{"nMixingAt0McTrue", 1, "number of triggers that are used for mc true mixing"};
    Configurable<int> nMixingAt0Hadron{"nMixingAt0Hadron", 16, "number of triggers that are used for hadron mixing"};
    Configurable<int> nMixingAt0Pipm{"nMixingAt0Pipm", 16, "number of triggers that are used for pipm mixing"};
    Configurable<int> nMixingAt0PhotonPcm{"nMixingAt0PhotonPcm", 256, "number of triggers that are saved for photonPcm mixing"};
    Configurable<int> nMixingAt0PhotonPcmPair{"nMixingAt0PhotonPcmPair", 1024, "number of triggers that are saved for photonPcmPair (pi0, eta) mixing"};

    Configurable<int> nNeighboursMixingPhotonPcmPair{"nNeighboursMixingPhotonPcmPair", 32, "number neighbours used for for photonPcm pair mixing"};

    Configurable<std::vector<double>> pi0PcmPeakMassRange{"pi0PcmPeakMassRange", {0.10, 0.15}, "photon-pair mass integration range for pi0Pcm"};
    Configurable<std::vector<double>> pi0PcmSideMassRange{"pi0PcmSideMassRange", {0.16, 0.22}, "photon-pair mass integration range outside pi0Pcm region"};
    Configurable<std::vector<double>> etaPcmPeakMassRange{"etaPcmPeakMassRange", {0.51, 0.56}, "photon-pair mass integration range for etaPcm"};
    Configurable<std::vector<double>> etaPcmLowSideMassRange{"etaPcmLowSideMassRange", {0.44, 0.48}, "photon-pair mass integration range below etaPcm region"};
    Configurable<std::vector<double>> etaPcmHighSideMassRange{"etaPcmHighSideMassRange", {0.58, 0.62}, "photon-pair mass integration range above etaPcm region"};

    Configurable<double> energyAsymmPhotonPairMax{"energyAsymmPhotonPairMax", 1.0, "maximum energy asymmetry ((E1 - E2) / (E1 + E2)) of photon pairs"};
  } config;

  // pcm cuts
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "pcmcut_group";
    o2::framework::Configurable<bool> cfg_require_v0_with_itstpc{"cfg_require_v0_with_itstpc", false, "flag to select V0s with ITS-TPC matched tracks"};
    o2::framework::Configurable<bool> cfg_require_v0_with_itsonly{"cfg_require_v0_with_itsonly", false, "flag to select V0s with ITSonly tracks"};
    o2::framework::Configurable<bool> cfg_require_v0_with_tpconly{"cfg_require_v0_with_tpconly", false, "flag to select V0s with TPConly tracks"};
    o2::framework::Configurable<float> cfg_min_pt_v0{"cfg_min_pt_v0", 0.1, "min pT for v0 photons at PV"};
    o2::framework::Configurable<float> cfg_max_pt_v0{"cfg_max_pt_v0", 1e+10, "max pT for v0 photons at PV"};
    o2::framework::Configurable<float> cfg_min_eta_v0{"cfg_min_eta_v0", -0.8, "min eta for v0 photons at PV"};
    o2::framework::Configurable<float> cfg_max_eta_v0{"cfg_max_eta_v0", 0.8, "max eta for v0 photons at PV"};
    o2::framework::Configurable<float> cfg_min_v0radius{"cfg_min_v0radius", 4.0, "min v0 radius"};
    o2::framework::Configurable<float> cfg_max_v0radius{"cfg_max_v0radius", 90.0, "max v0 radius"};
    o2::framework::Configurable<float> cfg_midL_v0radius{"cfg_midL_v0radius", -1.0, "middle low v0 radius for rejection if >0"};
    o2::framework::Configurable<float> cfg_midH_v0radius{"cfg_midH_v0radius", -1.0, "middle high v0 radius for rejection if >0"};
    o2::framework::Configurable<float> cfg_max_alpha_ap{"cfg_max_alpha_ap", 0.95, "max alpha for AP cut"};
    o2::framework::Configurable<float> cfg_max_qt_ap{"cfg_max_qt_ap", 0.01, "max qT for AP cut"};
    o2::framework::Configurable<float> cfg_min_cospa{"cfg_min_cospa", 0.999, "min V0 CosPA"};
    o2::framework::Configurable<float> cfg_max_pca{"cfg_max_pca", 1.5, "max distance btween 2 legs"};
    o2::framework::Configurable<float> cfg_max_chi2kf{"cfg_max_chi2kf", 1e+10, "max chi2/ndf with KF"};
    o2::framework::Configurable<bool> cfg_reject_v0_on_itsib{"cfg_reject_v0_on_itsib", true, "flag to reject V0s on ITSib"};
    o2::framework::Configurable<bool> cfg_apply_cuts_from_prefilter_derived{"cfg_apply_cuts_from_prefilter_derived", false, "flag to apply prefilter to V0"};

    o2::framework::Configurable<int> cfg_min_ncluster_tpc{"cfg_min_ncluster_tpc", 0, "min ncluster tpc"};
    o2::framework::Configurable<int> cfg_min_ncrossedrows{"cfg_min_ncrossedrows", 40, "min ncrossed rows"};
    o2::framework::Configurable<float> cfg_max_frac_shared_clusters_tpc{"cfg_max_frac_shared_clusters_tpc", 999.f, "max fraction of shared clusters in TPC"};
    o2::framework::Configurable<float> cfg_max_chi2tpc{"cfg_max_chi2tpc", 4.0, "max chi2/NclsTPC"};
    o2::framework::Configurable<float> cfg_max_chi2its{"cfg_max_chi2its", 36.0, "max chi2/NclsITS"};
    o2::framework::Configurable<float> cfg_min_TPCNsigmaEl{"cfg_min_TPCNsigmaEl", -3.0, "min. TPC n sigma for electron"};
    o2::framework::Configurable<float> cfg_max_TPCNsigmaEl{"cfg_max_TPCNsigmaEl", +3.0, "max. TPC n sigma for electron"};
    o2::framework::Configurable<bool> cfg_disable_itsonly_track{"cfg_disable_itsonly_track", false, "flag to disable ITSonly tracks"};
    o2::framework::Configurable<bool> cfg_disable_tpconly_track{"cfg_disable_tpconly_track", false, "flag to disable TPConly tracks"};
  } configPcmCuts;

  // histogram config
  struct : ConfigurableGroup {
    Configurable<int> nBinsZPvPlain{"nBinsZPvPlain", 28, "number zPv bins in histos"};
    Configurable<int> nBinsMultPlain{"nBinsMultPlain", 64, "number multiplicity bins in histos"};

    Configurable<int> nBinsPhi{"nBinsPhi", 72, "number phi bins"};
    Configurable<int> nBinsEta{"nBinsEta", 32, "number eta bins"};
    Configurable<int> nBinsDcaZ{"nBinsDcaZ", 100, "number dcaZ bins"};
    Configurable<int> nBinsMgg{"nBinsMgg", 240, "number mass-photon-pair bins"};

    Configurable<std::vector<double>> binsPtTrig{"binsPtTrig", {5, 10, 25, 100}, "correlation ptTrig bins"};
    Configurable<std::vector<double>> binsPtAssoc{"binsPtAssoc",
                                                  {0.2, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10,
                                                  12.5, 15, 17.5, 20, 25, 35, 50, 100},
                                                  "correlation ptAssoc bins"};
    Configurable<std::vector<double>> binsDPhi{"binsDPhi",
                                               {(0.00 - 1. / 4) * constants::math::TwoPI,
                                                (0.05 - 1. / 4) * constants::math::TwoPI, (0.10 - 1. / 4) * constants::math::TwoPI,
                                                (0.14 - 1. / 4) * constants::math::TwoPI, (0.17 - 1. / 4) * constants::math::TwoPI,
                                                (0.20 - 1. / 4) * constants::math::TwoPI,
                                                (0.22 - 1. / 4) * constants::math::TwoPI, (0.24 - 1. / 4) * constants::math::TwoPI,
                                                (0.26 - 1. / 4) * constants::math::TwoPI, (0.28 - 1. / 4) * constants::math::TwoPI,
                                                (0.30 - 1. / 4) * constants::math::TwoPI,
                                                (0.33 - 1. / 4) * constants::math::TwoPI, (0.36 - 1. / 4) * constants::math::TwoPI,
                                                (0.40 - 1. / 4) * constants::math::TwoPI, (0.45 - 1. / 4) * constants::math::TwoPI,
                                                (0.50 - 1. / 4) * constants::math::TwoPI,
                                                (0.55 - 1. / 4) * constants::math::TwoPI, (0.60 - 1. / 4) * constants::math::TwoPI,
                                                (0.64 - 1. / 4) * constants::math::TwoPI, (0.68 - 1. / 4) * constants::math::TwoPI,
                                                (0.71 - 1. / 4) * constants::math::TwoPI,
                                                (0.74 - 1. / 4) * constants::math::TwoPI, (0.76 - 1. / 4) * constants::math::TwoPI,
                                                (0.79 - 1. / 4) * constants::math::TwoPI,
                                                (0.82 - 1. / 4) * constants::math::TwoPI, (0.86 - 1. / 4) * constants::math::TwoPI,
                                                (0.90 - 1. / 4) * constants::math::TwoPI, (0.95 - 1. / 4) * constants::math::TwoPI,
                                                (1.00 - 1. / 4) * constants::math::TwoPI},
                                               "correlation bins DeltaPhi"};
    Configurable<std::vector<double>> binsDEta{"binsDEta",
                                               {(0 / 32. * 4 - 2) * absEtaMaxDefault,
                                                (2 / 32. * 4 - 2) * absEtaMaxDefault, (4 / 32. * 4 - 2) * absEtaMaxDefault,
                                                (6 / 32. * 4 - 2) * absEtaMaxDefault, (8 / 32. * 4 - 2) * absEtaMaxDefault,
                                                (9.5 / 32. * 4 - 2) * absEtaMaxDefault, (11 / 32. * 4 - 2) * absEtaMaxDefault,
                                                (12.5 / 32. * 4 - 2) * absEtaMaxDefault, (14 / 32. * 4 - 2) * absEtaMaxDefault,
                                                (15.5 / 32. * 4 - 2) * absEtaMaxDefault, (16.5 / 32. * 4 - 2) * absEtaMaxDefault,
                                                (18 / 32. * 4 - 2) * absEtaMaxDefault,
                                                (19.5 / 32. * 4 - 2) * absEtaMaxDefault, (21 / 32. * 4 - 2) * absEtaMaxDefault,
                                                (22.5 / 32. * 4 - 2) * absEtaMaxDefault, (24 / 32. * 4 - 2) * absEtaMaxDefault,
                                                (26 / 32. * 4 - 2) * absEtaMaxDefault, (28 / 32. * 4 - 2) * absEtaMaxDefault,
                                                (30 / 32. * 4 - 2) * absEtaMaxDefault, (32 / 32. * 4 - 2) * absEtaMaxDefault},
                                               "correlation bins DeltaEta"};
    Configurable<std::vector<double>> binsZPvBinning{"binsZPvBinning",
                                                     {-7, -5, -3, -1, 1, 3, 5, 7},
                                                     "zPv mixing bins"};
    Configurable<std::vector<double>> binsMultBinning{"binsMultBinning",
                                                      {-0.5, 10.5, 15.5, 20.5, 27.5, 42.5},
                                                      "multiplicity mixing bins"};
    Configurable<std::vector<double>> binsZPvBinningMcTrue{"binsZPvBinningMcTrue",
                                                           {-10000, 10000},
                                                           "zPv mixing bins for mc true"};
    Configurable<std::vector<double>> binsMultBinningMcTrue{"binsMultBinningMcTrue",
                                                            {-0.5, 16.5, 24.5, 31.5, 39.5, 64.5},
                                                            "multiplicity mixing bins for mc true"};
  } configHist;

  // further variables ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  enum class HistFillMode {
    PlainPt,
    PlainPhiEta,
    PlainZPvMult,
    Plain,
    Correlation,
    Mixing
  };

  enum class TriggerMemoryMode {
    Reco,
    True,
    TrueAssocEv,
    TrueAssocEvRecoPtTrig
  };

  struct HistPathEntry {
    HistFillMode fillMode;
    TruthLevel truthLevel;
    AssocType assocType;
    char const* path;
  };

  struct HistFillInfoPlainPt {
    double ptMax, ptAssoc;
    double weight = 1;
  };
  struct HistFillInfoPlainPhiEtaZPvMult {
    double ptMax, ptAssoc, phi, eta, zPv, mult;
    double weight = 1;
  };
  struct HistFillInfoCorr {
    double dPhi, dEta, ptTrig, ptAssoc, zPv, mult;
    double weight = 1;
  };

  std::mt19937 randomMt{3141u};
  int randomIntInInterval(int const lowEdgeInclusive, int const upEdgeInclusive)
  {
    return lowEdgeInclusive + (randomMt() % (upEdgeInclusive - lowEdgeInclusive + 1)); // speeds up task by factor > 2; modulo bias should not be relevant here
    // return std::uniform_int_distribution<int>(lowEdgeInclusive, upEdgeInclusive)(randomMt);
  }

  // objects to hold histograms
  HistogramRegistry histos{"histogramRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  // ccdb calls
  Service<ccdb::BasicCCDBManager> ccdb{};
  // for mc
  Service<framework::O2DatabasePDG> pdg{};

  // partitions
  SliceCache cache;

  // prepare for per collision slicing
  Preslice<aod::JetTracks> perColTracks = aod::jtrack::collisionId;
  Preslice<aod::Triggers> perColTriggers = aod::corr_common::corrCollisionId;
  Preslice<aod::Hadrons> perColHadrons = aod::corr_common::corrCollisionId;
  Preslice<aod::Pipms> perColPipms = aod::corr_common::corrCollisionId;
  Preslice<aod::PhotonPcms> perColPhotonPcms = aod::corr_common::corrCollisionId;
  Preslice<aod::PhotonPcmPairs> perColPhotonPcmPairs = aod::corr_common::corrCollisionId;
  Preslice<aod::JetParticles> perColMcParticles = aod::jmcparticle::mcCollisionId;
  Preslice<aod::TrueTriggers> perColTrueTriggers = aod::corr_common::jetMcCollisionId;

  // declare analysis variables

  // correction histograms
  TH1D* h1PtCorrectionTrigger = nullptr;
  TH1D* h1PtScalingK0Short = nullptr;

  // mixing trigger memory
  struct MixingTrigger {
    MixingTrigger(float pt, float phi, float eta)
      : fPt(pt), fPhi(phi), fEta(eta) {}
    float fPt, fPhi, fEta;
    [[nodiscard]] float pt() const { return fPt; }
    [[nodiscard]] float phi() const { return fPhi; }
    [[nodiscard]] float eta() const { return fEta; }
  };
  // class to handle trigger info from previous collisions (beyond single dataframe)
  // organised as zPv- and mult-bin matrix of deque (pt, phi, eta) to save triggers in buffer
  // last mult bin open to account for overflow and adjusted zVtx limits to account for rounding errors -> all events accounted for
  // (possibly replace by some advanced derived data method and O2 event mixing in future?)
  class MixingTriggerMemory
  {
   public:
    // finds bin that value belongs to (assumes ordered bin edges) (starts at 0; includes underflow (return -1) and overflow (return binEdges.size() - 1))
    // should be faster than some std binary search due to small number of bins (zPv, mult)
    static int findIntervalBin(double value, const std::vector<double>& binEdges)
    {
      if (binEdges.empty()) {
        throw std::invalid_argument("binEdges in findIntervalBin");
      }
      const int n = binEdges.size() - 1;
      if (value < binEdges[0]) {
        return -1; // underflow
      }
      for (int i_bin = 0; i_bin < n; i_bin++) {
        if (value < binEdges[i_bin + 1]) {
          return i_bin;
        }
      }
      return n; // overflow
    }

    MixingTriggerMemory(std::vector<double> const& binEdgesZPvIn, std::vector<double> const& binEdgesMultIn, size_t const nTriggerPerBinLimitIn)
      : binEdgesZPv(binEdgesZPvIn),
        nBinsZPv(binEdgesZPvIn.size() - 1),
        binEdgesMult(binEdgesMultIn),
        nBinsMult(binEdgesMultIn.size() - 1),
        nTriggerPerBinLimit(nTriggerPerBinLimitIn)
    {
      int const nEdgesMin = 2;
      if (binEdgesZPv.size() < nEdgesMin || binEdgesMult.size() < nEdgesMin) {
        throw std::invalid_argument("too few bin edges to define zPv or mult bins");
      }
      if (!std::is_sorted(binEdgesZPv.begin(), binEdgesZPv.end()) || !std::is_sorted(binEdgesMult.begin(), binEdgesMult.end())) {
        throw std::invalid_argument("bin edges must be sorted");
      }
      // prevent zPv over/underflow due to rounding errors
      binEdgesZPv.front() = binEdgesZPv.front() < 0 ? binEdgesZPv.front() * zPvRoundingErrorAdjust : binEdgesZPv.front() / zPvRoundingErrorAdjust;
      binEdgesZPv.back() = binEdgesZPv.back() > 0 ? binEdgesZPv.back() * zPvRoundingErrorAdjust : binEdgesZPv.back() / zPvRoundingErrorAdjust;
      // init correct size of zPv-mult matrix
      savedTriggersZPvMult.resize(nBinsZPv);
      for (size_t i_zPv = 0; i_zPv < nBinsZPv; i_zPv++) {
        savedTriggersZPvMult[i_zPv].resize(nBinsMult);
      }
    }

    // save trigger for mixing
    // up to nTriggerPerBinLimit stored in rolling buffer for each zPv-mult bin
    void saveTrigger(float const pt, float const phi, float const eta, double const zPv, double const mult)
    {
      int const iBinCorrZPv = getZPvBin(zPv);
      int const iBinCorrMult = getMultBin(mult);
      // fill
      std::deque<MixingTrigger>& triggerBin = savedTriggersZPvMult[iBinCorrZPv][iBinCorrMult];
      triggerBin.push_front(MixingTrigger{pt, phi, eta});
      if (triggerBin.size() > nTriggerPerBinLimit) {
        triggerBin.pop_back();
        return; // skip nTriggerBinMin update when unnecessary
      }
      // update nTriggerBinMin
      nTriggerBinMin = nTriggerPerBinLimit;
      for (size_t i_zPv = 0; i_zPv < nBinsZPv; i_zPv++) {
        for (size_t i_mult = 0; i_mult < nBinsMult; i_mult++) {
          nTriggerBinMin = std::min(nTriggerBinMin, savedTriggersZPvMult[i_zPv][i_mult].size());
        }
      }
    }

    // return deques of trigger pt, phi, eta in the given zPv/mult bin
    [[nodiscard]] std::deque<MixingTrigger> const& getTriggers(double const zPv, double const mult) const
    {
      int const iBinCorrZPv = getZPvBin(zPv);
      int const iBinCorrMult = getMultBin(mult);
      return savedTriggersZPvMult[iBinCorrZPv][iBinCorrMult];
    }

    [[nodiscard]] size_t getBinSizeMin() const
    {
      return nTriggerBinMin;
    }

   private:
    double zPvRoundingErrorAdjust = 1.0001;
    std::vector<double> binEdgesZPv;
    size_t nBinsZPv;
    std::vector<double> binEdgesMult;
    size_t nBinsMult;
    std::vector<std::vector<std::deque<MixingTrigger>>> savedTriggersZPvMult;
    size_t nTriggerPerBinLimit;
    size_t nTriggerBinMin = 0;

    [[nodiscard]] size_t getZPvBin(double const zPv) const
    {
      int const iBinInit = findIntervalBin(zPv, binEdgesZPv);
      if (iBinInit == -1 || iBinInit == static_cast<int>(nBinsZPv)) {
        throw std::runtime_error("zPv underflow or overflow in MixingTriggerMemory");
      }
      return iBinInit;
    }
    [[nodiscard]] size_t getMultBin(double const mult) const
    {
      int const iBinInit = findIntervalBin(mult, binEdgesMult);
      if (iBinInit == static_cast<int>(nBinsMult)) {
        return iBinInit - 1;
      }
      if (iBinInit == -1) {
        throw std::runtime_error("mult underflow in MixingTriggerMemory");
      }
      return iBinInit;
    }
  };

  // different buffers for mixing modes
  std::unique_ptr<MixingTriggerMemory> triggerMemoryReco;
  std::unique_ptr<MixingTriggerMemory> triggerMemoryTrue;
  std::unique_ptr<MixingTriggerMemory> triggerMemoryTrueAssocEv;
  std::unique_ptr<MixingTriggerMemory> triggerMemoryTrueAssocEvRecoPtTrig;

  // pcm cut class
  V0PhotonCut fV0PhotonCut;

  // functions ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ==========================================================================================================================================================================================

  // general (kenobi) /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  static std::vector<double> prependVariableWidth(std::vector<double> const& vec)
  {
    std::vector<double> result{VARIABLE_WIDTH};
    result.insert(result.end(), vec.begin(), vec.end());
    return result;
  }

  // get histograms from ccdb
  // save corrections from ccdb in histogram registry
  void initCcdbHistograms(int64_t const noLaterThanCcdb)
  {
    // trigger
    if (config.applyTriggerPtCorrection.value) {
      h1PtCorrectionTrigger = ccdb->getForTimeStamp<TH1D>(config.pathCcdbAnalysis.value + "/correction/trigger", noLaterThanCcdb);

      const double* bins = h1PtCorrectionTrigger->GetXaxis()->GetXbins()->GetArray();
      const AxisSpec axis{std::vector<double>(bins, bins + h1PtCorrectionTrigger->GetNbinsX() + 1), "#it{p}_{T}"};
      histos.add("usedObjectsCcdb/h1_pt_correction_trigger_ccdb", "h1_pt_correction_trigger_ccdb", kTH1D, {axis}, true);
      for (int iBin = 0; iBin <= h1PtCorrectionTrigger->GetNbinsX() + 1; iBin++) {
        histos.get<TH1>(HIST("usedObjectsCcdb/h1_pt_correction_trigger_ccdb"))->SetBinContent(iBin, h1PtCorrectionTrigger->GetBinContent(iBin));
        histos.get<TH1>(HIST("usedObjectsCcdb/h1_pt_correction_trigger_ccdb"))->SetBinError(iBin, h1PtCorrectionTrigger->GetBinError(iBin));
      }
    }
    // K0Short
    if (config.applyScalingK0ShortPt.value) {
      h1PtScalingK0Short = ccdb->getForTimeStamp<TH1D>(config.pathCcdbAnalysis.value + "/scaling/K0Short", noLaterThanCcdb);

      const double* bins = h1PtScalingK0Short->GetXaxis()->GetXbins()->GetArray();
      const AxisSpec axis{std::vector<double>(bins, bins + h1PtScalingK0Short->GetNbinsX() + 1), "#it{p}_{T}"};
      histos.add("usedObjectsCcdb/h1_pt_scaling_K0Short_ccdb", "h1_pt_scaling_K0Short_ccdb", kTH1D, {axis}, true);
      for (int iBin = 0; iBin <= h1PtScalingK0Short->GetNbinsX() + 1; iBin++) {
        histos.get<TH1>(HIST("usedObjectsCcdb/h1_pt_scaling_K0Short_ccdb"))->SetBinContent(iBin, h1PtScalingK0Short->GetBinContent(iBin));
        histos.get<TH1>(HIST("usedObjectsCcdb/h1_pt_scaling_K0Short_ccdb"))->SetBinError(iBin, h1PtScalingK0Short->GetBinError(iBin));
      }
    }
  }

  // create histograms
  void initHistograms()
  {
    // define axes
    const AxisSpec axisN{1, 0., 1., "#it{N}_{something}"};
    const AxisSpec axisCategories{32, 0., 32., "categories"};

    const AxisSpec axisZPv{configHist.nBinsZPvPlain, -7, 7, "#it{z}_{pv}"};
    const AxisSpec axisMult{configHist.nBinsMultPlain + 1, -0.5, configHist.nBinsMultPlain + 0.5, "multiplicity"};

    const AxisSpec axisPhi{configHist.nBinsPhi, 0, constants::math::TwoPI, "#it{#varphi}"};
    const AxisSpec axisEta{configHist.nBinsEta, -config.absEtaMax.value, config.absEtaMax.value, "#it{#eta}"};
    const AxisSpec axisDcaZ{configHist.nBinsDcaZ, -5, 5, "DCA_{z}"};
    const AxisSpec axisMgg{configHist.nBinsMgg, 0, 0.8, "#it{m}_{#gamma#gamma}"};

    const AxisSpec axisPtTrig{configHist.binsPtTrig, "#it{p}_{T}^{trig}"};
    const AxisSpec axisPtAssoc{configHist.binsPtAssoc, "#it{p}_{T}^{assoc}"};
    const AxisSpec axisDPhi{configHist.binsDPhi, "#Delta#it{#varphi}"};
    const AxisSpec axisDEta{configHist.binsDEta, "#Delta#it{#eta}"};
    const AxisSpec axisZPvBinning{configHist.binsZPvBinning, "#it{z}_{pv} correlation binning"};
    const AxisSpec axisMultBinning{configHist.binsMultBinning, "multiplicity correlation binning"};
    const AxisSpec axisZPvBinningMcTrue{configHist.binsZPvBinningMcTrue, "#it{z}_{pv} correlation binning for mc true"};
    const AxisSpec axisMultBinningMcTrue{configHist.binsMultBinningMcTrue, "multiplicity correlation binning for mc true"};

    // main info for assocType (spectra and correlation)
    auto const addMainAssocHists =
      [&](std::string const& recoOrMc, std::string const& nameTruthLevel, std::string const& nameAssoc, bool const resolvedInfo, bool const trueColInfo) {
        if (resolvedInfo) {
          histos.add(std::format("{}/plain/h4_ptTrigPtAssocPhiEta_{}_{}", recoOrMc, nameTruthLevel, nameAssoc).data(), std::format("h4_ptTrigPtAssocPhiEta_{}_{}", nameTruthLevel, nameAssoc).data(),
                     kTHnSparseD, {axisPtTrig, axisPtAssoc, axisPhi, axisEta}, true);
          histos.add(std::format("{}/plain/h4_ptTrigPtAssocZPvMult_{}_{}", recoOrMc, nameTruthLevel, nameAssoc).data(), std::format("h4_ptTrigPtAssocZPvMult_{}_{}", nameTruthLevel, nameAssoc).data(),
                     kTHnSparseD, {axisPtTrig, axisPtAssoc, axisZPv, axisMult}, true);
        } else {
          histos.add(std::format("{}/plain/h2_ptTrigPtAssoc_{}_{}", recoOrMc, nameTruthLevel, nameAssoc).data(), std::format("h2_ptTrigPtAssoc_{}_{}", nameTruthLevel, nameAssoc).data(),
                     kTHnSparseD, {axisPtTrig, axisPtAssoc}, true);
        }
        if (!trueColInfo) {
          histos.add(std::format("{}/corr/h6_corr_{}_{}", recoOrMc, nameTruthLevel, nameAssoc).data(), std::format("h6_corr_{}_{}", nameTruthLevel, nameAssoc).data(),
                    kTHnSparseD, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
          histos.add(std::format("{}/corr/h6_mix_{}_{}", recoOrMc, nameTruthLevel, nameAssoc).data(), std::format("h6_mix_{}_{}", nameTruthLevel, nameAssoc).data(),
                    kTHnSparseD, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
        } else {
          histos.add(std::format("{}/corr/h6_corr_{}_{}", recoOrMc, nameTruthLevel, nameAssoc).data(), std::format("h6_corr_{}_{}", nameTruthLevel, nameAssoc).data(),
                    kTHnSparseD, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinningMcTrue, axisMultBinningMcTrue}, true);
          histos.add(std::format("{}/corr/h6_mix_{}_{}", recoOrMc, nameTruthLevel, nameAssoc).data(), std::format("h6_mix_{}_{}", nameTruthLevel, nameAssoc).data(),
                    kTHnSparseD, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinningMcTrue, axisMultBinningMcTrue}, true);
        }
      };

    // reco info
    histos.add("reco/info/h1_nEvents", "h1_nEvents", kTH1D, {axisCategories}, true);
    histos.add("reco/info/h3_nEv_ptTrigZPvMult_reco", "h3_nEv_ptTrigZPvMult_reco", kTHnSparseD, {axisPtTrig, axisZPv, axisMult}, true);

    // reco (correlation) analysis
    histos.add("reco/corr/h3_ptPhiEta_trig", "h3_ptPhiEta_trig", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
    for (std::string const& nameAssoc : std::initializer_list<std::string>{"hadron", "pipm", "photonPcm"}) {
      addMainAssocHists("reco", "reco", nameAssoc, true, false);
    }
    for (std::string const& nameAssoc : std::initializer_list<std::string>{"pi0PcmPeak", "pi0PcmSide", "etaPcmPeak", "etaPcmSide"}) {
      addMainAssocHists("reco", "reco", nameAssoc, false, false);
    }

    // photonPcm pairs
    // mgg
    histos.add("reco/plain/h5_ptTrigPtAssocMggZPvMult_reco_photonPcmPair", "h5_ptTrigPtAssocMggZPvMult_reco_photonPcmPair",
               kTHnSparseD, {axisPtTrig, axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);
    histos.add("reco/corr/h5_ptTrigPtAssocMggZPvMult_reco_photonPcmPair", "h5_ptTrigPtAssocMggZPvMult_reco_photonPcmPair",
               kTHnSparseD, {axisPtTrig, axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);
    // event mixing for photon pairs
    histos.add("reco/plain/h4_ptMggZPvMult_photonPcmPair_evMix", "h4_ptMggZPvMult_photonPcmPair_evMix", kTHnSparseD, {axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);
    // pseudo mgg
    // pi0Pcm
    histos.add("mc/plain/h5_ptTrigPtAssocMggZPvMult_pseudoReco_pi0Pcm", "h5_ptTrigPtAssocMggZPvMult_pseudoReco_pi0Pcm",
               kTHnSparseD, {axisPtTrig, axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);
    histos.add("mc/corr/h5_ptTrigPtAssocMggZPvMult_pseudoReco_pi0Pcm", "h5_ptTrigPtAssocMggZPvMult_pseudoReco_pi0Pcm",
               kTHnSparseD, {axisPtTrig, axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);
    histos.add("mc/plain/h5_ptTrigPtAssocMggZPvMult_pseudoReco_etaPcm", "h5_ptTrigPtAssocMggZPvMult_pseudoReco_etaPcm",
               kTHnSparseD, {axisPtTrig, axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);
    histos.add("mc/corr/h5_ptTrigPtAssocMggZPvMult_pseudoReco_etaPcm", "h5_ptTrigPtAssocMggZPvMult_pseudoReco_etaPcm",
               kTHnSparseD, {axisPtTrig, axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);

    // mc info
    histos.add("mc/info/h3_nEv_ptTrigZPvMult_true", "h3_nEv_ptTrigZPvMult_true", kTHnSparseD, {axisPtTrig, axisZPv, axisMult}, true);
    histos.add("mc/info/h3_nEv_ptTrigZPvMult_trueAssocEv", "h3_nEv_ptTrigZPvMult_trueAssocEv", kTHnSparseD, {axisPtTrig, axisZPv, axisMult}, true);
    histos.add("mc/info/h3_nEv_ptTrigZPvMult_trueAssocEvRecoPtTrig", "h3_nEv_ptTrigZPvMult_trueAssocEvRecoPtTrig", kTHnSparseD, {axisPtTrig, axisZPv, axisMult}, true);

    // true++ (event selection)
    for (std::string const& nameTruthLevel : std::initializer_list<std::string>{"true", "trueAssocEv", "trueAssocEvRecoPtTrig"}) {
      histos.add(std::format("mc/corr/h3_ptPhiEta_trig_{}", nameTruthLevel).data(), std::format("h3_ptPhiEta_trig_{}", nameTruthLevel).data(),
                 kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
      for (std::string const& nameAssoc : std::initializer_list<std::string>{"hadron", "pipm", "photon"}) {
        addMainAssocHists("mc", nameTruthLevel, nameAssoc, true, true);
      }
      for (std::string const& nameAssoc : std::initializer_list<std::string>{"pi0", "eta"}) {
        addMainAssocHists("mc", nameTruthLevel, nameAssoc, false, true);
      }
    }
    // extra truthLevels for photonPcmPair
    for (std::string const& nameTruthLevel : std::initializer_list<std::string>{"trueDecaySelfTrig", "trueMeasDecay", "trueGeoAcc"}) {
      for (std::string const& nameAssoc : std::initializer_list<std::string>{"pi0", "eta"}) {
        addMainAssocHists("mc", nameTruthLevel, nameAssoc, false, true);
      }
    }
    // assocTypes just on true level
    for (std::string const& nameAssoc : std::initializer_list<std::string>{"photonDecay", "photonDirect", "photonPi0", "photonEta",
                                                                           "omega", "photonOmega", "etaPrime", "photonEtaPrime", "photonOtherMother"}) {
      addMainAssocHists("mc", "true", nameAssoc, false, true);
    }

    // reco++ (MC info)
    // std assoc
    for (std::string const& nameTruthLevel : std::initializer_list<std::string>{"recoAssocEv", "recoPure", "recoPureTruePtAssoc"}) {
      for (std::string const& nameAssoc : std::initializer_list<std::string>{"hadron", "pipm", "photonPcm"}) {
        addMainAssocHists("mc", nameTruthLevel, nameAssoc, true, false);
      }
    }
    // photonPcmPair assoc
    for (std::string const& nameTruthLevel : std::initializer_list<std::string>{"pseudoReco", "pseudoRecoAssocEv", "pseudoRecoPure", "pseudoRecoPureTruePtAssoc"}) {
      for (std::string const& nameAssoc : std::initializer_list<std::string>{"pi0Pcm", "etaPcm"}) {
        addMainAssocHists("mc", nameTruthLevel, nameAssoc, false, false);
      }
    }
    // impurity details
    histos.add("mc/plain/h3_ptTrigPtAssocCategory_recoImpurities_photonPcm", "h3_ptTrigPtAssocCategory_recoImpurities_photonPcm",
               kTHnSparseD, {axisPtTrig, axisPtAssoc, axisCategories}, true);
    histos.add("mc/plain/h3_ptTrigPtAssocCategory_pseudoRecoImpurities_pi0Pcm", "h3_ptTrigPtAssocCategory_pseudoRecoImpurities_pi0Pcm",
               kTHnSparseD, {axisPtTrig, axisPtAssoc, axisCategories}, true);

    // check reco-mc matching
    histos.add("mc/info/h1_pt_noMatch_event", "h1_pt_noMatch_event", kTH1D, {axisN}, true);
    histos.add("mc/info/h1_pt_noMatch_trigger", "h1_pt_noMatch_trigger", kTH1D, {axisPtTrig}, true);

    // pcm dcaZ
    histos.add("reco/plain/h5_ptTrigPtAssocDcaZZPvMult_photonPcm", "h5_ptTrigPtAssocDcaZZPvMult_photonPcm",
               kTHnSparseD, {axisPtTrig, axisPtAssoc, axisDcaZ, axisZPvBinning, axisMultBinning}, true);
    histos.add("mc/plain/h5_ptTrigPtAssocDcaZZPvMult_photonPcm", "h5_ptTrigPtAssocDcaZZPvMult_photonPcm",
               kTHnSparseD, {axisPtTrig, axisPtAssoc, axisDcaZ, axisZPvBinning, axisMultBinning}, true);
  }

  // init pcm cut class
  void initV0PhotonCut()
  {
    fV0PhotonCut = V0PhotonCut("fV0PhotonCut", "fV0PhotonCut");

    // for v0
    fV0PhotonCut.SetV0PtRange(configPcmCuts.cfg_min_pt_v0, configPcmCuts.cfg_max_pt_v0);
    fV0PhotonCut.SetV0EtaRange(configPcmCuts.cfg_min_eta_v0, configPcmCuts.cfg_max_eta_v0);
    fV0PhotonCut.SetMinCosPA(configPcmCuts.cfg_min_cospa);
    fV0PhotonCut.SetMaxPCA(configPcmCuts.cfg_max_pca);
    fV0PhotonCut.SetMaxChi2KF(configPcmCuts.cfg_max_chi2kf);
    fV0PhotonCut.SetRxyRange(configPcmCuts.cfg_min_v0radius, configPcmCuts.cfg_max_v0radius, configPcmCuts.cfg_midL_v0radius, configPcmCuts.cfg_midH_v0radius);
    fV0PhotonCut.SetAPRange(configPcmCuts.cfg_max_alpha_ap, configPcmCuts.cfg_max_qt_ap);
    fV0PhotonCut.RejectITSib(configPcmCuts.cfg_reject_v0_on_itsib);

    // for track
    fV0PhotonCut.SetMinNClustersTPC(configPcmCuts.cfg_min_ncluster_tpc);
    fV0PhotonCut.SetMinNCrossedRowsTPC(configPcmCuts.cfg_min_ncrossedrows);
    fV0PhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fV0PhotonCut.SetMaxFracSharedClustersTPC(configPcmCuts.cfg_max_frac_shared_clusters_tpc);
    fV0PhotonCut.SetChi2PerClusterTPC(0.0, configPcmCuts.cfg_max_chi2tpc);
    fV0PhotonCut.SetTPCNsigmaElRange(configPcmCuts.cfg_min_TPCNsigmaEl, configPcmCuts.cfg_max_TPCNsigmaEl);
    fV0PhotonCut.SetChi2PerClusterITS(-1e+10, configPcmCuts.cfg_max_chi2its);
    fV0PhotonCut.SetNClustersITS(0, 7);
    fV0PhotonCut.SetMeanClusterSizeITSob(0.0, 16.0);
    fV0PhotonCut.SetDisableITSonly(configPcmCuts.cfg_disable_itsonly_track);
    fV0PhotonCut.SetDisableTPConly(configPcmCuts.cfg_disable_tpconly_track);
    fV0PhotonCut.SetRequireITSTPC(configPcmCuts.cfg_require_v0_with_itstpc);
    fV0PhotonCut.SetRequireITSonly(configPcmCuts.cfg_require_v0_with_itsonly);
    fV0PhotonCut.SetRequireTPConly(configPcmCuts.cfg_require_v0_with_tpconly);
  }

  // selections ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // mc split event selection based in the thrid decimal of mc-true posZ
  bool checkSplitMcEventSelection(auto const& collision)
  {
    // select based on configurables
    return collision.thirdDecimalTruePosZ() % config.splitMcEvN.value == config.splitMcEvSelect.value - 1;
  }

  // total event selection (basic selections from producer and analysis level selections)
  bool totalEvSel(auto const& collision)
  {
    return collision.selEv() && checkSplitMcEventSelection(collision);
  }

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

  // checks if mcParticle daughters are two photons
  bool isDecayToGG(auto const& mcParticle)
  {
    auto const& daughters = mcParticle.template daughters_as<aod::JetParticles>();
    constexpr int NDaughtersToGG = 2;
    if (daughters.size() != NDaughtersToGG) {
      return false;
    }
    auto daughter = daughters.begin();
    auto daughterEnd = daughters.end();
    while (daughter != daughterEnd) {
      if ((*daughter).pdgCode() != PDG_t::kGamma) {
        return false;
      }
      daughter++;
    }
    return true;
  }
  // check if particle has associated mother in parent tree
  // includes associated particle checks
  bool checkForAssocDecay(auto const& mcParticle, AssocType const assocType)
  {
    int const pdgCode = [&]() {
      if (assocType == AssocType::Pi0Pcm) {
        return static_cast<int>(PDG_t::kPi0);
      } else if (assocType == AssocType::EtaPcm) {
        return static_cast<int>(constants::physics::Pdg::kEta);
      } else {
        throw std::runtime_error("invalid assocType for mother search");
      }
    }();

    if (!mcParticle.has_mothers()) {
      return false;
    }
    auto const mothers = mcParticle.template mothers_as<aod::JetParticles>();
    auto mother = mothers.begin();
    auto motherEnd = mothers.end();
    while (mother != motherEnd) {
      // exclude pre-hadonisation
      if ((*mother).producedByGenerator() && std::abs((*mother).getGenStatusCode()) < 80) {
        return false;
      }
      // check valid associated particle
      if ((*mother).pdgCode() == pdgCode || (*mother).pdgCode() == -pdgCode) {
        if (std::abs((*mother).eta()) < config.absEtaMax.value && isDecayPrimary(*mother)) {
          return true;
        }
      }
      if (checkForAssocDecay(*mother, assocType)) {
        return true;
      }
      mother++;
    }
    return false;
  }
  // get plain self trigger associated map for event
  template <TruthLevel TruthLevelV>
  std::map<AssocType, bool> getDecaySelfTrigPlainMap(auto const& triggersThisEvent)
  {
    std::map<AssocType, bool> decaySelfTrigPlainMap;

    for (AssocType const assocType : {AssocType::Pi0Pcm, AssocType::EtaPcm}) {
      decaySelfTrigPlainMap[assocType] = [&]() {
        auto trigger = triggersThisEvent.begin();
        auto triggerEnd = triggersThisEvent.end();
        // false for no triggers
        if (trigger == triggerEnd) {
          return false;
        }
        // check every trigger for self trigger and return false if unrelated (event triggers independently)
        while (trigger != triggerEnd) {
          if constexpr (TruthLevelV == TruthLevel::TrueAssocEvRecoPtTrig) {
            if (!isMcMatchedTrigger(*trigger)) {
              return false;
            }
            if (!checkForAssocDecay((*trigger).template jetTrack_as<aod::JetTracksMCD>().mcParticle(), assocType)) {
              return false;
            }
          } else if constexpr (TruthLevelV == TruthLevel::TrueAssocEv || TruthLevelV == TruthLevel::True) {
            if (!checkForAssocDecay((*trigger).jetMcParticle(), assocType)) {
              return false;
            }
          } else {
            static_assert(AlwaysFalseValue<TruthLevelV>, "truthLevel not found");
          }
          trigger++;
        }
        return true;
      }();
    }

    return decaySelfTrigPlainMap;
  }
  // recursive helper to exclude loopDeLoops
  template <typename T_mcParticle>
  bool producesTrigger(T_mcParticle const& mcParticle, T_mcParticle const& triggerParticle, std::unordered_set<int64_t>& visited)
  {
    if (!visited.insert(mcParticle.globalIndex()).second) {
      return false;
    }
    if (!mcParticle.has_daughters()) {
      return false;
    }
    auto const daughters = mcParticle.template daughters_as<aod::JetParticles>();
    auto daughter = daughters.begin();
    auto daughterEnd = daughters.end();
    while (daughter != daughterEnd) {
      if ((*daughter).globalIndex() == triggerParticle.globalIndex()) {
        return true;
      }
      if (producesTrigger(*daughter, triggerParticle, visited)) {
        return true;
      }
      daughter++;
    }
    return false;
  }
  // check if particle has trigger in daughter tree
  template <typename T_mcParticle>
  bool producesTrigger(T_mcParticle const& mcParticle, T_mcParticle const& triggerParticle)
  {
    std::unordered_set<int64_t> visited;
    return producesTrigger(mcParticle, triggerParticle, visited);
  }

  // checks if daughters are in acceptance
  bool checkDaughtersInAcceptance(auto const& mcParticle, std::pair<double, double> const etaRange, std::pair<double, double> const phiRange = {0, constants::math::TwoPI})
  {
    auto const& daughters = mcParticle.template daughters_as<aod::JetParticles>();
    auto daughter = daughters.begin();
    auto daughterEnd = daughters.end();
    while (daughter != daughterEnd) {
      if ((*daughter).eta() < etaRange.first || (*daughter).eta() > etaRange.second) {
        return false;
      }
      if ((*daughter).phi() < phiRange.first || (*daughter).phi() > phiRange.second) {
        return false;
      }
      daughter++;
    }
    return true;
  }

  // assocType dependent checks ---------------------------------------------------------------------------------------------------------------------------------------------------------------

  // check trivial self correlation for reconstruction
  bool isRecoSelfTrigger(auto const& trigger, auto const& associated) {
    for (auto const& intersectTriggerId : associated.intersectTriggerIds()) {
      if (intersectTriggerId == trigger.globalIndex()) {
        return true;
      }
    }
    return false;
  };

  // check reco-mc match
  bool isMcMatchedCollision(auto const& collision) {
    return collision.has_jetMcCollision();
  };
  bool isMcMatchedTrigger(auto const& trigger) {
    return trigger.template jetTrack_as<aod::JetTracksMCD>().has_mcParticle();
  };

  // check reco truthLevel reached
  template <TruthLevel TruthLevelV, AssocType AssocTypeV>
  bool checkTruthLevelReco(auto const& associated) {
    // std
    static constexpr std::array truthLevelsRecoStdRanked{
      TruthLevel::Reco,
      TruthLevel::RecoMcMatched,
      TruthLevel::RecoAssocEv,
      TruthLevel::RecoPure
    };
    // diPhoton
    static constexpr std::array truthLevelsRecoDiPhotonRanked{
      TruthLevel::Reco,
      TruthLevel::RecoMcMatched,
      TruthLevel::PseudoReco,
      TruthLevel::PseudoRecoAssocEv,
      TruthLevel::PseudoRecoPure
    };
    // position of truthLevel in chain, -1 if not in ranked truthLevels
    auto constexpr getRank = [](TruthLevel const truthLevel, auto const& truthLevelsRanked) constexpr -> int {
      for (size_t iRank = 0; iRank < truthLevelsRanked.size(); iRank++) {
        if (truthLevelsRanked[iRank] == truthLevel) {
          return static_cast<int>(iRank);
        }
      }
      return -1;
    };
    // TruthLevelV compared to the stored maximum level
    auto const isReached = [&getRank](TruthLevel const truthLevelRecoMax, auto const& truthLevelsRanked) {
      int const rankRecoToCheck = getRank(TruthLevelV, truthLevelsRanked);
      int const rankRecoMax = getRank(truthLevelRecoMax, truthLevelsRanked);
      if (rankRecoMax < 0) {
        throw std::runtime_error("rankRecoMax not found");
      }
      return rankRecoToCheck <= rankRecoMax;
    };
    // assocType cases
    if constexpr (AssocTypeV == AssocType::Hadron) {
      static_assert(getRank(TruthLevelV, truthLevelsRecoStdRanked) >= 0, "truthLevel not in std reco chain");
      TruthLevel const truthLevelRecoMax = static_cast<TruthLevel>(associated.truthLevelRecoHadron());
      return isReached(truthLevelRecoMax, truthLevelsRecoStdRanked);
    } else if constexpr (AssocTypeV == AssocType::Pipm) {
      static_assert(getRank(TruthLevelV, truthLevelsRecoStdRanked) >= 0, "truthLevel not in std reco chain");
      TruthLevel const truthLevelRecoMax = static_cast<TruthLevel>(associated.truthLevelRecoPipm());
      return isReached(truthLevelRecoMax, truthLevelsRecoStdRanked);
    } else if constexpr (AssocTypeV == AssocType::PhotonPcm) {
      static_assert(getRank(TruthLevelV, truthLevelsRecoStdRanked) >= 0, "truthLevel not in std reco chain");
      TruthLevel const truthLevelRecoMax = static_cast<TruthLevel>(associated.truthLevelRecoPhotonPcm());
      return isReached(truthLevelRecoMax, truthLevelsRecoStdRanked);
    } else if constexpr (AssocTypeV == AssocType::Pi0Pcm) {
      static_assert(getRank(TruthLevelV, truthLevelsRecoDiPhotonRanked) >= 0, "truthLevel not in diPhoton reco chain");
      TruthLevel const truthLevelRecoMax = static_cast<TruthLevel>(associated.truthLevelRecoPi0Pcm());
      return isReached(truthLevelRecoMax, truthLevelsRecoDiPhotonRanked);
    } else if constexpr (AssocTypeV == AssocType::EtaPcm) {
      static_assert(getRank(TruthLevelV, truthLevelsRecoDiPhotonRanked) >= 0, "truthLevel not in diPhoton reco chain");
      TruthLevel const truthLevelRecoMax = static_cast<TruthLevel>(associated.truthLevelRecoEtaPcm());
      return isReached(truthLevelRecoMax, truthLevelsRecoDiPhotonRanked);
    } else {
      static_assert(AlwaysFalseValue<AssocTypeV>, "assocType not found");
    }
  }

  // reco/systematic cuts ---------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // energy asymmetry for photon pair
  bool checkEnergyAsymDiPhoton(auto const& v0Photon1, auto const& v0Photon2)
  {
    double const e1 = v0Photon1.pt() * std::cosh(v0Photon1.eta());
    double const e2 = v0Photon2.pt() * std::cosh(v0Photon2.eta());
    double const energyAsym = std::abs((e1 - e2) / (e1 + e2));
    return energyAsym < config.energyAsymmPhotonPairMax.value;
  }

  // add pmeventId() due to unguarded call in V0PhotonCut.h qa; selection unaffected
  template <typename T_V0PhotonType>
  struct V0PhotonNoEvent : T_V0PhotonType {
    explicit V0PhotonNoEvent(T_V0PhotonType const& v0PhotonType) : T_V0PhotonType(v0PhotonType) {}
    [[nodiscard]] int pmeventId() const { return -1; }
  };
  bool isSelectedV0Photon(auto const& v0Photon)
  {
    using V0PhotonType = std::decay_t<decltype(v0Photon)>;
    return fV0PhotonCut.template IsSelected<V0PhotonNoEvent<V0PhotonType>, aod::V0Legs>(V0PhotonNoEvent<V0PhotonType>{v0Photon});
  }

  bool isSelectedPhotonPcm(auto const& photonPcm)
  {
    return isSelectedV0Photon(photonPcm.v0PhotonKF());
  }
  bool isSelectedPhotonPcmPair(auto const& photonPcmPair)
  {
    if (!isSelectedV0Photon(photonPcmPair.v0PhotonKF1()) || !isSelectedV0Photon(photonPcmPair.v0PhotonKF2())) {
      return false;
    }
    if (!checkEnergyAsymDiPhoton(photonPcmPair.v0PhotonKF1(), photonPcmPair.v0PhotonKF2())) {
      return false;
    }
    return true;
  }

  template <AssocType AssocTypeV>
  bool isSelectedAssoc(auto const& associated) {
    if constexpr (AssocTypeV == AssocType::Hadron) {
      return true;
    } else if constexpr (AssocTypeV == AssocType::Pipm) {
      return true;
    } else if constexpr (AssocTypeV == AssocType::PhotonPcm) {
      return isSelectedPhotonPcm(associated);
    } else if constexpr (AssocTypeV == AssocType::PhotonPcmPair) {
      return isSelectedPhotonPcmPair(associated);
    } else {
      static_assert(AlwaysFalseValue<AssocTypeV>, "assocType not found");
    }
  }

  // analysis helpers /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // select standard/main histogram compile-time path from array
  static consteval char const* getHistPathMain(auto const& histPaths, HistFillMode const fillMode, TruthLevel const truthLevel, AssocType const assocType) {
    for (auto const& entry : histPaths) {
      if (entry.fillMode == fillMode && entry.truthLevel == truthLevel && entry.assocType == assocType) {
        return entry.path;
      }
    }
    throw "path not found";
  }

  // fill main histogram
  template <auto const& HistPaths, HistFillMode HistFillModeV, TruthLevel TruthLevelV, AssocType AssocTypeV>
  void fillHistMain(auto const& histFillInfo) {
    if constexpr (HistFillModeV == HistFillMode::Plain) {
      if constexpr (AssocTypeV == AssocType::Hadron || AssocTypeV == AssocType::Pipm || AssocTypeV == AssocType::PhotonPcm) {
        histos.fill(HIST(getHistPathMain(HistPaths, HistFillMode::PlainPhiEta, TruthLevelV, AssocTypeV)),
                    histFillInfo.ptMax, histFillInfo.ptAssoc, histFillInfo.phi, histFillInfo.eta,
                    histFillInfo.weight);
        histos.fill(HIST(getHistPathMain(HistPaths, HistFillMode::PlainZPvMult, TruthLevelV, AssocTypeV)),
                    histFillInfo.ptMax, histFillInfo.ptAssoc, histFillInfo.zPv, histFillInfo.mult,
                    histFillInfo.weight);
      } else {
        histos.fill(HIST(getHistPathMain(HistPaths, HistFillMode::PlainPt, TruthLevelV, AssocTypeV)),
                    histFillInfo.ptMax, histFillInfo.ptAssoc,
                    histFillInfo.weight);
      }
    } else if constexpr (HistFillModeV == HistFillMode::Correlation) {
      histos.fill(HIST(getHistPathMain(HistPaths, HistFillMode::Correlation, TruthLevelV, AssocTypeV)),
                  histFillInfo.dPhi, histFillInfo.dEta, histFillInfo.ptTrig, histFillInfo.ptAssoc, histFillInfo.zPv, histFillInfo.mult,
                  histFillInfo.weight);
    } else if constexpr (HistFillModeV == HistFillMode::Mixing) {
      histos.fill(HIST(getHistPathMain(HistPaths, HistFillMode::Mixing, TruthLevelV, AssocTypeV)),
                  histFillInfo.dPhi, histFillInfo.dEta, histFillInfo.ptTrig, histFillInfo.ptAssoc, histFillInfo.zPv, histFillInfo.mult,
                  histFillInfo.weight);
    } else {
      static_assert(AlwaysFalseValue<HistFillModeV>, "fillMode not valid");
    }
  }

  // tcm component: exponential for low pt
  // parameter shape: [m0_h, A_e, T_e]
  double tcmInvariantSoftExp(double const* pt, double const* par) {
    double const ET_kin = std::sqrt(std::pow(pt[0], 2) + std::pow(par[0], 2)) - par[0];
    return std::abs(par[1]) * std::exp(-1 * ET_kin / par[2]);
  }
  // tcm component: power law for high pt
  // parameter shape: [A_p, T_p, N_p]
  double tcmInvariantHardPower(double const* pt, double const* par) {
    return std::abs(par[0]) / std::pow(1 + std::pow(pt[0], 2) / (std::pow(par[1], 2) * par[2]), par[2]);
  }
  // fit function for invariant pT spectrum
  // parameter shape: [m0_h, A_e, T_e, A_p, T_p, N_p]
  double tcmInvariantFit(double const* pt, double const* par) {
    int const nParSoftExp = 3;
    return tcmInvariantSoftExp(pt, par) + tcmInvariantHardPower(pt, par + nParSoftExp);
  }
  // tcm function (non-invariant)
  double tcmParam(double const pt, std::vector<double> const& par) {
    return pt * tcmInvariantFit(&pt, par.data());
  }

  // pt weight
  double ptWeightTcms(double const pt, std::vector<double> const& parNum, std::vector<double> const& parDenom) {
    return tcmParam(pt, parNum) / tcmParam(pt, parDenom);
  }

  // rescale of pt spectra for mc

  template <AssocType AssocTypeV>
  double getMcAssocPtScalingFromPt(double const ptTrue) {
    auto const getWeight = [this](double const pt, std::vector<double> const& ptRange, std::vector<double> const& parNum, std::vector<double> const& parDenom) {
      if (pt < ptRange.at(0)) {
        return ptWeightTcms(ptRange.at(0), parNum, parDenom);
      }
      if (pt > ptRange.at(1)) {
        return ptWeightTcms(ptRange.at(1), parNum, parDenom);
      }
      return ptWeightTcms(pt, parNum, parDenom);
    };
    if constexpr (AssocTypeV == AssocType::PhotonPcm) {
      return config.applyScalingPhotonPt.value ? getWeight(ptTrue, config.ptRangeScalingPhoton.value, config.parScalingPhotonPtData.value, config.parScalingPhotonPtMc.value) : 1;
    } else if constexpr (AssocTypeV == AssocType::Pi0Pcm) {
      return config.applyScalingPi0Pt.value ? getWeight(ptTrue, config.ptRangeScalingPi0.value, config.parScalingPi0PtData.value, config.parScalingPi0PtMc.value) : 1;
    } else if constexpr (AssocTypeV == AssocType::EtaPcm) {
      return config.applyScalingEtaPt.value ? getWeight(ptTrue, config.ptRangeScalingEta.value, config.parScalingEtaPtData.value, config.parScalingEtaPtMc.value) : 1;
    } else {
      static_assert(AlwaysFalseValue<AssocTypeV>, "assocType not found");
    }
  }
  template <AssocType AssocTypeV>
  double getMcAssocPtScaling(auto const& associated) {
    if constexpr (AssocTypeV == AssocType::Hadron || AssocTypeV == AssocType::Pipm) {
      return 1;
    } else if constexpr (AssocTypeV == AssocType::PhotonPcm) {
      if (!config.applyScalingPhotonPt.value) return 1;
      if (checkTruthLevelReco<TruthLevel::RecoPure, AssocTypeV>(associated)) {
        return getMcAssocPtScalingFromPt<AssocTypeV>(associated.ptTruePure());
      } else {
        return 1;
      }
    } else if constexpr (AssocTypeV == AssocType::PhotonPcmPair) {
      if (!config.applyScalingPi0Pt.value && !config.applyScalingEtaPt.value) return 1;
      if (checkTruthLevelReco<TruthLevel::PseudoRecoPure, AssocType::Pi0Pcm>(associated)) {
        return getMcAssocPtScalingFromPt<AssocType::Pi0Pcm>(associated.ptTruePure());
      } else if (checkTruthLevelReco<TruthLevel::PseudoRecoPure, AssocType::EtaPcm>(associated)) {
        return getMcAssocPtScalingFromPt<AssocType::EtaPcm>(associated.ptTruePure());
      } else {
        return 1;
      }
    } else {
      static_assert(AlwaysFalseValue<AssocTypeV>, "assocType not found");
    }
  }

  // get histogram y-quantity at x-value
  double getH1ValueAt(auto const*const h1, double const value)
  {
    return h1->GetBinContent(h1->FindFixBin(value));
  }

  // K0Short mc scaling
  double getMcK0ShortPtScalingFromPt(double const ptTrue) {
    return config.applyScalingK0ShortPt.value ? getH1ValueAt(h1PtScalingK0Short, ptTrue) : 1;
  }
  template <AssocType AssocTypeV>
  double getMcK0ShortPtScaling (auto const& associated) {
    if constexpr (AssocTypeV == AssocType::Hadron || AssocTypeV == AssocType::Pipm) {
      return 1;
    } else if constexpr (AssocTypeV == AssocType::PhotonPcm) {
      if (associated.impurityPhotonPcm() == static_cast<int>(PhotonImpurity::K0Short)) {
        return getMcK0ShortPtScalingFromPt(associated.ptTrueOrigin());
      } else {
        return 1;
      }
    } else if constexpr (AssocTypeV == AssocType::PhotonPcmPair) {
      if (associated.impurityPi0Pcm() == static_cast<int>(PhotonImpurity::K0Short)) {
        return getMcK0ShortPtScalingFromPt(associated.ptTrueOrigin());
      } else {
        return 1;
      }
    } else {
      static_assert(AlwaysFalseValue<AssocTypeV>, "assocType not found");
    }
  }

  // correction trigger as function of pt
  double getTriggerPtCorrection(double const pt)
  {
    return config.applyTriggerPtCorrection.value ? getH1ValueAt(h1PtCorrectionTrigger, pt) : 1;
  }

  // performs 'phi1 - phi2' and pushes it into the interval [-pi/2, 3pi/2]
  double getDeltaPhi(double const phi1, double const phi2)
  {
    return RecoDecay::constrainAngle(phi1 - phi2, -1 * constants::math::PIHalf);
  }

  // check if invariant mass range
  template <AssocType AssocTypeV>
  bool isInMassRange(double const mgg)
  {
    if constexpr (AssocTypeV == AssocType::Pi0PcmPeak) {
      return mgg > config.pi0PcmPeakMassRange.value.at(0) && mgg < config.pi0PcmPeakMassRange.value.at(1);
    } else if constexpr (AssocTypeV == AssocType::Pi0PcmSide) {
      return mgg > config.pi0PcmSideMassRange.value.at(0) && mgg < config.pi0PcmSideMassRange.value.at(1);
    } else if constexpr (AssocTypeV == AssocType::EtaPcmPeak) {
      return mgg > config.etaPcmPeakMassRange.value.at(0) && mgg < config.etaPcmPeakMassRange.value.at(1);
    } else if constexpr (AssocTypeV == AssocType::EtaPcmSide) {
      return (mgg > config.etaPcmLowSideMassRange.value.at(0) && mgg < config.etaPcmLowSideMassRange.value.at(1)) ||
             (mgg > config.etaPcmHighSideMassRange.value.at(0) && mgg < config.etaPcmHighSideMassRange.value.at(1));
    } else {
      static_assert(AlwaysFalseValue<AssocTypeV>, "assocType not found");
    }
  }

  // analysis /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // generalised correlation functions
  // per collision

  // plain info
  void corrProcessPlain(auto const& collision, auto const& associatedThisEvent,
                        auto const& funcPlain)
  {
    // normal spectra (per event - not per trigger)
    for (auto const& associated : associatedThisEvent) {
      funcPlain(collision, associated);
    }
  }

  // correlation
  void corrProcessCorrelation(auto const& collision, auto const& triggersThisEvent, auto const& associatedThisEvent,
                              auto const& funcCorrelation)
  {
    // correlation combinations
    for (auto const& [trigger, associated] : soa::combinations(soa::CombinationsFullIndexPolicy(triggersThisEvent, associatedThisEvent))) {
      funcCorrelation(collision, trigger, associated);
    }
  }

  // mixing-number pt scaling with power law
  size_t nMixingPtPowerScaling(size_t const n0, double const pt)
  {
    double const scale = 1 + config.nTriggerScaleCoefficient.value * std::pow(pt, config.nTriggerScaleExponent.value);
    return n0 * scale;
  }

  // mixing
  struct ConditionMixAcceptAll {
    constexpr bool operator()(auto const&) const
    {
      return true;
    }
  };
  template <TriggerMemoryMode TriggerMemoryModeV, typename T_AssocCondition = ConditionMixAcceptAll>
  void corrProcessMixing(auto const& collision, auto const& associatedThisEvent,
                         auto const& funcMixing, size_t const nTriggerMixingAt0, T_AssocCondition const& funcConditionAssocOnly = {})
  {
    size_t triggerMemoryBinSizeMin = 0;

    // mixing loops (more dynamic than O2 mixing)
    std::deque<MixingTrigger> const& savedTriggers =
      [&]() -> std::deque<MixingTrigger> const& {
      if constexpr (TriggerMemoryModeV == TriggerMemoryMode::Reco) {
        triggerMemoryBinSizeMin = triggerMemoryReco->getBinSizeMin();
        return triggerMemoryReco->getTriggers(collision.posZ(), collision.nGlobalTracks());
      } else if constexpr (TriggerMemoryModeV == TriggerMemoryMode::True) {
        triggerMemoryBinSizeMin = triggerMemoryTrue->getBinSizeMin();
        return triggerMemoryTrue->getTriggers(collision.posZ(), collision.nChargedInEtaRange());
      } else if constexpr (TriggerMemoryModeV == TriggerMemoryMode::TrueAssocEv) {
        triggerMemoryBinSizeMin = triggerMemoryTrueAssocEv->getBinSizeMin();
        return triggerMemoryTrueAssocEv->getTriggers(collision.posZ(), collision.nChargedInEtaRange());
      } else if constexpr (TriggerMemoryModeV == TriggerMemoryMode::TrueAssocEvRecoPtTrig) {
        triggerMemoryBinSizeMin = triggerMemoryTrueAssocEvRecoPtTrig->getBinSizeMin();
        return triggerMemoryTrueAssocEvRecoPtTrig->getTriggers(collision.posZ(), collision.nChargedInEtaRange());
      } else {
        static_assert(AlwaysFalseValue<TriggerMemoryModeV>, "triggerMemoryMode not found");
      }
    }();
    // require enough mixing statistics
    if (triggerMemoryBinSizeMin < static_cast<size_t>(config.nTriggerBinMinThreshold.value)) {
      return;
    }
    const size_t nSavedTriggers = savedTriggers.size();
    const size_t iStartTrigger = randomIntInInterval(0, config.iTriggerStartFraction.value * nSavedTriggers);
    // associated loop
    for (auto const& associated : associatedThisEvent) {
      if (!funcConditionAssocOnly(associated)) {
        continue;
      }
      // trigger settings
      const size_t nTriggerMixing = nMixingPtPowerScaling(nTriggerMixingAt0, associated.pt());
      const size_t mixUpToTriggerN = std::min(nSavedTriggers, iStartTrigger + nTriggerMixing);
      const float perTriggerWeight = 1. / (mixUpToTriggerN - iStartTrigger); // mixUpToTriggerN <= iStartTrigger caught by break statement of loop
      // trigger loop
      for (size_t i_mixingTrigger = iStartTrigger; i_mixingTrigger < mixUpToTriggerN; i_mixingTrigger++) {
        MixingTrigger const& mixingTrigger = savedTriggers[i_mixingTrigger];
        funcMixing(collision, mixingTrigger, associated, perTriggerWeight);
      }
    }
  }

  // analysis tasks ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // =/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/

  void init(InitContext const&)
  {
    // analysis info
    int64_t const noLaterThanCcdb =
      config.noLaterThanCcdbConfig.value == -1 ? std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() : config.noLaterThanCcdbConfig.value;
    ccdb->setURL(config.urlCcdb.value);
    // enabling object caching (otherwise each call goes to CCDB server)
    ccdb->setCaching(true);
    // ccdb->setLocalObjectValidityChecking();
    // not later than (avoids replacing objects while a train is running)
    ccdb->setCreatedNotAfter(noLaterThanCcdb);

    // init analysis variables

    // check coherent config
    double const doublePrecisionDeltaCheck = 1e-5;
    if (std::abs(config.ptTrigMin.value - configHist.binsPtTrig.value.at(0)) > doublePrecisionDeltaCheck) {
      throw std::runtime_error("ptTrigMin does not match trigger axis");
    }
    if (config.applyScalingPhotonPt.value && (config.parScalingPhotonPtData.value.size() != 6 || config.parScalingPhotonPtMc.value.size() != 6)) {
      throw std::runtime_error("wrong number of mc scaling params for photon");
    }
    if (config.applyScalingPi0Pt.value && (config.parScalingPi0PtData.value.size() != 6 || config.parScalingPi0PtMc.value.size() != 6)) {
      throw std::runtime_error("wrong number of mc scaling params for pi0");
    }
    if (config.applyScalingEtaPt.value && (config.parScalingEtaPtData.value.size() != 6 || config.parScalingEtaPtMc.value.size() != 6)) {
      throw std::runtime_error("wrong number of mc scaling params for eta");
    }

    // histograms from ccdb
    initCcdbHistograms(noLaterThanCcdb);

    // create analysis histograms
    initHistograms();

    // mixing trigger memory
    triggerMemoryReco = std::make_unique<MixingTriggerMemory>(configHist.binsZPvBinning.value, configHist.binsMultBinning.value, config.nTriggerSavedForMixing.value);
    triggerMemoryTrue = std::make_unique<MixingTriggerMemory>(configHist.binsZPvBinningMcTrue.value, configHist.binsMultBinningMcTrue.value, config.nTriggerSavedForMixing.value);
    triggerMemoryTrueAssocEv = std::make_unique<MixingTriggerMemory>(configHist.binsZPvBinningMcTrue.value, configHist.binsMultBinningMcTrue.value, config.nTriggerSavedForMixing.value);
    triggerMemoryTrueAssocEvRecoPtTrig = std::make_unique<MixingTriggerMemory>(configHist.binsZPvBinningMcTrue.value, configHist.binsMultBinningMcTrue.value, config.nTriggerSavedForMixing.value);

    // init EM pcm cut class
    initV0PhotonCut();
  }

  // reco only ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ==========================================================================================================================================================================================

  // correlations /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // standard reco histogram compile-time paths
  static constexpr auto HistPathsRecoStd = []() {
    using enum HistFillMode;
    using enum TruthLevel;
    using enum AssocType;

    return std::to_array<HistPathEntry>({
      {PlainPhiEta,  Reco, Hadron,    "reco/plain/h4_ptTrigPtAssocPhiEta_reco_hadron"},
      {PlainZPvMult, Reco, Hadron,    "reco/plain/h4_ptTrigPtAssocZPvMult_reco_hadron"},
      {Correlation,  Reco, Hadron,    "reco/corr/h6_corr_reco_hadron"},
      {Mixing,       Reco, Hadron,    "reco/corr/h6_mix_reco_hadron"},
      {PlainPhiEta,  Reco, Pipm,      "reco/plain/h4_ptTrigPtAssocPhiEta_reco_pipm"},
      {PlainZPvMult, Reco, Pipm,      "reco/plain/h4_ptTrigPtAssocZPvMult_reco_pipm"},
      {Correlation,  Reco, Pipm,      "reco/corr/h6_corr_reco_pipm"},
      {Mixing,       Reco, Pipm,      "reco/corr/h6_mix_reco_pipm"},
      {PlainPhiEta,  Reco, PhotonPcm, "reco/plain/h4_ptTrigPtAssocPhiEta_reco_photonPcm"},
      {PlainZPvMult, Reco, PhotonPcm, "reco/plain/h4_ptTrigPtAssocZPvMult_reco_photonPcm"},
      {Correlation,  Reco, PhotonPcm, "reco/corr/h6_corr_reco_photonPcm"},
      {Mixing,       Reco, PhotonPcm, "reco/corr/h6_mix_reco_photonPcm"}
    });
  }();

  template <AssocType AssocTypeV>
  void fillProcessRecoStd(auto const& collisions, auto const& triggers, auto const& assocs)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      // group collision ++
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const [assocsThisEvent, nMixingAt0] = [&]() {
        if constexpr (AssocTypeV == AssocType::Hadron) {
          return std::tuple{assocs.sliceBy(perColHadrons, collision.globalIndex()), config.nMixingAt0Hadron.value};
        } else if constexpr (AssocTypeV == AssocType::Pipm) {
          return std::tuple{assocs.sliceBy(perColPipms, collision.globalIndex()), config.nMixingAt0Pipm.value};
        } else if constexpr (AssocTypeV == AssocType::PhotonPcm) {
          return std::tuple{assocs.sliceBy(perColPhotonPcms, collision.globalIndex()), config.nMixingAt0PhotonPcm.value};
        } else {
          static_assert(AlwaysFalseValue<AssocTypeV>, "assocType not supported");
        }
      }();

      auto const funcPlain = [this](auto const& collision, auto const& associated) {
        // reco assoc selection
        if (!isSelectedAssoc<AssocTypeV>(associated)) {
          return;
        }

        HistFillInfoPlainPhiEtaZPvMult histFillInfo = {collision.ptMax(), associated.pt(), associated.phi(), associated.eta(), collision.posZ(), static_cast<double>(collision.nGlobalTracks())};
        histFillInfo.weight = getMcAssocPtScaling<AssocTypeV>(associated) * getMcK0ShortPtScaling<AssocTypeV>(associated);
        fillHistMain<HistPathsRecoStd, HistFillMode::Plain, TruthLevel::Reco, AssocTypeV>(histFillInfo);
      };
      corrProcessPlain(collision, assocsThisEvent, funcPlain);

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        // reco assoc selection
        if (!isSelectedAssoc<AssocTypeV>(associated)) {
          return;
        }
        // exclude self correlation
        if (isRecoSelfTrigger(trigger, associated)) {
          return;
        }

        HistFillInfoCorr histFillInfo = {getDeltaPhi(trigger.phi(), associated.phi()), trigger.eta() - associated.eta(),
                                         trigger.pt(), associated.pt(), collision.posZ(), static_cast<double>(collision.nGlobalTracks())};
        histFillInfo.weight = getMcAssocPtScaling<AssocTypeV>(associated) * getMcK0ShortPtScaling<AssocTypeV>(associated) * getTriggerPtCorrection(trigger.pt());
        fillHistMain<HistPathsRecoStd, HistFillMode::Correlation, TruthLevel::Reco, AssocTypeV>(histFillInfo);
      };
      corrProcessCorrelation(collision, triggersThisEvent, assocsThisEvent, funcCorrelation);

      // select mixing events
      if (config.doTrigEvMixing.value && !collision.trigEv()) {
        continue;
      }

      auto const conditionMixAssocOnly = [this](auto const& associated) {
        // reco assoc selection
        return isSelectedAssoc<AssocTypeV>(associated);
      };
      auto const funcMixing = [this](auto const& collision, auto const& mixingTrigger, auto const& associated, auto const perTriggerWeight) {
        HistFillInfoCorr histFillInfo = {getDeltaPhi(mixingTrigger.phi(), associated.phi()), mixingTrigger.eta() - associated.eta(),
                                         mixingTrigger.pt(), associated.pt(), collision.posZ(), static_cast<double>(collision.nGlobalTracks())};
        histFillInfo.weight = perTriggerWeight * getMcAssocPtScaling<AssocTypeV>(associated) * getMcK0ShortPtScaling<AssocTypeV>(associated) * getTriggerPtCorrection(mixingTrigger.pt());
        fillHistMain<HistPathsRecoStd, HistFillMode::Mixing, TruthLevel::Reco, AssocTypeV>(histFillInfo);
      };
      corrProcessMixing<TriggerMemoryMode::Reco>(collision, assocsThisEvent, funcMixing, nMixingAt0, conditionMixAssocOnly);
    }
  }

  void processRecoHadron(CorrCollisions const& collisions, aod::Triggers const& triggers, aod::Hadrons const& hadrons)
  {
    fillProcessRecoStd<AssocType::Hadron>(collisions, triggers, hadrons);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processRecoHadron, "process standard for associated hardon", false);

  void processRecoPipm(CorrCollisions const& collisions, aod::Triggers const& triggers, aod::Pipms const& pipms)
  {
    fillProcessRecoStd<AssocType::Pipm>(collisions, triggers, pipms);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processRecoPipm, "process standard for associated pipm", false);

  void processRecoPhotonPcm(CorrCollisions const& collisions, aod::Triggers const& triggers, PhotonPcmsWithV0 const& photonPcms,
                            aod::V0PhotonsKF const&, aod::V0Legs const&)
  {
    fillProcessRecoStd<AssocType::PhotonPcm>(collisions, triggers, photonPcms);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processRecoPhotonPcm, "process standard for associated photonPcm", false);

  // photonPcmPair reco histogram compile-time paths
  static constexpr auto HistPathsRecoPhotonPcmPair = []() {
    using enum HistFillMode;
    using enum TruthLevel;
    using enum AssocType;

    return std::to_array<HistPathEntry>({
      {PlainPt,      Reco, Pi0PcmPeak, "reco/plain/h2_ptTrigPtAssoc_reco_pi0PcmPeak"},
      {Correlation,  Reco, Pi0PcmPeak, "reco/corr/h6_corr_reco_pi0PcmPeak"},
      {Mixing,       Reco, Pi0PcmPeak, "reco/corr/h6_mix_reco_pi0PcmPeak"},
      {PlainPt,      Reco, Pi0PcmSide, "reco/plain/h2_ptTrigPtAssoc_reco_pi0PcmSide"},
      {Correlation,  Reco, Pi0PcmSide, "reco/corr/h6_corr_reco_pi0PcmSide"},
      {Mixing,       Reco, Pi0PcmSide, "reco/corr/h6_mix_reco_pi0PcmSide"},
      {PlainPt,      Reco, EtaPcmPeak, "reco/plain/h2_ptTrigPtAssoc_reco_etaPcmPeak"},
      {Correlation,  Reco, EtaPcmPeak, "reco/corr/h6_corr_reco_etaPcmPeak"},
      {Mixing,       Reco, EtaPcmPeak, "reco/corr/h6_mix_reco_etaPcmPeak"},
      {PlainPt,      Reco, EtaPcmSide, "reco/plain/h2_ptTrigPtAssoc_reco_etaPcmSide"},
      {Correlation,  Reco, EtaPcmSide, "reco/corr/h6_corr_reco_etaPcmSide"},
      {Mixing,       Reco, EtaPcmSide, "reco/corr/h6_mix_reco_etaPcmSide"}
    });
  }();

  template <HistFillMode HistFillModeV>
  void fillAssocsRecoPhotonPcmPair(auto const& histFillInfo, double const mgg) {
    // pi0
    if (isInMassRange<AssocType::Pi0PcmPeak>(mgg)) {
      fillHistMain<HistPathsRecoPhotonPcmPair, HistFillModeV, TruthLevel::Reco, AssocType::Pi0PcmPeak>(histFillInfo);
      return;
    }
    if (isInMassRange<AssocType::Pi0PcmSide>(mgg)) {
      fillHistMain<HistPathsRecoPhotonPcmPair, HistFillModeV, TruthLevel::Reco, AssocType::Pi0PcmSide>(histFillInfo);
      return;
    }
    // eta
    if (isInMassRange<AssocType::EtaPcmPeak>(mgg)) {
      fillHistMain<HistPathsRecoPhotonPcmPair, HistFillModeV, TruthLevel::Reco, AssocType::EtaPcmPeak>(histFillInfo);
      return;
    }
    if (isInMassRange<AssocType::EtaPcmSide>(mgg)) {
      fillHistMain<HistPathsRecoPhotonPcmPair, HistFillModeV, TruthLevel::Reco, AssocType::EtaPcmSide>(histFillInfo);
      return;
    }
  }

  void processRecoPhotonPcmPair(CorrCollisions const& collisions, aod::Triggers const& triggers, PhotonPcmPairsWithV0 const& photonPcmPairs,
                                aod::V0PhotonsKF const&, aod::V0Legs const&)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const photonPcmPairsThisEvent = photonPcmPairs.sliceBy(perColPhotonPcmPairs, collision.globalIndex());

      auto const funcPlain = [this](auto const& collision, auto const& associated) {
        // reco assoc selection
        if (!isSelectedAssoc<AssocType::PhotonPcmPair>(associated)) {
          return;
        }

        double const weight = getMcAssocPtScaling<AssocType::PhotonPcmPair>(associated) * getMcK0ShortPtScaling<AssocType::PhotonPcmPair>(associated);

        // mgg distribution
        histos.fill(HIST("reco/plain/h5_ptTrigPtAssocMggZPvMult_reco_photonPcmPair"),
                    collision.ptMax(), associated.pt(), associated.mgg(), collision.posZ(), collision.nGlobalTracks(),
                    weight);

        HistFillInfoPlainPt histFillInfo = {collision.ptMax(), associated.pt()};
        histFillInfo.weight = weight;
        fillAssocsRecoPhotonPcmPair<HistFillMode::Plain>(histFillInfo, associated.mgg());
      };
      corrProcessPlain(collision, photonPcmPairsThisEvent, funcPlain);

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        // reco assoc selection
        if (!isSelectedAssoc<AssocType::PhotonPcmPair>(associated)) {
          return;
        }
        // exclude self correlation
        if (isRecoSelfTrigger(trigger, associated)) {
          return;
        }

        double const weight = getMcAssocPtScaling<AssocType::PhotonPcmPair>(associated) * getMcK0ShortPtScaling<AssocType::PhotonPcmPair>(associated) * getTriggerPtCorrection(trigger.pt());

        // mgg distribution
        histos.fill(HIST("reco/corr/h5_ptTrigPtAssocMggZPvMult_reco_photonPcmPair"),
                    trigger.pt(), associated.pt(), associated.mgg(), collision.posZ(), collision.nGlobalTracks(),
                    weight);

        HistFillInfoCorr histFillInfo = {getDeltaPhi(trigger.phi(), associated.phi()), trigger.eta() - associated.eta(),
                                         trigger.pt(), associated.pt(), collision.posZ(), static_cast<double>(collision.nGlobalTracks())};
        histFillInfo.weight = weight;
        fillAssocsRecoPhotonPcmPair<HistFillMode::Correlation>(histFillInfo, associated.mgg());
      };
      corrProcessCorrelation(collision, triggersThisEvent, photonPcmPairsThisEvent, funcCorrelation);

      // select mixing events
      if (config.doTrigEvMixing.value && !collision.trigEv()) {
        continue;
      }

      auto const conditionMixAssocOnly = [this](auto const& associated) {
        // reco assoc selection
        return isSelectedAssoc<AssocType::PhotonPcmPair>(associated);
      };
      auto const funcMixing = [this](auto const& collision, auto const& mixingTrigger, auto const& associated, auto const perTriggerWeight) {
        HistFillInfoCorr histFillInfo = {getDeltaPhi(mixingTrigger.phi(), associated.phi()), mixingTrigger.eta() - associated.eta(),
                                         mixingTrigger.pt(), associated.pt(), collision.posZ(), static_cast<double>(collision.nGlobalTracks())};
        histFillInfo.weight = perTriggerWeight * getMcAssocPtScaling<AssocType::PhotonPcmPair>(associated) * getMcK0ShortPtScaling<AssocType::PhotonPcmPair>(associated) * getTriggerPtCorrection(mixingTrigger.pt());
        fillAssocsRecoPhotonPcmPair<HistFillMode::Mixing>(histFillInfo, associated.mgg());
      };
      corrProcessMixing<TriggerMemoryMode::Reco>(collision, photonPcmPairsThisEvent, funcMixing, config.nMixingAt0PhotonPcmPair.value, conditionMixAssocOnly);
    }
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processRecoPhotonPcmPair, "process standard for associated photonPcmPair", false);

  void processPhotonPcmPairMix(CorrCollisions const& collisions, PhotonPcmsWithV0 const& photonPcms,
                               aod::V0PhotonsKF const&, aod::V0Legs const&)
  {
    BinningZPvMult const binningZPvMult{{prependVariableWidth(configHist.binsZPvBinning.value), prependVariableWidth(configHist.binsMultBinning.value)}, true};
    auto photonPcmsTuple = std::make_tuple(photonPcms);
    SameKindPair<CorrCollisions, PhotonPcmsWithV0, BinningZPvMult> pairs{binningZPvMult, config.nNeighboursMixingPhotonPcmPair.value, -1, collisions, photonPcmsTuple, &cache};

    // mixed events
    for (auto pair = pairs.begin(); pair != pairs.end(); pair++) {
      auto const& [collision1, photonPcms1, collision2, photonPcms2] = *pair;
      // event selection
      if (!totalEvSel(collision1) || !totalEvSel(collision2)) {
        continue;
      }
      // mixing loop
      for (auto const& [photonPcm1, photonPcm2] : soa::combinations(soa::CombinationsFullIndexPolicy(photonPcms1, photonPcms2))) {
        ROOT::Math::PtEtaPhiMVector const p4photonPcm1(photonPcm1.pt(), photonPcm1.eta(), photonPcm1.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector const p4photonPcm2(photonPcm2.pt(), photonPcm2.eta(), photonPcm2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector const p4photonPcmPair = p4photonPcm1 + p4photonPcm2;

        if (std::abs(p4photonPcmPair.eta()) > config.absEtaMax.value) {
          continue;
        }
        // reco assoc selection
        if (!isSelectedAssoc<AssocType::PhotonPcm>(photonPcm1) || !isSelectedAssoc<AssocType::PhotonPcm>(photonPcm2)) {
          continue;
        }
        if (!checkEnergyAsymDiPhoton(photonPcm1.v0PhotonKF(), photonPcm2.v0PhotonKF())) {
          continue;
        }

        histos.fill(HIST("reco/plain/h4_ptMggZPvMult_photonPcmPair_evMix"),
                    p4photonPcmPair.pt(), p4photonPcmPair.M(), collision1.posZ(), collision1.nGlobalTracks());
      }
    }
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPhotonPcmPairMix, "process gamma-gamma mixing for photonPcm", false);

  // other ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  void processPcmDcaZ(CorrCollisions::iterator const& collision, PhotonPcmsWithV0 const& photonPcms,
                      aod::V0PhotonsKF const&, aod::V0Legs const&)
  {
    // event selection
    if (!totalEvSel(collision)) {
      return;
    }

    for (auto const& photonPcm : photonPcms) {
      // reco assoc selection
      if (!isSelectedAssoc<AssocType::PhotonPcm>(photonPcm)) {
        continue;
      }

      double const weight = getMcAssocPtScaling<AssocType::PhotonPcm>(photonPcm) * getMcK0ShortPtScaling<AssocType::PhotonPcm>(photonPcm);

      histos.fill(HIST("reco/plain/h5_ptTrigPtAssocDcaZZPvMult_photonPcm"),
                  collision.ptMax(), photonPcm.pt(), photonPcm.v0PhotonKF().dcaZtopv(), collision.posZ(), collision.nGlobalTracks(),
                  weight);
    }
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPcmDcaZ, "process to test dcaZ distribution of photons", false);

  // MC ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ==========================================================================================================================================================================================

  // correlations /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // mc true++ histogram compile-time paths
  static constexpr auto HistPathsMcTrue = []() {
    using enum HistFillMode;
    using enum TruthLevel;
    using enum AssocType;

    return std::to_array<HistPathEntry>({
      {PlainPhiEta,  True, Hadron,            "mc/plain/h4_ptTrigPtAssocPhiEta_true_hadron"},
      {PlainZPvMult, True, Hadron,            "mc/plain/h4_ptTrigPtAssocZPvMult_true_hadron"},
      {Correlation,  True, Hadron,            "mc/corr/h6_corr_true_hadron"},
      {Mixing,       True, Hadron,            "mc/corr/h6_mix_true_hadron"},
      {PlainPhiEta,  True, Pipm,              "mc/plain/h4_ptTrigPtAssocPhiEta_true_pipm"},
      {PlainZPvMult, True, Pipm,              "mc/plain/h4_ptTrigPtAssocZPvMult_true_pipm"},
      {Correlation,  True, Pipm,              "mc/corr/h6_corr_true_pipm"},
      {Mixing,       True, Pipm,              "mc/corr/h6_mix_true_pipm"},
      {PlainPhiEta,  True, PhotonPcm,         "mc/plain/h4_ptTrigPtAssocPhiEta_true_photon"},
      {PlainZPvMult, True, PhotonPcm,         "mc/plain/h4_ptTrigPtAssocZPvMult_true_photon"},
      {Correlation,  True, PhotonPcm,         "mc/corr/h6_corr_true_photon"},
      {Mixing,       True, PhotonPcm,         "mc/corr/h6_mix_true_photon"},
      {PlainPt,      True, PhotonDecay,       "mc/plain/h2_ptTrigPtAssoc_true_photonDecay"},
      {Correlation,  True, PhotonDecay,       "mc/corr/h6_corr_true_photonDecay"},
      {Mixing,       True, PhotonDecay,       "mc/corr/h6_mix_true_photonDecay"},
      {PlainPt,      True, PhotonDirect,      "mc/plain/h2_ptTrigPtAssoc_true_photonDirect"},
      {Correlation,  True, PhotonDirect,      "mc/corr/h6_corr_true_photonDirect"},
      {Mixing,       True, PhotonDirect,      "mc/corr/h6_mix_true_photonDirect"},
      {PlainPt,      True, Pi0Pcm,            "mc/plain/h2_ptTrigPtAssoc_true_pi0"},
      {Correlation,  True, Pi0Pcm,            "mc/corr/h6_corr_true_pi0"},
      {Mixing,       True, Pi0Pcm,            "mc/corr/h6_mix_true_pi0"},
      {PlainPt,      True, PhotonPi0,         "mc/plain/h2_ptTrigPtAssoc_true_photonPi0"},
      {Correlation,  True, PhotonPi0,         "mc/corr/h6_corr_true_photonPi0"},
      {Mixing,       True, PhotonPi0,         "mc/corr/h6_mix_true_photonPi0"},
      {PlainPt,      True, EtaPcm,            "mc/plain/h2_ptTrigPtAssoc_true_eta"},
      {Correlation,  True, EtaPcm,            "mc/corr/h6_corr_true_eta"},
      {Mixing,       True, EtaPcm,            "mc/corr/h6_mix_true_eta"},
      {PlainPt,      True, PhotonEta,         "mc/plain/h2_ptTrigPtAssoc_true_photonEta"},
      {Correlation,  True, PhotonEta,         "mc/corr/h6_corr_true_photonEta"},
      {Mixing,       True, PhotonEta,         "mc/corr/h6_mix_true_photonEta"},
      {PlainPt,      True, Omega,             "mc/plain/h2_ptTrigPtAssoc_true_omega"},
      {Correlation,  True, Omega,             "mc/corr/h6_corr_true_omega"},
      {Mixing,       True, Omega,             "mc/corr/h6_mix_true_omega"},
      {PlainPt,      True, PhotonOmega,       "mc/plain/h2_ptTrigPtAssoc_true_photonOmega"},
      {Correlation,  True, PhotonOmega,       "mc/corr/h6_corr_true_photonOmega"},
      {Mixing,       True, PhotonOmega,       "mc/corr/h6_mix_true_photonOmega"},
      {PlainPt,      True, EtaPrime,          "mc/plain/h2_ptTrigPtAssoc_true_etaPrime"},
      {Correlation,  True, EtaPrime,          "mc/corr/h6_corr_true_etaPrime"},
      {Mixing,       True, EtaPrime,          "mc/corr/h6_mix_true_etaPrime"},
      {PlainPt,      True, PhotonEtaPrime,    "mc/plain/h2_ptTrigPtAssoc_true_photonEtaPrime"},
      {Correlation,  True, PhotonEtaPrime,    "mc/corr/h6_corr_true_photonEtaPrime"},
      {Mixing,       True, PhotonEtaPrime,    "mc/corr/h6_mix_true_photonEtaPrime"},
      {PlainPt,      True, PhotonOtherMother, "mc/plain/h2_ptTrigPtAssoc_true_photonOtherMother"},
      {Correlation,  True, PhotonOtherMother, "mc/corr/h6_corr_true_photonOtherMother"},
      {Mixing,       True, PhotonOtherMother, "mc/corr/h6_mix_true_photonOtherMother"},

      {PlainPhiEta,  TrueAssocEv, Hadron,    "mc/plain/h4_ptTrigPtAssocPhiEta_trueAssocEv_hadron"},
      {PlainZPvMult, TrueAssocEv, Hadron,    "mc/plain/h4_ptTrigPtAssocZPvMult_trueAssocEv_hadron"},
      {Correlation,  TrueAssocEv, Hadron,    "mc/corr/h6_corr_trueAssocEv_hadron"},
      {Mixing,       TrueAssocEv, Hadron,    "mc/corr/h6_mix_trueAssocEv_hadron"},
      {PlainPhiEta,  TrueAssocEv, Pipm,      "mc/plain/h4_ptTrigPtAssocPhiEta_trueAssocEv_pipm"},
      {PlainZPvMult, TrueAssocEv, Pipm,      "mc/plain/h4_ptTrigPtAssocZPvMult_trueAssocEv_pipm"},
      {Correlation,  TrueAssocEv, Pipm,      "mc/corr/h6_corr_trueAssocEv_pipm"},
      {Mixing,       TrueAssocEv, Pipm,      "mc/corr/h6_mix_trueAssocEv_pipm"},
      {PlainPhiEta,  TrueAssocEv, PhotonPcm, "mc/plain/h4_ptTrigPtAssocPhiEta_trueAssocEv_photon"},
      {PlainZPvMult, TrueAssocEv, PhotonPcm, "mc/plain/h4_ptTrigPtAssocZPvMult_trueAssocEv_photon"},
      {Correlation,  TrueAssocEv, PhotonPcm, "mc/corr/h6_corr_trueAssocEv_photon"},
      {Mixing,       TrueAssocEv, PhotonPcm, "mc/corr/h6_mix_trueAssocEv_photon"},
      {PlainPt,      TrueAssocEv, Pi0Pcm,    "mc/plain/h2_ptTrigPtAssoc_trueAssocEv_pi0"},
      {Correlation,  TrueAssocEv, Pi0Pcm,    "mc/corr/h6_corr_trueAssocEv_pi0"},
      {Mixing,       TrueAssocEv, Pi0Pcm,    "mc/corr/h6_mix_trueAssocEv_pi0"},
      {PlainPt,      TrueAssocEv, EtaPcm,    "mc/plain/h2_ptTrigPtAssoc_trueAssocEv_eta"},
      {Correlation,  TrueAssocEv, EtaPcm,    "mc/corr/h6_corr_trueAssocEv_eta"},
      {Mixing,       TrueAssocEv, EtaPcm,    "mc/corr/h6_mix_trueAssocEv_eta"},

      {PlainPhiEta,  TrueAssocEvRecoPtTrig, Hadron,    "mc/plain/h4_ptTrigPtAssocPhiEta_trueAssocEvRecoPtTrig_hadron"},
      {PlainZPvMult, TrueAssocEvRecoPtTrig, Hadron,    "mc/plain/h4_ptTrigPtAssocZPvMult_trueAssocEvRecoPtTrig_hadron"},
      {Correlation,  TrueAssocEvRecoPtTrig, Hadron,    "mc/corr/h6_corr_trueAssocEvRecoPtTrig_hadron"},
      {Mixing,       TrueAssocEvRecoPtTrig, Hadron,    "mc/corr/h6_mix_trueAssocEvRecoPtTrig_hadron"},
      {PlainPhiEta,  TrueAssocEvRecoPtTrig, Pipm,      "mc/plain/h4_ptTrigPtAssocPhiEta_trueAssocEvRecoPtTrig_pipm"},
      {PlainZPvMult, TrueAssocEvRecoPtTrig, Pipm,      "mc/plain/h4_ptTrigPtAssocZPvMult_trueAssocEvRecoPtTrig_pipm"},
      {Correlation,  TrueAssocEvRecoPtTrig, Pipm,      "mc/corr/h6_corr_trueAssocEvRecoPtTrig_pipm"},
      {Mixing,       TrueAssocEvRecoPtTrig, Pipm,      "mc/corr/h6_mix_trueAssocEvRecoPtTrig_pipm"},
      {PlainPhiEta,  TrueAssocEvRecoPtTrig, PhotonPcm, "mc/plain/h4_ptTrigPtAssocPhiEta_trueAssocEvRecoPtTrig_photon"},
      {PlainZPvMult, TrueAssocEvRecoPtTrig, PhotonPcm, "mc/plain/h4_ptTrigPtAssocZPvMult_trueAssocEvRecoPtTrig_photon"},
      {Correlation,  TrueAssocEvRecoPtTrig, PhotonPcm, "mc/corr/h6_corr_trueAssocEvRecoPtTrig_photon"},
      {Mixing,       TrueAssocEvRecoPtTrig, PhotonPcm, "mc/corr/h6_mix_trueAssocEvRecoPtTrig_photon"},
      {PlainPt,      TrueAssocEvRecoPtTrig, Pi0Pcm,    "mc/plain/h2_ptTrigPtAssoc_trueAssocEvRecoPtTrig_pi0"},
      {Correlation,  TrueAssocEvRecoPtTrig, Pi0Pcm,    "mc/corr/h6_corr_trueAssocEvRecoPtTrig_pi0"},
      {Mixing,       TrueAssocEvRecoPtTrig, Pi0Pcm,    "mc/corr/h6_mix_trueAssocEvRecoPtTrig_pi0"},
      {PlainPt,      TrueAssocEvRecoPtTrig, EtaPcm,    "mc/plain/h2_ptTrigPtAssoc_trueAssocEvRecoPtTrig_eta"},
      {Correlation,  TrueAssocEvRecoPtTrig, EtaPcm,    "mc/corr/h6_corr_trueAssocEvRecoPtTrig_eta"},
      {Mixing,       TrueAssocEvRecoPtTrig, EtaPcm,    "mc/corr/h6_mix_trueAssocEvRecoPtTrig_eta"},

      {PlainPt,      TrueDecaySelfTrig, Pi0Pcm, "mc/plain/h2_ptTrigPtAssoc_trueDecaySelfTrig_pi0"},
      {Correlation,  TrueDecaySelfTrig, Pi0Pcm, "mc/corr/h6_corr_trueDecaySelfTrig_pi0"},
      {Mixing,       TrueDecaySelfTrig, Pi0Pcm, "mc/corr/h6_mix_trueDecaySelfTrig_pi0"},
      {PlainPt,      TrueDecaySelfTrig, EtaPcm, "mc/plain/h2_ptTrigPtAssoc_trueDecaySelfTrig_eta"},
      {Correlation,  TrueDecaySelfTrig, EtaPcm, "mc/corr/h6_corr_trueDecaySelfTrig_eta"},
      {Mixing,       TrueDecaySelfTrig, EtaPcm, "mc/corr/h6_mix_trueDecaySelfTrig_eta"},

      {PlainPt,      TrueMeasDecay, Pi0Pcm, "mc/plain/h2_ptTrigPtAssoc_trueMeasDecay_pi0"},
      {Correlation,  TrueMeasDecay, Pi0Pcm, "mc/corr/h6_corr_trueMeasDecay_pi0"},
      {Mixing,       TrueMeasDecay, Pi0Pcm, "mc/corr/h6_mix_trueMeasDecay_pi0"},
      {PlainPt,      TrueMeasDecay, EtaPcm, "mc/plain/h2_ptTrigPtAssoc_trueMeasDecay_eta"},
      {Correlation,  TrueMeasDecay, EtaPcm, "mc/corr/h6_corr_trueMeasDecay_eta"},
      {Mixing,       TrueMeasDecay, EtaPcm, "mc/corr/h6_mix_trueMeasDecay_eta"},

      {PlainPt,      TrueGeoAcc, Pi0Pcm, "mc/plain/h2_ptTrigPtAssoc_trueGeoAcc_pi0"},
      {Correlation,  TrueGeoAcc, Pi0Pcm, "mc/corr/h6_corr_trueGeoAcc_pi0"},
      {Mixing,       TrueGeoAcc, Pi0Pcm, "mc/corr/h6_mix_trueGeoAcc_pi0"},
      {PlainPt,      TrueGeoAcc, EtaPcm, "mc/plain/h2_ptTrigPtAssoc_trueGeoAcc_eta"},
      {Correlation,  TrueGeoAcc, EtaPcm, "mc/corr/h6_corr_trueGeoAcc_eta"},
      {Mixing,       TrueGeoAcc, EtaPcm, "mc/corr/h6_mix_trueGeoAcc_eta"}
    });
  }();

  template <HistFillMode HistFillModeV, TruthLevel TruthLevelV, typename T_DecaySelfTrigFunc = decltype([](){})>
  void fillAssocsMcTrue(auto histFillInfo, auto const& assocParticle, T_DecaySelfTrigFunc const& isDecaySelfTrig = {}) {
    // main eta cut
    if (std::abs(assocParticle.eta()) > config.absEtaMax.value) {
      return;
    }

    // standard particles (marked physical primary)
    if (assocParticle.isPhysicalPrimary()) {
      // hadron (charged primary)
      if (isCharged(assocParticle)) {
        fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevelV, AssocType::Hadron>(histFillInfo);
      }
      // pipm
      if (std::abs(assocParticle.pdgCode()) == PDG_t::kPiPlus) {
        fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevelV, AssocType::Pipm>(histFillInfo);
        return;
      }
      // photon
      if (assocParticle.pdgCode() == PDG_t::kGamma) {
        // mc-to-data scaling
        histFillInfo.weight *= getMcAssocPtScalingFromPt<AssocType::PhotonPcm>(assocParticle.pt());

        fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevelV, AssocType::PhotonPcm>(histFillInfo);

        // extra info for decay correlation only for total true level
        if constexpr (TruthLevelV == TruthLevel::True) {
          // decay and direct
          int const statusDecayLow = 91;
          int const statusDecayUp = 99;
          int const statusDirect = 62;

          if (std::abs(assocParticle.getGenStatusCode()) >= statusDecayLow && std::abs(assocParticle.getGenStatusCode()) <= statusDecayUp) {
            fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevelV, AssocType::PhotonDecay>(histFillInfo);
            // decays from different mothers
            int const pdgMother = assocParticle.template mothers_as<aod::JetParticles>().begin()->pdgCode();
            switch (pdgMother) {
              case PDG_t::kPi0:
                fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevelV, AssocType::PhotonPi0>(histFillInfo);
                break;
              case constants::physics::Pdg::kEta:
                fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevelV, AssocType::PhotonEta>(histFillInfo);
                break;
              case constants::physics::Pdg::kOmega:
                fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevelV, AssocType::PhotonOmega>(histFillInfo);
                break;
              case constants::physics::Pdg::kEtaPrime:
                fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevelV, AssocType::PhotonEtaPrime>(histFillInfo);
                break;
              default:
                fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevelV, AssocType::PhotonOtherMother>(histFillInfo);
                break;
            }
          } else {
            fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevelV, AssocType::PhotonDirect>(histFillInfo);
            // log unexpected status
            if (assocParticle.getGenStatusCode() != statusDirect) {
              LOGF(info, "filled direct photon with status: %i", assocParticle.getGenStatusCode());
            }
          }
        }
        return;
      }
      return;
    }
    // decaying particles (not marked physical primary)
    if (!isDecayPrimary(assocParticle)) {
      return;
    }
    // pi0
    if (assocParticle.pdgCode() == PDG_t::kPi0) {
      // mc-to-data scaling
      histFillInfo.weight *= getMcAssocPtScalingFromPt<AssocType::Pi0Pcm>(assocParticle.pt());

      if constexpr (TruthLevelV == TruthLevel::True) {
        fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevel::TrueDecaySelfTrig, AssocType::Pi0Pcm>(histFillInfo);
      }
      // self trigger
      if constexpr (HistFillModeV == HistFillMode::Plain) {
        if (isDecaySelfTrig(AssocType::Pi0Pcm)) {
          return;
        }
      } else if constexpr (HistFillModeV == HistFillMode::Correlation) {
        if (isDecaySelfTrig()) {
          return;
        }
      } else if constexpr (HistFillModeV == HistFillMode::Mixing) {
      } else {
        static_assert(AlwaysFalseValue<HistFillModeV>, "histFillMode not found");
      }
      fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevelV, AssocType::Pi0Pcm>(histFillInfo);
      // decaying particles extra info
      if constexpr (TruthLevelV == TruthLevel::TrueAssocEvRecoPtTrig) {
        // chosen decay
        if (!isDecayToGG(assocParticle)) {
          return;
        }
        fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevel::TrueMeasDecay, AssocType::Pi0Pcm>(histFillInfo);
        // daughters in acceptance
        if (!checkDaughtersInAcceptance(assocParticle, {-config.absEtaMax.value, config.absEtaMax.value})) {
          return;
        }
        fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevel::TrueGeoAcc, AssocType::Pi0Pcm>(histFillInfo);
      }
      return;
    }
    // eta
    if (assocParticle.pdgCode() == constants::physics::Pdg::kEta) {
      // mc-to-data scaling
      histFillInfo.weight *= getMcAssocPtScalingFromPt<AssocType::EtaPcm>(assocParticle.pt());

      if constexpr (TruthLevelV == TruthLevel::True) {
        fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevel::TrueDecaySelfTrig, AssocType::EtaPcm>(histFillInfo);
      }
      // self trigger
      if constexpr (HistFillModeV == HistFillMode::Plain) {
        if (isDecaySelfTrig(AssocType::EtaPcm)) {
          return;
        }
      } else if constexpr (HistFillModeV == HistFillMode::Correlation) {
        if (isDecaySelfTrig()) {
          return;
        }
      } else if constexpr (HistFillModeV == HistFillMode::Mixing) {
      } else {
        static_assert(AlwaysFalseValue<HistFillModeV>, "histFillMode not found");
      }
      fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevelV, AssocType::EtaPcm>(histFillInfo);
      // decaying particles extra info
      if constexpr (TruthLevelV == TruthLevel::TrueAssocEvRecoPtTrig) {
        // chosen decay
        if (!isDecayToGG(assocParticle)) {
          return;
        }
        fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevel::TrueMeasDecay, AssocType::EtaPcm>(histFillInfo);
        // daughters in acceptance
        if (!checkDaughtersInAcceptance(assocParticle, {-config.absEtaMax.value, config.absEtaMax.value})) {
          return;
        }
        fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevel::TrueGeoAcc, AssocType::EtaPcm>(histFillInfo);
      }
      return;
    }

    // extra info for decay correlation only for total true level
    if constexpr (TruthLevelV == TruthLevel::True) {
      // omega
      if (assocParticle.pdgCode() == constants::physics::Pdg::kOmega) {
        fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevelV, AssocType::Omega>(histFillInfo);
        return;
      }
      // etaPrime
      if (assocParticle.pdgCode() == constants::physics::Pdg::kEtaPrime) {
        fillHistMain<HistPathsMcTrue, HistFillModeV, TruthLevelV, AssocType::EtaPrime>(histFillInfo);
        return;
      }
    }
  }

  void processMcTrue(CorrMcCollisions const& mcCollisions, aod::TrueTriggers const& trueTriggers, aod::JetParticles const& mcParticles)
  {
    for (auto const& mcCollision : mcCollisions) {
      // event selection
      if (!totalEvSel(mcCollision)) {
        continue;
      }

      // group collision
      auto const trueTriggersThisEvent = trueTriggers.sliceBy(perColTrueTriggers, mcCollision.globalIndex());
      auto const mcParticlesThisEvent = mcParticles.sliceBy(perColMcParticles, mcCollision.globalIndex());

      // self correlation of events (plain)
      std::map<AssocType, bool> const decaySelfTrigPlainMap = getDecaySelfTrigPlainMap<TruthLevel::True>(trueTriggersThisEvent);

      auto const funcPlain = [this, &decaySelfTrigPlainMap](auto const& mcCollision, auto const& assocParticle) {
        // self trigger info for decaying particles
        auto const isDecaySelfTrig = [&decaySelfTrigPlainMap](AssocType const assocType) {
          return decaySelfTrigPlainMap.at(assocType);
        };
        HistFillInfoPlainPhiEtaZPvMult histFillInfo = {mcCollision.ptMax(), assocParticle.pt(), assocParticle.phi(), assocParticle.eta(), mcCollision.posZ(), static_cast<double>(mcCollision.nChargedInEtaRange())};
        fillAssocsMcTrue<HistFillMode::Plain, TruthLevel::True>(histFillInfo, assocParticle, isDecaySelfTrig);
      };
      corrProcessPlain(mcCollision, mcParticlesThisEvent, funcPlain);

      auto const funcCorrelation = [this](auto const& mcCollision, auto const& trueTrigger, auto const& assocParticle) {
        auto const& triggerParticle = trueTrigger.jetMcParticle();
        // self trigger info for decaying particles
        auto const isDecaySelfTrig = [this, &assocParticle, &triggerParticle]() {
          return producesTrigger(assocParticle, triggerParticle);
        };
        // exclude self correlation
        if (trueTrigger.jetMcParticleId() == assocParticle.globalIndex()) {
          return;
        }

        HistFillInfoCorr histFillInfo = {getDeltaPhi(trueTrigger.phi(), assocParticle.phi()), trueTrigger.eta() - assocParticle.eta(),
                                         trueTrigger.pt(), assocParticle.pt(), mcCollision.posZ(), static_cast<double>(mcCollision.nChargedInEtaRange())};
        fillAssocsMcTrue<HistFillMode::Correlation, TruthLevel::True>(histFillInfo, assocParticle, isDecaySelfTrig);
      };
      corrProcessCorrelation(mcCollision, trueTriggersThisEvent, mcParticlesThisEvent, funcCorrelation);

      // select mixing events
      if (config.doTrigEvMixing.value && !mcCollision.trigEv()) {
        continue;
      }

      auto const funcMixing = [this](auto const& mcCollision, auto const& mixingTrigger, auto const& assocParticle, auto const perTriggerWeight) {
        HistFillInfoCorr histFillInfo = {getDeltaPhi(mixingTrigger.phi(), assocParticle.phi()), mixingTrigger.eta() - assocParticle.eta(),
                                         mixingTrigger.pt(), assocParticle.pt(), mcCollision.posZ(), static_cast<double>(mcCollision.nChargedInEtaRange())};
        histFillInfo.weight = perTriggerWeight;
        fillAssocsMcTrue<HistFillMode::Mixing, TruthLevel::True>(histFillInfo, assocParticle);
      };
      corrProcessMixing<TriggerMemoryMode::True>(mcCollision, mcParticlesThisEvent, funcMixing, config.nMixingAt0McTrue.value);
    }
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processMcTrue, "process mc-true (all collisions) for multiple associated particles", false);

  void processMcTrueAssocEv(CorrCollisionsMcD const& collisions, CorrMcCollisions const&, aod::TrueTriggers const& trueTriggers, aod::JetParticles const& mcParticles)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      if (!isMcMatchedCollision(collision)) {
        continue;
      }

      auto const& mcCollision = collision.jetMcCollision_as<CorrMcCollisions>();

      // group collision
      auto const trueTriggersThisEvent = trueTriggers.sliceBy(perColTrueTriggers, collision.jetMcCollisionId());
      auto const mcParticlesThisEvent = mcParticles.sliceBy(perColMcParticles, collision.jetMcCollisionId());

      // self correlation of events (plain)
      std::map<AssocType, bool> const decaySelfTrigPlainMap = getDecaySelfTrigPlainMap<TruthLevel::TrueAssocEv>(trueTriggersThisEvent);

      auto const funcPlain = [this, &decaySelfTrigPlainMap](auto const& mcCollision, auto const& assocParticle) {
        // self trigger info for decaying particles
        auto const isDecaySelfTrig = [&decaySelfTrigPlainMap](AssocType const assocType) {
          return decaySelfTrigPlainMap.at(assocType);
        };
        HistFillInfoPlainPhiEtaZPvMult histFillInfo = {mcCollision.ptMax(), assocParticle.pt(), assocParticle.phi(), assocParticle.eta(), mcCollision.posZ(), static_cast<double>(mcCollision.nChargedInEtaRange())};
        fillAssocsMcTrue<HistFillMode::Plain, TruthLevel::TrueAssocEv>(histFillInfo, assocParticle, isDecaySelfTrig);
      };
      corrProcessPlain(mcCollision, mcParticlesThisEvent, funcPlain);

      auto const funcCorrelation = [this](auto const& mcCollision, auto const& trueTrigger, auto const& assocParticle) {
        auto const& triggerParticle = trueTrigger.jetMcParticle();
        // self trigger info for decaying particles
        auto const isDecaySelfTrig = [this, &assocParticle, &triggerParticle]() {
          return producesTrigger(assocParticle, triggerParticle);
        };
        // exclude self correlation
        if (trueTrigger.jetMcParticleId() == assocParticle.globalIndex()) {
          return;
        }

        HistFillInfoCorr histFillInfo = {getDeltaPhi(trueTrigger.phi(), assocParticle.phi()), trueTrigger.eta() - assocParticle.eta(),
                                         trueTrigger.pt(), assocParticle.pt(), mcCollision.posZ(), static_cast<double>(mcCollision.nChargedInEtaRange())};
        fillAssocsMcTrue<HistFillMode::Correlation, TruthLevel::TrueAssocEv>(histFillInfo, assocParticle, isDecaySelfTrig);
      };
      corrProcessCorrelation(mcCollision, trueTriggersThisEvent, mcParticlesThisEvent, funcCorrelation);

      // select mixing events
      if (config.doTrigEvMixing.value && !collision.trigEv()) {
        continue;
      }

      auto const funcMixing = [this](auto const& mcCollision, auto const& mixingTrigger, auto const& assocParticle, auto const perTriggerWeight) {
        HistFillInfoCorr histFillInfo = {getDeltaPhi(mixingTrigger.phi(), assocParticle.phi()), mixingTrigger.eta() - assocParticle.eta(),
                                         mixingTrigger.pt(), assocParticle.pt(), mcCollision.posZ(), static_cast<double>(mcCollision.nChargedInEtaRange())};
        histFillInfo.weight = perTriggerWeight;
        fillAssocsMcTrue<HistFillMode::Mixing, TruthLevel::TrueAssocEv>(histFillInfo, assocParticle);
      };
      corrProcessMixing<TriggerMemoryMode::TrueAssocEv>(mcCollision, mcParticlesThisEvent, funcMixing, config.nMixingAt0McTrue.value);
    }
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processMcTrueAssocEv, "process mc-true (reco collisions) correlation for multiple associated particles", false);

  void processMcTrueAssocEvRecoPtTrig(CorrCollisionsMcD const& collisions, soa::Join<aod::Triggers, aod::TriggerTrackIds> const& triggers, aod::JetTracksMCD const&,
                                      CorrMcCollisions const&, aod::JetParticles const& mcParticles)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      if (!isMcMatchedCollision(collision)) {
        continue;
      }

      auto const& mcCollision = collision.jetMcCollision_as<CorrMcCollisions>();

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const mcParticlesThisEvent = mcParticles.sliceBy(perColMcParticles, collision.jetMcCollisionId());

      // self correlation of events (plain)
      std::map<AssocType, bool> const decaySelfTrigPlainMap = getDecaySelfTrigPlainMap<TruthLevel::TrueAssocEvRecoPtTrig>(triggersThisEvent);

      auto const funcPlain = [this, &decaySelfTrigPlainMap](auto const& collision, auto const& assocParticle) {
        auto const& mcCollision = collision.template jetMcCollision_as<CorrMcCollisions>();
        // self trigger info for decaying particles
        auto const isDecaySelfTrig = [&decaySelfTrigPlainMap](AssocType const assocType) {
          return decaySelfTrigPlainMap.at(assocType);
        };
        HistFillInfoPlainPhiEtaZPvMult histFillInfo = {collision.ptMax(), assocParticle.pt(), assocParticle.phi(), assocParticle.eta(), mcCollision.posZ(), static_cast<double>(mcCollision.nChargedInEtaRange())};
        fillAssocsMcTrue<HistFillMode::Plain, TruthLevel::TrueAssocEvRecoPtTrig>(histFillInfo, assocParticle, isDecaySelfTrig);
      };
      corrProcessPlain(collision, mcParticlesThisEvent, funcPlain);

      auto const funcCorrelation = [this](auto const& mcCollision, auto const& trigger, auto const& assocParticle) {
        if (!isMcMatchedTrigger(trigger)) {
         return;
        }
        auto const& triggerParticle = trigger.template jetTrack_as<aod::JetTracksMCD>().mcParticle();
        // self trigger info for decaying particles
        auto const isDecaySelfTrig = [this, &assocParticle, &triggerParticle]() {
          return producesTrigger(assocParticle, triggerParticle);
        };
        // exclude self correlation
        if (triggerParticle.globalIndex() == assocParticle.globalIndex()) {
          return;
        }
        // match reco self-correlation exclusions
        // trigger from conversion leg of the photon (pi0,eta legs handled implicily in isDecaySelfTrig)
        if (assocParticle.pdgCode() == PDG_t::kGamma &&
            triggerParticle.has_mothers() && triggerParticle.mothersIds().front() == assocParticle.globalIndex() &&
            triggerParticle.getProcess() == TMCProcess::kPPair) {
          return;
        }

        HistFillInfoCorr histFillInfo = {getDeltaPhi(triggerParticle.phi(), assocParticle.phi()), triggerParticle.eta() - assocParticle.eta(),
                                         trigger.pt(), assocParticle.pt(), mcCollision.posZ(), static_cast<double>(mcCollision.nChargedInEtaRange())};
        fillAssocsMcTrue<HistFillMode::Correlation, TruthLevel::TrueAssocEvRecoPtTrig>(histFillInfo, assocParticle, isDecaySelfTrig);
      };
      corrProcessCorrelation(mcCollision, triggersThisEvent, mcParticlesThisEvent, funcCorrelation);

      // select mixing events
      if (config.doTrigEvMixing.value && !collision.trigEv()) {
        continue;
      }

      auto const funcMixing = [this](auto const& mcCollision, auto const& mixingTrigger, auto const& assocParticle, auto const perTriggerWeight) {
        HistFillInfoCorr histFillInfo = {getDeltaPhi(mixingTrigger.phi(), assocParticle.phi()), mixingTrigger.eta() - assocParticle.eta(),
                                         mixingTrigger.pt(), assocParticle.pt(), mcCollision.posZ(), static_cast<double>(mcCollision.nChargedInEtaRange())};
        histFillInfo.weight = perTriggerWeight;
        fillAssocsMcTrue<HistFillMode::Mixing, TruthLevel::TrueAssocEvRecoPtTrig>(histFillInfo, assocParticle, [](){ return false; });
      };
      corrProcessMixing<TriggerMemoryMode::TrueAssocEvRecoPtTrig>(mcCollision, mcParticlesThisEvent, funcMixing, config.nMixingAt0McTrue.value);
    }
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processMcTrueAssocEvRecoPtTrig, "process mc-true (reco collisions with reco ptTrig) correlation for multiple associated particles", false);

  // standard mc reco++ histogram compile-time paths
  static constexpr auto HistPathsMcRecoStd = []() {
    using enum HistFillMode;
    using enum TruthLevel;
    using enum AssocType;

    return std::to_array<HistPathEntry>({
      {PlainPhiEta,  RecoAssocEv, Hadron,    "mc/plain/h4_ptTrigPtAssocPhiEta_recoAssocEv_hadron"},
      {PlainZPvMult, RecoAssocEv, Hadron,    "mc/plain/h4_ptTrigPtAssocZPvMult_recoAssocEv_hadron"},
      {Correlation,  RecoAssocEv, Hadron,    "mc/corr/h6_corr_recoAssocEv_hadron"},
      {Mixing,       RecoAssocEv, Hadron,    "mc/corr/h6_mix_recoAssocEv_hadron"},
      {PlainPhiEta,  RecoAssocEv, Pipm,      "mc/plain/h4_ptTrigPtAssocPhiEta_recoAssocEv_pipm"},
      {PlainZPvMult, RecoAssocEv, Pipm,      "mc/plain/h4_ptTrigPtAssocZPvMult_recoAssocEv_pipm"},
      {Correlation,  RecoAssocEv, Pipm,      "mc/corr/h6_corr_recoAssocEv_pipm"},
      {Mixing,       RecoAssocEv, Pipm,      "mc/corr/h6_mix_recoAssocEv_pipm"},
      {PlainPhiEta,  RecoAssocEv, PhotonPcm, "mc/plain/h4_ptTrigPtAssocPhiEta_recoAssocEv_photonPcm"},
      {PlainZPvMult, RecoAssocEv, PhotonPcm, "mc/plain/h4_ptTrigPtAssocZPvMult_recoAssocEv_photonPcm"},
      {Correlation,  RecoAssocEv, PhotonPcm, "mc/corr/h6_corr_recoAssocEv_photonPcm"},
      {Mixing,       RecoAssocEv, PhotonPcm, "mc/corr/h6_mix_recoAssocEv_photonPcm"},

      {PlainPhiEta,  RecoPure, Hadron,    "mc/plain/h4_ptTrigPtAssocPhiEta_recoPure_hadron"},
      {PlainZPvMult, RecoPure, Hadron,    "mc/plain/h4_ptTrigPtAssocZPvMult_recoPure_hadron"},
      {Correlation,  RecoPure, Hadron,    "mc/corr/h6_corr_recoPure_hadron"},
      {Mixing,       RecoPure, Hadron,    "mc/corr/h6_mix_recoPure_hadron"},
      {PlainPhiEta,  RecoPure, Pipm,      "mc/plain/h4_ptTrigPtAssocPhiEta_recoPure_pipm"},
      {PlainZPvMult, RecoPure, Pipm,      "mc/plain/h4_ptTrigPtAssocZPvMult_recoPure_pipm"},
      {Correlation,  RecoPure, Pipm,      "mc/corr/h6_corr_recoPure_pipm"},
      {Mixing,       RecoPure, Pipm,      "mc/corr/h6_mix_recoPure_pipm"},
      {PlainPhiEta,  RecoPure, PhotonPcm, "mc/plain/h4_ptTrigPtAssocPhiEta_recoPure_photonPcm"},
      {PlainZPvMult, RecoPure, PhotonPcm, "mc/plain/h4_ptTrigPtAssocZPvMult_recoPure_photonPcm"},
      {Correlation,  RecoPure, PhotonPcm, "mc/corr/h6_corr_recoPure_photonPcm"},
      {Mixing,       RecoPure, PhotonPcm, "mc/corr/h6_mix_recoPure_photonPcm"},

      {PlainPhiEta,  RecoPureTruePtAssoc, Hadron,    "mc/plain/h4_ptTrigPtAssocPhiEta_recoPureTruePtAssoc_hadron"},
      {PlainZPvMult, RecoPureTruePtAssoc, Hadron,    "mc/plain/h4_ptTrigPtAssocZPvMult_recoPureTruePtAssoc_hadron"},
      {Correlation,  RecoPureTruePtAssoc, Hadron,    "mc/corr/h6_corr_recoPureTruePtAssoc_hadron"},
      {Mixing,       RecoPureTruePtAssoc, Hadron,    "mc/corr/h6_mix_recoPureTruePtAssoc_hadron"},
      {PlainPhiEta,  RecoPureTruePtAssoc, Pipm,      "mc/plain/h4_ptTrigPtAssocPhiEta_recoPureTruePtAssoc_pipm"},
      {PlainZPvMult, RecoPureTruePtAssoc, Pipm,      "mc/plain/h4_ptTrigPtAssocZPvMult_recoPureTruePtAssoc_pipm"},
      {Correlation,  RecoPureTruePtAssoc, Pipm,      "mc/corr/h6_corr_recoPureTruePtAssoc_pipm"},
      {Mixing,       RecoPureTruePtAssoc, Pipm,      "mc/corr/h6_mix_recoPureTruePtAssoc_pipm"},
      {PlainPhiEta,  RecoPureTruePtAssoc, PhotonPcm, "mc/plain/h4_ptTrigPtAssocPhiEta_recoPureTruePtAssoc_photonPcm"},
      {PlainZPvMult, RecoPureTruePtAssoc, PhotonPcm, "mc/plain/h4_ptTrigPtAssocZPvMult_recoPureTruePtAssoc_photonPcm"},
      {Correlation,  RecoPureTruePtAssoc, PhotonPcm, "mc/corr/h6_corr_recoPureTruePtAssoc_photonPcm"},
      {Mixing,       RecoPureTruePtAssoc, PhotonPcm, "mc/corr/h6_mix_recoPureTruePtAssoc_photonPcm"}
    });
  }();

  template <HistFillMode HistFillModeV, AssocType AssocTypeV>
  void fillTruthLevelsMcRecoStd(auto histFillInfo, auto const& collision, auto const& associated) {
    // photon impurity
    if constexpr (HistFillModeV == HistFillMode::Plain && AssocTypeV == AssocType::PhotonPcm) {
      PhotonImpurity const impurity = static_cast<PhotonImpurity>(associated.impurityPhotonPcm());
      histos.fill(HIST("mc/plain/h3_ptTrigPtAssocCategory_recoImpurities_photonPcm"),
                  collision.ptMax(), associated.pt(), binValuePhotonImpurity(impurity),
                  histFillInfo.weight);
    }
    // collision association
    if (!checkTruthLevelReco<TruthLevel::RecoAssocEv, AssocTypeV>(associated)) {
      return;
    }
    fillHistMain<HistPathsMcRecoStd, HistFillModeV, TruthLevel::RecoAssocEv, AssocTypeV>(histFillInfo);
    // purity
    if (!checkTruthLevelReco<TruthLevel::RecoPure, AssocTypeV>(associated)) {
      return;
    }
    fillHistMain<HistPathsMcRecoStd, HistFillModeV, TruthLevel::RecoPure, AssocTypeV>(histFillInfo);
    // true ptAssoc
    histFillInfo.ptAssoc = associated.ptTruePure();
    fillHistMain<HistPathsMcRecoStd, HistFillModeV, TruthLevel::RecoPureTruePtAssoc, AssocTypeV>(histFillInfo);
  }

  template <AssocType AssocTypeV>
  void fillProcessMcRecoStd(auto const& collisions, auto const& triggers, auto const& assocs)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      // group collision ++
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const [assocsThisEvent, nMixingAt0] = [&]() {
        if constexpr (AssocTypeV == AssocType::Hadron) {
          return std::tuple{assocs.sliceBy(perColHadrons, collision.globalIndex()), config.nMixingAt0Hadron.value};
        } else if constexpr (AssocTypeV == AssocType::Pipm) {
          return std::tuple{assocs.sliceBy(perColPipms, collision.globalIndex()), config.nMixingAt0Pipm.value};
        } else if constexpr (AssocTypeV == AssocType::PhotonPcm) {
          return std::tuple{assocs.sliceBy(perColPhotonPcms, collision.globalIndex()), config.nMixingAt0PhotonPcm.value};
        } else {
          static_assert(AlwaysFalseValue<AssocTypeV>, "assocType not supported");
        }
      }();

      auto const funcPlain = [this](auto const& collision, auto const& associated) {
        // reco assoc selection
        if (!isSelectedAssoc<AssocTypeV>(associated)) {
          return;
        }

        HistFillInfoPlainPhiEtaZPvMult histFillInfo = {collision.ptMax(), associated.pt(), associated.phi(), associated.eta(), collision.posZ(), static_cast<double>(collision.nGlobalTracks())};
        histFillInfo.weight = getMcAssocPtScaling<AssocTypeV>(associated) * getMcK0ShortPtScaling<AssocTypeV>(associated);
        fillTruthLevelsMcRecoStd<HistFillMode::Plain, AssocTypeV>(histFillInfo, collision, associated);
      };
      corrProcessPlain(collision, assocsThisEvent, funcPlain);

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        // reco assoc selection
        if (!isSelectedAssoc<AssocTypeV>(associated)) {
          return;
        }
        // exclude self correlation
        if (isRecoSelfTrigger(trigger, associated)) {
          return;
        }

        HistFillInfoCorr histFillInfo = {getDeltaPhi(trigger.phi(), associated.phi()), trigger.eta() - associated.eta(),
                                         trigger.pt(), associated.pt(), collision.posZ(), static_cast<double>(collision.nGlobalTracks())};
        histFillInfo.weight = getMcAssocPtScaling<AssocTypeV>(associated) * getMcK0ShortPtScaling<AssocTypeV>(associated) * getTriggerPtCorrection(trigger.pt());
        fillTruthLevelsMcRecoStd<HistFillMode::Correlation, AssocTypeV>(histFillInfo, collision, associated);
      };
      corrProcessCorrelation(collision, triggersThisEvent, assocsThisEvent, funcCorrelation);

      // select mixing events
      if (config.doTrigEvMixing.value && !collision.trigEv()) {
        continue;
      }

      auto const conditionMixAssocOnly = [this](auto const& associated) {
        // reco assoc selection
        return isSelectedAssoc<AssocTypeV>(associated);
      };
      auto const funcMixing = [this](auto const& collision, auto const& mixingTrigger, auto const& associated, auto const perTriggerWeight) {
        HistFillInfoCorr histFillInfo = {getDeltaPhi(mixingTrigger.phi(), associated.phi()), mixingTrigger.eta() - associated.eta(),
                                         mixingTrigger.pt(), associated.pt(), collision.posZ(), static_cast<double>(collision.nGlobalTracks())};
        histFillInfo.weight = perTriggerWeight * getMcAssocPtScaling<AssocTypeV>(associated) * getMcK0ShortPtScaling<AssocTypeV>(associated) * getTriggerPtCorrection(mixingTrigger.pt());
        fillTruthLevelsMcRecoStd<HistFillMode::Mixing, AssocTypeV>(histFillInfo, collision, associated);
      };
      corrProcessMixing<TriggerMemoryMode::Reco>(collision, assocsThisEvent, funcMixing, nMixingAt0, conditionMixAssocOnly);
    }
  }

  void processMcRecoHadron(CorrCollisions const& collisions, aod::Triggers const& triggers, aod::Hadrons const& hadrons)
  {
    fillProcessMcRecoStd<AssocType::Hadron>(collisions, triggers, hadrons);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processMcRecoHadron, "process standard for associated hardon with additional mc information", false);

  void processMcRecoPipm(CorrCollisions const& collisions, aod::Triggers const& triggers, aod::Pipms const& pipms)
  {
    fillProcessMcRecoStd<AssocType::Pipm>(collisions, triggers, pipms);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processMcRecoPipm, "process standard for associated pipm with additional mc information", false);

  void processMcRecoPhotonPcm(CorrCollisions const& collisions, aod::Triggers const& triggers, PhotonPcmsWithV0 const& photonPcms,
                              aod::V0PhotonsKF const&, aod::V0Legs const&)
  {
    fillProcessMcRecoStd<AssocType::PhotonPcm>(collisions, triggers, photonPcms);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processMcRecoPhotonPcm, "process standard for associated photonPcm with additional mc information", false);

  // photonPcmPair mc reco++ histogram compile-time paths
  static constexpr auto HistPathsMcRecoPhotonPcmPair = []() {
    using enum HistFillMode;
    using enum TruthLevel;
    using enum AssocType;

    return std::to_array<HistPathEntry>({
      {PlainPt,     PseudoReco, Pi0Pcm, "mc/plain/h2_ptTrigPtAssoc_pseudoReco_pi0Pcm"},
      {Correlation, PseudoReco, Pi0Pcm, "mc/corr/h6_corr_pseudoReco_pi0Pcm"},
      {Mixing,      PseudoReco, Pi0Pcm, "mc/corr/h6_mix_pseudoReco_pi0Pcm"},
      {PlainPt,     PseudoReco, EtaPcm, "mc/plain/h2_ptTrigPtAssoc_pseudoReco_etaPcm"},
      {Correlation, PseudoReco, EtaPcm, "mc/corr/h6_corr_pseudoReco_etaPcm"},
      {Mixing,      PseudoReco, EtaPcm, "mc/corr/h6_mix_pseudoReco_etaPcm"},

      {PlainPt,     PseudoRecoAssocEv, Pi0Pcm, "mc/plain/h2_ptTrigPtAssoc_pseudoRecoAssocEv_pi0Pcm"},
      {Correlation, PseudoRecoAssocEv, Pi0Pcm, "mc/corr/h6_corr_pseudoRecoAssocEv_pi0Pcm"},
      {Mixing,      PseudoRecoAssocEv, Pi0Pcm, "mc/corr/h6_mix_pseudoRecoAssocEv_pi0Pcm"},
      {PlainPt,     PseudoRecoAssocEv, EtaPcm, "mc/plain/h2_ptTrigPtAssoc_pseudoRecoAssocEv_etaPcm"},
      {Correlation, PseudoRecoAssocEv, EtaPcm, "mc/corr/h6_corr_pseudoRecoAssocEv_etaPcm"},
      {Mixing,      PseudoRecoAssocEv, EtaPcm, "mc/corr/h6_mix_pseudoRecoAssocEv_etaPcm"},

      {PlainPt,     PseudoRecoPure, Pi0Pcm, "mc/plain/h2_ptTrigPtAssoc_pseudoRecoPure_pi0Pcm"},
      {Correlation, PseudoRecoPure, Pi0Pcm, "mc/corr/h6_corr_pseudoRecoPure_pi0Pcm"},
      {Mixing,      PseudoRecoPure, Pi0Pcm, "mc/corr/h6_mix_pseudoRecoPure_pi0Pcm"},
      {PlainPt,     PseudoRecoPure, EtaPcm, "mc/plain/h2_ptTrigPtAssoc_pseudoRecoPure_etaPcm"},
      {Correlation, PseudoRecoPure, EtaPcm, "mc/corr/h6_corr_pseudoRecoPure_etaPcm"},
      {Mixing,      PseudoRecoPure, EtaPcm, "mc/corr/h6_mix_pseudoRecoPure_etaPcm"},

      {PlainPt,     PseudoRecoPureTruePtAssoc, Pi0Pcm, "mc/plain/h2_ptTrigPtAssoc_pseudoRecoPureTruePtAssoc_pi0Pcm"},
      {Correlation, PseudoRecoPureTruePtAssoc, Pi0Pcm, "mc/corr/h6_corr_pseudoRecoPureTruePtAssoc_pi0Pcm"},
      {Mixing,      PseudoRecoPureTruePtAssoc, Pi0Pcm, "mc/corr/h6_mix_pseudoRecoPureTruePtAssoc_pi0Pcm"},
      {PlainPt,     PseudoRecoPureTruePtAssoc, EtaPcm, "mc/plain/h2_ptTrigPtAssoc_pseudoRecoPureTruePtAssoc_etaPcm"},
      {Correlation, PseudoRecoPureTruePtAssoc, EtaPcm, "mc/corr/h6_corr_pseudoRecoPureTruePtAssoc_etaPcm"},
      {Mixing,      PseudoRecoPureTruePtAssoc, EtaPcm, "mc/corr/h6_mix_pseudoRecoPureTruePtAssoc_etaPcm"}
    });
  }();

  template <HistFillMode HistFillModeV, AssocType AssocTypeV>
  void fillTruthLevelsMcRecoPhotonPcmPair(auto histFillInfo, auto const& collision, auto const& associated) {
    // pi0 impurity
    if constexpr (HistFillModeV == HistFillMode::Plain && AssocTypeV == AssocType::Pi0Pcm) {
      PhotonImpurity const impurity = static_cast<PhotonImpurity>(associated.impurityPi0Pcm());
      histos.fill(HIST("mc/plain/h3_ptTrigPtAssocCategory_pseudoRecoImpurities_pi0Pcm"),
                  collision.ptMax(), associated.pt(), binValuePhotonImpurity(impurity),
                  histFillInfo.weight);
    }
    // collision association
    if (!checkTruthLevelReco<TruthLevel::PseudoRecoAssocEv, AssocTypeV>(associated)) {
      return;
    }
    fillHistMain<HistPathsMcRecoPhotonPcmPair, HistFillModeV, TruthLevel::PseudoRecoAssocEv, AssocTypeV>(histFillInfo);
    // purity
    if (!checkTruthLevelReco<TruthLevel::PseudoRecoPure, AssocTypeV>(associated)) {
      return;
    }
    fillHistMain<HistPathsMcRecoPhotonPcmPair, HistFillModeV, TruthLevel::PseudoRecoPure, AssocTypeV>(histFillInfo);
    // true ptAssoc
    histFillInfo.ptAssoc = associated.ptTruePure();
    fillHistMain<HistPathsMcRecoPhotonPcmPair, HistFillModeV, TruthLevel::PseudoRecoPureTruePtAssoc, AssocTypeV>(histFillInfo);
  }

  template <HistFillMode HistFillModeV>
  void fillAssocsMcRecoPhotonPcmPair(auto histFillInfo, auto const& collision, auto const& associated) {
    // pi0Pcm
    if (checkTruthLevelReco<TruthLevel::PseudoReco, AssocType::Pi0Pcm>(associated)) {
      // pseudo yield

      // mgg distribution
      if constexpr (HistFillModeV == HistFillMode::Plain) {
        histos.fill(HIST("mc/plain/h5_ptTrigPtAssocMggZPvMult_pseudoReco_pi0Pcm"),
                    histFillInfo.ptMax, histFillInfo.ptAssoc, associated.mgg(), collision.posZ(), collision.nGlobalTracks(),
                    histFillInfo.weight);
      } else if constexpr (HistFillModeV == HistFillMode::Correlation) {
        histos.fill(HIST("mc/corr/h5_ptTrigPtAssocMggZPvMult_pseudoReco_pi0Pcm"),
                    histFillInfo.ptTrig, histFillInfo.ptAssoc, associated.mgg(), collision.posZ(), collision.nGlobalTracks(),
                    histFillInfo.weight);
      }

      // mgg distribution

      fillHistMain<HistPathsMcRecoPhotonPcmPair, HistFillModeV, TruthLevel::PseudoReco, AssocType::Pi0Pcm>(histFillInfo);

      fillTruthLevelsMcRecoPhotonPcmPair<HistFillModeV, AssocType::Pi0Pcm>(histFillInfo, collision, associated);
      return;
    }
    // etaPcm
    if (checkTruthLevelReco<TruthLevel::PseudoReco, AssocType::EtaPcm>(associated)) {
      // pseudo yield

      // mgg distribution
      if constexpr (HistFillModeV == HistFillMode::Plain) {
        histos.fill(HIST("mc/plain/h5_ptTrigPtAssocMggZPvMult_pseudoReco_etaPcm"),
                    histFillInfo.ptMax, histFillInfo.ptAssoc, associated.mgg(), collision.posZ(), collision.nGlobalTracks(),
                    histFillInfo.weight);
      } else if constexpr (HistFillModeV == HistFillMode::Correlation) {
        histos.fill(HIST("mc/corr/h5_ptTrigPtAssocMggZPvMult_pseudoReco_etaPcm"),
                    histFillInfo.ptTrig, histFillInfo.ptAssoc, associated.mgg(), collision.posZ(), collision.nGlobalTracks(),
                    histFillInfo.weight);
      }

      fillHistMain<HistPathsMcRecoPhotonPcmPair, HistFillModeV, TruthLevel::PseudoReco, AssocType::EtaPcm>(histFillInfo);

      fillTruthLevelsMcRecoPhotonPcmPair<HistFillModeV, AssocType::EtaPcm>(histFillInfo, collision, associated);
      return;
    }
  }

  void fillProcesMcRecoPhotonPcmPair(auto const& collisions, auto const& triggers, auto const& photonPcmPairs) {
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const photonPcmPairsThisEvent = photonPcmPairs.sliceBy(perColPhotonPcmPairs, collision.globalIndex());

      auto const funcPlain = [this](auto const& collision, auto const& associated) {
        // reco assoc selection
        if (!isSelectedAssoc<AssocType::PhotonPcmPair>(associated)) {
          return;
        }

        HistFillInfoPlainPt histFillInfo = {collision.ptMax(), associated.pt()};
        histFillInfo.weight = getMcAssocPtScaling<AssocType::PhotonPcmPair>(associated) * getMcK0ShortPtScaling<AssocType::PhotonPcmPair>(associated);
        fillAssocsMcRecoPhotonPcmPair<HistFillMode::Plain>(histFillInfo, collision, associated);
      };
      corrProcessPlain(collision, photonPcmPairsThisEvent, funcPlain);

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        // reco assoc selection
        if (!isSelectedAssoc<AssocType::PhotonPcmPair>(associated)) {
          return;
        }
        // exclude self correlation
        if (isRecoSelfTrigger(trigger, associated)) {
          return;
        }

        HistFillInfoCorr histFillInfo = {getDeltaPhi(trigger.phi(), associated.phi()), trigger.eta() - associated.eta(),
                                         trigger.pt(), associated.pt(), collision.posZ(), static_cast<double>(collision.nGlobalTracks())};
        histFillInfo.weight = getMcAssocPtScaling<AssocType::PhotonPcmPair>(associated) * getMcK0ShortPtScaling<AssocType::PhotonPcmPair>(associated) * getTriggerPtCorrection(trigger.pt());
        fillAssocsMcRecoPhotonPcmPair<HistFillMode::Correlation>(histFillInfo, collision, associated);
      };
      corrProcessCorrelation(collision, triggersThisEvent, photonPcmPairsThisEvent, funcCorrelation);

      // select mixing events
      if (config.doTrigEvMixing.value && !collision.trigEv()) {
        continue;
      }

      auto const conditionMixAssocOnly = [this](auto const& associated) {
        // reco assoc selection
        return isSelectedAssoc<AssocType::PhotonPcmPair>(associated);
      };
      auto const funcMixing = [this](auto const& collision, auto const& mixingTrigger, auto const& associated, auto const perTriggerWeight) {
        HistFillInfoCorr histFillInfo = {getDeltaPhi(mixingTrigger.phi(), associated.phi()), mixingTrigger.eta() - associated.eta(),
                                         mixingTrigger.pt(), associated.pt(), collision.posZ(), static_cast<double>(collision.nGlobalTracks())};
        histFillInfo.weight = perTriggerWeight * getMcAssocPtScaling<AssocType::PhotonPcmPair>(associated) * getMcK0ShortPtScaling<AssocType::PhotonPcmPair>(associated) * getTriggerPtCorrection(mixingTrigger.pt());
        fillAssocsMcRecoPhotonPcmPair<HistFillMode::Mixing>(histFillInfo, collision, associated);
      };
      corrProcessMixing<TriggerMemoryMode::Reco>(collision, photonPcmPairsThisEvent, funcMixing, config.nMixingAt0PhotonPcmPair.value, conditionMixAssocOnly);
    }
  }

  void processMcRecoPhotonPcmPair(CorrCollisions const& collisions, aod::Triggers const& triggers, PhotonPcmPairsWithV0 const& photonPcmPairs,
                                  aod::V0PhotonsKF const&, aod::V0Legs const&)
  {
    fillProcesMcRecoPhotonPcmPair(collisions, triggers, photonPcmPairs);
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processMcRecoPhotonPcmPair, "process correlation for associated photonPcmPair with additional mc information", false);

  // other ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  void processPcmDcaZTrue(CorrCollisionsMcD::iterator const& collision, soa::Join<PhotonPcmsWithV0, aod::PhotonPcmTrackIds> const& photonPcms,
                          aod::V0PhotonsKF const&, aod::V0Legs const&, aod::JetTracksMCD const&, CorrMcCollisions const&, aod::JetParticles const&)
  {
    // event selection
    if (!totalEvSel(collision)) {
      return;
    }

    for (auto const& photonPcm : photonPcms) {
      // reco assoc selection
      if (!isSelectedAssoc<AssocType::PhotonPcm>(photonPcm)) {
        continue;
      }
      if (!checkTruthLevelReco<TruthLevel::RecoMcMatched, AssocType::PhotonPcm>(photonPcm)) {
        continue;
      }
      // check mc
      auto const& posTrack = photonPcm.posJetTrack_as<aod::JetTracksMCD>();
      auto const& negTrack = photonPcm.negJetTrack_as<aod::JetTracksMCD>();
      // collision association
      if (posTrack.mcParticle().mcCollisionId() != collision.jetMcCollisionId() || negTrack.mcParticle().mcCollisionId() != collision.jetMcCollisionId()) {
        continue;
      }

      double const weight = getMcAssocPtScaling<AssocType::PhotonPcm>(photonPcm) * getMcK0ShortPtScaling<AssocType::PhotonPcm>(photonPcm);

      histos.fill(HIST("mc/plain/h5_ptTrigPtAssocDcaZZPvMult_photonPcm"),
                  collision.ptMax(), photonPcm.pt(), photonPcm.v0PhotonKF().dcaZtopv(), collision.posZ(), collision.nGlobalTracks(),
                  weight);
    }
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processPcmDcaZTrue, "process to test true dcaZ distribution of photons in trigger collisions", false);

  // last /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  void processRecoLast(CorrCollisions const& collisions, aod::Triggers const& triggers)
  {
    // do at end of each data frame (after other reco correlation process functions)
    // (PROCESS_SWITCH of this process has to be declared last)

    for (auto const& collision : collisions) {
      // all events
      histos.fill(HIST("reco/info/h1_nEvents"), 0.5);

      // mc split
      if (!checkSplitMcEventSelection(collision)) {
        continue;
      }
      histos.fill(HIST("reco/info/h1_nEvents"), 1.5);

      // standard event selection
      if (!collision.selEv()) {
        continue;
      }
      histos.fill(HIST("reco/info/h1_nEvents"), 2.5);

      histos.fill(HIST("reco/info/h3_nEv_ptTrigZPvMult_reco"), collision.ptMax(), collision.posZ(), collision.nGlobalTracks());

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());

      // trigger loop
      for (auto const& trigger : triggersThisEvent) {
        // trigger info
        histos.fill(HIST("reco/corr/h3_ptPhiEta_trig"), trigger.pt(), trigger.phi(), trigger.eta(),
                    getTriggerPtCorrection(trigger.pt()));

        // save trigger for mixing
        triggerMemoryReco->saveTrigger(trigger.pt(), trigger.phi(), trigger.eta(), collision.posZ(), collision.nGlobalTracks());
      }
    }
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processRecoLast, "process general info for reconstructed collisions and later correlations last", true);

  void processMcLast(CorrCollisionsMcD const& collisions, soa::Join<aod::Triggers, aod::TriggerTrackIds> const& triggers, aod::JetTracksMCD const&,
                     CorrMcCollisions const& mcCollisions, aod::TrueTriggers const& trueTriggers, aod::JetParticles const&)
  {
    // do at end of each data frame (after other mc correlation process functions)
    // (PROCESS_SWITCH of this process has to be declared last)

    for (auto const& mcCollision : mcCollisions) {
      // event selection
      if (!totalEvSel(mcCollision)) {
        continue;
      }

      // group collision
      auto const trueTriggersThisEvent = trueTriggers.sliceBy(perColTrueTriggers, mcCollision.globalIndex());

      histos.fill(HIST("mc/info/h3_nEv_ptTrigZPvMult_true"), mcCollision.ptMax(), mcCollision.posZ(), mcCollision.nChargedInEtaRange());

      // trigger loop
      for (auto const& trueTrigger : trueTriggersThisEvent) {
        // trigger info
        histos.fill(HIST("mc/corr/h3_ptPhiEta_trig_true"), trueTrigger.pt(), trueTrigger.phi(), trueTrigger.eta());

        // save trigger for mixing
        triggerMemoryTrue->saveTrigger(trueTrigger.pt(), trueTrigger.phi(), trueTrigger.eta(), mcCollision.posZ(), mcCollision.nChargedInEtaRange());
      }
    }
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      if (!isMcMatchedCollision(collision)) {
        histos.fill(HIST("mc/info/h1_pt_noMatch_event"), 0.5);
        continue;
      }

      auto const& mcCollision = collision.jetMcCollision_as<CorrMcCollisions>();

      // group collision
      auto const trueTriggersThisEvent = trueTriggers.sliceBy(perColTrueTriggers, mcCollision.globalIndex());
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());

      histos.fill(HIST("mc/info/h3_nEv_ptTrigZPvMult_trueAssocEv"), mcCollision.ptMax(), mcCollision.posZ(), mcCollision.nChargedInEtaRange());
      histos.fill(HIST("mc/info/h3_nEv_ptTrigZPvMult_trueAssocEvRecoPtTrig"), collision.ptMax(), mcCollision.posZ(), mcCollision.nChargedInEtaRange());

      // trigger loops
      for (auto const& trueTrigger : trueTriggersThisEvent) {
        // trigger info
        histos.fill(HIST("mc/corr/h3_ptPhiEta_trig_trueAssocEv"), trueTrigger.pt(), trueTrigger.phi(), trueTrigger.eta());

        // save trigger for mixing
        triggerMemoryTrueAssocEv->saveTrigger(trueTrigger.pt(), trueTrigger.phi(), trueTrigger.eta(), mcCollision.posZ(), mcCollision.nChargedInEtaRange());
      }
      for (auto const& trigger : triggersThisEvent) {
        if (!isMcMatchedTrigger(trigger)) {
          histos.fill(HIST("mc/info/h1_pt_noMatch_trigger"), trigger.pt());
          continue;
        }

        auto const& triggerParticle = trigger.jetTrack_as<aod::JetTracksMCD>().mcParticle();

        // trigger info
        histos.fill(HIST("mc/corr/h3_ptPhiEta_trig_trueAssocEvRecoPtTrig"), trigger.pt(), triggerParticle.phi(), triggerParticle.eta());

        // save trigger for mixing
        triggerMemoryTrueAssocEvRecoPtTrig->saveTrigger(trigger.pt(), triggerParticle.phi(), triggerParticle.eta(), mcCollision.posZ(), mcCollision.nChargedInEtaRange());
      }
    }
  }
  PROCESS_SWITCH(HadronPhotonCorrelation, processMcLast, "process general MC info for collisions and later correlations last", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& configContext)
{
  return WorkflowSpec{adaptAnalysisTask<HadronPhotonCorrelation>(configContext)};
}
