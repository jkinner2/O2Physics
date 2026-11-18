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

///
/// \brief table definitions for hadron-photon correlation analyses
///
/// \author Julius Kinner
/// \file HadronPhotonCorrelation.h

#ifndef PWGJE_DATAMODEL_HADRONPHOTONCORRELATION_H_
#define PWGJE_DATAMODEL_HADRONPHOTONCORRELATION_H_

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include <Framework/AnalysisDataModel.h>

namespace hpc {

enum class AssocType {
  Hadron,
  Pipm,
  PhotonPcm,
  PhotonDecay,
  PhotonDirect,
  PhotonPcmPair,
  Pi0Pcm,
  Pi0PcmPeak,
  Pi0PcmSide,
  PhotonPi0,
  EtaPcm,
  EtaPcmPeak,
  EtaPcmSide,
  PhotonEta,
  Omega,
  PhotonOmega,
  EtaPrime,
  PhotonEtaPrime,
  PhotonOtherMother
};

enum class TruthLevel : int {
  Unknown = -1,
  Reco = 0,
  RecoMcMatched = 1,
  RecoAssocEv = 2,
  RecoPure = 3,
  RecoPureTruePtAssoc = 4,
  PseudoReco = 5,
  PseudoRecoAssocEv = 6,
  PseudoRecoPure = 7,
  PseudoRecoPureTruePtAssoc = 8,
  TrueGeoAcc = 9,
  TrueMeasDecay = 10,
  TrueRecoEvRecoPtTrig = 11,
  TrueRecoEv = 12,
  True = 13,
  TrueDecaySelfTrig = 14
};
// true level: associated particles taken from mc truth
constexpr bool isTrueLevel(TruthLevel const truthLevel)
{
  switch (truthLevel) {
    case TruthLevel::Unknown:
    case TruthLevel::Reco:
    case TruthLevel::RecoMcMatched:
    case TruthLevel::RecoAssocEv:
    case TruthLevel::RecoPure:
    case TruthLevel::RecoPureTruePtAssoc:
    case TruthLevel::PseudoReco:
    case TruthLevel::PseudoRecoAssocEv:
    case TruthLevel::PseudoRecoPure:
    case TruthLevel::PseudoRecoPureTruePtAssoc:
      return false;
    case TruthLevel::TrueGeoAcc:
    case TruthLevel::TrueMeasDecay:
    case TruthLevel::TrueRecoEvRecoPtTrig:
    case TruthLevel::TrueRecoEv:
    case TruthLevel::True:
    case TruthLevel::TrueDecaySelfTrig:
      return true;
  }
  return false;
}

enum class PhotonImpurity : int {
  Unknown = -1,
  Pure = 0,
  WrongEv = 1,
  MisId = 2,
  K0Short = 3,
  K0Long = 4,
  Lambda0 = 5,
  Other = 6
};
constexpr double binValuePhotonImpurity(PhotonImpurity const impurity)
{
  return static_cast<int>(impurity) + 0.5;
}

template <typename>
bool constexpr AlwaysFalseType = false;
template <auto>
bool constexpr AlwaysFalseValue = false;

} // namespace hpc

namespace o2::aod
{

// common correlation columns
namespace hpc_common
{
DECLARE_SOA_INDEX_COLUMN_FULL(JetMcCollision, jetMcCollision, int, JMcCollisions, "");
DECLARE_SOA_COLUMN(SelEv, selEv, bool);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(PtTruePure, ptTruePure, float);
DECLARE_SOA_COLUMN(PtTrueOrigin, ptTrueOrigin, float); // and the blind forest
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
} // namespace hpc_common

// reco

// correlation collision
namespace hpc_collision
{
// plain jet collision column for derived data
DECLARE_SOA_COLUMN(JetCollisionId, jetCollisionId, int);
DECLARE_SOA_COLUMN(ThirdDecimalTruePosZ, thirdDecimalTruePosZ, int);
DECLARE_SOA_COLUMN(TrigEv, trigEv, bool);
DECLARE_SOA_COLUMN(PtMax, ptMax, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(NGlobalTracks, nGlobalTracks, int);
} // namespace hpc_collision
DECLARE_SOA_TABLE(HpcCollisions, "AOD", "HPCCOL",
                  o2::soa::Index<>,
                  hpc_common::SelEv, hpc_collision::ThirdDecimalTruePosZ, hpc_collision::TrigEv, hpc_collision::PtMax, hpc_collision::PosZ, hpc_collision::NGlobalTracks);
DECLARE_SOA_TABLE(HpcCollisionIds, "AOD", "HPCCOLID",
                  hpc_collision::JetCollisionId);
DECLARE_SOA_TABLE(HpcCollisionMcIds, "AOD", "HPCCOLMCID",
                  hpc_common::JetMcCollisionId);

// add to common columns
namespace hpc_common
{
DECLARE_SOA_INDEX_COLUMN_FULL(CorrCollision, corrCollision, int, HpcCollisions, "");
} // namespace hpc_common

// trigger
namespace hpc_trigger
{
DECLARE_SOA_INDEX_COLUMN_FULL(JetTrack, jetTrack, int, JetTracks, "");
} // namespace hpc_trigger
DECLARE_SOA_TABLE(HpcTriggers, "AOD", "HPCTRIG",
                  o2::soa::Index<>, hpc_common::CorrCollisionId,
                  hpc_common::Pt, hpc_common::Phi, hpc_common::Eta);
DECLARE_SOA_TABLE(HpcTriggerTrackIds, "AOD", "HPCTRIGTRID",
                  hpc_trigger::JetTrackId);

// add to common columns
namespace hpc_common
{
DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(IntersectTrigger, intersectTrigger, int, HpcTriggers, "");
} // namespace hpc_common

// hadron (global tracks)
namespace hpc_hadron
{
DECLARE_SOA_INDEX_COLUMN_FULL(JetTrack, jetTrack, int, JetTracks, "");
DECLARE_SOA_COLUMN(TruthLevelRecoHadron, truthLevelRecoHadron, int);
} // namespace hpc_hadron
DECLARE_SOA_TABLE(HpcHadrons, "AOD", "HPCHAD",
                  o2::soa::Index<>, hpc_common::CorrCollisionId, hpc_common::IntersectTriggerIds,
                  hpc_common::Pt, hpc_common::Phi, hpc_common::Eta,
                  hpc_hadron::TruthLevelRecoHadron, hpc_common::PtTruePure);
DECLARE_SOA_TABLE(HpcHadronTrackIds, "AOD", "HPCHADTRID",
                  hpc_hadron::JetTrackId);

// pipm
namespace hpc_pipm
{
DECLARE_SOA_INDEX_COLUMN_FULL(JetTrack, jetTrack, int, JetTracks, "");
DECLARE_SOA_COLUMN(TruthLevelRecoPipm, truthLevelRecoPipm, int);
} // namespace hpc_pipm
DECLARE_SOA_TABLE(HpcPipms, "AOD", "HPCPIPM",
                  o2::soa::Index<>, hpc_common::CorrCollisionId, hpc_common::IntersectTriggerIds,
                  hpc_common::Pt, hpc_common::Phi, hpc_common::Eta,
                  hpc_pipm::TruthLevelRecoPipm, hpc_common::PtTruePure);
DECLARE_SOA_TABLE(HpcPipmTrackIds, "AOD", "HPCPIPMTRID",
                  hpc_pipm::JetTrackId);

// photonPcm
namespace hpc_photon_pcm
{
DECLARE_SOA_INDEX_COLUMN_FULL(V0PhotonKF, v0PhotonKF, int, V0PhotonsKF, "");
DECLARE_SOA_INDEX_COLUMN_FULL(PosJetTrack, posJetTrack, int, JetTracks, "_pos");
DECLARE_SOA_INDEX_COLUMN_FULL(NegJetTrack, negJetTrack, int, JetTracks, "_neg");
DECLARE_SOA_COLUMN(TruthLevelRecoPhotonPcm, truthLevelRecoPhotonPcm, int);
DECLARE_SOA_COLUMN(ImpurityPhotonPcm, impurityPhotonPcm, int);
} // namespace hpc_photon_pcm
DECLARE_SOA_TABLE(HpcPhotonPcms, "AOD", "HPCPHPCM",
                  o2::soa::Index<>, hpc_common::CorrCollisionId, hpc_common::IntersectTriggerIds,
                  hpc_common::Pt, hpc_common::Phi, hpc_common::Eta,
                  hpc_photon_pcm::TruthLevelRecoPhotonPcm, hpc_common::PtTruePure,
                  hpc_photon_pcm::ImpurityPhotonPcm, hpc_common::PtTrueOrigin);
DECLARE_SOA_TABLE(HpcPhotonPcmTrackIds, "AOD", "HPCPHPCMTRID",
                  hpc_photon_pcm::PosJetTrackId, hpc_photon_pcm::NegJetTrackId);
DECLARE_SOA_TABLE(HpcPhotonPcmV0PhotonIds, "AOD", "HPCPHPCMV0ID",
                  hpc_photon_pcm::V0PhotonKFId);

// photonPcm pairs
namespace hpc_photon_pcm_pair
{
DECLARE_SOA_INDEX_COLUMN_FULL(V0PhotonKF1, v0PhotonKF1, int, V0PhotonsKF, "_1");
DECLARE_SOA_INDEX_COLUMN_FULL(V0PhotonKF2, v0PhotonKF2, int, V0PhotonsKF, "_2");
DECLARE_SOA_INDEX_COLUMN_FULL(PosJetTrack1, posJetTrack1, int, JetTracks, "_pos1");
DECLARE_SOA_INDEX_COLUMN_FULL(NegJetTrack1, negJetTrack1, int, JetTracks, "_neg1");
DECLARE_SOA_INDEX_COLUMN_FULL(PosJetTrack2, posJetTrack2, int, JetTracks, "_pos2");
DECLARE_SOA_INDEX_COLUMN_FULL(NegJetTrack2, negJetTrack2, int, JetTracks, "_neg2");
DECLARE_SOA_COLUMN(Mgg, mgg, float);
DECLARE_SOA_COLUMN(TruthLevelRecoPi0Pcm, truthLevelRecoPi0Pcm, int);
DECLARE_SOA_COLUMN(TruthLevelRecoEtaPcm, truthLevelRecoEtaPcm, int);
DECLARE_SOA_COLUMN(ImpurityPi0Pcm, impurityPi0Pcm, int);
} // namespace hpc_photon_pcm_pair
DECLARE_SOA_TABLE(HpcPhotonPcmPairs, "AOD", "HPCPCMPAIR",
                  o2::soa::Index<>, hpc_common::CorrCollisionId, hpc_common::IntersectTriggerIds,
                  hpc_common::Pt, hpc_common::Phi, hpc_common::Eta, hpc_photon_pcm_pair::Mgg,
                  hpc_photon_pcm_pair::TruthLevelRecoPi0Pcm, hpc_photon_pcm_pair::TruthLevelRecoEtaPcm, hpc_common::PtTruePure,
                  hpc_photon_pcm_pair::ImpurityPi0Pcm, hpc_common::PtTrueOrigin);
DECLARE_SOA_TABLE(HpcPhotonPcmPairTrackIds, "AOD", "HPCPCMPAIRTRID",
                  hpc_photon_pcm_pair::PosJetTrack1Id, hpc_photon_pcm_pair::NegJetTrack1Id, hpc_photon_pcm_pair::PosJetTrack2Id, hpc_photon_pcm_pair::NegJetTrack2Id);
DECLARE_SOA_TABLE(HpcPhotonPcmPairV0PhotonIds, "AOD", "HPCPCMPAIRV0ID",
                  hpc_photon_pcm_pair::V0PhotonKF1Id, hpc_photon_pcm_pair::V0PhotonKF2Id);

// mc

// mcCollision extension
namespace hpc_mc_collision_extra
{
DECLARE_SOA_COLUMN(ThirdDecimalTruePosZ, thirdDecimalTruePosZ, int);
DECLARE_SOA_COLUMN(TrigEv, trigEv, bool);
DECLARE_SOA_COLUMN(PtMax, ptMax, float);
DECLARE_SOA_COLUMN(NChargedInEtaRange, nChargedInEtaRange, int);
} // namespace hpc_mc_collision_extra
DECLARE_SOA_TABLE(HpcMcCollisionsExtra, "AOD", "HPCMCCOLEXTRA",
                  hpc_common::SelEv, hpc_mc_collision_extra::ThirdDecimalTruePosZ, hpc_mc_collision_extra::TrigEv, hpc_mc_collision_extra::PtMax, hpc_mc_collision_extra::NChargedInEtaRange);

// trigger
namespace hpc_true_trigger
{
DECLARE_SOA_INDEX_COLUMN_FULL(JetMcParticle, jetMcParticle, int, JetParticles, "");
} // namespace hpc_true_trigger
DECLARE_SOA_TABLE(HpcTrueTriggers, "AOD", "HPCTRUETRIG",
                  o2::soa::Index<>, hpc_common::JetMcCollisionId, hpc_true_trigger::JetMcParticleId,
                  hpc_common::Pt, hpc_common::Phi, hpc_common::Eta);

} // namespace o2::aod

#endif // PWGJE_DATAMODEL_HADRONPHOTONCORRELATION_H_
