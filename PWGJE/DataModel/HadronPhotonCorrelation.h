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
  None = -1,
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
  TrueAssocEvRecoPtTrig = 11,
  TrueAssocEv = 12,
  True = 13,
  TrueDecaySelfTrig = 14
};

template <typename>
bool constexpr AlwaysFalseType = false;
template <auto>
bool constexpr AlwaysFalseValue = false;

namespace o2::aod
{

// common correlation columns
namespace corr_common
{
DECLARE_SOA_INDEX_COLUMN_FULL(JetMcCollision, jetMcCollision, int, JMcCollisions, "");
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(PtTrue, ptTrue, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
} // namespace corr_common

// reco

// correlation collision
namespace corr_collision
{
// plain jet collision column for derived data
DECLARE_SOA_COLUMN(JetCollisionId, jetCollisionId, int);
DECLARE_SOA_COLUMN(SelEv, selEv, bool);
DECLARE_SOA_COLUMN(ThirdDecimalTruePosZ, thirdDecimalTruePosZ, int);
DECLARE_SOA_COLUMN(TrigEv, trigEv, bool);
DECLARE_SOA_COLUMN(PtMax, ptMax, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(NGlobalTracks, nGlobalTracks, int);
} // namespace corr_collision
DECLARE_SOA_TABLE(CorrCollisions, "AOD", "CORRCOLLISION",
                  o2::soa::Index<>,
                  corr_collision::SelEv, corr_collision::ThirdDecimalTruePosZ, corr_collision::TrigEv, corr_collision::PtMax, corr_collision::PosZ, corr_collision::NGlobalTracks);
DECLARE_SOA_TABLE(CorrCollisionIds, "AOD", "CORRCOLLISIONID",
                  corr_collision::JetCollisionId);
DECLARE_SOA_TABLE(CorrCollisionMcIds, "AOD", "CORRCOLLISIONMCID",
                  corr_common::JetMcCollisionId);

// add to common columns
namespace corr_common
{
DECLARE_SOA_INDEX_COLUMN_FULL(CorrCollision, corrCollision, int, CorrCollisions, "");
} // namespace corr_common

// trigger
namespace trigger
{
DECLARE_SOA_INDEX_COLUMN_FULL(JetTrack, jetTrack, int, JetTracks, "");
} // namespace trigger
DECLARE_SOA_TABLE(Triggers, "AOD", "TRIGGER",
                  o2::soa::Index<>, corr_common::CorrCollisionId,
                  corr_common::Pt, corr_common::Phi, corr_common::Eta);
DECLARE_SOA_TABLE(TriggerTrackIds, "AOD", "TRIGGERTRACKID",
                  trigger::JetTrackId);

// add to common columns
namespace corr_common
{
DECLARE_SOA_ARRAY_INDEX_COLUMN_FULL(IntersectTrigger, intersectTrigger, int, Triggers, "");
} // namespace corr_common

// hadrons (global tracks)
namespace hadron
{
DECLARE_SOA_INDEX_COLUMN_FULL(JetTrack, jetTrack, int, JetTracks, "");
DECLARE_SOA_COLUMN(TruthLevelRecoHadron, truthLevelRecoHadron, int);
} // namespace hadron
DECLARE_SOA_TABLE(Hadrons, "AOD", "HADRON",
                  o2::soa::Index<>, corr_common::CorrCollisionId, corr_common::IntersectTriggerIds,
                  corr_common::Pt, corr_common::Phi, corr_common::Eta,
                  hadron::TruthLevelRecoHadron, corr_common::PtTrue);
DECLARE_SOA_TABLE(HadronTrackIds, "AOD", "HADRONTRACKID",
                  hadron::JetTrackId);

// pipm
namespace pipm
{
DECLARE_SOA_INDEX_COLUMN_FULL(JetTrack, jetTrack, int, JetTracks, "");
DECLARE_SOA_COLUMN(TruthLevelRecoPipm, truthLevelRecoPipm, int);
} // namespace pipm
DECLARE_SOA_TABLE(Pipms, "AOD", "PIPM",
                  o2::soa::Index<>, corr_common::CorrCollisionId, corr_common::IntersectTriggerIds,
                  corr_common::Pt, corr_common::Phi, corr_common::Eta,
                  pipm::TruthLevelRecoPipm, corr_common::PtTrue);
DECLARE_SOA_TABLE(PipmTrackIds, "AOD", "PIPMTRACKID",
                  pipm::JetTrackId);

// photonPcm
namespace photon_pcm
{
DECLARE_SOA_INDEX_COLUMN_FULL(V0PhotonKF, v0PhotonKF, int, V0PhotonsKF, "");
DECLARE_SOA_INDEX_COLUMN_FULL(PosJetTrack, posJetTrack, int, JetTracks, "_pos");
DECLARE_SOA_INDEX_COLUMN_FULL(NegJetTrack, negJetTrack, int, JetTracks, "_neg");
DECLARE_SOA_COLUMN(TruthLevelRecoPhotonPcm, truthLevelRecoPhotonPcm, int);
} // namespace photon_pcm
DECLARE_SOA_TABLE(PhotonPcms, "AOD", "PHOTONPCM",
                  o2::soa::Index<>, corr_common::CorrCollisionId, corr_common::IntersectTriggerIds,
                  corr_common::Pt, corr_common::Phi, corr_common::Eta,
                  photon_pcm::TruthLevelRecoPhotonPcm, corr_common::PtTrue);
DECLARE_SOA_TABLE(PhotonPcmTrackIds, "AOD", "PHOTONPCMTRACKID",
                  photon_pcm::PosJetTrackId, photon_pcm::NegJetTrackId);
DECLARE_SOA_TABLE(PhotonPcmV0PhotonIds, "AOD", "PHOTONPCMV0PHOTONID",
                  photon_pcm::V0PhotonKFId);

// photonPcm pairs
namespace photon_pcm_pair
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
} // namespace photon_pcm_pair
DECLARE_SOA_TABLE(PhotonPcmPairs, "AOD", "PHOTONPCMPAIR",
                  o2::soa::Index<>, corr_common::CorrCollisionId, corr_common::IntersectTriggerIds,
                  corr_common::Pt, corr_common::Phi, corr_common::Eta, photon_pcm_pair::Mgg,
                  photon_pcm_pair::TruthLevelRecoPi0Pcm, photon_pcm_pair::TruthLevelRecoEtaPcm, corr_common::PtTrue);
DECLARE_SOA_TABLE(PhotonPcmPairTrackIds, "AOD", "PHOTONPCMPAIRTRACKID",
                  photon_pcm_pair::PosJetTrack1Id, photon_pcm_pair::NegJetTrack1Id, photon_pcm_pair::PosJetTrack2Id, photon_pcm_pair::NegJetTrack2Id);
DECLARE_SOA_TABLE(PhotonPcmPairV0PhotonIds, "AOD", "PHOTONPCMPAIRV0PHOTONID",
                  photon_pcm_pair::V0PhotonKF1Id, photon_pcm_pair::V0PhotonKF2Id);

// mc

// mcCollision extension
namespace mc_collision_extra_corr
{
DECLARE_SOA_COLUMN(ThirdDecimalTruePosZ, thirdDecimalTruePosZ, int);
DECLARE_SOA_COLUMN(TrigEv, trigEv, bool);
DECLARE_SOA_COLUMN(PtMax, ptMax, float);
DECLARE_SOA_COLUMN(NChargedInEtaRange, nChargedInEtaRange, int);
} // namespace mc_collision_extra_corr
DECLARE_SOA_TABLE(McCollisionsExtraCorr, "AOD", "MCCOLLISIONEXTRACORR",
                  mc_collision_extra_corr::ThirdDecimalTruePosZ, mc_collision_extra_corr::TrigEv, mc_collision_extra_corr::PtMax, mc_collision_extra_corr::NChargedInEtaRange);

// trigger
namespace true_trigger
{
DECLARE_SOA_INDEX_COLUMN_FULL(JetMcParticle, jetMcParticle, int, JetParticles, "");
} // namespace true_trigger
DECLARE_SOA_TABLE(TrueTriggers, "AOD", "TRUETRIGGER",
                  o2::soa::Index<>, corr_common::JetMcCollisionId, true_trigger::JetMcParticleId,
                  corr_common::Pt, corr_common::Phi, corr_common::Eta);
} // namespace o2::aod

#endif // PWGJE_DATAMODEL_HADRONPHOTONCORRELATION_H_
