#ifndef TRACKREFMAP_H
#define TRACKREFMAP_H

#include <map>
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

namespace DVAna{
  typedef std::pair<reco::CandidatePtr, reco::TrackRef> TrackRefPair;
  typedef std::map<reco::CandidatePtr, reco::TrackRef> UnpackedCandidateTracksMap;
  typedef std::map<reco::TrackRef, reco::TrackRef> RescaledTrackMap;
}

//typedef std::pair<edm::Ptr<reco::Candidate>,edm::Ref<std::vector<reco::Track>,reco::Track,edm::refhelper::FindUsingAdvance<std::vector<reco::Track>,reco::Track> > > TrackRefPair;
//typedef std::map<edm::Ptr<reco::Candidate>,edm::Ref<std::vector<reco::Track>,reco::Track,edm::refhelper::FindUsingAdvance<std::vector<reco::Track>,reco::Track> > > UnpackedCandidateTracksMap;

#endif
