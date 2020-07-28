#ifndef DVAnalyzer_Formats_TrackAssociation_H
#define DVAnalyzer_Formats_TrackAssociation_H

#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

namespace DVAna {
  typedef edm::AssociationMap<edm::OneToMany<reco::VertexCollection, reco::TrackCollection> > TrackAssociation;
}


#endif
