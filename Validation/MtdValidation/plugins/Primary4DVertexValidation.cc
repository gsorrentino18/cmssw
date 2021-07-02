#include <numeric>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Common/interface/ValidHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"

// reco track and vertex
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

// TrackingParticle
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

// pile-up
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// associator
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/VertexAssociation/interface/calculateVertexSharedTracks.h"

// vertexing
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"

// simulated vertex
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociator.h"

// DQM
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#define verbose_ false
#define use_only_charged_tracks_ true
#define NOT_MATCHED 66666
#define debug_ true

//class declaration
class Primary4DVertexValidation : public DQMEDAnalyzer {
  typedef math::XYZTLorentzVector LorentzVector;

  //enum SignalVertexKind { HIGHEST_PT = 0, IS_ASSOC2FIRST_RECO = 1, IS_ASSOC2ANY_RECO = 2 };

  // auxiliary class holding simulated vertices
  struct simPrimaryVertex {
    simPrimaryVertex(double x1, double y1, double z1, double t1)
        : x(x1),
          y(y1),
          z(z1),
          t(t1),
          ptsq(0),
          closest_vertex_distance_z(-1.),
          nGenTrk(0),
          num_matched_reco_tracks(0),
          average_match_quality(0.0) {
      ptot.setPx(0);
      ptot.setPy(0);
      ptot.setPz(0);
      ptot.setE(0);
      p4 = LorentzVector(0, 0, 0, 0);
      r = sqrt(x * x + y * y);
    };
    double x, y, z, r, t;
    HepMC::FourVector ptot;
    LorentzVector p4;
    double ptsq;
    double closest_vertex_distance_z;
    int nGenTrk;
    int num_matched_reco_tracks;
    float average_match_quality;
    EncodedEventId eventId;
    TrackingVertexRef sim_vertex;

    unsigned int nwosmatch = 0;  // number of recvertices dominated by this simevt (by wos)
    unsigned int nwntmatch = 0;  // number of recvertices dominated by this simevt  (by nt)
    std::vector<unsigned int> wos_dominated_recv;  // list of dominated recv (by wos, size==nwosmatch)

    std::map<unsigned int, double> wnt;     // weighted number of tracks in recvtx (by index)
    std::map<unsigned int, double> wos;     // sum of wos in recvtx (by index) // oops -> this was int before 04-22
    double sumwos;                          // sum of wos in any recvtx
    double sumwnt;                           // sum of weighted tracks
    unsigned int rec;                    // best match  (NO_MATCH if not matched)
    unsigned int matchQuality;           // quality flag


    void addTrack(unsigned int irecv, double twos, double twt) {
      sumwnt += twt;
      if (wnt.find(irecv) == wnt.end()) {
        wnt[irecv] = twt;
      } else {
        wnt[irecv] += twt;
      }

      sumwos += twos;
      if (wos.find(irecv) == wos.end()) {
        wos[irecv] = twos;
      } else {
        wos[irecv] += twos;
      }
    }

  };

  // auxiliary class holding reconstructed vertices
  struct recoPrimaryVertex {
    //enum VertexProperties { NONE = 0, MATCHED = 1, DUPLICATE = 2, MERGED = 4 };
    recoPrimaryVertex(double x1, double y1, double z1)
        : x(x1),
          y(y1),
          z(z1),
          pt(0),
          ptsq(0),
          closest_vertex_distance_z(-1.),
          purity(-1.),
          nRecoTrk(0),
          num_matched_sim_tracks(0),
          kind_of_vertex(0),
          recVtx(nullptr) {
      r = sqrt(x * x + y * y);
    };
    double x, y, z, r;
    double pt;
    double ptsq;
    double closest_vertex_distance_z;
    double purity;  // calculated and assigned in calculatePurityAndFillHistograms
    int nRecoTrk;
    int num_matched_sim_tracks;
    int kind_of_vertex;
    const reco::Vertex *recVtx;
    reco::VertexBaseRef recVtxRef;
    
    std::map<unsigned int, double> wos;  // simevent -> wos
    std::map<unsigned int, double> wnt;  // simevent -> weighted number of truth matched tracks
    unsigned int wosmatch;               // index of the simevent providing the largest contribution to wos
    unsigned int wntmatch;               // index of the simevent providing the highest number of tracks
    double sumwos;                       // total sum of wos of all truth matched tracks
    double sumwnt;                       // total weighted number of truth matchted tracks
    double maxwos;                       // largest wos sum from one sim event (wosmatch)
    double maxwnt;                       // largest weighted  number of tracks from one sim event (ntmatch)
    int maxwosnt;                        // number of tracks from the simevt with highest wos
    unsigned int sim;      // best match  (NO_MATCH if not matched)
    unsigned int matchQuality;           // quality flag

    void addTrack(unsigned int iev, double twos, double wt) {
      sumwnt += wt;
      if (wnt.find(iev) == wnt.end()) {
        wnt[iev] = wt;
      } else {
        wnt[iev] += wt;
      }

      sumwos += twos;
      if (wos.find(iev) == wos.end()) {
        wos[iev] = twos;
      } else {
        wos[iev] += twos;
      }
    }
  };

public:
  explicit Primary4DVertexValidation(const edm::ParameterSet &);
  ~Primary4DVertexValidation() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void analyze(const edm::Event &, const edm::EventSetup &) override;
  void bookHistograms(DQMStore::IBooker &i, edm::Run const &, edm::EventSetup const &) override;

private:
  void matchReco2Sim(std::vector<recoPrimaryVertex> &,
                             std::vector<simPrimaryVertex> &,
                             const edm::ValueMap<float> &,
                             const edm::Handle<reco::BeamSpot>&);
  bool matchRecoTrack2SimSignal(const reco::TrackBaseRef &);
  bool matchRecoTrack2Sim(const reco::TrackBaseRef &, const TrackingVertexRef &);
  std::vector<Primary4DVertexValidation::simPrimaryVertex> getSimPVs(
      const edm::Handle<TrackingVertexCollection> &);
  std::vector<Primary4DVertexValidation::recoPrimaryVertex> getRecoPVs(
      const edm::Handle<edm::View<reco::Vertex>> &);

  // ----------member data ---------------------------
 
  const std::string folder_; 
  const double simUnit_; //sim time in s while reco time in ns
  const reco::RecoToSimCollection *r2s_;
  const reco::SimToRecoCollection *s2r_;

  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> vecPileupSummaryInfoToken_;

  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleCollectionToken_;
  edm::EDGetTokenT<TrackingVertexCollection> trackingVertexCollectionToken_;
  edm::EDGetTokenT<reco::SimToRecoCollection> simToRecoAssociationToken_;
  edm::EDGetTokenT<reco::RecoToSimCollection> recoToSimAssociationToken_;
  edm::EDGetTokenT<reco::BeamSpot> RecBeamSpotToken_;
  edm::EDGetTokenT<edm::View<reco::Vertex>> Rec4DVerToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0SafePidToken_;
 
  //histogram declaration  
  MonitorElement* meTimeRes_;
  MonitorElement* meTimePull_;
};

// constructors and destructor
Primary4DVertexValidation::Primary4DVertexValidation(const edm::ParameterSet& iConfig) 
    :  folder_(iConfig.getParameter<std::string>("folder")),
       simUnit_(iConfig.getParameter<double>("simUnit")) {
      vecPileupSummaryInfoToken_ = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag(std::string("addPileupInfo"))),
      trackingParticleCollectionToken_ = consumes<TrackingParticleCollection>(
          iConfig.getParameter<edm::InputTag>("SimTag")),
      trackingVertexCollectionToken_ = 
          consumes<TrackingVertexCollection>(iConfig.getParameter<edm::InputTag>("SimTag")),
      simToRecoAssociationToken_ = 
          consumes<reco::SimToRecoCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc")),
      recoToSimAssociationToken_ = 
          consumes<reco::RecoToSimCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc")),
      RecBeamSpotToken_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("offlineBS")),
      Rec4DVerToken_ = consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("offline4DPV")),
      Sigmat0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0SafePID"));
}

Primary4DVertexValidation::~Primary4DVertexValidation() {}

//
// member functions
//
void Primary4DVertexValidation::bookHistograms(DQMStore::IBooker& ibook,
                                                     edm::Run const& iRun,
                                                     edm::EventSetup const& iSetup) {

   ibook.setCurrentFolder(folder_);
   // --- histograms booking
   //
   //
   meTimeRes_ = ibook.book1D("TimeRes", "t_{rec} - t_{sim}; ", 200, -1., 1.);
   meTimePull_ = ibook.book1D("TimePull", "t_{rec} - t_{sim}/#sigma_{t rec}; ", 200, -10., 10.);
}


bool Primary4DVertexValidation::matchRecoTrack2SimSignal(const reco::TrackBaseRef& recoTrack) {
  auto found = r2s_->find(recoTrack);

  // reco track not matched to any TP
  if (found == r2s_->end())
    return false;

  //// reco track matched to some TP from signal vertex
  for (const auto& tp : found->val) {
    if (tp.first->eventId().bunchCrossing() == 0 && tp.first->eventId().event() == 0)
      return true;
  }

  // reco track not matched to any TP from signal vertex
  return false;
}

bool Primary4DVertexValidation::matchRecoTrack2Sim(const reco::TrackBaseRef& recoTrack, const TrackingVertexRef &vsim) {
  auto found = r2s_->find(recoTrack);

  // reco track not matched to any TP
  if (found == r2s_->end())
    return false;

  //matched TP equal to any TP of sim vertex 
  for (const auto& tp : found->val) {
    if (std::find_if(vsim->daughterTracks_begin(), vsim->daughterTracks_end(), [&](const TrackingParticleRef &vtp) {return tp.first == vtp;}) != vsim->daughterTracks_end())
    return true;
  }

  // reco track not matched to any TP from signal vertex
  return false;
}

/* Extract information form TrackingParticles/TrackingVertex and fill
 * the helper class simPrimaryVertex with proper generation-level
 * information */
std::vector<Primary4DVertexValidation::simPrimaryVertex> Primary4DVertexValidation::getSimPVs(
    const edm::Handle<TrackingVertexCollection>& tVC) {

  std::vector<Primary4DVertexValidation::simPrimaryVertex> simpv;
  int current_event = -1;

  for (TrackingVertexCollection::const_iterator v = tVC->begin(); v != tVC->end(); ++v) {
    //We keep only the first vertex from all the events at BX=0.
    if (v->eventId().bunchCrossing() != 0)
      continue;
    if (v->eventId().event() != current_event) {
      current_event = v->eventId().event();
    } else {
      continue;
    }
    // TODO(rovere) is this really necessary?
    if (fabs(v->position().z()) > 1000)
      continue;  // skip funny junk vertices

    // could be a new vertex, check  all primaries found so far to avoid
    // multiple entries
    simPrimaryVertex sv(v->position().x(), v->position().y(), v->position().z(), v->position().t());
    sv.eventId = v->eventId();
    sv.sim_vertex = TrackingVertexRef(tVC, std::distance(tVC->begin(), v));

    for (TrackingParticleRefVector::iterator iTrack = v->daughterTracks_begin(); iTrack != v->daughterTracks_end();
         ++iTrack) {
      // TODO(rovere) isn't it always the case? Is it really worth
      // checking this out?
      // sv.eventId = (**iTrack).eventId();
      assert((**iTrack).eventId().bunchCrossing() == 0);
    }
    // TODO(rovere) maybe get rid of this old logic completely ... ?
    simPrimaryVertex* vp = nullptr;  // will become non-NULL if a vertex
                                     // is found and then point to it
    for (std::vector<simPrimaryVertex>::iterator v0 = simpv.begin(); v0 != simpv.end(); v0++) {
      if ((sv.eventId == v0->eventId) && (fabs(sv.x - v0->x) < 1e-5) && (fabs(sv.y - v0->y) < 1e-5) &&
          (fabs(sv.z - v0->z) < 1e-5)) {
        vp = &(*v0);
        break;
      }
    }
    if (!vp) {
      // this is a new vertex, add it to the list of sim-vertices
      simpv.push_back(sv);
      vp = &simpv.back();
    }

    // Loop over daughter track(s) as Tracking Particles
    for (TrackingVertex::tp_iterator iTP = v->daughterTracks_begin(); iTP != v->daughterTracks_end(); ++iTP) {
      auto momentum = (*(*iTP)).momentum();
      const reco::Track* matched_best_reco_track = nullptr;
      double match_quality = -1;
      if (use_only_charged_tracks_ && (**iTP).charge() == 0)
        continue;
      if (s2r_->find(*iTP) != s2r_->end()) {
        matched_best_reco_track = (*s2r_)[*iTP][0].first.get();
        match_quality = (*s2r_)[*iTP][0].second;
      }
      if (verbose_) {
        std::cout << "  Daughter momentum:      " << momentum;
        std::cout << "  Daughter type     " << (*(*iTP)).pdgId();
        std::cout << "  matched: " << (matched_best_reco_track != nullptr);
        std::cout << "  match-quality: " << match_quality;
        std::cout << std::endl;
      }
      vp->ptot.setPx(vp->ptot.x() + momentum.x());
      vp->ptot.setPy(vp->ptot.y() + momentum.y());
      vp->ptot.setPz(vp->ptot.z() + momentum.z());
      vp->ptot.setE(vp->ptot.e() + (**iTP).energy());
      vp->ptsq += ((**iTP).pt() * (**iTP).pt());
      // TODO(rovere) only select charged sim-particles? If so, maybe
      // put it as a configuration parameter?
      if (matched_best_reco_track) {
        vp->num_matched_reco_tracks++;
        vp->average_match_quality += match_quality;
      }
      // TODO(rovere) get rid of cuts on sim-tracks
      // TODO(rovere) be consistent between simulated tracks and
      // reconstructed tracks selection
      // count relevant particles
      //if (((**iTP).pt() > 0.2) && (fabs((**iTP).eta()) < 2.5) && (**iTP).charge() != 0) {
        //vp->nGenTrk++;
      //}
    }  // End of for loop on daughters sim-particles
    if (vp->num_matched_reco_tracks)
      vp->average_match_quality /= static_cast<float>(vp->num_matched_reco_tracks);
    if (verbose_) {
      std::cout << "average number of associated tracks: "
                << vp->num_matched_reco_tracks / static_cast<float>(vp->nGenTrk)
                << " with average quality: " << vp->average_match_quality << std::endl;
    }
  }  // End of for loop on tracking vertices

  // In case of no simulated vertices, break here
  if (simpv.empty())
    return simpv;

  // Now compute the closest distance in z between all simulated vertex
  // first initialize
  auto prev_z = simpv.back().z;
  for (simPrimaryVertex& vsim : simpv) {
    vsim.closest_vertex_distance_z = std::abs(vsim.z - prev_z);
    prev_z = vsim.z;
  }
  // then calculate
  for (std::vector<simPrimaryVertex>::iterator vsim = simpv.begin(); vsim != simpv.end(); vsim++) {
    std::vector<simPrimaryVertex>::iterator vsim2 = vsim;
    vsim2++;
    for (; vsim2 != simpv.end(); vsim2++) {
      double distance = std::abs(vsim->z - vsim2->z);
      // need both to be complete
      vsim->closest_vertex_distance_z = std::min(vsim->closest_vertex_distance_z, distance);
      vsim2->closest_vertex_distance_z = std::min(vsim2->closest_vertex_distance_z, distance);
    }
  }
  return simpv;
}

/* Extract information form recoVertex and fill the helper class
 * recoPrimaryVertex with proper reco-level information */
std::vector<Primary4DVertexValidation::recoPrimaryVertex> Primary4DVertexValidation::getRecoPVs(
    const edm::Handle<edm::View<reco::Vertex>>& tVC) {


  std::vector<Primary4DVertexValidation::recoPrimaryVertex> recopv;
  for (auto v = tVC->begin(); v != tVC->end(); ++v) {
    // Skip junk vertices
    if (fabs(v->z()) > 1000)
      continue;
    if (v->isFake() || !v->isValid())
      continue;

    recoPrimaryVertex sv(v->position().x(), v->position().y(), v->position().z());
    sv.recVtx = &(*v);
    sv.recVtxRef = reco::VertexBaseRef(tVC, std::distance(tVC->begin(), v));

    // this is a new vertex, add it to the list of reco-vertices
    recopv.push_back(sv);
    Primary4DVertexValidation::recoPrimaryVertex* vp = &recopv.back();

    // Loop over daughter track(s)
    for (auto iTrack = v->tracks_begin(); iTrack != v->tracks_end(); ++iTrack) {
      auto momentum = (*(*iTrack)).innerMomentum();
      // TODO(rovere) better handle the pixelVertices, whose tracks
      // do not have the innerMomentum defined. This is a temporary
      // hack to overcome this problem.
      if (momentum.mag2() == 0)
        momentum = (*(*iTrack)).momentum();
      vp->pt += std::sqrt(momentum.perp2());
      vp->ptsq += (momentum.perp2());
      vp->nRecoTrk++;

      auto matched = r2s_->find(*iTrack);
      if (matched != r2s_->end()) {
        vp->num_matched_sim_tracks++;
      }

    }  // End of for loop on daughters reconstructed tracks
  }    // End of for loop on tracking vertices

  // In case of no reco vertices, break here
  if (recopv.empty())
    return recopv;

  // Now compute the closest distance in z between all reconstructed vertex
  // first initialize
  auto prev_z = recopv.back().z;
  for (recoPrimaryVertex& vreco : recopv) {
    vreco.closest_vertex_distance_z = std::abs(vreco.z - prev_z);
    prev_z = vreco.z;
  }
  for (std::vector<recoPrimaryVertex>::iterator vreco = recopv.begin(); vreco != recopv.end(); vreco++) {
    std::vector<recoPrimaryVertex>::iterator vreco2 = vreco;
    vreco2++;
    for (; vreco2 != recopv.end(); vreco2++) {
      double distance = std::abs(vreco->z - vreco2->z);
      // need both to be complete
      vreco->closest_vertex_distance_z = std::min(vreco->closest_vertex_distance_z, distance);
      vreco2->closest_vertex_distance_z = std::min(vreco2->closest_vertex_distance_z, distance);
    }
  }
  return recopv;
}

// ------------ method called to produce the data  ------------
void Primary4DVertexValidation::matchReco2Sim(std::vector<recoPrimaryVertex>& recopv,
                                                            std::vector<simPrimaryVertex>& simpv,
                                                            const edm::ValueMap<float>& ValueMap,
                                                            const edm::Handle<reco::BeamSpot>& BS) {


  for (auto vv: simpv) {
    vv.wnt.clear();
    vv.wos.clear();
  }
  for (auto rv: recopv) {
    rv.wnt.clear();
    rv.wos.clear();
  }

  for (unsigned int iv = 0; iv < recopv.size(); iv++) {
    const reco::Vertex* vertex = recopv.at(iv).recVtx;
 
    for (unsigned int iev = 0; iev < simpv.size(); iev++) {
      
      double wnt=0;
      double wos=0;
      double evwnt=0;
      double evwos=0;
      double evnt=0;

      for (auto iTrack = vertex->tracks_begin(); iTrack != vertex->tracks_end(); ++iTrack) {
        double pt = (*iTrack)->pt();
        if(vertex->trackWeight(*iTrack) < 0.5) continue;

        if (matchRecoTrack2Sim(*iTrack, simpv.at(iev).sim_vertex)) {
           
           double dz2_beam = pow((*BS).BeamWidthX() * cos((*iTrack)->phi()) / tan((*iTrack)->theta()), 2) + pow((*BS).BeamWidthY() * sin((*iTrack)->phi()) / tan((*iTrack)->theta()), 2);
           double dz2 = pow((*iTrack)->dz(), 2) + dz2_beam + pow(0.0020, 2); // added 20 um, some tracks have crazy small resolutions
           wos = vertex->trackWeight(*iTrack)/dz2;
           wnt = vertex->trackWeight(*iTrack) * std::min(pt, 1.0);

           //std::cout << "Sigmat0Safe[*iTrack]: " << ValueMap[(*iTrack)] << std::endl;
           //if ((*iTrack)->covt0t0() > 0) { //FIXME
           if (ValueMap[(*iTrack)] > 0) {
             double sigmaZ = (*BS).sigmaZ();
             double sigmaT = sigmaZ / 2.998e1;  // c in cm/ns
             wos = wos / erf(ValueMap[(*iTrack)]/sigmaT);
           }
           simpv.at(iev).addTrack(iv, wos, wnt);
           recopv.at(iv).addTrack(iev, wos, wnt);
           evwos += wos;
           evwnt += wnt;
           evnt++;
        }
      } //RecoTracks loop

      // require 2 tracks for a wos-match
      if ((evwos > 0) && (evwos > recopv.at(iv).maxwos) && (evnt > 1)) {
        recopv.at(iv).wosmatch = iev;
        recopv.at(iv).maxwos = evwos;
        recopv.at(iv).maxwosnt = evnt;

        simpv.at(iev).wos_dominated_recv.push_back(iv);
        simpv.at(iev).nwosmatch++;
      }

      // weighted track counting match, require at least one track
      if ((evnt > 0) && (evwnt > recopv.at(iv).maxwnt)) {
        recopv.at(iv).wntmatch = iev;
        recopv.at(iv).maxwnt = evwnt;
      }
    } //TrackingVertex loop

  } //RecoPrimaryVertex


  //after filling infos, goes for the sim-reco match 
  double zWosMatchMax = 1;
  for (auto& vrec : recopv) {
    vrec.sim = NOT_MATCHED;
    vrec.matchQuality = 0;
  }
  unsigned int iev = 0;
  for (auto& vv : simpv) {
    if (debug_) {
      std::cout << "iev: " << iev << std::endl;
      std::cout << "wos_dominated_recv.size: " << vv.wos_dominated_recv.size() << std::endl;
    }
    for (unsigned int i=0; i < vv.wos_dominated_recv.size(); i++) {
       auto recov = vv.wos_dominated_recv.at(i);
       if (debug_) {std::cout << "index of reco vertex: " << recov << " that has a wos: " << vv.wos.at(recov) << " at position " << i <<std::endl;}
    }
    vv.rec = NOT_MATCHED;
    vv.matchQuality = 0;
    iev++;
  }
  //this tries a one-to-one match, taking simPV with highest wos if there are > 1 simPV candidates
  for (unsigned int rank = 1; rank < 8; rank++)  {
      std::cout << "RANK: " << rank << std::endl;

      for (unsigned int iev = 0; iev < simpv.size(); iev++) { //loop on SimPV
          std::cout << "*****iev: " << iev << std::endl;
          std::cout << "simpv.at(iev).rec: " << simpv.at(iev).rec << std::endl;
	  if (simpv.at(iev).rec != NOT_MATCHED) continue;
	  if (simpv.at(iev).nwosmatch == 0) continue;
	  if (simpv.at(iev).nwosmatch > rank) continue;
	  unsigned int iv = NOT_MATCHED;
	  for (unsigned int k = 0; k < simpv.at(iev).wos_dominated_recv.size(); k++) {
	      unsigned int rec = simpv.at(iev).wos_dominated_recv.at(k);
              auto vrec = recopv.at(rec);
	      if (vrec.sim != NOT_MATCHED) continue; // already matched
	      if (fabs(simpv.at(iev).z - vrec.z) > zWosMatchMax) continue;// insanely far away
	      if ( (iv == NOT_MATCHED) || simpv.at(iev).wos.at(rec) > simpv.at(iev).wos.at(iv)) {
                iv = rec;
              }
	  }
	  if (iv != NOT_MATCHED) { //if the rec vertex has already been associated is possible that iv remains NOT_MATCHED at this point 
	      recopv.at(iv).sim = iev;
	      simpv.at(iev).rec = iv;
              recopv.at(iv).matchQuality = rank;
	      simpv.at(iev).matchQuality = rank;
	  }
      }
  }
  //give vertices a chance that have a lot of overlap, but are still recognizably
  //caused by a specific simvertex (without being classified as dominating)
  //like a small peak sitting on the flank of a larger nearby peak
  unsigned int ntry = 0;
  while (ntry++ < 10)
    {
      unsigned nmatch = 0;
      for (unsigned int iev = 0; iev < simpv.size(); iev++)
	{
	  if ((simpv.at(iev).rec != NOT_MATCHED) || (simpv.at(iev).wos.size() == 0))
	    continue;
	  // ok, single simvertex sim, who is your your favorite rec vertex?
	  unsigned int rec = NOT_MATCHED;
	  for (auto rv : simpv.at(iev).wos) {
	    if ((rec == NOT_MATCHED) || (rv.second > simpv.at(iev).wos.at(rec))) {
	      rec = rv.first;
	    }
	  }

	  if (rec == NOT_MATCHED){ //try with wnt match --> siamo sicuri di doverlo fare????
	    for (auto rv : simpv.at(iev).wnt) { 
	      if ((rec == NOT_MATCHED) || (rv.second > simpv.at(iev).wnt.at(rec))) {
		rec = rv.first;
	      }
	    }
	  }
	
	if (rec == NOT_MATCHED)
	  continue;  // should not happen
	if (recopv.at(rec).sim != NOT_MATCHED)
	  continue;  // already gone
	
	// do you, recvertex rec, take this simvertex sim as your lawful wedded truthmatch?
        unsigned int rec2sim = NOT_MATCHED;
	for (auto sv : recopv.at(rec).wos) {
	  if (simpv.at(sv.first).rec != NOT_MATCHED)
	    continue;  // already used
	  if ((rec2sim == NOT_MATCHED) || (sv.second > recopv.at(rec).wos.at(rec2sim))) {
	    rec2sim = sv.first;
	  }
	}
	if (iev == rec2sim) {
	  // I do
	  recopv.at(rec).sim = iev;
	  recopv.at(rec).matchQuality = 8;
	  simpv.at(iev).rec = rec;
	  simpv.at(iev).matchQuality = 8;
	  nmatch++;
	}
      }  //sim loop
      if (nmatch == 0) {
	break;
      }
    }  // ntry 
}

void Primary4DVertexValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using edm::Handle;
  using edm::View;
  using std::cout;
  using std::endl;
  using std::vector;
  using namespace reco;

  std::vector<float> pileUpInfo_z;

  // get the pileup information
  edm::Handle<std::vector<PileupSummaryInfo>> puinfoH;
  if (iEvent.getByToken(vecPileupSummaryInfoToken_, puinfoH)) {
    for (auto const& pu_info : *puinfoH.product()) {
      if (pu_info.getBunchCrossing() == 0) {
        pileUpInfo_z = pu_info.getPU_zpositions();
        if (verbose_) {
          for (auto const& p : pileUpInfo_z) {
            std::cout << "PileUpInfo on Z vertex: " << p << std::endl;
          }
        }
        break;
      }
    }
  }

  edm::Handle<TrackingParticleCollection> TPCollectionH;
  iEvent.getByToken(trackingParticleCollectionToken_, TPCollectionH);
  if (!TPCollectionH.isValid())
    edm::LogWarning("Primary4DVertexValidation") << "TPCollectionH is not valid";

  edm::Handle<TrackingVertexCollection> TVCollectionH;
  iEvent.getByToken(trackingVertexCollectionToken_, TVCollectionH);
  if (!TVCollectionH.isValid())
    edm::LogWarning("Primary4DVertexValidation") << "TVCollectionH is not valid";

  edm::Handle<reco::SimToRecoCollection> simToRecoH;
  iEvent.getByToken(simToRecoAssociationToken_, simToRecoH);
  if (simToRecoH.isValid())
    s2r_ = simToRecoH.product();
  else
    edm::LogWarning("Primary4DVertexValidation") << "simToRecoH is not valid";

  edm::Handle<reco::RecoToSimCollection> recoToSimH;
  iEvent.getByToken(recoToSimAssociationToken_, recoToSimH);
  if (recoToSimH.isValid())
    r2s_ = recoToSimH.product();
  else
    edm::LogWarning("Primary4DVertexValidation") << "recoToSimH is not valid";

  edm::Handle<reco::BeamSpot> BeamSpotH;
  iEvent.getByToken(RecBeamSpotToken_, BeamSpotH);
  if (!BeamSpotH.isValid())
    edm::LogWarning("Primary4DVertexValidation") << "BeamSpotH is not valid";

  std::vector<simPrimaryVertex> simpv;  // a list of simulated primary MC vertices
  simpv = getSimPVs(TVCollectionH);
  // TODO(rovere) 1 vertex is not, by definition, pileup, and should
  // probably be subtracted?
  int num_pileup_vertices = simpv.size();
  //this bool check if first vertex in that with highest pT (which means is the signal vertex?)
  bool signal_is_highest_pt =
      std::max_element(simpv.begin(), simpv.end(), [](const simPrimaryVertex& lhs, const simPrimaryVertex& rhs) {
        return lhs.ptsq < rhs.ptsq;
      }) == simpv.begin();

  std::vector<recoPrimaryVertex> recopv;  // a list of reconstructed primary MC vertices
  edm::Handle<edm::View<reco::Vertex>> recVtxs;
  iEvent.getByToken(Rec4DVerToken_, recVtxs);
  if (!recVtxs.isValid())
    edm::LogWarning("Primary4DVertexValidation") << "recVtxs is not valid";
  recopv = getRecoPVs(recVtxs);

  const auto& Sigmat0Safe = iEvent.get(Sigmat0SafePidToken_);

  // check upfront that refs to track are (likely) to be valid
  /*bool ok = true;
  for (const auto& v : *recVtxs) {
    if (v.tracksSize() > 0) {
      const auto& ref = v.trackRefAt(0);
      if (ref.isNull() || !ref.isAvailable()) {
        ok = false;
      }
    }
  }
  if (!ok)  continue; //use this bool to skip tracks whose ref is not valid
  */

  //I have simPV and recoPV collections 
  matchReco2Sim(recopv, simpv, Sigmat0Safe, BeamSpotH);

  //fill histograms here in a new loop
  
  if (debug_) {
    std::cout << "simpv.size: " << simpv.size() << std::endl;
    std::cout << "recopv.size: " << recopv.size() << std::endl;
  }
  for (unsigned int is=0; is<simpv.size(); is++) {
     if (simpv.at(is).rec == NOT_MATCHED) {
       if (debug_) {std::cout << "sim vertex: " << is << " is not matched with any reco" << std::endl;}
       continue;
     }
     for (unsigned int ir=0; ir<recopv.size(); ir++) {
        if (simpv.at(is).rec == ir) {
        }
        if (recopv.at(ir).sim == is && simpv.at(is).rec == ir) {
           meTimeRes_->Fill(recopv.at(ir).recVtx->t()-simpv.at(is).t * simUnit_);
           meTimePull_->Fill((recopv.at(ir).recVtx->t()-simpv.at(is).t * simUnit_)/recopv.at(ir).recVtx->tError());

           if (debug_) {
             std::cout << "*** Matching RECO: " << ir << "with SIM: " << is << " ***"<<std::endl;
             std::cout << "Pull is: " << (recopv.at(ir).recVtx->t()-simpv.at(is).t * simUnit_)/recopv.at(ir).recVtx->tError() << std:: endl;
             std::cout << "Match Quality is " << recopv.at(ir).matchQuality << std::endl;
             //std::cout << "Reco Vertex time: " << recopv.at(ir).recVtx->t() << "+- " << recopv.at(ir).recVtx->tError() << std::endl;
             //std::cout << "Sim Vertex time: " << simpv.at(is).t * simUnit_ << std::endl;
             std::cout << "****" << std::endl;
           }
        }  
     } 
  }


}  // end of analyze

void Primary4DVertexValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/dummy");
  desc.add<edm::InputTag>("TPtoRecoTrackAssoc", edm::InputTag("trackingParticleRecoTrackAsssociation"));
  desc.add<edm::InputTag>("SimTag", edm::InputTag("mix","MergedTrackTruth"));
  desc.add<edm::InputTag>("offlineBS", edm::InputTag("offlineBeamSpot"));
  desc.add<edm::InputTag>("offline4DPV", edm::InputTag("offlinePrimaryVertices4D"));
  desc.add<edm::InputTag>("sigmat0SafePID", edm::InputTag("tofPID:sigmat0safe"));

  desc.add<double>("simUnit", 1e9); 
  descriptions.add("vertices4D", desc);
}

DEFINE_FWK_MODULE(Primary4DVertexValidation);
