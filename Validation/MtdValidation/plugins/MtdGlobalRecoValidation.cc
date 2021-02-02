#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/Common/interface/ValidHandle.h"
#include "DataFormats/Math/interface/GeantUnits.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/Records/interface/MTDTopologyRcd.h"
#include "Geometry/MTDNumberingBuilder/interface/MTDTopology.h"
#include "Geometry/MTDCommonData/interface/MTDTopologyMode.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

class MtdGlobalRecoValidation : public DQMEDAnalyzer {
public:
  explicit MtdGlobalRecoValidation(const edm::ParameterSet&);
  ~MtdGlobalRecoValidation() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;

  // ------------ member data ------------

  const std::string folder_;
  const float trackMinEnergy_;
  const float trackMinEta_;
  const float trackMaxEta_;

  edm::EDGetTokenT<reco::TrackCollection> GenRecTrackToken_;
  edm::EDGetTokenT<reco::TrackCollection> MTDRecTrackToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> RecVertexToken_;

  edm::EDGetTokenT<edm::ValueMap<int>> trackAssocToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> t0Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmat0Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0safeToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> sigmat0safeToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probPToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> mtdQualMVAToken_;

  MonitorElement* meBTLTrackRPTime_;
  MonitorElement* meBTLTrackRPTimeErr_;
  MonitorElement* meBTLTrackRPBestTimeErr_;
  MonitorElement* meBTLTrackEffEtaTot_;
  MonitorElement* meBTLTrackEffPhiTot_;
  MonitorElement* meBTLTrackEffPtTot_;
  MonitorElement* meBTLTrackEffEtaMtd_;
  MonitorElement* meBTLTrackEffPhiMtd_;
  MonitorElement* meBTLTrackEffPtMtd_;

  MonitorElement* meETLTrackRPTime_[4];
  MonitorElement* meETLTrackRPTimeErr_[4];
  MonitorElement* meETLTrackRPBestTimeErr_[4];
  MonitorElement* meETLTrackNumHits_[4];
  MonitorElement* meETLTrackEffEtaTot_[2];
  MonitorElement* meETLTrackEffPhiTot_[2];
  MonitorElement* meETLTrackEffPtTot_[2];
  MonitorElement* meETLTrackEffEtaMtd_[2];
  MonitorElement* meETLTrackEffPhiMtd_[2];
  MonitorElement* meETLTrackEffPtMtd_[2];

  MonitorElement* meTrackNumHits_;

  MonitorElement* meVerNumber_;
  MonitorElement* meVerZ_;
  MonitorElement* meVerTime_;
};

// ------------ constructor and destructor --------------
MtdGlobalRecoValidation::MtdGlobalRecoValidation(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")),
      trackMinEnergy_(iConfig.getParameter<double>("trackMinimumEnergy")),
      trackMinEta_(iConfig.getParameter<double>("trackMinimumEta")),
      trackMaxEta_(iConfig.getParameter<double>("trackMaximumEta")) {
  GenRecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagG"));
  MTDRecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagT"));
  RecVertexToken_ = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("inputTagV"));
  trackAssocToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("trackAssocSrc"));
  t0Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0Src"));
  sigmat0Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0Src"));
  t0safeToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0safeSrc"));
  sigmat0safeToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0safeSrc"));
  probPiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probPiSrc"));
  probKToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probKSrc"));
  probPToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probPSrc"));
  mtdQualMVAToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("mtdQualMVASrc"));
}

MtdGlobalRecoValidation::~MtdGlobalRecoValidation() {}

// ------------ method called for each event  ------------
void MtdGlobalRecoValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace geant_units::operators;
  using namespace std;

  edm::ESHandle<MTDTopology> topologyHandle;
  iSetup.get<MTDTopologyRcd>().get(topologyHandle);
  const MTDTopology* topology = topologyHandle.product();

  bool topo1Dis = false;
  bool topo2Dis = false;
  if (topology->getMTDTopologyMode() <= static_cast<int>(MTDTopologyMode::Mode::barphiflat)) {
    topo1Dis = true;
  }
  if (topology->getMTDTopologyMode() > static_cast<int>(MTDTopologyMode::Mode::barphiflat)) {
    topo2Dis = true;
  }

  auto GenRecTrackHandle = makeValid(iEvent.getHandle(GenRecTrackToken_));
  //auto MTDRecTrackHandle = makeValid(iEvent.getHandle(MTDRecTrackToken_));
  auto RecVertexHandle = makeValid(iEvent.getHandle(RecVertexToken_));

  const auto& trackAssoc = iEvent.get(trackAssocToken_);

  const auto& t0 = iEvent.get(t0Token_);
  const auto& sigmat0 = iEvent.get(sigmat0Token_);
  const auto& t0safe = iEvent.get(t0safeToken_);
  const auto& sigmat0safe = iEvent.get(sigmat0safeToken_);
  const auto& probPi = iEvent.get(probPiToken_);
  const auto& probK = iEvent.get(probKToken_);
  const auto& probP = iEvent.get(probPToken_);
  const auto& mtdQualMVA = iEvent.get(mtdQualMVAToken_);

  unsigned int index = 0;
  // --- Loop over all RECO tracks ---
  for (const auto& trackGen : *GenRecTrackHandle) {
    const reco::TrackRef trackref(iEvent.getHandle(GenRecTrackToken_), index);
    index++;

    if (trackAssoc[trackref] == -1) {
      LogWarning("globalReco") << "Extended track not associated";
      continue;
    }

    const reco::TrackRef mtdTrackref = reco::TrackRef(iEvent.getHandle(MTDRecTrackToken_), trackAssoc[trackref]);
    const reco::Track track = *mtdTrackref;

    if (track.pt() < trackMinEnergy_)
      continue;

    if (fabs(track.eta()) < trackMinEta_) {
      // --- all BTL tracks (with and without hit in MTD) ---
      meBTLTrackEffEtaTot_->Fill(track.eta());
      meBTLTrackEffPhiTot_->Fill(track.phi());
      meBTLTrackEffPtTot_->Fill(track.pt());

      bool MTDBtl = false;
      int numMTDBtlvalidhits = 0;
      for (const auto hit : track.recHits()) {
        if (hit->isValid() == false)
          continue;
        MTDDetId Hit = hit->geographicalId();
        if ((Hit.det() == 6) && (Hit.subdetId() == 1) && (Hit.mtdSubDetector() == 1)) {
          MTDBtl = true;
          numMTDBtlvalidhits++;
        }
      }
      meTrackNumHits_->Fill(numMTDBtlvalidhits);

      // --- keeping only tracks with last hit in MTD ---
      if (MTDBtl == true) {
        if (sigmat0safe[trackref] > 0.0251) {
          edm::LogWarning("globalReco") << "pt " << track.pt() << " genpt " << trackGen.pt() << " mva "
                                        << mtdQualMVA[trackref] << " t0/st0/t0safe/st0safe " << t0[trackref] << " "
                                        << sigmat0[trackref] << " " << t0safe[trackref] << " " << sigmat0safe[trackref]
                                        << " prob(pi/K/p) " << probPi[trackref] << " " << probK[trackref] << " "
                                        << probP[trackref];
        }
        meBTLTrackEffEtaMtd_->Fill(track.eta());
        meBTLTrackEffPhiMtd_->Fill(track.phi());
        meBTLTrackEffPtMtd_->Fill(track.pt());
        meBTLTrackRPTime_->Fill(track.t0());
        meBTLTrackRPTimeErr_->Fill(track.t0Error());
        meBTLTrackRPBestTimeErr_->Fill(sigmat0safe[trackref]);
      }
    }

    else {
      // --- all ETL tracks (with and without hit in MTD) ---
      if ((track.eta() < -trackMinEta_) && (track.eta() > -trackMaxEta_)) {
        meETLTrackEffEtaTot_[0]->Fill(track.eta());
        meETLTrackEffPhiTot_[0]->Fill(track.phi());
        meETLTrackEffPtTot_[0]->Fill(track.pt());
      }

      if ((track.eta() > trackMinEta_) && (track.eta() < trackMaxEta_)) {
        meETLTrackEffEtaTot_[1]->Fill(track.eta());
        meETLTrackEffPhiTot_[1]->Fill(track.phi());
        meETLTrackEffPtTot_[1]->Fill(track.pt());
      }

      bool MTDEtlZnegD1 = false;
      bool MTDEtlZnegD2 = false;
      bool MTDEtlZposD1 = false;
      bool MTDEtlZposD2 = false;
      int numMTDEtlvalidhits = 0;
      for (const auto hit : track.recHits()) {
        if (hit->isValid() == false)
          continue;
        MTDDetId Hit = hit->geographicalId();
        if ((Hit.det() == 6) && (Hit.subdetId() == 1) && (Hit.mtdSubDetector() == 2)) {
          ETLDetId ETLHit = hit->geographicalId();

          if (topo2Dis) {
            if ((ETLHit.zside() == -1) && (ETLHit.nDisc() == 1)) {
              MTDEtlZnegD1 = true;
              meETLTrackRPTime_[0]->Fill(track.t0());
              meETLTrackRPTimeErr_[0]->Fill(track.t0Error());
              meETLTrackRPBestTimeErr_[0]->Fill(sigmat0safe[trackref]);
              numMTDEtlvalidhits++;
            }
            if ((ETLHit.zside() == -1) && (ETLHit.nDisc() == 2)) {
              MTDEtlZnegD2 = true;
              meETLTrackRPTime_[1]->Fill(track.t0());
              meETLTrackRPTimeErr_[1]->Fill(track.t0Error());
              meETLTrackRPBestTimeErr_[1]->Fill(sigmat0safe[trackref]);
              numMTDEtlvalidhits++;
            }
            if ((ETLHit.zside() == 1) && (ETLHit.nDisc() == 1)) {
              MTDEtlZposD1 = true;
              meETLTrackRPTime_[2]->Fill(track.t0());
              meETLTrackRPTimeErr_[2]->Fill(track.t0Error());
              meETLTrackRPBestTimeErr_[2]->Fill(sigmat0safe[trackref]);
              numMTDEtlvalidhits++;
            }
            if ((ETLHit.zside() == 1) && (ETLHit.nDisc() == 2)) {
              MTDEtlZposD2 = true;
              meETLTrackRPTime_[3]->Fill(track.t0());
              meETLTrackRPTimeErr_[3]->Fill(track.t0Error());
              meETLTrackRPBestTimeErr_[3]->Fill(sigmat0safe[trackref]);
              numMTDEtlvalidhits++;
            }
          }

          if (topo1Dis) {
            if (ETLHit.zside() == -1) {
              MTDEtlZnegD1 = true;
              meETLTrackRPTime_[0]->Fill(track.t0());
              meETLTrackRPTimeErr_[0]->Fill(track.t0Error());
              numMTDEtlvalidhits++;
            }
            if (ETLHit.zside() == 1) {
              MTDEtlZposD1 = true;
              meETLTrackRPTime_[2]->Fill(track.t0());
              meETLTrackRPTimeErr_[2]->Fill(track.t0Error());
              numMTDEtlvalidhits++;
            }
          }
        }
      }
      meTrackNumHits_->Fill(-numMTDEtlvalidhits);

      // --- keeping only tracks with last hit in MTD ---
      if ((track.eta() < -trackMinEta_) && (track.eta() > -trackMaxEta_)) {
        if ((MTDEtlZnegD1 == true) || (MTDEtlZnegD2 == true)) {
          if (sigmat0safe[trackref] > 0.0251) {
            edm::LogWarning("globalReco")
                << "pt " << track.pt() << " genpt " << trackGen.pt() << " mva " << mtdQualMVA[trackref]
                << " t0/st0/t0safe/st0safe " << t0[trackref] << " " << sigmat0[trackref] << " " << t0safe[trackref]
                << " " << sigmat0safe[trackref] << " prob(pi/K/p) " << probPi[trackref] << " " << probK[trackref] << " "
                << probP[trackref];
          }
          meETLTrackEffEtaMtd_[0]->Fill(track.eta());
          meETLTrackEffPhiMtd_[0]->Fill(track.phi());
          meETLTrackEffPtMtd_[0]->Fill(track.pt());
        }
      }
      if ((track.eta() > trackMinEta_) && (track.eta() < trackMaxEta_)) {
        if ((MTDEtlZposD1 == true) || (MTDEtlZposD2 == true)) {
          if (sigmat0safe[trackref] > 0.0251) {
            edm::LogWarning("globalReco")
                << "pt " << track.pt() << " genpt " << trackGen.pt() << " mva " << mtdQualMVA[trackref]
                << " t0/st0/t0safe/st0safe " << t0[trackref] << " " << sigmat0[trackref] << " " << t0safe[trackref]
                << " " << sigmat0safe[trackref] << " prob(pi/K/p) " << probPi[trackref] << " " << probK[trackref] << " "
                << probP[trackref];
          }
          meETLTrackEffEtaMtd_[1]->Fill(track.eta());
          meETLTrackEffPhiMtd_[1]->Fill(track.phi());
          meETLTrackEffPtMtd_[1]->Fill(track.pt());
        }
      }
    }
  }  //RECO tracks loop

  // --- Loop over the RECO vertices ---
  int nv = 0;
  for (const auto& v : *RecVertexHandle) {
    if (v.isValid()) {
      meVerZ_->Fill(v.z());
      meVerTime_->Fill(v.t());
      nv++;
    } else
      cout << "The vertex is not valid" << endl;
  }
  meVerNumber_->Fill(nv);
}

// ------------ method for histogram booking ------------
void MtdGlobalRecoValidation::bookHistograms(DQMStore::IBooker& ibook,
                                             edm::Run const& run,
                                             edm::EventSetup const& iSetup) {
  ibook.setCurrentFolder(folder_);

  // histogram booking
  meBTLTrackRPTime_ = ibook.book1D("TrackBTLRPTime", "Track t0 with respect to R.P.;t0 [ns]", 100, -1, 3);
  meBTLTrackRPTimeErr_ = ibook.book1D("TrackBTLRPTimeErr", "Track t0Error with respect to R.P.;t0 [ns]", 100, 0, 0.5);
  meBTLTrackRPBestTimeErr_ =
      ibook.book1D("TrackBTLRPBestTimeErr", "Track sigmat0safe with respect to R.P.;t0 [ns]", 100, 0, 0.5);
  meBTLTrackEffEtaTot_ = ibook.book1D("TrackBTLEffEtaTot", "Track efficiency vs eta (Tot);#eta_{RECO}", 100, -1.6, 1.6);
  meBTLTrackEffPhiTot_ =
      ibook.book1D("TrackBTLEffPhiTot", "Track efficiency vs phi (Tot);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meBTLTrackEffPtTot_ = ibook.book1D("TrackBTLEffPtTot", "Track efficiency vs pt (Tot);pt_{RECO} [GeV]", 50, 0, 10);
  meBTLTrackEffEtaMtd_ = ibook.book1D("TrackBTLEffEtaMtd", "Track efficiency vs eta (Mtd);#eta_{RECO}", 100, -1.6, 1.6);
  meBTLTrackEffPhiMtd_ =
      ibook.book1D("TrackBTLEffPhiMtd", "Track efficiency vs phi (Mtd);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meBTLTrackEffPtMtd_ = ibook.book1D("TrackBTLEffPtMtd", "Track efficiency vs pt (Mtd);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackRPTime_[0] =
      ibook.book1D("TrackETLRPTimeZnegD1", "Track t0 with respect to R.P. (-Z, First Disk);t0 [ns]", 100, -1, 3);
  meETLTrackRPTime_[1] =
      ibook.book1D("TrackETLRPTimeZnegD2", "Track t0 with respect to R.P. (-Z, Second Disk);t0 [ns]", 100, -1, 3);
  meETLTrackRPTime_[2] =
      ibook.book1D("TrackETLRPTimeZposD1", "Track t0 with respect to R.P. (+Z, First Disk);t0 [ns]", 100, -1, 3);
  meETLTrackRPTime_[3] =
      ibook.book1D("TrackETLRPTimeZposD2", "Track t0 with respect to R.P. (+Z, Second Disk);t0 [ns]", 100, -1, 3);
  meETLTrackRPTimeErr_[0] = ibook.book1D(
      "TrackETLRPTimeErrZnegD1", "Track t0Error with respect to R.P. (-Z, First Disk);t0 [ns]", 100, 0, 0.5);
  meETLTrackRPTimeErr_[1] = ibook.book1D(
      "TrackETLRPTimeErrZnegD2", "Track t0Error with respect to R.P. (-Z, Second Disk);t0 [ns]", 100, 0, 0.5);
  meETLTrackRPTimeErr_[2] = ibook.book1D(
      "TrackETLRPTimeErrZposD1", "Track t0Error with respect to R.P. (+Z, First Disk);t0 [ns]", 100, 0, 0.5);
  meETLTrackRPTimeErr_[3] = ibook.book1D(
      "TrackETLRPTimeErrZposD2", "Track t0Error with respect to R.P. (+Z, Second Disk);t0 [ns]", 100, 0, 0.5);
  meETLTrackRPBestTimeErr_[0] = ibook.book1D(
      "TrackETLRPBestTimeErrZnegD1", "Track sigmat0safe with respect to R.P. (-Z, First Disk);t0 [ns]", 100, 0, 0.5);
  meETLTrackRPBestTimeErr_[1] = ibook.book1D(
      "TrackETLRPBestTimeErrZnegD2", "Track sigmat0safe with respect to R.P. (-Z, Second Disk);t0 [ns]", 100, 0, 0.5);
  meETLTrackRPBestTimeErr_[2] = ibook.book1D(
      "TrackETLRPBestTimeErrZposD1", "Track sigmat0safe with respect to R.P. (+Z, First Disk);t0 [ns]", 100, 0, 0.5);
  meETLTrackRPBestTimeErr_[3] = ibook.book1D(
      "TrackETLRPBestTimeErrZposD2", "Track sigmat0safe with respect to R.P. (+Z, Second Disk);t0 [ns]", 100, 0, 0.5);
  meETLTrackEffEtaTot_[0] =
      ibook.book1D("TrackETLEffEtaTotZneg", "Track efficiency vs eta (Tot) (-Z);#eta_{RECO}", 100, -3.2, -1.4);
  meETLTrackEffEtaTot_[1] =
      ibook.book1D("TrackETLEffEtaTotZpos", "Track efficiency vs eta (Tot) (+Z);#eta_{RECO}", 100, 1.4, 3.2);
  meETLTrackEffPhiTot_[0] =
      ibook.book1D("TrackETLEffPhiTotZneg", "Track efficiency vs phi (Tot) (-Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPhiTot_[1] =
      ibook.book1D("TrackETLEffPhiTotZpos", "Track efficiency vs phi (Tot) (+Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPtTot_[0] =
      ibook.book1D("TrackETLEffPtTotZneg", "Track efficiency vs pt (Tot) (-Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackEffPtTot_[1] =
      ibook.book1D("TrackETLEffPtTotZpos", "Track efficiency vs pt (Tot) (+Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackEffEtaMtd_[0] =
      ibook.book1D("TrackETLEffEtaMtdZneg", "Track efficiency vs eta (Mtd) (-Z);#eta_{RECO}", 100, -3.2, -1.4);
  meETLTrackEffEtaMtd_[1] =
      ibook.book1D("TrackETLEffEtaMtdZpos", "Track efficiency vs eta (Mtd) (+Z);#eta_{RECO}", 100, 1.4, 3.2);
  meETLTrackEffPhiMtd_[0] =
      ibook.book1D("TrackETLEffPhiMtdZneg", "Track efficiency vs phi (Mtd) (-Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPhiMtd_[1] =
      ibook.book1D("TrackETLEffPhiMtdZpos", "Track efficiency vs phi (Mtd) (+Z);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackEffPtMtd_[0] =
      ibook.book1D("TrackETLEffPtMtdZneg", "Track efficiency vs pt (Mtd) (-Z);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackEffPtMtd_[1] =
      ibook.book1D("TrackETLEffPtMtdZpos", "Track efficiency vs pt (Mtd) (+Z);pt_{RECO} [GeV]", 50, 0, 10);
  meTrackNumHits_ = ibook.book1D("TrackNumHits", "Number of valid MTD hits per track ; Number of hits", 10, -5, 5);
  meVerZ_ = ibook.book1D("VerZ", "RECO Vertex Z;Z_{RECO} [cm]", 180, -18, 18);
  meVerTime_ = ibook.book1D("VerTime", "RECO Vertex Time;t0 [ns]", 100, -1, 1);
  meVerNumber_ = ibook.book1D("VerNumber", "RECO Vertex Number: Number of vertices", 100, 0, 500);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void MtdGlobalRecoValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/GlobalReco");
  desc.add<edm::InputTag>("inputTagG", edm::InputTag("generalTracks", ""));
  desc.add<edm::InputTag>("inputTagT", edm::InputTag("trackExtenderWithMTD", ""));
  desc.add<edm::InputTag>("inputTagV", edm::InputTag("offlinePrimaryVertices4D", ""));
  desc.add<edm::InputTag>("trackAssocSrc", edm::InputTag("trackExtenderWithMTD", "generalTrackassoc"))
      ->setComment("Association between General and MTD Extended tracks");
  desc.add<double>("trackMinimumEnergy", 1.0);  // [GeV]
  desc.add<double>("trackMinimumEta", 1.5);
  desc.add<double>("trackMaximumEta", 3.2);
  desc.add<edm::InputTag>("t0Src", edm::InputTag("tofPID", "t0"));
  desc.add<edm::InputTag>("sigmat0Src", edm::InputTag("tofPID", "sigmat0"));
  desc.add<edm::InputTag>("t0safeSrc", edm::InputTag("tofPID", "t0safe"));
  desc.add<edm::InputTag>("sigmat0safeSrc", edm::InputTag("tofPID", "sigmat0safe"));
  desc.add<edm::InputTag>("probPiSrc", edm::InputTag("tofPID", "probPi"));
  desc.add<edm::InputTag>("probKSrc", edm::InputTag("tofPID", "probK"));
  desc.add<edm::InputTag>("probPSrc", edm::InputTag("tofPID", "probP"));
  desc.add<edm::InputTag>("mtdQualMVASrc", edm::InputTag("mtdTrackQualityMVA", "mtdQualMVA"));

  descriptions.add("globalReco", desc);
}

DEFINE_FWK_MODULE(MtdGlobalRecoValidation);
