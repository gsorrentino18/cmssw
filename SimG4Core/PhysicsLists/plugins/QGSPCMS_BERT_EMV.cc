#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "QGSPCMS_BERT_EMV.h"

#include "G4DecayPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronicProcessStore.hh"
#include "G4IonPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4StoppingPhysics.hh"

#include "G4DataQuestionaire.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"

QGSPCMS_BERT_EMV::QGSPCMS_BERT_EMV(const edm::ParameterSet &p)
    : PhysicsList(p) {

  G4DataQuestionaire it(photon);

  int ver = p.getUntrackedParameter<int>("Verbosity", 0);
  bool emPhys = p.getUntrackedParameter<bool>("EMPhysics", true);
  bool hadPhys = p.getUntrackedParameter<bool>("HadPhysics", true);
  bool tracking = p.getParameter<bool>("TrackingCut");
  edm::LogInfo("PhysicsList") << "You are using the simulation engine: "
                              << "QGSP_BERT_EMV \n Flags for EM Physics "
                              << emPhys << ", for Hadronic Physics " << hadPhys
                              << " and tracking cut " << tracking;

  if (emPhys) {
    // EM Physics
    RegisterPhysics(new G4EmStandardPhysics_option1(ver));

    // Synchroton Radiation & GN Physics
    G4EmExtraPhysics *gn = new G4EmExtraPhysics(ver);
    RegisterPhysics(gn);
  }

  // Decays
  this->RegisterPhysics(new G4DecayPhysics(ver));

  if (hadPhys) {
    G4HadronicProcessStore::Instance()->SetVerbose(ver);

    // Hadron Elastic scattering
    RegisterPhysics(new G4HadronElasticPhysics(ver));

    // Hadron Physics
    RegisterPhysics(new G4HadronPhysicsQGSP_BERT(ver));

    // Stopping Physics
    RegisterPhysics(new G4StoppingPhysics(ver));

    // Ion Physics
    RegisterPhysics(new G4IonPhysics(ver));

    // Neutron tracking cut
    if (tracking) {
      RegisterPhysics(new G4NeutronTrackingCut(ver));
    }
  }
}
