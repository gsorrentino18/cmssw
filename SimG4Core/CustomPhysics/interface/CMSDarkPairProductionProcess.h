//--------------------------------------------------------
// File name:     CMSDarkPairProductionProcess
//
//  Author:        Dustin Stolp (dostolp@ucdavis.edu)
//                 Sushil S. Chauhan (schauhan@cern.ch)
// --------------------------------------------------------
#ifndef SimG4Core_CustomPhysics_CMSDarkPairProductionProcess_h
#define SimG4Core_CustomPhysics_CMSDarkPairProductionProcess_h

#include "G4Gamma.hh"
#include "G4VEmProcess.hh"
#include "SimG4Core/CustomPhysics/interface/CMSDarkPairProduction.h"
#include "globals.hh"

class G4ParticleDefinition;
class G4VEmModel;
class G4MaterialCutsCouple;
class G4DynamicParticle;

class CMSDarkPairProductionProcess : public G4VEmProcess

{
public: // with description
  CMSDarkPairProductionProcess(G4double df = 1E0,
                               const G4String &processName = "conv",
                               G4ProcessType type = fElectromagnetic);

  ~CMSDarkPairProductionProcess() override;

  // true for Gamma only.
  G4bool IsApplicable(const G4ParticleDefinition &) override;

  G4double MinPrimaryEnergy(const G4ParticleDefinition *,
                            const G4Material *) override;

  // Print few lines of informations about the process: validity range,
  void PrintInfo() override;

protected:
  void InitialiseProcess(const G4ParticleDefinition *) override;

private:
  G4bool isInitialised;
  G4double darkFactor;
};

#endif
