// Thermal-neutron capture on one isotopically pure lead target with NuDEX enabled (QGSP_BERT_HP).
// Counts capture steps and the gammas each capture emits, and prints the yield per capture of every
// gamma line (1 keV bins) above a threshold.
//
// Usage: nudex_lines A N_EVENTS SEED [MIN_YIELD]
//   A         target mass number: 204, 206 or 207 (any stable Pb isotope works)
//   N_EVENTS  number of 0.0253 eV neutrons started at the centre of a 6 m lead cube
//   SEED      seed for the CLHEP engine (MixMaxRng)
//   MIN_YIELD print only bins with at least this many gammas per capture (default 0.003)
// Requires G4NUDEXLIBDATA to point to an extracted G4NUDEXLIB data set.

#include "G4Box.hh"
#include "G4Element.hh"
#include "G4Gamma.hh"
#include "G4HadronicParameters.hh"
#include "G4Isotope.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Neutron.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4RunManagerFactory.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4UserSteppingAction.hh"
#include "G4VProcess.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "QGSP_BERT_HP.hh"
#include "Randomize.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>

namespace {

int gA = 207;
long gCaptures = 0;
long gCapturesOther = 0;  // captures on any nucleus other than the target isotope (should stay 0)
std::map<int, long> gGammaPerKeV;  // gamma count per 1 keV bin of kinetic energy

class Detector : public G4VUserDetectorConstruction {
 public:
  G4VPhysicalVolume* Construct() override {
    auto* iso = new G4Isotope("Pb" + std::to_string(gA), 82, gA);
    auto* el = new G4Element("PbEnriched" + std::to_string(gA), "Pb", 1);
    el->AddIsotope(iso, 1.0);
    auto* mat = new G4Material("Pb" + std::to_string(gA) + "_pure", 11.35 * g / cm3, 1);
    mat->AddElement(el, 1.0);
    auto* box = new G4Box("World", 3 * m, 3 * m, 3 * m);
    auto* lv = new G4LogicalVolume(box, mat, "World");
    return new G4PVPlacement(nullptr, {}, lv, "World", nullptr, false, 0);
  }
};

class Gun : public G4VUserPrimaryGeneratorAction {
 public:
  Gun() : fGun(1) {
    fGun.SetParticleDefinition(G4Neutron::Definition());
    fGun.SetParticleEnergy(0.0253 * eV);
    fGun.SetParticlePosition({});
  }
  void GeneratePrimaries(G4Event* ev) override {
    const double cost = 2 * G4UniformRand() - 1, phi = CLHEP::twopi * G4UniformRand();
    const double sint = std::sqrt(1 - cost * cost);
    fGun.SetParticleMomentumDirection({sint * std::cos(phi), sint * std::sin(phi), cost});
    fGun.GeneratePrimaryVertex(ev);
  }

 private:
  G4ParticleGun fGun;
};

class Stepping : public G4UserSteppingAction {
 public:
  void UserSteppingAction(const G4Step* step) override {
    const G4VProcess* p = step->GetPostStepPoint()->GetProcessDefinedStep();
    // QGSP_BERT_HP names the capture process "nCaptureHP"; other lists use "nCapture".
    if (!p || (p->GetProcessName() != "nCaptureHP" && p->GetProcessName() != "nCapture")) return;
    const auto* sec = step->GetSecondaryInCurrentStep();
    bool product = false;
    for (const G4Track* t : *sec) {
      const auto* d = t->GetDefinition();
      if (d->GetAtomicNumber() == 82 && d->GetAtomicMass() == gA + 1) product = true;
    }
    if (!product) { ++gCapturesOther; return; }
    ++gCaptures;
    for (const G4Track* t : *sec)
      if (t->GetDefinition() == G4Gamma::Definition())
        ++gGammaPerKeV[static_cast<int>(std::floor(t->GetKineticEnergy() / keV))];
  }
};

class Actions : public G4VUserActionInitialization {
 public:
  void Build() const override {
    SetUserAction(new Gun);
    SetUserAction(new Stepping);
  }
};

}  // namespace

int main(int argc, char** argv) {
  if (argc < 4) {
    std::fprintf(stderr, "usage: %s A N_EVENTS SEED [MIN_YIELD]\n", argv[0]);
    return 2;
  }
  gA = std::atoi(argv[1]);
  const long n = std::atol(argv[2]);
  const long seed = std::atol(argv[3]);
  const double minYield = argc > 4 ? std::atof(argv[4]) : 0.003;

  G4Random::setTheEngine(new CLHEP::MixMaxRng);
  G4Random::setTheSeed(seed);

  G4HadronicParameters::Instance()->SetEnableNUDEX(true);  // must precede the physics list
  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  rm->SetUserInitialization(new Detector);
  rm->SetUserInitialization(new QGSP_BERT_HP(0));
  rm->SetUserInitialization(new Actions);
  rm->Initialize();
  rm->BeamOn(static_cast<G4int>(n));

  std::printf("# target Pb-%d  events %ld  seed %ld  captures %ld  other-captures %ld\n", gA, n, seed,
              gCaptures, gCapturesOther);
  std::printf("# gammas per capture in 1 keV bins with yield >= %g (bin low edge in keV)\n", minYield);
  std::printf("bin_keV\tcount\tper_capture\tstat_err\n");
  for (const auto& [bin, c] : gGammaPerKeV) {
    const double y = gCaptures ? double(c) / gCaptures : 0;
    if (y >= minYield) std::printf("%d\t%ld\t%.5f\t%.5f\n", bin, c, y, std::sqrt(double(c)) / gCaptures);
  }
  delete rm;
  return 0;
}
