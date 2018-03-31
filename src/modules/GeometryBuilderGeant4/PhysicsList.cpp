#include "PhysicsList.hpp"

// General
#include "G4ParticleDefinition.hh"
#include "G4PhysicsListHelper.hh"

// Particles
#include "G4Gamma.hh"
#include "G4Alpha.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4GenericIon.hh"
#include "G4IonConstructor.hh"

// Decay Physics
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecay.hh"
#include "G4RadioactiveDecayPhysics.hh"

// EM Processes
#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4RayleighScattering.hh"
#include "G4EmProcessOptions.hh"

// Deexcitation
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"

PhysicsList::PhysicsList(bool enableDecay) :
	fEnableDecay(enableDecay) {
}

PhysicsList::~PhysicsList() {

}

void PhysicsList::ConstructParticle() {
	G4Gamma::GammaDefinition();
	G4Alpha::AlphaDefinition();
	G4Electron::ElectronDefinition();
	G4Positron::PositronDefinition();

	G4IonConstructor iConstructor;
	iConstructor.ConstructParticle();
}

void PhysicsList::ConstructProcess() {
	AddTransportation();
	ConstructEM();
	if (fEnableDecay) {
		ConstructDecay();
		ConstructDeexcite();
	}
}

void PhysicsList::ConstructEM() {
	G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
	G4ParticleDefinition* particle = G4Gamma::GammaDefinition();

	// Construct and register processes for gamma
	ph->RegisterProcess(new G4PhotoElectricEffect(), particle);
	ph->RegisterProcess(new G4ComptonScattering(), particle);
	ph->RegisterProcess(new G4GammaConversion(), particle);
	ph->RegisterProcess(new G4RayleighScattering(), particle);

	// EM Options
	G4EmProcessOptions emOptions;
	emOptions.SetMscStepLimitation(fUseSafety);
}

void PhysicsList::ConstructDecay() {
	G4RadioactiveDecay* radioactiveDecay = new G4RadioactiveDecay();
	radioactiveDecay->SetHLThreshold(1.);

	radioactiveDecay->SetICM(false);
	radioactiveDecay->SetARM(false
		);
	radioactiveDecay->SetAnalogueMonteCarlo(true);

	G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
	ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());
}

void PhysicsList::ConstructDeexcite() {
	G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
	de->SetFluo(true);
	de->SetAuger(false);
	de->SetPIXE(false);
	de->InitialiseAtomicDeexcitation();
	G4LossTableManager::Instance()->SetAtomDeexcitation(de);
}

void PhysicsList::SetCuts() {
	// SetCutValue(1*m, "gamma");
}
