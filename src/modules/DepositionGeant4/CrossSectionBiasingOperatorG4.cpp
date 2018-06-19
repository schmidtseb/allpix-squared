#include "CrossSectionBiasingOperatorG4.hpp"

#include <limits>
#include <memory>

#include "G4BiasingProcessInterface.hh"
#include "G4BOptnChangeCrossSection.hh"
#include "G4VProcess.hh"
#include "G4InteractionLawPhysical.hh"

#include "core/config/exceptions.h"
#include "core/utils/log.h"
#include "tools/geant4.h"

using namespace allpix;

CrossSectionBiasingOperatorG4::CrossSectionBiasingOperatorG4(G4ParticleDefinition* particleToBias, std::string name) : 
	G4VBiasingOperator(name),
	fParticleToBias(particleToBias),
	fSetup(true) {
}

CrossSectionBiasingOperatorG4::~CrossSectionBiasingOperatorG4() {
}

void CrossSectionBiasingOperatorG4::StartRun() {
	if(fSetup) {
		const G4ProcessManager* processManager = fParticleToBias->GetProcessManager();
		const G4BiasingProcessSharedData* sharedData =
			G4BiasingProcessInterface::GetSharedData( processManager );
		if( sharedData ) {
			for(auto& processInterfaces : sharedData->GetPhysicsBiasingProcessInterfaces()) {
				const G4BiasingProcessInterface* wrapperProcess = processInterfaces;
				std::string operationName = "XSchange-" + wrapperProcess->GetWrappedProcess()->GetProcessName();
				fChangeCrossSectionOperations[wrapperProcess] = new G4BOptnChangeCrossSection(operationName);
			}
		}
		fSetup = false;
	}
}

G4VBiasingOperation*
CrossSectionBiasingOperatorG4::ProposeOccurenceBiasingOperation(const G4Track* track, const G4BiasingProcessInterface* callingProcess) {
	// Check if current particle type is the one to bias
	if(track->GetDefinition() != fParticleToBias) {
		return nullptr;
	}

	double analogInteractionLength = callingProcess->GetWrappedProcess()->GetCurrentInteractionLength();
	if(analogInteractionLength > DBL_MAX/10.) {
		return nullptr;
	}

	double analogXS = 1./analogInteractionLength;
	double XStransformation = 10.;

	G4BOptnChangeCrossSection* operation = fChangeCrossSectionOperations[callingProcess];
	G4VBiasingOperation* previousOperation = callingProcess->GetPreviousOccurenceBiasingOperation();

	if(previousOperation == nullptr) {
		operation->SetBiasedCrossSection( XStransformation* analogXS );
		operation->Sample();
	} else {
		if(previousOperation != operation) {
			return nullptr;
		}
		if(operation->GetInteractionOccured()) {
			operation->SetBiasedCrossSection(XStransformation* analogXS);
			operation->Sample();
		} else {
			operation->UpdateForStep(callingProcess->GetPreviousStepSize());
			operation->SetBiasedCrossSection(XStransformation* analogXS);
			operation->UpdateForStep(0.);
		}
	}

	return operation;
}

void CrossSectionBiasingOperatorG4::
	OperationApplied(const G4BiasingProcessInterface*callingProcess, 
	    G4BiasingAppliedCase,
	    G4VBiasingOperation* occurenceOperationApplied,
	    G4double,
	    G4VBiasingOperation*,
	    const G4VParticleChange*) {
	G4BOptnChangeCrossSection* operation =
		fChangeCrossSectionOperations[callingProcess];
	if(operation == occurenceOperationApplied) {
		operation->SetInteractionOccured();
	}
}
