#include "ForceCollisionBiasingOperatorG4.hpp"

#include <limits>
#include <memory>

#include "G4BiasingProcessInterface.hh"
#include "G4BOptrForceCollision.hh"
#include "G4ParticleDefinition.hh"

#include "core/config/exceptions.h"
#include "core/utils/log.h"
#include "tools/geant4.h"

using namespace allpix;

ForceCollisionBiasingOperatorG4::ForceCollisionBiasingOperatorG4(G4ParticleDefinition* particleToBias) :
  G4VBiasingOperator("ForceCollision"),
  fParticleToBias(nullptr) {
    particleToBias = particleToBias;
    // Get name of particle to be biased
    std::string particleName = particleToBias->GetParticleName();

    // Create ForceCollision operator
    optr = new G4BOptrForceCollision(particleName,
        "ForceCollisionFor" + particleName);
}

ForceCollisionBiasingOperatorG4::~ForceCollisionBiasingOperatorG4() {
}


G4VBiasingOperation* 
ForceCollisionBiasingOperatorG4::ProposeOccurenceBiasingOperation(const G4Track* track, 
    const G4BiasingProcessInterface* callingProcess) {
    if(fCurrentOperator) {
        return fCurrentOperator->GetProposedOccurenceBiasingOperation(track, callingProcess);
    } else {
        return nullptr;
    }
}

G4VBiasingOperation*
ForceCollisionBiasingOperatorG4::ProposeNonPhysicsBiasingOperation(const G4Track* track,
    const G4BiasingProcessInterface* callingProcess) {
    if(fCurrentOperator) {
        return fCurrentOperator->GetProposedNonPhysicsBiasingOperation(track, callingProcess);
    } else {
        return nullptr;
    }
}

G4VBiasingOperation*
ForceCollisionBiasingOperatorG4::ProposeFinalStateBiasingOperation(const G4Track* track,
    const G4BiasingProcessInterface* callingProcess) {
    if(fCurrentOperator) {
        return fCurrentOperator->GetProposedFinalStateBiasingOperation(track, callingProcess);
    } else {
        return nullptr;
    }
}

void ForceCollisionBiasingOperatorG4::StartTracking(const G4Track*) {
    fCurrentOperator = optr;
}

void ForceCollisionBiasingOperatorG4::
OperationApplied( const G4BiasingProcessInterface*         callingProcess,
                  G4BiasingAppliedCase                        biasingCase,
                  G4VBiasingOperation*                   operationApplied,
                  const G4VParticleChange*         particleChangeProduced ) {
    if(fCurrentOperator) {
        fCurrentOperator->ReportOperationApplied( callingProcess,
            biasingCase,
            operationApplied,
            particleChangeProduced );
    }
}

void ForceCollisionBiasingOperatorG4::
OperationApplied( const G4BiasingProcessInterface*        callingProcess,
                  G4BiasingAppliedCase                       biasingCase,
                  G4VBiasingOperation*         occurenceOperationApplied,
                  G4double                 weightForOccurenceInteraction,
                  G4VBiasingOperation*        finalStateOperationApplied, 
                  const G4VParticleChange*        particleChangeProduced ) {
    if(fCurrentOperator) {
        fCurrentOperator->ReportOperationApplied( callingProcess,
            biasingCase,
            occurenceOperationApplied,
            weightForOccurenceInteraction,
            finalStateOperationApplied, 
            particleChangeProduced );
    }
}

void ForceCollisionBiasingOperatorG4::
ExitBiasing( const G4Track*                           track,
             const G4BiasingProcessInterface* callingProcess ) {
    if(fCurrentOperator) {
        fCurrentOperator->ExitingBiasing(track, callingProcess);
    }
}
