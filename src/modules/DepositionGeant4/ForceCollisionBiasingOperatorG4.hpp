#ifndef ALLPIX_FORCE_COLLISION_BIASING_OPERATOR_G4_H
#define ALLPIX_FORCE_COLLISION_BIASING_OPERATOR_G4_H

#include "G4VBiasingOperator.hh"
#include <map>

#include "core/config/Configuration.hpp"

class G4BOptrForceCollision;
class G4ParticleDefinition;

namespace allpix {
    class ForceCollisionBiasingOperatorG4 : public G4VBiasingOperator {
    public:
        explicit ForceCollisionBiasingOperatorG4(G4ParticleDefinition* particleToBias);
        ~ForceCollisionBiasingOperatorG4();

        virtual void StartTracking( const G4Track* track ) final;

    private:

        virtual G4VBiasingOperation*
        ProposeNonPhysicsBiasingOperation(const G4Track* track,
                                        const G4BiasingProcessInterface* callingProcess) final;
        virtual G4VBiasingOperation* 
        ProposeOccurenceBiasingOperation(const G4Track* track,
                                        const G4BiasingProcessInterface* callingProcess) final;
        virtual G4VBiasingOperation*
        ProposeFinalStateBiasingOperation(const G4Track* track,
                                        const G4BiasingProcessInterface* callingProcess) final;

        // Needed implementions to forward calls to the underneath
        // G4BOptrForceCollision biasing operator
        void OperationApplied( const G4BiasingProcessInterface*         callingProcess,
                                G4BiasingAppliedCase                        biasingCase,
                                G4VBiasingOperation*                   operationApplied,
                                const G4VParticleChange*         particleChangeProduced ) final;
        void OperationApplied( const G4BiasingProcessInterface*         callingProcess,
                                G4BiasingAppliedCase                        biasingCase,
                                G4VBiasingOperation*          occurenceOperationApplied,
                                G4double                  weightForOccurenceInteraction,
                                G4VBiasingOperation*         finalStateOperationApplied, 
                                const G4VParticleChange*         particleChangeProduced ) final;
          
        void ExitBiasing( const G4Track*, const G4BiasingProcessInterface* ) final;

        G4ParticleDefinition* fParticleToBias;
        G4BOptrForceCollision* optr;
        G4BOptrForceCollision* fCurrentOperator;
    };
} // namespace allpix

#endif
