#ifndef ALLPIX_CROSS_SECTION_BIASING_OPERATOR_G4_H
#define ALLPIX_CROSS_SECTION_BIASING_OPERATOR_G4_H

#include "G4VBiasingOperator.hh"
#include <map>

#include "core/config/Configuration.hpp"

class G4BOptnChangeCrossSection;
class G4ParticleDefinition;

namespace allpix {
    class CrossSectionBiasingOperatorG4 : public G4VBiasingOperator {
    public:
        explicit CrossSectionBiasingOperatorG4(G4ParticleDefinition* particleToBias, std::string name);
        ~CrossSectionBiasingOperatorG4();

        // Method called at the beginning of a run
        virtual void StartRun();

    private:
        virtual G4VBiasingOperation*
        ProposeOccurenceBiasingOperation(const G4Track*                            track,
            const G4BiasingProcessInterface* callingProcess);

        // -- Methods not used:
        virtual G4VBiasingOperation*
        ProposeFinalStateBiasingOperation(const G4Track*, const G4BiasingProcessInterface*)
        {return nullptr;}
        virtual G4VBiasingOperation*
        ProposeNonPhysicsBiasingOperation(const G4Track*, const G4BiasingProcessInterface*)
        {return nullptr;}

        using G4VBiasingOperator::OperationApplied;
        virtual void OperationApplied( const G4BiasingProcessInterface*                callingProcess,
            G4BiasingAppliedCase biasingCase,
            G4VBiasingOperation* occurenceOperationApplied,
            G4double weightForOccurenceInteraction,
            G4VBiasingOperation* finalStateOperationApplied, 
            const G4VParticleChange* particleChangeProduced );

        // List of associations between processes and biasing operations:
        std::map<const G4BiasingProcessInterface*, G4BOptnChangeCrossSection*> fChangeCrossSectionOperations;
        const G4ParticleDefinition* fParticleToBias;
        bool fSetup;
    };
} // namespace allpix

#endif
