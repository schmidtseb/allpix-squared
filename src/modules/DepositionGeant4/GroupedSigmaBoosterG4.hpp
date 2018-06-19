#ifndef GROUPED_SIGMA_BOOSTER_H
#define GROUPED_SIGMA_BOOSTER_H 1

#include "globals.hh"
#include "G4VProcess.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4VDiscreteProcess.hh"
#include "G4Track.hh"

class GroupedSigmaBooster : public G4VContinuousDiscreteProcess {
    public:
        GroupedSigmaBooster(const std::string processName = "SigmaBooster",
            G4VProcess *wrapProc = nullptr) : 
        G4VContinuousDiscreteProcess(processName + ":" + wrapProc->GetProcessName()),
        lastStepBoost(1.0) {
            proc = dynamic_cast<G4VDiscreteProcess*>(wrapProc);
            // if (!proc)
                // G4Exception("Cannont wrap a non-G4VDiscreteProcess with sigma booster: "+wrapProc->GetProcessName());
            pParticleChange = new G4ParticleChange;
        }

        ~GroupedSigmaBooster() {
            if(pParticleChange) delete pParticleChange;
            pParticleChange = nullptr;
        }

        virtual double PostStepGetPhysicalInteractionLength(const G4Track& aTrack, double prev, G4ForceCondition* cond);
        virtual double AlongStepGetPhysicalInteractionLength(const G4Track&, double, double, double&, G4GPILSelection*) {
            return DBL_MAX;
        }
        virtual double GetContinuousStepLimit(const G4Track&, double, double, double&) {
            return DBL_MAX;
        }
        virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
        virtual G4VParticleChange* AlongStepDoIt(const G4Track& aTrack, const G4Step& aStep);

        virtual void ResetNumberOfInteractionLengthLeft() {
            lastStepBoost = 1.0;
            G4VContinuousDiscreteProcess::ResetNumberOfInteractionLengthLeft();
            proc->ResetNumberOfInteractionLengthLeft();
        }

        virtual void StartTracking(G4Track* aTrack) {
            G4VContinuousDiscreteProcess::StartTracking(aTrack);
            proc->StartTracking(aTrack);
        }

        virtual void EndTracking() {
            G4VContinuousDiscreteProcess::EndTracking();
            proc->EndTracking();
        }

        virtual void SetProcessManager(const G4ProcessManager* procMan) {
            proc->SetProcessManager(procMan);
        }
        virtual const G4ProcessManager* GetProcessManager() {
            return proc->GetProcessManager();
        }
        virtual bool IsApplicable(const G4ParticleDefinition& aParticleType) {
            return proc->IsApplicable(aParticleType);
        }
        virtual void BuildPhysicsTable(const G4ParticleDefinition& aParticleType) {
            proc->BuildPhysicsTable(aParticleType);
        }
        virtual void PreparePhysicsTable(const G4ParticleDefinition& aParticleType) {
            proc->PreparePhysicsTable(aParticleType);
        }
        virtual bool StorePhysicsTable(const G4ParticleDefinition* ptcl, const G4String& fname, bool ascii = false) {
            return proc->StorePhysicsTable(ptcl, fname, ascii);
        } 
        using G4VContinuousDiscreteProcess::RetrievePhysicsTable;
        virtual bool RetrievePhysicsTable(const G4ParticleDefinition* ptcl, const std::string fname, bool ascii=false) {
            return proc->RetrievePhysicsTable(ptcl, fname, ascii);
        }

        static void SetCrossSectionBias(G4double BiasFactor=1) {
            biasFactor = BiasFactor;
        }
        static void SetPrimaryOnly(bool flag) {
            primary_only = flag;
        }
        static void SetUseWeighting(bool flag) {
            useWeighting = flag;
        }
        static void SetEnableProbabilityConservation(bool flag) {
            conserve_probability = flag;
        }

        G4VDiscreteProcess *GetProcess() {
            return proc;
        }

        static double GetAccumulatedWeight() {
            return accumulatedWeight;
        }

        static void ResetAccumulatedWeight() {
            accumulatedWeight = 1;
            weightInitialized = true;
        }

    protected:
        using G4VContinuousDiscreteProcess::GetMeanFreePath;
        virtual G4double GetMeanFreePath(const G4Track& aTrack, G4double prev, G4ForceCondition* cond) {
            return proc->PostStepGetPhysicalInteractionLength(aTrack, prev*lastStepBoost, cond);
        }
        G4VDiscreteProcess *proc;
        double lastStepBoost, lastInteractionLengths;

        static double biasFactor;
        static bool primary_only, useWeighting, conserve_probability;

        static double accumulatedWeight;
        static bool weightInitialized;

    private:
};

#endif
