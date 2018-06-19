#include "GroupedSigmaBoosterG4.hpp"

double GroupedSigmaBooster::biasFactor = 1.0;
bool GroupedSigmaBooster::primary_only = 0;
bool GroupedSigmaBooster::useWeighting = 0;
bool GroupedSigmaBooster::conserve_probability = true;
double GroupedSigmaBooster::accumulatedWeight = 0;
bool GroupedSigmaBooster::weightInitialized = false;

double GroupedSigmaBooster::PostStepGetPhysicalInteractionLength(const G4Track& aTrack, G4double prev, G4ForceCondition* cond) {
	if(biasFactor == 0.0) {
		lastStepBoost = 0.0;
		return DBL_MAX;
	}
	double pil = proc->PostStepGetPhysicalInteractionLength(aTrack, prev*lastStepBoost, cond);

	if(!primary_only || aTrack.GetTrackID() == 1) {
		pil /= biasFactor;
		lastStepBoost = biasFactor;
	} else lastStepBoost = 1.0;

	return pil;
}

G4VParticleChange* GroupedSigmaBooster::AlongStepDoIt(const G4Track& aTrack, const G4Step& aStep) {
	pParticleChange->Initialize(aTrack);

	if(conserve_probability && lastStepBoost != 1.0) {
		G4ParticleChange* aParticleChange_ = static_cast<G4ParticleChange*>(pParticleChange);
		aParticleChange_->SetParentWeightByProcess(false);
		double adjust = std::exp((lastStepBoost - 1)*(aStep.GetStepLength()/proc->GetCurrentInteractionLength()));
		double newWeight = aTrack.GetWeight()*adjust;
		aParticleChange_->ProposeWeight(newWeight);
		if(weightInitialized) accumulatedWeight *= adjust;
	}

	return pParticleChange;
}

G4VParticleChange* GroupedSigmaBooster::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) {
	G4VParticleChange* vParticleChange = proc->PostStepDoIt(aTrack, aStep);
	G4ParticleChange* aParticleChange_;
	if(useWeighting && lastStepBoost != 1.0 && (aParticleChange_ = dynamic_cast<G4ParticleChange*>(vParticleChange)) != nullptr) {
		int nSec = aParticleChange_->GetNumberOfSecondaries();
		double newWeight = aTrack.GetWeight() / lastStepBoost;
		for(int i = 0; i < nSec; i++) {
			aParticleChange_->GetSecondary(i)->SetWeight(newWeight);
		}
		aParticleChange_->SetParentWeightByProcess(false);
		aParticleChange_->ProposeWeight(newWeight);

		if(weightInitialized) accumulatedWeight /= lastStepBoost;
	}

	return vParticleChange;
}
