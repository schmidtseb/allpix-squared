#ifndef PHYSICSLIST_HPP
#define PHYSICSLIST_HPP

#include "G4VModularPhysicsList.hh"

class PhysicsList: public G4VModularPhysicsList {
	public:
		PhysicsList(bool);
		~PhysicsList();

	private:
		void ConstructParticle();
		void ConstructProcess();
		void ConstructEM();
		void ConstructDecay();
		void ConstructDeexcite();
		void SetCuts();

		bool fEnableDecay;
};

#endif
