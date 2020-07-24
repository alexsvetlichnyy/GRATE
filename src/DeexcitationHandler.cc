#include "DeexcitationHandler.hh"

DeexcitationHandler::DeexcitationHandler(){
}

DeexcitationHandler::~DeexcitationHandler() = default;

G4ReactionProductVector *DeexcitationHandler::BreakUp(const G4Fragment &theInitialFragment) {
    if(theInitialFragment.GetA()>MaxAforFermiBreakUpForPureNeutronFragments && theInitialFragment.GetZ()==0){return BreakUpPureNeutrons(theInitialFragment);}
    else { return BreakItUp(theInitialFragment);} //parent class method
} 

G4ReactionProductVector *DeexcitationHandler::BreakUpPureNeutrons(const G4Fragment &theInitialFragment) {
   G4ParticleDefinition* neutron = G4Neutron::NeutronDefinition();
    G4ReactionProductVector* outVec = new G4ReactionProductVector();
    std::vector<G4double> mr;
    mr.reserve(theInitialFragment.GetA_asInt());
    G4double M = theInitialFragment.GetMomentum().m();
    for(G4int k = 0; k<theInitialFragment.GetA_asInt(); k++){
        mr.push_back(mn);
    }
    std::vector<G4LorentzVector*>* mom = PhaseSpaceDecay.Decay(M,mr);
    for(G4int k = 0; k<theInitialFragment.GetA_asInt(); k++){
        G4ReactionProduct* New = new G4ReactionProduct(neutron);
        New->SetMomentum(mom->at(k)->vect());
        New->SetTotalEnergy(mn);
        outVec->push_back(New);
    }

    return outVec;
};

