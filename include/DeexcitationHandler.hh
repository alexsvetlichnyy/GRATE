
#include <list>

 #include "G4ExcitationHandler.hh"
 #include "G4SystemOfUnits.hh"
 #include "G4LorentzVector.hh"
 #include "G4NistManager.hh"
 #include "G4ParticleTable.hh"
 #include "G4ParticleTypes.hh"
 #include "G4Ions.hh"
 #include "G4Evaporation.hh"
 #include "G4StatMF.hh"
 #include "G4PhotonEvaporation.hh"
 #include "G4Pow.hh"
 #include "G4FermiPhaseSpaceDecay.hh"

 class DeexcitationHandler: public G4ExcitationHandler {
 public:
     DeexcitationHandler();
    ~DeexcitationHandler();

    G4ReactionProductVector* BreakUp(const G4Fragment &theInitialFragment);
    G4ReactionProductVector* BreakUpPureNeutrons(const G4Fragment &theInitialFragment);
    inline void SetMaxAforPureNeutronFragments(G4int in_A) {MaxAforFermiBreakUpForPureNeutronFragments = in_A;};

 private:
    G4int MaxAforFermiBreakUpForPureNeutronFragments = 18;
    G4double mn = 939.5731*MeV;
    G4FermiPhaseSpaceDecay PhaseSpaceDecay;
};
