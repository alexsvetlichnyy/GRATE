#include "InitialConditions.hh"

InitialConditions::InitialConditions() {

}

InitialConditions::~InitialConditions() {

}
InitialConditions::InitialConditions(G4double KinEn_in, G4String SysA_in, G4String SysB_in, G4bool IsCollider_in) {

    if(SysA_in =="Pb"){SysA+="*"; sourceA = 208; sourceZ = 82; SysA = SysA_in;}
    else if(SysA_in == "Pbpnrw"){sourceA = 208; sourceZ = 82; SysA = SysA_in;}
    else if(SysA_in == "Cu"){SysA+="2";sourceA = 64; sourceZ = 29; SysA = SysA_in;}
    else if(SysA_in == "O"){sourceA = 16; sourceZ = 8; SysA = SysA_in;}
    else if(SysA_in == "Au"){sourceA = 197; sourceZ = 79; SysA = SysA_in;}
    else if(SysA_in == "Ag"){sourceA = 109; sourceZ = 47; SysA = SysA_in;}
    else if(SysA_in == "Br"){sourceA = 79; sourceZ = 35; SysA = SysA_in;}
    else if(SysA_in == "Xe"){sourceA = 129; sourceZ = 54; SysA = SysA_in;}
    else if(SysA_in == "C"){sourceA = 12; sourceZ = 6; SysA = SysA_in;}
    else if(SysA_in == "Al"){sourceA = 27; sourceZ = 13; SysA = SysA_in;}
    else if(SysA_in == "U"){sourceA = 238; sourceZ = 92; SysA = SysA_in;}
    else if(SysA_in == "U2"){sourceA = 238; sourceZ = 92; SysA = SysA_in;}
    else{ G4Exception("Nucleus input in GRATE", "GRATE-0", FatalErrorInArgument, "There is no matched nucleus in GRATE");
    }

    if(SysB_in == "Pb"){SysB+="*"; sourceAb = 208; sourceZb = 82; SysB = SysB_in;}
    else if(SysB_in == "Pbpnrw"){sourceAb = 208; sourceZb = 82; SysB = SysB_in;}
    else if(SysB_in == "Cu"){SysB+="2";sourceAb = 64; sourceZb = 29; SysB = SysB_in;}
    else if(SysB_in == "O") {sourceAb = 16; sourceZb = 8; SysB = SysB_in;}
    else if(SysB_in == "Au"){sourceAb = 197; sourceZb = 79; SysB = SysB_in;}
    else if(SysB_in == "Ag"){sourceAb = 109; sourceZb = 47; SysB = SysB_in;}
    else if(SysB_in == "Br"){sourceAb = 79; sourceZb = 35; SysB = SysB_in;}
    else if(SysB_in == "Xe"){sourceAb = 129; sourceZb = 54; SysB = SysB_in;}
    else if(SysB_in == "C") {sourceAb = 12; sourceZb = 6; SysB = SysB_in;}
    else if(SysB_in == "Al"){sourceAb = 27; sourceZb = 13; SysB = SysB_in;}
    else if(SysB_in == "U") {sourceAb = 238; sourceZb = 92; SysB = SysB_in;}
    else if(SysB_in == "U2"){sourceAb = 238; sourceZb = 92; SysB = SysB_in;}
    else{ G4Exception("Nucleus input in GRATE", "GRATE-0", FatalErrorInArgument, "There is no matched nucleus in GRATE");
    }

    IsCollider = IsCollider_in;
    if(IsCollider){
        PzA =    pow(KinEn_in*KinEn_in*0.25 - nucleonAverMass*nucleonAverMass,0.5);
        PzB = -1*pow(KinEn_in*KinEn_in*0.25 - nucleonAverMass*nucleonAverMass,0.5);
    } else{
        PzA = pow(KinEn_in*(KinEn_in+2*nucleonAverMass),0.5);
        PzB = 0;
    }
    PzA*= sourceA*GeV;
    PzB*= sourceAb*GeV;
    KinEn = KinEn_in*sourceA*GeV;
}

G4bool InitialConditions::SetSysA(G4String SysA_in) {
    if(SysA_in == "Pb"){ sourceA = 208; sourceZ = 82; SysA = SysA_in; SysA+="*";}
    else if(SysA_in == "Pbpnrw"){sourceA = 208; sourceZ = 82; SysA = SysA_in;}
    else if(SysA_in == "Cu"){sourceA = 64; sourceZ = 29; SysA = SysA_in; SysA+="2";}
    else if(SysA_in == "O"){sourceA = 16; sourceZ = 8; SysA = SysA_in;}
    else if(SysA_in == "Au"){sourceA = 197; sourceZ = 79; SysA = SysA_in;}
    else if(SysA_in == "Ag"){sourceA = 109; sourceZ = 47; SysA = SysA_in;}
    else if(SysA_in == "Br"){sourceA = 79; sourceZ = 35; SysA = SysA_in;}
    else if(SysA_in == "Xe"){sourceA = 129; sourceZ = 54; SysA = SysA_in;}
    else if(SysA_in == "C"){sourceA = 12; sourceZ = 6; SysA = SysA_in;}
    else if(SysA_in == "Al"){sourceA = 27; sourceZ = 13; SysA = SysA_in;}
    else if(SysA_in == "U"){sourceA = 238; sourceZ = 92; SysA = SysA_in;}
    else if(SysA_in == "U2"){sourceA = 238; sourceZ = 92; SysA = SysA_in;}
    else{ G4Exception("Nucleus input in GRATE", "GRATE-0", JustWarning, "There is no matched nucleus in GRATE");
        return 0;
    }

    return 1;
}

G4bool InitialConditions::SetSysB(G4String SysB_in) {
    if(SysB_in == "Pb"){sourceAb = 208; sourceZb = 82; SysB = SysB_in; SysB+="*";}
    else if(SysB_in == "Pbpnrw"){sourceAb = 208; sourceZb = 82; SysB = SysB_in;}
    else if(SysB_in == "Cu"){sourceAb = 64; sourceZb = 29; SysB = SysB_in; SysB+="2";}
    else if(SysB_in == "O") {sourceAb = 16; sourceZb = 8; SysB = SysB_in;}
    else if(SysB_in == "Au"){sourceAb = 197; sourceZb = 79; SysB = SysB_in;}
    else if(SysB_in == "Ag"){sourceAb = 109; sourceZb = 47; SysB = SysB_in;}
    else if(SysB_in == "Br"){sourceAb = 79; sourceZb = 35; SysB = SysB_in;}
    else if(SysB_in == "Xe"){sourceAb = 129; sourceZb = 54; SysB = SysB_in;}
    else if(SysB_in == "C") {sourceAb = 12; sourceZb = 6; SysB = SysB_in;}
    else if(SysB_in == "Al"){sourceAb = 27; sourceZb = 13; SysB = SysB_in;}
    else if(SysB_in == "U") {sourceAb = 238; sourceZb = 92; SysB = SysB_in;}
    else if(SysB_in == "U2"){sourceAb = 238; sourceZb = 92; SysB = SysB_in;}
    else{ G4Exception("Nucleus input in GRATE", "GRATE-0", JustWarning, "There is no matched nucleus in GRATE");
        return 0;
    }

    return 1;
}

void InitialConditions::SetKinEn(G4double KinEn_in) {

    if(IsCollider){
        PzA =    pow(KinEn_in*KinEn_in*0.25 - nucleonAverMass*nucleonAverMass,0.5);
        PzB = -1*pow(KinEn_in*KinEn_in*0.25 - nucleonAverMass*nucleonAverMass,0.5);
    } else{
        PzA = pow(KinEn_in*(KinEn_in+2*nucleonAverMass),0.5);
        PzB = 0;
    }
    PzA*= sourceA*GeV;
    PzB*= sourceAb*GeV;
    KinEn = KinEn_in*sourceA*GeV;
}


