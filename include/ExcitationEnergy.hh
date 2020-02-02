#ifndef ExcitationEnergy_h;
#define ExcitationEnergy_h 1

#include "globals.hh"
#include "Randomize.hh"
#include "G4ExceptionHandler.hh"
#include "cmath"
#include <fstream>

class ExcitationEnergy {
public:
    ExcitationEnergy(G4int ExEnLabel_in, G4int initA_in);

    ~ExcitationEnergy();

public:
    G4double GetEnergy(G4int A);

    G4double GetEnergyEricson(G4int A);

    G4double GetEnergyGaimardSchmidt(G4int A);

    G4double GetEnergyALADIN(G4int A);

    G4double GetEnergyCorrectedALADIN(G4int A);

    void SetInitNuclMass(G4int initA_in);

    void SetParametersEricson(G4double g0_in);

    void SetParametersGaimardSchmidt(G4double g0_in, G4double g1_in);

    void SetParametersALADIN(G4double e0_in, G4double sigma0_in, G4double b0_in);

    void SetParametersCorrectedALADIN(G4double e0_in,G4double c0_in, G4double sigma0_in, G4double b0_in, G4double b1_in);

    void SetParametersCorrectedALADINFromFile();

private:
     G4double g0;
     G4double g1;
     G4double e0;
     G4double sigma0;
     G4double b0;
     G4double c0;
     G4double b1;

     G4int initA;
     G4int ExEnLabel;

     G4double LowExEn;
     G4double UpExEn;
     G4double Ebound;

     std::ifstream ParamFile;
};
#endif