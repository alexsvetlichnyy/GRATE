
// Igor Pshenichnov 06.08.2008

#ifndef GRATEmanager_h
#define GRATEmanager_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "TTree.h"
#include "TParticle.h"
#include <fstream>
#include "math.h"
class TFile;
class TH1D;
class TH2D;

class GRATEmanager
{
  public:

  GRATEmanager();
   ~GRATEmanager();

  public:

  TH1D* GetHisto(G4int id) {return histo[id];};
  TTree* GetTree() {return Glauber;};
  TH2D* GetHisto2(G4int id) {return histo2[id];};
 
       
  void BookHisto();
  void CalcXsectNN();
  void CleanHisto();
  void CalcHyppGeomHisto();
  void CalcHyppGeomArray(G4int RestNucleons);
  void CalcHyppGeomArrayB(G4int RestNucleons);
  void FillConditionsTree(G4double Xsect);
  
  inline G4String GetSysA() {return SysA;}
  inline G4String GetSysB() {return SysB;}
  inline G4int GetSourceZ() {return sourceZ;}
  inline G4int GetSourceA() {return sourceA;}
  inline G4int GetSourceZb() {return sourceZb;}
  inline G4int GetSourceAb() {return sourceAb;}
  inline G4int GetStatType() {return StatisticsLabel;}
  inline G4int GetIterations()  {return iterations;};
  inline G4double GetXsectNN() {return XsectNN;}
  inline G4double GetKinEn() {return KinEn;};
  inline G4double GetLowEn() {return lowLimitExEn;};
  inline G4double GetUpEn() {return upperLimitExEn;};
  inline G4double GetLowEnB() {return lowLimitExEnB;};
  inline G4double GetUpEnB() {return upperLimitExEnB;};
  inline G4double GetHyppGeomArray(G4int id) {return HyppGeomArray[id];};
  inline G4double GetHyppGeomArrayB(G4int id) {return HyppGeomArrayB[id];};


  
  private:

  
    TFile* fFile;
    TFile* IsotopMap;
    TH1D*  histo[20];
    TH2D*  HyppGeomHisto2;
    TH2D*  HyppGeomHisto2b;
    TH2D*  histo2[10];
    TTree* Glauber;
    TTree* modelingCo;
 
        
    G4double HyppGeomArray[100];
    G4double HyppGeomArrayB[100];
    G4int sourceZ;
    G4int sourceA;
    G4int sourceZb;
    G4int sourceAb;
    G4int iterations; 
    G4int StatisticsLabel;   

    G4String fileName;
    G4String fileType;
    G4String fileOpenPath;

    G4String     fileFullName;
    G4String     SysA;
    G4String     SysB;
    G4int        compressionFactor;

    G4double XsectNN;

    G4double KinEn;	
    G4double lowLimitExEn;
    G4double upperLimitExEn;
    G4double lowLimitExEnB;
    G4double upperLimitExEnB;
    G4int    binsExEn;
    G4int    eventsPerBin;      

    std::ifstream XsectFile;

};

#endif
