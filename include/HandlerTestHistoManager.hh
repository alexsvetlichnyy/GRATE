
// Igor Pshenichnov 06.08.2008

#ifndef HandlerTestHistoManager_h
#define HandlerTestHistoManager_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
class TFile;
class TH1D;
class TH2D;

class HandlerTestHistoManager
{
  public:

  HandlerTestHistoManager();
   ~HandlerTestHistoManager();

  public:

  TH1D* GetHisto(G4int id) {return histo[id];};
  TH2D* GetHisto2(G4int id) {return histo2[id];};
  TH2D* GetIsotopMap() {return IsotopMapHisto;};
       
  void BookHisto();
  void NormalizeHisto();
  void CleanHisto();
  void OpenHisto();
  void CloseMap();
  
  inline G4int GetSourceZ() {return sourceZ;}
  inline G4int GetSourceA() {return sourceA;}
  inline G4int GetStatType() {return StatisticsLabel;}
  inline G4int GetZ()    {return Z;};
  inline G4int GetA()    {return A;};
  inline G4int GetIterations()  {return iterations;};
  inline G4double GetLowEn() {return lowLimitExEn;};
  inline G4double GetUpEn() {return upperLimitExEn;};


  
  private:

  
    TFile* fFile;
    TFile* IsotopMap;
    TH1D*  histo[20];
    TH2D*  histo2[10];
    TH2D* IsotopMapHisto;
    
    G4int sourceZ;
    G4int sourceA;
    G4int Z;
    G4int A;
    G4int iterations; 
    G4int StatisticsLabel;   

    G4String fileName;
    G4String fileType;
    G4String fileOpenPath;

    G4String     fileFullName;
    G4String     IsotopMapPath;
    G4int        compressionFactor;

    G4double lowLimitExEn;
    G4double upperLimitExEn;
    G4int    binsExEn;
    G4int    eventsPerBin;      

};

#endif
