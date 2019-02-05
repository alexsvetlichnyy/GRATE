

#include "HandlerTestHistoManager.hh"

 #include "TFile.h"
 #include "TH1D.h"
 #include "TH2D.h" 


HandlerTestHistoManager::HandlerTestHistoManager()
  : Z(0), A(0),sourceZ(0),sourceA(0), lowLimitExEn(-1.), upperLimitExEn(-1.), binsExEn(1), eventsPerBin(1), StatisticsLabel(-1)
{  
  std::cout << "######### Abrasion-Ablation model using Geant4" << G4endl;
  
  while ( (sourceZ<1) || (sourceZ>100) ){
  std::cout << "Please enter Z of colliding nuclei from 1 to 100: ";
  std::cin >> sourceZ;
  }
  
  while ( (sourceA<1) || (sourceA>300) ){
  std::cout << "Please enter A of colliding nuclei from 1 to 300: ";
  std::cin >> sourceA;
  }

  
  while ( (lowLimitExEn<0.) || (lowLimitExEn>100.) ) {
  std::cout << "Please choose the low limit for the range of excitation energy (per nucleon in MeV) : from 0 to 100: ";
  std::cin >> lowLimitExEn;
  }

  while ( (upperLimitExEn<0.) || (upperLimitExEn>100.) || (upperLimitExEn<lowLimitExEn) ) {
  std::cout << "Please choose the upper limit for the range of excitation energy (per nucleon in MeV) : from 0 to 100: ";
  std::cin >> upperLimitExEn;
  }
  lowLimitExEn *= sourceA;
  lowLimitExEn *= MeV;
  upperLimitExEn *=sourceA;
  upperLimitExEn *=MeV;

 while ( (StatisticsLabel<0) || (StatisticsLabel>3) || (upperLimitExEn<lowLimitExEn) ) {
  std::cout << "Please choose the level density function to be used: 1 - Ericson, 2 - Gaimard-Schmidt, 3 - ALADIN parametrization : ";
  std::cin >> StatisticsLabel;
  }


  while ( (iterations<10) || (iterations>10000000) ) {
    std::cout<<"Please enter number of iterations: ";
    std::cin >> iterations;

  } 

  std::cout << "Please enter the file name to write histograms (.root will be supplied): ";
  std::cin >> fileName;
  std::cout << "Please enter the file name to open with the path: ";
  std::cin>> IsotopMapPath;  
 


  for (G4int j=0; j<20; j++) histo[j] = 0;
  for (G4int l=0; l<10; l++) histo2[l] = 0;

}


HandlerTestHistoManager::~HandlerTestHistoManager()
{
  
}

void HandlerTestHistoManager::OpenHisto()
{
 TFile *IsotopMap= new TFile(IsotopMapPath);
 IsotopMap->ls();
 IsotopMapHisto = (TH2D*)IsotopMap->Get("h2"); 
 std::cout<<"Histogram graped"<<std::endl;
}

void HandlerTestHistoManager::BookHisto()
{

// Open a file to keep histograms inside it ... 

 if ( fileName == "") fileName = "Dexcitation";
 fileType = "root";
 fileFullName = fileName+"."+fileType;
 compressionFactor = 8;
 fFile = new TFile(fileFullName, "RECREATE", fileName, compressionFactor);
 G4cout << "Histograms will be written to " << fileFullName << G4endl;

// Book all histograms there ...

 histo[0] =  new TH1D("multip","Average multiplicity",binsExEn,lowLimitExEn/A,upperLimitExEn/A);

 histo[2] =  new TH1D("M distr"," ;M;entries",100, -0.5, 100+0.5);

 histo[6] =  new TH1D("Charge distruibution"," ;Z;entries",sourceZ+1,-0.5, sourceZ+0.5);

 histo[8] =  new TH1D("Charge distruibution for M>3"," ;Z;entries",sourceZ+1,-0.5, sourceZ+0.5); 

 histo[7] =  new TH1D("Mass distribution", " ;A,entries",sourceA+1,-0.5, sourceA+0.5);
 
 histo2[1] = new TH2D("Ex En distribution"," ;E*/A;A_pf/A",1000, 0, 13, sourceA+1, 0, 1);
 
 histo2[2] = new TH2D("Mass and Charge distribution"," ;Z;A",sourceZ+1, -0.5, sourceZ+0.5, sourceA+1, -0.5, sourceA+0.5);

}

void HandlerTestHistoManager::CloseMap(){
  IsotopMap->Close();
}


void HandlerTestHistoManager::CleanHisto()
{
  fFile->Write();
  G4cout << "\n----> Histograms were written into the file " << fileFullName << G4endl;
  //  delete [] histo;
  //  delete [] histo2;
  delete fFile;
}

