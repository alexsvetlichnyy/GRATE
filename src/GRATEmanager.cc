

#include "GRATEmanager.hh"

 #include "TFile.h"
 #include "TH1D.h"
 #include "TH2D.h" 
 #include "TMath.h"


GRATEmanager::GRATEmanager()
  : sourceZ(0), sourceA(0), KinEn(-1.), XsectNN(-1.), lowLimitExEn( 0.), upperLimitExEn( 100.), binsExEn(1), eventsPerBin(1), StatisticsLabel(-1), iterations(-1), wM(0), wP(0), NucleusInputLabel(0)
{  
 std::cout << "######### Abrasion-Ablation model using Glauber Monte Carlo and Geant4" <<std::endl;
 
while(!NucleusInputLabel){
  NucleusInputLabel = 1;  
  std::cout << "Please enter colliding nucleus name (side A). U, Pb, Pbpnrw(with neutron skin) Au, Xe, Ag, Br, Cu, Al, O, Oh(Oxygen parametrized with HO), C is available : ";
  std::cin >> SysA;
  if(SysA =="Pb"){SysA+="*"; sourceA = 208; sourceZ = 82;}
  else if(SysA == "Pbpnrw"){sourceA = 208; sourceZ = 82;}
  else if(SysA == "Cu"){SysA+="2";sourceA = 64; sourceZ = 29;}
  else if(SysA == "O"){sourceA = 16; sourceZ = 8;}
  else if(SysA == "Oh"){sourceA = 16; sourceZ = 8;}
  else if(SysA == "Au"){sourceA = 197; sourceZ = 79;}
  else if(SysA == "Ag"){sourceA = 109; sourceZ = 47;}
  else if(SysA == "Br"){sourceA = 79; sourceZ = 35;}
  else if(SysA == "Xe"){sourceA = 129; sourceZ = 54;}
  else if(SysA == "C"){sourceA = 12; sourceZ = 6;}
  else if(SysA == "Al"){sourceA = 27; sourceZ = 13;}
  else if(SysA == "U"){sourceA = 238; sourceZ = 92;}
  else if(SysA == "U2"){sourceA = 238; sourceZ = 92;}
  else{ G4Exception("Nucleus input in GRATE", "GRATE-0", JustWarning, "There is no matched nucleus in GRATE");
	NucleusInputLabel = 0;
      }
}
NucleusInputLabel = 0;

while(!NucleusInputLabel){
  NucleusInputLabel = 1;  
   
 std::cout << "Please enter colliding nucleus name (side B). U, Pb, Pbpnrw(with neutron skin) Au, Xe, Ag, Br, Cu, Al, O, Oh(Oxygen parametrized with HO), C is available : ";
  std::cin >> SysB;
 
  if(SysB == "Pb"){SysB+="*"; sourceAb = 208; sourceZb = 82;}
  else if(SysB == "Pbpnrw"){sourceAb = 208; sourceZb = 82;}
  else if(SysB == "Cu"){sourceAb = 64; sourceZb = 29;}
  else if(SysB == "O"){sourceAb = 16; sourceZb = 8;}
  else if(SysB == "Oh"){sourceA = 16; sourceZ = 8;}
  else if(SysB == "Au"){sourceAb = 197; sourceZb = 79;}
  else if(SysB == "Ag"){sourceAb = 109; sourceZb = 47;}
  else if(SysB == "Br"){sourceAb = 79; sourceZb = 35;}
  else if(SysB == "Xe"){sourceAb = 129; sourceZb = 54;}
  else if(SysB == "C"){sourceAb = 12; sourceZb = 6;}
  else if(SysB == "Al"){sourceAb = 27; sourceZb = 13;}
  else if(SysB == "U"){sourceAb = 238; sourceZb = 92;}
  else if(SysB == "U2"){sourceAb = 238; sourceZb = 92;}
  else{ G4Exception("Nucleus input in GRATE", "GRATE-0", JustWarning, "There is no matched nucleus in GRATE");
	NucleusInputLabel = 0;
      }
}

 while ( KinEn<0. ) {
  std::cout << "Please enter kinetic enegy of projectile nucleus (per nucleon in GeV) : ";
  std::cin >> KinEn;
  }

  pZ = pow(KinEn*(KinEn+2*0.983),0.5)*GeV;
  KinEn*= sourceA*GeV;
  
  lowLimitExEnB = lowLimitExEn;
  upperLimitExEnB = upperLimitExEn;
  lowLimitExEn *= sourceA;
  lowLimitExEn *= MeV;
  upperLimitExEn *=sourceA;
  upperLimitExEn *=MeV;
  lowLimitExEnB *= sourceAb;
  lowLimitExEnB *= MeV;
  upperLimitExEnB *=sourceAb;
  upperLimitExEnB *=MeV;

 while ( (StatisticsLabel<0) || (StatisticsLabel>3) || (upperLimitExEn<lowLimitExEn) ) {
  std::cout << "Please choose the level density function to be used: 1 - Ericson, 2 - Gaimard-Schmidt, 3 - ALADIN parametrization : ";
  std::cin >> StatisticsLabel;
  }

 std::cout<<"Write momentum of a fragments? ";
 std::cin>>wM;


 std::cout<<"Write pseudorapidity of a fragments? ";
 std::cin>>wP;

 while ( (iterations<0) || (iterations>10000000) ) {
    std::cout<<"Please enter number of events to be generated: ";
    std::cin >> iterations;

  } 

  std::cout << "Please enter the file name to write histograms (.root will be supplied): ";
  std::cin >> fileName; 
 


  for (G4int j=0; j<20; j++) histo[j] = 0;
  for (G4int l=0; l<10; l++) histo2[l] = 0;

  

}


GRATEmanager::~GRATEmanager()
{
  
}

void GRATEmanager::BookHisto()
{

// Open a file to keep histograms inside it 

 if ( fileName == "") fileName = "GRATE_"+SysA+SysB+"_"+std::to_string(KinEn/GeV)+"_GeV_"+std::to_string(iterations)+"_events";
 fileType = "root";
 fileFullName = fileName+"."+fileType;
 compressionFactor = 9;
 fFile = new TFile(fileFullName, "RECREATE", fileName, compressionFactor);
 G4cout << "Histograms will be written to " << fileFullName << G4endl;
//Creating a Trees
 Glauber = new TTree("Glauber","Events from glauber modeling");
 modelingCo = new TTree("Conditions","preconditions for modeling");
// Book all histograms there ...
 histo[0] =  new TH1D("Charge distruibution for side B"," ;Z;entries",sourceZb+1,-0.5, sourceZb+0.5); 

 histo[1] =  new TH1D("M distr"," ;M;entries",100, -0.5, 100+0.5);

 histo[2] =  new TH1D("pz for neutrons",";pz;",1000,0, pZ+5000);
 histo[3] =  new TH1D("pz for protons"," ;pz;",1000, 0, pZ+5000);
 histo[4] =  new TH1D("pz for IMF"," ;pz;",1000, -6*1e4, 6*1e4);
 histo[5] =  new TH1D("pz for heavy fragments"," ;pz;",1000, -6*1e4, 6*1e4);


 histo[6] =  new TH1D("Charge distruibution"," ;Z;entries",sourceZ+1,-0.5, sourceZ+0.5);

 histo[7] =  new TH1D("Mass distribution", " ;A,entries",sourceA+1,-0.5, sourceA+0.5);
 
 histo2[1] = new TH2D("Ex En distribution"," ;E*/A;A_{pf}/A",300, 0, 13, sourceA+1, 0, 1);
 
 histo2[2] = new TH2D("Mass and Charge distribution"," ;Z;A",sourceZ+1, -0.5, sourceZ+0.5, sourceA+1, -0.5, sourceA+0.5);

 histo2[3] = new TH2D("px vs py for neutrons", ";px;py", 100,-200,200,100,-200,200);
 histo2[4] = new TH2D("px vs py for protons", ";px;py", 100,-200,200,100,-200,200);
 histo2[5] = new TH2D("px vs py for IMF", ";px;py", 100,-200,200,100,-200,200);
 histo2[6] = new TH2D("px vs py for heavy fragments", ";px;py", 100,-200,200,100,-200,200);
 
}


void GRATEmanager::CalcXsectNN()
{
if(KinEn/G4double(sourceA) < 425*GeV){
	G4double Tkin[2];
	G4double xsect[2];
	XsectFile.open("../src/bystricky.dat");

	while (Tkin[0]*GeV < KinEn/G4double(sourceA)) {
		Tkin[0] = Tkin[1];	
		xsect[0] = xsect[1];
		XsectFile >> Tkin[1] >> xsect[1];
		if (!XsectFile.good()) break;
	}
	G4double a = (xsect[1]-xsect[0])/(Tkin[1]-Tkin[0]);
	G4double b = xsect[1] - a*Tkin[1];
	XsectNN = a*KinEn/(G4double(sourceA)*GeV)+b;
   }
else{
   //G4double s = KinEn*KinEn;
   G4double S = 4*0.938272*GeV*0.938272*GeV+2*(KinEn/G4double(sourceA))*0.938272*GeV;
   XsectNN = 25.0+0.146*pow(log(S/(GeV*GeV)),2);
   }

}


void GRATEmanager::CleanHisto()
{ 
  
   

  fFile->Write();
  G4cout << "\n----> Histograms were written into the file " << fileFullName << G4endl;
  //  delete [] histo;
  //  delete [] histo2;
  delete fFile;
}


void GRATEmanager::FillConditionsTree(G4double Xsect){

G4double XsectTot = 0;
G4double KineticEnergy = 0;
G4double pz = 0;
G4double Mass_on_A = 0;
G4double Mass_on_B = 0;
G4double Charge_on_A = 0;
G4double Charge_on_B = 0;

modelingCo->Branch("Xsect_total", &XsectTot,"Xsect_total/d");
modelingCo->Branch("Kinetic_energy_per_nucleon_of_projectile_in_GeV", &KineticEnergy,"Kinetic_energy_of_per_nucleon_projectile_in_GeV/d");
modelingCo->Branch("pZ_in_MeV", &pz,"pZ_in_MeV/d");
modelingCo->Branch("Mass_on_A", &Mass_on_A,"Mass_on_A/d");
modelingCo->Branch("Mass_on_B", &Mass_on_B,"Mass_on_B/d");
modelingCo->Branch("Charge_on_A", &Charge_on_A,"Charge_on_A/d");
modelingCo->Branch("Charge_on_B", &Charge_on_B,"Charge_on_B/d");

XsectTot = Xsect;
KineticEnergy = KinEn/(sourceA*GeV);
pz = pZ/MeV;
Mass_on_A = sourceA;
Mass_on_B = sourceAb;
Charge_on_A = sourceZ;
Charge_on_B = sourceZb;

modelingCo->Fill();

}

