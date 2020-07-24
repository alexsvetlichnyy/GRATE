
// Igor Pshenichnov 06.08.2008 + Alexandr Svetlichnyi at present
// ROOT 6 or higher and Geant4 10 or higher installation is reqired

#include "G4RunManager.hh"
#include "G4StateManager.hh"

#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4SystemOfUnits.hh"
	
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4ExcitationHandler.hh"
#include "G4NucleiProperties.hh"
#include "G4Evaporation.hh"
#include "GRATEmanager.hh"
#include "InitialConditions.hh"
#include "ExcitationEnergy.hh"
#include "DeexcitationHandler.hh"
#include "GRATEPhysicsList.hh"
#include "TGlauber/TGlauberMC.hh"
#include "TGlauber/TGlauNucleon.hh"
#include "TGlauber/TGlauNucleus.hh"
#include "TVector3.h"
#include "TObjArray.h"
#include "TObject.h"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4Threading.hh"

#include "G4UImanager.hh"
#include "G4IonTable.hh"
#include "G4GenericIon.hh"
#include "G4Ions.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclearLevelData.hh"


#include <fstream>

 #include "TFile.h"
 #include "TH1D.h"
 #include "TH2D.h" 


int main()
{
   //Seting parameters for Deexcitation
  G4NuclearLevelData* fLevelData = G4NuclearLevelData::GetInstance(); 
  G4DeexPrecoParameters* fParam = fLevelData->GetParameters();
  fParam->SetMinExPerNucleounForMF(3*MeV);

  G4StateManager* fStateManager = G4StateManager::GetStateManager();
   
  G4Random::setTheEngine(new CLHEP::RanluxEngine());

  auto seed1 = 0x7ffce81312f0;
  CLHEP::HepRandom::setTheSeeds(CLHEP::HepRandom::getTheSeeds());
   //CLHEP::HepRandom::setTheSeed(seed1);

  CLHEP::HepRandom::setTheEngine(new CLHEP::RanluxEngine);

  G4RunManager * runManager = new G4RunManager;
  runManager->SetUserInitialization(new GRATEPhysicsList);
  G4BosonConstructor pCBos;
  pCBos.ConstructParticle();
 
  G4LeptonConstructor pCLept;
  pCLept.ConstructParticle();

  G4MesonConstructor pCMes;
  pCMes.ConstructParticle();

  G4BaryonConstructor pCBar;
  pCBar.ConstructParticle();

  G4IonConstructor pCIon;
  pCIon.ConstructParticle();

  G4GenericIon* gion = G4GenericIon::GenericIon();
  gion->SetProcessManager(new G4ProcessManager(gion));

  G4StateManager::GetStateManager()->SetNewState(G4State_Init); // To let create ions
  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  G4IonTable* ions = partTable->GetIonTable();
  partTable->SetReadiness();
  ions->CreateAllIon();
  ions->CreateAllIsomer();

  //Definition of decay processes 
  G4Evaporation * theEvaporation   = new G4Evaporation();	


  
  //Definition of level density functions
  
  G4double GaimardSchmidt(G4double, G4double, G4int ,G4int);

  G4double Ericson(G4double, G4double, G4int ,G4int);

  // The user will be asked for the nuclear name to simulate it's collisions. 

  GRATEmanager histoManager;
 
  //arrays for tree creation
  std::vector<G4float> MassOnSideA;
  std::vector<G4float> MassOnSideB;
  std::vector<G4float> ChargeOnSideA;
  std::vector<G4float> ChargeOnSideB;
  std::vector<G4double> pXonSideA;
  std::vector<G4double> pYonSideA;
  std::vector<G4double> pZonSideA;
  std::vector<G4double> pXonSideB;
  std::vector<G4double> pYonSideB;
  std::vector<G4double> pZonSideB;
  std::vector<G4double> pseudorapidity_A;
  std::vector<G4double> pseudorapidity_B;
  G4float b;
  G4float ExEn;
  G4int id;
  G4int Ncoll;
  G4int Ncollpp;
  G4int Ncollpn;
  G4int Ncollnn;

  // Histograms will be booked now.
  histoManager.BookHisto();
  histoManager.GetTree()->Branch("id", &id, "id/i");
  histoManager.GetTree()->Branch("A_on_A", "std::vector" ,&MassOnSideA);
  histoManager.GetTree()->Branch("A_on_B", "std::vector" ,&MassOnSideB);
  histoManager.GetTree()->Branch("Z_on_A", "std::vector" ,&ChargeOnSideA);
  histoManager.GetTree()->Branch("Z_on_B", "std::vector" ,&ChargeOnSideB);
  histoManager.GetTree()->Branch("Ncoll", &Ncoll, "Ncoll/I");
  histoManager.GetTree()->Branch("Ncollpp", &Ncollpp, "Ncollpp/I");
  histoManager.GetTree()->Branch("Ncollpn", &Ncollpn, "Ncollpn/I");
  histoManager.GetTree()->Branch("Ncollnn", &Ncollnn, "Ncollnn/I");

if(histoManager.WritePseudorapidity()){
  histoManager.GetTree()->Branch("pseudorapidity_on_A", "std::vector", &pseudorapidity_A);
  histoManager.GetTree()->Branch("pseudorapidity_on_B", "std::vector", &pseudorapidity_B);
  }
if(histoManager.WriteMomentum()){
  histoManager.GetTree()->Branch("pX_on_A", "std::vector" ,&pXonSideA,128000,1);
  histoManager.GetTree()->Branch("pY_on_A", "std::vector" ,&pYonSideA,128000,1);
  histoManager.GetTree()->Branch("pZ_on_A", "std::vector" ,&pZonSideA,128000,1);
  histoManager.GetTree()->Branch("pX_on_B", "std::vector" ,&pXonSideB,128000,1);
  histoManager.GetTree()->Branch("pY_on_B", "std::vector" ,&pYonSideB,128000,1);
  histoManager.GetTree()->Branch("pZ_on_B", "std::vector" ,&pZonSideB,128000,1);
  }

  histoManager.GetTree()->Branch("impact_parameter", &b, "impact_parameter/f");
  histoManager.GetTree()->Branch("Ex_En_per_nucleon", &ExEn, "Ex_En_per_nucleon/f");
  

  G4bool Verbose = 0;
 
  //Get Z and A of nuclei
    G4int sourceA = histoManager.GetInitialContidions().GetSourceA();
    G4int sourceAb = histoManager.GetInitialContidions().GetSourceAb();
    G4double rA_plus_rB = 1.3*(pow(G4double(sourceA),0.3333333)+pow(G4double(sourceAb),0.3333333));
  //Get nuclear mass
	  G4double init_nucl_mass_A = G4NucleiProperties::GetNuclearMass(histoManager.GetInitialContidions().GetSourceA(),histoManager.GetInitialContidions().GetSourceZ());
	  G4double init_nucl_mass_B = G4NucleiProperties::GetNuclearMass(histoManager.GetInitialContidions().GetSourceAb(),histoManager.GetInitialContidions().GetSourceZb());


//Starting of modeling
//###########################################################################################s#######################################
//Setting up ExcitationHandler
  G4ExcitationHandler* handler = new G4ExcitationHandler();
  DeexcitationHandler* handlerNew = new DeexcitationHandler(); //bad_alloc while initializing with new
  handler->SetEvaporation(theEvaporation);
  handler->SetMinEForMultiFrag(3*MeV);
  handler->SetMaxAandZForFermiBreakUp(19, 9); //problems when 12 and 6
    /*
    handlerNew->SetEvaporation(theEvaporation);
    handlerNew->SetMinEForMultiFrag(3*MeV);
    handlerNew->SetMaxAandZForFermiBreakUp(19, 8);*/
 
//Setting up Glauber code
  histoManager.CalcXsectNN();
  G4float omega = -1;
  G4float signn = histoManager.GetXsectNN();
  auto seed = static_cast<unsigned long int>(*CLHEP::HepRandom::getTheSeeds()); //setting the same seed to TGlauber
  //const unsigned long int seed = 0x7ffce81312f0;

  TGlauberMC *mcg=new TGlauberMC(histoManager.GetInitialContidions().GetSysA(),histoManager.GetInitialContidions().GetSysB(),signn,omega,seed);
  mcg->SetMinDistance(0);
  mcg->SetNodeDistance(0);
  mcg->SetCalcLength(0);
  mcg->SetCalcArea(0);
  mcg->SetCalcCore(0);
  mcg->SetDetail(99);

  if((histoManager.GetUpB() > 0 && histoManager.GetLowB() >= 0) && (histoManager.GetUpB() > histoManager.GetLowB()) && histoManager.GetUpB() < rA_plus_rB+7.5) {
      mcg->SetBmin(histoManager.GetLowB());
      mcg->SetBmax(histoManager.GetUpB());
     } 
  else {std::cout << " \n------>Calculation will be MB \n"<<std::endl; }
   
    //Setting Excitation Energy
    ExcitationEnergy* ExEnA = new ExcitationEnergy(histoManager.GetStatType(), sourceA);
    ExcitationEnergy* ExEnB = new ExcitationEnergy(histoManager.GetStatType(), sourceAb);
	
    if(histoManager.GetStatType()> 2){
    //Parameters for ALADIN parametrization
    G4double e_0=8*MeV;//MeV
    G4double sigma0 = 2	;//should fit our results
    G4double c0 = 2; // From Bondorf 1995
    G4double sigmaE0 = 1*MeV;
    G4double b0 = 0.1;
    G4double Pe = 24*MeV;
    G4double Pm = 0.2;
    ExEnA->SetParametersALADIN(e_0,sigma0,b0);
    ExEnB->SetParametersALADIN(e_0,sigma0,b0);
    ExEnA->SetParametersParabolicApproximation(Pe, Pm, sigma0, b0, 0.01);
    ExEnB->SetParametersParabolicApproximation(Pe, Pm, sigma0, b0, 0.01);
    //ExEnA->SetParametersCorrectedALADIN(0.01,1000,sigma0,b0,0);
    //ExEnB->SetParametersCorrectedALADIN(0.01,1000,sigma0,b0,0);
    }

    //Goldhaber model parameter
    G4double GoldhaberDev0 = 90*MeV; //Model parameter



for(G4int count=0;count<histoManager.GetIterations() ;count++){

      //auto remFile = remove("last_event.txt");
      //ofstream last_event("last_event.txt");
 
 id = count;

 //An event generated by GlauberMC is here.

  mcg->Run(1);

  histoManager.CalcNucleonDensity(mcg->GetNucleons(), mcg->GetB());

  //Side A$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  TGlauNucleus *nucA   = mcg->GetNucleusA(); 
  G4int NpartA = mcg->GetNpartA(); 
  G4int NpartB = mcg->GetNpartB();
  b = mcg->GetB();
  Ncoll = mcg->GetNcoll();
  Ncollpp = mcg->GetNcollpp();
  Ncollpn = mcg->GetNcollpn();
  Ncollnn = mcg->GetNcollnn();
  TObjArray* nucleons=mcg->GetNucleons();

  G4int A = 0;
  G4int Z = 0;
  G4int Ab = 0;
  G4int Zb = 0;
  G4int totBarNumA = 0;
  G4int totChargeNumA = 0;
  G4int totBarNumB = 0;
  G4int totChargeNumB = 0;

  for(G4int iArray = 0; iArray < nucleons->GetEntries(); iArray++){

  TGlauNucleon *nucleon=(TGlauNucleon*)(nucleons->At(iArray));
  if(nucleon->IsSpectator() && nucleon->IsInNucleusA()){A+=1;}
  if(nucleon->IsSpectator() && nucleon->IsInNucleusA() && nucleon->IsProton()){Z+=1;} 
  if(nucleon->IsWounded()   && nucleon->IsInNucleusA() && nucleon->IsProton()){totChargeNumA+=1;} 
  if(nucleon->IsSpectator() && nucleon->IsInNucleusB()){Ab+=1;}
  if(nucleon->IsSpectator() && nucleon->IsInNucleusB() && nucleon->IsProton()){Zb+=1;}
  if(nucleon->IsWounded()   && nucleon->IsInNucleusB() && nucleon->IsProton()){totChargeNumB+=1;} 
  }
  totBarNumA +=NpartA;
  totBarNumB +=NpartB;


 if(!(A == 0 && Ab ==0)) {//should be solved in some other way

    // last_event<<"A = "<<A<<", Z = "<<Z<<std::endl;
    // last_event<<"Ab = "<<Ab<<", Zb = "<<Zb<<std::endl;
    // last_event<<"b = "<<b<<std::endl;

//~~~~~~~~~~~~~~GoldhaberModel~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
     G4int thisEventNumFragments = 0;


     std::cout.setf(std::ios::scientific, std::ios::floatfield);

     //Excitation energy computing
     CLHEP::RandGauss randGauss(0, 1);
     CLHEP::RandFlat randFlat(new CLHEP::RanluxEngine);

     G4double energy = ExEnA->GetEnergy(A);

     G4double GoldhaberDev = GoldhaberDev0 * pow(G4double(A * NpartA) / G4double(sourceA - 1), 0.5);
     CLHEP::RandGauss randGoldhaber(new CLHEP::RanluxEngine, 0, GoldhaberDev);

     histoManager.GetHisto2(1)->Fill(energy / G4double(A), G4double(A) / G4double(sourceA));
     ExEn = energy / G4double(A);

     G4double p = randGauss.shoot() * GoldhaberDev;
     G4double theta = randFlat.shoot(1) * 2 * 3.141592;
     G4double px = p * sin(theta);
     G4double py = p * cos(theta);
     G4double pz = histoManager.GetInitialContidions().GetPzA()* G4double(A)/G4double(sourceA);
   
     G4double NuclearMass = G4NucleiProperties::GetNuclearMass(A, Z) + energy;
     G4LorentzVector p4(px, py, pz, pow(pow(NuclearMass, 2) + pow(pz, 2) + pow(p, 2), 0.5));
     G4Fragment aFragment(A, Z, p4);

     G4double eta_A = 0;


     G4ReactionProductVector *theProduct = handlerNew->BreakUp(aFragment);

     thisEventNumFragments = theProduct->size();

     histoManager.GetHisto(1)->Fill(thisEventNumFragments);

     if (Verbose) {
         std::cout << "### event  at " << energy / A << " MeV/nucleon"
                   << " with " << thisEventNumFragments << " particles  #####\n";
     }
    

     for (G4ReactionProductVector::iterator iVector = theProduct->begin(); iVector != theProduct->end(); ++iVector) {
         G4double thisFragmentZ = 0;
         G4double thisFragmentA = 0;

         const G4ParticleDefinition *pd = (*iVector)->GetDefinition();

         G4String particleEmitted = pd->GetParticleName();

         if (particleEmitted != "gamma" && particleEmitted != "e-" && particleEmitted != "e+") {
             thisFragmentZ = pd->GetAtomicNumber();
             thisFragmentA = pd->GetAtomicMass();
             if (pd->GetAtomicMass() == 0) { G4cout << "ERROR, pn = " << pd->GetParticleName() << G4endl; }
             MassOnSideA.push_back(thisFragmentA);
             ChargeOnSideA.push_back(thisFragmentZ);

             G4double eeA = (*iVector)->GetTotalEnergy();
             G4LorentzVector product_p4((*iVector)->GetMomentum().x(), (*iVector)->GetMomentum().y(),
                                        (*iVector)->GetMomentum().z(), eeA);
             G4double pXonA = product_p4.x() / MeV;
             G4double pYonA = product_p4.y() / MeV;
             G4double pZonA = product_p4.z() / MeV;
             p4 = p4 - product_p4;

             eta_A = 0.5 * log((std::sqrt(pXonA * pXonA + pYonA * pYonA + pZonA * pZonA) + pZonA) /
                               (std::sqrt(pXonA * pXonA + pYonA * pYonA + pZonA * pZonA) - pZonA));

             pXonSideA.push_back(pXonA);
             pYonSideA.push_back(pYonA);
             pZonSideA.push_back(pZonA);
             pseudorapidity_A.push_back(eta_A);


             if (thisFragmentZ == 0) {
                 histoManager.GetHisto2(3)->Fill(pXonA, pYonA);
                 histoManager.GetHisto(2)->Fill(pZonA);
             } else if (thisFragmentZ == 1 && thisFragmentA == 1) {
                 histoManager.GetHisto2(4)->Fill(pXonA, pYonA);
                 histoManager.GetHisto(3)->Fill(pZonA);
             } else if (thisFragmentZ < 20 && thisFragmentZ > 2) {
                 histoManager.GetHisto2(5)->Fill(pXonA, pYonA);
                 histoManager.GetHisto(4)->Fill(pZonA);
             } else {
                 histoManager.GetHisto2(6)->Fill(pXonA, pYonA);
                 histoManager.GetHisto(5)->Fill(pZonA);
             }
         }

         histoManager.GetHisto(6)->Fill(thisFragmentZ);
         histoManager.GetHisto(7)->Fill(thisFragmentA);
         histoManager.GetHisto2(2)->Fill(thisFragmentZ, thisFragmentA);
         totBarNumA += thisFragmentA;
         totChargeNumA += thisFragmentZ;
         delete (*iVector);
     }

     delete theProduct;

//$$$$$$$$$$$$$$$//Side B//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     TGlauNucleus *nucB = mcg->GetNucleusB();
     G4double energyB = 0;
     energyB = ExEnB->GetEnergy(Ab);

     G4double pxB = -px;
     G4double pyB = -py;
     G4double pzB = histoManager.GetInitialContidions().GetPzB()* G4double(Ab) /G4double(sourceAb);
   
     G4double NuclearMassB = G4NucleiProperties::GetNuclearMass(Ab, Zb) + energyB;
     G4LorentzVector p4b(pxB, pyB, pzB, pow(pow(NuclearMassB, 2) + pow(pxB, 2) + pow(pyB, 2) + pow(pzB, 2), 0.5));
     G4Fragment aFragmentB(Ab, Zb, p4b);
     G4double beta_z = std::sqrt(pow(pzB / NuclearMassB, 2) / (1 + pow(pzB / NuclearMassB, 2)));
     G4double beta_y = std::sqrt(pow(pyB / NuclearMassB, 2) / (1 + pow(pyB / NuclearMassB, 2)));
     G4double beta_x = std::sqrt(pow(pxB / NuclearMassB, 2) / (1 + pow(pxB / NuclearMassB, 2)));
     G4float eta_B = 0;


     G4ReactionProductVector *theProductB = handlerNew->BreakUp(aFragment);

     for (G4ReactionProductVector::iterator kVector = theProductB->begin(); kVector != theProductB->end(); ++kVector) {
         G4int thisFragmentZb = 0;
         G4int thisFragmentAb = 0;

         const G4ParticleDefinition *pdB = (*kVector)->GetDefinition();

         G4String particleEmittedB = pdB->GetParticleName();

         if (particleEmittedB != "gamma" && particleEmittedB != "e-" && particleEmittedB != "e+") {
             thisFragmentZb = pdB->GetAtomicNumber();
             thisFragmentAb = pdB->GetAtomicMass();
             MassOnSideB.push_back(thisFragmentAb);
             ChargeOnSideB.push_back(thisFragmentZb);

             G4double eeB = (*kVector)->GetTotalEnergy();
             G4LorentzVector product_p4b((*kVector)->GetMomentum().x(), (*kVector)->GetMomentum().y(),
                                         (*kVector)->GetMomentum().z(), eeB);
             G4double pXonB = product_p4b.x() / MeV;
             G4double pYonB = product_p4b.y() / MeV;
             G4double pZonB = product_p4b.z() / MeV;
             p4b = p4b - product_p4b;

             eta_B = 0.5 * log((std::sqrt(pXonB * pXonB + pYonB * pYonB + pZonB * pZonB) + pZonB) /
                               (std::sqrt(pXonB * pXonB + pYonB * pYonB + pZonB * pZonB) - pZonB));
             pXonSideB.push_back(pXonB);
             pYonSideB.push_back(pYonB);
             pZonSideB.push_back(pZonB);
             pseudorapidity_B.push_back(eta_B);

             totBarNumB += thisFragmentAb;
             totChargeNumB += thisFragmentZb;

         }


         histoManager.GetHisto(0)->Fill(thisFragmentZb);
         delete (*kVector);
     }
    
     delete theProductB;
     //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

     histoManager.GetTree()->Fill();
     MassOnSideA.clear();
     MassOnSideB.clear();
     ChargeOnSideB.clear();
     ChargeOnSideA.clear();
     pXonSideA.clear();
     pXonSideB.clear();
     pYonSideA.clear();
     pYonSideB.clear();
     pZonSideA.clear();
     pZonSideB.clear();
     pseudorapidity_A.clear();
     pseudorapidity_B.clear();

     if (!G4bool(count % 10)) { G4cout << "Program is working," << count << " events calculated    \r" << std::flush; }

     auto seed_now = static_cast<unsigned long int>(*CLHEP::HepRandom::getTheSeeds());
    // last_event << "Current seed = " << CLHEP::HepRandom::getTheSeeds() << std::endl;
    // last_event.close();
    }
  }
G4cout<<"----> collided "<<histoManager.GetIterations()<<" nuclei "<<histoManager.GetSysA()<< " with " << histoManager.GetSysB() <<" at N-N x-section "<<signn<<" mb"<<G4endl;
G4cout<<"----> total x-sect = "<<mcg->GetTotXSect()<< " +- " << mcg->GetTotXSectErr() <<" b";
      
  histoManager.FillConditionsTree(mcg->GetTotXSect());
  histoManager.CleanHisto();

  delete runManager;
  delete theEvaporation;
  delete handler;
  delete handlerNew;
  delete mcg;
  delete ExEnA;
  delete ExEnB;
  return 0;
}
