

#include "jpsi.h"
//#include "caloheader.h"
//#include "clusterizer.cxx"

void jpsi(
	TString inFile            	= "../rootfiles/jpsi_r1.root",
	TString inFileGeometry      = "../rootfiles/geometry_ideal_CEMCproj_EEMC_LFHCAL_corrected.root",
	bool do_reclus              = true,
    unsigned short primaryTrackSource = 0

){

	// load tree
    TChain *const tt_event = new TChain("event_tree");
    if (inFile.EndsWith(".root")) {                     // are we loading a single root tree?
        std::cout << "loading a single root file" << std::endl;
        tt_event->AddFile(inFile);
    }
    else {                                              // or are we loading a bunch?
        std::cout << "loading a list of files" << std::endl;
        std::ifstream files(inFile);
        std::string filePath;

        while (std::getline(files, filePath)) {
            tt_event->AddFile(filePath.c_str());
        }
        files.close();
    }
    if(!tt_event){ std::cout << "tree not found... returning!"<< std::endl; return;}

     // // load geometry tree
    /*tt_geometry =  (TTree *) (new TFile(inFileGeometry.Data(), "READ"))->Get("geometry_tree");
    if(!tt_geometry){ cout << "geometry tree not found... returning!"<< endl; return;}*/

    Long64_t nEntriesTree                 = tt_event->GetEntries();
    std::cout << "Number of events in tree: " << nEntriesTree << std::endl;

    SetBranchAddressesTree(tt_event);
    //SetBranchAddressesGeometryTree(tt_geometry);
    //SetGeometryIndices();

    double elect_energy = 10;
    double proton_energy = 100;

    TFile *MyFile = new TFile("jpsi_test.root","RECREATE");
	TTree *EvTree = new TTree( "T", "T" );

	double JpsimassMC, JpsithetaMC, JpsiphiMC, JpsietaMC, JpsipMC, JpsiEMC;
	double epthetaMC, epphiMC, eppTMC, epetaMC, eppMC, epEMC;
	double emthetaMC, emphiMC, empTMC, emetaMC, empMC, emEMC;
	double pthetaMC, pphiMC, ppTMC, petaMC, ppMC, pEMC;
	double elthetaMC, elphiMC, elpTMC, eletaMC, elpMC, elEMC;
	double Q2MClepton, xbjMClepton, nuMClepton, thetaelMClepton, yMClepton, sMClepton;
	double Q2MChadron, xbjMChadron, nuMChadron, thetaelMChadron, yMChadron, sMChadron;
	double Jpsi_track3_M1_etap1, Jpsi_track3_M1_etap2, Jpsi_track3_M2_etap1, Jpsi_track3_M2_etap2;
	double Jpsi_track3_M3_etap1, Jpsi_track3_M3_etap2;
	double Jpsi_track3_M1_E1, Jpsi_track3_M1_E2;

	int nTracks;
	int n_mass;
	int trackID[20];
	

	double Jpsi_track2_M, Jpsi_track3_M1, Jpsi_track3_M2, Jpsi_track3_M3;
	double Jpsi_track3_eta1, Jpsi_track3_eta2, Jpsi_track3_eta3;
	double Jpsi_track2_etap1, Jpsi_track2_etap2;
	double Q2_track3_lepton, xbj_track3_lepton;
	double Jpsi_track[10];

	EvTree->Branch("JpsimassMC",&JpsimassMC,"JpsimassMC/D");
	EvTree->Branch("JpsithetaMC",&JpsithetaMC,"JpsithetaMC/D");
	EvTree->Branch("JpsiphiMC",&JpsiphiMC,"JpsiphiMC/D");
	EvTree->Branch("JpsietaMC",&JpsietaMC,"JpsietaMC/D");
	EvTree->Branch("JpsipMC",&JpsipMC,"JpsipMC/D");
	EvTree->Branch("JpsiEMC",&JpsiEMC,"JpsiEEMC/D");

	EvTree->Branch("epthetaMC",&epthetaMC,"epthetaMC/D");
	EvTree->Branch("epphiMC",&epphiMC,"epphiMC/D");
	EvTree->Branch("epetaMC",&epetaMC,"epetaMC/D");
	EvTree->Branch("eppTMC",&eppTMC,"eppTMC/D");
	EvTree->Branch("eppMC",&eppMC,"eppMC/D");
	EvTree->Branch("epEMC",&epEMC,"epEMC/D");

	EvTree->Branch("emthetaMC",&emthetaMC,"emthetaMC/D");
	EvTree->Branch("emphiMC",&emphiMC,"emphiMC/D");
	EvTree->Branch("emetaMC",&emetaMC,"emetaMC/D");
	EvTree->Branch("empTMC",&empTMC,"empTMC/D");
	EvTree->Branch("empMC",&empMC,"empMC/D");
	EvTree->Branch("emEMC",&emEMC,"emEMC/D");

	EvTree->Branch("pthetaMC",&pthetaMC,"pthetaMC/D");
	EvTree->Branch("pphiMC",  &pphiMC,  "pphiMC/D");
	EvTree->Branch("petaMC",  &petaMC,  "petaMC/D");
	EvTree->Branch("ppTMC",   &ppTMC,   "ppTMC/D");
	EvTree->Branch("ppMC",    &ppMC,    "ppMC/D");
	EvTree->Branch("pEMC",    &pEMC,    "pEMC/D");

	EvTree->Branch("elthetaMC",&elthetaMC,"elthetaMC/D");
	EvTree->Branch("elphiMC",  &elphiMC,  "elphiMC/D");
	EvTree->Branch("eletaMC",  &eletaMC,  "eletaMC/D");
	EvTree->Branch("elpTMC",   &elpTMC,   "elpTMC/D");
	EvTree->Branch("elpMC",    &elpMC,    "elpMC/D");
	EvTree->Branch("elEMC",    &elEMC,    "elEMC/D");

	EvTree->Branch("Q2MClepton",&Q2MClepton,"Q2MClepton/D");
	EvTree->Branch("xbjMClepton",&xbjMClepton,"xbjMClepton/D");
	EvTree->Branch("nuMClepton",&nuMClepton,"nuMClepton/D");
	EvTree->Branch("thetaelMClepton",&thetaelMClepton,"thetaelMClepton/D");
	EvTree->Branch("yMClepton",&yMClepton,"yMClepton/D");
	EvTree->Branch("sMClepton",&sMClepton,"sMClepton/D");

	EvTree->Branch("Q2MChadron",&Q2MChadron,"Q2MChadron/D");
	EvTree->Branch("xbjMChadron",&xbjMChadron,"xbjMChadron/D");
	EvTree->Branch("nuMChadron",&nuMChadron,"nuMChadron/D");
	EvTree->Branch("thetaelMChadron",&thetaelMChadron,"thetaelMChadron/D");
	EvTree->Branch("yMChadron",&yMChadron,"yMChadron/D");
	EvTree->Branch("sMChadron",&sMChadron,"sMChadron/D");

	EvTree->Branch("nTracks",&nTracks,"nTracks/I");
    EvTree->Branch("trackID",trackID,"trackID[nTracks]/I");
	EvTree->Branch("Jpsi_track2_M",&Jpsi_track2_M,"Jpsi_track2_M/D");
	
	EvTree->Branch("Jpsi_track2_etap1",&Jpsi_track2_etap1,"Jpsi_track2_etap1/D");
	EvTree->Branch("Jpsi_track2_etap2",&Jpsi_track2_etap2,"Jpsi_track2_etap2/D");


	EvTree->Branch("Jpsi_track3_M1",&Jpsi_track3_M1,"Jpsi_track3_M1/D");
	EvTree->Branch("Jpsi_track3_M2",&Jpsi_track3_M2,"Jpsi_track3_M2/D");
	EvTree->Branch("Jpsi_track3_M3",&Jpsi_track3_M3,"Jpsi_track3_M3/D");

	EvTree->Branch("Jpsi_track3_M1_etap1",&Jpsi_track3_M1_etap1,"Jpsi_track3_M1_etap1/D");
	EvTree->Branch("Jpsi_track3_M1_etap2",&Jpsi_track3_M1_etap2,"Jpsi_track3_M1_etap2/D");
	EvTree->Branch("Jpsi_track3_M1_E1",&Jpsi_track3_M1_E1,"Jpsi_track3_M1_E1/D");
	EvTree->Branch("Jpsi_track3_M1_E2",&Jpsi_track3_M1_E2,"Jpsi_track3_M1_E2/D");

	EvTree->Branch("Jpsi_track3_M2_etap1",&Jpsi_track3_M2_etap1,"Jpsi_track3_M2_etap1/D");
	EvTree->Branch("Jpsi_track3_M2_etap2",&Jpsi_track3_M2_etap2,"Jpsi_track3_M2_etap2/D");
	EvTree->Branch("Jpsi_track3_M3_etap1",&Jpsi_track3_M3_etap1,"Jpsi_track3_M3_etap1/D");
	EvTree->Branch("Jpsi_track3_M3_etap2",&Jpsi_track3_M3_etap2,"Jpsi_track3_M3_etap2/D");
	
	EvTree->Branch("Q2_track3_lepton",&Q2_track3_lepton,"Q2_track3_lepton/D");

   


	double CEMC_cluster_E_ep, CEMC_cluster_seed_ep, CEMC_cluster_Eta_ep, CEMC_cluster_Phi_ep; 
	double CEMC_cluster_Z_ep, CEMC_cluster_X_ep, CEMC_cluster_Y_ep, CEMC_cluster_M02_ep;
	double CEMC_cluster_M20_ep, CEMC_cluster_isMatched_ep, CEMC_cluster_NTowers_ep, CEMC_cluster_trueID_ep;
	double CEMC_cluster_NtrueID_ep;
	double CEMC_cluster_E_em, CEMC_cluster_seed_em, CEMC_cluster_Eta_em, CEMC_cluster_Phi_em; 
	double CEMC_cluster_Z_em, CEMC_cluster_X_em, CEMC_cluster_Y_em, CEMC_cluster_M02_em;
	double CEMC_cluster_M20_em, CEMC_cluster_isMatched_em, CEMC_cluster_NTowers_em, CEMC_cluster_trueID_em;
	double CEMC_cluster_NtrueID_em;
	double CEMC_cluster_E_p, CEMC_cluster_seed_p, CEMC_cluster_Eta_p, CEMC_cluster_Phi_p; 
	double CEMC_cluster_Z_p, CEMC_cluster_X_p, CEMC_cluster_Y_p, CEMC_cluster_M02_p;
	double CEMC_cluster_M20_p, CEMC_cluster_isMatched_p, CEMC_cluster_NTowers_p, CEMC_cluster_trueID_p;
	double CEMC_cluster_NtrueID_p;
/*
	if(kCEMC){
		EvTree->Branch("CEMC_cluster_E_ep",&CEMC_cluster_E_ep,"CEMC_cluster_E_ep/D");
		EvTree->Branch("CEMC_cluster_seed_ep",&CEMC_cluster_seed_ep,"CEMC_cluster_seed_ep/D");
		EvTree->Branch("CEMC_cluster_Eta_ep",&CEMC_cluster_Eta_ep,"CEMC_cluster_Eta_ep/D");
		EvTree->Branch("CEMC_cluster_Phi_ep",&CEMC_cluster_Phi_ep,"CEMC_cluster_Phi_ep/D");
		EvTree->Branch("CEMC_cluster_Z_ep",&CEMC_cluster_Z_ep,"CEMC_cluster_Z_ep/D");
		EvTree->Branch("CEMC_cluster_X_ep",&CEMC_cluster_X_ep,"CEMC_cluster_X_ep/D");
		EvTree->Branch("CEMC_cluster_Y_ep",&CEMC_cluster_Y_ep,"CEMC_cluster_Y_ep/D");
		EvTree->Branch("CEMC_cluster_M02_ep",&CEMC_cluster_M02_ep,"CEMC_cluster_M02_ep/D");
		EvTree->Branch("CEMC_cluster_M20_ep",&CEMC_cluster_M20_ep,"CEMC_cluster_M20_ep/D");
		EvTree->Branch("CEMC_cluster_isMatched_ep",&CEMC_cluster_isMatched_ep,"CEMC_cluster_isMatched_ep/D");
		EvTree->Branch("CEMC_cluster_NTowers_ep",&CEMC_cluster_NTowers_ep,"CEMC_cluster_NTowers_ep/D");
		EvTree->Branch("CEMC_cluster_trueID_ep",&CEMC_cluster_trueID_ep,"CEMC_cluster_trueID_ep/D");
		EvTree->Branch("CEMC_cluster_NtrueID_ep",&CEMC_cluster_NtrueID_ep,"CEMC_cluster_NtrueID_ep/D");


		EvTree->Branch("CEMC_cluster_E_em",&CEMC_cluster_E_em,"CEMC_cluster_E_em/D");
		EvTree->Branch("CEMC_cluster_seed_em",&CEMC_cluster_seed_em,"CEMC_cluster_seed_em/D");
		EvTree->Branch("CEMC_cluster_Eta_em",&CEMC_cluster_Eta_em,"CEMC_cluster_Eta_em/D");
		EvTree->Branch("CEMC_cluster_Phi_em",&CEMC_cluster_Phi_em,"CEMC_cluster_Phi_em/D");
		EvTree->Branch("CEMC_cluster_Z_em",&CEMC_cluster_Z_em,"CEMC_cluster_Z_em/D");
		EvTree->Branch("CEMC_cluster_X_em",&CEMC_cluster_X_em,"CEMC_cluster_X_em/D");
		EvTree->Branch("CEMC_cluster_Y_em",&CEMC_cluster_Y_em,"CEMC_cluster_Y_em/D");
		EvTree->Branch("CEMC_cluster_M02_em",&CEMC_cluster_M02_em,"CEMC_cluster_M02_em/D");
		EvTree->Branch("CEMC_cluster_M20_em",&CEMC_cluster_M20_em,"CEMC_cluster_M20_em/D");
		EvTree->Branch("CEMC_cluster_isMatched_em",&CEMC_cluster_isMatched_em,"CEMC_cluster_isMatched_em/D");
		EvTree->Branch("CEMC_cluster_NTowers_em",&CEMC_cluster_NTowers_em,"CEMC_cluster_NTowers_em/D");
		EvTree->Branch("CEMC_cluster_trueID_em",&CEMC_cluster_trueID_em,"CEMC_cluster_trueID_em/D");
		EvTree->Branch("CEMC_cluster_NtrueID_em",&CEMC_cluster_NtrueID_em,"CEMC_cluster_NtrueID_em/D");

		EvTree->Branch("CEMC_cluster_E_p",&CEMC_cluster_E_p,"CEMC_cluster_E_p/D");
		EvTree->Branch("CEMC_cluster_seed_p",&CEMC_cluster_seed_p,"CEMC_cluster_seed_p/D");
		EvTree->Branch("CEMC_cluster_Eta_p",&CEMC_cluster_Eta_p,"CEMC_cluster_Eta_p/D");
		EvTree->Branch("CEMC_cluster_Phi_p",&CEMC_cluster_Phi_p,"CEMC_cluster_Phi_p/D");
		EvTree->Branch("CEMC_cluster_Z_p",&CEMC_cluster_Z_p,"CEMC_cluster_Z_p/D");
		EvTree->Branch("CEMC_cluster_X_p",&CEMC_cluster_X_p,"CEMC_cluster_X_p/D");
		EvTree->Branch("CEMC_cluster_Y_p",&CEMC_cluster_Y_p,"CEMC_cluster_Y_p/D");
		EvTree->Branch("CEMC_cluster_M02_p",&CEMC_cluster_M02_p,"CEMC_cluster_M02_p/D");
		EvTree->Branch("CEMC_cluster_M20_p",&CEMC_cluster_M20_p,"CEMC_cluster_M20_p/D");
		EvTree->Branch("CEMC_cluster_isMatched_p",&CEMC_cluster_isMatched_p,"CEMC_cluster_isMatched_p/D");
		EvTree->Branch("CEMC_cluster_NTowers_p",&CEMC_cluster_NTowers_p,"CEMC_cluster_NTowers_p/D");
		EvTree->Branch("CEMC_cluster_trueID_p",&CEMC_cluster_trueID_p,"CEMC_cluster_trueID_p/D");
		EvTree->Branch("CEMC_cluster_NtrueID_p",&CEMC_cluster_NtrueID_p,"CEMC_cluster_NtrueID_p/D");

	}*/

   	_nEventsTree=0;

    // main event loop
    for (Long64_t i=0; i<nEntriesTree;i++) {

    	Jpsi_track2_M = Jpsi_track3_M1 = Jpsi_track3_M2 = Jpsi_track3_M3 = -100;
    	Jpsi_track3_eta1  = Jpsi_track3_eta2  = Jpsi_track3_eta3 = -100;
    	Jpsi_track2_etap1 = Jpsi_track2_etap2 = -100;
    	Q2_track3_lepton  = xbj_track3_lepton = -100;
    	Jpsi_track3_M1_E1 = Jpsi_track3_M1_E2 = -100;

    	// load current event
        tt_event->GetEntry(i);
        if(_nMCPart!=4) continue;
        //if(_nMCPart!=2) continue;
        TLorentzVector epMC;
        TLorentzVector emMC;
        TLorentzVector pMC;
        TLorentzVector eMC;

        for(int imc=0; imc<_nMCPart; imc++){
        //if (_mcpart_PDG[imc] == -11){
        	 //===== J/Psi True
        	 if(_mcpart_PDG[imc] == 11 && _mcpart_BCID[imc]==10009)  {
        	 	epMC[0] = _mcpart_px[imc];
        	 	epMC[1]	= _mcpart_py[imc];
        	 	epMC[2]	= _mcpart_pz[imc];
        	 	epMC[3]	= _mcpart_E[imc];
        	 }
        	 else if(_mcpart_PDG[imc] == -11 && _mcpart_BCID[imc]==10008)  {
        	 	emMC[0] = _mcpart_px[imc];
        	 	emMC[1]	= _mcpart_py[imc];
        	 	emMC[2]	= _mcpart_pz[imc];
        	 	emMC[3]	= _mcpart_E[imc];
        	 }
        	 else if(_mcpart_PDG[imc] == 2212 && _mcpart_BCID[imc]==10007)  {
        	 	pMC[0] = _mcpart_px[imc];
        	 	pMC[1]	= _mcpart_py[imc];
        	 	pMC[2]	= _mcpart_pz[imc];
        	 	pMC[3]	= _mcpart_E[imc];
        	 }
        	 else if(_mcpart_PDG[imc] == 11 && _mcpart_BCID[imc]==10002)  {
        	 	eMC[0] = _mcpart_px[imc];
        	 	eMC[1]	= _mcpart_py[imc];
        	 	eMC[2]	= _mcpart_pz[imc];
        	 	eMC[3]	= _mcpart_E[imc];
        	 }
        	 else
        	 {
        	 	continue;
        	 }
        }


       	TLorentzVector JpsiMC = epMC + emMC;
        JpsimassMC  = JpsiMC.M();
      	JpsithetaMC = JpsiMC.Theta()*180/TMath::Pi();
      	JpsiphiMC   = JpsiMC.Phi()*180/TMath::Pi();
        JpsietaMC   = -TMath::Log(TMath::Tan(JpsiMC.Theta()/2));
        JpsipMC 	=  JpsiMC.Vect().Mag();
        JpsiEMC     =  JpsiMC.E();
        //====== ep
        epphiMC     = epMC.Phi()*180/TMath::Pi();
      	epthetaMC   = epMC.Theta()*180/TMath::Pi();
        epetaMC 	= -TMath::Log(TMath::Tan(epMC.Theta()/2));
        eppTMC		= epMC.Perp();
        eppMC       = epMC.Vect().Mag();
        epEMC       = epMC.E();

        //====== em
       	emphiMC     = emMC.Phi()*180/TMath::Pi();
      	emthetaMC   = emMC.Theta()*180/TMath::Pi();
      	emetaMC 	 = -TMath::Log(TMath::Tan(emMC.Theta()/2));
        empTMC		 = emMC.Perp();
        empMC       = emMC.Vect().Mag();
        emEMC       = emMC.E();

        //====== p
        pphiMC     = pMC.Phi()*180/TMath::Pi();
        pthetaMC   = pMC.Theta()*180/TMath::Pi();
        petaMC 	= -TMath::Log(TMath::Tan(pMC.Theta()/2));
        ppTMC		= pMC.Perp();
        ppMC       = pMC.Vect().Mag();
        pEMC       = pMC.E();

        //====== e
      	elphiMC     = eMC.Phi()*180/TMath::Pi();
      	elthetaMC   = eMC.Theta()*180/TMath::Pi();
      	eletaMC 	= -TMath::Log(TMath::Tan(eMC.Theta()/2));
      	elpTMC		=  eMC.Perp();
      	elpMC       = eMC.Vect().Mag();
        elEMC       = eMC.E();

        //============Rotating to the head on collisions
        TLorentzVector JpsiMC_HO = rotate_reco(JpsiMC); 
        TLorentzVector pMC_HO = rotate_reco(pMC);
        TLorentzVector eMC_HO = rotate_reco(eMC);

        TLorentzVector init_e(0,0,-elect_energy,elect_energy);
		TLorentzVector init_p(0,0,proton_energy,proton_energy);        	 

		//============ Kinematic Variables
		//============ Only Lepton Information
		nuMClepton        = init_e.E() - eMC_HO.E();
		thetaelMClepton   = init_e.Vect().Angle(eMC_HO.Vect());
		Q2MClepton  = -(init_e-eMC_HO)*(init_e-eMC_HO);
		//xbjMClepton = Q2MClepton/(2*init_e*(init_e-eMC_HO));
		yMClepton   = 1 - ((eMC_HO.E()*(1-pow(cos(thetaelMClepton/2),2)))/(2*init_e.E()));	 
		sMClepton   = sqrt((init_e + init_p)*(init_e + init_p));
		xbjMClepton = Q2MClepton/(sMClepton*yMClepton);


		//============ JB method
		Double_t num = (JpsiMC_HO.E() - JpsiMC_HO.Pz()) + (pMC_HO.E() - pMC_HO.Pz()) + (init_p.E() - init_p.Pz());

		yMChadron 	  = num/(2*init_e.E());
			 
		TLorentzVector PT  = JpsiMC_HO + pMC_HO + init_p;

		Q2MChadron   = PT.Perp()/(1-yMChadron);
			 
		sMChadron    = sqrt((JpsiMC_HO + pMC_HO + eMC_HO)*(JpsiMC_HO + pMC_HO + eMC_HO));
			 
		xbjMChadron  = Q2MChadron/(sMChadron*yMChadron);

		nTracks = _nTracks;
		n_mass = 0;

		for (int l = 0; l<nTracks; l++){
			trackID[l] = _mcpart_BCID[l];
		}
		Jpsi_track2_M -100;
		if(_nTracks == 2){

			TVector3 v1(_track_px[0],_track_px[0],_track_pz[0]);
			TVector3 v2(_track_px[1],_track_px[1],_track_pz[1]);
			//======= test if are J/Psi
			TLorentzVector hyp1;
			TLorentzVector hyp2;

			hyp1.SetXYZM(_track_px[0],_track_px[0],_track_pz[0], me);
			hyp2.SetXYZM(_track_px[1],_track_px[1],_track_pz[1], me);
			TLorentzVector JPsihyp = hyp1 + hyp2;

			Jpsi_track2_M = JPsihyp.M();
			Jpsi_track2_etap1 = JPsihyp.Eta();
			Jpsi_track2_etap2 = JPsihyp.Eta();
		}
		else if(_nTracks ==3){

			TLorentzVector hyp1;
			TLorentzVector hyp2;
			TLorentzVector hyp3;

			hyp1.SetXYZM(_track_px[0],_track_px[0],_track_pz[0], me);
			hyp2.SetXYZM(_track_px[1],_track_px[1],_track_pz[1], me);
			hyp3.SetXYZM(_track_px[2],_track_px[2],_track_pz[2], me);

			TLorentzVector JPsihyp1 = hyp1 + hyp2;
			TLorentzVector JPsihyp2 = hyp1 + hyp3;
			TLorentzVector JPsihyp3 = hyp2 + hyp3;

	
			Jpsi_track3_M1 = JPsihyp1.M();
			Jpsi_track3_eta1 = JPsihyp1.Eta();
			Jpsi_track3_M1_etap1 = hyp1.Eta();
			Jpsi_track3_M1_etap2 = hyp2.Eta();
			Jpsi_track3_M1_E1    = hyp1.E();
			Jpsi_track3_M1_E2    = hyp2.E();

			Jpsi_track3_M2 = JPsihyp2.M();
			Jpsi_track3_eta2 = JPsihyp2.Eta();
			Jpsi_track3_M2_etap1 = hyp1.Eta();
			Jpsi_track3_M2_etap2 = hyp3.Eta();

			Jpsi_track3_M3 = JPsihyp3.M();
			Jpsi_track3_M3_etap1 = hyp2.Eta();
			Jpsi_track3_M3_etap2 = hyp3.Eta();

			//========== Calculating kinematics
			TLorentzVector hyp3_HO =  rotate_reco(hyp3);
			
			Q2_track3_lepton  = -(init_e-hyp3_HO)*(init_e-hyp3_HO);
			//xbj_track3_lepton = Q2MClepton/(2*init_e*(init_e-hyp3_HO));

			nuMClepton        = init_e.E() - eMC_HO.E();
			
			double thetael_track3_lepton   = init_e.Vect().Angle(hyp3_HO.Vect());
		
			double y_track3_lepton   = 1 - ((hyp3_HO.E()*(1-pow(cos(thetael_track3_lepton/2),2)))/(2*init_e.E()));	 
		    double s_track3_lepton   = sqrt((init_e + init_p)*(init_e + init_p));
		
			xbj_track3_lepton = Q2_track3_lepton/(s_track3_lepton*y_track3_lepton);
			

	

		}
 
        	 
    	
		if (Q2MClepton<5.) continue;


    	EvTree->Fill();
    	//if (i>10000) break;


    }

    EvTree->Write();
 	MyFile->Close();
	
}

