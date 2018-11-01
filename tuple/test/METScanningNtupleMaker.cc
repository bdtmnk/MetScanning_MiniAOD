#include "MetScanning/tuple/test/METScanningNtupleMaker.h"
#include <iostream>
//User
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//Include Generator Level information
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <DataFormats/HLTReco/interface/TriggerTypeDefs.h>
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
//#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

//#include "DesyTauAnalyses/NTupleMaker/interface/idAlgos.h"


METScanningNtupleMaker::METScanningNtupleMaker(const edm::ParameterSet& iConfig):
    cdata(iConfig.getUntrackedParameter<bool>("IsData", false)),
    crecmuon(iConfig.getUntrackedParameter<bool>("RecMuon", false)),
    crecpuppimet(iConfig.getUntrackedParameter<bool>("RecPuppiMet", false)),
    crecpfmet(iConfig.getUntrackedParameter<bool>("RecPFMet", false)),
    crecprimvertex(iConfig.getUntrackedParameter<bool>("RecPrimVertex", false)),
    cFlags(iConfig.getUntrackedParameter<vector<string> >("Flags")),
    cFlagsProcesses(iConfig.getUntrackedParameter<vector<string> >("FlagsProcesses")),

    PuppiMetCollectionToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("PuppiMetCollectionTag"))),
    MetCollectionToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("MetCollectionTag"))),
    cHLTriggerPaths(iConfig.getUntrackedParameter<vector<string> >("HLTriggerPaths")),
    cTriggerProcess(iConfig.getUntrackedParameter<string>("TriggerProcess", "HLT")),
    MuonCollectionToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("MuonCollectionTag"))),
    PVToken_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("PVCollectionTag")))

{
 
  //TO MINIAOD
  //LD
  consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", cTriggerProcess));
  for(std::vector<string>::iterator it = cFlagsProcesses.begin();
      it != cFlagsProcesses.end(); it++){
    consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", it->data()));
  }



}

void METScanningNtupleMaker::beginJob(){

  // The root tuple
  tf1 = new TFile("tuple.root", "RECREATE");
  tree = new TTree("tree","tree");

  //basic informations ==========================
  tree->Branch("run",&run,"run/l");
  tree->Branch("lumi",&lumiBlock,"lumi/l");
  tree->Branch("event",&event,"event/l");
  tree->Branch("time",&time,"time/l");

  flags_ = new std::map<std::string, int>();
  tree->Branch("flags", "std::map<std::string, int>", &flags_);


  //MET Tokens:
	//PF MET
  tree->Branch("pfmet_ex", &pfmet_ex, "pfmet_ex/F");
    tree->Branch("pfmet_ey", &pfmet_ey, "pfmet_ey/F");
    tree->Branch("pfmet_ez", &pfmet_ey, "pfmet_ez/F");
    tree->Branch("pfmet_pt", &pfmet_pt, "pfmet_pt/F");
    tree->Branch("pfmet_phi", &pfmet_phi, "pfmet_phi/F");
    tree->Branch("pfmet_sigxx", &pfmet_sigxx, "pfmet_sigxx/F");
    tree->Branch("pfmet_sigxy", &pfmet_sigxy, "pfmet_sigxy/F");
    tree->Branch("pfmet_sigyx", &pfmet_sigyx, "pfmet_sigyx/F");
    tree->Branch("pfmet_sigyy", &pfmet_sigyy, "pfmet_sigyy/F");
    tree->Branch("pfmet_sig", &pfmet_sig, "pfmet_sig/F");




  //Puppi MET:
  tree->Branch("puppimet_ex", &puppimet_ex, "puppimet_ex/F");
  tree->Branch("puppimet_ey", &puppimet_ey, "puppimet_ey/F");
  tree->Branch("puppimet_ez", &puppimet_ey, "puppimet_ez/F");
  tree->Branch("puppimet_pt", &puppimet_pt, "puppimet_pt/F");
  tree->Branch("puppimet_phi", &puppimet_phi, "puppimet_phi/F");
  tree->Branch("puppimet_sigxx", &puppimet_sigxx, "puppimet_sigxx/F");
  tree->Branch("puppimet_sigxy", &puppimet_sigxy, "puppimet_sigxy/F");
  tree->Branch("puppimet_sigyx", &puppimet_sigyx, "puppimet_sigyx/F");
  tree->Branch("puppimet_sigyy", &puppimet_sigyy, "puppimet_sigyy/F");


  //MUONS:
  tree->Branch("muon_count", &muon_count, "muon_count/i");
  tree->Branch("muon_px", muon_px, "muon_px[muon_count]/F");
  tree->Branch("muon_py", muon_py, "muon_py[muon_count]/F");
  tree->Branch("muon_pz", muon_pz, "muon_pz[muon_count]/F");
  tree->Branch("muon_pt", muon_pt, "muon_pt[muon_count]/F");
  tree->Branch("muon_eta", muon_eta, "muon_eta[muon_count]/F");
  tree->Branch("muon_phi", muon_phi, "muon_phi[muon_count]/F");
  tree->Branch("muon_pterror", muon_pterror, "muon_pterror[muon_count]/F");
  //tree->Branch("muon_segmentComp", muon_segmentComp, "muon_segmentComp[muon_count]/F");
  //tree->Branch("muon_isPF",muon_isPF,"muon_isPF[muon_count]/O");
  //tree->Branch("muon_isGlobal",muon_isGlobal,"muon_isGlobal[muon_count]/O");
  tree->Branch("muon_isTracker",muon_isTracker,"muon_isTracker[muon_count]/O");
  //tree->Branch("muon_isTight",muon_isTight,"muon_isTight[muon_count]/O");
  //tree->Branch("muon_isLoose",muon_isLoose,"muon_isLoose[muon_count]/O");
  //tree->Branch("muon_isMedium",muon_isMedium,"muon_isMedium[muon_count]/O");
  //tree->Branch("muon_isICHEP",muon_isICHEP,"muon_isICHEP[muon_count]/O");
  //tree->Branch("muon_genmatch", muon_genmatch, "muon_genmatch[muon_count]/I");
  //tree->Branch("muon_isDuplicate",muon_isDuplicate,"muon_isDuplicate[muon_count]/O");
  //tree->Branch("muon_isBad",muon_isBad,"muon_isBad[muon_count]/O");

  //N Vertex:
  tree->Branch("primvertex_count", &primvertex_count, "primvertex_count/i");
  tree->Branch("goodprimvertex_count", &goodprimvertex_count, "goodprimvertex_count/i");
  tree->Branch("primvertex_x", &primvertex_x, "primvertex_x/F");
  tree->Branch("primvertex_y", &primvertex_y, "primvertex_y/F");
  tree->Branch("primvertex_z", &primvertex_z, "primvertex_z/F");
  tree->Branch("primvertex_chi2", &primvertex_chi2, "primvertex_chi2/F");
  tree->Branch("primvertex_ndof", &primvertex_ndof, "primvertex_ndof/F");
  tree->Branch("primvertex_ptq", &primvertex_ptq, "primvertex_pdf/F");
  tree->Branch("primvertex_ntracks", &primvertex_ntracks, "primvertex_ntracks/I");
  tree->Branch("primvertex_cov", primvertex_cov, "primvertex_cov[6]/F");
  tree->Branch("primvertex_mindz", &primvertex_mindz, "primvertex_mindz/F");


}


METScanningNtupleMaker::~METScanningNtupleMaker() { 
  
  tf1->cd();
  tree->Write();
  tf1->Write();
  tf1->Close();  
}


void METScanningNtupleMaker::beginRun(const edm::Run& run, 
				 const edm::EventSetup & es) { }

bool METScanningNtupleMaker::AddFlags(const edm::Event& iEvent, const char* module, const char* label, const char* process) {

  iEvent.getByLabel(edm::InputTag( module, label, process), Flags);
  if (!Flags.isValid())
    return false;

  const edm::TriggerNames& FlagNames_ = iEvent.triggerNames(*Flags);
  for(unsigned i = 0 ; i < Flags->size(); i++){
    if(!Flags->wasrun(i) )continue;
    std::string flagName=FlagNames_.triggerName(i);
    if(cFlags.size() > 0){
      for(size_t ip = 0; ip < cFlags.size(); ip++){
if(flagName.find(cFlags[ip]) != string::npos){

  flags_->insert(std::pair<string, int>(flagName, Flags->accept(i)));
  TString TriggerName(flagName);
  //std::cout << flagName << " : " << Flags->accept(i) << std::endl;
         }
      }
     }
    }

  return true;
  }


unsigned int METScanningNtupleMaker::AddMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //edm::ESHandle<TransientTrackBuilder> builder;
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",builder);
  //const TransientTrackBuilder * transientTrackBuilder = builder.product();

  edm::Handle<pat::MuonCollection> Muons;
  iEvent.getByToken(MuonCollectionToken_, Muons);

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken( PackedCantidateCollectionToken_, pfcands);
  /*
  edm::Handle<edm::PtrVector<reco::Muon>> BadDuplicateMuons;
  iEvent.getByToken(BadDuplicateMuonsToken_, BadDuplicateMuons);
  edm::Handle<edm::PtrVector<reco::Muon>> BadGlobalMuons;
  iEvent.getByToken(BadGlobalMuonsToken_, BadGlobalMuons);
  */

  if(Muons.isValid())
    {
      for(unsigned i = 0 ; i < Muons->size() ; i++){


        muon_px[muon_count] = (*Muons)[i].px();
        muon_py[muon_count] = (*Muons)[i].py();
        muon_pz[muon_count] = (*Muons)[i].pz();
        muon_pt[muon_count] = (*Muons)[i].pt();
        muon_eta[muon_count] = (*Muons)[i].eta();
        muon_phi[muon_count] = (*Muons)[i].phi();
        muon_charge[muon_count] = (*Muons)[i].charge();
        muon_vx[muon_count] = (*Muons)[i].vx(); // gives the same as (*Muons)[i].muonBestTrack()->vx()
        muon_vy[muon_count] = (*Muons)[i].vy();
        muon_vz[muon_count] = (*Muons)[i].vz();

        const pat::Muon &lep = (*Muons)[i];

        if((*Muons)[i].globalTrack().isNonnull())
          {
            muon_globalTrack[muon_count] = true;
            muon_pterror[muon_count] = (*Muons)[i].globalTrack()->ptError();
            muon_chi2[muon_count] = (*Muons)[i].globalTrack()->chi2();
            muon_ndof[muon_count] = (*Muons)[i].globalTrack()->ndof();
            muon_nMuonHits[muon_count] = (*Muons)[i].globalTrack()->hitPattern().numberOfValidMuonHits();
            muon_normChi2[muon_count] = (*Muons)[i].globalTrack()->normalizedChi2();
          }
        else
          {
            muon_globalTrack[muon_count] = false;
            muon_pterror[muon_count] = -1.;
            muon_chi2[muon_count] = -1.;
            muon_ndof[muon_count] = 0;
            muon_nMuonHits[muon_count] = 0;
            muon_normChi2[muon_count] = -1;
          }

        //	std::cout << "  chi2 = " << muon_chi2[muon_count] << "  ndof = " << muon_ndof[muon_count] << std::endl;

        muon_nMuonStations[muon_count] = (*Muons)[i].numberOfMatchedStations();

        muon_isTracker[muon_count] = (*Muons)[i].isTrackerMuon();
        muon_isPF[muon_count] = (*Muons)[i].isPFMuon();
        //muon_isTight[muon_count] = (*Muons)[i].isTightMuon(primvertex);
        //muon_isLoose[muon_count] = (*Muons)[i].isLooseMuon();
        //muon_isGlobal[muon_count] = (*Muons)[i].isGlobalMuon();
        //muon_isMedium[muon_count] = (*Muons)[i].isMediumMuon();
        //muon_isICHEP[muon_count] = idAlgos::isICHEPMuon((*Muons)[i]);

        muon_chargedHadIso[muon_count] = (*Muons)[i].chargedHadronIso();
        muon_neutralHadIso[muon_count] = (*Muons)[i].neutralHadronIso();
        muon_photonIso[muon_count] = (*Muons)[i].photonIso();
        muon_puIso[muon_count] = (*Muons)[i].puChargedHadronIso();

        muon_r03_sumChargedHadronPt[muon_count] = (*Muons)[i].pfIsolationR03().sumChargedHadronPt;
        muon_r03_sumChargedParticlePt[muon_count] = (*Muons)[i].pfIsolationR03().sumChargedParticlePt;
        muon_r03_sumNeutralHadronEt[muon_count] = (*Muons)[i].pfIsolationR03().sumNeutralHadronEt;
        muon_r03_sumPhotonEt[muon_count] = (*Muons)[i].pfIsolationR03().sumPhotonEt;
        muon_r03_sumNeutralHadronEtHighThreshold[muon_count] = (*Muons)[i].pfIsolationR03().sumNeutralHadronEtHighThreshold;
        muon_r03_sumPhotonEtHighThreshold[muon_count] = (*Muons)[i].pfIsolationR03().sumPhotonEtHighThreshold;
        muon_r03_sumPUPt[muon_count] = (*Muons)[i].pfIsolationR03().sumPUPt;

        muon_r04_sumChargedHadronPt[muon_count] = (*Muons)[i].pfIsolationR04().sumChargedHadronPt;
        muon_r04_sumChargedParticlePt[muon_count] = (*Muons)[i].pfIsolationR04().sumChargedParticlePt;
        muon_r04_sumNeutralHadronEt[muon_count] = (*Muons)[i].pfIsolationR04().sumNeutralHadronEt;
        muon_r04_sumPhotonEt[muon_count] = (*Muons)[i].pfIsolationR04().sumPhotonEt;
        muon_r04_sumNeutralHadronEtHighThreshold[muon_count] = (*Muons)[i].pfIsolationR04().sumNeutralHadronEtHighThreshold;
        muon_r04_sumPhotonEtHighThreshold[muon_count] = (*Muons)[i].pfIsolationR04().sumPhotonEtHighThreshold;
        muon_r04_sumPUPt[muon_count] = (*Muons)[i].pfIsolationR04().sumPUPt;

        TrackRef innertrack = (*Muons)[i].innerTrack();
        TrackRef bestTrack  = (*Muons)[i].muonBestTrack();

        muon_combQ_trkKink[muon_count] = (*Muons)[i].combinedQuality().trkKink;
        muon_combQ_chi2LocalPosition[muon_count] = (*Muons)[i].combinedQuality().chi2LocalPosition;
        //muon_segmentComp[muon_count] = (*Muons)[i].segmentCompatibility();

        if (bestTrack.isNonnull()) {
          muon_dxy[muon_count]    = bestTrack->dxy(pv_position);
          muon_dz[muon_count]     = bestTrack->dz(pv_position);
          muon_dxyerr[muon_count]    = bestTrack->dxyError();
          muon_dzerr[muon_count]     = bestTrack->dzError();
        }
        else {
          muon_dxy[muon_count]    = -9999;
          muon_dz[muon_count]     = -9999;
          muon_dxyerr[muon_count]    = -9999;
          muon_dzerr[muon_count]     = -9999;
        }

        if(innertrack.isNonnull())
          {
            muon_innerTrack[muon_count] = true;
            muon_nPixelHits[muon_count] = innertrack->hitPattern().numberOfValidPixelHits();
            muon_nTrackerHits[muon_count] = innertrack->hitPattern().trackerLayersWithMeasurement();
            muon_validFraction[muon_count] = innertrack->validFraction();
          }
        else
          {
            muon_innerTrack[muon_count] = false;
            muon_nPixelHits[muon_count] = 0;
            muon_nTrackerHits[muon_count] = 0;
            muon_validFraction[muon_count] = 0;
          }


        // Dimuons
        if( !(*Muons)[i].innerTrack().isNull()){
          for(unsigned j = i+1 ; j < Muons->size() ; j++){
            if (dimuon_count==M_muonmaxcount*(M_muonmaxcount-1)/2) {
              cerr << "number of dimuons > M_muonmaxcount*(M_muonmaxcount-1)/2. They are missing." << endl; errors |= 1<<1;
              break;
            }

            //if ((*Muons)[j].pt() < cMuPtMin) continue;
            //if (fabs(((*Muons)[j].eta()))>cMuEtaMax) continue;
            if( (*Muons)[j].innerTrack().isNull()) continue;

            dimuon_leading[dimuon_count] = i;
            dimuon_trailing[dimuon_count] = j;
            if ((*Muons)[i].pt() < (*Muons)[j].pt()){
              dimuon_leading[dimuon_count] = j;
              dimuon_trailing[dimuon_count] = i;
            }

            if (fabs( (*Muons)[i].pt() - (*Muons)[j].pt()) < 1.e-4){
              std::cout<<"WTF!!"<<std::endl;
            }

            dimuon_dist2D[dimuon_count] = -1.;
            dimuon_dist2DE[dimuon_count] = -1.;
            dimuon_dist3D[dimuon_count] = -1.;
            dimuon_dist3DE[dimuon_count] = -1.;



            dimuon_count++;
          }
        }

        muon_count++;

      }
    }
  return muon_count;
}


void METScanningNtupleMaker::analyze(const Event& iEvent,
				const EventSetup& iSetup) {
  
   //Counters:
   track_count = 0;
   goodprimvertex_count = 0;
   primvertex_count = 0;
   muon_count = 0;
   dimuon_count = 0;
   tau_count = 0;
   l1muon_count = 0;
   l1egamma_count = 0;
   l1tau_count = 0;
   l1isotau_count = 0;
   gentau_count = 0;
   pfjet_count = 0;
   electron_count = 0;
   photon_count = 0;
   genparticles_count = 0;
   genjets_count = 0;
   errors = 0;
   trigobject_count = 0;
   mvamet_count = 0;


  edm::ESHandle<EcalSeverityLevelAlgo> sevlv;
  iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv);

  //Run, Event, LumiBlock, Time:
  irun  = iEvent.id().run();
  ievent  = iEvent.id().event();
  ilumiBlock = iEvent.id().luminosityBlock();
  itime = iEvent.time();

  run = (size_t)irun;
  event = (size_t)ievent;
  lumiBlock = (size_t)ilumiBlock;
  time = (size_t)((iEvent.time().value())>>32);

  //Triggers:
  flags_->clear();
  for(std::vector<string>::iterator it = cFlagsProcesses.begin(); it != cFlagsProcesses.end(); it++){

    AddFlags(iEvent,"TriggerResults", "", it->data());
  }

   //Puppi MET
  if(crecpuppimet==true){
  edm::Handle<pat::METCollection> patMet;
  iEvent.getByToken(PuppiMetCollectionToken_, patMet);

  assert(patMet->size() > 0);
  puppimet_ex = (*patMet)[0].px();
  puppimet_ey = (*patMet)[0].py();
  puppimet_ez = (*patMet)[0].pz();
  puppimet_pt = (*patMet)[0].pt();
  puppimet_phi = (*patMet)[0].phi();
  puppimet_sigxx = (*patMet)[0].getSignificanceMatrix()(0,0);
  puppimet_sigxy = (*patMet)[0].getSignificanceMatrix()(0,1);
  puppimet_sigyx = (*patMet)[0].getSignificanceMatrix()(1,0);
  puppimet_sigyy = (*patMet)[0].getSignificanceMatrix()(1,1);

  }

   //PF MET
  //if(crecpfmet==true)
 // {
      edm::Handle<pat::METCollection> patMet;
      iEvent.getByToken(MetCollectionToken_, patMet);

      assert(patMet->size() > 0);
      pfmet_ex = (*patMet)[0].px();
      pfmet_ey = (*patMet)[0].py();
      pfmet_ez = (*patMet)[0].pz();
      pfmet_pt = (*patMet)[0].pt();
      pfmet_phi = (*patMet)[0].phi();

      pfmet_sigxx = (*patMet)[0].getSignificanceMatrix()(0,0);
      pfmet_sigxy = (*patMet)[0].getSignificanceMatrix()(0,1);
      pfmet_sigyx = (*patMet)[0].getSignificanceMatrix()(1,0);
      pfmet_sigyy = (*patMet)[0].getSignificanceMatrix()(1,1);
      pfmet_sig   = (*patMet)[0].significance();



  //N Vertex:
  if(crecprimvertex==true)
    {
      edm::Handle<VertexCollection> Vertex;
      iEvent.getByToken(PVToken_, Vertex);
      if(Vertex.isValid()) {
        primvertex_mindz = 999;
        for(unsigned i = 0 ; i < Vertex->size(); i++) {
          primvertex_count++;
          if(i == 0) {
            primvertex_x = (*Vertex)[i].x();
            primvertex_y = (*Vertex)[i].y();
            primvertex_z = (*Vertex)[i].z();
            primvertex_chi2 = (*Vertex)[i].chi2();
            primvertex_ndof = (*Vertex)[i].ndof();
            primvertex_ntracks = (*Vertex)[i].tracksSize();
            primvertex_cov[0] = (*Vertex)[i].covariance(0,0); // xError()
            primvertex_cov[1] = (*Vertex)[i].covariance(0,1);
            primvertex_cov[2] = (*Vertex)[i].covariance(0,2);
            primvertex_cov[3] = (*Vertex)[i].covariance(1,1); // yError()
            primvertex_cov[4] = (*Vertex)[i].covariance(1,2);
            primvertex_cov[5] = (*Vertex)[i].covariance(2,2); // zError()
            Float_t ptq = 0.;
            for(Vertex::trackRef_iterator it = (*Vertex)[i].tracks_begin() ; it != (*Vertex)[i].tracks_end() ; ++it)
              {
                ptq += (*it)->pt() * (*it)->pt();
              }
            primvertex_ptq = ptq;

            pv_position = (*Vertex)[i].position();
            primvertex = (*Vertex)[i];
          } else {
            if(std::abs((*Vertex)[i].z()-(*Vertex)[0].z()) < primvertex_mindz)
              primvertex_mindz = std::abs((*Vertex)[i].z()-(*Vertex)[0].z()); //minimal longitudinal distance between the PV and other vertex
          }
          if((*Vertex)[i].isValid() && !(*Vertex)[i].isFake() && (*Vertex)[i].ndof() >= 4 && (*Vertex)[i].z() > -24 && (*Vertex)[i].z() < 24 && (*Vertex)[i].position().Rho() < 2.)
            goodprimvertex_count++;
        }
      }
    }


  //Muons:
  if(crecmuon==true){
    int numberOfMuons = int(AddMuons(iEvent, iSetup));
  }


  //Tracks:
  //TODO Insert the track collections;


  //tree filling ===========================
  tree->Fill();

}



DEFINE_FWK_MODULE(METScanningNtupleMaker);






