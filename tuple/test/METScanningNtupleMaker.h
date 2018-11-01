#ifndef RecoParticleFlow_PFPatProducer_METScan_
#define RecoParticleFlow_PFPatProducer_METScan_

// system include files
#include <memory>
#include <string>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFClusterMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFClusterMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"

//ld Trigger
//#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

//ld
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/PatCandidates/interface/PATTauDiscriminator.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"




#define M_muonmaxcount 1000
#define M_trackmaxcount 1000

#include <math.h>

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TVector3.h>

using namespace std;
using namespace edm;
using namespace reco;


class METScanningNtupleMaker : public edm::EDAnalyzer {
 public:

  explicit METScanningNtupleMaker(const edm::ParameterSet&);

  ~METScanningNtupleMaker();
   virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  unsigned int AddMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  bool AddFlags(const edm::Event& iEvent, const char* module, const char* label, const char* process);
  virtual void beginRun(const edm::Run & r, const edm::EventSetup & c);

  //Configuration
   bool cdata;
   bool crecmuon;
   bool crecpuppimet;
   bool crecpfmet;
   bool crecprimvertex;

   vector<string> cFlags;
   vector<string> cFlagsProcesses;

     edm::EDGetTokenT<bool> BadChCandFilterToken_;
   edm::EDGetTokenT<bool> BadPFMuonFilterToken_;
   edm::EDGetTokenT<pat::METCollection> PuppiMetCollectionToken_;
   edm::EDGetTokenT<pat::METCollection> MetCollectionToken_;
   vector<string> cHLTriggerPaths;
   string cTriggerProcess;
   edm::EDGetTokenT<pat::MuonCollection> MuonCollectionToken_;


   edm::EDGetTokenT<VertexCollection> PVToken_;

   //Useless cuts inherited from MINIAOD Code:
   double cMuPtMin = 10000.0;
   double cMuEtaMax = 10000.0;

   double getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
                           const reco::Candidate* ptcl,
                           double r_iso_min, double r_iso_max, double kt_scale,
                           bool charged_only) {

       if (ptcl->pt()<5.) return 99999.;

       double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
       if(ptcl->isElectron()) {
         if (fabs(ptcl->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
       } else if(ptcl->isMuon()) {
         deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;
       } else {
         //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
       }

       double iso_nh(0.); double iso_ch(0.);
       double iso_ph(0.); double iso_pu(0.);
       double ptThresh(0.5);
       if(ptcl->isElectron()) ptThresh = 0;
       double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/ptcl->pt()));
       for (const pat::PackedCandidate &pfc : *pfcands) {
         if (abs(pfc.pdgId())<7) continue;

         double dr = deltaR(pfc, *ptcl);
         if (dr > r_iso) continue;

         //////////////////  NEUTRALS  /////////////////////////
         if (pfc.charge()==0){
           if (pfc.pt()>ptThresh) {
             /////////// PHOTONS ////////////
             if (abs(pfc.pdgId())==22) {
               if(dr < deadcone_ph) continue;
               iso_ph += pfc.pt();
               /////////// NEUTRAL HADRONS ////////////
             } else if (abs(pfc.pdgId())==130) {
               if(dr < deadcone_nh) continue;
               iso_nh += pfc.pt();
             }
           }
           //////////////////  CHARGED from PV  /////////////////////////
         } else if (pfc.fromPV()>1){
           if (abs(pfc.pdgId())==211) {
             if(dr < deadcone_ch) continue;
             iso_ch += pfc.pt();
           }
           //////////////////  CHARGED from PU  /////////////////////////
         } else {
           if (pfc.pt()>ptThresh){
             if(dr < deadcone_pu) continue;
             iso_pu += pfc.pt();
           }
         }
       }
       double iso(0.);
       if (charged_only){
         iso = iso_ch;
       } else {
         iso = iso_ph + iso_nh;
         iso -= 0.5*iso_pu;
         if (iso>0) iso += iso_ch;
         else iso = iso_ch;
       }
       iso = iso/ptcl->pt();

       return iso;
     }


   string cMyTriggerProcess;



   edm::EDGetTokenT<reco::GenParticleCollection> GenParticleCollectionToken_;
   edm::EDGetTokenT<pat::PackedCandidateCollection> PackedCantidateCollectionToken_;
   edm::ESHandle<TransientTrackBuilder>  TTrackBuilder;

   // pat muons
    UInt_t muon_count;
    Float_t muon_px[M_muonmaxcount];
    Float_t muon_py[M_muonmaxcount];
    Float_t muon_pz[M_muonmaxcount];
    Float_t muon_pt[M_muonmaxcount];
    Float_t muon_eta[M_muonmaxcount];
    Float_t muon_phi[M_muonmaxcount];
    Float_t muon_pterror[M_muonmaxcount];
    Float_t muon_chi2[M_muonmaxcount];
    Float_t muon_normChi2[M_muonmaxcount];
    Float_t muon_ndof[M_muonmaxcount];
    Float_t muon_charge[M_muonmaxcount];
    Float_t muon_miniISO[M_muonmaxcount];
    Bool_t muon_isPF[M_muonmaxcount];
    Bool_t muon_isGlobal[M_muonmaxcount];
    Bool_t muon_isTracker[M_muonmaxcount];
    Bool_t muon_isTight[M_muonmaxcount];
    Bool_t muon_isLoose[M_muonmaxcount];
    Bool_t muon_isMedium[M_muonmaxcount];
    Bool_t muon_isICHEP[M_muonmaxcount];
    Bool_t muon_isDuplicate[M_muonmaxcount];
    Bool_t muon_isBad[M_muonmaxcount];
    Bool_t muon_globalTrack[M_muonmaxcount];
    Bool_t muon_innerTrack[M_muonmaxcount];
    Int_t muon_genmatch[M_muonmaxcount];
    UInt_t muon_nMuonStations[M_muonmaxcount];
    UInt_t muon_nMuonHits[M_muonmaxcount];
    UInt_t muon_nPixelHits[M_muonmaxcount];
    UInt_t muon_nTrackerHits[M_muonmaxcount];
    Float_t muon_chargedHadIso[M_muonmaxcount];
    Float_t muon_neutralHadIso[M_muonmaxcount];
    Float_t muon_photonIso[M_muonmaxcount];
    Float_t muon_puIso[M_muonmaxcount];
    Float_t muon_dz[M_muonmaxcount];
    Float_t muon_dzerr[M_muonmaxcount];
    Float_t muon_dxy[M_muonmaxcount];
    Float_t muon_dxyerr[M_muonmaxcount];
    Float_t muon_vx[M_muonmaxcount];
    Float_t muon_vy[M_muonmaxcount];
    Float_t muon_vz[M_muonmaxcount];

    Float_t muon_r03_sumChargedHadronPt[M_muonmaxcount];
    Float_t muon_r03_sumChargedParticlePt[M_muonmaxcount];
    Float_t muon_r03_sumNeutralHadronEt[M_muonmaxcount];
    Float_t muon_r03_sumPhotonEt[M_muonmaxcount];
    Float_t muon_r03_sumNeutralHadronEtHighThreshold[M_muonmaxcount];
    Float_t muon_r03_sumPhotonEtHighThreshold[M_muonmaxcount];
    Float_t muon_r03_sumPUPt[M_muonmaxcount];

    Float_t muon_combQ_chi2LocalPosition[M_muonmaxcount];
     Float_t muon_combQ_trkKink[M_muonmaxcount];
     Float_t muon_validFraction[M_muonmaxcount];
     Float_t muon_segmentComp[M_muonmaxcount];

     Float_t muon_r04_sumChargedHadronPt[M_muonmaxcount];
     Float_t muon_r04_sumChargedParticlePt[M_muonmaxcount];
     Float_t muon_r04_sumNeutralHadronEt[M_muonmaxcount];
     Float_t muon_r04_sumPhotonEt[M_muonmaxcount];
     Float_t muon_r04_sumNeutralHadronEtHighThreshold[M_muonmaxcount];
     Float_t muon_r04_sumPhotonEtHighThreshold[M_muonmaxcount];
     Float_t muon_r04_sumPUPt[M_muonmaxcount];

     //tracks
     UInt_t errors;
     UInt_t track_count;
     UInt_t photon_count;
     UInt_t genparticles_count;
     UInt_t genjets_count;
     UInt_t trigobject_count;
     UInt_t mvamet_count;
     UInt_t electron_count;
     UInt_t gentau_count;
     UInt_t l1isotau_count;
     UInt_t tau_count;
     UInt_t  l1muon_count;
     UInt_t  l1egamma_count;
     UInt_t l1tau_count;
     UInt_t pfjet_count;

     UInt_t dimuon_count;
     UInt_t dimuon_leading[M_muonmaxcount*(M_muonmaxcount - 1)/2];
     UInt_t dimuon_trailing[M_muonmaxcount*(M_muonmaxcount - 1)/2];
     Float_t dimuon_dist2D[M_muonmaxcount*(M_muonmaxcount - 1)/2];
     Float_t dimuon_dist2DE[M_muonmaxcount*(M_muonmaxcount - 1)/2];
     Float_t dimuon_dist3D[M_muonmaxcount*(M_muonmaxcount - 1)/2];
     Float_t dimuon_dist3DE[M_muonmaxcount*(M_muonmaxcount - 1)/2];


     Float_t track_px[M_trackmaxcount];
     Float_t track_py[M_trackmaxcount];
     Float_t track_pz[M_trackmaxcount];
     Float_t track_pt[M_trackmaxcount];
     Float_t track_eta[M_trackmaxcount];
     Float_t track_phi[M_trackmaxcount];
     Float_t track_mass[M_trackmaxcount];
     Float_t track_charge[M_trackmaxcount];
     Float_t track_outerx[M_trackmaxcount];
     Float_t track_outery[M_trackmaxcount];
     Float_t track_outerz[M_trackmaxcount];
     Float_t track_closestpointx[M_trackmaxcount];
     Float_t track_closestpointy[M_trackmaxcount];
     Float_t track_closestpointz[M_trackmaxcount];
     Float_t track_chi2[M_trackmaxcount];
     Float_t track_ndof[M_trackmaxcount];
     Float_t track_dxy[M_trackmaxcount];
     Float_t track_dxyerr[M_trackmaxcount];
     Float_t track_dz[M_trackmaxcount];
     Float_t track_dzerr[M_trackmaxcount];
     Float_t track_dedxharmonic2[M_trackmaxcount];
     UChar_t track_nhits[M_trackmaxcount];
     UChar_t track_nmissinghits[M_trackmaxcount];
     UChar_t track_npixelhits[M_trackmaxcount];
     UChar_t track_npixellayers[M_trackmaxcount];
     UChar_t track_nstriplayers[M_trackmaxcount];
     Int_t track_ID[M_trackmaxcount];
     Bool_t track_highPurity[M_trackmaxcount];
     Float_t track_vx[M_trackmaxcount];
     Float_t track_vy[M_trackmaxcount];
     Float_t track_vz[M_trackmaxcount];


    // primary vertex
    UInt_t  primvertex_count;
    UInt_t  goodprimvertex_count;
    Float_t primvertex_x;
    Float_t primvertex_y;
    Float_t primvertex_z;
    Float_t primvertex_chi2;
    Float_t primvertex_ndof;
    Float_t primvertex_ptq;
    Int_t   primvertex_ntracks;
    Float_t primvertex_cov[6];
    Float_t primvertex_mindz;

    // met
    Float_t pfmet_ex;
    Float_t pfmet_ey;
    Float_t pfmet_ez;
    Float_t pfmet_pt;
    Float_t pfmet_phi;
    Float_t pfmet_sumet;

    Float_t pfmet_sig;
    Float_t pfmet_sigxx;
    Float_t pfmet_sigxy;
    Float_t pfmet_sigyx;
    Float_t pfmet_sigyy;

    Float_t pfmet_ex_JetEnUp;
    Float_t pfmet_ey_JetEnUp;

    Float_t pfmet_ex_JetEnDown;
    Float_t pfmet_ey_JetEnDown;

    Float_t pfmet_ex_UnclusteredEnUp;
    Float_t pfmet_ey_UnclusteredEnUp;

    Float_t pfmet_ex_UnclusteredEnDown;
    Float_t pfmet_ey_UnclusteredEnDown;

    //Puppi Met
    Float_t puppimet_ex;
    Float_t puppimet_ey;
    Float_t puppimet_ez;
    Float_t puppimet_pt;
    Float_t puppimet_phi;
    Float_t puppimet_sumet;

    Float_t puppimet_ex_JetEnUp;
    Float_t puppimet_ey_JetEnUp;

    Float_t puppimet_ex_JetEnDown;
    Float_t puppimet_ey_JetEnDown;

    Float_t puppimet_ex_UnclusteredEnUp;
    Float_t puppimet_ey_UnclusteredEnUp;

    Float_t puppimet_ex_UnclusteredEnDown;
    Float_t puppimet_ey_UnclusteredEnDown;

    Float_t puppimet_ex_JetResUp;
    Float_t puppimet_ey_JetResUp;

    Float_t puppimet_ex_JetResDown;
    Float_t puppimet_ey_JetResDown;

    Float_t puppimet_sigxx;
    Float_t puppimet_sigxy;
    Float_t puppimet_sigyx;
    Float_t puppimet_sigyy;

   //Filters Flags:
   edm::Handle<edm::TriggerResults> Flags;
   std::map<std::string, int>* flags_;

    //Vertex Token
    edm::EDGetTokenT<vector<reco::Vertex> >  vertex_token;


    bool isReco_token;
    bool cgen;
    bool csusyinfo;
    bool ctrigger;
    bool cbeamspot;
    bool crectrack;
    bool crecelectron;
    bool crectau;
    bool cl1isotau;
    bool cl1objects;
    bool crecphoton;
    bool crecpfjet;
    bool crecpfmetcorr;
    bool crecmvamet;


  size_t run,event,lumiBlock,time;
  edm::RunNumber_t irun;
  edm::EventNumber_t ievent;
  edm::LuminosityBlockNumber_t ilumiBlock;
  edm::Timestamp itime;

  size_t nVtx;
  math::XYZPoint pv_position;
  Vertex primvertex;

  Float_t genmet_ex;
  Float_t genmet_ey;

  //tree stuff
  std::string outputfile_;
  TFile* tf1;
  TTree* tree;

 private:




};

#endif
