// -*- C++ -*-
//
// Package:    VNAnalyzer
// Class:      VNAnalyzer
// 


// system include files
#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/Vector3D.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/DataRecord/interface/HeavyIonRPRcd.h"
#include "CondFormats/DataRecord/interface/HeavyIonRcd.h"
#include "CondFormats/HIObjects/interface/CentralityTable.h"
#include "CondFormats/HIObjects/interface/RPFlatParams.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "HeavyIonsAnalysis/VNAnalysis/interface/TrackEfficiency.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TH1I.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include <time.h>
#include <cstdlib>
	
#include <vector>
using std::vector;
using std::rand;
using namespace std;
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/LoadEPDB.h"
using namespace hi;
using namespace edm;

static const int ntrkbins = 25;
static const  double trkBins[]={0, 10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 135, 150, 160, 185, 210, 230, 250, 270, 300, 330, 350, 370, 390, 420, 500};


static const int ncentbins = 14;
static const  double centbins[]={0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100};

static const int nptbins = 28;
static const float ptbins[]={0.3, 0.4, 0.5,  0.6,  0.8,  1.0,  1.25,  1.50,  2.0,
			      2.5,  3.0,  3.5,  4.0,  5.0,  6.0,  7.0, 8.0, 10., 12.0, 14.0, 16.0,  20.0, 26.0, 35.0, 45.0, 60.0, 80.0, 100., 200.};

static const int MaxTracks = 50;

static const int netabinsDefault = 12;
static const float etabinsDefault[]={-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};


//
// class declaration
//

class VNAnalyzer : public edm::EDAnalyzer {
public:
  explicit VNAnalyzer(const edm::ParameterSet&);
  ~VNAnalyzer();
      
private:
  int NtrkToBin(int ntrk){
    for(int i = 0; i<=ntrkbins; i++) {
      if(ntrk < trkBins[i]) return i;
    }
    return ntrkbins;
  }
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  bool CaloMatch(const reco::Track & track, const edm::Event & iEvent, unsigned int idx);
  // ----------member data ---------------------------
  int eporder_;


  std::string centralityVariable_;
  std::string centralityLabel_;
  std::string centralityMC_;

  edm::InputTag centralityBinTag_;
  edm::EDGetTokenT<int> centralityBinToken;
  edm::Handle<int> cbin_;
  edm::EDGetTokenT<int> tag_;

  edm::InputTag centralityTag_;
  edm::EDGetTokenT<reco::Centrality> centralityToken;
  edm::Handle<reco::Centrality> centrality_;


  edm::InputTag pfTag_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfToken_;


  edm::InputTag vertexTag_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken;
  edm::Handle<std::vector<reco::Vertex>> vertex_;

  edm::InputTag trackTag_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken;
  edm::Handle<reco::TrackCollection> trackCollection_;

  edm::InputTag inputPlanesTag_;
  edm::EDGetTokenT<reco::EvtPlaneCollection> inputPlanesToken;
  edm::Handle<reco::EvtPlaneCollection> inputPlanes_;

  edm::Service<TFileService> fs;
  TFile * frecenter;
  string offsetFileName;

  double caloCentRef_;
  double caloCentRefWidth_;
  int caloCentRefMinBin_;
  int caloCentRefMaxBin_;

  double nCentBins_;
  bool useNtrk_;

  int vs_sell;   // vertex collection size
  float vzr_sell;
  float vzErr_sell;
  TH1D * hcent;
  TH1D * hcentbins;
  TH1D * hcentres;
  TH1D * hptNtrk;
  TH1D * hptNtrkGood;
  TH1I * hNtrkRet;
  TH2D * hEff[ntrkbins];
  double centval;
  int ntrkval;
  double vtx;
  int Noff;
  double reso_;
  bool bCaloMatching_;
  int nvtx_;
  double minvz_;
  double maxvz_;
  double dzerr_;
  double chi2_;

  Double_t epang[NumEPNames];
  Double_t eporig[NumEPNames];
  Double_t epsin[NumEPNames];
  Double_t epcos[NumEPNames];

  Double_t qx[NumEPNames];
  Double_t qy[NumEPNames];
  Double_t q[NumEPNames];
  Double_t epmult[NumEPNames];
  Double_t sumw[NumEPNames];
  Double_t sumw2[NumEPNames];
  Double_t vn[NumEPNames];

  Double_t rescor[NumEPNames];
  Double_t rescorErr[NumEPNames];
  TH1D * hPsi[NumEPNames];
  TH1D * hPsiOffset[NumEPNames];
  TH1D * hPsiFlat[NumEPNames];


  unsigned int runno_;

  TH1D * hNtrkoff;
  int nEtaBins;
  TH1I * hrun;
  string rpnames[NumEPNames];
  string effTable_;
  TTree * tree;
  TrackEfficiency *teff;
  int FlatOrder_;
  int NumFlatBins_;
  int CentBinCompression_;
  int Noffmin_;
  int Noffmax_;
  int EPOrder_;
  int EPLevel_;
  TH2D * qxtrk1;
  TH2D * qytrk1;
  TH2D * qxtrk2;
  TH2D * qytrk2;
  TH2D * qxtrk3;
  TH2D * qytrk3;
  TH2D * qxtrk4;
  TH2D * qytrk4;
  TH2D * qxtrk5;
  TH2D * qytrk5;
  TH2D * qxtrk6;
  TH2D * qytrk6;
  TH2D * qxtrk7;
  TH2D * qytrk7;
  TH2D * qcnt;
  TH2D * avpt;
  HiEvtPlaneFlatten * flat[NumEPNames];
  bool loadDB_;
  bool useNtrkBins_; 
  bool bypassCentrality_;
  bool FirstEvent_;
  bool MB_;
  bool Recenter_;
  int minrun_;
  int maxrun_;
  TH2D * wqxtrkRef[7][40];
  TH2D * wqytrkRef[7][40];


  int ntrack;
  float sppt[MaxTracks];
  float spphi[MaxTracks];
  float speta[MaxTracks];

  enum    TrackCut {trackUndefine = 0, ppReco = 1, HIReco, Pixel};
  TrackCut sTrackQuality;
  double  dzdzerror_;
  double  d0d0error_;
  double  pterrorpt_;
  double  dzdzerror_pix_;
  bool TrackQuality_ppReco(const reco::TrackCollection::const_iterator&, const reco::VertexCollection&);
  bool TrackQuality_HIReco(const reco::TrackCollection::const_iterator&, const reco::VertexCollection&);
  bool TrackQuality_Pixel(const reco::TrackCollection::const_iterator&, const reco::VertexCollection&);

  int getNoff(const edm::Event& iEvent, const edm::EventSetup& iSetup, int bin)
  {
    int Noff = 0;
    using namespace edm;
    using namespace reco;
    qxtrk1->Reset();
    qytrk1->Reset();
    qxtrk2->Reset();
    qytrk2->Reset();
    qxtrk3->Reset();
    qytrk3->Reset();
    qxtrk4->Reset();
    qytrk4->Reset();
    qxtrk5->Reset();
    qytrk5->Reset();
    qxtrk6->Reset();
    qytrk6->Reset();
    qxtrk7->Reset();
    qytrk7->Reset();
    qcnt->Reset();
    avpt->Reset();
    iEvent.getByToken(vertexToken,vertex_);
    VertexCollection recoVertices = *vertex_;
    if ( recoVertices.size() > 100 ) return -1;
    sort(recoVertices.begin(), recoVertices.end(), [](const reco::Vertex &a, const reco::Vertex &b){
	if ( a.tracksSize() == b.tracksSize() ) return a.chi2() < b.chi2();
	return a.tracksSize() > b.tracksSize();
      });
   
    int primaryvtx = 0;
   
    double vz = recoVertices[primaryvtx].z();
    if (fabs(vz) < -15 || fabs(vz) > 15) {
      return -1;
    }
 
    math::XYZPoint v1( recoVertices[primaryvtx].position().x(), recoVertices[primaryvtx].position().y(), recoVertices[primaryvtx].position().z() );
    

    iEvent.getByLabel(trackTag_,trackCollection_);

    for(TrackCollection::const_iterator itTrack = trackCollection_->begin(); itTrack != trackCollection_->end(); ++itTrack) {    

      if(fabs(itTrack->eta())<1) hptNtrk->Fill(itTrack->pt());
 
      if ( sTrackQuality == HIReco and not TrackQuality_HIReco(itTrack, recoVertices) ) continue;
      else if ( sTrackQuality == ppReco and not TrackQuality_ppReco(itTrack, recoVertices) ) continue;
      else if ( sTrackQuality == Pixel  and not TrackQuality_Pixel (itTrack, recoVertices) ) continue;


      if(fabs(itTrack->eta())<1) hptNtrkGood->Fill(itTrack->pt());
      if(itTrack->pt()>0.4) ++Noff;
      int ipt = qxtrk2->GetXaxis()->FindBin(itTrack->pt());
      int ieta = qxtrk2->GetYaxis()->FindBin(itTrack->eta());
      qxtrk1->Fill(itTrack->pt(), itTrack->eta(), TMath::Cos(itTrack->phi()) - wqxtrkRef[0][bin]->GetBinContent(ipt,ieta));
      qytrk1->Fill(itTrack->pt(), itTrack->eta(), TMath::Sin(itTrack->phi()) - wqytrkRef[0][bin]->GetBinContent(ipt,ieta));
      qxtrk2->Fill(itTrack->pt(), itTrack->eta(), TMath::Cos(2.*itTrack->phi()) - wqxtrkRef[1][bin]->GetBinContent(ipt,ieta));
      qytrk2->Fill(itTrack->pt(), itTrack->eta(), TMath::Sin(2.*itTrack->phi()) - wqytrkRef[1][bin]->GetBinContent(ipt,ieta));
      qxtrk3->Fill(itTrack->pt(), itTrack->eta(), TMath::Cos(3.*itTrack->phi()) - wqxtrkRef[2][bin]->GetBinContent(ipt,ieta));
      qytrk3->Fill(itTrack->pt(), itTrack->eta(), TMath::Sin(3.*itTrack->phi()) - wqytrkRef[2][bin]->GetBinContent(ipt,ieta));
      qxtrk4->Fill(itTrack->pt(), itTrack->eta(), TMath::Cos(4.*itTrack->phi()) - wqxtrkRef[3][bin]->GetBinContent(ipt,ieta));
      qytrk4->Fill(itTrack->pt(), itTrack->eta(), TMath::Sin(4.*itTrack->phi()) - wqytrkRef[3][bin]->GetBinContent(ipt,ieta));
      qxtrk5->Fill(itTrack->pt(), itTrack->eta(), TMath::Cos(5.*itTrack->phi()) - wqxtrkRef[4][bin]->GetBinContent(ipt,ieta));
      qytrk5->Fill(itTrack->pt(), itTrack->eta(), TMath::Sin(5.*itTrack->phi()) - wqytrkRef[4][bin]->GetBinContent(ipt,ieta));
      qxtrk6->Fill(itTrack->pt(), itTrack->eta(), TMath::Cos(6.*itTrack->phi()) - wqxtrkRef[5][bin]->GetBinContent(ipt,ieta));
      qytrk6->Fill(itTrack->pt(), itTrack->eta(), TMath::Sin(6.*itTrack->phi()) - wqytrkRef[5][bin]->GetBinContent(ipt,ieta));
      qxtrk7->Fill(itTrack->pt(), itTrack->eta(), TMath::Cos(7.*itTrack->phi()) - wqxtrkRef[6][bin]->GetBinContent(ipt,ieta));
      qytrk7->Fill(itTrack->pt(), itTrack->eta(), TMath::Sin(7.*itTrack->phi()) - wqytrkRef[6][bin]->GetBinContent(ipt,ieta));

      qcnt->Fill(itTrack->pt(), itTrack->eta());
      avpt->Fill(itTrack->pt(), itTrack->eta(), itTrack->pt());
      
      if( itTrack->pt() < 0.2 ) continue;
      hEff[bin]->Fill(itTrack->phi(),itTrack->eta());
    }
    return Noff;
  }


};


//
// constructors and destructor
//
VNAnalyzer::VNAnalyzer(const edm::ParameterSet& iConfig):runno_(0)
  
{
  runno_ = 0;
  loadDB_ = kTRUE;
  FirstEvent_ = kTRUE;
  for(int i = 0; i<NumEPNames; i++) {
    epang[i] = -10;
    eporig[i] = -10;
    epsin[i] = 0;
    epcos[i] = 0;
    qx[i] = 0;
    qy[i] = 0;
    q[i] = 0;
    epmult[i] = 0;
    rescor[i] = 0;
    rescorErr[i] = 0;
  }
  EPOrder_ = iConfig.getUntrackedParameter<int>("EPOrder_",2);
  FlatOrder_ = iConfig.getUntrackedParameter<int>("FlatOrder_", 9);
  NumFlatBins_ = iConfig.getUntrackedParameter<int>("NumFlatBins_",20);
  caloCentRef_ = iConfig.getUntrackedParameter<double>("caloCentRef_",80.);
  caloCentRefWidth_ = iConfig.getUntrackedParameter<double>("caloCentRefWidth_",5.);
  CentBinCompression_ = iConfig.getUntrackedParameter<int>("CentBinCompression_",5);
  Noffmin_ = iConfig.getUntrackedParameter<int>("Noffmin_", 0);
  Noffmax_ = iConfig.getUntrackedParameter<int>("Noffmax_", 50000);	
  minrun_ = iConfig.getUntrackedParameter<int>("minrun_", 0);
  maxrun_ = iConfig.getUntrackedParameter<int>("maxrun_", 50000);	
  effTable_ = iConfig.getParameter<std::string>("effTable_");

  centralityVariable_ = iConfig.getParameter<std::string>("centralityVariable");
  if(iConfig.exists("nonDefaultGlauberModel")){
    centralityMC_ = iConfig.getParameter<std::string>("nonDefaultGlauberModel");
  }
  centralityLabel_ = centralityVariable_+centralityMC_;

  centralityBinTag_ = iConfig.getParameter<edm::InputTag>("centralityBinTag_");
  centralityBinToken = consumes<int>(centralityBinTag_);

  centralityTag_ = iConfig.getParameter<edm::InputTag>("centralityTag_");
  centralityToken = consumes<reco::Centrality>(centralityTag_);
  if(centralityToken.isUninitialized()) {
    std::cout<<"centralityToken is uninitialized."<<std::endl;
  }
  vertexTag_  = iConfig.getParameter<edm::InputTag>("vertexTag_");
  vertexToken = consumes<std::vector<reco::Vertex>>(vertexTag_);
  if(vertexToken.isUninitialized()) {
    std::cout<<"vertexToken is uninitialized."<<std::endl;
  }

  trackTag_ = iConfig.getParameter<edm::InputTag>("trackTag_");
  trackToken = consumes<reco::TrackCollection>(trackTag_);
  if(trackToken.isUninitialized()) {
    std::cout<<"trackToken is uninitialized."<<std::endl;
  }
  if ( trackTag_.label() == "hiGeneralTracks" ) {
    sTrackQuality = HIReco;
    cout<<"hiGeneralTracks"<<endl;
  } else if ( trackTag_.label() == "generalTracks" ) {
    sTrackQuality = ppReco;
    cout<<"generalTracks"<<endl;
  } else if ( trackTag_.label() == "hiGeneralAndPixelTracks" ) {
    sTrackQuality = Pixel;
    cout<<"hiGeneralAndPixelTracks"<<endl;
  } else {
    sTrackQuality = trackUndefine;
  }
  useNtrk_ = iConfig.getUntrackedParameter<bool>("useNtrk",false);
  if(useNtrk_) {
    NumFlatBins_ = ntrkbins;
    CentBinCompression_ = 1;
  }

  inputPlanesTag_ = iConfig.getParameter<edm::InputTag>("inputPlanesTag_");
  inputPlanesToken = consumes<reco::EvtPlaneCollection>(inputPlanesTag_);
  if(inputPlanesToken.isUninitialized()) {
    std::cout<<"inputPlanesToken is uninitialized."<<std::endl;
  }
  tag_ = consumes<int>(iConfig.getParameter<edm::InputTag>("BinLabel"));

  bCaloMatching_ = iConfig.getUntrackedParameter<bool>("bCaloMaching", false);
  MB_ = iConfig.getUntrackedParameter<bool>("MB_",true);
  Recenter_ = iConfig.getUntrackedParameter<bool>("Recenter",true);
  EPLevel_ = iConfig.getUntrackedParameter<int>("EPLevel",2);

  nvtx_ = iConfig.getUntrackedParameter<int>("nvtx_", 100);
  reso_ = iConfig.getUntrackedParameter<double>("reso", 0.2);
  if(reso_<0.01) bCaloMatching_ = false;
  if(bCaloMatching_) {
    pfTag_ = iConfig.getUntrackedParameter<edm::InputTag>("pfTag");
    pfToken_ = consumes<reco::PFCandidateCollection>(pfTag_);
  }
  dzdzerror_ = iConfig.getUntrackedParameter<double>("dzdzerror_", 3.);
  d0d0error_ = iConfig.getUntrackedParameter<double>("d0d0error_", 3.);
  pterrorpt_ = iConfig.getUntrackedParameter<double>("pterrorpt_",0.1);
  dzdzerror_pix_ = iConfig.getUntrackedParameter<double>("dzdzerror_pix", 8.); // pixel: nominal 8., tight 6., loose 10
  chi2_  = iConfig.getUntrackedParameter<double>("chi2",12.) ;
  teff = 0;
  if(!effTable_.empty()) teff = new TrackEfficiency(effTable_.data());
  minvz_ = iConfig.getUntrackedParameter<double>("minvz_", -15.);
  maxvz_ = iConfig.getUntrackedParameter<double>("maxvz_", 15.);
  offsetFileName = iConfig.getUntrackedParameter<std::string>("offsetFile");
  frecenter = new TFile(offsetFileName.data(),"read");
  int mx = ntrkbins;
  if(!useNtrk_) {
    mx = ncentbins;
  }
  for(int i = 0; i<mx; i++) {
    for(int j = 1; j<=7; j++){
      wqxtrkRef[j-1][i] = (TH2D *) frecenter->Get(Form("wqxtrk_%d_%d",j,i));
      wqytrkRef[j-1][i] = (TH2D *) frecenter->Get(Form("wqytrk_%d_%d",j,i));
      if(!Recenter_) {
	wqxtrkRef[j-1][i]->Reset();
	wqytrkRef[j-1][i]->Reset();
      }
    }
  }
  
  std::cout<<"==============================================="<<std::endl;
  std::cout<<"centralityBinTag_           "<<centralityBinTag_.encode()<<std::endl;
  std::cout<<"centralityTag_              "<<centralityTag_.encode()<<std::endl;
  std::cout<<"vertexTag_                  "<<vertexTag_.encode()<<std::endl;
  std::cout<<"trackTag_                   "<<trackTag_.encode()<<std::endl;
  std::cout<<"inputPlanesTag_             "<<inputPlanesTag_.encode()<<std::endl;
  std::cout<<"EPOrder_                    "<<EPOrder_<<std::endl;
  std::cout<<"FlatOrder_                  "<<FlatOrder_<<std::endl;
  std::cout<<"NumFlatBins_                "<<NumFlatBins_<<std::endl;
  std::cout<<"caloCentRef_                "<<caloCentRef_<<std::endl;
  std::cout<<"caloCentRefWidth_           "<<caloCentRefWidth_<<std::endl;
  std::cout<<"CentBinCompression_         "<<CentBinCompression_<<std::endl;
  std::cout<<"Noffmin_                    "<<Noffmin_<<std::endl;
  std::cout<<"Noffmax_                    "<<Noffmax_<<std::endl;
  std::cout<<"minrun_                     "<<minrun_<<std::endl;
  std::cout<<"maxrun_                     "<<maxrun_<<std::endl;
  std::cout<<"effTable_                   "<<effTable_<<std::endl;
  std::cout<<"dzerror_                    "<<dzdzerror_<<endl;
  std::cout<<"d0error_                    "<<d0d0error_<<endl;
  std::cout<<"pterrorpt_                    "<<pterrorpt_<<endl;
  std::cout<<"nvtx_                       "<<nvtx_<<endl;
  if(bCaloMatching_) { 
    std::cout<<"bCaloMatching_              true"<<std::endl;
    std::cout<<"reso_                     "<<reso_<<std::endl;   
  }
  std::cout<<"dzdzerror_pix_               "<<dzdzerror_pix_<<std::endl;
  std::cout<<"chi2_                        "<<chi2_<<std::endl;
  std::cout<<"==============================================="<<std::endl;

  TDirectory * save = gDirectory;
  TFileDirectory conddir = fs->mkdir("Conditions");
  conddir.make<TH1I>(centralityBinTag_.label().data(),centralityBinTag_.label().data(),1,0,1);
  conddir.make<TH1I>(centralityTag_.label().data(), centralityTag_.label().data(),1,0,1);
  conddir.make<TH1I>(vertexTag_.label().data(), vertexTag_.label().data(),1,0,1);
  conddir.make<TH1I>(trackTag_.label().data(), trackTag_.label().data(),1,0,1);
  conddir.make<TH1I>(inputPlanesTag_.label().data(), inputPlanesTag_.label().data(),1,0,1);
  string etable = Form("EffTable_%s",effTable_.data());
  conddir.make<TH1I>(etable.data(), etable.data(),1,0,1);
  string note_EPOrder = Form("EPOrder_%d",EPOrder_);
  conddir.make<TH1I>(note_EPOrder.data(), note_EPOrder.data(),1,0,1);
  string note_EPLevel = Form("EPLevel_%d",EPLevel_);
  conddir.make<TH1I>(note_EPLevel.data(), note_EPLevel.data(),1,0,1);
  string note_FlatOrder = Form("FlatOrder_%d",FlatOrder_);
  conddir.make<TH1I>(note_FlatOrder.data(), note_FlatOrder.data(),1,0,1);
  string note_NumFlatBins = Form("NumFlatBins_%d",NumFlatBins_);
  conddir.make<TH1I>(note_NumFlatBins.data(), note_NumFlatBins.data(),1,0,1);
  string note_caloCentRef = Form("caloCentRef_%d",(int)caloCentRef_);
  conddir.make<TH1I>(note_caloCentRef.data(), note_caloCentRef.data(),1,0,1);
  string note_caloCentRefWidth = Form("caloCentRefWidth_%d",(int)caloCentRefWidth_);
  conddir.make<TH1I>(note_caloCentRefWidth.data(), note_caloCentRefWidth.data(),1,0,1);
  string note_dzdzerror = Form("dzdzerror_%07.2f",dzdzerror_);
  conddir.make<TH1I>(note_dzdzerror.data(), note_dzdzerror.data(),1,0,1);
  string note_d0d0error = Form("d0d0error_%07.2f",d0d0error_);
  conddir.make<TH1I>(note_d0d0error.data(), note_d0d0error.data(),1,0,1);
  string note_pterrorpt = Form("dterrorpt_%07.2f",pterrorpt_);
  conddir.make<TH1I>(note_pterrorpt.data(), note_pterrorpt.data(),1,0,1);
  string note_dzdzerror_pix = Form("dzdzerror_pix_%07.2f",dzdzerror_pix_);
  conddir.make<TH1I>(note_dzdzerror_pix.data(), note_dzdzerror_pix.data(),1,0,1);
  string note_chi2 = Form("chi2_%07.2f",chi2_);
  conddir.make<TH1I>(note_chi2.data(), note_chi2.data(),1,0,1);
  if(Recenter_) {
    conddir.make<TH1I>("RecenterTracks", "RecenterTracks",1,0,1);
  } else {
    conddir.make<TH1I>("DoNotRecenterTracks", "DoNotRecenterTracks",1,0,1);
  }
  if(MB_) {
    conddir.make<TH1I>("MB_Set_True", "MB_Set_True",1,0,1);
  } else {
    conddir.make<TH1I>("MB_Set_False", "MB_Set_False",1,0,1);
  }
  if(useNtrk_) {
    conddir.make<TH1I>("useNtrk_Set_True", "useNtrk_Set_True",1,0,1);
  } else {
    conddir.make<TH1I>("useNtrk_Set_False", "useNtrk_Set_False",1,0,1);
  }

  save->cd();
  hNtrkoff = fs->make<TH1D>("Ntrkoff","Ntrkoff",1001,0,3000);
  int npt = nptbins;
  qxtrk1 = fs->make<TH2D>("qxtrk1","qxtrk1",npt,ptbins, netabinsDefault, etabinsDefault);
  qytrk1 = fs->make<TH2D>("qytrk1","qytrk1",npt,ptbins, netabinsDefault, etabinsDefault);
  qxtrk2 = fs->make<TH2D>("qxtrk2","qxtrk2",npt,ptbins, netabinsDefault, etabinsDefault);
  qytrk2 = fs->make<TH2D>("qytrk2","qytrk2",npt,ptbins, netabinsDefault, etabinsDefault);
  qxtrk3 = fs->make<TH2D>("qxtrk3","qxtrk3",npt,ptbins, netabinsDefault, etabinsDefault);
  qytrk3 = fs->make<TH2D>("qytrk3","qytrk3",npt,ptbins, netabinsDefault, etabinsDefault);
  qxtrk4 = fs->make<TH2D>("qxtrk4","qxtrk4",npt,ptbins, netabinsDefault, etabinsDefault);
  qytrk4 = fs->make<TH2D>("qytrk4","qytrk4",npt,ptbins, netabinsDefault, etabinsDefault);
  qxtrk5 = fs->make<TH2D>("qxtrk5","qxtrk5",npt,ptbins, netabinsDefault, etabinsDefault);
  qytrk5 = fs->make<TH2D>("qytrk5","qytrk5",npt,ptbins, netabinsDefault, etabinsDefault);
  qxtrk6 = fs->make<TH2D>("qxtrk6","qxtrk6",npt,ptbins, netabinsDefault, etabinsDefault);
  qytrk6 = fs->make<TH2D>("qytrk6","qytrk6",npt,ptbins, netabinsDefault, etabinsDefault);
  qxtrk7 = fs->make<TH2D>("qxtrk7","qxtrk7",npt,ptbins, netabinsDefault, etabinsDefault);
  qytrk7 = fs->make<TH2D>("qytrk7","qytrk7",npt,ptbins, netabinsDefault, etabinsDefault);
  qcnt =  fs->make<TH2D>("qcnt", "qcnt",npt,ptbins, netabinsDefault, etabinsDefault);
  avpt =  fs->make<TH2D>("avpt","avpt",npt,ptbins, netabinsDefault, etabinsDefault);
  qxtrk1->SetOption("colz");
  qytrk1->SetOption("colz");
  qxtrk2->SetOption("colz");
  qytrk2->SetOption("colz");
  qxtrk3->SetOption("colz");
  qytrk3->SetOption("colz");
  qxtrk4->SetOption("colz");
  qytrk4->SetOption("colz");
  qxtrk5->SetOption("colz");
  qytrk5->SetOption("colz");
  qxtrk6->SetOption("colz");
  qytrk6->SetOption("colz");
  qxtrk7->SetOption("colz");
  qytrk7->SetOption("colz");
  qcnt->SetOption("colz");
  avpt->SetOption("colz");
  qxtrk1->Sumw2();
  qytrk1->Sumw2();
  qxtrk2->Sumw2();
  qytrk2->Sumw2();
  qxtrk3->Sumw2();
  qytrk3->Sumw2();
  qxtrk4->Sumw2();
  qytrk4->Sumw2();
  qxtrk5->Sumw2();
  qytrk5->Sumw2();
  qxtrk6->Sumw2();
  qytrk6->Sumw2();
  qxtrk7->Sumw2();
  qytrk7->Sumw2();
  qcnt->Sumw2();
  avpt->Sumw2();



  hcent = fs->make<TH1D>("cent","cent",220,-10,110);
  hcentbins = fs->make<TH1D>("centbins","centbins",201,0,200);
  hcentres = fs->make<TH1D>("centres","centres",ncentbins,centbins);
  hrun = fs->make<TH1I>("runs","runs",maxrun_-minrun_+1,minrun_,maxrun_);
  hptNtrk = fs->make<TH1D>("ptNtrk","ptNtrk",npt,ptbins);
  hptNtrk->SetXTitle("p_{T} (GeV/c)");
  hptNtrk->SetYTitle("Ntrks (|#eta|<1; 0-5)");
  hptNtrkGood = fs->make<TH1D>("ptNtrkGood","ptNtrkGood",npt,ptbins);
  hptNtrkGood->SetXTitle("p_{T} (GeV/c)");
  hptNtrkGood->SetYTitle("Ntrks (Good) (|#eta|<1; 0-5)");
  hNtrkRet = fs->make<TH1I>("NtrkRet","NtrkRet", 10,0,10);
  for(int i = 0; i<mx; i++) {
    TString hn = Form("Eff_%d_%d",(int)trkBins[i],(int)trkBins[i+1]);
    if(!useNtrk_) {
      hn = Form("Eff_%d_%d",(int)centbins[i],(int)centbins[i+1]);
    }
    hEff[i] = fs->make<TH2D>(hn.Data(),hn.Data(),50,-TMath::Pi(),TMath::Pi(),50,-2.4,2.4);
    hEff[i]->Sumw2();
    hEff[i]->SetXTitle("#phi (radians)");
    hEff[i]->SetYTitle("#eta");
    hEff[i]->SetOption("colz");
  }
  TString epnames = EPNames[0].data();
  epnames = epnames+"/D";
  for(int i = 0; i<NumEPNames; i++) {
    if(i>0) epnames = epnames + ":" + EPNames[i].data() + "/D";
    TFileDirectory subdir = fs->mkdir(Form("%s",EPNames[i].data()));
    flat[i] = new HiEvtPlaneFlatten();
    flat[i]->init(FlatOrder_,NumFlatBins_,EPNames[i],EPOrder[i]);
    Double_t psirange = 4;
    if(EPOrder[i]==1 ) psirange = 3.5;
    if(EPOrder[i]==2 ) psirange = 2;
    if(EPOrder[i]==3 ) psirange = 1.5;
    if(EPOrder[i]==4 ) psirange = 1;
    if(EPOrder[i]==5) psirange = 0.8;
    if(EPOrder[i]==6) psirange = 0.6;
    if(EPOrder[i]==7) psirange = 0.5;

    hPsi[i] = subdir.make<TH1D>("psi","psi",800,-psirange,psirange);
    hPsi[i]->SetXTitle("#Psi");
    hPsi[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));
    
    hPsiOffset[i] = subdir.make<TH1D>("psiOffset","psiOffset",800,-psirange,psirange);
    hPsiOffset[i]->SetXTitle("#Psi");
    hPsiOffset[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));

    
    hPsiFlat[i] = subdir.make<TH1D>("psiFlat","psiFlat",800,-psirange,psirange);
    hPsiFlat[i]->SetXTitle("#Psi");
    hPsiFlat[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));

  }
  std::cout<<"Done with flat init"<<std::endl;

  tree = fs->make<TTree>("tree","EP tree");

  tree->Branch("Cent",&centval,"cent/D");
  tree->Branch("NtrkOff",&Noff,"Noff/I");
  tree->Branch("ntrkflat",&ntrkval,"nofftrak/I");
  tree->Branch("Vtx",&vtx,"vtx/D");
  tree->Branch("epang",&epang, epnames.Data());
  tree->Branch("eporig",&eporig, epnames.Data());
  tree->Branch("qx",      &qx,       epnames.Data());
  tree->Branch("qy",      &qy,       epnames.Data());
  tree->Branch("q",       &q,       epnames.Data());
  tree->Branch("vn", &vn, epnames.Data());
  tree->Branch("mult",    &epmult,  epnames.Data());
  tree->Branch("sumw",    &sumw,  epnames.Data());
  tree->Branch("sumw2",    &sumw2,  epnames.Data());
  tree->Branch("Run",     &runno_,   "run/i");
  tree->Branch("Rescor",  &rescor,   epnames.Data());
  tree->Branch("RescorErr",  &rescorErr,   epnames.Data());
  tree->Branch("qxtrk1",  "TH2D",  &qxtrk1, 128000, 0);
  tree->Branch("qytrk1",  "TH2D",  &qytrk1, 128000, 0);
  tree->Branch("qxtrk2",  "TH2D",  &qxtrk2, 128000, 0);
  tree->Branch("qytrk2",  "TH2D",  &qytrk2, 128000, 0);
  tree->Branch("qxtrk3",  "TH2D",  &qxtrk3, 128000, 0);
  tree->Branch("qytrk3",  "TH2D",  &qytrk3, 128000, 0);
  tree->Branch("qxtrk4",  "TH2D",  &qxtrk4, 128000, 0);
  tree->Branch("qytrk4",  "TH2D",  &qytrk4, 128000, 0);
  tree->Branch("qxtrk5",  "TH2D",  &qxtrk5, 128000, 0);
  tree->Branch("qytrk5",  "TH2D",  &qytrk5, 128000, 0);
  tree->Branch("qxtrk6",  "TH2D",  &qxtrk6, 128000, 0);
  tree->Branch("qytrk6",  "TH2D",  &qytrk6, 128000, 0);
  tree->Branch("qxtrk7",  "TH2D",  &qxtrk7, 128000, 0);
  tree->Branch("qytrk7",  "TH2D",  &qytrk7, 128000, 0);
  tree->Branch("qcnt",    "TH2D",  &qcnt, 128000, 0);
  tree->Branch("avpt",    "TH2D",  &avpt, 128000, 0);
}



VNAnalyzer::~VNAnalyzer()
{
  frecenter->Close();  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VNAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  Bool_t newrun = kFALSE;
  if(runno_ != iEvent.id().run()) newrun = kTRUE;
  runno_ = iEvent.id().run();
  hrun->Fill(runno_);
  if(FirstEvent_ || newrun) {
    FirstEvent_ = kFALSE;
    newrun = kFALSE;
    //
    //Get Size of Centrality Table
    //
    if(!useNtrk_) {
      edm::ESHandle<CentralityTable> centDB_;
      iSetup.get<HeavyIonRcd>().get(centralityLabel_,centDB_);
      nCentBins_ = (int) centDB_->m_table.size();
      for(int i = 0; i<NumEPNames; i++) {
	flat[i]->setCaloCentRefBins(-1,-1);
	if(caloCentRef_>0) {
	  int minbin = (caloCentRef_-caloCentRefWidth_/2.)*nCentBins_/100.;
	  int maxbin = (caloCentRef_+caloCentRefWidth_/2.)*nCentBins_/100.;
	  minbin/=CentBinCompression_;
	  maxbin/=CentBinCompression_;
	  if(minbin>0 && maxbin>=minbin) {
	    if(EPDet[i]==HF || EPDet[i]==Castor) flat[i]->setCaloCentRefBins(minbin,maxbin);
	  }
	}
      }
    }
    //
    //Get flattening parameter file.  
    //
    edm::ESHandle<RPFlatParams> flatparmsDB_;
    iSetup.get<HeavyIonRPRcd>().get(flatparmsDB_);
    LoadEPDB * db = new LoadEPDB(flatparmsDB_,flat);
    if(!db->IsSuccess()) {
 	std::cout<<"Flattening db inconsistancy, will continue without: "<<std::endl;
     loadDB_ = kFALSE;
    }        
  } //First event
  
  
  // //
  // //Get Centrality
  // //

   int Noff=0;

  int bin = 0;
  int centres = 0;
  if(!useNtrk_) {
    int ntrkval=0;
    if(Noffmin_>=0) {
      iEvent.getByToken(centralityToken, centrality_);
      ntrkval = centrality_->Ntracks();
      if ( (ntrkval < Noffmin_) || (ntrkval >= Noffmax_) ) {
	return;
      }
    }

   iEvent.getByToken(centralityBinToken, cbin_);
   int cbin = *cbin_;
   bin = cbin/CentBinCompression_; 
   double cscale = 100./nCentBins_;
   centval = cscale*cbin;
   centres = hcentres->FindBin(centval)-1;
   hcentres->Fill(centval);
  } else {
    iEvent.getByToken(tag_,cbin_);
    ntrkval = *cbin_;
    hNtrkoff->Fill(ntrkval);
    bin = NtrkToBin(ntrkval)-1;
    centval = bin;
    centres = bin;
  }

  hcent->Fill(centval);
  hcentbins->Fill(bin);
  // //
  // //Get Vertex
  // //
  iEvent.getByToken(vertexToken,vertex_);
  if(!vertex_.isValid()) {
    std::cout<<"Error! Can't get vertex!"<<std::endl;
    return;
  }

  iEvent.getByToken(vertexToken, vertex_);
  VertexCollection recoVertices = *vertex_;
  if ( recoVertices.size() > 100 ) return;
  sort(recoVertices.begin(), recoVertices.end(), [](const reco::Vertex &a, const reco::Vertex &b){
      if ( a.tracksSize() == b.tracksSize() ) return a.chi2() < b.chi2();
      return a.tracksSize() > b.tracksSize();
    });
  
  int primaryvtx = 0;
  
  double vz = recoVertices[primaryvtx].z();
  if (fabs(vz) < -15 || fabs(vz) > 15) {
    return;
  }
  vtx = vz; 
  // //
  // //Get Event Planes
  // //
  iEvent.getByToken(inputPlanesToken,inputPlanes_);
  
  if(!inputPlanes_.isValid()){
     cout << "Error (VNAnalyzer)! Can't get hiEvtPlaneFlat product!" << endl;
     return ;
   }
  
   Int_t indx = 0;
   for(int i = 0; i<NumEPNames; i++) {
     epang[i] = -10;
     epsin[i] = 0;
     epcos[i] = 0;
     qx[i] = 0;
     qy[i] = 0;
     q[i] = 0;
     vn[i] = 0;
     epmult[i] = 0;
     sumw[i] = 0;
     sumw2[i] = 0;
   }
   for (EvtPlaneCollection::const_iterator rp = inputPlanes_->begin();rp !=inputPlanes_->end(); rp++) {
     if(indx != rp->indx() ) {
       cout<<"EP inconsistency found. Abort."<<endl;
       return;
     }
     if((rp->sumSin()!=0 || rp->sumCos()!=0) && rp->angle(0)>-4) {
       epang[indx]=rp->angle(EPLevel_);
       eporig[indx]=rp->angle(0);
       epsin[indx] = rp->sumSin(EPLevel_);
       epcos[indx] = rp->sumCos(EPLevel_);
       if(rp->mult()>3 && fabs(vtx)<15) {
	 hPsi[indx]->Fill(rp->angle(0));
	 hPsiOffset[indx]->Fill(rp->angle(1));
	 hPsiFlat[indx]->Fill(rp->angle(2));
       }
      
       qx[indx]  = rp->qx(EPLevel_); 
       qy[indx]  = rp->qy(EPLevel_);
       q[indx]   = rp->q(EPLevel_);
       vn[indx] = rp->vn(0);
       epmult[indx] = (double) rp->mult();
       sumw[indx] = rp->sumw();
       sumw2[indx] = rp->sumw2();
       
       rescor[indx] = flat[indx]->getCentRes1((int) centval);
       rescorErr[indx] = flat[indx]->getCentResErr1((int) centval);
     }
     ++indx; 
   }


  ntrkval = Noff;
  if ( Noff == -2 ) {
    return;
  }
  ntrkval=getNoff(iEvent, iSetup, centres);
  hNtrkoff->Fill(ntrkval);
  
  tree->Fill(); 
}



// ------------ method called once each job just before starting event loop  ------------
void 
VNAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VNAnalyzer::endJob() {
}
bool
VNAnalyzer::CaloMatch(const reco::Track & track, const edm::Event & iEvent, unsigned int idx)
{
  if ( !bCaloMatching_ ) return true;
  edm::Handle<reco::PFCandidateCollection> pfCand;
  iEvent.getByToken( pfToken_, pfCand );
  double energy = 0;
  for ( reco::PFCandidateCollection::const_iterator it = pfCand->begin(); it != pfCand->end(); ++it ) {
    reco::TrackRef trackRef = it->trackRef();
    if ( !((it->particleId() != reco::PFCandidate::h) ||
	 (it->particleId() != reco::PFCandidate::e) ||
	   (it->particleId() != reco::PFCandidate::mu) )) continue;
    if ( idx == trackRef.key() ) {
      energy = it->ecalEnergy() + it->hcalEnergy();
      break;
    }
  }
  
  if( track.pt() < 20 || ( energy/( track.pt()*TMath::CosH(track.eta() ) ) > reso_ && (energy)/(TMath::CosH(track.eta())) > (track.pt() - 80.0) )  ) return true;
  else {
    return false;
  }
}
///
bool
VNAnalyzer::TrackQuality_ppReco(const reco::TrackCollection::const_iterator& itTrack, const reco::VertexCollection& recoVertices)
{
        if ( itTrack->charge() == 0 ) return false;
        if ( !itTrack->quality(reco::TrackBase::highPurity) ) return false;
        if ( itTrack->ptError()/itTrack->pt() > pterrorpt_ ) return false;
	int primaryvtx = 0;
	math::XYZPoint v1( recoVertices[primaryvtx].position().x(), recoVertices[primaryvtx].position().y(), recoVertices[primaryvtx].position().z() );
	double vxError = recoVertices[primaryvtx].xError();
	double vyError = recoVertices[primaryvtx].yError();
	double vzError = recoVertices[primaryvtx].zError();
        double d0 = -1.* itTrack->dxy(v1);
        double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
        if ( fabs( d0/derror ) > d0d0error_ ) {
                return false;
        }
        double dz=itTrack->dz(v1);
        double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
        if ( fabs( dz/dzerror ) > dzdzerror_ ) {
                return false;
        }
        return true;
}

///
bool
VNAnalyzer::TrackQuality_HIReco(const reco::TrackCollection::const_iterator& itTrack, const reco::VertexCollection& recoVertices)
{
	if ( itTrack->charge() == 0 ) return false;
	if ( !itTrack->quality(reco::TrackBase::highPurity) ) return false;
	if ( itTrack->numberOfValidHits() < 11 ) return false;
	if ( itTrack->normalizedChi2() / itTrack->hitPattern().trackerLayersWithMeasurement() > 0.15 ) {
		return false;
	}
	if ( itTrack->ptError()/itTrack->pt() > pterrorpt_ ) {
		return false;
	}
	if (
		itTrack->originalAlgo() != 4 and
		itTrack->originalAlgo() != 5 and
		itTrack->originalAlgo() != 6 and
		itTrack->originalAlgo() != 7
	) {
		return false;
	}

	int primaryvtx = 0;
	math::XYZPoint v1( recoVertices[primaryvtx].position().x(), recoVertices[primaryvtx].position().y(), recoVertices[primaryvtx].position().z() );
	double vxError = recoVertices[primaryvtx].xError();
	double vyError = recoVertices[primaryvtx].yError();
	double vzError = recoVertices[primaryvtx].zError();
	double d0 = -1.* itTrack->dxy(v1);
	double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
	if ( fabs( d0/derror ) > d0d0error_ ) {
		return false;
	}

	double dz=itTrack->dz(v1);
	double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
	if ( fabs( dz/dzerror ) > dzdzerror_ ) {
		return false;
	}
	return true;
}

///
bool
VNAnalyzer::TrackQuality_Pixel(const reco::TrackCollection::const_iterator& itTrack, const reco::VertexCollection& recoVertices)
{
	if ( itTrack->charge() == 0 ) return false;
	if ( !itTrack->quality(reco::TrackBase::highPurity) ) return false;
	bool bPix = false;
	int nHits = itTrack->numberOfValidHits();

	int primaryvtx = 0;
	math::XYZPoint v1( recoVertices[primaryvtx].position().x(), recoVertices[primaryvtx].position().y(), recoVertices[primaryvtx].position().z() );
	double vxError = recoVertices[primaryvtx].xError();
	double vyError = recoVertices[primaryvtx].yError();
	double vzError = recoVertices[primaryvtx].zError();
	double d0 = -1.* itTrack->dxy(v1);

	double dz=itTrack->dz(v1);
	double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
//	std::cout << __LINE__ << "\tnHits = " << nHits << std::endl;
	if ( itTrack->pt() < 2.4 and (nHits==3 or nHits==4 or nHits==5 or nHits==6) ) bPix = true;
	if ( not bPix ) {
		if ( nHits < 11 ) return false;
		if ( itTrack->normalizedChi2() / itTrack->hitPattern().trackerLayersWithMeasurement() > 0.15 ) {
			return false;
		}
		if ( itTrack->ptError()/itTrack->pt() > pterrorpt_ ) {
			return false;
		}
		if (
			itTrack->pt() > 2.4 and
			itTrack->originalAlgo() != 4 and
			itTrack->originalAlgo() != 5 and
			itTrack->originalAlgo() != 6 and
			itTrack->originalAlgo() != 7
		) {
			return false;
		}

		double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
		if ( fabs( d0/derror ) > d0d0error_ ) {
			return false;
		}

		if ( fabs( dz/dzerror ) > dzdzerror_ ) {
			return false;
		}
	} else {
		if ( itTrack->normalizedChi2() / itTrack->hitPattern().trackerLayersWithMeasurement() > chi2_ ) return false;
		if ( fabs( dz/dzerror ) > dzdzerror_pix_ ) {
			return false;
		}
	}
	return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(VNAnalyzer);

