//
// Simplied analyzer for source calibration analysis - collects CaloHit and CaloMCHit
// Author: Sophie Middleton 2022, 2024
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"

#include "Offline/CaloCluster/inc/ClusterUtils.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"

#include "TTree.h"
namespace
{
  constexpr int ntupLen = 16384;
}

int Contains(std::vector<int> v, int x)
{
  return std::count(v.begin(), v.end(), x);
}
namespace mu2e {

    class SourceCalibAna : public art::EDAnalyzer {

      public:
      struct Config
      {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<art::InputTag> caloHitCollection{Name("CaloHitCollection"),Comment("Calo Hit collection name")};
        fhicl::Atom<art::InputTag> caloHitTruthCollection{Name("CaloHitMCCollection"),Comment("CaloHit truth name")};
        fhicl::Atom<int> diagLevel{Name("diagLevel"),Comment("Diag Level"),0};
      };
      explicit SourceCalibAna(const art::EDAnalyzer::Table<Config>& config);
      virtual ~SourceCalibAna() {}

      virtual void beginJob();
      virtual void endJob() {};
      virtual void analyze(const art::Event& e);


      private:

        art::InputTag         caloHitTag_;
        art::InputTag         caloHitTruthTag_;

        int                   diagLevel_;
        int                   nProcess_;

        TTree* Ntup_;
        unsigned int _evt,_run;

        unsigned int ncalhitHit_, nCrystals_;
        float truetotalEnergyDep_;
        int cryId_[ntupLen],crySectionId_[ntupLen],crycalhitMCIdx_[ntupLen],crySimLen_[ntupLen];
        float cryEtot_,calhitRecoTime_[ntupLen],calhitRecoEdep_[ntupLen],calhitRecoEdepErr_[ntupLen],calhitRecoPosX_[ntupLen],calhitRecoPosY_[ntupLen],calhitRecoPosZ_[ntupLen],_cryLeak[ntupLen];

        unsigned int ncalhitMCHit_,crycalhitMCId_[ntupLen];
        
        int crycalhitMCPdgId_[ntupLen],crycalhitMCCrCode_[ntupLen],crycalhitMCGenIdx_[ntupLen],calhitConv_[ntupLen];
        float crycalhitMCMom_[ntupLen],crycalhitMCX_[ntupLen],crycalhitMCY_[ntupLen],crycalhitMCZ_[ntupLen],crycalhitMCStartT_[ntupLen];
        float crycalhitMCEndX_[ntupLen],crycalhitMCEndY_[ntupLen],crycalhitMCEndZ_[ntupLen],crycalhitMCEndT_[ntupLen];
        float crycalhitTime_[ntupLen],crycalhitMCEdep_[ntupLen],calhitRecoTimeErr_[ntupLen],calhitRecoT1_[ntupLen],calhitRecoT2_[ntupLen],calhitRecoT1Err_[ntupLen],calhitRecoT2Err_[ntupLen];

        };


    SourceCalibAna::SourceCalibAna(const art::EDAnalyzer::Table<Config>& config):
      EDAnalyzer{config},
      caloHitTag_         (config().caloHitCollection()),
      caloHitTruthTag_    (config().caloHitTruthCollection()),
      diagLevel_          (config().diagLevel()),
      nProcess_(0),
      Ntup_(0)
      {}

    void SourceCalibAna::beginJob(){

      art::ServiceHandle<art::TFileService> tfs;

      Ntup_  = tfs->make<TTree>("SourceCalibAna","SourceCalibAna");

      Ntup_->Branch("evt",          &_evt ,         "evt/I");
      Ntup_->Branch("run",          &_run ,         "run/I");

      // Reconstructed carystal hit info (from CaloHitCollection):
      Ntup_->Branch("calhitRecoEtot",      &cryEtot_ ,     "calhitRecoEtot/F");
      Ntup_->Branch("ncrystals",         &ncalhitHit_ ,       "ncalhitHit/I");
      Ntup_->Branch("cryId",        &cryId_ ,       "cryId[ncalhitHit]/I");
      Ntup_->Branch("crySectionId", &crySectionId_, "crySectionId[ncalhitHit]/I");
      Ntup_->Branch("calhitRecoPosX",      &calhitRecoPosX_ ,     "calhitRecoPosX[ncalhitHit]/F");
      Ntup_->Branch("calhitRecoPosY",      &calhitRecoPosY_ ,     "calhitRecoPosY[ncalhitHit]/F");
      Ntup_->Branch("calhitRecoPosZ",      &calhitRecoPosZ_ ,     "calhitRecoPosZ[ncalhitHit]/F");
      Ntup_->Branch("calhitRecoEdep",      &calhitRecoEdep_ ,     "calhitRecoEdep[ncalhitHit]/F");
      Ntup_->Branch("calhitRecoEdepErr",   &calhitRecoEdepErr_ ,  "calhitRecoEdepErr[ncalhitHit]/F");
      Ntup_->Branch("calhitRecoTime",      &calhitRecoTime_ ,     "calhitRecoTime[ncalhitHit]/F");
      Ntup_->Branch("calhitRecoTimeErr",   &calhitRecoTimeErr_ ,  "calhitRecoTimeErr[ncalhitHit]/F");
      Ntup_->Branch("calhitRecoT1",        &calhitRecoT1_ ,       "calhitRecoT1[ncalhitHit]/F");
      Ntup_->Branch("calhitRecoT2",        &calhitRecoT2_ ,       "calhitRecoT2[ncalhitHit]/F");
      Ntup_->Branch("calhitRecoT1Err",     &calhitRecoT1Err_ ,    "calhitRecoT1Err[ncalhitHit]/F");
      Ntup_->Branch("calhitRecoT2Err",     &calhitRecoT2Err_ ,    "calhitRecoT2Err[ncalhitHit]/F");
      Ntup_->Branch("calhitConv",      &calhitConv_ ,     "calhitConv[ncalhitHit]/I");

      // Truth crystal hit info (from CaloHitMCCollection):
      Ntup_->Branch("calhitMCEtot",      &truetotalEnergyDep_ ,     "calhitMCEtot/F");
      Ntup_->Branch("ncalhitMC",         &ncalhitMCHit_ ,     "ncalhitMC/I");
      Ntup_->Branch("crycalhitMCIdx",    &crycalhitMCIdx_ ,   "crycalhitMCIdx[ncalhitMC]/I");
      Ntup_->Branch("crySimLen",    &crySimLen_ ,   "crySimLen[ncalhitMC]/I");
      Ntup_->Branch("calhitMCId",        &crycalhitMCId_ ,    "calhitMCId[ncalhitMC]/I");
      Ntup_->Branch("calhitMCPdgId",     &crycalhitMCPdgId_ , "calhitMCPdgId[ncalhitMC]/I");
      Ntup_->Branch("calhitMCCrCode",    &crycalhitMCCrCode_ ,"calhitMCCrCode[ncalhitMC]/I");
      Ntup_->Branch("calhitMCMom",       &crycalhitMCMom_ ,   "calhitMCMom[ncalhitMC]/F");
      Ntup_->Branch("calhitMCX",    &crycalhitMCX_ ,"calhitMCX[ncalhitMC]/F");
      Ntup_->Branch("calhitMCY",    &crycalhitMCY_ ,"calhitMCY[ncalhitMC]/F");
      Ntup_->Branch("calhitMCZ",    &crycalhitMCZ_ ,"calhitMCZ[ncalhitMC]/F");
      Ntup_->Branch("calhitMCStartT",    &crycalhitMCStartT_ ,"calhitMCStartT[ncalhitMC]/F");
      Ntup_->Branch("calhitMCEndX",    &crycalhitMCEndX_ ,"calhitMCEndX[ncalhitMC]/F");
      Ntup_->Branch("calhitMCEndY",    &crycalhitMCEndY_ ,"calhitMCEndY[ncalhitMC]/F");
      Ntup_->Branch("calhitMCEndZ",    &crycalhitMCEndZ_ ,"calhitMCEndZ[ncalhitMC]/F");
      Ntup_->Branch("calhitMCEndT",    &crycalhitMCEndT_ ,"calhitMCEndT[ncalhitMC]/F");
      Ntup_->Branch("calhitTime",      &crycalhitTime_ ,  "calhitTime[ncalhitMC]/F");
      Ntup_->Branch("calhitMCEdep",      &crycalhitMCEdep_ ,  "calhitMCEdep[ncalhitMC]/F");
      Ntup_->Branch("calhitMCGenIdx",    &crycalhitMCGenIdx_ ,"calhitMCGenIdx[ncalhitMC]/I");
    }

    void SourceCalibAna::analyze(const art::Event& event){
        ++nProcess_;
        if (nProcess_%10==0 && diagLevel_ > 0) std::cout<<"Processing event from SourceCalibAna =  "<<nProcess_ <<std::endl;

        //Handle to the calorimeter
        art::ServiceHandle<GeometryService> geom;
        if (!geom->hasElement<Calorimeter>() ) return;
        const Calorimeter& cal = *(GeomHandle<Calorimeter>());

        //Calorimeter crystal hits (average from readouts)
        art::Handle<CaloHitCollection> CaloHitsHandle;
        event.getByLabel(caloHitTag_, CaloHitsHandle);
        const CaloHitCollection& CaloHits(*CaloHitsHandle);

        //Calo digi truth assignment
        art::Handle<CaloHitMCCollection> caloHitMCHandle;
        event.getByLabel(caloHitTruthTag_, caloHitMCHandle);
        const CaloHitMCCollection& caloHitTruthCollection(*caloHitMCHandle);

        _evt = event.id().event();
        _run = event.run();

        if (diagLevel_ == 3){std::cout << "processing event in calo_example " << nProcess_ << " run and event  = " << _run << " " << _evt << std::endl;}

        //--------------------------  Do calorimeter hits --------------------------------
        ncalhitHit_ = ncalhitMCHit_ = 0;
        cryEtot_ = 0.0;
        truetotalEnergyDep_ = 0.0;

        std::vector<int> crystalsHit;
        for (unsigned int ic=0; ic<CaloHits.size();++ic)
        {

          const CaloHit& hit            = CaloHits.at(ic);
          int diskId                    = cal.crystal(hit.crystalID()).diskID();
          CLHEP::Hep3Vector crystalPos  = cal.mu2eToDiskFF(diskId,cal.crystal(hit.crystalID()).position());

          const auto eDepMCs =  caloHitTruthCollection[ic].energyDeposits();


          constexpr float invalid(999.0);
          float calhitRecoT1(invalid),calhitRecoT2(invalid),calhitRecoT1Err(invalid),calhitRecoT2Err(invalid);
          if (hit.recoCaloDigis().size()>1)
          {
            int idx0 = CaloSiPMId(hit.recoCaloDigis().at(0)->SiPMID()).SiPMLocalId();
            int idx1 = CaloSiPMId(hit.recoCaloDigis().at(1)->SiPMID()).SiPMLocalId();
            calhitRecoT1    = hit.recoCaloDigis().at(idx0)->time();
            calhitRecoT2    = hit.recoCaloDigis().at(idx1)->time();
            calhitRecoT1Err = hit.recoCaloDigis().at(idx0)->timeErr();
            calhitRecoT2Err = hit.recoCaloDigis().at(idx1)->timeErr();
          }

          cryId_[ncalhitHit_]        = hit.crystalID();
          if(ic == 0) crystalsHit.push_back(hit.crystalID());
          else if(Contains(crystalsHit, hit.crystalID()) == 0) crystalsHit.push_back(hit.crystalID());
          crySectionId_[ncalhitHit_] = diskId;
          calhitRecoEdep_[ncalhitHit_]      = hit.energyDep();
          calhitRecoEdepErr_[ncalhitHit_]   = hit.energyDepErr();
          calhitRecoTime_[ncalhitHit_]      = hit.time();
          calhitRecoTimeErr_[ncalhitHit_]   = hit.timeErr();
          calhitRecoT1_[ncalhitHit_]        = calhitRecoT1;
          calhitRecoT2_[ncalhitHit_]        = calhitRecoT2;
          calhitRecoT1Err_[ncalhitHit_]     = calhitRecoT1Err;
          calhitRecoT2Err_[ncalhitHit_]     = calhitRecoT2Err;
          GeomHandle<DetectorSystem> det;
          CLHEP::Hep3Vector Mu2ePos = det->toMu2e(crystalPos); // in mu2e coordinates for comparison
          calhitRecoPosX_[ncalhitHit_]      = Mu2ePos.x();
          calhitRecoPosY_[ncalhitHit_]      = Mu2ePos.y();
          calhitRecoPosZ_[ncalhitHit_]      = Mu2ePos.z();

          cryEtot_             += hit.energyDep();

          crycalhitMCIdx_[ncalhitHit_]    = ncalhitMCHit_;
          crySimLen_[ncalhitHit_]    = eDepMCs.size();

          double sumEdepMC(0),edepTime(0);
          for (unsigned i=0;i< eDepMCs.size();++i)
          {
            const auto& eDepMC = eDepMCs[i];

            auto parent(eDepMC.sim());
            while (parent->hasParent()) parent = parent->parent();
            int genId=-1;
            if (parent->genParticle()) genId = parent->genParticle()->generatorId().id();

            crycalhitMCId_[ncalhitMCHit_]      = eDepMC.sim()->id().asInt();
            crycalhitMCPdgId_[ncalhitMCHit_]   = eDepMC.sim()->pdgId();
            crycalhitMCCrCode_[ncalhitMCHit_]  = eDepMC.sim()->creationCode();
            crycalhitTime_[ncalhitMCHit_]    = eDepMC.time();
            crycalhitMCEdep_[ncalhitMCHit_]    = eDepMC.energyDep();
            crycalhitMCMom_[ncalhitMCHit_]     = eDepMC.momentumIn();

            crycalhitMCX_[ncalhitMCHit_]  = parent->startPosition().x();
            crycalhitMCY_[ncalhitMCHit_]  = parent->startPosition().y();
            crycalhitMCZ_[ncalhitMCHit_]  = parent->startPosition().z();
            crycalhitMCStartT_[ncalhitMCHit_]  = parent->startGlobalTime();

            crycalhitMCEndX_[ncalhitMCHit_]    = parent->endPosition().x();
            crycalhitMCEndY_[ncalhitMCHit_]    = parent->endPosition().y();
            crycalhitMCEndZ_[ncalhitMCHit_]    = parent->endPosition().z();
            crycalhitMCEndT_[ncalhitMCHit_]  = parent->endGlobalTime();

            crycalhitMCGenIdx_[ncalhitMCHit_]  = genId;
            ++ncalhitMCHit_;

            sumEdepMC += eDepMC.energyDep();
            truetotalEnergyDep_ += sumEdepMC;
            if (edepTime<1) edepTime = eDepMC.time();
          }
          ++ncalhitHit_;

        }
              nCrystals_ = crystalsHit.size();

        Ntup_->Fill();
  }
}
DEFINE_ART_MODULE(mu2e::SourceCalibAna)
