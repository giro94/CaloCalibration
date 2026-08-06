// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <list>
#include <deque>
//Mu2e Geom:
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "TrackerGeom/inc/Tracker.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeometryService/inc/DetectorSystem.hh"

//Mu2e Data Prods:
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

// Mu2e Utilities
#include "GeometryService/inc/GeomHandle.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"


// ROOT incldues
#include "TLegend.h"
#include "TLatex.h"
#include "TTree.h"
#include "TH2D.h"
#include "TF1.h"
#include "TH3D.h"
#include "Rtypes.h"
#include "TApplication.h"
#include "TArc.h"
#include "TTUBE.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLine.h"
#include "TNtuple.h"
#include "TPolyMarker.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TStyle.h"
#include "TText.h"
#include "TRotMatrix.h"


using namespace std; 

namespace mu2e 
{
  class CaloEventDisplays : public art::EDAnalyzer {
    public:
	    struct Config{
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<art::InputTag> calocrysTag{Name("CaloCrystalHitCollection"),Comment("cal reco crystal hit info")};
        fhicl::Atom<art::InputTag> caloclusterTag{Name("CaloClusterCollection"),Comment("cal reco cluster info")};
        fhicl::Atom<bool> doDisplay{Name("doDisplay"),Comment("use display"), false};
        fhicl::Atom<bool> clickToAdvance{Name("clickToAdvance"),Comment("next event"), false};
	     };
      typedef art::EDAnalyzer::Table<Config> Parameters;

      explicit CaloEventDisplays(const Parameters& conf);

      virtual ~CaloEventDisplays();
      virtual void beginJob();
      virtual void analyze(const art::Event& e) override;
    private: 
	    Config _conf;

	    
	    //For Event Displays:
	    TApplication* application_;
	    TDirectory*   directory_ = nullptr;
	    TCanvas*      canvas_ = nullptr;
	    TH2D* _display = nullptr;
	    TNtuple* _ntTrack = nullptr;
	    TNtuple* _ntHit = nullptr;

	    
      art::InputTag _calocrysTag;
      art::InputTag _caloclusterTag;

      const CaloCrystalHitCollection*  _calcryhitcol;
      const CaloClusterCollection* _calclustercol;

    };

    CaloEventDisplays::CaloEventDisplays(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _calocrysTag(conf().calocrysTag()),
    _caloclusterTag(conf().caloclusterTag()),
    doDisplay_ (conf().doDisplay()),
    clickToAdvance_ (conf().clickToAdvance())
    {}
   
    CaloEventDisplays::~CaloEventDisplays(){}

    void CaloEventDisplays::beginJob() {
      
      art::ServiceHandle<art::TFileService> tfs;
      directory_ = gDirectory;

      TString name  = "test";
      TString title = "test";
      int window_size_x(1300);
      int window_size_y(600);
      canvas_ = tfs->make<TCanvas>(name,title,window_size_x,window_size_y);
      canvas_->Divide(1,2);
   
  }

  void CaloEventDisplays::analyze(const art::Event& event) {
	 
  art::ServiceHandle<GeometryService> geom;
  if( ! geom->hasElement<Calorimeter>() ) return;
  Calorimeter const & cal = *(GeomHandle<Calorimeter>());

    TText  text;     
    TArc arcOut, arcIn;
    TBox box;
    canvas_->SetTitle("foo title");
    auto pad = canvas_->cd();
    pad->Clear();
    canvas_->SetTitle("bar title");

    auto xyplot = pad->DrawFrame(-1000,-1000, 1000,1000);
    xyplot->GetYaxis()->SetTitleOffset(1.25);
    xyplot->SetTitle( "View of Calo Disk 1 in YZ Plane; Z(mm);Y(mm)");

    art::ServiceHandle<mu2e::GeometryService>   geom;

    mu2e::GeomHandle<mu2e::DetectorSystem>      ds;
    mu2e::GeomHandle<mu2e::VirtualDetector>     vdet;

    Hep3Vector vd_tt_back = ds->toDetector(vdet->getGlobal(mu2e::VirtualDetectorId::TT_Back));

    Disk const & disk =  cal.disk(1);
    double outerR = disk.outerRadius();
    double innerR= disk.innerRadius();

    arcOut.SetFillColor(kGray);
    arcIn.SetFillColor(kWhite);
    arcIn.SetLineColor(kGray+1);
    arcOut.SetLineColor(kGray+1);
    arcOut.DrawArc(0.,0., outerR);
    arcIn.DrawArc(0.,0., innerR);

    unsigned int i =0;
    for(int i=0;i<674;i++){
        int crystalID(i);
	      Crystal const &crystal = cal.crystal(i);
        double crystalXLen = crystal.size().x();
        double crystalYLen = crystal.size().y();
        CLHEP::Hep3Vector crystalPos   = cal.mu2eToDiskFF(1,crystal.position());
        box.SetLineColor(kGray+1);
        box.DrawBox(crystalPos.x()-crystalXLen/2, crystalPos.y()-crystalYLen/2,crystalPos.x()+crystalXLen/2, crystalPos.y()+crystalYLen/2);
  int diskId     = cal.crystal(hit->id()).diskId();

        CLHEP::Hep3Vector crystalPos   = cal.mu2eToDiskFF(diskId,cal.crystal(hit->id()).position());

        box.DrawBox(crystalPos.x()-crystalXLen/2, crystalPos.y()-crystalYLen/2,crystalPos.x()+crystalXLen/2, crystalPos.y()+crystalYLen/2);
     
        TLatex latex;
        stringstream crys;
        crys<<round(crystalID);
        const char* str_crys = crys.str().c_str();
        latex.SetTextSize(0.02);
        latex.DrawLatex(crystalPos.x()-crystalXLen/2, crystalPos.y()-crystalYLen/2,str_crys);
      }
  }


}  // end namespace mu2e 

using mu2e::CaloEventDisplays;
DEFINE_ART_MODULE(CaloEventDisplays);
 /*   
TBox box;
	TArc arc;
	auto xyplot = canvas_->DrawFrame(-1000,-1000, 1000,1000);
	xyplot->GetYaxis()->SetTitleOffset(1.25);
	xyplot->SetTitle( "XY; X(mm);Y(mm)");
	art::ServiceHandle<GeometryService> geom;
        if( ! geom->hasElement<Calorimeter>() ) return;
        Calorimeter const & cal = *(GeomHandle<Calorimeter>());
	canvas_->SetTitle("foo title");
	Disk const & disk =  cal.disk(1);
	double disk_outR = disk.outerRadius();
	double disk_inR= disk.innerRadius();
	arc.SetFillStyle(0);
      	arc.DrawArc(0.,0., disk_inR);
      	arc.DrawArc(0.,0., disk_outR);
	for(int i=0;i<NCrystals;i++){
	        Crystal const &crystal = cal.crystal(i);
	   	double crystalXLen = crystal.size().x();
		double crystalYLen = crystal.size().y();
		CLHEP::Hep3Vector crystalPos   = cal.mu2eToDiskFF(1,crystal.position());

		box.DrawBox(crystalPos.x()-crystalXLen/2, crystalPos.y()-crystalYLen/2,crystalPos.x()+crystalXLen/2, crystalPos.y()+crystalYLen/2);
		
	}
	canvas_->Update();
	canvas_->SaveAs("AllCrystalsInAllEvents.root");
	
}
*/
