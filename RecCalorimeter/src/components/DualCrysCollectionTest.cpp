#include "DualCrysCollectionTest.h"



#include <Gaudi/PluginServiceV2.h>
#include <GaudiKernel/ISvcLocator.h>

#include <cctype>
#include <cstdlib>  // abs
#include <edm4hep/TimeSeriesCollection.h>
#include <k4FWCore/Transformer.h>
#include <tuple>
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "GaudiKernel/MsgStream.h"
#include "edm4hep/CalorimeterHit.h"
#include "DualCrysCalorimeterHit.h"
#include "edm4hep/Constants.h"

#include "DualCrysSiPMConstants.h"

DECLARE_COMPONENT(DualCrysCollectionTest)






DualCrysCollectionTest::DualCrysCollectionTest(const std::string &name, ISvcLocator *svcLoc)
: MultiTransformer(name, svcLoc, 
		   {
		     KeyValues("CALCollection", {"EcalHitCollection"}),
		     KeyValues("HeaderName", {"EventHeader"}),
		   },
		   {
		     KeyValues("ScintillationSiPMWaveforms", {"ScintWaveforms"}),
		     KeyValues("CherenkovSiPMWaveforms", {"CherenkovWaveforms"})}
		   ) {
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
  }

  m_geoSvc = serviceLocator()->service("GeoSvc");  // important to initialize m_geoSvc
}



StatusCode DualCrysCollectionTest::initialize() {

    // Initialize random services
  m_randSvc = service("RndmGenSvc", false);

  if (!m_randSvc) {
    error() << "Couldn't get RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_rndmUniform.initialize(m_randSvc, Rndm::Flat(0., 1.)).isFailure()) {
    error() << "Couldn't initialize RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::TimeSeriesCollection, 
	   edm4hep::TimeSeriesCollection>
DualCrysCollectionTest::operator()(const edm4hep::SimCalorimeterHitCollection &simCaloHits,
				   const edm4hep::EventHeaderCollection& headers) const {

  
  debug() << " process event : " << headers[0].getEventNumber() << " - run  " << headers[0].getRunNumber()
          << endmsg;  // headers[0].getRunNumber(),headers[0].getEventNumber()

  edm4hep::TimeSeriesCollection cherenkovWaveforms;
  edm4hep::TimeSeriesCollection scintillationWaveforms;

  auto *detector = m_geoSvc->getDetector();
  auto &constants = detector->constants();
  debug() << "Count of constants:" << constants.size() << endmsg;
  debug() << "Hit Count:" << simCaloHits.size() << endmsg; 
  for (auto &[k, v] : constants) {
    debug() << "Detector Constant:" << k << endmsg;
  }

  std::string initString;
  initString =       "system:3,ix:-7,iy:-7,slice:3,layer:3,wc:3";
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);  // check!


  std::vector<photon> sPhotons;
  std::vector<photon> cPhotons; 
  
  for (const auto& hit : simCaloHits)
    {
      auto photonTuple = processHit(hit);
      sPhotons.insert(end(sPhotons), begin(std::get<0>(photonTuple)), end(std::get<0>(photonTuple))); 
      cPhotons.insert(end(cPhotons), begin(std::get<1>(photonTuple)), end(std::get<1>(photonTuple))); 

    }

  // Now we have all the photons for this event
  info() << "Found " << sPhotons.size() << " Scintillation photons." << endmsg;
  info() << "Found " << cPhotons.size() << " Cherenkov photons." << endmsg;


  auto wvlfilter = [](photon &p) { return p.wavelength <= 300.0 || p.wavelength >= 1000.0;}; 
  
  std::erase_if(sPhotons, wvlfilter);
  std::erase_if(cPhotons, wvlfilter); 
  info() << "After Erasing " << sPhotons.size() << " Scintillation photons." << endmsg;
  info() << "After Erasing " << cPhotons.size() << " Cherenkov photons." << endmsg;


  auto respfilter = [&](photon &p) {
      double wvl, response; 
      std::tie(wvl,response) = findnearest(SiPM_Type::RGB, p.wavelength);
      double  randval = m_rndmUniform.shoot();
      if (randval > response)
	return true;
      return false;
  }; 


  std::erase_if(sPhotons, respfilter);
  std::erase_if(cPhotons, respfilter);

  // double  randval;
  // double wvl, response; 
  // for (auto it = sPhotons.begin(); it != sPhotons.end(); ++it) {
  //   std::tie(wvl,response) = findnearest(SiPM_Type::RGB, (*it).wavelength);
  //   randval =  m_rndmUniform.shoot();
  //   if (randval > response)
  //     sPhotons.erase(it);
  // }
  // for (auto it = cPhotons.begin(); it != cPhotons.end(); ++it) {
  //   std::tie(wvl,response) = findnearest(SiPM_Type::RGB, (*it).wavelength);
  //   randval =  m_rndmUniform.shoot();
  //   if (randval > response)
  //     cPhotons.erase(it);
  // }

  info() << "After Response Cut " << sPhotons.size() << " Scintillation photons." << endmsg;
  info() << "After Response Cut " << cPhotons.size() << " Cherenkov photons." << endmsg;


  std::vector<double> xs(1024);
  std::vector<double> scintSignal(1024);
  std::vector<double> cherenkovSignal(1024);
  double dt = 0.2; // sampling time in ns 
  for (size_t i= 0; i< 1024; i++) {
    xs[i] = dt*i;
  }

  
  auto fillWaveform = [&](edm4hep::MutableTimeSeries &ts,
			  std::vector<double> &wv, int ix, int iy, int layer) {
    ts.setInterval(dt);
    ts.setTime(0);
    ts.setCellID((ix<<3)|(iy<<10)|(layer<<20));

    for (size_t i = 0; i < 1024; i++) {
      ts.addToAmplitude(wv[i]);
    }

  }; 

  std::map<key,std::vector<double>> CherenWaveforms;
  std::map<key,std::vector<double>> ScintWaveforms;

  // create waveform positions 
  for (auto &p : sPhotons) {
    int idx = int(round(p.time/dt));
    key k{.ix = p.ix,
      .iy = p.iy,
      .layer = p.layer}; 

    if (ScintWaveforms.count(k) == 0) {
      auto &vec = ScintWaveforms[k];
      vec.resize(1024);
    }
    auto &pvec = ScintWaveforms[k]; 
    
    for (; idx < 1024; idx++) {
      double offset = xs[idx]-p.time;
      if (offset < 0.0)
	offset = 0; 
      pvec[idx] += DESY_SPR(offset); 
    }


  }
  for (auto &p : cPhotons) {
    int idx = int(round(p.time/dt));
    key k{};
    k.ix = p.ix;
    k.iy = p.iy;
    k.layer = p.layer; 
    if (CherenWaveforms.count(k) == 0) {
      auto &vec = CherenWaveforms[k];
      vec.resize(1024);
    }
    auto &cvec = CherenWaveforms[k]; 
    
    for (; idx < 1024; idx++) {
      double offset = xs[idx]-p.time;
      if (offset < 0.0)
	offset = 0; 
      cvec[idx] += DESY_SPR(offset); 
    }

    
  }

  for (auto &[k,v] : ScintWaveforms) {
    int ix = k.ix;
    int iy = k.iy;
    int layer = k.layer; 
    
    auto wv = scintillationWaveforms.create();

    fillWaveform(wv, v, ix, iy,layer); 
  }

  for (auto &[k,v] : CherenWaveforms) {

    int ix = k.ix;
    int iy = k.iy;
    int layer = k.layer; 
    auto wv = cherenkovWaveforms.create();

    fillWaveform(wv, v, ix, iy,layer); 
  }
  
  
  
  return std::make_tuple(std::move(cherenkovWaveforms), std::move(scintillationWaveforms)); 
  
}


bool DualCrysCollectionTest::useLayer(CHT::Layout caloLayout, unsigned int layer) const {
  switch (caloLayout) {
    case CHT::ecal:
      if (layer > m_useLayersEcalVec.size() || m_useLayersEcalVec.size() == 0)
        return true;
      return m_useLayersEcalVec[layer];  //break not needed, because of return
    case CHT::hcal:
      if (layer > m_useLayersHcalVec.size() || m_useLayersHcalVec.size() == 0)
        return true;
      return m_useLayersHcalVec[layer];  //break not needed, because of return
      //For all other cases, always keep the hit
    default:
      return true;
  }
}  //useLayer


typedef std::vector<DualCrysCollectionTest::photon> PhotonVector;
std::tuple<PhotonVector, PhotonVector> 
DualCrysCollectionTest::processHit(const edm4hep::SimCalorimeterHit &hit) const {

  PhotonVector sPhotons;
  PhotonVector cPhotons; 
  const int cellID = hit.getCellID();
  int slice_id = ((0x7<<17&cellID)>>17);
  int layer_id = ((0x7<<20&cellID)>>20);
  int ix = (0x7f<<3&cellID)>>3;
  int iy = (0x7f<<10&cellID)>>10;
  if (ix > 63)
    ix -= 64;
  if (iy > 63)
    iy -= 64;
  

  
  bool first_pd = (slice_id == 4) && (layer_id == 1);
  bool second_pd = (slice_id == 1) && (layer_id == 0);
  if (first_pd || second_pd) {
    if (hit.isAvailable()) {
      for (auto step = hit.contributions_begin();
	   step != hit.contributions_end(); step++) {
	
	float energy = step->getEnergy()/CLHEP::eV;
	double wavelength = 1239.84187 / (1000*energy); 
	photon p {
	  .wavelength = wavelength,
	  .ix = ix,
	  .iy = iy,
	  .layer = layer_id,
	  .time = step->getTime(),
		.photon_type = step->getPDG()
	} ;
	if (step->getPDG() == -22)
	  sPhotons.push_back(p);
	else
	  cPhotons.push_back(p); 
      }
    }
  }
  return std::make_tuple(sPhotons, cPhotons); 
  

}
