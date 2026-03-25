#include "DualCrysSiPMAlgo.h"
#include "DualCrysSiPMConstants.h"
#include <CLHEP/Units/SystemOfUnits.h>
#include <GaudiKernel/MsgStream.h>
#include <Math/Interpolator.h>
#include <algorithm>
#include <edm4hep/CalorimeterHit.h>
#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/MutableCalorimeterHit.h>

using namespace calvision;


DECLARE_COMPONENT(DualCrysSiPMAlgo)

DualCrysSiPMAlgo::DualCrysSiPMAlgo(const std::string& name, ISvcLocator* svclocator)
: Gaudi::Algorithm(name, svclocator)
{




}



StatusCode DualCrysSiPMAlgo::initialize()
{

  StatusCode sc = Gaudi::Algorithm::initialize();

  if (sc.isFailure())
    return sc;

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

  if (!calvision::filterInit)
    init_filters(); 

  
  info() << "Dual Crystal SiPM Algorithm Initialized" << endmsg; 

  return StatusCode::SUCCESS;

}


// This builds and returns waveforms for a given set of hits 
StatusCode DualCrysSiPMAlgo::execute(const EventContext&) const
{
  
  const edm4hep::CalorimeterHitCollection* simHits = m_simHits.get();
  //const edm4hep::CaloHitSimCaloHitLinkCollection* linkCollection = m_links.get();

  // Output
  edm4hep::TimeSeriesCollection* waveforms = m_waveforms.createAndPut();
  edm4hep::TimeSeriesCollection* cherenwaveforms = m_cherenwaveforms.createAndPut();
  edm4hep::TimeSeriesCollection* scintwaveforms = m_scintwaveforms.createAndPut();

  edm4hep::CalorimeterHitCollection* killedCherenkovPhts = m_killedCherenPhotons.createAndPut();
  edm4hep::CalorimeterHitCollection* killedScintPhts = m_killedScintPhotons.createAndPut();

  edm4hep::CalorimeterHitCollection* passedScintPhts = m_passedScintPhotons.createAndPut();

  edm4hep::CalorimeterHitCollection* passedCherenkovPhts = 
    m_passedCherenPhotons.createAndPut(); 
  
  info() << "Sim Hit Size:" << simHits->size() << " :: ";
  //info() << "Link Size:" << linkCollection->size() << endmsg; 
  

  std::vector<photon> totalPhotons; 
  
  for (size_t i = 0; i < simHits->size(); i++) {

    const auto &hit = simHits->at(i); 
    const int cellID = hit.getCellID();
    const int slice_id = ((0x7<<17&cellID)>>17);
    const int layer_id = ((0x7<<20&cellID)>>20);
    int ix = (0x7f<<3&cellID)>>3;
    int iy = (0x7f<<10&cellID)>>10;
    if (ix > 63)
      ix -= 64;
    if (iy > 63)
      iy -= 64;
    info() << "ix " << ix << ", iy " << iy;
    info() << ", slice " << slice_id << ", layer" << layer_id;
    info() << endmsg; 
    // not using the slice or layer at the moment,
    // but the upstream feeder should only give us particles from
    // the SiPM tower slice
    // 1.2398 / energy = wavelength hc/e ~ 1.239e-6 eV*m
    // I think the value we get is in GeV, have to check the calculation
    // but it looks ~ right in the results but needs verifying / rewriting
    if (hit.getType() == -22) {
      // scint
      float energy = hit.getEnergy()/CLHEP::eV;
      double wavelength = 1239.84187 / (1000*energy); 
      info() << "Scint @ " << hit.getTime() << " Wavelength:" << wavelength;
      info() << " Energy: " << energy << " eV?"; 
      info() << endmsg;
      photon sphoton {
	.wavelength = wavelength,
	.ix = ix,
	.iy = iy,
	.layer = layer_id,
	.time = hit.getTime(),
	.photon_type = -22,
	.hitidx = i
      };
      totalPhotons.push_back(sphoton); 

    }
    else if ( hit.getType() == -44) {
      // cherenkov
      float energy = hit.getEnergy()/CLHEP::eV;
      double wavelength = 1239.84187 / (1000*energy); 
      info() << "Cherenk @ " << hit.getTime() << " Wavelength:" << wavelength;
      info() << " Energy: " << energy << " eV?"; 
      info() << endmsg;
      photon cphoton {
	.wavelength = wavelength,
	.ix = ix,
	.iy = iy,
	.layer = layer_id,
	.time = hit.getTime(),
	.photon_type = -44,
	.hitidx = i
      };
      totalPhotons.push_back(cphoton); 
      
    }

  }

  std::vector<double> xs(1024);
  std::vector<double> scintSignal(1024);
  std::vector<double> cherenkovSignal(1024);
  double dt = 0.2; // sampling time in ns 
  for (size_t i= 0; i< 1024; i++) {
    xs[i] = dt*i;
  }

  std::map<key,std::vector<double>> totalWaveforms;
  std::map<key,std::vector<double>> CherenWaveforms;
  std::map<key,std::vector<double>> ScintWaveforms;
  Filter_Type ftype = Filter_Type::NONE;
  if (m_U330_Filter.value() && m_O58_Filter.value()) {
    info() << "Error! Both filters set active! Using none" << endmsg;
    ftype = Filter_Type::NONE; 
  }
  else { 
    if (m_U330_Filter.value()) { 
      ftype = Filter_Type::U330;
    }
    if (m_O58_Filter.value()) {
      ftype = Filter_Type::O58;
    }
  }


  switch (ftype) {
  case (Filter_Type::NONE): {
    info() << "Using no filter on crystal" << endmsg;
    break;
  }
  case (Filter_Type::O58): {
    info() << "Using O58 filter" << endmsg;
    break;
  }
  case (Filter_Type::U330): {
    info() << "Using U330 Filter on crystal" << endmsg;
    break;
  }
  }


  auto storeKilledHit = [&](photon &p) {
    auto hit = simHits->at(p.hitidx);
    edm4hep::MutableCalorimeterHit calHit;
    calHit = hit.clone();

      if (p.photon_type == -44) {
	// Cherenkov Photon
	killedCherenkovPhts->push_back(calHit); 
      }
      else if (p.photon_type == -22) {
	// Scintillation Photon
	killedScintPhts->push_back(calHit); 
      }
  };

  auto storePassedHit = [&](photon &p) {
    auto hit = simHits->at(p.hitidx);
    edm4hep::MutableCalorimeterHit calHit;
    calHit = hit.clone();
    if (p.photon_type == -44) {
      // Cherenkov Photon
      passedCherenkovPhts->push_back(calHit); 
    }
    else if (p.photon_type == -22) {
      // Scintillation Photon
	passedScintPhts->push_back(calHit); 
    }
    
  };

  
  for (auto &p : totalPhotons) {

    if ((p.wavelength <= 300.0) || (p.wavelength >= 1000.0)) {
      info() << "Skipping photons <= 300 nm or >= 1000 nm" << endmsg;
      storeKilledHit(p); 
      continue;
    }
    
    double response;
    double filterresponse;
    switch (ftype) {
    case calvision::Filter_Type::U330: { 
      //filterresponse = u330_filter.Eval(p.wavelength);
      //std::lock_guard<std::mutex> lg(calvision::guard);
      filterresponse = calvision::u330_filter.Eval(p.wavelength);
      break;
    }
    case calvision::Filter_Type::O58: {
      //filterresponse = o58_filter.Eval(p.wavelength);
      //std::lock_guard<std::mutex> lg(calvision::guard); 
      filterresponse = calvision::o58_filter.Eval(p.wavelength);
      break; 
    }
    case calvision::Filter_Type::NONE: {
      filterresponse = 0.;
      break; 
    }
    }

    // Note: Add Choice for UV vs RGB
    //response = rgb_sipm_filter.Eval(p.wavelength);
    {
      //std::lock_guard<std::mutex> lg(calvision::guard);
      response = calvision::rgb_sipm_filter.Eval(p.wavelength);
    }
    
    // filter cut
    double  randval =  m_rndmUniform.shoot()*100;
    if (ftype != Filter_Type::NONE) { 
      if (randval > filterresponse ) {
	storeKilledHit(p); 
	info() << "Skipping " << p.wavelength << " nm photon.";
	info() << " filter resp: " << filterresponse << " randval " << randval << endmsg; 
	//info() << "Skipping this photon, random val > filter resp" << endmsg;
	continue;
      }
    }
    randval =  m_rndmUniform.shoot();
    if (randval > response ) {
      storeKilledHit(p); 
      info() << "Skipping " << p.wavelength << " nm photon.";
      info() << " sipm resp: " << response << " randval: " << randval << endmsg; 
	//info() << "Skipping this photon, random val > sipm response" << endmsg;
      continue;
    }

    storePassedHit(p); 
    
    int idx = int(round(p.time/dt));
    key k{};
    k.ix = p.ix;
    k.iy = p.iy;
    k.layer = p.layer; 
    if (totalWaveforms.count(k)== 0) {
      auto &vec = totalWaveforms[k];
      vec.resize(1024);
    }
    
    std::vector<double> &wave = totalWaveforms[k];
    std::vector<double> &pwave = p.photon_type == -22 ? ScintWaveforms[k]  :
      CherenWaveforms[k]; 
    if (pwave.size() != 1024)
      pwave.resize(1024);
    
    if (wave.size() != 1024) // first time we get here 
      wave.resize(1024);
    for (; idx < 1024; idx++) {
      double offset = xs[idx]-p.time;
      // in case our rounding puts us in a higher bin 
      if (offset < 0.0)
	offset = 0; 
      wave[idx] += DESY_SPR(offset);
      pwave[idx] += DESY_SPR(offset); 
    
    }
  }


  //Create and store the waveform created by Cherenkov photons, scintillation photons
  // and their combination 
  


  auto fillWaveform = [&](edm4hep::MutableTimeSeries &ts,
			  std::vector<double> &wv, int ix, int iy, int layer) {
    ts.setInterval(dt);
    ts.setTime(0);
    ts.setCellID((ix<<3)|(iy<<10)|(layer<<20));

    for (size_t i = 0; i < 1024; i++) {
      ts.addToAmplitude(wv[i]);
    }

  }; 


  // store the group of waveforms based on ix, iy positions
  for (auto &[k,v] : totalWaveforms) {
    int ix = k.ix;
    int iy = k.iy;
    int layer = k.layer;
    auto wv = waveforms->create();

    fillWaveform(wv, v, ix, iy, layer); 

  }

  
  for (auto &[k,v] : ScintWaveforms) {
    int ix = k.ix;
    int iy = k.iy;
    int layer = k.layer; 
    
    auto wv = scintwaveforms->create();

    fillWaveform(wv, v, ix, iy,layer); 
  }

  
  for (auto &[k,v] : CherenWaveforms) {

    int ix = k.ix;
    int iy = k.iy;
    int layer = k.layer; 
    auto wv = cherenwaveforms->create();

    fillWaveform(wv, v, ix, iy,layer); 
  }

  
  
  return StatusCode::SUCCESS;

}

StatusCode DualCrysSiPMAlgo::finalize()
{

  return Gaudi::Algorithm::finalize();

}
  
 
  
