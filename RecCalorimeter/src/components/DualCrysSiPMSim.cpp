#include "DualCrysSiPMSim.h"
#include <CLHEP/Units/SystemOfUnits.h>
#include <GaudiKernel/MsgStream.h>
#include <edm4hep/MutableTimeSeries.h>
#include <sipm/SiPMAnalogSignal.h>
#include <sipm/SiPMProperties.h>
#include <sipm/SiPMSensor.h>
#include <algorithm>
#include <edm4hep/CalorimeterHit.h>
#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/MutableCalorimeterHit.h>
#include <memory>

DECLARE_COMPONENT(DualCrysSiPMSim)

DualCrysSiPMSim::DualCrysSiPMSim(const std::string& name, ISvcLocator* svclocator)
: Gaudi::Algorithm(name, svclocator)
{




}

struct key {
  int ix;
  int iy;
  int layer;
  };

// Defined in DualCrysSiPMAlgo, I'll change this around
// eventually so it has a local version 
bool operator<(const key &lhs, const key &rhs);

    


StatusCode DualCrysSiPMSim::initialize()
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


  // setup sipm properties
  sipmProp.setSignalLength(204.8);
  sipmProp.setSampling(0.2);

  //sipmSensor = std::make_unique<sipm::SiPMSensor>(sipmProp); 


   
  // 
  
  info() << "Dual Crystal SiPM Algorithm Initialized" << endmsg; 

  return StatusCode::SUCCESS;

}

// This builds and returns waveforms for a given set of hits 
StatusCode DualCrysSiPMSim::execute(const EventContext&) const
{
  
  const edm4hep::CalorimeterHitCollection* simHits = m_simHits.get();
  //const edm4hep::CaloHitSimCaloHitLinkCollection* linkCollection = m_links.get();

  
  // Output
  edm4hep::TimeSeriesCollection* totalCollection = m_waveforms.createAndPut();
  edm4hep::TimeSeriesCollection* cherenCollection = m_cherenwaveforms.createAndPut();
  edm4hep::TimeSeriesCollection* scintCollection = m_scintwaveforms.createAndPut();


  edm4hep::CalorimeterHitCollection* killedCherenkovPhts = m_killedCherenPhotons.createAndPut();
  edm4hep::CalorimeterHitCollection* killedScintPhts = m_killedScintPhotons.createAndPut();

  edm4hep::CalorimeterHitCollection* passedScintPhts = m_passedScintPhotons.createAndPut();

  edm4hep::CalorimeterHitCollection* passedCherenkovPhts = 
    m_passedCherenPhotons.createAndPut(); 

  
  info() << "Sim Hit Size:" << simHits->size() << " :: ";

  std::vector<photon> totalPhotons;
  // find and separate out Cherenkov and Scintillation Photons 
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
      // cerenkov
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


  
  std::map<key,sipm::SiPMSensor> totalWaveforms;
  std::map<key,sipm::SiPMSensor> CherenWaveforms;
  std::map<key,sipm::SiPMSensor> ScintWaveforms;

  size_t idx = 0; 
  for (auto &p : totalPhotons) {

    if ((p.wavelength <= 300.0) || (p.wavelength >= 1000.0)) {
      info() << "Skipping photons <= 300 nm or >= 1000 nm" << endmsg;
      storeKilledHit(p); 
      continue;
    }
    // Filter cuts and SiPM Cuts go here 



    storePassedHit(p); 
    
    key k{};
    k.ix = p.ix;
    k.iy = p.iy;
    k.layer = p.layer;

    if (totalWaveforms.count(k) == 0) {
      totalWaveforms[k] = sipm::SiPMSensor(sipmProp);
      totalWaveforms[k].resetState(); 
    }
    if (p.photon_type == -22) {
      if (ScintWaveforms.count(k) == 0) {
	ScintWaveforms[k] =  sipm::SiPMSensor(sipmProp);
      }
    }
    else if (p.photon_type == -44) {
      if (CherenWaveforms.count(k) == 0) {
	CherenWaveforms[k] = sipm::SiPMSensor(sipmProp);
      }
    }
    auto &sensor = totalWaveforms[k];
    auto &tsensor = p.photon_type == -22 ? ScintWaveforms[k] : CherenWaveforms[k];

    sensor.addPhoton(p.time);
    tsensor.addPhoton(p.time); 
    if ((idx %10) == 0) {
      info() << "Processing photon " << idx << endmsg;
    }
    idx++; 
  }

  idx = 0; 


  auto fillWaveform = [&](  key k, sipm::SiPMSensor &sensor, edm4hep::MutableTimeSeries &waveform) {

      waveform.setInterval(sensor.properties().sampling()); 
      waveform.setTime(0);
      uint64_t cID = (k.ix << 3)|(k.iy << 10)|(k.layer<<20); 
      waveform.setCellID(cID);
      sensor.runEvent();
      sipm::SiPMAnalogSignal signal = sensor.signal();
    
      for (auto pt : signal.waveform()) {
	waveform.addToAmplitude(pt);
      }
    }; 


  for (auto &[k,sensor] : ScintWaveforms) {
    auto wvfrm = scintCollection->create();
    fillWaveform(k, sensor, wvfrm);
  }
  for (auto &[k,sensor] : CherenWaveforms) {
    auto wvfrm = cherenCollection->create();
    fillWaveform(k, sensor, wvfrm);
  }
  for (auto &[k,sensor] : totalWaveforms) {
    auto wvfrm = totalCollection->create();
    fillWaveform(k, sensor, wvfrm);
  }


  return StatusCode::SUCCESS;
}


StatusCode DualCrysSiPMSim::finalize()
{

  return Gaudi::Algorithm::finalize();

}
