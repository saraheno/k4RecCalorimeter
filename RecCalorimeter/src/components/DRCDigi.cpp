#include "DRCDigi.h"



#include <Gaudi/PluginServiceV2.h>
#include <GaudiKernel/ISvcLocator.h>

#include <Math/Interpolator.h>
#include <cctype>
#include <cstdint>
#include <cstdlib>  // abs
#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/TimeSeriesCollection.h>
#include <exception>
#include <k4FWCore/Transformer.h>
#include <memory>
#include <mutex>
#include <sipm/SiPMAnalogSignal.h>
#include <sipm/SiPMProperties.h>
#include <sipm/SiPMSensor.h>
#include <stdexcept>
#include <tuple>
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "GaudiKernel/MsgStream.h"
#include "edm4hep/CalorimeterHit.h"
#include "DualCrysCalorimeterHit.h"
#include "edm4hep/Constants.h"

#include "DualCrysSiPMConstants.h"


using namespace calvision; 

DECLARE_COMPONENT(DRCDigi)






DRCDigi::DRCDigi(const std::string &name, ISvcLocator *svcLoc)
: MultiTransformer(name, svcLoc, 
		   {
		     KeyValues("CALCollection", {"EcalHitCollection"}),
		     KeyValues("HeaderName", {"EventHeader"}),
		   },
		   {
		     KeyValues("ScintillationSiPMWaveforms", {"ScintWaveforms"}),
		     KeyValues("CherenkovSiPMWaveforms", {"CherenkovWaveforms"}),
		     KeyValues("passedCherenkovHits", {"passedCherenkovHits"}),
		     KeyValues("passedScintillationHits", {"passedScintillationHits"}),
		     KeyValues("killedCherenkovHits", {"killedCherenkovHits"}),
		     KeyValues("killedScintillationHits", {"killedScintillationHits"})
		   })
		    {
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
  }

  m_geoSvc = serviceLocator()->service("GeoSvc");  // important to initialize m_geoSvc
}



StatusCode DRCDigi::initialize() {

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

  try {
    sipmAlgo = sipmAlgoMap.at(m_SiPMAlgorithm.toString());
  }
  catch(const std::out_of_range& ex) {
    error() << m_SiPMAlgorithm.toString() << " not found in algorithms" << endmsg;
    error() << "Using FNAL2023 as default." << endmsg;
    sipmAlgo = SiPM_Algorithm::FNAL2023; 
  }
  
  if (sipmAlgo == SiPM_Algorithm::SIM_SIPM)
    {
      // setup sipm properties

      //sipmProp.setSignalLength(204.8);
      //    sipmProp.setSampling(0.2);
      if (m_samplerate > 0.) 
	sipmProp.setSampling(m_samplerate);
      else
	sipmProp.setSampling(0.2);
      
      sipmProp.setSignalLength(sipmProp.sampling()*m_samples);
      sipmProp.setSize(m_sipmSize);
      sipmProp.setDcr(m_Dcr);
      sipmProp.setXt(m_Xt);
      sipmProp.setSampling(m_samplerate);
      sipmProp.setRecoveryTime(m_recovery);
      sipmProp.setPitch(m_cellPitch);
      sipmProp.setAp(m_afterpulse);
      sipmProp.setFallTimeFast(m_falltimeFast);
      sipmProp.setRiseTime(m_risetime);
      sipmProp.setSnr(m_snr);
      // Set the PDE type to spectrum PDE
      sipmProp.setPdeType(sipm::SiPMProperties::PdeType::kSpectrumPde);
      // Set the PDE spectrum
      sipmProp.setPdeSpectrum(m_wavelen.value(), m_sipmEff.value());
    }

  if (!calvision::filterInit) 
    calvision::init_filters(); 

  try {
    sipmAlgo = sipmAlgoMap.at(m_SiPMAlgorithm.toString());
  }
  catch(const std::out_of_range& ex) {
    error() << m_SiPMAlgorithm.toString() << " not found in algorithms" << endmsg;
    error() << "Using FNAL2023 as default." << endmsg;
    sipmAlgo = SiPM_Algorithm::FNAL2023; 
  }

  
  try { 
    sipmType = sipmTypeMap.at(m_sipmType.toString());
  }
  catch (const std::out_of_range  &e) {
    error() << m_sipmType.toString() << ":" << e.what() << endmsg;
    error() << "Using RGB as default" << endmsg;
    sipmType = calvision::SiPM_Type::RGB; 
  }

  try {
    crystal_filter = filterTypeMap.at(m_filter_type.toString());
  }
  catch ( const std::out_of_range &e) {
    error() << "Bad choice for crystal filter:" << e.what();
    error() << "Applying No crystal filter as default" << endmsg;
    crystal_filter = calvision::Filter_Type::NONE;
  }

  
  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::TimeSeriesCollection, 
	   edm4hep::TimeSeriesCollection,
	   edm4hep::CalorimeterHitCollection,
	   edm4hep::CalorimeterHitCollection,
	   edm4hep::CalorimeterHitCollection,
	   edm4hep::CalorimeterHitCollection>
DRCDigi::operator()(const edm4hep::SimCalorimeterHitCollection &simCaloHits,
				   const edm4hep::EventHeaderCollection& headers) const {

  
  debug() << " process event : " << headers[0].getEventNumber() << " - run  " << headers[0].getRunNumber()
          << endmsg;  // headers[0].getRunNumber(),headers[0].getEventNumber()

  edm4hep::TimeSeriesCollection cherenkovWaveforms;
  edm4hep::TimeSeriesCollection scintillationWaveforms;

  edm4hep::CalorimeterHitCollection passedCherenkovHits;
  edm4hep::CalorimeterHitCollection passedScintillationHits;
  edm4hep::CalorimeterHitCollection killedCherenkovHits;
  edm4hep::CalorimeterHitCollection killedScintillationHits;
  

  

  // Pull the cellID Description from the geometry if available
  // If not, pull from the bitField property
  auto *detector = m_geoSvc->getDetector();
  auto &constants = detector->constants();
  std::string cellid_definition = m_bitField; 
  if (constants.empty()) {
    debug() << "No Constants....did you load the geometry?" << endmsg; 
  }
  else {

    auto readout = detector->readout(m_detectorNameEcal);
    info() << "Readout name:" << readout.name() << endmsg;
    auto cidDesc = readout.segmentation().segmentation()->fieldDescription(); 
    info() << "Readout Segmentation Field Description:" << cidDesc << endmsg; 
    cellid_definition = cidDesc; 

  }
  
  debug() << "Hit Count:" << simCaloHits.size() << endmsg;
  switch (sipmAlgo) { 
  case SiPM_Algorithm::SIM_SIPM: {
    debug() << "Using SIM_SIPM for Waveform construction" << endmsg;
    break;
  }
  case SiPM_Algorithm::DESY: {
    debug() << "Using DESY for Waveform construction" << endmsg;
    break;
  }
  case calvision::SiPM_Algorithm::FNAL2023: {
    debug() << "Using FNAL 2023 SPR for Waveform construction" << endmsg;
    break; 
  }  
  };

  
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(cellid_definition);  


  std::vector<photon> sPhotons;
  std::vector<photon> cPhotons; 
  std::vector<std::function<bool (const photon &p)>> filterFns;

  ROOT::Math::Interpolator *sipm_response = &calvision::rgb_sipm_filter;


  switch(sipmType) {
  case calvision::SiPM_Type::RGB: {
    sipm_response = &calvision::rgb_sipm_filter;
    break; 
  }
  case calvision::SiPM_Type::UV: {
    sipm_response = &calvision::uv_sipm_filter;
    break;
  }
  case calvision::SiPM_Type::BROADCOM: {
    sipm_response = &calvision::broadcom_2x1_sipm_filter;
    break;
  }
  };

  

  
  auto wvlfilter = [](const photon &p) {
    return p.wavelength <= 300.0 || p.wavelength >= 1000.0;
  };

  //  std::function<bool (const photon &p)> respFilter;

  
  auto respfilter = [&](const photon &p) {
      double response; 
      //response = rgb_sipm_filter.Eval(p.wavelength);
      //std::lock_guard<std::mutex> lg(calvision::guard); 
      response = sipm_response->Eval(p.wavelength);
      if (std::isnan(response)) {
	info() << p.wavelength << " nm is past SiPM's response curve." << std::endl;
	response = 0.0;
      }

      double  randval = m_rndmUniform.shoot();
      if (randval > response)
	return true;
      return false;
  }; 


  filterFns.push_back(wvlfilter);

  // Only use our extra response filtering with the DESY / FNAL Algorithm
  if (sipmAlgo == SiPM_Algorithm::DESY || sipmAlgo == SiPM_Algorithm::FNAL2023) 
    filterFns.push_back(respfilter);

  
  auto processPhoton = [](const std::vector<photon> &ps,
			   std::vector<photon> &livePhotons,
			   edm4hep::CalorimeterHitCollection &pass,
			   edm4hep::CalorimeterHitCollection &kill,
			   std::vector<std::function<bool (const photon &p)>> &filters
			   ) {

    for (const photon &p : ps) {
      bool passPhoton = true;
      for (auto &f : filters) {
	if (f(p))
	  passPhoton = false; 
      }
      
      
      if (!passPhoton) { 
	auto hit = kill.create();
	hit.setCellID(p.cellid);
	hit.setEnergy(p.energy);
      }
      else {
	auto hit = pass.create();
	hit.setCellID(p.cellid);
	hit.setEnergy(p.energy); 
	livePhotons.push_back(p);
      }

    }
  }; 
  
  
  for (const auto& hit : simCaloHits)
    {

      auto [sphts, cphts] = processHit(hit);

      processPhoton(sphts, sPhotons, passedScintillationHits, 
		    killedScintillationHits,filterFns); 
      processPhoton(cphts, cPhotons, 
		    passedCherenkovHits, killedCherenkovHits,filterFns); 

    }

  // Now we have all the photons for this event
  info() << "Found " << sPhotons.size() << " Scintillation photons." << endmsg;
  info() << "Found " << cPhotons.size() << " Cherenkov photons." << endmsg;


 

  std::vector<double> xs(m_samples);
  std::vector<double> scintSignal(m_samples);
  std::vector<double> cherenkovSignal(m_samples);

  double dt;
  if (m_samplerate < 0. || m_samplerate == 0.) 
    dt = 0.2; 
  else
    dt = m_samplerate; // sampling time in ns 




  
  if (sipmAlgo == SiPM_Algorithm::DESY) { 

    std::map<key,std::vector<double>> CherenWaveforms;
    std::map<key,std::vector<double>> ScintWaveforms;
  
  
    for (size_t i= 0; i< m_samples; i++) {
      xs[i] = dt*i;
    }

    auto fillWaveform = [](edm4hep::MutableTimeSeries &ts,
			   std::vector<double> &wv, 
			   int ix, int iy, int layer,
			   double sampleRate,
			   size_t samples) {
      ts.setInterval(sampleRate);
      ts.setTime(0);
      ts.setCellID((ix<<3)|(iy<<10)|(layer<<20));

      for (size_t i = 0; i < samples; i++) {
	ts.addToAmplitude(wv[i]);
      }
    }; 



    // create waveform positions 
    generate_waveform_positions(sPhotons, xs, ScintWaveforms, dt, m_samples);
    generate_waveform_positions(cPhotons, xs, CherenWaveforms, dt, m_samples);

    auto filler = [&](std::map<key,std::vector<double>> &waves,
		      edm4hep::TimeSeriesCollection &waveCollection)
    {
      for (auto &[k,v] : waves) {
      int ix = k.ix;
      int iy = k.iy;
      int layer = k.layer; 
      auto wv = waveCollection.create();
      fillWaveform(wv, v, ix, iy,layer,dt,m_samples);
      }
    };
      
    filler(ScintWaveforms, scintillationWaveforms);
    filler(CherenWaveforms, cherenkovWaveforms); 
  
  }

  if (sipmAlgo == SiPM_Algorithm::SIM_SIPM) {

    std::map<key,sipm::SiPMSensor> CherenWaveforms;
    std::map<key,sipm::SiPMSensor> ScintWaveforms;

    

    generate_waveform_positions(sPhotons, ScintWaveforms, sipmProp);
    generate_waveform_positions(cPhotons, CherenWaveforms, sipmProp);

    auto fillWaveform = []( key k, sipm::SiPMSensor &sensor,
			     edm4hep::MutableTimeSeries &ts)
    {

      ts.setInterval(sensor.properties().sampling());
      ts.setTime(0);
      uint64_t cID = (k.ix << 3)|(k.iy << 10)|(k.layer<<20);
      ts.setCellID(cID);
      sensor.runEvent();
      sipm::SiPMAnalogSignal signal = sensor.signal();

      for (auto pt: signal.waveform()) {
	ts.addToAmplitude(pt);
      }

    }; 

    for (auto &[k,sensor] : ScintWaveforms) {
      auto wvfrm = scintillationWaveforms.create();
      fillWaveform(k, ScintWaveforms[k], wvfrm);
    }
    for (auto &[k,sensor] : CherenWaveforms) {
      auto wvfrm = cherenkovWaveforms.create();
      fillWaveform(k, CherenWaveforms[k], wvfrm);
    }

  }
  
  return std::make_tuple(std::move(cherenkovWaveforms), 
			 std::move(scintillationWaveforms),
			 std::move(passedCherenkovHits),
			 std::move(passedScintillationHits),
			 std::move(killedCherenkovHits),
			 std::move(killedScintillationHits)
			 ); 
  
}


bool DRCDigi::useLayer(CHT::Layout caloLayout, unsigned int layer) const {
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


typedef std::vector<photon> PhotonVector;
std::tuple<PhotonVector, PhotonVector> 
DRCDigi::processHit(const edm4hep::SimCalorimeterHit &hit) const {

  PhotonVector sPhotons;
  PhotonVector cPhotons; 
  const uint64_t cellID = hit.getCellID();

  unpackedcellID cid = unpackCellID(cellID);
  

  
  bool first_pd = (cid.slice == 4) && (cid.layer == 1);
  bool second_pd = (cid.slice == 1) && (cid.layer == 0);
  if (first_pd || second_pd) {
    if (hit.isAvailable()) {
      for (auto step = hit.contributions_begin();
	   step != hit.contributions_end(); step++) {
	
	float energy = step->getEnergy()/CLHEP::eV;
	double wavelength = 1239.84187 / (1000*energy); 
	photon p {
	  .wavelength = wavelength,
	  .energy = energy,
	  .ix = cid.ix,
	  .iy = cid.iy,
	  .layer = cid.layer,
	  .time = step->getTime(),
	  .photon_type = step->getPDG(),
	  .cellid = cellID
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


void DRCDigi::generate_waveform_positions(const std::vector<photon> &photons,
							 std::map<key,sipm::SiPMSensor> &waveMap,
							 const sipm::SiPMProperties &sipmprops) const
{
      for (const auto &p : photons ) {

	key k{};
	k.ix = p.ix;
	k.iy = p.iy;
	k.layer = p.layer;
	if (waveMap.count(k) == 0) {
	  waveMap[k] = sipm::SiPMSensor(sipmprops);
	}
	auto &sensor = waveMap[k];
	sensor.addPhoton(p.time, p.wavelength);
      }
}



void DRCDigi::generate_waveform_positions(const std::vector<photon> &photons,
							 const std::vector<double> &timev,
							 std::map<key,std::vector<double>> &waveMap,
							 double sampleRate,
							 size_t samples) const
{

  std::function<double(double)> sprFn = DESY_SPR;

  if (sipmAlgo == SiPM_Algorithm::FNAL2023)
    sprFn = FNAL2023_SPR;
  
  for (const auto &p : photons) {
    size_t idx = int(round(p.time/sampleRate));
    key k{.ix = p.ix,
	  .iy = p.iy,
	  .layer = p.layer
    }; 
    
    if (waveMap.count(k) == 0) {
      auto &vec = waveMap[k];
      vec.resize(samples);
    }
    auto &pvec = waveMap[k]; 
    
    for (; idx < samples; idx++) {
      double offset = timev[idx]-p.time;
      if (offset < 0.0)
	offset = 0; 
      pvec[idx] += sprFn(offset); 
    }
  }
  
}
