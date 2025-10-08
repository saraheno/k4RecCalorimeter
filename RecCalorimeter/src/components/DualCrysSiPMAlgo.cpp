#include "DualCrysSiPMAlgo.h"
#include <CLHEP/Units/SystemOfUnits.h>

DECLARE_COMPONENT(DualCrysSiPMAlgo)

DualCrysSiPMAlgo::DualCrysSiPMAlgo(const std::string& name, ISvcLocator* svclocator)
: Gaudi::Algorithm(name, svclocator)
{




}



/* Find the nearest wavelength from our fixed set of wavelengths
 This works by rounding the wavelength off to the nearest integer
 then finding the lowest value greater than given rounded wavelength
 and then checking the nearest neighbors for the closest value and returning
 the pair of nearest wavelength and PDE for a given SiPM type
*/
const std::pair<double,double> DualCrysSiPMAlgo::findnearest(SiPM_Type stype, double wavelength) const {
  //info() << "FindNearest Wavelength:" << wavelength << endmsg; 
  if (wavelength < 0.)
    return std::make_pair(0.0,0.0); 

  auto lb = UV_Wavelengths.lower_bound(round(wavelength)); 
  if (stype == SiPM_Type::RGB) {
    lb = RGB_Wavelengths.lower_bound(round(wavelength));
  }


  
  int start = *lb;
  int delta = abs(wavelength - *lb);
  lb--;
  int lowerdelta = abs(wavelength - *lb);
  if (lowerdelta < delta) {
    start = *lb;
    delta = lowerdelta;
  }
  lb++;
  lb++;
  int upperdelta = abs(wavelength - *lb);
  if (upperdelta < delta) {
    start = *lb;
  }

  if (stype == SiPM_Type::RGB)
    return RGB_Map.at(start);
  else if (stype == SiPM_Type::UV)
    return UV_Map.at(start); 
  else
    return std::make_pair(0.0,0.0); 

}


// pulled straight from Resolution.C, need to modify it to use the class parameters

double DualCrysSiPMAlgo::SPR(double now) const
{

  //  double tMin_  = 0.0;
  //  double tMax_  = 1000.0;

  double tRise       = 0.853;
  double tDecay      = 6.538;
  double tUnderShoot = 101.7;
  //  double norm        = 0.111051;


  double a = 1./ tRise;
  double b = 1./ tDecay;
  double A = -a * b / (a - b);
  double B = -A;
  double result = A * exp(-a*now) + B * exp(-b*now);
    
  double g = 1./ tUnderShoot;
  double Atmp = -A * g / ( a - g);
  double Btmp = -B * g / ( b - g);
  double G = - Atmp - Btmp ;
  A = Atmp;
  B = Btmp;
  result -= A * exp(-a*now) + B * exp(-b*now) + G * exp(-g*now);
    
  return result;

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

  
  info() << "Dual Crystal SiPM Algorithm Initialized" << endmsg; 

  return StatusCode::SUCCESS;

}


// This builds and returns waveforms for a given set of hits 
StatusCode DualCrysSiPMAlgo::execute(const EventContext&) const
{
  
  const edm4hep::CalorimeterHitCollection* simHits = m_simHits.get();
  const edm4hep::CaloHitSimCaloHitLinkCollection* linkCollection = m_links.get();

  // Output
  edm4hep::TimeSeriesCollection* waveforms = m_waveforms.createAndPut();
  edm4hep::TimeSeriesCollection* cerenwaveforms = m_cerenwaveforms.createAndPut();
  edm4hep::TimeSeriesCollection* scintwaveforms = m_scintwaveforms.createAndPut();
  
  info() << "Sim Hit Size:" << simHits->size() << " :: ";
  info() << "Link Size:" << linkCollection->size() << endmsg; 
  std::vector<double> scintPhotons;
  std::vector<double> scintPhotonWavelength;
  std::vector<double> cerenkPhotons;
  std::vector<double> cerenkPhotonWavelength;
  for (size_t i = 0; i < simHits->size(); i++) {

    const auto &hit = simHits->at(i); 
    const int cellID = hit.getCellID();
    const int slice_id = ((0x7<<17&cellID)>>17);
    const int layer_id = ((0x7<<20&cellID)>>20);
    // not using the slice or layer at the moment,
    // but the upstream feeder should only give us particles from
    // the SiPM tower slice
    // 1.2398 / energy = wavelength hc/e ~ 1.239e-6 eV*m
    // I think the value we get is in GeV, have to check the calculation
    // but it looks ~ right in the results but needs verifying / rewriting
    if (hit.getType() == -22) {
      // scint
      scintPhotons.push_back(hit.getTime());
      float energy = hit.getEnergy()/CLHEP::eV;

      scintPhotonWavelength.push_back(1239.84187 / (1000*energy));
      info() << "Scint @ " << hit.getTime() << " Wavelength:" << 1239.84187/(1000*energy);
      info() << " Energy: " << energy << " eV?"; 
      info() << endmsg; 
    }
    else if ( hit.getType() == -44) {
      // cerenkov
      cerenkPhotons.push_back(hit.getTime());
      float energy = hit.getEnergy()/CLHEP::eV;
      cerenkPhotonWavelength.push_back(1239.84187 / (1000*energy));
      info() << "Cerenk @ " << hit.getTime() << " Wavelength:" << 1239.84187/(1000*energy);
      info() << " Energy: " << energy << " eV?"; 
      info() << endmsg; 
    }

  }

  std::vector<double> xs(1024);
  std::vector<double> scintSignal(1024);
  std::vector<double> cerenkovSignal(1024);
  double dt = 0.2; // sampling time in ns 
  for (size_t i= 0; i< 1024; i++) {
    xs[i] = dt*i;
  }

  // Use RGB SiPM for now 
  for (size_t i = 0; i < scintPhotons.size(); i++) {
    double phot = scintPhotons[i];
    double wvl,response;
    std::tie(wvl,response) = findnearest(SiPM_Type::RGB, scintPhotonWavelength[i]); 
    double randval =  m_rndmUniform.shoot();
    info() << "Rand value: " << randval << " WaveLen:" << wvl;
    info() << " PDE:" << response << endmsg; 
    if (randval > response ) {
      info() << "Skipping this scint. photon, random val > resp" << endmsg;
      continue;
    }
    int idx = int(round(phot/dt));
    for (; idx < 1024; idx++) {
      double offset = xs[idx]-phot;
      // in case our rounding puts us in a higher bin 
      if (offset < 0.0)
	offset = 0; 
	
      scintSignal[idx] += SPR(offset); 
    }
  }
  for (size_t i = 0; i < cerenkPhotons.size(); i++) {
    double phot = cerenkPhotons[i];
    double wvl,response;
    std::tie(wvl,response) = findnearest(SiPM_Type::RGB, cerenkPhotonWavelength[i]); 
    double randval =  m_rndmUniform.shoot();
    info() << "Rand value: " << randval << " WaveLen:" << wvl;
    info() << " PDE:" << response << endmsg; 
    if (randval > response ) {
      info() << "Skipping this cerenk. photon, random val > resp" << endmsg;
      continue;
    }
    int idx = int(round(phot/dt));
    for (; idx < 1024; idx++) {
      double offset = xs[idx]-phot;
      // in case our rounding puts us in a higher bin 
      if (offset < 0.0)
	offset = 0; 
	
      cerenkovSignal[idx] += SPR(offset); 
    }
  }

  //Create and store the waveform created by Cerenkov photons, scintillation photons
  // and their combination 
  
  info() << "Filling scint waveform with " << scintPhotons.size() << " photons.";
  info() << endmsg; 
  auto scintwaveform = scintwaveforms->create();
  scintwaveform.setInterval(dt);
  scintwaveform.setTime(0);
  scintwaveform.setCellID(simHits->at(0).getCellID());

  for (double sample : scintSignal) {
    scintwaveform.addToAmplitude(sample);
  }

  info() << "Filling cerenkov waveform with " << cerenkPhotons.size() << " photons.";
  info() << endmsg; 

  auto cerenkwaveform = cerenwaveforms->create();
  cerenkwaveform.setInterval(dt);
  cerenkwaveform.setTime(0);
  cerenkwaveform.setCellID(simHits->at(0).getCellID());

  for (double sample : cerenkovSignal) {
    cerenkwaveform.addToAmplitude(sample);
  }

  auto waveform = waveforms->create();
  waveform.setInterval(dt);
  waveform.setTime(0);
  waveform.setCellID(simHits->at(0).getCellID());
  for (size_t i = 0; i < 1024; i++) {
    waveform.addToAmplitude(cerenkovSignal[i]+scintSignal[i]);
  }
    
  info() << "Scintillation Photons:" << scintPhotons.size() << endmsg;
  info() << "Cerenkov Photons:" << cerenkPhotons.size() << endmsg;
  
  return StatusCode::SUCCESS;

}

StatusCode DualCrysSiPMAlgo::finalize()
{

  return Gaudi::Algorithm::finalize();

}
  
 
  
