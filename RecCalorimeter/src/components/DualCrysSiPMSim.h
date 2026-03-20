#pragma once



// EDM4HEP includes
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/RawTimeSeriesCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/TimeSeriesCollection.h"

// k4FWCore includes
#include "k4FWCore/DataHandle.h"

// Gaudi includes
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ToolHandle.h"
#include <Gaudi/Property.h>
#include <GaudiKernel/DataHandle.h>


#include <sipm/SiPMProperties.h>
#include <sipm/SiPMSensor.h>
#include <memory>


#include "DualCrysSiPMConstants.h"

/** @class DualCrysSiPMAlgo

    Algorithm for generating waveforms using the SimSiPM library 

    @author Thomas Anderson
    @date 2026-1-20


*/



class DualCrysSiPMSim : public Gaudi::Algorithm {
 public: 
  DualCrysSiPMSim(const std::string& name, ISvcLocator* svcLoc);
  virtual ~DualCrysSiPMSim() {};

  StatusCode initialize() override;
  StatusCode execute(const EventContext&) const override;
  StatusCode finalize() override;





  struct photon {
    double wavelength;
    int ix;
    int iy;
    int layer; 
    double time;
    int photon_type;
    size_t hitidx; 
  };

  

 private:

  // Template Properties
  sipm::SiPMProperties sipmProp;
  //std::unique_ptr<sipm::SiPMSensor> sipmSensor;

  // Properties pulled from SimulatedSiPMwithOpticalPhoton.hamamatsu 
  // SiPM properties (defaults set based on Hamamatsu S14160-1310PS)
  // Signal properties
  Gaudi::Property<double> m_sigLength{this, "signalLength", 200., "Signal length in ns"};
  Gaudi::Property<double> m_sampling{this, "sampling", 0.1, "SiPM sampling rate in ns"};
  Gaudi::Property<double> m_risetime{this, "risetime", 1., "Signal rise time in ns"};
  Gaudi::Property<double> m_falltimeFast{this, "falltimeFast", 6.5, "Signal fast component decay time in ns"};


  // SiPM physical properties
  Gaudi::Property<double> m_sipmSize{this, "SiPMsize", 1.3, "Width of photosensitive area in mm"};
  Gaudi::Property<double> m_cellPitch{this, "cellpitch", 10., "SiPM cell size in um"};
  Gaudi::Property<double> m_recovery{this, "recovery", 10., "SiPM cell recovery time in ns"};

  // Noise parameters
  Gaudi::Property<double> m_Dcr{this, "DCR", 120e3, "SiPM dark count rate in Hz"};
  Gaudi::Property<double> m_Xt{this, "Xtalk", 0.01, "SiPM optical crosstalk probability"};
  Gaudi::Property<double> m_afterpulse{this, "afterpulse", 0.03, "Afterpulse probability"};
  Gaudi::Property<double> m_snr{this, "SNR", 20., "Signal-to-noise ratio in dB"};

  // Integration parameters
  Gaudi::Property<double> m_gateStart{this, "gateStart", 5., "Integration gate starting time in ns"};
  Gaudi::Property<double> m_gateL{this, "gateLength", 95., "Integration gate length in ns"};
  Gaudi::Property<double> m_thres{this, "threshold", 1.5, "Integration threshold in photoelectrons"};

  // other parameters (attention, will override above parameters if set)
  Gaudi::Property<std::map<std::string,double>> m_params{this, "params", {}, "optional parameters"};

  // SiPM efficiency
  Gaudi::Property<std::vector<double>> m_wavelen{
      this, "wavelength", {1000., 100.}, "wavelength vector in nm (decreasing order)"};
  Gaudi::Property<std::vector<double>> m_sipmEff{this, "sipmEfficiency", {0.1, 0.1}, "SiPM efficiency vs wavelength"};



  
  Gaudi::Property<std::string> m_hitCollection{this, "inputHitCollection", "CalorimeterHit",
      "Input containing the collection of photon hits"};
  Gaudi::Property<std::string> m_outTimeColl{this, "outputTimeStructCollection", "CalvisionSiPMDigiWaveform",
                                             "calvision waveform collection name"};
  Gaudi::Property<std::string> m_outScintTimeColl{this, "scintoutputTimeStructCollection", "CalvisionSiPMScintWaveform",
                                             "calvision scint waveform collection name"};
  Gaudi::Property<std::string> m_outCherenTimeColl{this, "cherenoutputTimeStructCollection", "CalvisionSiPMCherenWaveform",
                                             "calvision ceren waveform collection name"};

  Gaudi::Property<std::string> m_killedScintCollection{this,
			       "killedScintPhotonCollection",
			       "killedScintPhotons",
			       "Output containing the collection of killed Scintillation photons"};
  Gaudi::Property<std::string> m_killedCherenCollection{this, 
				"killedCherenPhotonCollection", "killedCherenPhotons",
				"Output containing the collection of killed Cherenkov photons"};

  Gaudi::Property<std::string> m_passedScintCollection{this,
			       "passedScintPhotonCollection", 
			       "passedScintPhotons",
			       "Output containing the collection of passed Scintillation photons"};

  Gaudi::Property<std::string> m_passedCherenCollection{this, 
			       "passedCherenPhotonCollection",
			       "passedCherenPhotons",
			       "Output containing the collection of passed Cherenkov photons"};


    // Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;
  Rndm::Numbers m_rndmUniform;

  // Input Collections
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_simHits{m_hitCollection,
      Gaudi::DataHandle::Reader,
      this};

  // Output

  mutable k4FWCore::DataHandle<edm4hep::TimeSeriesCollection> m_cherenwaveforms{m_outCherenTimeColl,
      Gaudi::DataHandle::Writer,
      this};

  mutable k4FWCore::DataHandle<edm4hep::TimeSeriesCollection> m_scintwaveforms{m_outScintTimeColl,
      Gaudi::DataHandle::Writer,
      this};

  mutable k4FWCore::DataHandle<edm4hep::TimeSeriesCollection> m_waveforms{m_outTimeColl,
      Gaudi::DataHandle::Writer,
      this};

  // photon hit collections to get wavelength info from, probably replacable with something more
  // efficient in the future
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_killedScintPhotons{m_killedScintCollection,
										       Gaudi::DataHandle::Writer,
										       this};
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_killedCherenPhotons{m_killedCherenCollection,
											Gaudi::DataHandle::Writer,
											this}; 
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_passedScintPhotons{m_passedScintCollection,
										       Gaudi::DataHandle::Writer,
										       this};
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_passedCherenPhotons{m_passedCherenCollection,
											Gaudi::DataHandle::Writer,
											this}; 
  

}; 

