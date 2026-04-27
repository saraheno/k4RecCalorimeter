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


#include "DualCrysSiPMConstants.h"
#include <map> 
#include <thread>
/** @class DualCrysSiPMAlgo

    Algorithm for generating waveforms from a template fitting method developed
    by Sasha Ledovsky

    @author Thomas Anderson
    @date 2025-10-07


*/


using calvision::SiPM_Algorithm; 

class DualCrysSiPMAlgo : public Gaudi::Algorithm {
 public:
  DualCrysSiPMAlgo(const std::string& name, ISvcLocator* svcLoc);
  virtual ~DualCrysSiPMAlgo() {};

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

  std::thread::id init_id; 
  
  // Template Properties
  Gaudi::Property<double> m_Rise{this, "Rise", 0.853, "SiPM SPR Rise Time in ns"};
  Gaudi::Property<double> m_Decay{this, "Decay", 6.538, "SiPM SPR Decay Time in ns"};
  Gaudi::Property<double> m_UnderShoot{this, "Undershoot", 101.7, "SPR Undershoot"};
  Gaudi::Property<double> m_norm{this, "norm", 0.111051, "Unsure what this does..."};


  calvision::Filter_Type crystal_filter; 
  Gaudi::Property<std::string> m_filter_type{this, "filter_type", "NONE",
					     "Applied Crystal Filter, u330, o58, none"}; 

  calvision::SiPM_Type sipmType; 
  Gaudi::Property<std::string> m_sipmType{this, "sipm_type","RGB",
					  "SiPM Response type, Broadcom-2x1, UV, or RGB"}; 
  Gaudi::Property<std::string> m_hitCollection{this, "inputHitCollection", "CalorimeterHit",
      "Input containing the collection of photon hits"};

  Gaudi::Property<std::string> m_killedScintCollection{this,
			       "killedScintPhotonCollection",
			       "killedScintPhotons",
			       "Output containing the collection of killed Scintillation photons"};
  Gaudi::Property<std::string> m_killedCherenCollection{this, 
				"killedCherenPhotonCollection",
				"killedCherenPhotons",
				"Output containing the collection of killed Cherenkov photons"};

  Gaudi::Property<std::string> m_passedScintCollection{this,
			       "passedScintPhotonCollection", 
			       "passedScintPhotons",
			       "Output containing the collection of passed Scintillation photons"};

  Gaudi::Property<std::string> m_passedCherenCollection{this, 
			       "passedCherenPhotonCollection",
			       "passedCherenPhotons",
			       "Output containing the collection of passed Cherenkov photons"};

  
  Gaudi::Property<std::string> m_outTimeColl{this, "outputTimeStructCollection", "CalvisionSiPMDigiWaveform",
                                             "calvision waveform collection name"};
  Gaudi::Property<std::string> m_outScintTimeColl{this, "scintoutputTimeStructCollection", "CalvisionSiPMScintWaveform",
                                             "calvision scint waveform collection name"};
  Gaudi::Property<std::string> m_outCherenTimeColl{this, "cherenoutputTimeStructCollection", "CalvisionSiPMCherenWaveform",
                                             "calvision cheren waveform collection name"};
  SiPM_Algorithm sipmAlgo; 
  Gaudi::Property<std::string> m_SiPMAlgorithm { this, "SiPMAlgorithm", "FNAL2023",
						 "The SiPM Algorithm used for waveform construction"}; 


  // Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;
  Rndm::Numbers m_rndmUniform;
    
  // Input Collections
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_simHits{m_hitCollection,
      Gaudi::DataHandle::Reader,
      this};

  /*mutable k4FWCore::DataHandle<edm4hep::CaloHitSimCaloHitLinkCollection> m_links{m_RelCollection,
      Gaudi::DataHandle::Reader,
      this};*/

  mutable k4FWCore::DataHandle<edm4hep::TimeSeriesCollection> m_waveforms{m_outTimeColl,
      Gaudi::DataHandle::Writer,
      this};
  
  mutable k4FWCore::DataHandle<edm4hep::TimeSeriesCollection> m_cherenwaveforms{m_outCherenTimeColl,
      Gaudi::DataHandle::Writer,
      this};

  mutable k4FWCore::DataHandle<edm4hep::TimeSeriesCollection> m_scintwaveforms{m_outScintTimeColl,
      Gaudi::DataHandle::Writer,
      this};
  
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
    
