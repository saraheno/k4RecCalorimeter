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


/** @class DualCrysSiPMAlgo

    Algorithm for generating waveforms from a template fitting method developed
    by Sasha Ledovsky

    @author Thomas Anderson
    @date 2025-10-07


*/



class DualCrysSiPMAlgo : public Gaudi::Algorithm {
 public:
  DualCrysSiPMAlgo(const std::string& name, ISvcLocator* svcLoc);
  virtual ~DualCrysSiPMAlgo() {};

  StatusCode initialize() override;
  StatusCode execute(const EventContext&) const override;
  StatusCode finalize() override;

  double SPR(double now) const; // FNAL Electronics result to single photoelectron 

  
 private:

  // Template Properties
  Gaudi::Property<double> m_Rise{this, "Rise", 0.853, "SiPM SPR Rise Time in ns"};
  Gaudi::Property<double> m_Decay{this, "Decay", 6.538, "SiPM SPR Decay Time in ns"};
  Gaudi::Property<double> m_UnderShoot{this, "Undershoot", 101.7, "SPR Undershoot"};
  Gaudi::Property<double> m_norm{this, "norm", 0.111051, "Unsure what this does..."};

  Gaudi::Property<std::string> m_hitCollection{this, "inputHitCollection", "CalorimeterHit",
      "Input containing the collection of photon hits"};
  Gaudi::Property<std::string> m_RelCollection{this, "RelationInputCollection", "RelationCalHit",
      "Input containing the relational hit collection"};
  
  Gaudi::Property<std::string> m_outTimeColl{this, "outputTimeStructCollection", "CalvisionSiPMDigiWaveform",
                                             "calvision waveform collection name"};
  Gaudi::Property<std::string> m_outScintTimeColl{this, "scintoutputTimeStructCollection", "CalvisionSiPMScintWaveform",
                                             "calvision scint waveform collection name"};
  Gaudi::Property<std::string> m_outCerenTimeColl{this, "cerenoutputTimeStructCollection", "CalvisionSiPMCerenWaveform",
                                             "calvision ceren waveform collection name"};
    
  // Input Collections
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_simHits{m_hitCollection,
      Gaudi::DataHandle::Reader,
      this};

  mutable k4FWCore::DataHandle<edm4hep::CaloHitSimCaloHitLinkCollection> m_links{m_RelCollection,
      Gaudi::DataHandle::Reader,
      this};

  mutable k4FWCore::DataHandle<edm4hep::TimeSeriesCollection> m_waveforms{m_outTimeColl,
      Gaudi::DataHandle::Writer,
      this};
  
  mutable k4FWCore::DataHandle<edm4hep::TimeSeriesCollection> m_cerenwaveforms{m_outCerenTimeColl,
      Gaudi::DataHandle::Writer,
      this};

  mutable k4FWCore::DataHandle<edm4hep::TimeSeriesCollection> m_scintwaveforms{m_outScintTimeColl,
      Gaudi::DataHandle::Writer,
      this};
  


};
    
