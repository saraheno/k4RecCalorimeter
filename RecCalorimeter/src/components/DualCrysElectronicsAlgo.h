#pragma once


#include <ngspice/sharedspice.h>



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

#include <memory>

/** @class DualCrysSiPMAlgo

    Algorithm for adding electronics modeling to SiPM waveforms using ngspice

    @author Thomas Anderson
    @date 2026-03-12


*/


class DualCrysElectronicsAlgo : public Gaudi::Algorithm {
 public:

  DualCrysElectronicsAlgo(const std::string &name, ISvcLocator* svcLoc);
  virtual ~DualCrysElectronicsAlgo() {};

  StatusCode initialize() override;
  StatusCode execute(const EventContext&) const override;
  StatusCode finalize() override;



  // spice callbacks, I may move this to a separate struct in the future

  static int send_char(char* str, int len, void* user);
  static int send_stat(char* str, int len, void* user);
  static int controlled_exit(int exit_status, bool immediate_exit, bool on_quit, int ident, void* user);
  static int send_data(pvecvaluesall data, int num_vectors, int ident, void* user);
  static int send_init_data(pvecinfoall init_data, int ident, void* user);
  static int bg_thread_running(bool running, int ident, void* user);
  // Helper function to cast userData to your class type
  static DualCrysElectronicsAlgo* get_callback_data(void* user) {
        return static_cast<DualCrysElectronicsAlgo*>(user);
    }

  
 private:
  Gaudi::Property<std::string> m_outScintElecColl{this, "scintoutputElecTimeStructCollection", 
      "ScintCircuitWaveform",
      "calvision scint electronics waveform collection name"};

  Gaudi::Property<std::string> m_outCherenElecColl{this, "cherenoutputTimeElecStructCollection", 
      "CherenkovCircuitWaveform",
      "calvision ceren electronics waveform collection name"};

    Gaudi::Property<std::string> m_inScintTimeColl{this, "ScintSiPMInputWaveforms", 
						 "ScintWaveforms",
						 "calvision scint waveform collection name"};
  Gaudi::Property<std::string> m_inCherenTimeColl{this, 
						  "CherenkovSiPMInputWaveforms",
						  "CherenkovWaveforms",
						  "calvision ceren waveform collection name"};
  Gaudi::Property<std::string> m_SiPMDigiTimeColl{this, 
						  "CalvisionSiPMDigiWaveform",
						  "CalvisionSiPMDigiWaveform",
						  "calvision combo waveform collection name"};

  Gaudi::Property<std::string> m_spiceFile{this, 
						  "spicefile",
						  "50ohm.cir",
						  "Spice circuit to use"};


						  /*						  

  Gaudi::Property<std::string> m_killedScint{this, 
						  "killedScintPhotons",
						  "killedScintPhotons",
						  "killed scint photons"};

  Gaudi::Property<std::string> m_passedScint{this, 
						  "passedScintPhotons",
						  "passedScintPhotons",
						  "passed scint photons"};

  Gaudi::Property<std::string> m_killedCheren{this, 
						  "killedCherenPhotons",
						  "killedCherenPhotons",
						  "killed cheren photons"};

  Gaudi::Property<std::string> m_passedCheren{this, 
						  "passedCherenPhotons",
						  "passedCherenPhotons",
						  "passed cheren photons"};
  */

  // Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;
  Rndm::Numbers m_rndmUniform;

  mutable k4FWCore::DataHandle<edm4hep::TimeSeriesCollection> m_in_cherenwaveforms{m_inCherenTimeColl,
      Gaudi::DataHandle::Reader,
      this};

  mutable k4FWCore::DataHandle<edm4hep::TimeSeriesCollection> m_in_scintwaveforms{m_inScintTimeColl,
      Gaudi::DataHandle::Reader,
      this};

  mutable k4FWCore::DataHandle<edm4hep::TimeSeriesCollection> m_in_combo{"CalvisionSiPMDigiWaveform",
      Gaudi::DataHandle::Reader,
      this};

  
  mutable std::vector<double> waveform; 
  

  // out
  mutable k4FWCore::DataHandle<edm4hep::TimeSeriesCollection> m_cherenwaveforms{m_outCherenElecColl,
      Gaudi::DataHandle::Writer,
      this};

  mutable k4FWCore::DataHandle<edm4hep::TimeSeriesCollection> m_scintwaveforms{m_outScintElecColl,
      Gaudi::DataHandle::Writer,
      this};
  

}; 
