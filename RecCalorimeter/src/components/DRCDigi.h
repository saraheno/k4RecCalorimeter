#pragma once

#include "Gaudi/Property.h"
#include "edm4hep/CaloHitContributionCollection.h"
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include "CalorimeterHitType.h"
#include "DualCrysCalorimeterHit.h"
#include "DDRec/SurfaceManager.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"


#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"

#include <Gaudi/PluginServiceV2.h>
#include <GaudiKernel/ISvcLocator.h>
#include <edm4hep/SimCalorimeterHit.h>
#include <edm4hep/TimeSeriesCollection.h>
#include <random>
#include <string>
#include <vector>

#include "DualCrysSiPMConstants.h"

#include <sipm/SiPMProperties.h>
#include <sipm/SiPMSensor.h>


using namespace calvision; 

struct DRCDigi final :

  k4FWCore::MultiTransformer<
    std::tuple<edm4hep::TimeSeriesCollection,
	       edm4hep::TimeSeriesCollection,
	       edm4hep::CalorimeterHitCollection,
	       edm4hep::CalorimeterHitCollection,
	       edm4hep::CalorimeterHitCollection,
	       edm4hep::CalorimeterHitCollection>
    (const edm4hep::SimCalorimeterHitCollection &,
     const edm4hep::EventHeaderCollection&)> {


  
  DRCDigi(const std::string &name, ISvcLocator *svcLoc);


    StatusCode initialize() override;
  // StatusCode finalize() override;

    std::tuple<edm4hep::TimeSeriesCollection, 
	       edm4hep::TimeSeriesCollection,
      edm4hep::CalorimeterHitCollection,
      edm4hep::CalorimeterHitCollection,
      edm4hep::CalorimeterHitCollection,
	       edm4hep::CalorimeterHitCollection>

    operator()(const edm4hep::SimCalorimeterHitCollection &simCaloHits,
	       const edm4hep::EventHeaderCollection& headers) const override;


  std::tuple<std::vector<photon>, std::vector<photon>> processHit(const edm4hep::SimCalorimeterHit &hit) const; 


  // Use DESY Model
  void generate_waveform_positions(const std::vector<photon> &photons,
			     std::map<key,sipm::SiPMSensor> &waveMap,
				   const sipm::SiPMProperties &sipmprops) const;

  // Use Sim SiPM Model 
  void generate_waveform_positions(const std::vector<photon> &photons,
				   const std::vector<double> &timev,
				   std::map<key,std::vector<double>> &waveMap,
				   double sampleRate,
				   size_t samples) const;

  
  
 private:

  bool useLayer(CHT::Layout caloLayout, unsigned int layer) const; 

  SiPM_Algorithm sipmAlgo; 
  Gaudi::Property<std::string> m_SiPMAlgorithm { this, "SiPMAlgorithm", "DESY",
						    "The SiPM Algorithm used for waveform construction"}; 
  
  Gaudi::Property<std::string> m_calCollections{this, "calCollections", "DRCNoSegment",
                                                "The input collection of calorimeters"};
  Gaudi::Property<std::string> outputCalCollection{this, "outputCalCollection", "outputCalCollection",
                                                   "The output collection of calorimeters"};
  Gaudi::Property<std::string> m_bitField{this, "bitField", "system:3,ix:-7,iy:-7,slice:3,layer:3,wc1:3,wc2:3,wc3:3",
                                                "The bitfield used with the ECAL detector"};
  Gaudi::Property<size_t> m_samples{this, "samples", 1024,
				    "SiPM Waveform sample count"};
  Gaudi::Property<double> m_samplerate{this, "samplerate", 0.2,
				    "SiPM Waveform sample rate"};


  Gaudi::Property<std::string> m_encodingStringVariable{
      this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
      "The name of the DD4hep constant that contains the encoding string for tracking detectors"};


  Gaudi::Property<float> m_calibrCoeffCal{this, "calibrationCoeffCal", 120000.0, "Calibration coefficient of calorimeters"};
  Gaudi::Property<std::string> m_detectorNameEcal{this, "detectorNameEcal", "DRCNoSegment", "Name of ECAL"};
  Gaudi::Property<std::string> m_detectorNameHcal{this, "detectorNameHcal", "DRFtubeFiber", "Name of HCAL"};
  Gaudi::Property<std::vector<bool>> m_useLayersEcalVec{this, "useLayersEcal", {}, "Enable/disable ECAL layers"};
  Gaudi::Property<std::vector<bool>> m_useLayersHcalVec{this, "useLayersHcal", {}, "Enable/disable HCAL layers"};

  std::string m_collName;

  SmartIF<IGeoSvc> m_geoSvc;
  SmartIF<IUniqueIDGenSvc> m_uidSvc;
  // Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;
  Rndm::Numbers m_rndmUniform;



  // Sim SiPM Properties
  sipm::SiPMProperties sipmProp;
  

  // Properties pulled from SimulatedSiPMwithOpticalPhoton.hamamatsu 
  // SiPM properties (defaults set based on Hamamatsu S14160-1310PS)
  // Signal properties
  //  Gaudi::Property<double> m_sigLength{this, "signalLength", 200., "Signal length in ns"};
  //  Gaudi::Property<double> m_sampling{this, "sampling", 0.1, "SiPM sampling rate in ns"};
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
  Gaudi::Property<std::vector<double>> m_sipmEff{this, "sipmEfficiency", {0.1, 0.1},
						 "SiPM efficiency vs wavelength"};



};


DECLARE_COMPONENT(DRCDigi)
