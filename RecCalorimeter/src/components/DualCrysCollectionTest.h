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

#include <Gaudi/PluginServiceV2.h>
#include <GaudiKernel/ISvcLocator.h>
#include <edm4hep/SimCalorimeterHit.h>
#include <edm4hep/TimeSeriesCollection.h>
#include <random>
#include <string>
#include <vector>




struct DualCrysCollectionTest final :

  k4FWCore::MultiTransformer<
    std::tuple<edm4hep::TimeSeriesCollection,
	       edm4hep::TimeSeriesCollection>
    (const edm4hep::SimCalorimeterHitCollection &,
     const edm4hep::EventHeaderCollection&)> {



  struct photon {
    double wavelength;
    int ix;
    int iy;
    int layer; 
    double time;
    int photon_type;
  };


  
  DualCrysCollectionTest(const std::string &name, ISvcLocator *svcLoc);


    StatusCode initialize() override;
  // StatusCode finalize() override;

    std::tuple<edm4hep::TimeSeriesCollection, 
      edm4hep::TimeSeriesCollection> operator()(const edm4hep::SimCalorimeterHitCollection &simCaloHits,
						const edm4hep::EventHeaderCollection& headers) const override;


  std::tuple<std::vector<photon>, std::vector<photon>> processHit(const edm4hep::SimCalorimeterHit &hit) const; 

 private:

  bool useLayer(CHT::Layout caloLayout, unsigned int layer) const; 


  Gaudi::Property<std::string> m_calCollections{this, "calCollections", "DRCNoSegment",
                                                "The input collection of calorimeters"};
  Gaudi::Property<std::string> outputCalCollection{this, "outputCalCollection", "outputCalCollection",
                                                   "The output collection of calorimeters"};


  Gaudi::Property<std::string> m_encodingStringVariable{
      this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
      "The name of the DD4hep constant that contains the encoding string for tracking detectors"};


  Gaudi::Property<float> m_calibrCoeffCal{this, "calibrationCoeffCal", 120000.0, "Calibration coefficient of calorimeters"};
  Gaudi::Property<std::string> m_detectorNameEcal{this, "detectorNameEcal", "DRCrystal", "Name of ECAL"};
  Gaudi::Property<std::string> m_detectorNameHcal{this, "detectorNameHcal", "DRFtubeFiber", "Name of HCAL"};
  Gaudi::Property<std::vector<bool>> m_useLayersEcalVec{this, "useLayersEcal", {}, "Enable/disable ECAL layers"};
  Gaudi::Property<std::vector<bool>> m_useLayersHcalVec{this, "useLayersHcal", {}, "Enable/disable HCAL layers"};

  std::string m_collName;

  SmartIF<IGeoSvc> m_geoSvc;
  SmartIF<IUniqueIDGenSvc> m_uidSvc;

    
};


DECLARE_COMPONENT(DualCrysCollectionTest)
