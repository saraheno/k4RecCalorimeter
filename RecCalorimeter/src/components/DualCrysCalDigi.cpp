/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <cctype>
#include <cstdlib>  // abs
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "GaudiKernel/MsgStream.h"
#include "edm4hep/CalorimeterHit.h"
#include "DualCrysCalorimeterHit.h"
#include "edm4hep/Constants.h"
#include "DualCrysCalDigi.h"

DECLARE_COMPONENT(DualCrysCalDigi)

DualCrysCalDigi::DualCrysCalDigi(const std::string& aName, ISvcLocator* aSvcLoc)
    : MultiTransformer(aName, aSvcLoc,
                       {
                           KeyValues("CALCollection", {"ECalEcalCollection"}),
                           KeyValues("HeaderName", {"EventHeader"}),
                       },
                       {
			 KeyValues("FrontPhotonCollection", {"FrontPhotons"}),
			 KeyValues("RearPhotonCollection", {"RearPhotons"}),
			 KeyValues("OtherPhotonCollection", {"OtherPhotons"}),
						 
		       }
		       ) {
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
  }

  m_geoSvc = serviceLocator()->service("GeoSvc");  // important to initialize m_geoSvc
}

StatusCode DualCrysCalDigi::initialize() {
  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::CalorimeterHitCollection,
	   edm4hep::CalorimeterHitCollection,
	   edm4hep::CalorimeterHitCollection> 
DualCrysCalDigi::operator()(const edm4hep::SimCalorimeterHitCollection& SimCaloHits, const edm4hep::EventHeaderCollection& headers) const {
  debug() << " process event : " << headers[0].getEventNumber() << " - run  " << headers[0].getRunNumber()
          << endmsg;  // headers[0].getRunNumber(),headers[0].getEventNumber()


  auto frontHits = edm4hep::CalorimeterHitCollection();
  auto rearHits = edm4hep::CalorimeterHitCollection();
  auto otherHits    = edm4hep::CalorimeterHitCollection();
  
  std::string initString;

  std::string colName    = m_calCollections;
  CHT::Layout caloLayout = layoutFromString(colName);

  // hack for now, hardcode cell id
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


  debug() << "Hit Count:" << SimCaloHits.size() << endmsg; 
  for (auto &[k, v] : constants) {
    debug() << "Detector Constant:" << k << endmsg;
  }

  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(cellid_definition); 


  for (const auto& hit : SimCaloHits)
    {
      const int cellID = hit.getCellID();
      
      //Get the layer number
      unsigned int layer = bitFieldCoder.get(cellID, "layer");
      //Check if we want to use this later, else go to the next hit
      if (!useLayer(caloLayout, layer))
	continue;

      int slice_id = ((0x7<<17&cellID)>>17);
      int layer_id = ((0x7<<20&cellID)>>20);

      // cut based on cell ID
      // we're (currently) looking for layer 0, slice 1
      // or layer 1, slice 4
      bool first_pd = (slice_id == 4) && (layer_id == 1);
      bool second_pd = (slice_id == 1) && (layer_id == 0);

      auto *collection = &frontHits;
      std::vector<std::string> cols = { "rear", "front", "other" };
      size_t colsel = 0; 
      if (first_pd) { 
	collection = &rearHits;
	colsel = 0;
      }
      else if (second_pd) { 
	collection = &frontHits;
	colsel = 1;
      }
      else {
	colsel = 2; 
	collection = &otherHits;
      }
      
      if (hit.isAvailable()) {
	debug() << "Cell ID " << cellID << " collection " << cols[colsel] << endmsg; 
	debug() << "Hit Slice ID " << slice_id << " ,Layer ID " << layer_id << endmsg;

	for (auto step = hit.contributions_begin(); 
	     step != hit.contributions_end(); step++) {
	  edm4hep::CaloHitContribution contrib = *step;
	  if (contrib.isAvailable()) {
	    auto pdg = contrib.getPDG();
	    // scint or cherenkov photon, maybe redundant as might get checked upstream 
	    if (pdg == -22 || pdg == -44) { 
	      edm4hep::MutableCalorimeterHit calHit = collection->create();
	      calHit.setCellID(cellID);
	      calHit.setEnergy(contrib.getEnergy());
	      calHit.setTime(contrib.getTime());
	      calHit.setPosition(hit.getPosition());
	      calHit.setType(contrib.getPDG());
	    }
	    
	    //debug() << contrib.getPDG() << " time:";
	    //debug() << contrib.getTime() << " ns." << endmsg;
	  }
	  else {
	    debug() << "Contrib not available" << endmsg;
	  }

	}
      }
    }
  return std::make_tuple(std::move(frontHits),std::move(rearHits),std::move(otherHits));
}

//StatusCode DualCrysCalDigi::finalize() { return StatusCode::SUCCESS; }

//If the vectors are empty, we are keeping everything
bool DualCrysCalDigi::useLayer(CHT::Layout caloLayout, unsigned int layer) const {
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

// Placeholder
CHT::Layout layoutFromString(const std::string& name) {
  
  return CHT::any; 

}
