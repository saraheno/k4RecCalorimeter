#include "DualCrysSiPMAlgo.h"


DECLARE_COMPONENT(DualCrysSiPMAlgo)

DualCrysSiPMAlgo::DualCrysSiPMAlgo(const std::string& name, ISvcLocator* svclocator)
: Gaudi::Algorithm(name, svclocator)
{




}


StatusCode DualCrysSiPMAlgo::initialize()
{

  StatusCode sc = Gaudi::Algorithm::initialize();

  if (sc.isFailure())
    return sc;

  info() << "Dual Crystal SiPM Algorithm Initialized" << endmsg; 

  return StatusCode::SUCCESS;

}

StatusCode DualCrysSiPMAlgo::execute(const EventContext&) const
{
  // Do nothing

  const edm4hep::SimCalorimeterHitCollection* simHits = m_simHits.get();
  const edm4hep::CaloHitSimCaloHitLinkCollection* linkCollection = m_links.get();
  info() << "Sim Hit Size:" << simHits->size() << " :: ";
  info() << "Link Size:" << linkCollection->size() << endmsg; 

  for (size_t i = 0; i < simHits->size(); i++) {

    const auto &hit = simHits->at(i); 
    const int cellID = hit.getCellID();
    const int slice_id = ((0x7<<17&cellID)>>17);
    const int layer_id = ((0x7<<20&cellID)>>20);

    std::vector<double> scintPhotons;
    std::vector<double> cerenkPhotons;

    for (auto step = hit.contributions_begin(); 
	 step != hit.contributions_end(); step++)
      {
	edm4hep::CaloHitContribution contrib = *step;
	if (contrib.isAvailable()) {
	  if (-22 == contrib.getPDG()) {
	    // scint
	    scintPhotons.push_back(contrib.getTime()); 
	  }
	  else if (-44 == contrib.getPDG()) {
	    // ceren
	    cerenkPhotons.push_back(contrib.getTime());
	  }
	}
      }

    info() << "Scintillation Photons:" << scintPhotons.size() << endmsg;
    info() << "Cerenkov Photons:" << cerenkPhotons.size() << endmsg;
    //info() << "Hit Slice ID " << slice_id << " ,Layer ID " << layer_id << endmsg;
    
    

  }
  
  return StatusCode::SUCCESS;

}

StatusCode DualCrysSiPMAlgo::finalize()
{

  return Gaudi::Algorithm::finalize();

}
  
 
  
