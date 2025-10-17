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



  struct photon {
    double wavelength;
    int ix;
    int iy;
    double time;
    int photon_type; 
  };


  
  enum class SiPM_Type {
    UV,
    RGB
  };

  /* We may want to move the SiPM Data into a file or other source eventually.

     I took the approach of converting the wavelength doubles into ints in sets and
     use those to figure out the "nearest" wavelength from the given wavelengths.
     We'll probably want to replace this with some kind of interpolation scheme
     would be my guess, but for now this will work.

     
   */
  const std::set<int> UV_Wavelengths = {361,365,380,388,397,406,412,427,437,456,477,496,518,548,573,599,615,649,671,705,756,774,798};

  const std::set<int> RGB_Wavelengths = {305,318,335,352,370,396,416,444,467,477,492,515,530,557,583,610,636,663,684,712,739,756,775,795,826,850,874,895,901};

  const std::map<int, std::pair<double,double> > UV_Map = { { 361, {361.161,0.770121} },{ 365, {364.766,0.787349} },{ 380, {379.794,0.879305} },{ 388, {387.614,0.94252} },{ 397, {396.624,0.982752} },{ 406, {406.226,1} },{ 412, {411.617,0.982752} },{ 427, {426.594,0.94252} },{ 437, {436.769,0.890797} },{ 456, {455.931,0.816089} },{ 477, {477.492,0.741381} },{ 496, {496.061,0.683901} },{ 518, {517.627,0.620686} },{ 548, {547.583,0.545978} },{ 573, {573.349,0.488498} },{ 599, {598.521,0.448266} },{ 615, {615.299,0.41379} },{ 649, {649.46,0.35633} },{ 671, {671.039,0.327591} },{ 705, {705.202,0.275867} },{ 756, {755.548,0.201139} },{ 774, {773.531,0.178155} },{ 798, {798.108,0.149416} }};

  const std::map<int, std::pair<double,double> > RGB_Map = { { 305, {305.28,0.0346782} },{ 318, {318.47,0.144499} },{ 335, {334.67,0.271679} },{ 352, {352.06,0.42775} },{ 370, {370.06,0.525999} },{ 396, {396.44,0.635839} },{ 416, {416.23,0.705196} },{ 444, {443.81,0.786125} },{ 467, {466.6,0.878607} },{ 477, {477.39,0.907518} },{ 492, {491.78,0.93641} },{ 515, {515.17,0.994214} },{ 530, {529.56,1} },{ 557, {556.53,0.976875} },{ 583, {582.91,0.942196} },{ 610, {610.49,0.901732} },{ 636, {636.26,0.849715} },{ 663, {663.24,0.780338} },{ 684, {684.22,0.734107} },{ 712, {712.39,0.664731} },{ 739, {738.76,0.583802} },{ 756, {755.55,0.520232} },{ 775, {774.73,0.485554} },{ 795, {795.11,0.42775} },{ 826, {825.68,0.364161} },{ 850, {850.26,0.289018} },{ 874, {874.23,0.225428} },{ 895, {894.61,0.167624} },{ 901, {900.61,0.144499} }};
  
  const std::pair<double,double> findnearest(SiPM_Type stype, double wavelength) const;

  
 private:

  // Template Properties
  Gaudi::Property<double> m_Rise{this, "Rise", 0.853, "SiPM SPR Rise Time in ns"};
  Gaudi::Property<double> m_Decay{this, "Decay", 6.538, "SiPM SPR Decay Time in ns"};
  Gaudi::Property<double> m_UnderShoot{this, "Undershoot", 101.7, "SPR Undershoot"};
  Gaudi::Property<double> m_norm{this, "norm", 0.111051, "Unsure what this does..."};

  Gaudi::Property<std::string> m_hitCollection{this, "inputHitCollection", "CalorimeterHit",
      "Input containing the collection of photon hits"};
  
  Gaudi::Property<std::string> m_outTimeColl{this, "outputTimeStructCollection", "CalvisionSiPMDigiWaveform",
                                             "calvision waveform collection name"};
  Gaudi::Property<std::string> m_outScintTimeColl{this, "scintoutputTimeStructCollection", "CalvisionSiPMScintWaveform",
                                             "calvision scint waveform collection name"};
  Gaudi::Property<std::string> m_outCerenTimeColl{this, "cerenoutputTimeStructCollection", "CalvisionSiPMCerenWaveform",
                                             "calvision ceren waveform collection name"};


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
  
  mutable k4FWCore::DataHandle<edm4hep::TimeSeriesCollection> m_cerenwaveforms{m_outCerenTimeColl,
      Gaudi::DataHandle::Writer,
      this};

  mutable k4FWCore::DataHandle<edm4hep::TimeSeriesCollection> m_scintwaveforms{m_outScintTimeColl,
      Gaudi::DataHandle::Writer,
      this};
  


};
    
