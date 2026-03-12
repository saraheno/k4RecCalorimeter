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
#include <GaudiKernel/DataHandle.h>


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
    int layer; 
    double time;
    int photon_type;
    size_t hitidx; 
  };


  enum class Filter_Type {
    NONE,
    U330,
    O58
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
const std::set<int> u330_filter_wavelengths = { 230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500,510,520,530,540,550,560,570,580,590,600,610,620,630,640,650,660,670,680,690,700,710,720,730,740,750,760,770,780,790,800,810,820,830,840,850,860,870,880,890,900,910,920,930,940,950,960,970,980,990,1000,1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,1110,1120,1130,1140,1150,1160,1170,1180,1190,1200};
  const std::map<int, double > u330_filterMap = {{ 230, 0.0 },{ 240, 15.0 },{ 250, 62.0 },{ 260, 76.0 },{ 270, 83.0 },{ 280, 86.0 },{ 290, 87.0 },{ 300, 87.5 },{ 310, 87.8 },{ 320, 87.8 },{ 330, 87.8 },{ 340, 87.5 },{ 350, 86.5 },{ 360, 85.0 },{ 370, 80.0 },{ 380, 65.0 },{ 390, 35.0 },{ 400, 10.0 },{ 410, 1.0 },{ 420, 0.0 },{ 430, 0.0 },{ 440, 0.0 },{ 450, 0.0 },{ 460, 0.0 },{ 470, 0.0 },{ 480, 0.0 },{ 490, 0.0 },{ 500, 0.0 },{ 510, 0.0 },{ 520, 0.0 },{ 530, 0.0 },{ 540, 0.0 },{ 550, 0.0 },{ 560, 0.0 },{ 570, 0.0 },{ 580, 0.0 },{ 590, 0.0 },{ 600, 0.0 },{ 610, 0.0 },{ 620, 0.0 },{ 630, 0.0 },{ 640, 0.0 },{ 650, 1.0 },{ 660, 3.0 },{ 670, 10.0 },{ 680, 22.0 },{ 690, 38.0 },{ 700, 46.0 },{ 710, 49.0 },{ 720, 49.5 },{ 730, 47.0 },{ 740, 43.0 },{ 750, 39.0 },{ 760, 35.0 },{ 770, 32.0 },{ 780, 29.0 },{ 790, 27.0 },{ 800, 25.0 },{ 810, 24.0 },{ 820, 23.0 },{ 830, 22.5 },{ 840, 22.0 },{ 850, 22.0 },{ 860, 22.0 },{ 870, 22.2 },{ 880, 22.5 },{ 890, 23.0 },{ 900, 23.8 },{ 910, 24.5 },{ 920, 25.2 },{ 930, 26.0 },{ 940, 27.0 },{ 950, 28.0 },{ 960, 29.0 },{ 970, 30.0 },{ 980, 31.0 },{ 990, 31.8 },{ 1000, 32.5 },{ 1010, 33.0 },{ 1020, 33.5 },{ 1030, 33.8 },{ 1040, 33.8 },{ 1050, 33.5 },{ 1060, 33.0 },{ 1070, 32.2 },{ 1080, 31.5 },{ 1090, 30.5 },{ 1100, 29.0 },{ 1110, 27.5 },{ 1120, 26.0 },{ 1130, 24.5 },{ 1140, 23.0 },{ 1150, 21.5 },{ 1160, 20.0 },{ 1170, 18.5 },{ 1180, 17.0 },{ 1190, 16.0 },{ 1200, 15.0 }};

  
  


const std::set<int> o58_filter_wavelengths = { 200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,500,510,520,530,540,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,641,651,661,671,681,691,701,711,721,731,741,751,761,771,781,791,801,811,821,831,841,851,861,871,881,891,901,911,921,931,941,951,961,971,981,991,1001,1011,1021,1031,1041,1051,1061,1071,1081,1091,1101,1111,1121,1131,1141,1151,1161,1171,1181,1191};
const std::map<int, double > o58_filterMap = {{ 200, 0.0 },{ 210, 0.0 },{ 220, 0.0 },{ 230, 0.0 },{ 240, 0.0 },{ 250, 0.0 },{ 260, 0.0 },{ 270, 0.0 },{ 280, 0.0 },{ 290, 0.0 },{ 300, 0.0 },{ 310, 0.0 },{ 320, 0.0 },{ 330, 0.0 },{ 340, 0.0 },{ 350, 0.0 },{ 360, 0.0 },{ 370, 0.0 },{ 380, 0.0 },{ 390, 0.0 },{ 400, 0.0 },{ 410, 0.0 },{ 420, 0.0 },{ 430, 0.0 },{ 440, 0.0 },{ 450, 0.0 },{ 460, 0.0 },{ 470, 0.0 },{ 480, 0.0 },{ 490, 0.0 },{ 500, 0.0 },{ 510, 0.0 },{ 520, 0.0 },{ 530, 0.0 },{ 540, 0.0 },{ 550, 0.0 },{ 551, 0.1 },{ 552, 0.2 },{ 553, 0.3 },{ 554, 0.5 },{ 555, 0.8 },{ 556, 1.2 },{ 557, 1.8 },{ 558, 2.5 },{ 559, 3.5 },{ 560, 4.8 },{ 561, 6.5 },{ 562, 8.5 },{ 563, 11.0 },{ 564, 14.0 },{ 565, 17.5 },{ 566, 21.5 },{ 567, 26.0 },{ 568, 31.0 },{ 569, 36.5 },{ 570, 42.0 },{ 571, 47.5 },{ 572, 53.0 },{ 573, 58.0 },{ 574, 62.5 },{ 575, 66.5 },{ 576, 70.0 },{ 577, 73.0 },{ 578, 75.5 },{ 579, 77.8 },{ 580, 79.5 },{ 581, 81.0 },{ 582, 82.2 },{ 583, 83.2 },{ 584, 84.1 },{ 585, 84.8 },{ 586, 85.5 },{ 587, 86.0 },{ 588, 86.5 },{ 589, 86.9 },{ 590, 87.2 },{ 591, 87.5 },{ 592, 87.8 },{ 593, 88.0 },{ 594, 88.2 },{ 595, 88.4 },{ 596, 88.6 },{ 597, 88.8 },{ 598, 88.9 },{ 599, 89.0 },{ 600, 89.1 },{ 601, 89.2 },{ 602, 89.3 },{ 603, 89.4 },{ 604, 89.5 },{ 605, 89.6 },{ 606, 89.7 },{ 607, 89.8 },{ 608, 89.8 },{ 609, 89.9 },{ 610, 90.0 },{ 611, 90.0 },{ 612, 90.1 },{ 613, 90.1 },{ 614, 90.2 },{ 615, 90.2 },{ 616, 90.2 },{ 617, 90.3 },{ 618, 90.3 },{ 619, 90.3 },{ 620, 90.4 },{ 621, 90.4 },{ 622, 90.4 },{ 623, 90.4 },{ 624, 90.5 },{ 625, 90.5 },{ 626, 90.5 },{ 627, 90.5 },{ 628, 90.5 },{ 629, 90.5 },{ 630, 90.6 },{ 631, 90.6 },{ 641, 90.60070298769772 },{ 651, 90.60140597539542 },{ 661, 90.60210896309314 },{ 671, 90.60281195079085 },{ 681, 90.60351493848857 },{ 691, 90.60421792618628 },{ 701, 90.604920913884 },{ 711, 90.60562390158172 },{ 721, 90.60632688927943 },{ 731, 90.60702987697715 },{ 741, 90.60773286467486 },{ 751, 90.60843585237258 },{ 761, 90.60913884007029 },{ 771, 90.60984182776801 },{ 781, 90.61054481546573 },{ 791, 90.61124780316344 },{ 801, 90.61195079086116 },{ 811, 90.61265377855887 },{ 821, 90.61335676625659 },{ 831, 90.6140597539543 },{ 841, 90.61476274165202 },{ 851, 90.61546572934974 },{ 861, 90.61616871704744 },{ 871, 90.61687170474517 },{ 881, 90.61757469244287 },{ 891, 90.6182776801406 },{ 901, 90.6189806678383 },{ 911, 90.61968365553602 },{ 921, 90.62038664323374 },{ 931, 90.62108963093145 },{ 941, 90.62179261862917 },{ 951, 90.62249560632688 },{ 961, 90.6231985940246 },{ 971, 90.62390158172231 },{ 981, 90.62460456942003 },{ 991, 90.62530755711775 },{ 1001, 90.62601054481546 },{ 1011, 90.62671353251318 },{ 1021, 90.62741652021089 },{ 1031, 90.62811950790861 },{ 1041, 90.62882249560631 },{ 1051, 90.62952548330404 },{ 1061, 90.63022847100176 },{ 1071, 90.63093145869946 },{ 1081, 90.63163444639719 },{ 1091, 90.6323374340949 },{ 1101, 90.63304042179261 },{ 1111, 90.63374340949034 },{ 1121, 90.63444639718804 },{ 1131, 90.63514938488576 },{ 1141, 90.63585237258347 },{ 1151, 90.63655536028119 },{ 1161, 90.6372583479789 },{ 1171, 90.63796133567662 },{ 1181, 90.63866432337434 },{ 1191, 90.63936731107205 }};
  
  const std::set<int> UV_Wavelengths = {361,365,380,388,397,406,412,427,437,456,477,496,518,548,573,599,615,649,671,705,756,774,798};

  const std::set<int> RGB_Wavelengths = {305,318,335,352,370,396,416,444,467,477,492,515,530,557,583,610,636,663,684,712,739,756,775,795,826,850,874,895,901};

  const std::map<int, std::pair<double,double> > UV_Map = { { 361, {361.161,0.770121} },{ 365, {364.766,0.787349} },{ 380, {379.794,0.879305} },{ 388, {387.614,0.94252} },{ 397, {396.624,0.982752} },{ 406, {406.226,1} },{ 412, {411.617,0.982752} },{ 427, {426.594,0.94252} },{ 437, {436.769,0.890797} },{ 456, {455.931,0.816089} },{ 477, {477.492,0.741381} },{ 496, {496.061,0.683901} },{ 518, {517.627,0.620686} },{ 548, {547.583,0.545978} },{ 573, {573.349,0.488498} },{ 599, {598.521,0.448266} },{ 615, {615.299,0.41379} },{ 649, {649.46,0.35633} },{ 671, {671.039,0.327591} },{ 705, {705.202,0.275867} },{ 756, {755.548,0.201139} },{ 774, {773.531,0.178155} },{ 798, {798.108,0.149416} }};

  const std::map<int, std::pair<double,double> > RGB_Map = { { 305, {305.28,0.0346782} },{ 318, {318.47,0.144499} },{ 335, {334.67,0.271679} },{ 352, {352.06,0.42775} },{ 370, {370.06,0.525999} },{ 396, {396.44,0.635839} },{ 416, {416.23,0.705196} },{ 444, {443.81,0.786125} },{ 467, {466.6,0.878607} },{ 477, {477.39,0.907518} },{ 492, {491.78,0.93641} },{ 515, {515.17,0.994214} },{ 530, {529.56,1} },{ 557, {556.53,0.976875} },{ 583, {582.91,0.942196} },{ 610, {610.49,0.901732} },{ 636, {636.26,0.849715} },{ 663, {663.24,0.780338} },{ 684, {684.22,0.734107} },{ 712, {712.39,0.664731} },{ 739, {738.76,0.583802} },{ 756, {755.55,0.520232} },{ 775, {774.73,0.485554} },{ 795, {795.11,0.42775} },{ 826, {825.68,0.364161} },{ 850, {850.26,0.289018} },{ 874, {874.23,0.225428} },{ 895, {894.61,0.167624} },{ 901, {900.61,0.144499} }};
  
  const std::pair<double,double> findnearest(SiPM_Type stype, double wavelength) const;
  const std::pair<double,double> findnearest(Filter_Type ftype, double wavelength) const;

  
 private:

  // Template Properties
  Gaudi::Property<double> m_Rise{this, "Rise", 0.853, "SiPM SPR Rise Time in ns"};
  Gaudi::Property<double> m_Decay{this, "Decay", 6.538, "SiPM SPR Decay Time in ns"};
  Gaudi::Property<double> m_UnderShoot{this, "Undershoot", 101.7, "SPR Undershoot"};
  Gaudi::Property<double> m_norm{this, "norm", 0.111051, "Unsure what this does..."};
  Gaudi::Property<bool> m_U330_Filter{this, "U330", false, "Use U330 Crystal Filter, default no"};
  Gaudi::Property<bool> m_O58_Filter{this, "O58", false, "Use O58 Crystal Filter, default no"};
  Gaudi::Property<std::string> m_hitCollection{this, "inputHitCollection", "CalorimeterHit",
      "Input containing the collection of photon hits"};

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

  
  Gaudi::Property<std::string> m_outTimeColl{this, "outputTimeStructCollection", "CalvisionSiPMDigiWaveform",
                                             "calvision waveform collection name"};
  Gaudi::Property<std::string> m_outScintTimeColl{this, "scintoutputTimeStructCollection", "CalvisionSiPMScintWaveform",
                                             "calvision scint waveform collection name"};
  Gaudi::Property<std::string> m_outCherenTimeColl{this, "cherenoutputTimeStructCollection", "CalvisionSiPMCherenWaveform",
                                             "calvision cheren waveform collection name"};


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
    
