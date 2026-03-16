#include "DualCrysElectronicsAlgo.h"


#include <GaudiKernel/MsgStream.h>
#include <GaudiKernel/StatusCode.h>
#include <cstdlib>
#include <edm4hep/CalorimeterHit.h>
#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/MutableCalorimeterHit.h>
#include <edm4hep/TimeSeriesCollection.h>
#include <fstream>
#include <memory>
#include <sstream>
#include <filesystem>


DECLARE_COMPONENT(DualCrysElectronicsAlgo)



DualCrysElectronicsAlgo::DualCrysElectronicsAlgo(const std::string& name, ISvcLocator* svclocator)
: Gaudi::Algorithm(name, svclocator)
{




}

// Spice Callbacks

int DualCrysElectronicsAlgo::send_char(char* str, int len, void* user) {
    DualCrysElectronicsAlgo* data = get_callback_data(user);
    data->info() << "send_char:" << str << std::endl; 
  
    //std::cout << "send_char:" << data->name << std::endl; 
    // Use your data here
    return 0;
}

int DualCrysElectronicsAlgo::send_stat(char* str, int len, void* user) {
    DualCrysElectronicsAlgo* data = get_callback_data(user);
    //std::cout << "send_stat:" << data->name << std::endl; 
    // Use your data here
    return 0;
}

int DualCrysElectronicsAlgo::controlled_exit(int exit_status, bool immediate_exit, bool on_quit, int ident, void* user) {

  DualCrysElectronicsAlgo* data = get_callback_data(user);
  //  std::cout << "controlled_exit:" << data->name << std::endl; 
  // Use your data here
  return 0;
  

}


int DualCrysElectronicsAlgo::send_data(pvecvaluesall data, int num_vectors, int ident, void* user) {

    DualCrysElectronicsAlgo* cb = get_callback_data(user);

    for (int i = 0; i < num_vectors; i++) {
      auto *pv = data->vecsa[i];

      std::string brname(pv->name);
      
      //cb->info() << pv->name << "," << data->veccount << ",";
      //cb->info() << data->vecindex << ","; 
      if (brname.find("vout") != std::string::npos) {
	if (pv->is_complex)
	  cb->info() << ",Complex,"; 
	//cb->info() << pv->creal << std::endl;
	cb->waveform.push_back(pv->creal); 
      }
      else if (brname.find("time") != std::string::npos) { 
	// time
	//cb->info() << "," << pv->creal << ",";
	//outvals->xs.push_back(pv->creal); 
      }
    }

    //std::cout << "send_data:" << cb->name << std::endl; 
    // Use your data here
    return 0;

}


int DualCrysElectronicsAlgo::send_init_data(pvecinfoall init_data, int ident, void* user) {
    DualCrysElectronicsAlgo* cb = get_callback_data(user);
    //    std::cout << "send_init_data:" << cb->name << std::endl; 
    // Use your data here
    return 0;



  
}

int DualCrysElectronicsAlgo::bg_thread_running(bool running, int ident, void* user) {
    DualCrysElectronicsAlgo* cb = get_callback_data(user);
    //std::cout << "bg_thread_running:" << cb->name << std::endl; 
    // Use your data here
    return 0;


}






StatusCode DualCrysElectronicsAlgo::initialize()
{

  StatusCode sc = Gaudi::Algorithm::initialize();

  if (sc.isFailure())
    return sc;

  // Initialize random services
  m_randSvc = service("RndmGenSvc", false);

  if (!m_randSvc) {
    error() << "Couldn't get RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_rndmUniform.initialize(m_randSvc, Rndm::Flat(0., 1.)).isFailure()) {
    error() << "Couldn't initialize RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }

  


  // try to init ngspice
  int r = ngSpice_Init(send_char,
		       send_stat,
		       controlled_exit,
		       send_data,
		       send_init_data,
		       bg_thread_running,
		       this); 

  if (r != 0) {
    error() << "Failed to init ngspice!" << std::endl; 
    return StatusCode::FAILURE;
  }

  std::string loadCmd = "source " + m_spiceFile.toString();
  info() << "Trying to load circuit file:" << loadCmd << endmsg;

  r = ngSpice_Command((char*) loadCmd.c_str());


  if (r != 0) {
    error() << "Failed to load circuit!" << std::endl; 
    return StatusCode::FAILURE;
  }
  
  
  info() << "Dual Crystal SiPM Electronics Algorithm Initialized" << endmsg;

  return StatusCode::SUCCESS;

}



StatusCode DualCrysElectronicsAlgo::execute(const EventContext&) const
{
  edm4hep::TimeSeriesCollection* cherenkovWaveforms;
  edm4hep::TimeSeriesCollection* scintillationWaveforms; 

  info() << "Starting Electron Sim" << std::endl;

  cherenkovWaveforms = m_cherenwaveforms.createAndPut();
  scintillationWaveforms = m_scintwaveforms.createAndPut();

  //  Input 
  const edm4hep::TimeSeriesCollection *cherenkovTS = m_in_cherenwaveforms.get();
  const edm4hep::TimeSeriesCollection *scintillationTS = m_in_scintwaveforms.get();

  info() << "Cherenkov Size:" << cherenkovTS->size() << std::endl;
  info() << "Scintillation Size:" << scintillationTS->size() << std::endl;




  for (size_t i = 0; i < cherenkovTS->size(); i++)  {
    std::stringstream s; 
    const auto &ts = cherenkovTS->at(i);
    const auto cellID = ts.getCellID();
    int slice_id = ((0x7<<17&cellID)>>17);
    int layer_id = ((0x7<<20&cellID)>>20);
    int ix = (0x7f<<3&cellID)>>3;
    int iy = (0x7f<<10&cellID)>>10;
    if (ix > 63)
      ix -= 64;
    if (iy > 63)
      iy -= 64;
    double dt = ts.getInterval();
    
    //    s << "alter @V1[pwl] =  [";

    for (size_t ypos = 0; ypos < ts.amplitude_size(); ypos++) {
	s << ypos*dt << "e-09" << "," << ts.getAmplitude(ypos) << std::endl;
      //s << dt*ypos << "ns," << ts.getAmplitude(ypos) << ","; 
    }
    //    s << "]";

    
    std::string pwl = s.str();
    //info() << pwl << endmsg;
    auto tmpPath = std::filesystem::temp_directory_path();
    std::stringstream fname;
    //    double rnd = round(50*m_rndmUniform.shoot());
    fname << "waveform" << ".pwl";
    tmpPath /= fname.str();
    info() << "Writing to temporary " << tmpPath << std::endl;
    //info() << pwl << endmsg; 
    std::fstream fout(tmpPath, std::ios_base::trunc | std::ios_base::out);
    fout << pwl;
    fout.close();
    

    std::vector<std::string> cmdlst;
    std::string loadCmd = "source " + m_spiceFile.toString(); 
    cmdlst.push_back(loadCmd.c_str()); // have to reload the source file to get updated data 
    cmdlst.push_back("tran 200ps 205ns 0 ");
    cmdlst.push_back("run");
    cmdlst.push_back("destroy all");
    cmdlst.push_back("reset"); 
    //cmdlst.push_back("write spice4qucs.tr1.plot v(out) ");
    for (auto &cmd : cmdlst) {
      info() << "Sending cmd:" << cmd << endmsg;  
      ngSpice_Command((char*) cmd.c_str()); 

    }
    auto wvf = cherenkovWaveforms->create();
    wvf.setInterval(ts.getInterval());
    wvf.setCellID(ts.getCellID()); 
    wvf.setTime(0);
    for (size_t pt = 0; pt < waveform.size(); pt++) {
      wvf.addToAmplitude(waveform[pt]);
    }

    waveform.clear(); 
    //ngSpice_Command((char*)");
    
    //if (i == 2)
    //break; 
  }

  for (size_t i = 0; i < scintillationTS->size(); i++)  {
    std::stringstream s; 
    const auto &ts = scintillationTS->at(i);
    const auto cellID = ts.getCellID();
    int slice_id = ((0x7<<17&cellID)>>17);
    int layer_id = ((0x7<<20&cellID)>>20);
    int ix = (0x7f<<3&cellID)>>3;
    int iy = (0x7f<<10&cellID)>>10;
    if (ix > 63)
      ix -= 64;
    if (iy > 63)
      iy -= 64;
    double dt = ts.getInterval();
    
    //    s << "alter @V1[pwl] =  [";

    for (size_t ypos = 0; ypos < ts.amplitude_size(); ypos++) {
	s << ypos*dt << "e-09" << "," << ts.getAmplitude(ypos) << std::endl;
      //s << dt*ypos << "ns," << ts.getAmplitude(ypos) << ","; 
    }
    //    s << "]";

    
    std::string pwl = s.str();
    //info() << pwl << endmsg;
    auto tmpPath = std::filesystem::temp_directory_path();
    std::stringstream fname;
    //    double rnd = round(50*m_rndmUniform.shoot());
    fname << "waveform" << ".pwl";
    tmpPath /= fname.str();
    info() << "Writing to temporary " << tmpPath << std::endl;
    //info() << pwl << endmsg; 
    std::fstream fout(tmpPath, std::ios_base::trunc | std::ios_base::out);
    fout << pwl;
    fout.close();
    

    std::vector<std::string> cmdlst;
    std::string loadCmd = "source " + m_spiceFile.toString(); 
    cmdlst.push_back(loadCmd.c_str()); // have to reload the source file to get updated data 
    cmdlst.push_back("tran 200ps 205ns 0 ");
    cmdlst.push_back("run");
    cmdlst.push_back("destroy all");
    cmdlst.push_back("reset"); 
    //cmdlst.push_back("write spice4qucs.tr1.plot v(out) ");
    for (auto &cmd : cmdlst) {
      info() << "Sending cmd:" << cmd << endmsg;  
      ngSpice_Command((char*) cmd.c_str()); 

    }
    auto wvf = scintillationWaveforms->create();
    wvf.setInterval(ts.getInterval());
    wvf.setCellID(ts.getCellID()); 
    wvf.setTime(0);
    for (size_t pt = 0; pt < waveform.size(); pt++) {
      wvf.addToAmplitude(waveform[pt]);
    }

    waveform.clear(); 
    //ngSpice_Command((char*)");
    
    //    if (i == 2)
    //break; 
  }

  return StatusCode::SUCCESS;

  
}

StatusCode DualCrysElectronicsAlgo::finalize()
{

  return Gaudi::Algorithm::finalize();

}
