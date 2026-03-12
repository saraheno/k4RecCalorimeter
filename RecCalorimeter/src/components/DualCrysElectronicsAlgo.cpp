#include "DualCrysElectronicsAlgo.h"


#include <GaudiKernel/StatusCode.h>
#include <edm4hep/CalorimeterHit.h>
#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/MutableCalorimeterHit.h>
#include <memory>

DECLARE_COMPONENT(DualCrysElectronicsAlgo)



DualCrysElectronicsAlgo::DualCrysElectronicsAlgo(const std::string& name, ISvcLocator* svclocator)
: Gaudi::Algorithm(name, svclocator)
{




}

// Spice Callbacks

int DualCrysElectronicsAlgo::send_char(char* str, int len, void* user) {
    DualCrysElectronicsAlgo* data = get_callback_data(user);
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

int DualCrysElectronicsAlgo::controlled_exit(int exit_status, NG_BOOL immediate_exit, NG_BOOL on_quit, int ident, void* user) {

  DualCrysElectronicsAlgo* data = get_callback_data(user);
  //  std::cout << "controlled_exit:" << data->name << std::endl; 
  // Use your data here
  return 0;
  

}


int DualCrysElectronicsAlgo::send_data(pvecvaluesall data, int num_vectors, int ident, void* user) {

    DualCrysElectronicsAlgo* cb = get_callback_data(user);
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

int DualCrysElectronicsAlgo::bg_thread_running(NG_BOOL running, int ident, void* user) {
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

  
  info() << "Dual Crystal SiPM Electronics Algorithm Initialized" << endmsg;

  // try to init ngspice
  int r = ngSpice_Init(send_char,
		       send_stat,
		       controlled_exit,
		       send_data,
		       send_init_data,
		       bg_thread_running,
		       this); 
  if (r != 0)
    return StatusCode::FAILURE; 


  return StatusCode::SUCCESS;

}



StatusCode DualCrysElectronicsAlgo::execute(const EventContext&) const
{


  edm4hep::TimeSeriesCollection* cherenkovWaveforms = m_cherenwaveforms.createAndPut();
  edm4hep::TimeSeriesCollection* scintillationWaveforms = m_scintwaveforms.createAndPut();

  return StatusCode::SUCCESS;

  
}

StatusCode DualCrysElectronicsAlgo::finalize()
{

  return Gaudi::Algorithm::finalize();

}
