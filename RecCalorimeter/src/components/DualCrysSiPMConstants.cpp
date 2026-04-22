#include "DualCrysSiPMConstants.h"
#include "math.h" 
#include <iostream> 
#include <ostream>
namespace calvision { 
  /* Find the nearest wavelength from our fixed set of wavelengths
     This works by rounding the wavelength off to the nearest integer
     then finding the lowest value greater than given rounded wavelength
     and then checking the nearest neighbors for the closest value and returning
     the pair of nearest wavelength and PDE for a given SiPM type
  */


  double DESY_SPR(double now) 
  {

    //  double tMin_  = 0.0;
    //  double tMax_  = 1000.0;

    double tRise       = 0.853;
    double tDecay      = 6.538;
    double tUnderShoot = 101.7;
    //  double norm        = 0.111051;


    double a = 1./ tRise;
    double b = 1./ tDecay;
    double A = -a * b / (a - b);
    double B = -A;
    double result = A * exp(-a*now) + B * exp(-b*now);
    
    double g = 1./ tUnderShoot;
    double Atmp = -A * g / ( a - g);
    double Btmp = -B * g / ( b - g);
    double G = - Atmp - Btmp ;
    A = Atmp;
    B = Btmp;
    result -= A * exp(-a*now) + B * exp(-b*now) + G * exp(-g*now);
    
    return result;

  }

  bool operator<(const key &lhs, const key &rhs) {
    int lvalue = (lhs.ix<<3|lhs.iy<<10|lhs.layer<<20); 
    int rvalue = (rhs.ix<<3|rhs.iy<<10|rhs.layer<<20); 
    if (lvalue < rvalue)
      return true;
    else
      return false;

  }

  unpackedcellID unpackCellID(uint64_t cellID) {

    unpackedcellID cid;

    int slice_id = ((0x7<<17&cellID)>>17);
    int layer_id = ((0x7<<20&cellID)>>20);
    int ix = (0x7f<<3&cellID)>>3;
    int iy = (0x7f<<10&cellID)>>10;
    if (ix > 63)
      ix -= 64;
    if (iy > 63)
      iy -= 64;


    cid.system = 0x7&cellID; 
    cid.ix = ix;
    cid.iy = iy;
    cid.slice = slice_id;
    cid.layer = layer_id;
    cid.wc = (0x7<<23&cellID)>>23;
    cid.wc2 = (0x7<<26&cellID)>>26;
    cid.wc3 = (0x7<<29&cellID)>>29;
    return cid; 

  }


  thread_local ROOT::Math::Interpolator u330_filter;
  thread_local ROOT::Math::Interpolator o58_filter;
  thread_local ROOT::Math::Interpolator broadcom_2x1_sipm_filter;
  thread_local ROOT::Math::Interpolator rgb_sipm_filter;
  thread_local ROOT::Math::Interpolator uv_sipm_filter;

  thread_local bool filterInit; 
  std::mutex guard;   
  
  bool init_filters() 
  {
    //std::lock_guard<std::mutex> lg(guard);
    if (calvision::filterInit)
      return true;
    else { 
      std::cout << "Init U330 Filter" << std::endl; 
      calvision::u330_filter.SetData(u330_wavelengths, u330_fltreff);
      std::cout << "Init O58 Filter" << std::endl; 
      calvision::o58_filter.SetData(o58_wavelengths, o58_fltreff);
      std::cout << "Init Broadcom SiPM Filter" << std::endl;
      calvision::broadcom_2x1_sipm_filter.SetData(Broadcom_2x1_Wvl, Broadcom_2x1_Eff);
      std::cout << "Init RGB SiPM Filter" << std::endl; 
      calvision::rgb_sipm_filter.SetData(RGB_Wvl, RGB_Eff); 
      std::cout << "Init UV SiPM Filter" << std::endl; 
      calvision::uv_sipm_filter.SetData(UV_Wvl, UV_Eff);
      calvision::filterInit = true;
    }
    return true; 
  }

  
}
