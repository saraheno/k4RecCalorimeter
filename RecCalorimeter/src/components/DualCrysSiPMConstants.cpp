#include "DualCrysSiPMConstants.h"
#include "math.h" 


/* Find the nearest wavelength from our fixed set of wavelengths
 This works by rounding the wavelength off to the nearest integer
 then finding the lowest value greater than given rounded wavelength
 and then checking the nearest neighbors for the closest value and returning
 the pair of nearest wavelength and PDE for a given SiPM type
*/

const std::pair<double,double> findnearest(Filter_Type ftype, double wavelength)  {
  if (wavelength < 0.) 
    return std::make_pair(0.0,0.0); 

  auto lb = u330_filter_wavelengths.lower_bound(round(wavelength)); 
  if (ftype == Filter_Type::O58) {
    lb = o58_filter_wavelengths.lower_bound(round(wavelength));
  }


  
  int start = *lb;
  int delta = abs(wavelength - *lb);
  lb--;
  int lowerdelta = abs(wavelength - *lb);
  if (lowerdelta < delta) {
    start = *lb;
    delta = lowerdelta;
  }
  lb++;
  lb++;
  int upperdelta = abs(wavelength - *lb);
  if (upperdelta < delta) {
    start = *lb;
  }

  if (ftype == Filter_Type::U330)
    return std::make_pair(start, u330_filterMap.at(start));
  else if (ftype == Filter_Type::O58)
    return std::make_pair(start, o58_filterMap.at(start));
  else
    return std::make_pair(0.0,0.0); 


}

const std::pair<double,double> findnearest(SiPM_Type stype, double wavelength)  {
  //info() << "FindNearest Wavelength:" << wavelength << endmsg;

  if (wavelength < 0.) 
    return std::make_pair(0.0,0.0); 

  auto lb = UV_Wavelengths.lower_bound(round(wavelength)); 
  if (stype == SiPM_Type::RGB) {
    lb = RGB_Wavelengths.lower_bound(round(wavelength));
  }


  
  int start = *lb;
  int delta = abs(wavelength - *lb);
  lb--;
  int lowerdelta = abs(wavelength - *lb);
  if (lowerdelta < delta) {
    start = *lb;
    delta = lowerdelta;
  }
  lb++;
  lb++;
  int upperdelta = abs(wavelength - *lb);
  if (upperdelta < delta) {
    start = *lb;
  }

  if (stype == SiPM_Type::RGB)
    return RGB_Map.at(start);
  else if (stype == SiPM_Type::UV)
    return UV_Map.at(start); 
  else
    return std::make_pair(0.0,0.0); 

}

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
