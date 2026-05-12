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


  double DESY_SPR(double now) {
    return 0.0;
  }
  

  double FNAL2023_SPR(double now) 
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



  PulseSpline::PulseSpline(PulseType ptype)
  {
    // start with default vaules
    par[0]=0.0;  // t0
    par[1]=1.0;  // amplituide
    par[2]=0.0;  // baseline

    switch(ptype) {
    case PulseType::DESY24_PbF2: {
      // DESY 2024, Typical PbF2 pulse shape
      nk_ = 11;
      double p[11][6] = {
	{     0 ,            0 ,            0 ,   0.00380779 ,   -0.0045358 ,   0.00194151 }  ,
	{   2.5 ,    0.0287671 ,    0.0553372 ,    0.0773282 ,    -0.015079 ,   0.00085068 }  ,
	{   7.5 ,     0.885459 ,     0.123035 ,   -0.0473158 ,   0.00754599 , -0.000414265 }  ,
	{  11.5 ,     0.997435 ,  0.000663746 ,  -0.00191597 , -7.03457e-05 ,  3.14136e-05 }  ,
	{  15.5 ,     0.972975 ,   -0.0099987 ,   0.00107499 , -8.06713e-05 ,  2.19843e-06 }  ,
	{    25 ,     0.923746 ,  -0.00387615 ,  8.77992e-05 , -1.04601e-05 ,  1.65228e-07 }  ,
	{    40 ,      0.85842 ,  -0.00607216 ,   2.9781e-06 ,  1.54456e-07 , -5.51846e-10 }  ,
	{   200 ,      0.23411 ,  -0.00229836 ,  1.94217e-05 , -1.36284e-07 ,  4.09358e-10 }  ,
	{   300 ,     0.103143 , -0.000865112 ,  4.47868e-06 , -1.43529e-08 ,  1.92479e-11 }  ,
	{   500 ,    0.0252406 , -0.000180058 ,  6.70644e-07 , -4.97031e-10 , -1.75271e-12 }  ,
	{   750 ,   0.00752879 ,  -0.00630559 ,            0 ,            0 ,            0 }
      };
      for(int i=0; i<nk_; i++){
	for(int j=0; j<np_; j++){
	  p_[i][j] = p[i][j];
	}
      }
      break; 
    }
    case PulseType::JLAB25_PbF2LG: {
      // JLAB 2025, Typical PbF2 pulse shape in Low Gain
      nk_ = 11;
      double p[11][6] = {
	{     0 ,            0 ,            0 ,    0.0198848 ,    0.0269754 ,    0.0347858 }  ,
	{   0.7 ,    0.0273482 ,     0.115219 ,      0.25472 ,   -0.0444562 ,  -0.00724739 }  ,
	{     3 ,     0.896112 ,     0.228697 ,    -0.195429 ,     0.130166 ,   -0.0567278 }  ,
	{   4.1 ,      1.00141 ,   -0.0307628 ,  -0.00211941 ,   0.00211809 , -0.000757375 }  ,
	{     7 ,     0.892459 ,   -0.0635024 ,  -0.00345521 ,  0.000854172 , -3.59311e-05 }  ,
	{    15 ,     0.453469 ,   -0.0283716 ,    0.0011344 , -2.81671e-05 ,  2.74625e-07 }  ,
	{    50 ,    0.0545451 ,  -0.00537969 ,  0.000177559 , -3.56658e-06 ,  2.80642e-08 }  ,
	{    80 ,   -0.0206085 ,  -0.00132501 ,  2.82877e-05 ,  -3.1703e-07 ,   1.6496e-09 }  ,
	{   130 ,   -0.0454583 , -4.91609e-05 ,  4.91927e-06 ,  -4.9982e-08 ,  1.82029e-10 }  ,
	{   200 ,   -0.0375684 ,  0.000154545 ,  2.48105e-07 , -2.74739e-09 ,   4.2938e-12 }  ,
	{   400 ,   -0.0118442 ,  -0.00519259 ,            0 ,            0 ,            0 }
      };
      for(int i=0; i<nk_; i++){
	for(int j=0; j<np_; j++){
	  p_[i][j] = p[i][j];
	}
      }
      break;
    }
    case PulseType::JLAB25_PbF2HG: {
      // JLAB 2025, Typical PbF2 pulse shape in High Gain
      nk_ = 8;
      double p[8][6] = {
	{     0 ,            0 ,            0 ,   0.00902277 ,    0.0124065 ,    0.0161493 }  ,
	{   0.7 ,     0.012554 ,    0.0530262 ,     0.221329 ,   -0.0191583 ,   -0.0105919 }  ,
	{     3 ,     0.775842 ,     0.251612 ,   -0.0936897 ,   0.00549861 ,   0.00244342 }  ,
	{   5.1 ,     0.989498 ,    0.0213759 ,  -0.00801679 , -0.000515797 ,  0.000153308 }  ,
	{   9.2 ,     0.950148 ,   -0.0281088 ,  0.000296773 , -1.55071e-06 ,  3.71353e-09 }  ,
	{    60 ,     0.109521 ,  -0.00801487 ,  0.000121226 ,  -9.4282e-07 ,  3.22049e-09 }  ,
	{   150 ,    -0.105909 ,  0.000286165 ,  1.09346e-06 , -5.53985e-09 ,  6.14741e-12 }  ,
	{   400 ,   -0.0285731 ,  -0.00624314 ,            0 ,            0 ,            0 }
      };
      for(int i=0; i<nk_; i++){
	for(int j=0; j<np_; j++){
	  p_[i][j] = p[i][j];
	}
      }
      break;
    }
    case PulseType::SLJan26_SPR: {
      // Broadcom + S.Los amp, Jan 2026. 1-pe pulse shape
      nk_ = 10;
      double p[10][6] = {
	{     0 ,            0 ,            0 ,       0.1854 ,   -0.0783658 ,     0.116088 }  ,
	{     1 ,     0.223123 ,     0.600056 ,     0.564185 ,    -0.568587 ,     0.102022 }  ,
	{     2 ,     0.920799 ,     0.430754 ,   -0.0756142 ,     -2.38348 ,      2.38137 }  ,
	{   2.5 ,     0.968173 ,    -0.241781 ,   -0.0140137 ,    0.0176933 ,  -0.00234401 }  ,
	{     5 ,     0.461029 ,    -0.126602 ,    0.0151215 , -0.000867766 ,  1.84769e-05 }  ,
	{    15 ,    0.0241648 ,   -0.0105938 ,   0.00125789 , -5.58816e-05 ,  8.16153e-07 }  ,
	{    40 ,  -0.00884069 ,  -0.00146784 ,   6.8088e-05 , -1.73965e-06 ,  2.62483e-08 }  ,
	{    60 ,   -0.0206798 ,  8.04551e-06 , -9.15471e-07 ,  5.06326e-08 , -4.17816e-10 }  ,
	{   130 ,   -0.0172672 ,  5.09341e-05 ,  9.46525e-07 ,  1.56585e-09 , -8.63305e-11 }  ,
	{   200 ,   -0.0105996 ,  -0.00830411 ,            0 ,            0 ,            0 }
      };
      for(int i=0; i<nk_; i++){
	for(int j=0; j<np_; j++){
	  p_[i][j] = p[i][j];
	}
      }
      break;
    }
    case PulseType::DESY24_SPR: {
      // DESY 2024 TB, 1-pe pulse shape
      nk_ = 7;
      double p[7][6] = {
	{     0 ,            0 ,            0 ,    0.0150194 ,  -0.00083076 , -4.80737e-06 }  ,
	{   5.5 ,     0.311721 ,     0.086623 ,    0.0138609 ,  -0.00266192 ,  0.000108897 }  ,
	{  15.5 ,     0.991088 ,  0.000852331 , -0.000359749 ,   1.2507e-05 , -1.60501e-07 }  ,
	{  50.5 ,     0.875611 ,  -0.00589292 ,   6.8182e-05 , -1.66198e-06 ,  1.43062e-08 }  ,
	{ 100.5 ,     0.633086 ,  -0.00438649 ,  6.59134e-06 ,  9.43806e-08 , -5.27449e-10 }  ,
	{ 200.5 ,     0.301986 ,   -0.0023466 ,  1.09794e-05 , -2.84432e-08 ,  2.71366e-11 }  ,
	{   400 ,     0.087965 ,  -0.00568492 ,            0 ,            0 ,            0 }
      };
      for(int i=0; i<nk_; i++){
	for(int j=0; j<np_; j++){
	  p_[i][j] = p[i][j];
	}
      }
      break;
    }
    }; 
    
  }

  PulseSpline::~PulseSpline() {}

  double PulseSpline::Evaluate(double *x, double *par)
  {
    double t = x[0] - par[0];
    double f = 0;
    if(t < p_[0][0]){
      f = 0;
    }else if( t >= p_[nk_-1][0] ){
      f = p_[nk_-1][1] * exp( p_[nk_-1][2] * (t - p_[nk_-1][0]));
    }else{
      for(int i=0; i<nk_-1; i++){
	if(t >= p_[i][0] && t < p_[i+1][0]){
	  for(int j=0; j<np_-1; j++){
	    f += p_[i][j+1] * pow(t - p_[i][0], j);
	  }
	}
      }
    }
    return f * par[1] + par[2];
  }

  double PulseSpline::Eval(double x)
  {
    return Evaluate(&x, par);
  }

  void PulseSpline::SetParameters(double *newpars)
  {
    for (int i=0; i<3; ++i) par[i]=newpars[i];
  }

  
}
