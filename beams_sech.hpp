//******************** FILE: BEAMS.HPP ********************
//
// Description: This header file contains the beam class definition as well as the definition of its derived class.
//
// Author: Vincent Marceau (vincent.marceau.2@ulaval.ca)
//
// Since: October 2012
// Last update: October 2012
//
//*********************************************************

#ifndef BEAMS_HPP_INCLUDED
#define BEAMS_HPP_INCLUDED

// Standard header files
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>

// GSL header files
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

// Project specific header files
#include "./constants.hpp"
#include "./general.hpp"


using namespace std;


//******************** BEAM CLASS DEFINITION ********************
// CLASS CBeam is an abstract base class. No instances of it can be created
class CBeam {
  protected:
    double E0; // Amplitude parameter [V/m]
    double lambda0; // Central wavelength [m]
    double k0; // Wavevector corresponding to lambda0 [1/m]
    double omega0; // Angular frequency corresponding to lambda0 [rad/s]
    double z0; // Confocal parameter [m]
    double w0; // Beat waist size [m]
    double phi0; // Field phase [rad]
    double T; // Pulse duration [s]
    double xi0; // Pulse duration parameter []
  public:
    double power; // Power carried by the beam [W]
    double charlength_r; // Characteristic length in the radial direction
    double charlength_z; // Characteristic length in the axial direction
    virtual vector<complex<double> > fields (double,double,double) =0; // Electromagnetic fields
    virtual void reset_phase(double) =0; // Reset the pulse phase phi0 to a given value
    virtual void printparameters (ostream&) =0; // Print beam parameters
    double errormap_fixedtime (double,double,double,double,double,double,double,double); // Compute Maxwell's equations error
                                                                                         // parameter map at fixed time (no output)
    double errormap_fixedtime (double,double,double,double,double,double,double,double,ostream&); // Compute Maxwell's equations
                                                                                                  // error parameter map at fixed
                                                                                                  // time (full output)
};


// MEMBER FUNCTION errormap_fixedtime (double,double,double,double,double,double,double,double) definition
// Compute the Maxwell's equation error parameter map at a given time (no output)
//
// Input
//   - t: time
//   - maxr: upper boundary for the r coordinate
//   - stepr: grid step for the r coordinate
//   - maxz: upper and lower boundaries for the z coordinate
//   - stepz: grid step for the z coordinate
//   - hr: step size for the computation of the r derivative
//   - hz: step size for the computation of the z derivative
//   - ht: step size for the computation of the t derivative
//   - out: valid ofstream object (output stream)
// Output
//   - epsilonMmax: maximum value of the error parameter
//
double CBeam::errormap_fixedtime (double t, double maxr, double stepr, double maxz, double stepz, double hhr, double hhz, double hht) {

  // Fields and derivatives
  double Er, Ez, Bphi;
  double dErdz, dErdt, dEzdr, dEzdt, dBphidr, dBphidz, dBphidt;

  // Make sure that ht is a representable number
  volatile double tempt = 1e-15*(t/1e-15 + hht/1e-15);
  double ht = 1e-15*(tempt/1e-15 - t/1e-15);

  // Define result container
  double Ermax = 0.0;
  double Ezmax = 0.0;
  double epsilonMmax = 0.0;
  int sizer = static_cast<int>(maxr/stepr);
  int sizez = static_cast<int>(2*maxz/stepz + 1);
  double epsilonM[sizer][sizez];

  // Iterate over map and calculate the unscaled error parameter
  for (int nr = 0; nr<sizer; nr++) {

    double r = (nr+1)*stepr; // Current r coordinate
    volatile double tempr = r + hhr; // Make sure that hr is a representable number
    double hr = tempr - r;

    for (int nz = 0; nz<sizez; nz++) {

      double z = nz*stepz - maxz; // Current z coordinate
      volatile double tempz = z + hhz; // Make sur that hz is a representable number
      double hz = tempz - z;


      // Compute the fields around the given point
      vector<complex<double> > theFields = this->fields(r,z,t);
      vector<complex<double> > theFieldsRback = this->fields(r-hr,z,t);
      vector<complex<double> > theFieldsRfront = this->fields(r+hr,z,t);
      vector<complex<double> > theFieldsZback = this->fields(r,z-hz,t);
      vector<complex<double> > theFieldsZfront = this->fields(r,z+hz,t);
      vector<complex<double> > theFieldsTback = this->fields(r,z,(t/1e-15-ht/1e-15)*1e-15);
      vector<complex<double> > theFieldsTfront = this->fields(r,z,(t/1e-15+ht/1e-15)*1e-15);

      // Compute the derivatives with three point rule
      dErdz = (real(theFieldsZfront[0]) - real(theFieldsZback[0]))/(2*hz);
      dErdt = (real(theFieldsTfront[0]) - real(theFieldsTback[0]))/(2*ht);
      dEzdr = (real(theFieldsRfront[1]) - real(theFieldsRback[1]))/(2*hr);
      dEzdt = (real(theFieldsTfront[1]) - real(theFieldsTback[1]))/(2*ht);
      dBphidr = MU0*(real(theFieldsRfront[2]) - real(theFieldsRback[2]))/(2*hr);
      dBphidz = MU0*(real(theFieldsZfront[2]) - real(theFieldsZback[2]))/(2*hz);
      dBphidt = MU0*(real(theFieldsTfront[2]) - real(theFieldsTback[2]))/(2*ht);

      // Compute field amplitudes
      Er = real(theFields[0]);
      Ez = real(theFields[1]);
      Bphi = MU0*real(theFields[2]);

      // Compute unscaled error parameter (SI units)
      epsilonM[nr][nz] = abs(dErdz-dEzdr+dBphidt) + C0*sqrt(pow(dBphidz+dErdt/pow(C0,2.0),2.0) + pow(Bphi/r+dBphidr-dEzdt/pow(C0,2.0),2.0));

      // Check for maximum field amplitude
      if (abs(Er) > Ermax)
        Ermax = Er;
      if (abs(Ez) > Ezmax)
        Ezmax = Ez;

      // Check for maximum error parameter
      if (epsilonM[nr][nz] > epsilonMmax)
        epsilonMmax = epsilonM[nr][nz];

    } // end for
  } // end for


  // Compute scaling factor
  double scale;
  if (Ermax > Ezmax)
    scale = lambda0/Ermax;
  else
    scale = lambda0/Ezmax;

  return epsilonMmax*scale;

} // end member function errormap_fixedtime (double,double,double,double,double,double,double,double) definition


// MEMBER FUNCTION errormap_fixedtime (double,double,double,double,double,double,double,double,ostream&) definition
// Compute the Maxwell's equation error parameter map at a given time (full output)
//
// Input
//   - t: time
//   - maxr: upper boundary for the r coordinate
//   - stepr: grid step for the r coordinate
//   - maxz: upper and lower boundaries for the z coordinate
//   - stepz: grid step for the z coordinate
//   - hr: step size for the computation of the r derivative
//   - hz: step size for the computation of the z derivative
//   - ht: step size for the computation of the t derivative
//   - out: valid ofstream object (output stream)
// Output
//   - epsilonMmax: maximum value of the error parameter
//
double CBeam::errormap_fixedtime (double t, double maxr, double stepr, double maxz, double stepz, double hhr, double hhz, double hht, ostream& out) {

  // Fields and derivatives
  double Er, Ez, Bphi;
  double dErdz, dErdt, dEzdr, dEzdt, dBphidr, dBphidz, dBphidt;

  // Make sure that ht is a representable number
  volatile double tempt = 1e-15*(t/1e-15 + hht/1e-15);
  double ht = 1e-15*(tempt/1e-15 - t/1e-15);


  // Define result container
  double Ermax = 0.0;
  double Ezmax = 0.0;
  double epsilonMmax = 0.0;
  int sizer = static_cast<int>(maxr/stepr);
  int sizez = static_cast<int>(2*maxz/stepz + 1);
  double epsilonM[sizer][sizez];

  // Iterate over map and calculate the unscaled error parameter
  for (int nr = 0; nr<sizer; nr++) {

    double r = (nr+1)*stepr; // Current r coordinate
    volatile double tempr = r + hhr; // Make sure that hr is a representable number
    double hr = tempr - r;

    for (int nz = 0; nz<sizez; nz++) {

      double z = nz*stepz - maxz; // Current z coordinate
      volatile double tempz = z + hhz; // Make sur that hz is a representable number
      double hz = tempz - z;

      // Compute the fields around the given point
      vector<complex<double> > theFields = this->fields(r,z,t);
      vector<complex<double> > theFieldsRback = this->fields(r-hr,z,t);
      vector<complex<double> > theFieldsRfront = this->fields(r+hr,z,t);
      vector<complex<double> > theFieldsZback = this->fields(r,z-hz,t);
      vector<complex<double> > theFieldsZfront = this->fields(r,z+hz,t);
      vector<complex<double> > theFieldsTback = this->fields(r,z,(t/1e-15-ht/1e-15)*1e-15);
      vector<complex<double> > theFieldsTfront = this->fields(r,z,(t/1e-15+ht/1e-15)*1e-15);

      // Compute the derivatives with three point rule
      dErdz = (real(theFieldsZfront[0]) - real(theFieldsZback[0]))/(2*hz);
      dErdt = (real(theFieldsTfront[0]) - real(theFieldsTback[0]))/(2*ht);
      dEzdr = (real(theFieldsRfront[1]) - real(theFieldsRback[1]))/(2*hr);
      dEzdt = (real(theFieldsTfront[1]) - real(theFieldsTback[1]))/(2*ht);
      dBphidr = MU0*(real(theFieldsRfront[2]) - real(theFieldsRback[2]))/(2*hr);
      dBphidz = MU0*(real(theFieldsZfront[2]) - real(theFieldsZback[2]))/(2*hz);
      dBphidt = MU0*(real(theFieldsTfront[2]) - real(theFieldsTback[2]))/(2*ht);

      // Compute field amplitudes
      Er = real(theFields[0]);
      Ez = real(theFields[1]);
      Bphi = MU0*real(theFields[2]);

      // Compute unscaled error parameter (SI units)
      epsilonM[nr][nz] = abs(dErdz-dEzdr+dBphidt) + C0*sqrt(pow(dBphidz+dErdt/pow(C0,2.0),2.0) + pow(Bphi/r+dBphidr-dEzdt/pow(C0,2.0),2.0));

      // Check for maximum field amplitude
      if (abs(Er) > Ermax)
        Ermax = Er;
      if (abs(Ez) > Ezmax)
        Ezmax = Ez;

      // Check for maximum error parameter
      if (epsilonM[nr][nz] > epsilonMmax)
        epsilonMmax = epsilonM[nr][nz];

    } // end for
  } // end for


  // Compute scaling factor
  double scale;
  if (Ermax > Ezmax)
    scale = lambda0/Ermax;
  else
    scale = lambda0/Ezmax;

  // Output results
  out << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out << "#" << endl;
  out << "# DESCRIPTION" << endl;
  out << "# Maxwell's equations error parameter map at t = " << t << endl;
  out << "#" << endl;
  this->printparameters(out);
  out << "#" << endl;
  out << "# FIELD PARAMETERS" << endl;
  out << "# Ermax = " << Ermax << " , Ezmax = " << Ezmax;
  out << " , scaling factor = " << scale << " , epsilonMmax = " << epsilonMmax*scale << endl;
  out << "#" << endl;
  out << "# r/lambda0  z/lambda0  epsilonM" << endl;


  for (int nr = 0; nr<sizer; nr++) {
    double r = (nr+1)*stepr; // Current r coordinate
    for (int nz = 0; nz<sizez; nz++) {
      double z = nz*stepz - maxz; // Current z coordinate
      out << r/lambda0 << "  " << z/lambda0 << "  " << epsilonM[nr][nz]*scale << endl;
    } // end for
    out << endl;
  } // end for

  return epsilonMmax*scale;

} // end member function errormap_fixedtime (double,double,double,double,double,double,double,double,ostream&) definition



//******************** EXACTTM01 CLASS DEFINITION ********************
// CLASS CExactTM01 is derived from class Cbeam
// Class whose instances represent exact TM01 pulsed beams
class CExactTM01: public CBeam {
  public:
    CExactTM01 (double,double,double,double,double,double); // Constructor
    vector<complex<double> > fields (double,double,double); // Compute the electromagnetic fields
    void reset_phase(double); // Reset the pulse phase phi0 to a given value
    void printparameters (ostream&); // Print the beam parameters to output stream
};


// CONSTRUCTOR CExactTM01 definition
// The beam's member attributes are set to the prescribed values.
//
// Input
//   - wavelength: value for the lambda0 member attribute
//   - avpower: value for the power member attribute
//   - amplitude: value of the E0 member attribute. Note: avpower and amplitude cannot be simultaneously nonzero
//   - k0z0: normalized confocal parameter
//   - duration: value for the T member attribute
//   - phase: value for the phi0 member attribute
//
CExactTM01::CExactTM01 (double wavelength, double avpower, double amplitude, double k0z0, double duration, double phase){

  // Check if avpower and amplitude aren't both nonzero
  if ((avpower != 0.0 && amplitude != 0.0) || (avpower == 0.0 && amplitude == 0.0)) {
    cout << "Error in class CExactTM01 constructor: average power and amplitude cannot be both zero or nonzero" << endl;
    exit (EXIT_FAILURE);
  } // end if

  // Set member attribute values
  lambda0 = wavelength;
  k0 = 2.0*PI/lambda0;
  omega0 = C0*k0;
  z0 = k0z0/k0;
  w0 = sqrt(2.0)*sqrt(sqrt(1.0+pow(k0z0,2.0))-1)/k0;
  phi0 = phase;
  T = duration;
  xi0 = omega0*T/1.657454454;

  if (avpower != 0.0) {
    power = avpower;
    E0 = sqrt(8.0*ETA0*pow(k0,2.0)*pow(k0*z0,3.0)*power/(PI*((k0*z0-0.5)-(k0*z0+0.5)*exp(-4.0*k0*z0)+(1.0-2.0*pow(k0*z0,2.0))*exp(-2.0*k0*z0))));
  } // end if
  else {
    E0 = amplitude;
    power = PI*pow(E0,2.0)*((k0*z0-0.5)-(k0*z0+0.5)*exp(-4.0*k0*z0)+(1.0-2.0*pow(k0*z0,2.0))*exp(-2.0*k0*z0))/(8.0*ETA0*pow(k0,2.0)*pow(k0*z0,3.0));
  } // end else

  // Set characteristic lengths
  charlength_r = w0;
  charlength_z = z0;

} // End constructor CExactTM01 definition


// MEMBER FUNCTION fields(double,double,double) definition
// The electromagnetic fields of the exact TM01 pulse are calculated according to the closed-form expressions presented
// in the PhD thesis of Alexandre April.
//
// Input
//   - r: radial coordinate
//   - z: longitudinal coordinate
//   - t: time
// Output
//   - TM01fields[0]: analytic (complex) value of the Er field component
//   - TM01fields[1]: analytic (complex) value of the Ez field component
//   - TM01fields[2]: analytic (complex) value of the Hphi field component
//
vector<complex<double> > CExactTM01::fields(double r,double z,double t)
{

  // Complex radius, complex trigonometric functions
  complex<double> Rtilde = I*z0*sqrt(1.0 - (pow(r,2.0)+pow(z,2.0)+2.0*I*z*z0)/pow(z0,2.0));
  complex<double> Rscaled = sqrt(1.0 - (pow(r,2.0)+pow(z,2.0)+2.0*I*z*z0)/pow(z0,2.0));
  complex<double> sintheta = r/(Rtilde);
  complex<double> costheta = (z+I*z0)/(Rtilde);

  // Pulse temporal dependency
  complex<double> pulse_temporal = exp(I*(omega0*t-phi0))/cosh((omega0*t-k0*z)/xi0);

  // Compute spherical Bessel functions times exp(-k0z0)
  complex<double> expSphericalBesselj0;
  complex<double> expSphericalBesselj1;
  complex<double> expSphericalBesselj2;
  if (abs(z) > 1e-14 || abs(r-z0) > 1e-12 ){
    complex<double> expSin = (exp(-k0*z0*(1.0+Rscaled)) - exp(-k0*z0*(1.0-Rscaled)))/(2.0*I);
    complex<double> expCos = (exp(-k0*z0*(1.0+Rscaled)) + exp(-k0*z0*(1.0-Rscaled)))/2.0;
    expSphericalBesselj0 = expSin/(k0*Rtilde);
    expSphericalBesselj1 = expSin/pow(k0*Rtilde,2.0) - expCos/(k0*Rtilde);
    expSphericalBesselj2 = expSin*(3.0/pow(k0*Rtilde,3.0)-1.0/(k0*Rtilde)) - 3.0*expCos/pow(k0*Rtilde,2.0);
  } // end if
  else{
    expSphericalBesselj0 = exp(-k0*z0)*complex<double>(1.0,0);
    expSphericalBesselj1 = complex<double>(0,0);
    expSphericalBesselj2 = complex<double>(0,0);
  } // end else

  // Radial electric field component Er
  complex<double> Er = -I*E0*expSphericalBesselj2*sintheta*costheta*pulse_temporal;

  // Longitudinal electric field component Ez
  complex<double> Ez = -(2.0/3.0)*I*E0*(expSphericalBesselj0 - expSphericalBesselj2*(1.0-3.0*pow(costheta,2.0))/2.0)*pulse_temporal;

  // Azimutal magnetic field component Hphi
  complex<double> Hphi = (E0/ETA0)*expSphericalBesselj1*sintheta*pulse_temporal;


  // Return the components of the electromagnetic field
  vector<complex<double> > TM01fields (3,complex<double>(0,0));
  TM01fields[0] = Er;
  TM01fields[1] = Ez;
  TM01fields[2] = Hphi;

  return TM01fields;

} // End member function fields(double,double,double) definition


// MEMBER FUNCTION reset_phase(double) definition
// Reset the pulse phase to a given value
//
// Input
//   - newphi0: new value of the pulse phase
// Output
//     none
//
void CExactTM01::reset_phase(double newphi0) {

  phi0 = newphi0;

} // end member function reset_phase(double) definition


// MEMBER FUNCTION printparameters(ofstream) definition
// Print the pulsed beam parameters to a given output stream
//
// Input
//   - out: valid ofstream object (output stream)
// Output
//     none
//
void CExactTM01::printparameters(ostream& out) {

  out << "# BEAM PARAMETERS" << endl;
  out << "# Type: Exact TM01 pulsed beam" << endl;
  out << "# Central wavelength: lambda0 = " << lambda0 << endl;
  out << "# Average power: P = " << power << endl;
  out << "# Amplitude parameter: E0 = " << E0 << endl;
  out << "# Confocal parameter: z0 = " << z0 << endl;
  out << "# Scaled confocal parameter: k0z0 = " << k0*z0 << endl;
  out << "# Pulse duration: T = " << T << endl;
  out << "# Pulse duration parameter: xi0 = " << xi0 << endl;
  out << "# Pulse phase: phi0 = " << phi0 << endl;
  out << "# Beam waist size at wavelength lambda0: w0 = " << w0 << endl;

} // end member function printparameters(ofstream) definition



//******************** PARAXIALTM01 CLASS DEFINITION ********************
// CLASS CParaxialTM01 is derived from class Cbeam
// Class whose instances represent paraxial TM01 pulsed beams
class CParaxialTM01: public CBeam {
  public:
    CParaxialTM01 (double,double,double,double,double,double); // Constructor
    vector<complex<double> > fields (double,double,double); // Compute the electromagnetic fields
    void reset_phase(double); // Reset the pulse phase phi0 to a given value
    void printparameters (ostream&); // Print the beam parameters to output stream
};


// CONSTRUCTOR CParaxialTM01 definition
// The beam's member attributes are set to the prescribed values.
//
// Input
//   - wavelength: value for the lambda0 member attribute
//   - avpower: value for the power member attribute
//   - amplitude: value of the E0 member attribute. Note: avpower and amplitude cannot be simultaneously nonzero
//   - k0z0: normalized confocal parameter
//   - duration: value for the T member attribute
//   - phase: value for the phi0 member attribute
//
CParaxialTM01::CParaxialTM01 (double wavelength, double avpower, double amplitude, double k0z0, double duration, double phase){

  // Check if avpower and amplitude aren't both nonzero
  if ((avpower != 0.0 && amplitude != 0.0) || (avpower == 0.0 && amplitude == 0.0)) {
    cout << "Error in class CParaxialTM01 constructor: average power and amplitude cannot be both zero or nonzero" << endl;
    exit (EXIT_FAILURE);
  } // end if

  // Set member attribute values
  lambda0 = wavelength;
  k0 = 2.0*PI/lambda0;
  omega0 = C0*k0;
  z0 = k0z0/k0;
  w0 = sqrt(2.0*z0/k0);
  phi0 = phase;
  T = duration;
  xi0 = omega0*T/1.657454454;

  if (avpower != 0.0) {
    power = avpower;
    E0 = sqrt(8.0*ETA0*pow(k0,2.0)*pow(k0*z0,3.0)*power/(PI*((k0*z0-0.5)-(k0*z0+0.5)*exp(-4.0*k0*z0)+(1.0-2.0*pow(k0*z0,2.0))*exp(-2.0*k0*z0))));
  } // end if
  else {
    E0 = amplitude;
    power = PI*pow(E0,2.0)*((k0*z0-0.5)-(k0*z0+0.5)*exp(-4.0*k0*z0)+(1.0-2.0*pow(k0*z0,2.0))*exp(-2.0*k0*z0))/(8.0*ETA0*pow(k0,2.0)*pow(k0*z0,3.0));
  } // end else

  // Set characteristic lengths
  charlength_r = w0;
  charlength_z = z0;

} // End constructor CParaxialTM01 definition


// MEMBER FUNCTION fields(double,double,double) definition
// The electromagnetic fields of the paraxial TM01 pulsed beam correspond to the first term of perturbative expansion of the exact fields
//
// Input
//   - r: radial coordinate
//   - z: longitudinal coordinate
//   - t: time
// Output
//   - TM01fields[0]: analytic (complex) value of the Er field component
//   - TM01fields[1]: analytic (complex) value of the Ez field component
//   - TM01fields[2]: analytic (complex) value of the Hphi field component
//
vector<complex<double> > CParaxialTM01::fields(double r,double z,double t)
{

  // Compute common quantities
  complex<double> qtilde = z+I*z0;
  complex<double> pulse_spatial = exp(-I*k0*pow(r,2.0)/(2.0*qtilde));
  complex<double> pulse_temporal = 1.0/cosh((omega0*t-k0*z)/xi0);
  complex<double> propagator = exp(I*(omega0*t-k0*z-phi0));
  complex<double> common_factor = pulse_spatial*pulse_temporal*propagator;

  // Radial electric field component Er
  complex<double> Er = -E0*r/(2.0*k0*pow(qtilde,2.0))*common_factor;

  // Longitudinal electric field component Ez
  complex<double> Ez = (I*E0/pow(k0*qtilde,2.0))*(1.0 - I*k0*pow(r,2.0)/(2.0*qtilde))*common_factor;

  // Azimutal magnetic field component Hphi
  complex<double> Hphi = Er/ETA0;

  // Return the components of the electromagnetic field
  vector<complex<double> > TM01fields (3,complex<double>(0,0));
  TM01fields[0] = Er;
  TM01fields[1] = Ez;
  TM01fields[2] = Hphi;

  return TM01fields;

} // End member function fields(double,double,double) definition


// MEMBER FUNCTION reset_phase(double) definition
// Reset the pulse phase to a given value
//
// Input
//   - newphi0: new value of the pulse phase
// Output
//     none
//
void CParaxialTM01::reset_phase(double newphi0) {

  phi0 = newphi0;

} // end member function reset_phase(double) definition


// MEMBER FUNCTION printparameters(ofstream) definition
// Print the pulsed beam parameters to a given output stream
//
// Input
//   - out: valid ofstream object (output stream)
// Output
//     none
//
void CParaxialTM01::printparameters(ostream& out) {

  out << "# BEAM PARAMETERS" << endl;
  out << "# Type: Paraxial TM01 pulsed beam" << endl;
  out << "# Central wavelength: lambda0 = " << lambda0 << endl;
  out << "# Average power: P = " << power << endl;
  out << "# Amplitude parameter: E0 = " << E0 << endl;
  out << "# Confocal parameter: z0 = " << z0 << endl;
  out << "# Scaled confocal parameter: k0z0 = " << k0*z0 << endl;
  out << "# Pulse duration: T = " << T << endl;
  out << "# Pulse duration parameter: xi0 = " << xi0 << endl;
  out << "# Pulse phase: phi0 = " << phi0 << endl;
  out << "# Beam waist size at wavelength lambda0: w0 = " << w0 << endl;

} // end member function printparameters(ofstream) definition



//******************** PARAXIAL1STCORRTM01 CLASS DEFINITION ********************
// CLASS CParaxial1stCorrTM01 is derived from class Cbeam
// Class whose instances represent paraxial TM01 pulsed beams with first non paraxial corrections to Er and Hphi
class CParaxial1stCorrTM01: public CBeam {
  public:
    CParaxial1stCorrTM01 (double,double,double,double,double,double); // Constructor
    vector<complex<double> > fields (double,double,double); // Compute the electromagnetic fields
    void reset_phase(double); // Reset the pulse phase phi0 to a given value
    void printparameters (ostream&); // Print the beam parameters to output stream
};


// CONSTRUCTOR CParaxial1stCorrTM01 definition
// The beam's member attributes are set to the prescribed values.
//
// Input
//   - wavelength: value for the lambda0 member attribute
//   - avpower: value for the power member attribute
//   - amplitude: value of the E0 member attribute. Note: avpower and amplitude cannot be simultaneously nonzero
//   - k0z0: normalized confocal parameter
//   - duration: value for the T member attribute
//   - phase: value for the phi0 member attribute
//
CParaxial1stCorrTM01::CParaxial1stCorrTM01 (double wavelength, double avpower, double amplitude, double k0z0, double duration, double phase){

  // Check if avpower and amplitude aren't both nonzero
  if ((avpower != 0.0 && amplitude != 0.0) || (avpower == 0.0 && amplitude == 0.0)) {
    cout << "Error in class CParaxial1stCorrTM01 constructor: average power and amplitude cannot be both zero or nonzero" << endl;
    exit (EXIT_FAILURE);
  } // end if

  // Set member attribute values
  lambda0 = wavelength;
  k0 = 2.0*PI/lambda0;
  omega0 = C0*k0;
  z0 = k0z0/k0;
  w0 = sqrt(2.0*z0/k0);
  phi0 = phase;
  T = duration;
  xi0 = omega0*T/1.657454454;

  if (avpower != 0.0) {
    power = avpower;
    E0 = sqrt(8.0*ETA0*pow(k0,2.0)*pow(k0*z0,3.0)*power/(PI*((k0*z0-0.5)-(k0*z0+0.5)*exp(-4.0*k0*z0)+(1.0-2.0*pow(k0*z0,2.0))*exp(-2.0*k0*z0))));
  } // end if
  else {
    E0 = amplitude;
    power = PI*pow(E0,2.0)*((k0*z0-0.5)-(k0*z0+0.5)*exp(-4.0*k0*z0)+(1.0-2.0*pow(k0*z0,2.0))*exp(-2.0*k0*z0))/(8.0*ETA0*pow(k0,2.0)*pow(k0*z0,3.0));
  } // end else

  // Set characteristic lengths
  charlength_r = w0;
  charlength_z = z0;

} // End constructor CExactTM01 definition


// MEMBER FUNCTION fields(double,double,double) definition
// The electromagnetic fields of the corrected paraxial TM01 pulse correspond to the first terms in the perturbative expansion of
// the fields plus the second terms for Er and Hphi components.
//
// Input
//   - r: radial coordinate
//   - z: longitudinal coordinate
//   - t: time
// Output
//   - TM01fields[0]: analytic (complex) value of the Er field component
//   - TM01fields[1]: analytic (complex) value of the Ez field component
//   - TM01fields[2]: analytic (complex) value of the Hphi field component
//
vector<complex<double> > CParaxial1stCorrTM01::fields(double r,double z,double t)
{

  // Compute common quantities
  complex<double> qtilde = z+I*z0;
  complex<double> pulse_spatial = exp(-I*k0*pow(r,2.0)/(2.0*qtilde));
  complex<double> pulse_temporal = 1.0/cosh((omega0*t-k0*z)/xi0);
  complex<double> propagator = exp(I*(omega0*t-k0*z-phi0));
  complex<double> common_factor = pulse_spatial*pulse_temporal*propagator;

  // Radial electric field component Er
  complex<double> Er = -(E0/(2.0*k0))*(r/pow(qtilde,2.0) - 3.0*pow(r,3.0)/(2.0*pow(qtilde,4.0)) + I*k0*pow(r,5.0)/(8.0*pow(qtilde,5.0)) - 3.0*I*r/(k0*pow(qtilde,3.0)))*common_factor;

  // Longitudinal electric field component Ez
  complex<double> Ez = (I*E0/pow(k0*qtilde,2.0))*(1.0 - I*k0*pow(r,2.0)/(2.0*qtilde))*common_factor;

  // Azimutal magnetic field component Hphi
  complex<double> Hphi = -(E0/(2.0*ETA0*k0))*(r/pow(qtilde,2.0) - pow(r,3.0)/pow(qtilde,4.0) + I*k0*pow(r,5.0)/(8.0*pow(qtilde,5.0)) - I*r/(k0*pow(qtilde,3.0)))*common_factor;


  // Return the components of the electromagnetic field
  vector<complex<double> > TM01fields (3,complex<double>(0,0));
  TM01fields[0] = Er;
  TM01fields[1] = Ez;
  TM01fields[2] = Hphi;

  return TM01fields;

} // End member function fields(double,double,double) definition


// MEMBER FUNCTION reset_phase(double) definition
// Reset the pulse phase to a given value
//
// Input
//   - newphi0: new value of the pulse phase
// Output
//     none
//
void CParaxial1stCorrTM01::reset_phase(double newphi0) {

  phi0 = newphi0;

} // end member function reset_phase(double) definition


// MEMBER FUNCTION printparameters(ofstream) definition
// Print the pulsed beam parameters to a given output stream
//
// Input
//   - out: valid ofstream object (output stream)
// Output
//     none
//
void CParaxial1stCorrTM01::printparameters(ostream& out) {

  out << "# BEAM PARAMETERS" << endl;
  out << "# Type: Corrected paraxial TM01 pulsed beam" << endl;
  out << "# Central wavelength: lambda0 = " << lambda0 << endl;
  out << "# Average power: P = " << power << endl;
  out << "# Amplitude parameter: E0 = " << E0 << endl;
  out << "# Confocal parameter: z0 = " << z0 << endl;
  out << "# Scaled confocal parameter: k0z0 = " << k0*z0 << endl;
  out << "# Pulse duration: T = " << T << endl;
  out << "# Pulse duration parameter: xi0 = " << xi0 << endl;
  out << "# Pulse phase: phi0 = " << phi0 << endl;
  out << "# Beam waist size at wavelength lambda0: w0 = " << w0 << endl;

} // end member function printparameters(ofstream) definition


#endif // BEAMS_HPP_INCLUDED
