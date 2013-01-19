//******************** FILE: 2DACC_BEAMCHARACTERIZATION.CPP ********************
//
// Description: Source file that contains algorithms to characterize the beam objects
//
// Author: Vincent Marceau (vincent.marceau.2@ulaval.ca)
//
// Since: October 2012
// Last update: October 2012
//
//******************************************************************************

// Standard header files
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>
#include <ctime>

// Project specific header files
#include "./beams.hpp"
#include "./constants.hpp"
#include "./general.hpp"


using namespace std;


// Main function
int main() {

  /*// FEATURE TO BE CHARACTERIZED
  // Density of electric energy versus radial coordinate r and time t at beam waist

  // Definition of the pulsed beams to be used
  double lambda0 = 800e-9;
  double power = 1.0e14;
  double duration = 10e-15;
  double phase = 0.0;
  // Exact beams
  CExactTM01 myExactBeam1(lambda0,power,1.0,duration,phase);
  CExactTM01 myExactBeam5(lambda0,power,5.0,duration,phase);
  CExactTM01 myExactBeam10(lambda0,power,10.0,duration,phase);
  CExactTM01 myExactBeam20(lambda0,power,20.0,duration,phase);
  CExactTM01 myExactBeam50(lambda0,power,50.0,duration,phase);
  CExactTM01 myExactBeam200(lambda0,power,200.0,duration,phase);
  // Paraxial beams
  CParaxialTM01 myParaxialBeam1(lambda0,power,1.0,duration,phase);
  CParaxialTM01 myParaxialBeam5(lambda0,power,5.0,duration,phase);
  CParaxialTM01 myParaxialBeam10(lambda0,power,10.0,duration,phase);
  CParaxialTM01 myParaxialBeam20(lambda0,power,20.0,duration,phase);
  CParaxialTM01 myParaxialBeam50(lambda0,power,50.0,duration,phase);
  CParaxialTM01 myParaxialBeam200(lambda0,power,200.0,duration,phase);
  // Corrected Paraxial beams
  CParaxial1stCorrTM01 myParaxialCorrBeam1(lambda0,power,1.0,duration,phase);
  CParaxial1stCorrTM01 myParaxialCorrBeam5(lambda0,power,5.0,duration,phase);
  CParaxial1stCorrTM01 myParaxialCorrBeam10(lambda0,power,10.0,duration,phase);
  CParaxial1stCorrTM01 myParaxialCorrBeam20(lambda0,power,20.0,duration,phase);
  CParaxial1stCorrTM01 myParaxialCorrBeam50(lambda0,power,50.0,duration,phase);
  CParaxial1stCorrTM01 myParaxialCorrBeam200(lambda0,power,200.0,duration,phase);

  // Calculation for EXACT BEAMS
  // Open output file
  ofstream out1;
  out1.open("./dat/test_beamcharacterization/ExactTM01_we_vs_r.dat", ios::out);
  out1 << "# 2DACCELERATION_SVEA PROJECT OUTPUT FILE" << endl;
  out1 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out1 << "#" << endl;
  out1 << "# DESCRIPTION" << endl;
  out1 << "# Electric field energy density at beam waist versus radial coordinate" << endl;
  out1 << "#" << endl;
  myExactBeam1.printparameters(out1);
  out1 << "#" << endl;
  out1 << "# r/lambda0  w_e(1)  w_e(5)  w_e(10)  w_e(20)  w_e(50)  w_e(200)" << endl;
  for (double r=0.0*lambda0; r<=6.0*lambda0; r+=lambda0/100.0) {
    vector<complex<double> > fields1 = myExactBeam1.fields(r,0.0,0.0);
    vector<complex<double> > fields5 = myExactBeam5.fields(r,0.0,0.0);
    vector<complex<double> > fields10 = myExactBeam10.fields(r,0.0,0.0);
    vector<complex<double> > fields20 = myExactBeam20.fields(r,0.0,0.0);
    vector<complex<double> > fields50 = myExactBeam50.fields(r,0.0,0.0);
    vector<complex<double> > fields200 = myExactBeam200.fields(r,0.0,0.0);
    double w_e1 = EPS0*(pow(abs(fields1[0]),2.0) + pow(abs(fields1[1]),2.0))/2.0;
    double w_e5 = EPS0*(pow(abs(fields5[0]),2.0) + pow(abs(fields5[1]),2.0))/2.0;
    double w_e10 = EPS0*(pow(abs(fields10[0]),2.0) + pow(abs(fields10[1]),2.0))/2.0;
    double w_e20 = EPS0*(pow(abs(fields20[0]),2.0) + pow(abs(fields20[1]),2.0))/2.0;
    double w_e50 = EPS0*(pow(abs(fields50[0]),2.0) + pow(abs(fields50[1]),2.0))/2.0;
    double w_e200 = EPS0*(pow(abs(fields200[0]),2.0) + pow(abs(fields200[1]),2.0))/2.0;
    out1 << r/lambda0 << "  " << w_e1 << "  " << w_e5 << "  " << w_e10;
    out1 << "  " << w_e20 << "  " << w_e50 << "  " << w_e200 << endl;
  } // end for
  // Close output file
  out1 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out1.close();

  // Calculation for PARAXIAL BEAMS
  // Open output file
  ofstream out2;
  out2.open("./dat/test_beamcharacterization/ParaxialTM01_we_vs_r.dat", ios::out);
  out2 << "# 2DACCELERATION_SVEA PROJECT OUTPUT FILE" << endl;
  out2 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out2 << "#" << endl;
  out2 << "# DESCRIPTION" << endl;
  out2 << "# Electric field energy density at beam waist versus radial coordinate" << endl;
  out2 << "#" << endl;
  myParaxialBeam1.printparameters(out2);
  out2 << "#" << endl;
  out2 << "# r/lambda0  w_e(1)  w_e(5)  w_e(10)  w_e(20)  w_e(50)  w_e(200)" << endl;
  for (double r=0.0*lambda0; r<=6.0*lambda0; r+=lambda0/100.0) {
    vector<complex<double> > fields1 = myParaxialBeam1.fields(r,0.0,0.0);
    vector<complex<double> > fields5 = myParaxialBeam5.fields(r,0.0,0.0);
    vector<complex<double> > fields10 = myParaxialBeam10.fields(r,0.0,0.0);
    vector<complex<double> > fields20 = myParaxialBeam20.fields(r,0.0,0.0);
    vector<complex<double> > fields50 = myParaxialBeam50.fields(r,0.0,0.0);
    vector<complex<double> > fields200 = myParaxialBeam200.fields(r,0.0,0.0);
    double w_e1 = EPS0*(pow(abs(fields1[0]),2.0) + pow(abs(fields1[1]),2.0))/2.0;
    double w_e5 = EPS0*(pow(abs(fields5[0]),2.0) + pow(abs(fields5[1]),2.0))/2.0;
    double w_e10 = EPS0*(pow(abs(fields10[0]),2.0) + pow(abs(fields10[1]),2.0))/2.0;
    double w_e20 = EPS0*(pow(abs(fields20[0]),2.0) + pow(abs(fields20[1]),2.0))/2.0;
    double w_e50 = EPS0*(pow(abs(fields50[0]),2.0) + pow(abs(fields50[1]),2.0))/2.0;
    double w_e200 = EPS0*(pow(abs(fields200[0]),2.0) + pow(abs(fields200[1]),2.0))/2.0;
    out2 << r/lambda0 << "  " << w_e1 << "  " << w_e5 << "  " << w_e10;
    out2 << "  " << w_e20 << "  " << w_e50 << "  " << w_e200 << endl;
  } // end for
  // Close output file
  out2 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out2.close();

  // Calculation for CORRECTED PARAXIAL BEAMS
  // Open output file
  ofstream out3;
  out3.open("./dat/test_beamcharacterization/ParaxialCorrTM01_we_vs_r.dat", ios::out);
  out3 << "# 2DACCELERATION_SVEA PROJECT OUTPUT FILE" << endl;
  out3 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out3 << "#" << endl;
  out3 << "# DESCRIPTION" << endl;
  out3 << "# Electric field energy density at beam waist versus radial coordinate" << endl;
  out3 << "#" << endl;
  myParaxialCorrBeam1.printparameters(out3);
  out3 << "#" << endl;
  out3 << "# r/lambda0  w_e(1)  w_e(5)  w_e(10)  w_e(20)  w_e(50)  w_e(200)" << endl;
  for (double r=0.0*lambda0; r<=6.0*lambda0; r+=lambda0/100.0) {
    vector<complex<double> > fields1 = myParaxialCorrBeam1.fields(r,0.0,0.0);
    vector<complex<double> > fields5 = myParaxialCorrBeam5.fields(r,0.0,0.0);
    vector<complex<double> > fields10 = myParaxialCorrBeam10.fields(r,0.0,0.0);
    vector<complex<double> > fields20 = myParaxialCorrBeam20.fields(r,0.0,0.0);
    vector<complex<double> > fields50 = myParaxialCorrBeam50.fields(r,0.0,0.0);
    vector<complex<double> > fields200 = myParaxialCorrBeam200.fields(r,0.0,0.0);
    double w_e1 = EPS0*(pow(abs(fields1[0]),2.0) + pow(abs(fields1[1]),2.0))/2.0;
    double w_e5 = EPS0*(pow(abs(fields5[0]),2.0) + pow(abs(fields5[1]),2.0))/2.0;
    double w_e10 = EPS0*(pow(abs(fields10[0]),2.0) + pow(abs(fields10[1]),2.0))/2.0;
    double w_e20 = EPS0*(pow(abs(fields20[0]),2.0) + pow(abs(fields20[1]),2.0))/2.0;
    double w_e50 = EPS0*(pow(abs(fields50[0]),2.0) + pow(abs(fields50[1]),2.0))/2.0;
    double w_e200 = EPS0*(pow(abs(fields200[0]),2.0) + pow(abs(fields200[1]),2.0))/2.0;
    out3 << r/lambda0 << "  " << w_e1 << "  " << w_e5 << "  " << w_e10;
    out3 << "  " << w_e20 << "  " << w_e50 << "  " << w_e200 << endl;
  } // end for
  // Close output file
  out3 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out3.close();*/


  // FEATURE TO BE CHARACTERIZED
  // Beam power verification

  // Definition of the pulsed beams to be used
  /*double lambda0 = 800e-9;
  double power = 1.0e14;
  double ka = 50.0;
  double duration = 10e-15;
  double phase = 0.0;
  CExactTM01 myExactBeam(lambda0,0.0,1e15,ka,duration,phase);
  CParaxialTM01 myParaxialBeam(lambda0,0.0,1e15,ka,duration,phase);
  CParaxial1stCorrTM01 myParaxialCorrBeam(lambda0,0.0,1e15,ka,duration,phase);

  // Calculate power by numerical integration
  double deltar = lambda0/10000.0;
  double sum1 = 0.0;
  double sum2 = 0.0;
  double sum3 = 0.0;
  for (double r=0.0*lambda0; r<=100.0*lambda0; r+=deltar) {
    vector<complex<double> > fields1 = myExactBeam.fields(r,0.0,0.0);
    vector<complex<double> > fields2 = myParaxialBeam.fields(r,0.0,0.0);
    vector<complex<double> > fields3 = myParaxialCorrBeam.fields(r,0.0,0.0);
    sum1 += PI*r*real(fields1[0]*conj(fields1[2]))*deltar;
    sum2 += PI*r*real(fields2[0]*conj(fields2[2]))*deltar;
    sum3 += PI*r*real(fields3[0]*conj(fields3[2]))*deltar;
  } // end for
  cout << "Average power of the exact beam: P = " << sum1 << " W" << endl;
  cout << "Average power of the paraxial beam: P = " << sum2 << " W" << endl;
  cout << "Average power of the corrected paraxial beam: P = " << sum3 << " W" << endl;*/

  /*// FEATURE TO BE CHARACTERIZED
  // Power radiated radially

  // Definition of the pulsed beams to be used
  double lambda0 = 800e-9;
  double power = 1.0e9;
  double ka = 1.0;
  double duration = 10e-15;
  double phase = 0.0;
  CExactTM01 myExactBeam(lambda0,power,0.0,ka,duration,phase);
  CParaxialTM01 myParaxialBeam(lambda0,power,0.0,ka,duration,phase);
  CParaxial1stCorrTM01 myParaxialCorrBeam(lambda0,power,0.0,ka,duration,phase);

  // Calculate power by numerical integration
  double r = 3*lambda0;
  double deltaz = lambda0/2000.0;
  double sum1 = 0.0;
  double sum2 = 0.0;
  double sum3 = 0.0;
  for (double z=-20.0*lambda0; z<=20.000001*lambda0; z+=deltaz) {
    vector<complex<double> > fields1 = myExactBeam.fields(r,z,0.0);
    vector<complex<double> > fields2 = myParaxialBeam.fields(r,z,0.0);
    vector<complex<double> > fields3 = myParaxialCorrBeam.fields(r,z,0.0);
    sum1 += PI*r*real(fields1[1]*conj(fields1[2]))*deltaz;
    sum2 += PI*r*real(fields2[1]*conj(fields2[2]))*deltaz;
    sum3 += PI*r*real(fields3[1]*conj(fields3[2]))*deltaz;
  } // end for
  cout << "Average radial power of the exact beam: P = " << sum1 << " W" << endl;
  cout << "Average radial power of the paraxial beam: P = " << sum2 << " W" << endl;
  cout << "Average radial power of the corrected paraxial beam: P = " << sum3 << " W" << endl;*/


  /*// FEATURE TO BE CHARACTERIZED
  // Electromagnetic field components at beam waist versus radial coordinate

  // Definition of the pulsed beams to be used
  double lambda0 = 800e-9;
  double w0 = 2e-6;
  double ka = 124.366;
  double zR = PI*pow(w0,2.0)/lambda0;
  CNonparaxialTM01 npTM01(lambda0,2.0e14,ka,277,0.0933*PI);
  CParaxialTM01 pTM01(lambda0,2.0e14,w0,10e-15,0.0933*PI);
  CBeam * beam1 = &npTM01;
  CBeam * beam2 = &pTM01;

  // Open output file
  ofstream out2;
  out2.open("./dat/test_beamcharacterization/pnpTM01_fields_vs_r.dat", ios::out);
  out2 << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out2 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out2 << "#" << endl;
  out2 << "# DESCRIPTION" << endl;
  out2 << "# Electromagnetic field components at beam waist at t=0" << endl;
  out2 << "#" << endl;
  beam1->printparameters(out2);
  out2 << "#" << endl;
  beam2->printparameters(out2);
  out2 << "#" << endl;
  out2 << "# r/lambda0  npTM01[Er,Ez,Hphi]  pTM01[Er,Ez,Hphi]" << endl;

  for (double r=-5.0*w0; r<=5.0*w0; r+=w0/10000.0) {
      vector<complex<double> > fields1 = beam1->fields(r,0.0,0.0);
      vector<complex<double> > fields2 = beam2->fields(r,0.0,0.0);
      out2 << r/w0 << "  " << real(fields1[0]) << "  " << real(fields1[1]) << "  " << real(fields1[2]);
      out2 << "  " << real(fields2[0]) << "  " << real(fields2[1]) << "  " << real(fields2[2]) << endl;
  } // end for

  // Close output file
  out2 << endl;
  out2 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out2.close();*/

  // FEATURE TO BE CHARACTERIZED
  // Electromagnetic field components at pulse peak versus r and z coordinates

  // Definition of the pulsed beams to be used
  double lambda0 = 800e-9;
  double power = 1.0e15;
  double ka = 500.0;
  double duration = 10e-15;
  double phase = PI;
  double k0 = 2.0*PI/lambda0;
  double zR = ka/k0;
  double w0 = sqrt(2.0*zR/k0);
  CExactTM01 myExactBeam(lambda0,power,0.0,ka,duration,phase);
  CParaxialTM01 myParaxialBeam(lambda0,power,0.0,ka,duration,phase);
  CParaxial1stCorrTM01 myParaxialCorrBeam(lambda0,power,0.0,ka,duration,phase);
  CBeam * beam1 = &myExactBeam;
  CBeam * beam2 = &myParaxialBeam;
  CBeam * beam3 = &myParaxialCorrBeam;

  // Open output file
  ofstream out3;
  out3.open("./dat/test_beamcharacterization/TM01_fields_rz.dat", ios::out);
  out3 << "# 2DACCELERATION_SVEA PROJECT OUTPUT FILE" << endl;
  out3 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out3 << "#" << endl;
  out3 << "# DESCRIPTION" << endl;
  out3 << "# Electromagnetic field components at pulse peak versus r and z coordinates" << endl;
  out3 << "#" << endl;
  beam1->printparameters(out3);
  out3 << "#" << endl;
  beam2->printparameters(out3);
  out3 << "#" << endl;
  beam3->printparameters(out3);
  out3 << "# z/zR  r/w0  exact[Er-eta0*Hphi,Ez]  par[Er-eta0*Hphi,Ez]  parcorr[Er-eta0*Hphi,Ez]" << endl;

  for (double z=-zR; z<=zR; z+=zR/1000.0) {
    for (double r=0; r<=2.0*w0; r+=w0/500.0) {
      vector<complex<double> > fields1 = beam1->fields(r,z,0);
      vector<complex<double> > fields2 = beam2->fields(r,z,0);
      vector<complex<double> > fields3 = beam3->fields(r,z,0);
      out3 << z/lambda0 << "  " << r/lambda0;
      out3 << "  " << -QE*real(fields1[0]-ETA0*fields1[2]) << "  " << -QE*real(fields1[1]);
      out3 << "  " << -QE*real(fields2[0]-ETA0*fields2[2]) << "  " << -QE*real(fields2[1]);
      out3 << "  " << -QE*real(fields3[0]-ETA0*fields3[2]) << "  " << -QE*real(fields3[1]) << endl;
    } // end for
    out3 << endl;
  } // end for

  // Close output file
  out3 << endl;
  out3 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out3.close();


  /*// FEATURE TO BE CHARACTERIZED
  // Electromagnetic field components versus r and z coordinates averaged over time

  // Definition of the pulsed beams to be used
  double lambda0 = 800e-9;
  double w0 = 2e-6;
  double ka = 124.366;
  double zR = PI*pow(w0,2.0)/lambda0;
  CNonparaxialTM01 npTM01(lambda0,1.0e12,ka,277,0.0);
  CParaxialTM01 pTM01(lambda0,1.0e12,w0,10e-15,0.0);
  CBeam * beam1 = &npTM01;
  CBeam * beam2 = &pTM01;

  // Open output file
  ofstream out4;
  out4.open("./dat/test_beamcharacterization/pnpTM01_fields_rz.dat", ios::out);
  out4 << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out4 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out4 << "#" << endl;
  out4 << "# DESCRIPTION" << endl;
  out4 << "# Electromagnetic field components versus r and z coordinates averaged over time" << endl;
  out4 << "#" << endl;
  beam1->printparameters(out4);
  out4 << "#" << endl;
  beam2->printparameters(out4);
  out4 << "#" << endl;
  out4 << "# r/w0  z/zR  npTM01[<E^2>]  pTM01[<E^2>]" << endl;

  // Results container
  int N = 201;
  vector<double> sample (N,0.0);
  vector<vector<double> > EsquaredP, EsquaredNP;
  EsquaredP.resize(N);
  EsquaredNP.resize(N);
  for (int n=0; n<N; n++) {
    EsquaredP[n] = sample;
    EsquaredNP[n] = sample;
  }

  for (double t=-100e-14; t<=100e-14; t+=1e-14) {
    cout << t << endl;
    for (int indexR=0; indexR<N; indexR+=1) {
      for (int indexZ=0; indexZ<N; indexZ+=1) {
        double r = indexR*w0/40.0;
        double z = -5.0*zR + indexZ*zR/20.0;
        vector<complex<double> > fields1 = beam1->fields(r,z,t);
        vector<complex<double> > fields2 = beam2->fields(r,z,t);
        EsquaredNP[indexR][indexZ] += pow(abs(fields1[0]),2.0) + pow(abs(fields1[1]),2.0);
        EsquaredP[indexR][indexZ] += pow(abs(fields2[0]),2.0) + pow(abs(fields2[1]),2.0);
      } // end for
    } // end for
  } // end for

  // Write results
  for (int indexR=0; indexR<N; indexR+=1) {
    for (int indexZ=0; indexZ<N; indexZ+=1) {
      double r = indexR*w0/40.0;
      double z = -5.0*zR + indexZ*zR/20.0;
      out4 << r/w0 << "  " << z/zR << "  " << EsquaredNP[indexR][indexZ] << "  " << EsquaredP[indexR][indexZ] << endl;
    } // end for
    out4 << endl;
  } // end for

  // Close output file
  out4 << endl;
  out4 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out4.close();*/


  /*// FEATURE TO BE CHARACTERIZED
  // Electromagnetic force components at pulse peak versus r and z coordinates

  // Definition of the pulsed beams to be used
  double lambda0 = 800e-9;
  double w0 = 2e-6;
  double ka = 124.366;
  double zR = PI*pow(w0,2.0)/lambda0;
  CNonparaxialTM01 npTM01(lambda0,2.0e14,ka,277,PI/2.0);
  CParaxialTM01 pTM01(lambda0,2.0e14,w0,10e-15,PI/2.0);
  CBeam * beam1 = &npTM01;
  CBeam * beam2 = &pTM01;

  // Open output file
  ofstream out5;
  out5.open("./dat/test_beamcharacterization/pnpTM01_force_rz.dat", ios::out);
  out5 << "# 2DACCELERATION PROJECT OUTPUT FILE" << endl;
  out5 << "# File created on " << GetDate() << " at " << GetTime() << endl;
  out5 << "#" << endl;
  out5 << "# DESCRIPTION" << endl;
  out5 << "# Electromagnetic force components at pulse peak versus r and z coordinates" << endl;
  out5 << "# Force computed for an electron with velocity 0.995c along the z axis" << endl;
  out5 << "#" << endl;
  beam1->printparameters(out5);
  out5 << "#" << endl;
  beam2->printparameters(out5);
  out5 << "#" << endl;
  out5 << "# z/zR  r/w0  npTM01[Fr,Fz]  pTM01[Fr,Fz,Hphi]" << endl;

  double vz = 0.999*C0;
  double gamma = pow(1.0-pow(vz/C0,2.0),-0.5);
  for (double z=-5.0*zR; z<=5.0*zR; z+=zR/20.0) {
    for (double r=0*w0; r<=5.0*w0; r+=w0/40.0) {
      vector<complex<double> > fields1 = beam1->fields(r,z,z/C0);
      vector<complex<double> > fields2 = beam2->fields(r,z,z/C0);
      double npFr = -QE*(real(fields1[0]) - vz*MU0*real(fields1[2]));
      double npFz = -QE*real(fields1[1]);
      double pFr = -QE*(real(fields2[0]) - vz*MU0*real(fields2[2]));
      double pFz = -QE*real(fields2[1]);
      out5 << z/zR << "  " << r/w0 << "  " << npFr << "  " << npFz << "  " << pFr << "  " << pFz << endl;
    } // end for
    out5 << endl;
  } // end for

  // Close output file
  out5 << endl;
  out5 << "# File closed on " << GetDate() << " at " << GetTime() << endl;
  out5.close();*/


  /*// FEATURE TO BE CHARACTERIZED
  // Maxwell's equation error parameter map

  // Definition of the pulsed beams to be used
  double lambda0 = 800e-9;
  double power = 1e12;
  double ka = 300;
  double duration = 15e-15;
  double phase = 0.0;
  CExactTM01 myExactBeam(lambda0,power,ka,duration,phase);
  CParaxialTM01 myParaxialBeam(lambda0,power,ka,duration,phase);
  CParaxial1stCorrTM01 myParaxialCorrBeam(lambda0,power,ka,duration,phase);
  CBeam * beam1 = &myExactBeam;
  CBeam * beam2 = &myParaxialBeam;
  CBeam * beam3 = &myParaxialCorrBeam;

  // Declaration of the output files
  ofstream out1,out2,out3;
  out1.open("./dat/test_beamcharacterization/ExactTM01_errormap_ka100.dat", ios::out);
  out2.open("./dat/test_beamcharacterization/ParaxialTM01_errormap_ka100.dat", ios::out);
  out3.open("./dat/test_beamcharacterization/ParaxialCorrTM01_errormap_ka100.dat", ios::out);
  double eps1 = beam1->errormap_fixedtime(0.0,10*lambda0,0.025*lambda0,10*lambda0,0.025*lambda0,1e-11,1e-11,1e-20,out1);
  double eps2 = beam2->errormap_fixedtime(0.0,10*lambda0,0.025*lambda0,10*lambda0,0.025*lambda0,1e-11,1e-11,1e-20,out2);
  double eps3 = beam3->errormap_fixedtime(0.0,10*lambda0,0.025*lambda0,10*lambda0,0.025*lambda0,1e-11,1e-11,1e-20,out3);*/

  return 0;
}
