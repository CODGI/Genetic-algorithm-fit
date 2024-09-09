#include "model_jacobi.h"
#include <boost/array.hpp>
#include <iostream>
using namespace std;

 model_jacobi::model_jacobi(vector<double> params) {
    this->T = params[0];
    this->g = params[1];
    this->n_trans = params[2];
    this->tau1 = params[3];
    this->P1 = params[4];
    this->T1 = params[5];
    this->d = params[6];
    this->L = params[7]; 
    this->tau_nonrad = params[8];
    this->tau_rad = params[9]; 
    this->F = params[10];
    this->G = params[11];
    this->f_p = params[12];
    this->beta = params[13];
    this->A = params[14];
    this->n_abs = params[15];
    this->A_auger = params[16];
 }

 double model_jacobi::pulses(double t) {
    return (-2*((t-(this->tau1))/pow(this->T1,2))*(this->P1)*exp(-pow((t-(this->tau1))/(this->T1),2)));
}

 double model_jacobi::pulses2(double t) {
    return (this->P1)*exp(-pow((t-(this->tau1))/(this->T1),2));
}

/*
 template< class State , class Matrix > 
 void model_jacobi::operator()( const State &x , Matrix &J , const double &t , State &dfdt ) {
    J(0,0) = -((1/tau_nonrad)+(G*F*beta/tau_rad) + ((3E8/16)*g*G*F*x[1]/x[0]));
    J(0,1) = -(3E8/16)*g*log(x[0]/n_trans)*G*F;
    J(0,2) = 0;
    J(1,0) = (G*F*beta/tau_rad) + ((3E8/16)*g*G*F*x[1]/x[0]);
    J(1,1) = (3E8/16)*g*log(x[0]/n_trans)*G*F - f_p;
    J(1,2) = 0;
    J(2,0) = 0;
    J(2,1) = 1;
    J(2,2) = 0;
    dfdt[0] = pulses(t)*1E10*L*d/(hbar*(2*pi*3E8/350E-9)*L*pi*(d/2)*(d/2));
    dfdt[1] = 0;
    dfdt[2] = 0;
    
    double P = pulses(t)*1E10;
    double P_n = P*L*d/(hbar*(2*pi*3E8/350E-9)*L*pi*(d/2)*(d/2));
    double n_nonrad = -x[0]/tau_nonrad;
    double n_sp = -G*F*beta*x[0]/tau_rad;
    double n_st = -(3E8/16)*g*log(x[0]/n_trans)*G*F*x[1];
    dxdt[0] = P_n+n_nonrad+n_sp+n_st;
    dxdt[1] = -n_sp - n_st - f_p*x[1];
    dxdt[2] = A*x[1];
    
}


*/

