#include "model.h"
#include <boost/array.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <iostream>
using namespace std;

 model::model(vector<double> params) {
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

 double model::pulses(double t) {
    return (this->P1)*exp(-pow((t-(this->tau1))/(this->T1),2));
}

/*
template< class State >
void model::operator()(const State &x, State &dxdt, const double t){
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
