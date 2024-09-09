#pragma once

#include <vector>
using namespace std;

/*
const double pi = 3.141592653589793;
const double hbar = 1.0545718001391127e-34;
*/

class model_jacobi 
{
    public:
        model_jacobi(vector<double>);
        template< class State , class Matrix > void operator()( const State & , Matrix & , const double & , State & );
        static const int Nparameters{17};
    private:
        double T;           //0
        double g;           //1
        double n_trans;     //2
        double tau1;        //3
        double P1;          //4
        double T1;          //5
        double d;           //10
        double L;           //11
        double tau_nonrad;  //12
        double tau_rad;     //13
        double F;           //14
        double G;           //15
        double f_p;         //16
        double beta;        //17
        double A;           //18
        double n_abs;       //19
        double A_auger;     //20
        double pulses(double);
        double pulses2(double);
};

 template< class State , class Matrix > 
 void model_jacobi::operator()( const State &x , Matrix &J , const double &t , State &dfdt ) {
    J(0,0) = -pulses2(t)*1E10*L*d/(n_abs*hbar*(2*pi*3E8/350E-9)*L*pi*(d/2)*(d/2))-((1/tau_nonrad)+(3*A_auger*pow(x[0],2))+(F/tau_rad) + ((3E8/16)*g*G*F*x[1]/(x[0])));
    J(0,1) = -(3E8/16)*g*G*F*log(x[0]/n_trans);
    J(0,2) = 0;
    J(1,0) = (F*beta/tau_rad) + (3E8/16)*g*G*F*x[1]/(x[0]);
    J(1,1) = (3E8/16)*g*G*F*log(x[0]/n_trans) - f_p;
    J(1,2) = 0;
    J(2,0) = 0;
    J(2,1) = A;
    J(2,2) = 0;
    dfdt[0] = pulses(t)*1E10*L*d/(hbar*(2*pi*3E8/350E-9)*L*pi*(d/2)*(d/2))*(n_abs-x[0])/n_abs;
    dfdt[1] = 0;
    dfdt[2] = 0;
    /*
    double P = pulses(t)*1E10;
    double P_n = P*L*d/(hbar*(2*pi*3E8/350E-9)*L*pi*(d/2)*(d/2));
    double n_nonrad = -x[0]/tau_nonrad;
    double n_sp = -G*F*beta*x[0]/tau_rad;
    double n_st = -(3E8/16)*g*log(x[0]/n_trans)*G*F*x[1];
    dxdt[0] = P_n+n_nonrad+n_sp+n_st;
    dxdt[1] = -n_sp - n_st - f_p*x[1];
    dxdt[2] = A*x[1];
    */
}