#pragma once

#include <vector>
#include <boost/array.hpp>
#include <boost/numeric/ublas/vector.hpp>
using namespace std;

typedef boost::numeric::ublas::vector< double> state_type;

const double pi = 3.141592653589793;
const double hbar = 1.0545718001391127e-34;


class model 
{
    public:
        model(vector<double>);
        template< class State > void operator()(const State&, State&, const double);
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
};

template< class State >
void model::operator()(const State &x, State &dxdt, const double t){
    //cout << t << endl;
    double P = pulses(t)*1E10;
    double P_n = P*L*d/(hbar*(2*pi*3E8/350E-9)*L*pi*(d/2)*(d/2))*(n_abs-x[0])/n_abs;
    double n_nonrad = -A_auger*pow(x[0],3)-x[0]/tau_nonrad;
    double n_sp = -F*x[0]/tau_rad;
    double n_st = -(3E8/16)*g*log(x[0]/n_trans)*G*F*x[1];
    dxdt[0] = P_n+n_nonrad+n_sp+n_st;
    dxdt[1] = -beta*n_sp - n_st - f_p*x[1];
    dxdt[2] = A*x[1];
}