#pragma once
#include <vector>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
using namespace std;
using namespace boost::numeric::odeint;

typedef boost::array< double , 2> state_type;
typedef runge_kutta_dopri5< state_type > stepper_type;


class integralObserver{
    public:
        double integral;
        calculator(int, vector<vector<double>>, double);
        void setParameters(vector<double>);
        double run();
    private:
        double T;
        vector<double> params;
        vector<vector<double>> data;
};