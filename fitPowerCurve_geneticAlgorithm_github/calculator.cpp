#include "calculator.h"
#include "model.h"
#include "model_jacobi.h"
#include <vector>
/*
#include <boost/array.hpp>
*/
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <math.h>
#include <float.h>

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::numeric::ublas::vector< double> state_type;
typedef rosenbrock4< state_type > stepper_type;

calculator::calculator(int _id, vector<vector<double>> _data, double _T) {
    id = _id;
    data = _data;
    params = vector<double>(model::Nparameters,0.0);
    T = _T;
}

void calculator::setParameters(vector<double> _params) {
    params = _params;
}


double calculator::run(bool showResult) {
    vector<double> errors;
	for (int i=0; i < data.size(); i++) {
		params[4] = data[i][0]/1000;
		model m(params);
		model_jacobi mj(params);
    	state_type x(3);x[0] = 1E22; x[1] = 0.0; x[2] = 0.0;
		//integrate_adaptive(make_controlled< rosenbrock4< double > >( 1.0e-12 , 1.0e-12 ), make_pair(m, mj), x, 0.0, T, 1E-14);
		integrate_const(rosenbrock4<double>(), make_pair(m, mj), x, 0.0, T, 10E-12);
		double value = x[2];
		//cout << "value" << value << endl;
		if (showResult) {
			cout << data[i][0] <<", "<< data[i][1] << ", " << value  << endl;
		}
		errors.push_back(pow(value-data[i][1],2)/pow(data[i][1],2));
	}
	double e = 0;
	for (int i=0; i < errors.size(); i++) {
		e += errors[i];
	}
	e = sqrt(e);
	if (std::isnan(e)) return DBL_MAX;
    return e;
}