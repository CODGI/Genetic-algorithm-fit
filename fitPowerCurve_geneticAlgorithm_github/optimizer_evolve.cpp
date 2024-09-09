#include <vector>
#include <numeric>
#include <algorithm>
#include "optimizer_evolve.h"
#include "calculator.h"
#include <math.h>
#include <random>
#include <float.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "calculator.h"
using namespace std;

vector<int> getLargestIndexes(vector<double> errors, int N) {
	vector<double> indices(errors.size());
	std::iota(indices.begin(), indices.end(), 0);
	std::partial_sort(indices.begin(), indices.begin()+N, indices.end(),
						[&errors](int i, int j) {
							return errors[i] > errors[j];
							});
	return vector<int>(indices.begin(), indices.begin()+N);
}

vector<int> getSmallestIndexes(vector<double> errors, int N) {
	vector<double> indices(errors.size());
	std::iota(indices.begin(), indices.end(), 0);
	std::partial_sort(indices.begin(), indices.begin()+N, indices.end(),
						[&errors](int i, int j) {
							return errors[i] < errors[j];
							});
	return vector<int>(indices.begin(), indices.begin()+N);
}

optimizer_evolve::optimizer_evolve(int _Nparameters, int _populationSize, int _survivors, double _mutationBound, vector<double> _paramsStart,
                                    vector<vector<double>> _data, double _T, double _stepLength, double _lr) {
    Nparameters = _Nparameters;
    populationSize = _populationSize;
    survivors = _survivors;
    mutationBound = _mutationBound;
    paramsStart = _paramsStart;
    allParams = vector<vector<double>>(populationSize,vector<double>(Nparameters,0.0));
    allErrors = vector<double>(populationSize,0.0);
    data = _data;
    T = _T;
    stepLength = _stepLength;
    lr = _lr;
}

void optimizer_evolve::gradientDescend(int numberOfSteps) {
	vector<int> lowestInd =  getSmallestIndexes(allErrors, 1);
    paramsStart = allParams[lowestInd[0]];
    cout << "Starting gradient descend" << endl;
    for (int i=0;i<numberOfSteps; i++) {
        cout << "Step " << i <<  endl;
        cout << "Learning rate : " << lr << endl;
        cout << "Step length: " << stepLength << endl;
        calculator c1(0, data, T);
        c1.setParameters(paramsStart);
        double e1 = c1.run(false);
        cout << "Error: "<< e1 << endl;
        vector<double> newParams;
        vector<double> newParams2;
        for(int j=0; j<Nparameters; j++) {
            if (std::find(dontVary.begin(),dontVary.end(),j) != std::end(dontVary))  {
			    newParams.push_back(paramsStart[j]);
                newParams2.push_back(paramsStart[j]);
            } else {
                calculator c2(0, data, T);
                vector<double> params2 = paramsStart;
                params2[j] = (1.0 + stepLength)*paramsStart[j];
                double step = stepLength*paramsStart[j];
                c2.setParameters(params2);
                double e2 = c2.run(false);
                newParams.push_back(paramsStart[j]-lr*((e2-e1)/step));
                newParams2.push_back(paramsStart[j]-(lr/2)*((e2-e1)/step));
                cout << "e2-e1: " << e2-e1 << endl;
            }
        }
        calculator c(0, data, T);
        c.setParameters(newParams);
        double errorNew = c.run(false);
        c.setParameters(newParams2);
        double errorNew2 = c.run(false);
        if (errorNew2 < errorNew) {
            newParams = newParams2;
            errorNew = errorNew2;
            if (lr>1E-8) lr /= 2;
        }
        if (errorNew<e1) {
            paramsStart = newParams;
        }
        else if (stepLength>1E-8) stepLength /= 2;
        
    }
    allParams[lowestInd[0]] = paramsStart;
    calculator c(0, data, T);
    c.setParameters(paramsStart);
    double errorNew = c.run(false);
    allErrors[lowestInd[0]] = errorNew;
}



void optimizer_evolve::setBounds(vector<int> indexes, vector<double> _lowerBounds, vector<double> _upperBounds) {
    indexesBound = indexes;
    lowerBounds = _lowerBounds;
    upperBounds = _upperBounds;
}

vector<double> optimizer_evolve::randomMutation(vector<double> params) {
    vector<double> mutatedParams;
	std::random_device rd;
	std::default_random_engine eng(rd());
	std::uniform_real_distribution<double> distr(-1, 1);
    srand(time(NULL));
	for (int i=0; i<params.size(); i++) {
		if (std::find(dontVary.begin(),dontVary.end(),i) != std::end(dontVary))  {
			mutatedParams.push_back(params[i]);
		}
		else {
            int mutate = rand()%2;
            double value = 0;
            if (mutate==1) {
                value = params[i]*(1+(distr(eng)/mutationBound));
            } else {
                value = params[i];
            }
            mutatedParams.push_back(value);
		}

	}
	return mutatedParams;
}

void optimizer_evolve::enforceBounds() {
    for (int i=0;i<indexesBound.size();i++){
        int boundIndex = indexesBound[i];
        for(int j=0;j<populationSize;j++) {
            if (allParams[j][boundIndex] < lowerBounds[i]) allParams[j][boundIndex] = lowerBounds[i];
            if (allParams[j][boundIndex] > upperBounds[i]) allParams[j][boundIndex] = upperBounds[i];
        }
    }
}

void optimizer_evolve::set_dontVary(vector<int> _dontVary) {
    dontVary = _dontVary;
}

vector<double> optimizer_evolve::pair(int j, int k) {
	vector<double> child(Nparameters,0.0);
    srand(time(NULL));
    for (int i=0; i<Nparameters; i++) {
        int number = rand()%2;
        if (number == 0) {
            child[i] = allParams[j][i];
        } else {
            child[i] = allParams[k][i];
        }
    }
    return child;
}

vector<vector<double>> optimizer_evolve::getAllParameters() {
    return allParams;
}

void optimizer_evolve::initialize() {
    allParams[0] = paramsStart;
    int k = 1;
    for (int i=1; i<populationSize; i++) {
        allParams[i] = randomMutation(paramsStart);
    }
}

void optimizer_evolve::set_Errors(vector<double> _allErrors) {
    allErrors = _allErrors;
}

void optimizer_evolve::outputAllErrors(){
    for (int i=0; i< populationSize; i++) {
        cout << allErrors[i] << endl;
    }
}

void optimizer_evolve::outputAllParameters(){
    for (int i=0; i< populationSize; i++) {
        cout << "Parameterset " << i << ", error: " << allErrors[i] << endl;
        for (int j=0; j < Nparameters; j++) {
            cout <<  allParams[i][j] << endl;
        }
    }
}

void optimizer_evolve::removeUnfit() {
    vector<int> lowestInd =  getSmallestIndexes(allErrors, survivors);
	vector<vector<double>> newParameters = vector<vector<double>>(populationSize,vector<double>(Nparameters,0.0));
	int lowestIndexNotInfinite = -1;
	for (int m=0;m<lowestInd.size();m++) {
		if (!(allErrors[lowestInd[m]] == DBL_MAX)) {
			lowestIndexNotInfinite = m;
		} else {
            break;
        }
	}
	for (int i=0; i<lowestInd.size(); i++) {
		if (allErrors[lowestInd[i]] == DBL_MAX) {
			if (lowestIndexNotInfinite>-1) {
				newParameters[i] = randomMutation(allParams[lowestInd[0]]);
			} else {
				newParameters[i] = randomMutation(paramsStart);
			}
				}
		else {
			newParameters[i] = allParams[lowestInd[i]];
		}
	}
    allParams = newParameters;
}


void optimizer_evolve::generateChildren() {
    srand(time(NULL));
    for (int i=survivors+1; i < populationSize; i++) {
        int ind1 = rand()%survivors;
        int ind2 = rand()%survivors;
        while (ind2 == ind1) {
            ind2 = rand()%survivors;
        }
        vector<double> child = pair(ind1, ind2);
        int mutate = rand()%2;
        if (mutate == 0) 
        {
            allParams[i] = randomMutation(child);
        } else {
            allParams[i] = child;
        }
    }
}

double optimizer_evolve::getLowestError() {
    vector<int> lowestInd =  getSmallestIndexes(allErrors, 1);
    return allErrors[lowestInd[0]];
}

double optimizer_evolve::getLargestError() {
    vector<int> largestInd =  getLargestIndexes(allErrors, 1);
    return allErrors[largestInd[0]];
}

void optimizer_evolve::printSurvivors(int N) {
    vector<int> bestPrameters =  getSmallestIndexes(allErrors, survivors);
		for(int i=0; i< N; i++) {
			cout << "-----" << endl;
            cout << "Survivor " << i << ", error: " << allErrors[bestPrameters[i]] << endl;
            calculator c(0, data, T);
            c.setParameters(allParams[bestPrameters[i]]);
            c.run(true);
			for (int j=0;j<Nparameters;j++) {
				cout << allParams[bestPrameters[i]][j] << ", " << endl;
			}
			cout << endl;
		}
}

void optimizer_evolve::increaseMutationBound(double upperBound) {
    if (mutationBound < upperBound) {
        mutationBound *= 2;
    }
    cout << "Mutation bound: " << mutationBound << endl;
}

void optimizer_evolve::decreaseMutationBound(double lowerBound) {
    if (mutationBound > lowerBound) {
        mutationBound /= 2;
    }
    cout << "Mutation bound: " << mutationBound << endl;
}

void optimizer_evolve::setMutationBound(double error) {
    if (error < 0.5) {
        mutationBound = 1/error;
    } else {
        mutationBound = 2;
    }
    cout << "Mutation bound: " << mutationBound << endl;
}


