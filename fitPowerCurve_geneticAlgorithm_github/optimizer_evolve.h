#pragma once

#include <vector>
using namespace std;

class optimizer_evolve {
    public:
        optimizer_evolve(int, int, int, double, vector<double>, vector<vector<double>>, double, double, double);
        void set_dontVary(vector<int>);
        void set_Errors(vector<double>);
        vector<vector<double>> getAllParameters();
        void initialize();
        void removeUnfit();
        void generateChildren();
        double getLowestError();
        double getLargestError();
        void printSurvivors(int);
        void gradientDescend(int);
        void newton(int);
        void increaseMutationBound(double);
        void decreaseMutationBound(double);
        void setMutationBound(double);
        void setBounds(vector<int>, vector<double>, vector<double>);
        void enforceBounds();
    private:
        double T;
        vector<vector<double>> data;
        vector<double> paramsStart;
        vector<vector<double>> allParams;
        vector<double> allErrors;
        int populationSize;
        int survivors;
        double mutationBound;
        int Nparameters;
        double stepLength;
        double lr;
        vector<int> dontVary;
        vector<int> indexesBound;
        vector<double> lowerBounds;
        vector<double> upperBounds;
        vector<double> randomMutation(vector<double>);
        vector<double> pair( int, int);
        void outputAllErrors();
        void outputAllParameters();
};