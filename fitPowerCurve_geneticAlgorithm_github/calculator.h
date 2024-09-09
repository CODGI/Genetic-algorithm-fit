#pragma once
#include <vector>
using namespace std;

class calculator{
    public:
        int id;
        calculator(int, vector<vector<double>>, double);
        void setParameters(vector<double>);
        double run(bool);
    private:
        double T;
        vector<double> params;
        vector<vector<double>> data;
};