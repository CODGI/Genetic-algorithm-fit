#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <mpi.h>
#include "calculator.h"
#include "optimizer_evolve.h"
#include "model.h"
#include <float.h>
using namespace std;

double T = 100E-9;
double g = 700000;
double n_trans = 2e+24;
double tau1 = 10E-9;
double P1 = 4;
double T1 = 7E-9;
double d = 200*1E-9;
double L = 6E-6;
double tau_nonrad = 700E-12;
double tau_rad = 700E-12;
double F =	45.8;
double G = 3.4;
double f_p = 1.42968e+13;
double beta = 0.94;
double A = 8.57383e-14;
double n_abs = 5.61131e+37;
double A_auger = 1E-57;
vector<double> paramsStart{T, g, n_trans, tau1, P1, T1, d, L, tau_nonrad, tau_rad, F, G, f_p, beta, A, n_abs, A_auger};


vector<vector<double> > readData(string filename) {
    vector<vector<double> > data;
	vector<double> row;
	string line, word;

    fstream file (filename, ios::in);
	if(file.is_open())
	{
		while(getline(file, line))
		{
			row.clear();
 
			stringstream str(line);
 
			while(getline(str, word, ' ')) {
				row.push_back(stod(word));
			}
			data.push_back(row);
		}
	}
	else
		cout<<"Could not open the file\n";

	return data;
}


int main() {
	MPI_Init(NULL, NULL);
	int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	vector<vector<double> > data = readData("/beegfs/ni86did/fitPowerCurve_geneticAlgorithm/Data_NW4_onAg_pumpPERP_S5.txt");

	int N_generations = 20000;

	if (world_rank == 0) {
		int numberOfSurvivors = 60;
		int numberToPrint = 10;
		double mutationBound = 2;	
		vector<int> dontVary{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 16};
		double stepLength = 0.1;
		double learningrate = 0.1;
		optimizer_evolve oe(model::Nparameters, world_size-1, numberOfSurvivors, mutationBound, paramsStart, data, T, stepLength, learningrate);
		vector<int> indBounded{13}; vector<double> low{0.1}; vector<double> high{0.99};
		oe.setBounds(indBounded, low, high);
		oe.set_dontVary(dontVary);
		oe.initialize();
		for (int n=0; n<N_generations; n++) {
			vector<vector<double>> allParams = oe.getAllParameters();
			for (int k=1; k<world_size;k++) {
				vector<double> sendLine = allParams[k-1];
				MPI_Send(&sendLine[0], model::Nparameters, MPI_DOUBLE, k, n, MPI_COMM_WORLD);
			}
			vector<double> allErrors(world_size-1,0.0);
			for (int k=1; k<world_size;k++) {
				MPI_Recv(&allErrors[k-1], 1, MPI_DOUBLE, k, n, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			oe.set_Errors(allErrors);
			double lowestError = oe.getLowestError();
			cout << "Lowest error " << lowestError << endl;
			double largestError = oe.getLargestError();
			cout << "Largest error " << largestError << endl;
			if (n%100 == 99) oe.printSurvivors(numberToPrint);
			if (lowestError < DBL_MAX) {
				oe.setMutationBound(lowestError);
				oe.removeUnfit();
				oe.generateChildren();
				oe.enforceBounds();
			} else {
				oe.initialize();
			}
		}
		oe.printSurvivors(numberOfSurvivors);
	} else  {
		calculator c(world_rank, data, T);
		vector<double> params = vector<double>(model::Nparameters,0.0);
		for(int n=0; n<N_generations; n++) {
			MPI_Recv(&params[0], model::Nparameters, MPI_DOUBLE, 0, n,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			c.setParameters(params);
			double error = c.run(false);
			MPI_Send(&error, 1, MPI_DOUBLE, 0, n, MPI_COMM_WORLD);
		}
	}

	
    MPI_Finalize();
	return 0;
}