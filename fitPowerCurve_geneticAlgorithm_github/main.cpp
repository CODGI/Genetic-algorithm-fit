#include "domainHandler.h"
#include "entity.h"
#include "material.h"
#include "position.h"
#include "source.h"
#include "fdtd.h"
#include <iostream>
#include <math.h>
#include "monitor.h"
#include <cstdio>
#include <string>

int main()
{

    double radiusWire = 0.1;
    double spacerWidth = 0.1;
    double metalWidth = 0.1;
    double Lx = 2;
    double Ly = 1;
    double dx = 0.001;
    double C = 0.5;
    double dt = C*dx;  
    double T = 5/3.33;
    int sample = 10;

    int numberOfPoints = 11;
    double width = 0.1;
    double centerWavelength = 0.38;

    domainHandler dh(Lx, Ly, dx, dx, dt);
    ellipse wire(radiusWire, radiusWire);
    rectangle spacer(Lx, spacerWidth);
    position pSource = position(0,0.1);
    position pWire = position(0,-0.1);
    position pSpacer = position(0,-0.1-radiusWire-spacerWidth/2);


    material Al = material("Al",1,1,0,0);
    Al.adddDrudePole(3.33E-15*2*3.1415*3800E12, 3.33E-15*1E15);


    //cout << "Placing Wire." << endl;
    //dh.addEntity(wire, pWire, material("ZnO",6.5,1,0,0), false, true);
    //cout << "Placing Spacer." << endl;
    //dh.addEntity(spacer, pWire, Al, false,true);
    //cout << "Creating dipole source." << endl;
    //dipoleSource s(Lx, Ly, dx, dx, p1, omega, orientation, amplitude, T/100);
    dh.addPML(30,0.2);
    //dh.printDomain();

    string direction = "x";
    double length = 1;
    double omega = 2*M_PI/0.38;
    double amplitude = 1;
    vector<int> orientation{0,0,1,1};
    gaussianSource s(Lx, Ly, dx, dx, pSource, length, direction, omega, orientation, amplitude, 10*T, 0.4, 0.1);
    //lineSource s(Lx, Ly, dx, dx, pSource, length, direction, omega, orientation, amplitude, 10*T);

    cout << "Creating monitor." << endl;
    freopen( "I_planeWave.csv", "w", stderr );
    monitor m(Lx, Ly); 

    cout << "Creating fdtd."<< endl;
    fdtd f(dh, dt, T, sample, &s, &m);

    f.run();
    
    return 0;
}
